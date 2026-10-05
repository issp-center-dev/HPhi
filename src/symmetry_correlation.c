#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "Common.h"
#include "symmetry_basis.h"
#include "symmetry_checked.h"
#include "symmetry_correlation.h"
#include "symmetry_directory.h"
#include "symmetry_matvec_plan.h"
#include "symmetry_memory_policy.h"
#include "symmetry_mpi_exchange.h"
#include "symmetry_terms.h"
#include "symmetry_vector_halo.h"
#ifdef MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

struct CorrelationOrbit {
  unsigned int factors;
  size_t member_count;
  size_t offdiag_count;
  int *members;            /* member_count * 4 * factors, sorted; members[0..] is the orbit key */
  unsigned char *offdiag;  /* member_count flags: 1 when the product can change the state */
};

struct CorrelationOrbitTable {
  size_t orbit_count;
  struct CorrelationOrbit *orbits;
  size_t *op_to_orbit;
  uint64_t offdiag_total;
  uint64_t bytes;          /* live allocation of the table after construction */
  uint64_t scratch_peak;   /* largest construction scratch, freed before return */
};

struct CorrelationTransition {
  uint32_t orbit;
  uint32_t reserved;
  unsigned long int key;        /* replicated: global index; distributed: representative state */
  double complex coefficient;   /* a_r * sign * phase / norm_r (replicated: also * norm_target) */
};

struct CorrelationContext {
  const struct BindStruct *X;
  const double complex *vec;
  const struct CorrelationOrbitTable *table;
  struct SymmetryRepresentativeDirectory *directory;   /* owned here; distributed layout with off-diagonal members only */
  uint64_t directory_bytes;
  uint64_t fixed_bytes;                                 /* table + directory + accumulators + thread partials */
  struct SymmetryMemoryPolicy policy;
  unsigned long int block_rows;
  unsigned int threads;
  int mpi_active;
  int replicated;
  int warning_emitted;
  struct SymmetryCorrelationStats stats;
};

/* ---- small helpers ----------------------------------------------------- */

static unsigned int team_size(void)
{
#ifdef _OPENMP
  int n = omp_get_max_threads();
  return n > 0 ? (unsigned int)n : 1U;
#else
  return 1U;
#endif
}

static int checked_add64(uint64_t a, uint64_t b, uint64_t *out)
{
  if (a > UINT64_MAX - b) return -1;
  *out = a + b;
  return 0;
}

static int checked_mul64(uint64_t a, uint64_t b, uint64_t *out)
{
  if (a != 0U && b > UINT64_MAX / a) return -1;
  *out = a * b;
  return 0;
}

static int compare_tuple(const int *a, const int *b, size_t len)
{
  size_t i;
  for (i = 0; i < len; i++) if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
  return 0;
}

/* Same criterion as emit() in symmetry_terms.c: the occupation change at
 * every out-orbital vanishes. */
static int term_is_diagonal(unsigned int factors, const int *index)
{
  unsigned int i, j;
  for (i = 0; i < factors; ++i) {
    int balance = 0;
    for (j = 0; j < factors; ++j) {
      if (index[4*i] == index[4*j] && index[4*i+1] == index[4*j+1]) ++balance;
      if (index[4*i] == index[4*j+2] && index[4*i+1] == index[4*j+3]) --balance;
    }
    if (balance != 0) return 0;
  }
  return 1;
}

static int validate_operator(const struct DefineList *def,
                             const struct SymmetryCorrelationOperator *op)
{
  unsigned int f, nint;
  int spin_max = def->iCalcModel == SpinlessFermion ? 0 : 1;
  if (op->factors < 1U || op->factors > UINT_MAX / 4U || op->index == NULL) return -1;
  nint = 4U * op->factors;
  for (f = 0; f < nint; f += 2U)
    if (op->index[f] < 0 || (unsigned int)op->index[f] >= def->Nsite ||
        op->index[f+1] < 0 || op->index[f+1] > spin_max) return -1;
  if (def->iCalcModel == Spin)
    for (f = 0; f < op->factors; ++f)
      if (op->index[4*f] != op->index[4*f+2]) return -1;
  return 0;
}

static void permute_tuple(const struct DefineList *def, unsigned int g,
                          unsigned int factors, const int *in, int *out)
{
  unsigned int f;
  for (f = 0; f < factors; ++f) {
    out[4*f]   = def->SymTrans[g][in[4*f]];
    out[4*f+1] = in[4*f+1];
    out[4*f+2] = def->SymTrans[g][in[4*f+2]];
    out[4*f+3] = in[4*f+3];
  }
}

static void free_orbit_table(struct CorrelationOrbitTable *table)
{
  size_t k;
  if (table->orbits != NULL)
    for (k = 0; k < table->orbit_count; k++) { free(table->orbits[k].members); free(table->orbits[k].offdiag); }
  free(table->orbits);
  free(table->op_to_orbit);
  memset(table, 0, sizeof(*table));
}

/* ---- orbit table ------------------------------------------------------- */

static int build_orbit_table(const struct DefineList *def,
                             const struct SymmetryCorrelationOperator *ops,
                             size_t count, struct CorrelationOrbitTable *table)
{
  size_t t, k;
  int *sorted = NULL, *tuple = NULL;
  size_t sorted_capacity = 0;
  int status = -1;
  memset(table, 0, sizeof(*table));
  if (def == NULL || ops == NULL || def->NSymTrans == 0U || def->SymTrans == NULL) return -1;
  table->orbits = (struct CorrelationOrbit *)calloc(count, sizeof(*table->orbits));
  table->op_to_orbit = (size_t *)calloc(count, sizeof(*table->op_to_orbit));
  if (table->orbits == NULL || table->op_to_orbit == NULL) goto done;
  table->bytes = (uint64_t)count * (sizeof(*table->orbits) + sizeof(*table->op_to_orbit));
  for (t = 0; t < count; t++) {
    unsigned int factors = ops[t].factors;
    size_t len, n = 0, pos, need;
    unsigned int g;
    uint64_t scratch;
    if (validate_operator(def, &ops[t]) != 0) goto done;
    len = 4U * (size_t)factors;
    if (SymmetryCheckedSizeMul(len, (size_t)def->NSymTrans, &need) != 0) goto done;
    if (need > sorted_capacity) {
      free(sorted); free(tuple);
      sorted = (int *)malloc(need * sizeof(int));
      tuple = (int *)malloc(len * sizeof(int));
      if (sorted == NULL || tuple == NULL) goto done;
      sorted_capacity = need;
    }
    scratch = ((uint64_t)sorted_capacity + (uint64_t)len) * sizeof(int);
    if (scratch > table->scratch_peak) table->scratch_peak = scratch;
    for (g = 0; g < def->NSymTrans; g++) {
      int c = 1;
      permute_tuple(def, g, factors, ops[t].index, tuple);
      for (pos = 0; pos < n; pos++) {
        c = compare_tuple(tuple, sorted + pos * len, len);
        if (c <= 0) break;
      }
      if (pos < n && c == 0) continue;               /* duplicate member */
      memmove(sorted + (pos + 1) * len, sorted + pos * len, (n - pos) * len * sizeof(int));
      memcpy(sorted + pos * len, tuple, len * sizeof(int));
      n++;
    }
    for (k = 0; k < table->orbit_count; k++)
      if (table->orbits[k].factors == factors &&
          compare_tuple(table->orbits[k].members, sorted, len) == 0) break;
    if (k == table->orbit_count) {
      struct CorrelationOrbit *orbit = &table->orbits[k];
      size_t m;
      orbit->factors = factors;
      orbit->member_count = n;
      orbit->members = (int *)malloc(n * len * sizeof(int));
      orbit->offdiag = (unsigned char *)calloc(n, sizeof(unsigned char));
      table->orbit_count++;
      if (orbit->members == NULL || orbit->offdiag == NULL) goto done;
      table->bytes += (uint64_t)n * len * sizeof(int) + (uint64_t)n;
      memcpy(orbit->members, sorted, n * len * sizeof(int));
      for (m = 0; m < n; m++) {
        orbit->offdiag[m] = term_is_diagonal(factors, orbit->members + m * len) ? 0U : 1U;
        orbit->offdiag_count += orbit->offdiag[m];
      }
      table->offdiag_total += orbit->offdiag_count;
    }
    table->op_to_orbit[t] = k;
  }
  status = 0;
done:
  free(sorted); free(tuple);
  if (status != 0) free_orbit_table(table);
  return status;
}

int SymmetryCorrelationOrbitCount(const struct DefineList *def,
                                  const struct SymmetryCorrelationOperator *ops,
                                  size_t count, size_t *orbit_count, size_t *member_total)
{
  struct CorrelationOrbitTable table;
  size_t k, members = 0;
  if (orbit_count == NULL || member_total == NULL) return -1;
  if (count == 0) { *orbit_count = 0; *member_total = 0; return 0; }
  if (build_orbit_table(def, ops, count, &table) != 0) return -1;
  for (k = 0; k < table.orbit_count; k++) members += table.orbits[k].member_count;
  *orbit_count = table.orbit_count;
  *member_total = members;
  free_orbit_table(&table);
  return 0;
}

/* ---- MPI helpers --------------------------------------------------------- */

static int agree_block_rows(int mpi_active, unsigned long int rows)
{
#ifdef MPI
  if (mpi_active != FALSE && nproc > 1) {
    unsigned long int rows_min, rows_max;
    if (MPI_Allreduce(&rows, &rows_min, 1, MPI_UNSIGNED_LONG, MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&rows, &rows_max, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS ||
        rows_min != rows_max) return -1;
  }
#else
  (void)mpi_active; (void)rows;
#endif
  return 0;
}

static int max_wave_count(int mpi_active, unsigned long long local, unsigned long long *max)
{
  *max = local;
#ifdef MPI
  if (mpi_active != FALSE && nproc > 1) {
    if (MPI_Allreduce(&local, max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS) return -1;
  }
#else
  (void)mpi_active;
#endif
  return 0;
}

static int max_fixed_bytes(int mpi_active, uint64_t local, uint64_t *max)
{
  *max = local;
#ifdef MPI
  if (mpi_active != FALSE && nproc > 1) {
    if (MPI_Allreduce(&local, max, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS) return -1;
  }
#else
  (void)mpi_active;
#endif
  return 0;
}

static int allreduce_sum(int mpi_active, double complex *acc, size_t n)
{
#ifdef MPI
  if (mpi_active != FALSE && nproc > 1) {
    if (n > (size_t)INT_MAX) return -1;
    if (MPI_Allreduce(MPI_IN_PLACE, acc, (int)n, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) return -1;
  }
#else
  (void)mpi_active; (void)acc; (void)n;
#endif
  return 0;
}

static int compare_ulong(const void *a, const void *b)
{
  unsigned long int x = *(const unsigned long int *)a, y = *(const unsigned long int *)b;
  return x < y ? -1 : (x > y ? 1 : 0);
}

static int find_ghost(const struct SymmetryVectorHaloPlan *halo, unsigned long int global_index, size_t *position)
{
  size_t left = 0U, right = halo->ghost_count;
  while (left < right) {
    size_t middle = left + (right - left) / 2U;
    if (halo->ghost_global_index[middle] < global_index) left = middle + 1U; else right = middle;
  }
  if (left >= halo->ghost_count || halo->ghost_global_index[left] != global_index) return -1;
  *position = left;
  return 0;
}

/* ---- memory accounting ----------------------------------------------------- */

/* Bytes that live simultaneously per block row, from the arrays allocated in
 * process_wave(): row_counts, and per off-diagonal member one transition, one
 * key, one resolved index, one norm, one fetched value and one halo column. */
static int bytes_per_row(uint64_t offdiag_total, uint64_t *out)
{
  uint64_t per_member = sizeof(struct CorrelationTransition) + sizeof(unsigned long int) +
                        sizeof(unsigned long int) + sizeof(double) + sizeof(double complex) +
                        sizeof(unsigned long int);
  uint64_t products;
  if (checked_mul64(offdiag_total, per_member, &products) != 0) return -1;
  return checked_add64(products, sizeof(size_t), out);
}

static int choose_block_rows(const struct CorrelationContext *ctx, unsigned long int *rows)
{
  uint64_t b = (uint64_t)HPHI_SYMMETRY_PLAN_LOCAL_ROWS_PER_BLOCK;
  uint64_t offdiag_total = ctx->table->offdiag_total;
  if (offdiag_total > 0U) {
    uint64_t cap = (uint64_t)HPHI_SYMMETRY_CORRELATION_BLOCK_TRANSITIONS / offdiag_total;
    if (cap == 0U) return -1;      /* one row does not fit the transition cap */
    if (cap < b) b = cap;
  }
  if (ctx->policy.hard_byte_limit != 0U) {
    uint64_t per_row, fit;
    if (bytes_per_row(offdiag_total, &per_row) != 0) return -1;
    if (ctx->fixed_bytes >= ctx->policy.hard_byte_limit) return -1;
    fit = (ctx->policy.hard_byte_limit - ctx->fixed_bytes) / per_row;
    if (fit == 0U) return -1;      /* one row does not fit the hard limit */
    if (fit < b) b = fit;
  }
  if (b > (uint64_t)ULONG_MAX) b = (uint64_t)ULONG_MAX;
  *rows = (unsigned long int)b;
  return 0;
}

/* Collective: compares the observed live bytes on this stage with the policy. */
static int check_memory(struct CorrelationContext *ctx, uint64_t observed, const char *stage)
{
  int status;
  if (observed > ctx->stats.peak_bytes) ctx->stats.peak_bytes = observed;
  status = SymmetryCheckMemoryPolicy(&ctx->policy, observed, ctx->mpi_active, myrank,
                                     "correlation", stage, ctx->warning_emitted == 0);
  if (status < 0) return -1;
  if (status > 0) ctx->warning_emitted = 1;
  return 0;
}

/* ---- apply within one block ------------------------------------------------ */

static int collect_block(const struct CorrelationContext *ctx,
                         unsigned long int row_begin, unsigned long int row_count,
                         struct CorrelationTransition *transitions, size_t *row_counts,
                         size_t *transition_count, double complex *acc)
{
  const struct BindStruct *X = ctx->X;
  const struct DefineList *def = &X->Def;
  const struct CorrelationOrbitTable *table = ctx->table;
  const size_t per_row = (size_t)table->offdiag_total;
  const size_t orbit_count = table->orbit_count;
  int error = 0;
  *transition_count = 0;
#pragma omp parallel
  {
    double complex *partial = (double complex *)calloc(orbit_count > 0U ? orbit_count : 1U, sizeof(*partial));
    int thread_error = partial == NULL;
    long i;
#pragma omp for schedule(dynamic, 16)
    for (i = 0; i < (long)row_count; i++) {
      unsigned long int local_index = row_begin + (unsigned long int)i + 1UL;
      const struct SymmetryBasisVector *entry = SymmetryBasisLocalEntry(X->Sym, local_index);
      double complex a;
      struct CorrelationTransition *slot = per_row > 0U ? transitions + (size_t)i * per_row : NULL;
      size_t written = 0, k, m;
      if (thread_error != 0) continue;
      if (entry == NULL || entry->norm == 0.0) { thread_error = 1; continue; }
      a = ctx->vec[local_index];
      for (k = 0; k < orbit_count && thread_error == 0; k++) {
        const struct CorrelationOrbit *orbit = &table->orbits[k];
        size_t len = 4U * (size_t)orbit->factors;
        for (m = 0; m < orbit->member_count; m++) {
          unsigned long int out = 0UL;
          double sign = 0.0;
          int status = ApplySymmetryFactors(def, orbit->factors, orbit->members + m * len, entry->rep_state, &out, &sign);
          if (status < 0) { thread_error = 1; break; }
          if (status == 0) continue;
          if (out == entry->rep_state) { partial[k] += sign * conj(a) * a; continue; }
          if (orbit->offdiag[m] == 0U || slot == NULL || written >= per_row) { thread_error = 1; break; }
          if (ctx->replicated) {
            struct SymmetryCanonicalResult res;
            const struct SymmetryBasisVector *target;
            if (SymmetryCanonicalizeState(X, out, &res) != 0) { thread_error = 1; break; }
            if (res.found != TRUE) continue;
            target = SymmetryBasisReplicatedGlobalEntry(X->Sym, res.basis_index);
            if (target == NULL) { thread_error = 1; break; }
            slot[written].key = res.basis_index;
            slot[written].coefficient = a * sign * res.phase * (target->norm / entry->norm);
          } else {
            struct SymmetryRepresentativeResult rep;
            if (SymmetryFindRepresentative(X, out, &rep) != 0) { thread_error = 1; break; }
            slot[written].key = rep.rep_state;
            slot[written].coefficient = a * sign * rep.phase / entry->norm;
          }
          slot[written].orbit = (uint32_t)k;
          slot[written].reserved = 0U;
          written++;
        }
      }
      row_counts[i] = written;
    }
#pragma omp critical
    {
      size_t k;
      if (thread_error != 0) error = 1;
      else for (k = 0; k < orbit_count; k++) acc[k] += partial[k];
    }
    free(partial);
  }
  if (error == 0) {
    size_t total = 0;
    unsigned long int r;
    for (r = 0; r < row_count; r++) {
      if (row_counts[r] > 0U && total != (size_t)r * per_row)
        memmove(transitions + total, transitions + (size_t)r * per_row, row_counts[r] * sizeof(*transitions));
      total += row_counts[r];
    }
    *transition_count = total;
  }
  return error ? -1 : 0;
}

/* ---- one wave: collect, resolve, fetch, accumulate ------------------------- */

static int process_wave(struct CorrelationContext *ctx,
                        unsigned long int row_begin, unsigned long int row_count,
                        double complex *acc)
{
  const struct BindStruct *X = ctx->X;
  struct CorrelationTransition *transitions = NULL;
  size_t *row_counts = NULL;
  unsigned long int *keys = NULL, *beta = NULL, *columns = NULL;
  double *norm = NULL;
  double complex *value = NULL;
  struct SymmetryVectorHaloPlan halo;
  struct SymmetryGlobalColumnSpan span;
  size_t capacity = 0, count = 0, unique = 0, found = 0, u, p, local_columns = 0, remote_columns = 0;
  uint64_t observed;
  int local_error = 0, halo_ready = 0;
  memset(&halo, 0, sizeof(halo));

  row_counts = (size_t *)calloc(row_count > 0UL ? row_count : 1UL, sizeof(*row_counts));
  if (row_counts == NULL) local_error = 1;
  else if (SymmetryCheckedSizeMul((size_t)row_count, (size_t)ctx->table->offdiag_total, &capacity) != 0) local_error = 1;
  else if (capacity > 0U) {
    transitions = (struct CorrelationTransition *)malloc(capacity * sizeof(*transitions));
    if (transitions == NULL) local_error = 1;
  }
  if (local_error == 0 && row_count > 0UL &&
      collect_block(ctx, row_begin, row_count, transitions, row_counts, &count, acc) != 0) local_error = 1;
  if (SymmetryMpiAgreeError(ctx->mpi_active, local_error) != 0) goto fail;
  ctx->stats.transitions += count;

  /* Diagonal-only requests: no keys, no directory, no halo. The condition is
   * global (offdiag_total is identical on every rank), so every rank skips. */
  if (ctx->table->offdiag_total == 0U) {
    observed = ctx->fixed_bytes + (uint64_t)row_count * sizeof(*row_counts);
    if (check_memory(ctx, observed, "diagonal-wave") != 0) goto fail;
    free(row_counts);
    return 0;
  }

  keys = (unsigned long int *)malloc((count > 0U ? count : 1U) * sizeof(*keys));
  if (keys == NULL) local_error = 1;
  else {
    for (p = 0; p < count; p++) keys[p] = transitions[p].key;
    qsort(keys, count, sizeof(*keys), compare_ulong);
    for (p = 0; p < count; p++) if (p == 0U || keys[p] != keys[unique - 1U]) keys[unique++] = keys[p];
  }
  if (local_error == 0) {
    size_t n = unique > 0U ? unique : 1U;
    beta = (unsigned long int *)calloc(n, sizeof(*beta));
    norm = (double *)calloc(n, sizeof(*norm));
    value = (double complex *)calloc(n, sizeof(*value));
    columns = (unsigned long int *)calloc(n, sizeof(*columns));
    if (beta == NULL || norm == NULL || value == NULL || columns == NULL) local_error = 1;
  }
  if (SymmetryMpiAgreeError(ctx->mpi_active, local_error) != 0) goto fail;
  ctx->stats.unique_keys += unique;

  if (ctx->replicated) {
    for (u = 0; u < unique; u++) { beta[u] = keys[u]; norm[u] = 1.0; }
  } else if (SymmetryResolveRepresentativeBatchWithOptions(
                 ctx->directory, keys, (uint64_t)unique, beta, norm, NULL) != 0) {
    local_error = 1;   /* collective; called with unique == 0 as well */
  }
  if (SymmetryMpiAgreeError(ctx->mpi_active, local_error) != 0) goto fail;

  for (u = 0; u < unique; u++) if (beta[u] != 0UL) columns[found++] = beta[u];
  span.columns = columns;
  span.count = found;
  if (BuildSymmetryVectorHaloPlan(&halo, X->Sym->dim, X->Sym->local_offset, X->Sym->local_dim,
                                  &span, 1U, nproc, myrank, &local_columns, &remote_columns) != 0) goto fail;
  halo_ready = 1;
  if (ExchangeSymmetryVectorHalo(&halo, ctx->vec) != 0) local_error = 1;
  if (SymmetryMpiAgreeError(ctx->mpi_active, local_error) != 0) goto fail;

  for (u = 0; u < unique && local_error == 0; u++) {
    unsigned long int g = beta[u];
    if (g == 0UL) continue;
    if (g > X->Sym->local_offset && g <= X->Sym->local_offset + X->Sym->local_dim) {
      value[u] = ctx->vec[g - X->Sym->local_offset];
    } else {
      size_t position;
      if (find_ghost(&halo, g, &position) != 0) local_error = 1;
      else value[u] = halo.ghost_values[position];
    }
  }
  for (p = 0; p < count && local_error == 0; p++) {
    const unsigned long int *hit = (const unsigned long int *)bsearch(
        &transitions[p].key, keys, unique, sizeof(*keys), compare_ulong);
    if (hit == NULL) { local_error = 1; break; }
    u = (size_t)(hit - keys);
    if (beta[u] == 0UL) continue;
    acc[transitions[p].orbit] += transitions[p].coefficient * norm[u] * conj(value[u]);
  }
  observed = ctx->fixed_bytes +
             (uint64_t)row_count * sizeof(*row_counts) +
             (uint64_t)capacity * sizeof(*transitions) +
             (uint64_t)count * sizeof(*keys) +
             (uint64_t)unique * (sizeof(*beta) + sizeof(*norm) + sizeof(*value) + sizeof(*columns)) +
             (uint64_t)halo.topology_scratch_bytes + (uint64_t)halo.schedule_bytes + (uint64_t)halo.runtime_buffer_bytes;
  if (check_memory(ctx, observed, "block-wave") != 0) local_error = 1;     /* collective */
  if (SymmetryMpiAgreeError(ctx->mpi_active, local_error) != 0) goto fail;

  FreeSymmetryVectorHaloPlan(&halo);
  free(transitions); free(row_counts); free(keys); free(beta); free(norm); free(value); free(columns);
  return 0;
fail:
  if (halo_ready) FreeSymmetryVectorHaloPlan(&halo);
  free(transitions); free(row_counts); free(keys); free(beta); free(norm); free(value); free(columns);
  return -1;
}

/* ---- driver ------------------------------------------------------------------ */

static int expectation_impl(const struct BindStruct *X, const double complex *vec,
                            const struct SymmetryCorrelationOperator *ops, size_t count,
                            double complex *values, struct SymmetryCorrelationStats *stats_out,
                            int verify)
{
  struct CorrelationOrbitTable table;
  struct CorrelationContext ctx;
  double complex *acc = NULL;
  unsigned long long local_waves = 0ULL, waves = 0ULL, wave;
  uint64_t setup_bytes = 0U, product = 0U;
  int local_error = 0;
  size_t t;
  memset(&table, 0, sizeof(table));
  memset(&ctx, 0, sizeof(ctx));
  if (count == 0U) return 0;
  ctx.mpi_active = SymmetryMpiCollectivesActive();
  if (X == NULL || X->Sym == NULL || X->Sym->enabled != TRUE || vec == NULL || ops == NULL || values == NULL ||
      build_orbit_table(&X->Def, ops, count, &table) != 0) local_error = 1;
  if (SymmetryMpiAgreeError(ctx.mpi_active, local_error) != 0) goto fail;
  ctx.X = X; ctx.vec = vec; ctx.table = &table;
  ctx.replicated = X->Sym->basis_layout == SYMMETRY_BASIS_REPLICATED;
  ctx.threads = team_size();
  ctx.stats.threads = ctx.threads;
  if (SymmetryLoadMemoryPolicy((uint64_t)HPHI_SYMMETRY_PLAN_BLOCK_MEMORY_BYTES, ctx.mpi_active, myrank,
                               "correlation", &ctx.policy) != 0) goto fail;   /* collective */

  /* Stage 1: orbit table (plus its construction scratch, already freed) under the policy. */
  if (checked_add64(table.bytes, table.scratch_peak, &setup_bytes) != 0) local_error = 1;
  if (SymmetryMpiAgreeError(ctx.mpi_active, local_error) != 0) goto fail;
  if (check_memory(&ctx, setup_bytes, "orbit-table") != 0) goto fail;

  /* Stage 2: temporary directory for distributed layouts with off-diagonal members.
   * The solver releases X->Sym->representative_directory after its plan is built, so
   * never rely on it; rebuild from the retained rank-local basis, as RebuildSymmetryTEPlan does. */
  if (!ctx.replicated && table.offdiag_total > 0U) {
    struct SymmetryRepresentativeDirectoryInfo info;
    struct SymmetryLocalRepresentativeIndexStats index_stats;
    if (X->Sym->local_basis == NULL || X->Sym->rank_offsets == NULL ||
        BuildSymmetryRepresentativeDirectory(X->Sym->local_basis, X->Sym->dim, X->Sym->local_dim,
                                             X->Sym->local_capacity, X->Sym->local_offset,
                                             X->Sym->rank_offsets, myrank, nproc, &ctx.directory) != 0 ||
        GetSymmetryRepresentativeDirectoryInfo(ctx.directory, &info) != 0 ||
        GetSymmetryRepresentativeDirectoryLocalIndexStats(ctx.directory, &index_stats) != 0) local_error = 1;
    else ctx.directory_bytes = (uint64_t)info.splitter_bytes + (uint64_t)index_stats.table_bytes;
    if (SymmetryMpiAgreeError(ctx.mpi_active, local_error) != 0) goto fail;
  }

  /* Stage 3: accumulators and per-thread partials; fixed bytes for block sizing. */
  acc = (double complex *)calloc(table.orbit_count, sizeof(*acc));
  if (acc == NULL) local_error = 1;
  else if (checked_mul64((uint64_t)table.orbit_count * sizeof(*acc), (uint64_t)ctx.threads + 1U, &product) != 0 ||
           checked_add64(table.bytes, ctx.directory_bytes, &ctx.fixed_bytes) != 0 ||
           checked_add64(ctx.fixed_bytes, product, &ctx.fixed_bytes) != 0) local_error = 1;
  if (SymmetryMpiAgreeError(ctx.mpi_active, local_error) != 0) goto fail;
  /* Directory tables and OpenMP team sizes may differ by rank. Size every
   * block against the largest fixed allocation so all ranks choose one safe
   * row count under a finite hard limit. */
  if (max_fixed_bytes(ctx.mpi_active, ctx.fixed_bytes, &ctx.fixed_bytes) != 0) goto fail;
  if (choose_block_rows(&ctx, &ctx.block_rows) != 0) local_error = 1;
  if (SymmetryMpiAgreeError(ctx.mpi_active, local_error) != 0) goto fail;
  if (agree_block_rows(ctx.mpi_active, ctx.block_rows) != 0) goto fail;
  if (check_memory(&ctx, ctx.fixed_bytes, "setup") != 0) goto fail;

  /* Stage 4: bounded waves, the same number on every rank. Diagonal-only
   * requests walk the same waves without directory or halo traffic. */
  local_waves = (X->Sym->local_dim + ctx.block_rows - 1UL) / ctx.block_rows;
  if (max_wave_count(ctx.mpi_active, local_waves, &waves) != 0) goto fail;
  ctx.stats.waves = waves;
  for (wave = 0ULL; wave < waves; wave++) {
    unsigned long int row_begin = 0UL, row_count = 0UL;
    if (wave < local_waves) {
      row_begin = (unsigned long int)wave * ctx.block_rows;
      row_count = X->Sym->local_dim - row_begin;
      if (row_count > ctx.block_rows) row_count = ctx.block_rows;
    }
    if (process_wave(&ctx, row_begin, row_count, acc) != 0) goto fail;
  }
  if (allreduce_sum(ctx.mpi_active, acc, table.orbit_count) != 0) goto fail;
  for (t = 0; t < count; t++) {
    const struct CorrelationOrbit *orbit = &table.orbits[table.op_to_orbit[t]];
    values[t] = acc[table.op_to_orbit[t]] / (double)orbit->member_count;
  }
  if (verify) {
    /* Verify the orbit grouping only: each row evaluated alone (without
     * verification, so this cannot recurse) must give the same value. */
    for (t = 0; t < count; t++) {
      double complex single = 0.0;
      if (expectation_impl(X, vec, &ops[t], 1U, &single, NULL, 0) != 0 ||
          cabs(single - values[t]) > 1.0e-10 * (1.0 + cabs(single))) {
        fprintf(stdoutMPI, "Error: symmetry correlation orbit grouping mismatch at row %zu.\n", t);
        goto fail;
      }
    }
  }
  if (stats_out != NULL) *stats_out = ctx.stats;
  fprintf(stdoutMPI,
          "Symmetry correlation: rows=%zu orbits=%zu waves=%llu transitions=%llu "
          "peak_bytes=%llu threads=%u layout=%s\n",
          count, table.orbit_count, ctx.stats.waves, ctx.stats.transitions,
          ctx.stats.peak_bytes, ctx.stats.threads, ctx.replicated ? "replicated" : "distributed");
  FreeSymmetryRepresentativeDirectory(ctx.directory);
  free(acc);
  free_orbit_table(&table);
  return 0;
fail:
  FreeSymmetryRepresentativeDirectory(ctx.directory);
  free(acc);
  free_orbit_table(&table);
  return -1;
}

int SymmetryCorrelationExpectationWithStats(const struct BindStruct *X, const double complex *vec,
                                            const struct SymmetryCorrelationOperator *ops, size_t count,
                                            double complex *values, struct SymmetryCorrelationStats *stats)
{
#ifdef HPHI_SYMMETRY_CORRELATION_VERIFY
  return expectation_impl(X, vec, ops, count, values, stats, 1);
#else
  return expectation_impl(X, vec, ops, count, values, stats, 0);
#endif
}

int SymmetryCorrelationExpectation(const struct BindStruct *X, const double complex *vec,
                                   const struct SymmetryCorrelationOperator *ops, size_t count,
                                   double complex *values)
{
  return SymmetryCorrelationExpectationWithStats(X, vec, ops, count, values, NULL);
}
