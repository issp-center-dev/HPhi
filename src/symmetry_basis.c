#include <limits.h>
#include <math.h>
#include <stdint.h>
#include "DefCommon.h"
#include "global.h"
#include "symmetry_basis.h"
#include "symmetry_matvec_plan.h"
#include "struct.h"
#include "CalcTime.h"
#include "wrapperMPI.h"

#define SYMMETRY_BASIS_DISTRIBUTION_CHUNK 1024UL

#ifdef MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

unsigned long int SymmetryApplyToSpinBits(unsigned long int state,
                                          const int *perm,
                                          unsigned int nsite)
{
  unsigned int site;
  unsigned long int out = 0;
  for (site = 0; site < nsite; site++) {
    if ((state & (1UL << site)) != 0UL) {
      out |= (1UL << (unsigned int)perm[site]);
    }
  }
  return out;
}

static int parity_sign_from_mapped_orbitals(const unsigned int *mapped,
                                            unsigned int count)
{
  unsigned int i, j;
  unsigned int inversions = 0;
  for (i = 0; i < count; i++) {
    for (j = i + 1U; j < count; j++) {
      if (mapped[i] > mapped[j]) inversions++;
    }
  }
  return (inversions % 2U == 0U) ? 1 : -1;
}

static int apply_fermion_site_permutation(unsigned long int state,
                                          const int *perm,
                                          unsigned int nsite,
                                          unsigned int orbitals_per_site,
                                          struct SymmetryTransformResult *result)
{
  const unsigned int max_bits = (unsigned int)(sizeof(unsigned long int) * CHAR_BIT);
  unsigned int norb = nsite * orbitals_per_site;
  unsigned int orb;
  unsigned int count = 0;
  unsigned int mapped_orbitals[sizeof(unsigned long int) * CHAR_BIT];
  unsigned long int out = 0UL;
  if (orbitals_per_site == 0U || nsite > max_bits / orbitals_per_site) return -1;
  if (norb > max_bits) return -1;
  for (orb = 0; orb < norb; orb++) {
    if ((state & (1UL << orb)) != 0UL) {
      unsigned int site = orb / orbitals_per_site;
      unsigned int spin = orb % orbitals_per_site;
      unsigned int target_site = (unsigned int)perm[site];
      unsigned int mapped = orbitals_per_site * target_site + spin;
      if (mapped >= max_bits) return -1;
      mapped_orbitals[count++] = mapped;
      out |= (1UL << mapped);
    }
  }
  result->state = out;
  result->amplitude = (double)parity_sign_from_mapped_orbitals(mapped_orbitals, count);
  return 0;
}

int SymmetryApplyToState(const struct DefineList *def,
                         unsigned long int state,
                         unsigned int op,
                         struct SymmetryTransformResult *result)
{
  if (def == NULL || result == NULL || op >= def->NSymTrans) return -1;
  memset(result, 0, sizeof(*result));
  switch (def->iCalcModel) {
  case Spin:
    result->state = SymmetryApplyToSpinBits(state, def->SymTrans[op], def->Nsite);
    result->amplitude = 1.0;
    return 0;
  case SpinlessFermion:
    return apply_fermion_site_permutation(state, def->SymTrans[op], def->Nsite, 1U, result);
  case Hubbard:
    return apply_fermion_site_permutation(state, def->SymTrans[op], def->Nsite, 2U, result);
  default:
    return -1;
  }
}

static int same_perm(const int *a, const int *b, unsigned int nsite)
{
  unsigned int i;
  for (i = 0; i < nsite; i++) {
    if (a[i] != b[i]) return FALSE;
  }
  return TRUE;
}

static int find_perm(const struct DefineList *def, const int *perm)
{
  unsigned int g;
  for (g = 0; g < def->NSymTrans; g++) {
    if (same_perm(def->SymTrans[g], perm, def->Nsite) == TRUE) return (int)g;
  }
  return -1;
}

static int is_identity_perm(const int *perm, unsigned int nsite)
{
  unsigned int site;
  for (site = 0; site < nsite; site++) {
    if (perm[site] != (int)site) return FALSE;
  }
  return TRUE;
}

static int validate_bijection_and_anti(const struct DefineList *def)
{
  unsigned int g, site;
  int *seen = (int *)malloc(sizeof(int) * def->Nsite);
  if (seen == NULL) return -1;
  for (g = 0; g < def->NSymTrans; g++) {
    for (site = 0; site < def->Nsite; site++) seen[site] = 0;
    for (site = 0; site < def->Nsite; site++) {
      int target = def->SymTrans[g][site];
      if (target < 0 || (unsigned int)target >= def->Nsite) {
        free(seen);
        fprintf(stdoutMPI, "Error: TransSym op %u maps site %u outside [0, Nsite).\n", g, site);
        return -1;
      }
      if (seen[target] != 0) {
        free(seen);
        fprintf(stdoutMPI, "Error: TransSym op %u is not a bijection.\n", g);
        return -1;
      }
      if (def->SymTransAnti[g][site] != 1) {
        free(seen);
        fprintf(stdoutMPI, "Error: TransSym Anti must be 1 in v1; op=%u site=%u anti=%d.\n",
                g, site, def->SymTransAnti[g][site]);
        return -1;
      }
      seen[target] = 1;
    }
  }
  free(seen);
  return 0;
}

int ValidateSymmetryGroupInput(const struct DefineList *def)
{
  const double eps_ch = 1.0e-10;
  unsigned int g, h, site;
  int identity = -1;
  int *composed;

  if (def->iFlgSymmetryBasis == FALSE) return 0;
  if (def->NSymTrans == 0) {
    fprintf(stdoutMPI, "Error: TransSym requires NQPTrans > 0.\n");
    return -1;
  }
  if (validate_bijection_and_anti(def) != 0) return -1;

  for (g = 0; g < def->NSymTrans; g++) {
    unsigned int prev;
    double norm = cabs(def->SymTransChar[g]);
    if (!isfinite(creal(def->SymTransChar[g])) ||
        !isfinite(cimag(def->SymTransChar[g]))) {
      fprintf(stdoutMPI, "Error: TransSym character must be finite; op=%u.\n", g);
      return -1;
    }
    if (fabs(norm - 1.0) > eps_ch) {
      fprintf(stdoutMPI, "Error: TransSym character must have unit norm; op=%u abs=% .16e.\n",
              g, norm);
      return -1;
    }
    for (prev = 0; prev < g; prev++) {
      if (same_perm(def->SymTrans[g], def->SymTrans[prev], def->Nsite) == TRUE) {
        fprintf(stdoutMPI, "Error: duplicate TransSym operation: op %u and op %u are identical.\n",
                prev, g);
        return -1;
      }
    }
    if (is_identity_perm(def->SymTrans[g], def->Nsite) == TRUE) {
      if (identity >= 0) {
        fprintf(stdoutMPI, "Error: TransSym contains duplicate identity operations.\n");
        return -1;
      }
      identity = (int)g;
    }
  }
  if (identity < 0) {
    fprintf(stdoutMPI, "Error: TransSym group must contain identity permutation.\n");
    return -1;
  }
  if (cabs(def->SymTransChar[identity] - 1.0) > eps_ch) {
    fprintf(stdoutMPI, "Error: TransSym identity character must be 1.\n");
    return -1;
  }

  composed = (int *)malloc(sizeof(int) * def->Nsite);
  if (composed == NULL) return -1;
  for (g = 0; g < def->NSymTrans; g++) {
    for (h = 0; h < def->NSymTrans; h++) {
      int gh;
      for (site = 0; site < def->Nsite; site++) {
        composed[site] = def->SymTrans[g][def->SymTrans[h][site]];
      }
      gh = find_perm(def, composed);
      if (gh < 0) {
        free(composed);
        fprintf(stdoutMPI, "Error: TransSym operations are not closed; op %u * op %u is missing.\n",
                g, h);
        return -1;
      }
      if (cabs(def->SymTransChar[gh] - def->SymTransChar[g] * def->SymTransChar[h]) > eps_ch) {
        free(composed);
        fprintf(stdoutMPI, "Error: TransSym character is not multiplicative for op %u * op %u.\n",
                g, h);
        return -1;
      }
    }
  }
  free(composed);
  return 0;
}

struct SymmetryBasisCollector {
  struct SymmetryBasisVector *entries;
  unsigned long int count;
  unsigned long int capacity;
  unsigned long long raw_states;
  unsigned long long representative_candidates;
  unsigned long long compatible_survivors;
  unsigned long long transform_calls;
  int error;
};

static unsigned long int symmetry_distribution_chunk(unsigned long int full_dim,
                                                     int nrank)
{
  unsigned long int unrank = (unsigned long int)nrank;
  if (full_dim == 0UL) return 1UL;
  if (full_dim / unrank < SYMMETRY_BASIS_DISTRIBUTION_CHUNK) {
    return full_dim / unrank + (full_dim % unrank != 0UL ? 1UL : 0UL);
  }
  return SYMMETRY_BASIS_DISTRIBUTION_CHUNK;
}

static unsigned long int symmetry_rank_raw_count(unsigned long int full_dim,
                                                 unsigned long int chunk,
                                                 int rank,
                                                 int nrank)
{
  unsigned long int full_chunks = full_dim / chunk;
  unsigned long int tail = full_dim % chunk;
  unsigned long int urank = (unsigned long int)rank;
  unsigned long int unrank = (unsigned long int)nrank;
  unsigned long int count = (full_chunks / unrank) * chunk;
  if (urank < full_chunks % unrank) {
    count += chunk;
  }
  if (tail > 0UL && urank == full_chunks % unrank) count += tail;
  return count;
}

static unsigned long int symmetry_rank_raw_index(unsigned long int local_index,
                                                 unsigned long int chunk,
                                                 int rank,
                                                 int nrank)
{
  unsigned long int local_chunk = local_index / chunk;
  unsigned long int within_chunk = local_index % chunk;
  unsigned long int global_chunk =
      local_chunk * (unsigned long int)nrank + (unsigned long int)rank;
  return global_chunk * chunk + within_chunk + 1UL;
}

static int ensure_collector_capacity(struct SymmetryBasisCollector *collector,
                                     unsigned long int needed)
{
  struct SymmetryBasisVector *next;
  unsigned long int next_capacity;
  size_t element_count;
  if (needed <= collector->capacity) return 0;
  next_capacity = (collector->capacity == 0UL) ? 16UL : collector->capacity;
  while (next_capacity < needed) {
    if (next_capacity > ULONG_MAX / 2UL) return -1;
    next_capacity *= 2UL;
  }
  if (next_capacity > SIZE_MAX / sizeof(*next)) return -1;
  element_count = (size_t)next_capacity;
  next = (struct SymmetryBasisVector *)realloc(collector->entries,
      sizeof(*next) * element_count);
  if (next == NULL) return -1;
  collector->entries = next;
  collector->capacity = next_capacity;
  return 0;
}

static int compare_basis_rep_state(const void *lhs, const void *rhs)
{
  const struct SymmetryBasisVector *a = (const struct SymmetryBasisVector *)lhs;
  const struct SymmetryBasisVector *b = (const struct SymmetryBasisVector *)rhs;
  if (a->rep_state < b->rep_state) return -1;
  if (a->rep_state > b->rep_state) return 1;
  return 0;
}

static int find_representative_state(const struct DefineList *def,
                                     unsigned long int state,
                                     unsigned long int *rep,
                                     unsigned long long *transform_calls)
{
  unsigned int g;
  if (rep == NULL) return -1;
  *rep = state;
  for (g = 0; g < def->NSymTrans; g++) {
    struct SymmetryTransformResult moved;
    if (SymmetryApplyToState(def, state, g, &moved) != 0) return -1;
    if (transform_calls != NULL) (*transform_calls)++;
    if (moved.state < *rep) *rep = moved.state;
  }
  return 0;
}

static int analyze_basis_candidate(const struct DefineList *def,
                                   unsigned long int state,
                                   int *is_representative,
                                   unsigned int *orbit_size,
                                   unsigned int *stabilizer_size,
                                   double complex *stabilizer_sum,
                                   unsigned long long *transform_calls)
{
  unsigned int g;
  *is_representative = FALSE;
  *stabilizer_size = 0;
  *stabilizer_sum = 0.0;
  for (g = 0; g < def->NSymTrans; g++) {
    struct SymmetryTransformResult moved;
    if (SymmetryApplyToState(def, state, g, &moved) != 0) return -1;
    if (transform_calls != NULL) (*transform_calls)++;
    if (moved.state < state) return 0;
    if (moved.state == state) {
      (*stabilizer_size)++;
      *stabilizer_sum += conj(def->SymTransChar[g]) * moved.amplitude;
    }
  }
  if (*stabilizer_size == 0U) return -1;
  *is_representative = TRUE;
  *orbit_size = def->NSymTrans / *stabilizer_size;
  return 0;
}

static unsigned long int rep_state_hash(unsigned long int state);
static unsigned long int next_power_of_two(unsigned long int value);
static int insert_rep_hash(struct SymmetryBasisRuntime *sym,
                           unsigned long int rep_state,
                           unsigned long int basis_index);
static int build_rep_hash(struct SymmetryBasisRuntime *sym);
static void symmetry_block_range(unsigned long int dim,
                                 int rank,
                                 int nrank,
                                 unsigned long int *offset,
                                 unsigned long int *count);

static unsigned long int find_basis_index_by_rep(const struct SymmetryBasisRuntime *sym,
                                                 unsigned long int rep_state)
{
  if (sym->rep_hash_size > 0UL && sym->rep_hash_values != NULL &&
      sym->rep_hash_keys != NULL) {
    unsigned long int mask = sym->rep_hash_size - 1UL;
    unsigned long int slot = rep_state_hash(rep_state) & mask;
    unsigned long int probes;
    for (probes = 0; probes < sym->rep_hash_size; probes++) {
      unsigned long int value = sym->rep_hash_values[slot];
      if (value == 0UL) return 0UL;
      if (sym->rep_hash_keys[slot] == rep_state) return value;
      slot = (slot + 1UL) & mask;
    }
    return 0UL;
  }

  unsigned long int lo = 1;
  unsigned long int hi = sym->dim;
  while (lo <= hi) {
    unsigned long int mid = lo + (hi - lo) / 2;
    if (sym->basis[mid].rep_state == rep_state) return mid;
    if (sym->basis[mid].rep_state < rep_state) {
      lo = mid + 1;
    } else {
      hi = mid - 1;
    }
  }
  return 0;
}

static unsigned long int rep_state_hash(unsigned long int state)
{
#if ULONG_MAX > 0xffffffffUL
  state ^= state >> 30;
  state *= 0xbf58476d1ce4e5b9UL;
  state ^= state >> 27;
  state *= 0x94d049bb133111ebUL;
  state ^= state >> 31;
#else
  state ^= state >> 16;
  state *= 0x7feb352dUL;
  state ^= state >> 15;
  state *= 0x846ca68bUL;
  state ^= state >> 16;
#endif
  return state;
}

static unsigned long int next_power_of_two(unsigned long int value)
{
  unsigned long int size = 1UL;
  while (size < value) {
    if (size > ULONG_MAX / 2UL) return 0UL;
    size *= 2UL;
  }
  return size;
}

static int insert_rep_hash(struct SymmetryBasisRuntime *sym,
                           unsigned long int rep_state,
                           unsigned long int basis_index)
{
  unsigned long int mask = sym->rep_hash_size - 1UL;
  unsigned long int slot = rep_state_hash(rep_state) & mask;
  unsigned long int probes;
  for (probes = 0; probes < sym->rep_hash_size; probes++) {
    if (sym->rep_hash_values[slot] == 0UL) {
      sym->rep_hash_keys[slot] = rep_state;
      sym->rep_hash_values[slot] = basis_index;
      return 0;
    }
    if (sym->rep_hash_keys[slot] == rep_state) return -1;
    slot = (slot + 1UL) & mask;
  }
  return -1;
}

static int build_rep_hash(struct SymmetryBasisRuntime *sym)
{
  unsigned long int i;
  unsigned long int target_size;
  if (sym->dim == 0UL) return 0;
  if (sym->dim > (ULONG_MAX - 1UL) / 2UL) return -1;
  target_size = next_power_of_two(sym->dim * 2UL + 1UL);
  if (target_size == 0UL) return -1;
  if (target_size < 4UL) target_size = 4UL;
  if (target_size > SIZE_MAX / sizeof(*sym->rep_hash_keys)) return -1;

  sym->rep_hash_keys = (unsigned long int *)calloc((size_t)target_size,
                                                   sizeof(*sym->rep_hash_keys));
  sym->rep_hash_values = (unsigned long int *)calloc((size_t)target_size,
                                                     sizeof(*sym->rep_hash_values));
  if (sym->rep_hash_keys == NULL || sym->rep_hash_values == NULL) return -1;
  sym->rep_hash_size = target_size;

  for (i = 1; i <= sym->dim; i++) {
    if (insert_rep_hash(sym, sym->basis[i].rep_state, i) != 0) return -1;
  }
  return 0;
}

static void symmetry_block_range(unsigned long int dim,
                                 int rank,
                                 int nrank,
                                 unsigned long int *offset,
                                 unsigned long int *count)
{
  unsigned long int base = dim / (unsigned long int)nrank;
  unsigned long int rem = dim % (unsigned long int)nrank;
  unsigned long int urank = (unsigned long int)rank;
  *count = base + (urank < rem ? 1UL : 0UL);
  *offset = base * urank + (urank < rem ? urank : rem);
}

static void initialize_basis_vector(struct SymmetryBasisVector *entry,
                                    unsigned long int rep_state,
                                    unsigned int orbit_size,
                                    unsigned int stabilizer_size,
                                    double complex stabilizer_sum,
                                    double diagonal)
{
  entry->rep_state = rep_state;
  entry->orbit_size = orbit_size;
  entry->stabilizer_size = stabilizer_size;
  entry->stabilizer_character_sum = stabilizer_sum;
  entry->norm = sqrt((double)orbit_size *
                     creal(conj(stabilizer_sum) * stabilizer_sum));
  entry->diagonal = diagonal;
}

static int append_basis_vector(struct SymmetryBasisCollector *collector,
                               unsigned long int rep_state,
                               unsigned int orbit_size,
                               unsigned int stabilizer_size,
                               double complex stabilizer_sum,
                               double diagonal)
{
  if (collector->count == ULONG_MAX ||
      ensure_collector_capacity(collector, collector->count + 1UL) != 0) {
    return -1;
  }
  initialize_basis_vector(&collector->entries[collector->count], rep_state,
                          orbit_size, stabilizer_size, stabilizer_sum, diagonal);
  collector->count++;
  return 0;
}

static void free_basis_collectors(struct SymmetryBasisCollector *collectors,
                                  int collector_count)
{
  int thread;
  if (collectors == NULL) return;
  for (thread = 0; thread < collector_count; thread++) {
    free(collectors[thread].entries);
  }
  free(collectors);
}

#ifdef MPI
static int create_symmetry_basis_vector_type(MPI_Datatype *vector_type)
{
  struct SymmetryBasisVector sample;
  int block_lengths[6] = {1, 1, 1, 1, 1, 1};
  MPI_Aint base;
  MPI_Aint displacements[6];
  MPI_Datatype member_types[6] = {
    MPI_UNSIGNED_LONG, MPI_UNSIGNED, MPI_UNSIGNED,
    MPI_DOUBLE, MPI_DOUBLE_COMPLEX, MPI_DOUBLE
  };
  MPI_Datatype packed_type;
  int ierr;

  MPI_Get_address(&sample, &base);
  MPI_Get_address(&sample.rep_state, &displacements[0]);
  MPI_Get_address(&sample.orbit_size, &displacements[1]);
  MPI_Get_address(&sample.stabilizer_size, &displacements[2]);
  MPI_Get_address(&sample.norm, &displacements[3]);
  MPI_Get_address(&sample.stabilizer_character_sum, &displacements[4]);
  MPI_Get_address(&sample.diagonal, &displacements[5]);
  for (ierr = 0; ierr < 6; ierr++) displacements[ierr] -= base;

  ierr = MPI_Type_create_struct(6, block_lengths, displacements, member_types,
                                &packed_type);
  if (ierr != MPI_SUCCESS) return -1;
  ierr = MPI_Type_create_resized(packed_type, 0,
                                 (MPI_Aint)sizeof(struct SymmetryBasisVector),
                                 vector_type);
  MPI_Type_free(&packed_type);
  if (ierr != MPI_SUCCESS) return -1;
  ierr = MPI_Type_commit(vector_type);
  if (ierr != MPI_SUCCESS) {
    MPI_Type_free(vector_type);
    return -1;
  }
  return 0;
}
#endif

static int gather_symmetry_basis(struct SymmetryBasisRuntime *sym)
{
  if (nproc <= 1) return 0;
#ifdef MPI
  int local_count;
  int total_count = 0;
  int local_error;
  int global_error;
  int ierr;
  int rank;
  int *counts = NULL;
  int *displacements = NULL;
  struct SymmetryBasisVector *global_basis = NULL;
  MPI_Datatype vector_type = MPI_DATATYPE_NULL;

  local_error = sym->dim > (unsigned long int)INT_MAX ? 1 : 0;
  counts = (int *)malloc((size_t)nproc * sizeof(*counts));
  displacements = (int *)malloc((size_t)nproc * sizeof(*displacements));
  if (counts == NULL || displacements == NULL) local_error = 1;
  global_error = SumMPI_i(local_error);
  if (global_error != 0) goto fail;

  local_count = (int)sym->dim;
  ierr = MPI_Allgather(&local_count, 1, MPI_INT,
                       counts, 1, MPI_INT, MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) goto fail;
  for (rank = 0; rank < nproc; rank++) {
    if (counts[rank] < 0 || counts[rank] > INT_MAX - total_count) {
      local_error = 1;
      break;
    }
    displacements[rank] = total_count;
    total_count += counts[rank];
  }
  global_error = SumMPI_i(local_error);
  if (global_error != 0) goto fail;

  if ((unsigned long int)total_count >
      SIZE_MAX / sizeof(*global_basis) - 1UL) {
    local_error = 1;
  } else {
    global_basis = (struct SymmetryBasisVector *)calloc(
        (size_t)total_count + 1U, sizeof(*global_basis));
    if (global_basis == NULL) local_error = 1;
  }
  if (create_symmetry_basis_vector_type(&vector_type) != 0) local_error = 1;
  global_error = SumMPI_i(local_error);
  if (global_error != 0) goto fail;

  ierr = MPI_Allgatherv(local_count > 0 ? sym->basis + 1 : sym->basis,
                        local_count, vector_type,
                        global_basis + 1, counts, displacements, vector_type,
                        MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) goto fail;

  MPI_Type_free(&vector_type);
  free(counts);
  free(displacements);
  free(sym->basis);
  sym->basis = global_basis;
  sym->dim = (unsigned long int)total_count;
  sym->capacity = sym->dim;
  sym->basis_gather_entries = (unsigned long long)total_count;
  sym->basis_gather_bytes =
      (unsigned long long)total_count * sizeof(*global_basis);
  return 0;

fail:
  if (vector_type != MPI_DATATYPE_NULL) MPI_Type_free(&vector_type);
  free(global_basis);
  free(counts);
  free(displacements);
  return -1;
#else
  (void)sym;
  return -1;
#endif
}

int BuildSymmetryBasis(struct BindStruct *X)
{
  unsigned long int raw, full_dim;
  unsigned long int local_raw_index;
  unsigned long int rank_raw_count;
  unsigned long int distribution_chunk;
  unsigned long int basis_offset;
  unsigned long int global_representative_candidates;
  int collector_count = 1;
  int actual_thread_count = 1;
  int local_error = 0;
  int global_error;
  int thread;
  struct SymmetryBasisCollector *collectors = NULL;
  struct SymmetryBasisRuntime *sym;

  if (X->Def.iFlgSymmetryBasis == FALSE) return 0;
  full_dim = X->Check.idim_max;
  sym = (struct SymmetryBasisRuntime *)calloc(1, sizeof(*sym));
  global_error = SumMPI_i(sym == NULL ? 1 : 0);
  if (global_error != 0) {
    free(sym);
    return -1;
  }

  sym->enabled = TRUE;
  sym->nsite = X->Def.Nsite;
  sym->group_order = X->Def.NSymTrans;
  sym->full_dim = full_dim;
  distribution_chunk = symmetry_distribution_chunk(full_dim, nproc);
  rank_raw_count = symmetry_rank_raw_count(full_dim, distribution_chunk,
                                           myrank, nproc);

#ifdef _OPENMP
  collector_count = omp_get_max_threads();
#endif
  if (collector_count < 1) goto fail;
  collectors = (struct SymmetryBasisCollector *)calloc(
      (size_t)collector_count, sizeof(*collectors));
  global_error = SumMPI_i(collectors == NULL ? 1 : 0);
  if (global_error != 0) goto fail;

  StartTimer(1110);
#ifdef _OPENMP
#pragma omp parallel shared(actual_thread_count, collectors, X, rank_raw_count, distribution_chunk)
#endif
  {
    int thread_id = 0;
    struct SymmetryBasisCollector *collector;
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#pragma omp single
    actual_thread_count = omp_get_num_threads();
#endif
    collector = &collectors[thread_id];
#ifdef _OPENMP
    /* Cyclic chunks spread representative-heavy regions while retaining
       locality in list_1 and list_Diagonal. */
#pragma omp for schedule(static, distribution_chunk)
#endif
    for (local_raw_index = 0UL; local_raw_index < rank_raw_count;
         local_raw_index++) {
      unsigned long int raw_index;
      unsigned long int state;
      int is_representative = FALSE;
      unsigned int orbit_size = 0;
      unsigned int stabilizer_size = 0;
      double complex stabilizer_sum = 0.0;
      double diagonal;
      if (collector->error != 0) continue;
      collector->raw_states++;
      raw_index = symmetry_rank_raw_index(local_raw_index, distribution_chunk,
                                          myrank, nproc);
      state = list_1[raw_index];
      diagonal = (list_Diagonal != NULL) ? list_Diagonal[raw_index] : 0.0;
      if (analyze_basis_candidate(&X->Def, state, &is_representative,
                                  &orbit_size, &stabilizer_size,
                                  &stabilizer_sum,
                                  &collector->transform_calls) != 0) {
        collector->error = 1;
        continue;
      }
      if (is_representative != TRUE) continue;
      collector->representative_candidates++;
      if (cabs(stabilizer_sum) < 0.5) continue;
      if (append_basis_vector(collector, state, orbit_size, stabilizer_size,
                              stabilizer_sum, diagonal) != 0) {
        collector->error = 1;
        continue;
      }
      collector->compatible_survivors++;
    }
  }

  for (thread = 0; thread < actual_thread_count; thread++) {
    struct SymmetryBasisCollector *collector = &collectors[thread];
    if (collector->error != 0 ||
        sym->dim > ULONG_MAX - collector->count) {
      local_error = 1;
      continue;
    }
    sym->dim += collector->count;
    sym->basis_raw_states += collector->raw_states;
    sym->basis_representative_candidates +=
        collector->representative_candidates;
    sym->basis_compatible_survivors += collector->compatible_survivors;
    sym->basis_transform_calls += collector->transform_calls;
    if (collector->raw_states > sym->basis_thread_raw_states_max)
      sym->basis_thread_raw_states_max = collector->raw_states;
    if (collector->representative_candidates >
        sym->basis_thread_representative_candidates_max)
      sym->basis_thread_representative_candidates_max =
          collector->representative_candidates;
    if (collector->compatible_survivors >
        sym->basis_thread_compatible_survivors_max)
      sym->basis_thread_compatible_survivors_max =
          collector->compatible_survivors;
    if (collector->transform_calls > sym->basis_thread_transform_calls_max)
      sym->basis_thread_transform_calls_max = collector->transform_calls;
  }
  global_error = SumMPI_i(local_error);
  if (global_error != 0) {
    StopTimer(1110);
    free_basis_collectors(collectors, collector_count);
    collectors = NULL;
    goto fail;
  }
  sym->basis_thread_count = (unsigned int)actual_thread_count;
  sym->basis_orbit_metadata_calls =
      sym->basis_representative_candidates;

  if (sym->dim > SIZE_MAX / sizeof(*sym->basis) - 1UL) {
    local_error = 1;
  } else {
    sym->basis = (struct SymmetryBasisVector *)calloc(
        (size_t)sym->dim + 1U, sizeof(*sym->basis));
    if (sym->basis == NULL) local_error = 1;
  }
  global_error = SumMPI_i(local_error);
  if (global_error != 0) {
    StopTimer(1110);
    free_basis_collectors(collectors, collector_count);
    collectors = NULL;
    goto fail;
  }
  sym->capacity = sym->dim;
  basis_offset = 1UL;
  for (thread = 0; thread < actual_thread_count; thread++) {
    struct SymmetryBasisCollector *collector = &collectors[thread];
    if (collector->count > 0UL) {
      memcpy(sym->basis + basis_offset, collector->entries,
             (size_t)collector->count * sizeof(*sym->basis));
      basis_offset += collector->count;
    }
  }
  free_basis_collectors(collectors, collector_count);
  collectors = NULL;
  if (gather_symmetry_basis(sym) != 0) {
    StopTimer(1110);
    goto fail;
  }
  global_representative_candidates = SumMPI_li(
      (unsigned long int)sym->basis_representative_candidates);
  StopTimer(1110);

  StartTimer(1111);
  if (sym->dim > 1) {
    qsort(sym->basis + 1, sym->dim, sizeof(struct SymmetryBasisVector),
          compare_basis_rep_state);
  }
  StopTimer(1111);
  StartTimer(1112);
  local_error = build_rep_hash(sym) != 0 ? 1 : 0;
  global_error = SumMPI_i(local_error);
  if (global_error != 0) {
    StopTimer(1112);
    goto fail;
  }
  StopTimer(1112);

  StartTimer(1113);
  if (sym->dim > SIZE_MAX / sizeof(*sym->sym_diagonal) - 1UL) {
    StopTimer(1113);
    goto fail;
  }
  sym->sym_diagonal = (double *)calloc((size_t)sym->dim + 1U,
                                      sizeof(*sym->sym_diagonal));
  local_error = sym->sym_diagonal == NULL ? 1 : 0;
  global_error = SumMPI_i(local_error);
  if (global_error != 0) {
    StopTimer(1113);
    goto fail;
  }
  for (raw = 1; raw <= sym->dim; raw++) {
    sym->sym_diagonal[raw] = sym->basis[raw].diagonal;
  }
  StopTimer(1113);

  fprintf(stdoutMPI, "Symmetry basis: raw_dim=%lu sector_dim=%lu group_order=%u\n",
          sym->full_dim, sym->dim, sym->group_order);
  fprintf(stdoutMPI,
          "Symmetry basis build: raw_states=%lu representative_candidates=%lu "
          "compatible_survivors=%lu\n",
          sym->full_dim,
          global_representative_candidates, sym->dim);
  if (sym->dim == 0) {
    fprintf(stdoutMPI, "Error: TransSym sector has zero basis dimension.\n");
    FreeSymmetryBasis(sym);
    return -1;
  }
  X->Sym = sym;
  return 0;

fail:
  free_basis_collectors(collectors, collector_count);
  FreeSymmetryBasis(sym);
  return -1;
}

int SymmetryCanonicalizeState(const struct BindStruct *X,
                              unsigned long int state,
                              struct SymmetryCanonicalResult *result)
{
  unsigned int g;
  unsigned long int rep_state;
  unsigned long int basis_index;
  if (result == NULL) return -1;
  memset(result, 0, sizeof(*result));
  if (X == NULL || X->Sym == NULL || X->Sym->enabled != TRUE) return -1;

  if (find_representative_state(&X->Def, state, &rep_state, NULL) != 0) return -1;
  basis_index = find_basis_index_by_rep(X->Sym, rep_state);
  if (basis_index == 0) return 0;

  for (g = 0; g < X->Def.NSymTrans; g++) {
    struct SymmetryTransformResult moved;
    if (SymmetryApplyToState(&X->Def, rep_state, g, &moved) != 0) return -1;
    if (moved.state == state) {
      result->found = TRUE;
      result->basis_index = basis_index;
      result->op_rep_to_state = g;
      result->phase = X->Def.SymTransChar[g] * moved.amplitude;
      return 0;
    }
  }
  return 0;
}

int SymmetryCanonicalizeSpinState(const struct BindStruct *X,
                                  unsigned long int state,
                                  struct SymmetryCanonicalResult *result)
{
  return SymmetryCanonicalizeState(X, state, result);
}

int ActivateSymmetryBasisDimension(struct BindStruct *X)
{
  if (X->Sym != NULL && X->Sym->enabled == TRUE) {
    int rank = 0;
    if (nproc < 1 || myrank < 0 || myrank >= nproc) return -1;
    FreeSymmetryMatvecPlan(X->Sym->matvec_plan);
    X->Sym->matvec_plan = NULL;
    symmetry_block_range(X->Sym->dim, myrank, nproc,
                         &X->Sym->local_offset, &X->Sym->local_dim);
#ifdef MPI
    free(X->Sym->mpi_recvcounts);
    free(X->Sym->mpi_displs);
    free(X->Sym->mpi_full_v1);
    X->Sym->mpi_recvcounts = NULL;
    X->Sym->mpi_displs = NULL;
    X->Sym->mpi_full_v1 = NULL;
    if (nproc > 1) {
      if (X->Sym->dim > (unsigned long int)INT_MAX) {
        fprintf(stdoutMPI,
                "Error: TransSym MPI sector dimension %lu exceeds MPI int count limit.\n",
                X->Sym->dim);
        return -1;
      }
      X->Sym->mpi_recvcounts = (int *)calloc((size_t)nproc, sizeof(int));
      X->Sym->mpi_displs = (int *)calloc((size_t)nproc, sizeof(int));
      X->Sym->mpi_full_v1 = (double complex *)calloc(X->Sym->dim + 1UL,
                                                     sizeof(double complex));
      if (X->Sym->mpi_recvcounts == NULL || X->Sym->mpi_displs == NULL ||
          X->Sym->mpi_full_v1 == NULL) {
        return -1;
      }
      for (rank = 0; rank < nproc; rank++) {
        unsigned long int offset;
        unsigned long int count;
        symmetry_block_range(X->Sym->dim, rank, nproc, &offset, &count);
        X->Sym->mpi_recvcounts[rank] = (int)count;
        X->Sym->mpi_displs[rank] = (int)offset;
      }
    }
#else
    (void)rank;
#endif
    X->Check.idim_max = X->Sym->local_dim;
    X->Check.idim_maxMPI = X->Sym->dim;
  }
  return 0;
}

int SymmetryBasisGlobalToLocal(const struct SymmetryBasisRuntime *sym,
                               unsigned long int global_index,
                               unsigned long int *local_index)
{
  if (sym == NULL || global_index == 0UL) return FALSE;
  if (global_index <= sym->local_offset ||
      global_index > sym->local_offset + sym->local_dim) {
    return FALSE;
  }
  if (local_index != NULL) *local_index = global_index - sym->local_offset;
  return TRUE;
}

int ValidateSymmetrySectorOptions(const struct BindStruct *X)
{
  if (X->Def.iFlgSymmetryBasis == FALSE) return 0;
  if (X->Sym == NULL || X->Sym->enabled != TRUE) return 0;
  if (X->Def.k_exct > X->Sym->dim) {
    fprintf(stdoutMPI,
            "Error: TransSym sector dimension %lu is smaller than exct=%u.\n",
            X->Sym->dim, X->Def.k_exct);
    return -1;
  }
  return 0;
}

void FreeSymmetryBasis(struct SymmetryBasisRuntime *sym)
{
  if (sym == NULL) return;
  FreeSymmetryMatvecPlan(sym->matvec_plan);
  free(sym->basis);
  free(sym->sym_diagonal);
  free(sym->rep_hash_keys);
  free(sym->rep_hash_values);
  free(sym->mpi_recvcounts);
  free(sym->mpi_displs);
  free(sym->mpi_full_v1);
  free(sym);
}
