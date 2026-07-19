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

static int ensure_basis_capacity(struct SymmetryBasisRuntime *sym,
                                 unsigned long int needed)
{
  struct SymmetryBasisVector *next;
  unsigned long int next_capacity;
  size_t element_count;
  if (needed <= sym->capacity) return 0;
  next_capacity = (sym->capacity == 0UL) ? 16UL : sym->capacity;
  while (next_capacity < needed) {
    if (next_capacity > ULONG_MAX / 2UL) return -1;
    next_capacity *= 2UL;
  }
  if (next_capacity > SIZE_MAX / sizeof(*next) - 1UL) return -1;
  element_count = (size_t)next_capacity + 1U;
  next = (struct SymmetryBasisVector *)realloc(sym->basis,
      sizeof(*next) * element_count);
  if (next == NULL) return -1;
  sym->basis = next;
  sym->capacity = next_capacity;
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
                                     unsigned long int *rep)
{
  unsigned int g;
  if (rep == NULL) return -1;
  *rep = state;
  for (g = 0; g < def->NSymTrans; g++) {
    struct SymmetryTransformResult moved;
    if (SymmetryApplyToState(def, state, g, &moved) != 0) return -1;
    if (moved.state < *rep) *rep = moved.state;
  }
  return 0;
}

static int compute_orbit_metadata(const struct DefineList *def,
                                  unsigned long int rep_state,
                                  unsigned int *orbit_size,
                                  unsigned int *stabilizer_size,
                                  double complex *stabilizer_sum)
{
  unsigned int g;
  *stabilizer_size = 0;
  *stabilizer_sum = 0.0;
  for (g = 0; g < def->NSymTrans; g++) {
    struct SymmetryTransformResult moved;
    if (SymmetryApplyToState(def, rep_state, g, &moved) != 0) return -1;
    if (moved.state == rep_state) {
      (*stabilizer_size)++;
      *stabilizer_sum += conj(def->SymTransChar[g]) * moved.amplitude;
    }
  }
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

static int store_basis_vector(struct SymmetryBasisRuntime *sym,
                              unsigned long int basis_id,
                              unsigned long int rep_state,
                              unsigned int orbit_size,
                              unsigned int stabilizer_size,
                              double complex stabilizer_sum,
                              double diagonal)
{
  sym->basis[basis_id].rep_state = rep_state;
  sym->basis[basis_id].orbit_size = orbit_size;
  sym->basis[basis_id].stabilizer_size = stabilizer_size;
  sym->basis[basis_id].stabilizer_character_sum = stabilizer_sum;
  sym->basis[basis_id].norm = sqrt((double)orbit_size * creal(conj(stabilizer_sum) * stabilizer_sum));
  sym->basis[basis_id].diagonal = diagonal;
  return 0;
}

int BuildSymmetryBasis(struct BindStruct *X)
{
  unsigned long int raw, full_dim;
  unsigned long int representative_candidates = 0UL;
  struct SymmetryBasisRuntime *sym;

  if (X->Def.iFlgSymmetryBasis == FALSE) return 0;
  full_dim = X->Check.idim_max;
  sym = (struct SymmetryBasisRuntime *)calloc(1, sizeof(*sym));
  if (sym == NULL) return -1;

  sym->enabled = TRUE;
  sym->nsite = X->Def.Nsite;
  sym->group_order = X->Def.NSymTrans;
  sym->full_dim = full_dim;

  StartTimer(1110);
  for (raw = 1; raw <= full_dim; raw++) {
    unsigned long int state = list_1[raw];
    unsigned long int rep_state = 0UL;
    unsigned int orbit_size = 0;
    unsigned int stabilizer_size = 0;
    double complex stabilizer_sum = 0.0;
    double diagonal = (list_Diagonal != NULL) ? list_Diagonal[raw] : 0.0;
    if (find_representative_state(&X->Def, state, &rep_state) != 0) {
      StopTimer(1110);
      goto fail;
    }
    if (state != rep_state) continue;
    representative_candidates++;
    if (compute_orbit_metadata(&X->Def, rep_state, &orbit_size, &stabilizer_size,
                               &stabilizer_sum) != 0) {
      StopTimer(1110);
      goto fail;
    }
    if (cabs(stabilizer_sum) < 0.5) continue;
    sym->dim++;
    if (ensure_basis_capacity(sym, sym->dim) != 0) {
      StopTimer(1110);
      goto fail;
    }
    if (store_basis_vector(sym, sym->dim, rep_state, orbit_size, stabilizer_size,
                           stabilizer_sum, diagonal) != 0) {
      StopTimer(1110);
      goto fail;
    }
  }
  StopTimer(1110);

  StartTimer(1111);
  if (sym->dim > 1) {
    qsort(sym->basis + 1, sym->dim, sizeof(struct SymmetryBasisVector),
          compare_basis_rep_state);
  }
  StopTimer(1111);
  StartTimer(1112);
  if (build_rep_hash(sym) != 0) {
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
  if (sym->sym_diagonal == NULL) {
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
          sym->full_dim, representative_candidates, sym->dim);
  if (sym->dim == 0) {
    fprintf(stdoutMPI, "Error: TransSym sector has zero basis dimension.\n");
    FreeSymmetryBasis(sym);
    return -1;
  }
  X->Sym = sym;
  return 0;

fail:
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

  if (find_representative_state(&X->Def, state, &rep_state) != 0) return -1;
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
