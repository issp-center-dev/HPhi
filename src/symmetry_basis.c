#include <limits.h>
#include <math.h>
#include "DefCommon.h"
#include "symmetry_basis.h"
#include "struct.h"
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
  if (needed <= sym->capacity) return 0;
  next_capacity = (sym->capacity == 0UL) ? 16UL : sym->capacity;
  while (next_capacity < needed) {
    if (next_capacity > ULONG_MAX / 2UL) return -1;
    next_capacity *= 2UL;
  }
  next = (struct SymmetryBasisVector *)realloc(sym->basis,
      sizeof(struct SymmetryBasisVector) * (next_capacity + 1UL));
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

static unsigned long int find_representative_spin_state(const struct DefineList *def,
                                                        unsigned long int state)
{
  unsigned int g;
  unsigned long int rep = state;
  for (g = 0; g < def->NSymTrans; g++) {
    unsigned long int moved = SymmetryApplyToSpinBits(state, def->SymTrans[g], def->Nsite);
    if (moved < rep) rep = moved;
  }
  return rep;
}

static void compute_orbit_metadata(const struct DefineList *def,
                                   unsigned long int rep_state,
                                   unsigned int *orbit_size,
                                   unsigned int *stabilizer_size,
                                   double complex *stabilizer_sum)
{
  unsigned int g;
  *stabilizer_size = 0;
  *stabilizer_sum = 0.0;
  for (g = 0; g < def->NSymTrans; g++) {
    unsigned long int moved = SymmetryApplyToSpinBits(rep_state, def->SymTrans[g], def->Nsite);
    if (moved == rep_state) {
      (*stabilizer_size)++;
      *stabilizer_sum += conj(def->SymTransChar[g]);
    }
  }
  *orbit_size = def->NSymTrans / *stabilizer_size;
}

static unsigned long int rep_state_hash(unsigned long int state);
static unsigned long int next_power_of_two(unsigned long int value);
static int insert_rep_hash(struct SymmetryBasisRuntime *sym,
                           unsigned long int rep_state,
                           unsigned long int basis_index);
static int build_rep_hash(struct SymmetryBasisRuntime *sym);

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
  state ^= state >> 16;
  state *= 0x7feb352dUL;
  state ^= state >> 15;
  state *= 0x846ca68bUL;
  state ^= state >> 16;
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

  sym->rep_hash_keys = (unsigned long int *)calloc(target_size, sizeof(unsigned long int));
  sym->rep_hash_values = (unsigned long int *)calloc(target_size, sizeof(unsigned long int));
  if (sym->rep_hash_keys == NULL || sym->rep_hash_values == NULL) return -1;
  sym->rep_hash_size = target_size;

  for (i = 1; i <= sym->dim; i++) {
    if (insert_rep_hash(sym, sym->basis[i].rep_state, i) != 0) return -1;
  }
  return 0;
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
  struct SymmetryBasisRuntime *sym;

  if (X->Def.iFlgSymmetryBasis == FALSE) return 0;
  full_dim = X->Check.idim_max;
  sym = (struct SymmetryBasisRuntime *)calloc(1, sizeof(*sym));
  if (sym == NULL) return -1;

  sym->enabled = TRUE;
  sym->nsite = X->Def.Nsite;
  sym->group_order = X->Def.NSymTrans;
  sym->full_dim = full_dim;

  for (raw = 1; raw <= full_dim; raw++) {
    unsigned long int state = list_1[raw];
    unsigned long int rep_state = find_representative_spin_state(&X->Def, state);
    unsigned int orbit_size = 0;
    unsigned int stabilizer_size = 0;
    double complex stabilizer_sum = 0.0;
    double diagonal = (list_Diagonal != NULL) ? list_Diagonal[raw] : 0.0;
    if (state != rep_state) continue;
    compute_orbit_metadata(&X->Def, rep_state, &orbit_size, &stabilizer_size, &stabilizer_sum);
    if (cabs(stabilizer_sum) < 0.5) continue;
    sym->dim++;
    if (ensure_basis_capacity(sym, sym->dim) != 0) goto fail;
    if (store_basis_vector(sym, sym->dim, rep_state, orbit_size, stabilizer_size,
                           stabilizer_sum, diagonal) != 0) goto fail;
  }

  if (sym->dim > 1) {
    qsort(sym->basis + 1, sym->dim, sizeof(struct SymmetryBasisVector),
          compare_basis_rep_state);
  }
  if (build_rep_hash(sym) != 0) goto fail;

  sym->sym_diagonal = (double *)calloc(sym->dim + 1, sizeof(double));
  if (sym->sym_diagonal == NULL) goto fail;
  for (raw = 1; raw <= sym->dim; raw++) {
    sym->sym_diagonal[raw] = sym->basis[raw].diagonal;
  }

  fprintf(stdoutMPI, "Symmetry basis: raw_dim=%lu sector_dim=%lu group_order=%u\n",
          sym->full_dim, sym->dim, sym->group_order);
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

int SymmetryCanonicalizeSpinState(const struct BindStruct *X,
                                  unsigned long int state,
                                  struct SymmetryCanonicalResult *result)
{
  unsigned int g;
  unsigned long int rep_state;
  unsigned long int basis_index;
  if (result == NULL) return -1;
  memset(result, 0, sizeof(*result));
  if (X == NULL || X->Sym == NULL || X->Sym->enabled != TRUE) return -1;

  rep_state = find_representative_spin_state(&X->Def, state);
  basis_index = find_basis_index_by_rep(X->Sym, rep_state);
  if (basis_index == 0) return 0;

  for (g = 0; g < X->Def.NSymTrans; g++) {
    unsigned long int moved = SymmetryApplyToSpinBits(rep_state, X->Def.SymTrans[g], X->Def.Nsite);
    if (moved == state) {
      result->found = TRUE;
      result->basis_index = basis_index;
      result->op_rep_to_state = g;
      result->phase = X->Def.SymTransChar[g];
      return 0;
    }
  }
  return 0;
}

void ActivateSymmetryBasisDimension(struct BindStruct *X)
{
  if (X->Sym != NULL && X->Sym->enabled == TRUE) {
    X->Check.idim_max = X->Sym->dim;
    X->Check.idim_maxMPI = X->Sym->dim;
  }
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
  free(sym->basis);
  free(sym->sym_diagonal);
  free(sym->rep_hash_keys);
  free(sym->rep_hash_values);
  free(sym);
}
