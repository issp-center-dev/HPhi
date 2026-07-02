#include <math.h>
#include <bitcalc.h>
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
    double norm = cabs(def->SymTransChar[g]);
    if (fabs(norm - 1.0) > eps_ch) {
      fprintf(stdoutMPI, "Error: TransSym character must have unit norm; op=%u abs=% .16e.\n",
              g, norm);
      return -1;
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
  next = (struct SymmetryBasisVector *)realloc(sym->basis,
      sizeof(struct SymmetryBasisVector) * (needed + 1));
  if (next == NULL) return -1;
  sym->basis = next;
  return 0;
}

static int raw_index_from_state(const struct BindStruct *X,
                                unsigned long int state,
                                unsigned long int *raw_index)
{
  return GetOffComp(list_2_1, list_2_2, state,
                    X->Large.irght, X->Large.ilft, X->Large.ihfbit,
                    raw_index);
}

static int store_basis_vector(struct SymmetryBasisRuntime *sym,
                              unsigned long int basis_id,
                              unsigned long int rep_state,
                              double complex *raw_coeff,
                              unsigned long int full_dim,
                              double norm)
{
  unsigned long int raw;
  unsigned int count = 0;
  unsigned int pos = 0;
  sym->basis[basis_id].rep_state = rep_state;
  for (raw = 1; raw <= full_dim; raw++) {
    if (cabs(raw_coeff[raw]) > 1.0e-12) count++;
  }
  sym->basis[basis_id].count = count;
  sym->basis[basis_id].raw_index = (unsigned long int *)malloc(sizeof(unsigned long int) * count);
  sym->basis[basis_id].coeff = (double complex *)malloc(sizeof(double complex) * count);
  if (sym->basis[basis_id].raw_index == NULL || sym->basis[basis_id].coeff == NULL) return -1;
  for (raw = 1; raw <= full_dim; raw++) {
    if (cabs(raw_coeff[raw]) > 1.0e-12) {
      double complex c = raw_coeff[raw] / norm;
      sym->basis[basis_id].raw_index[pos] = raw;
      sym->basis[basis_id].coeff[pos] = c;
      sym->raw_to_sym[raw] = basis_id;
      sym->raw_to_coeff[raw] = c;
      pos++;
    }
  }
  return 0;
}

int BuildSymmetryBasis(struct BindStruct *X)
{
  unsigned long int raw, full_dim;
  int *visited;
  double complex *acc;
  struct SymmetryBasisRuntime *sym;

  if (X->Def.iFlgSymmetryBasis == FALSE) return 0;
  full_dim = X->Check.idim_max;
  sym = (struct SymmetryBasisRuntime *)calloc(1, sizeof(*sym));
  visited = (int *)calloc(full_dim + 1, sizeof(int));
  acc = (double complex *)calloc(full_dim + 1, sizeof(double complex));
  if (sym == NULL || visited == NULL || acc == NULL) return -1;

  sym->enabled = TRUE;
  sym->nsite = X->Def.Nsite;
  sym->group_order = X->Def.NSymTrans;
  sym->full_dim = full_dim;
  sym->raw_to_sym = (unsigned long int *)calloc(full_dim + 1, sizeof(unsigned long int));
  sym->raw_to_coeff = (double complex *)calloc(full_dim + 1, sizeof(double complex));
  if (sym->raw_to_sym == NULL || sym->raw_to_coeff == NULL) return -1;

  for (raw = 1; raw <= full_dim; raw++) {
    unsigned int g;
    double norm2 = 0.0;
    if (visited[raw] != 0) continue;
    memset(acc, 0, sizeof(double complex) * (full_dim + 1));
    for (g = 0; g < X->Def.NSymTrans; g++) {
      unsigned long int moved_state = SymmetryApplyToSpinBits(list_1[raw], X->Def.SymTrans[g], X->Def.Nsite);
      unsigned long int moved_raw = 0;
      if (raw_index_from_state(X, moved_state, &moved_raw) != TRUE) {
        fprintf(stdoutMPI, "Error: TransSym maps a fixed-Sz state outside the fixed-Sz basis.\n");
        return -1;
      }
      acc[moved_raw] += conj(X->Def.SymTransChar[g]);
    }
    for (g = 0; g < X->Def.NSymTrans; g++) {
      unsigned long int moved_state = SymmetryApplyToSpinBits(list_1[raw], X->Def.SymTrans[g], X->Def.Nsite);
      unsigned long int moved_raw = 0;
      if (raw_index_from_state(X, moved_state, &moved_raw) != TRUE) {
        fprintf(stdoutMPI, "Error: TransSym maps a fixed-Sz state outside the fixed-Sz basis.\n");
        return -1;
      }
      if (visited[moved_raw] == 0) {
        visited[moved_raw] = 1;
      }
    }
    {
      unsigned long int idx;
      for (idx = 1; idx <= full_dim; idx++) norm2 += creal(conj(acc[idx]) * acc[idx]);
    }
    if (norm2 <= 1.0e-20) continue;
    sym->dim++;
    if (ensure_basis_capacity(sym, sym->dim) != 0) return -1;
    if (store_basis_vector(sym, sym->dim, list_1[raw], acc, full_dim, sqrt(norm2)) != 0) return -1;
  }

  sym->sym_diagonal = (double *)calloc(sym->dim + 1, sizeof(double));
  if (sym->sym_diagonal == NULL) return -1;
  for (raw = 1; raw <= full_dim; raw++) {
    unsigned long int b = sym->raw_to_sym[raw];
    if (b != 0) {
      double weight = creal(conj(sym->raw_to_coeff[raw]) * sym->raw_to_coeff[raw]);
      sym->sym_diagonal[b] += weight * list_Diagonal[raw];
    }
  }

  X->Sym = sym;
  free(visited);
  free(acc);
  fprintf(stdoutMPI, "Symmetry basis: raw_dim=%lu sector_dim=%lu group_order=%u\n",
          sym->full_dim, sym->dim, sym->group_order);
  if (sym->dim == 0) {
    fprintf(stdoutMPI, "Error: TransSym sector has zero basis dimension.\n");
    return -1;
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

void FreeSymmetryBasis(struct SymmetryBasisRuntime *sym)
{
  unsigned long int b;
  if (sym == NULL) return;
  if (sym->basis != NULL) {
    for (b = 1; b <= sym->dim; b++) {
      free(sym->basis[b].raw_index);
      free(sym->basis[b].coeff);
    }
  }
  free(sym->basis);
  free(sym->raw_to_sym);
  free(sym->raw_to_coeff);
  free(sym->sym_diagonal);
  free(sym);
}
