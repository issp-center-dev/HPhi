#include <math.h>
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

int BuildSymmetryBasis(struct BindStruct *X)
{
  (void)X;
  fprintf(stdoutMPI, "Error: internal symmetry basis builder is not wired yet.\n");
  return -1;
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
  if (sym == NULL) return;
  free(sym);
}
