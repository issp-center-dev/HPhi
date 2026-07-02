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

int ValidateSymmetryGroupInput(const struct DefineList *def)
{
  if (def->iFlgSymmetryBasis == FALSE) return 0;
  if (def->NSymTrans == 0) {
    fprintf(stdoutMPI, "Error: TransSym requires NQPTrans > 0.\n");
    return -1;
  }
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
