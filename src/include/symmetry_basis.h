#ifndef HPHI_SYMMETRY_BASIS_H
#define HPHI_SYMMETRY_BASIS_H

#include "Common.h"

struct BindStruct;
struct DefineList;

struct SymmetryBasisVector {
  unsigned long int rep_state;
  unsigned int count;
  unsigned long int *raw_index;
  double complex *coeff;
};

struct SymmetryBasisRuntime {
  int enabled;
  unsigned int nsite;
  unsigned int group_order;
  unsigned long int full_dim;
  unsigned long int dim;
  unsigned long int capacity;
  struct SymmetryBasisVector *basis;
  unsigned long int *raw_to_sym;
  double complex *raw_to_coeff;
  double *sym_diagonal;
};

int ValidateSymmetryGroupInput(const struct DefineList *def);
unsigned long int SymmetryApplyToSpinBits(unsigned long int state,
                                          const int *perm,
                                          unsigned int nsite);
int BuildSymmetryBasis(struct BindStruct *X);
void ActivateSymmetryBasisDimension(struct BindStruct *X);
int ValidateSymmetrySectorOptions(const struct BindStruct *X);
void FreeSymmetryBasis(struct SymmetryBasisRuntime *sym);

#endif /* HPHI_SYMMETRY_BASIS_H */
