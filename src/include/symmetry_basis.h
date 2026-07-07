#ifndef HPHI_SYMMETRY_BASIS_H
#define HPHI_SYMMETRY_BASIS_H

#include "Common.h"

struct BindStruct;
struct DefineList;

struct SymmetryBasisVector {
  unsigned long int rep_state;
  unsigned int orbit_size;
  unsigned int stabilizer_size;
  double norm;
  double complex stabilizer_character_sum;
  double diagonal;
};

struct SymmetryCanonicalResult {
  int found;
  unsigned long int basis_index;
  unsigned int op_rep_to_state;
  double complex phase;
};

struct SymmetryBasisRuntime {
  int enabled;
  unsigned int nsite;
  unsigned int group_order;
  unsigned long int full_dim;
  unsigned long int dim;
  unsigned long int capacity;
  struct SymmetryBasisVector *basis;
  double *sym_diagonal;
  unsigned long int rep_hash_size;
  unsigned long int *rep_hash_keys;
  unsigned long int *rep_hash_values;
};

int ValidateSymmetryGroupInput(const struct DefineList *def);
unsigned long int SymmetryApplyToSpinBits(unsigned long int state,
                                          const int *perm,
                                          unsigned int nsite);
int BuildSymmetryBasis(struct BindStruct *X);
int SymmetryCanonicalizeSpinState(const struct BindStruct *X,
                                  unsigned long int state,
                                  struct SymmetryCanonicalResult *result);
void ActivateSymmetryBasisDimension(struct BindStruct *X);
int ValidateSymmetrySectorOptions(const struct BindStruct *X);
void FreeSymmetryBasis(struct SymmetryBasisRuntime *sym);

#endif /* HPHI_SYMMETRY_BASIS_H */
