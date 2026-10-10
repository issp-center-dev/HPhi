#ifndef HPHI_SYMMETRY_KONDO_H
#define HPHI_SYMMETRY_KONDO_H

#include <stdint.h>
struct DefineList;
struct SymmetryKondoIdentity {
  uint64_t local_site_mask, fixed_flags, nup, ndown, ne, phase;
};
int GetSymmetryKondoIdentity(const struct DefineList *def,
                             struct SymmetryKondoIdentity *identity);

int IsSymmetryKondoModel(int model);
int NormalizeSymmetryKondoQuantumNumbers(struct DefineList *def,
                                         int has_ncond, int has_sz,
                                         int has_nup, int has_ndown);
int ValidateSymmetryKondoSpace(const struct DefineList *def);
int SymmetryKondoLocalMask(const struct DefineList *def,
                           unsigned long *mask);
int SymmetryKondoStateIsPhysical(const struct DefineList *def,
                                 unsigned long state);
int SymmetryKondoPermutationSign(const struct DefineList *def,
                                 const int *permutation, int *sign);

#endif /* HPHI_SYMMETRY_KONDO_H */
