#ifndef HPHI_SYMMETRY_SECTOR_H
#define HPHI_SYMMETRY_SECTOR_H

#include <stdint.h>
struct DefineList;
struct BindStruct;
struct SymmetryBasisRuntime;

/* A multiset fingerprint, not a vector ordering or Hamiltonian identity. */
struct SymmetrySectorDigest {
  uint64_t count;
  uint64_t xor_hash;
  uint64_t sum_hash;
};

int ComputeSymmetryKondoSpaceDigest(const struct DefineList *def, uint64_t *digest);
int ComputeSymmetryGroupDigest(const struct DefineList *def, uint64_t *digest);
int ComputeSymmetryHamiltonianDigest(const struct DefineList *def, uint64_t *digest);
/* Collective: all ranks must call, including ranks with no local entries. */
int ComputeSymmetrySectorDigest(const struct SymmetryBasisRuntime *sym,
                                struct SymmetrySectorDigest *digest);
int WriteSymmetrySectorManifest(const struct BindStruct *X);
/* Append the validated cTPQ schedule before producing any thermal data. */
int WriteSymmetryCanonicalTPQSchedule(const struct BindStruct *X, int rows,
                                     const double *beta, const int *orders);

#endif
