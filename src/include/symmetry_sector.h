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

/* FNV-1a 64-bit step over the low `bytes` (4 or 8) little-endian bytes. */
uint64_t SymmetryHashInteger(uint64_t hash, uint64_t value, unsigned int bytes);
int ComputeSymmetryKondoSpaceDigest(const struct DefineList *def, uint64_t *digest);
/* HubbardGC/SpinlessFermionGC fixed-quantity flags: neither N (bit0) nor Sz
 * (bit1) is fixed; bit2 marks a ParticleParity and bit3 the odd one, so all
 * N, even N and odd N are 0, 4 and 12. Other models return -1. */
int GetSymmetryFermionGCFixedFlags(const struct DefineList *def, uint64_t *flags);
int ComputeSymmetryGroupDigest(const struct DefineList *def, uint64_t *digest);
/* "hphi-group-fnv1a64-v1", or the SpinFlipParity domain whose digest also
 * covers the Px flag of every operation and p. */
const char *SymmetryGroupDigestDomain(const struct DefineList *def);
int ComputeSymmetryHamiltonianDigest(const struct DefineList *def, uint64_t *digest);
/* Collective: all ranks must call, including ranks with no local entries. */
int ComputeSymmetrySectorDigest(const struct SymmetryBasisRuntime *sym,
                                struct SymmetrySectorDigest *digest);
int WriteSymmetrySectorManifest(const struct BindStruct *X);
/* Sector FullDiag outputs. The mask helper is pure: OutputHam alone, or
 * eigenvalues plus each requested correlation/eigenvector output. */
enum SymmetryFullDiagOutput {
  SYM_FD_EIGENVALUES = 1u, SYM_FD_CORRELATIONS = 2u,
  SYM_FD_EIGENVECTORS = 4u, SYM_FD_HAMILTONIAN = 8u
};
unsigned int SymmetryFullDiagOutputMask(const struct DefineList *def);
/* Collective: appends the matrix digest, the completed outputs and
 * fulldiag_status=complete to the manifest after every output succeeded. */
int WriteSymmetryFullDiagCompletion(const struct BindStruct *X,
    int from_file, uint64_t matrix_digest);
/* Append the validated cTPQ schedule before producing any thermal data. */
int WriteSymmetryCanonicalTPQSchedule(const struct BindStruct *X, int rows,
                                     const double *beta, const int *orders);

#endif
