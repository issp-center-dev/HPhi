#ifndef HPHI_SYMMETRY_HAMILTONIAN_IO_H
#define HPHI_SYMMETRY_HAMILTONIAN_IO_H

#include <stdint.h>
struct BindStruct;

/* Basis identity of a sector matrix: the checkpoint normalization of model
 * and quantum numbers, plus the group/sector digests and the global basis
 * order. MPI rank count and dense backend are deliberately not part of it. */
struct SymmetryHamiltonianIdentity {
  uint64_t model, nsite, nup, ndown, ne, local_site_mask, fixed_flags;
  uint64_t phase, full_dim, sector_dim;
  uint64_t group_digest, sector_count, sector_xor, sector_sum, basis_order;
};

/* matrix_digest is hphi-sector-matrix-fnv1a64-v1 of the canonical matrix
 * handed to the dense solver. */
struct SymmetryHamiltonianInfo {
  struct SymmetryHamiltonianIdentity identity;
  uint64_t matrix_digest;
  int from_file; /* 0=generated, 1=input */
};

/* Both are collective and return 0 on every rank or -1 on every rank. The
 * output struct is written only after all ranks succeeded. Prepare checks
 * the stored H (replicated Ham or the ELPA column panel) for finite,
 * Hermitian coefficients, rebuilds the upper triangle from the lower one,
 * zeroes the diagonal imaginary part and hashes the result. It does not
 * change the solver dimensions, the symmetry basis or any eigenvector. */
int BuildSymmetryHamiltonianIdentity(const struct BindStruct *X,
    struct SymmetryHamiltonianIdentity *identity);
int PrepareSymmetryHamiltonian(const struct BindStruct *X, int from_file,
    struct SymmetryHamiltonianInfo *info);

/* <head>_ham_sector.mtx in the output directory: a Matrix Market
 * "coordinate complex hermitian" lower triangle with 17 significant digits,
 * preceded by one "% HPHI key=value" line for format, scalar and each
 * identity field. Both are collective and return 0 on every rank or -1 on
 * every rank.
 *
 * Input zeroes the stored H, reads and validates the file on the root,
 * sends the records in chunks to every rank (a column panel keeps only its
 * own columns, replicated storage every record), detects duplicates by
 * column owner and ends with Prepare(X, 1, info).
 * Output writes the prepared matrix through <name>.part and renames it
 * only after every write and the close succeeded. */
int InputSymmetryHamiltonian(const struct BindStruct *X,
    struct SymmetryHamiltonianInfo *info);
int OutputSymmetryHamiltonian(const struct BindStruct *X,
    const struct SymmetryHamiltonianInfo *info);

#endif
