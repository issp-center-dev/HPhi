#ifndef HPHI_SYMMETRY_CHECKPOINT_H
#define HPHI_SYMMETRY_CHECKPOINT_H

#include <complex.h>
#include <stdint.h>
#include <stdio.h>
struct BindStruct;

/* What the HAMILTONIAN word identifies. CG/TimeEvolution store the digest
 * of the parsed definition-file H; FullDiag stores the canonical sector
 * matrix digest. The two are never compared with each other. */
enum SymmetryHamiltonianDigestKind {
  SYMMETRY_HAMILTONIAN_PARSED = 0,
  SYMMETRY_HAMILTONIAN_MATRIX_V1 = 1
};

struct SymmetryCheckpointInfo {
  uint64_t source_method;
  uint64_t state_index;
  uint64_t step;
  double time;
  uint64_t hamiltonian_digest;
  /* Appended: existing five-field initializers keep PARSED (0). A reader
   * derives it from METHOD on success. */
  uint64_t hamiltonian_digest_kind;
};

/* Collective, including ranks with no vector elements. Paths are relative to
 * HPhi's output directory and include the rank suffix. Vectors are one-based.
 * A load is a same-sector initial-state import: a changed H is allowed. It is
 * not a solver restart.  This low-level reader requires the saved rank count,
 * local ownership, index order and phase to match; layout 0/1 may differ.
 * RestartImportVector handles rank-count redistribution for valid v1/v2 sets.
 * Non-Kondo: v1, 28 little-endian uint64 words, phase 1. Kondo only: v2,
 * 30 words, phase 2; words 28/29 hold local-site mask/fixed-quantity flags.
 * Both append local_dim complex binary64 pairs, without vector index zero. */
int WriteSymmetryCheckpoint(const struct BindStruct *X, const char *name,
                            const double complex *vector,
                            const struct SymmetryCheckpointInfo *info);
int ReadSymmetryCheckpoint(const struct BindStruct *X, const char *name,
                           double complex *vector,
                           struct SymmetryCheckpointInfo *info);

/* Collective identity for words 0..29; preserves v1/v2 encoding. */
int SymmetryCheckpointIdentityWords(const struct BindStruct *X, uint64_t *words);

/* Reproduce one writer-rank contribution to the partition-dependent DATA
 * checksum used by sector v1/v2 files.  The stream must be positioned at the
 * first payload word. */
int SymmetryCheckpointLegacyFileHash(FILE *fp, uint64_t rank, uint64_t offset,
                                     uint64_t local_dim, uint64_t *hash,
                                     double *norm2);

#endif
