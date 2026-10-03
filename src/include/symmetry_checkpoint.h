#ifndef HPHI_SYMMETRY_CHECKPOINT_H
#define HPHI_SYMMETRY_CHECKPOINT_H

#include <complex.h>
#include <stdint.h>
struct BindStruct;

struct SymmetryCheckpointInfo {
  uint64_t source_method;
  uint64_t state_index;
  uint64_t step;
  double time;
  uint64_t hamiltonian_digest;
};

/* Collective, including ranks with no vector elements. Paths are relative to
 * HPhi's output directory and include the rank suffix. Vectors are one-based.
 * A load is a same-sector initial-state import: a changed H is allowed. It is
 * not a solver restart. Layout, rank count, index order and phase must match. */
int WriteSymmetryCheckpoint(const struct BindStruct *X, const char *name,
                            const double complex *vector,
                            const struct SymmetryCheckpointInfo *info);
int ReadSymmetryCheckpoint(const struct BindStruct *X, const char *name,
                           double complex *vector,
                           struct SymmetryCheckpointInfo *info);

#endif
