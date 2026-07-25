#ifndef HPHI_SYMMETRY_DISTRIBUTION_H
#define HPHI_SYMMETRY_DISTRIBUTION_H

#include <stddef.h>
#include <stdint.h>

struct BindStruct;
struct SymmetryBasisRuntime;
struct SymmetryBasisVector;

/*
 * One-origin rank-local basis storage. Capacity is an element count, not a
 * byte count, and includes entries[0], which is a reserved sentinel.
 *
 * A populated run has entries != NULL and capacity >= count + 1.  A valid
 * empty run is either {NULL, 0, 0} or sentinel-only storage with capacity
 * at least one.
 */
struct SymmetryBasisRun {
  struct SymmetryBasisVector *entries;
  unsigned long int count;
  unsigned long int capacity;
};

struct SymmetryBasisOwnership {
  unsigned long int dim;
  unsigned long int local_offset;
  unsigned long int local_dim;
  unsigned long int *rank_offsets;
};

struct SymmetryBasisDistributionStats {
  uint64_t local_survivor_entries;
  uint64_t sample_send_entries;
  uint64_t sample_recv_entries;
  uint64_t range_entries;
  uint64_t rebalance_send_entries;
  uint64_t rebalance_recv_entries;
  size_t sample_temporary_peak_bytes;
  size_t rebalance_temporary_peak_bytes;
};

/*
 * The caller must pass an initially empty {NULL, 0, 0} run.  On success the
 * run owns [1..count] (including sentinel-only storage when count is zero)
 * and must be released with FreeSymmetryBasisRun().  On failure the run
 * remains {NULL, 0, 0}; no caller-owned storage is consumed.
 */
int BuildRankLocalSymmetryBasisRun(
    const struct BindStruct *X,
    struct SymmetryBasisRuntime *sym,
    struct SymmetryBasisRun *run);

void FreeSymmetryBasisRun(struct SymmetryBasisRun *run);
void FreeSymmetryBasisOwnership(
    struct SymmetryBasisOwnership *ownership);

#endif /* HPHI_SYMMETRY_DISTRIBUTION_H */
