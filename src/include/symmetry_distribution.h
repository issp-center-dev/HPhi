#ifndef HPHI_SYMMETRY_DISTRIBUTION_H
#define HPHI_SYMMETRY_DISTRIBUTION_H

#include <stddef.h>
#include <stdint.h>
#include "symmetry_memory_policy.h"

#ifndef HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES
#ifdef HPHI_SYMMETRY_SAMPLE_SORT_MEMORY_BYTES
#define HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES \
  HPHI_SYMMETRY_SAMPLE_SORT_MEMORY_BYTES
#else
#define HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES \
  HPHI_SYMMETRY_MEMORY_LIMIT_BYTES
#endif
#endif

/*
 * Compatibility alias for callers that configured the C3-only name before
 * exact rebalance was placed under the same distribution memory contract.
 */
#ifndef HPHI_SYMMETRY_SAMPLE_SORT_MEMORY_BYTES
#define HPHI_SYMMETRY_SAMPLE_SORT_MEMORY_BYTES \
  HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES
#endif

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
  uint64_t nonempty_rank_count;
  uint64_t samples_per_nonempty_rank;
  uint64_t local_sample_entries;
  uint64_t global_sample_entries;
  int sample_gather_used_chunked;
  uint64_t sample_gather_message_byte_limit;
  uint64_t sample_gather_max_message_bytes;
  /* C3 range redistribution after splitter selection. */
  uint64_t range_send_entries;
  uint64_t range_recv_entries;
  int range_exchange_used_chunked;
  uint64_t range_exchange_message_byte_limit;
  uint64_t range_exchange_max_message_bytes;
  uint64_t local_sample_gap;
  uint64_t global_sample_gap_max;
  uint64_t global_sample_gap_sum;
  uint64_t bucket_sample_entries;
  uint64_t bucket_entry_upper_bound;
  uint64_t global_entries;
  uint64_t range_entries;
  uint64_t range_max_entries;
  /* C4 exact block rebalance; left zero by C3. */
  uint64_t rebalance_send_entries;
  uint64_t rebalance_recv_entries;
  int rebalance_exchange_used_chunked;
  uint64_t rebalance_exchange_message_byte_limit;
  uint64_t rebalance_exchange_max_message_bytes;
  uint64_t distribution_memory_warning_byte_threshold;
  uint64_t distribution_memory_byte_limit;
  uint64_t splitter_digest;
  uint64_t range_digest;
  double range_max_over_mean;
  /* Complete C3 temporary peak, including sampling and range merge. */
  size_t sort_temporary_peak_bytes;
  /* C4 peak including the live input run; left zero by C3. */
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

/* Shared qsort comparator for replicated and distributed basis ordering. */
int SymmetryCompareBasisRepState(const void *lhs, const void *rhs);

/*
 * Sorts a valid run locally and redistributes it into deterministic
 * rep_state ranges.  On success the old storage is freed and run owns the
 * new sentinel-based range.  On failure run retains its storage and count,
 * but its entries may have been locally sorted. Collective invariant
 * failures emit one diagnostic from rank zero; memory-limit diagnostics
 * include the maximum observed bytes and the configured limit.
 */
int SymmetrySampleSortBasisRun(
    struct SymmetryBasisRun *run,
    int rank,
    int nrank,
    struct SymmetryBasisDistributionStats *stats);

int SymmetryBlockRange(
    unsigned long int dim,
    int rank,
    int nrank,
    unsigned long int *offset,
    unsigned long int *count);

/*
 * Redistributes a globally range-sorted run into the exact block layout.
 * The caller must pass an empty ownership object. On success run owns the
 * sentinel-based exact block and ownership owns rank_offsets[0..nrank].
 * The configured distribution memory limit is enforced collectively before
 * the data exchange. On failure run, ownership, and stats retain their input
 * ownership/values. On success C4 updates only the rebalance_* fields; other
 * stats fields retain their C3-provided or caller-provided values. Collective
 * invariant failures use the same rank-zero diagnostic contract as C3.
 */
int SymmetryExactRebalanceBasisRun(
    struct SymmetryBasisRun *run,
    int rank,
    int nrank,
    struct SymmetryBasisOwnership *ownership,
    struct SymmetryBasisDistributionStats *stats);

void FreeSymmetryBasisRun(struct SymmetryBasisRun *run);
void FreeSymmetryBasisOwnership(
    struct SymmetryBasisOwnership *ownership);

#endif /* HPHI_SYMMETRY_DISTRIBUTION_H */
