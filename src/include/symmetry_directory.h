#ifndef HPHI_SYMMETRY_DIRECTORY_H
#define HPHI_SYMMETRY_DIRECTORY_H

#include <stddef.h>
#include <stdint.h>

#ifndef HPHI_SYMMETRY_DIRECTORY_MEMORY_BYTES
#define HPHI_SYMMETRY_DIRECTORY_MEMORY_BYTES UINT64_C(1073741824)
#endif

struct SymmetryBasisVector;
struct SymmetryLocalRepresentativeIndex;
struct SymmetryRepresentativeDirectory;

/**
 * @brief Construction statistics for one rank-local representative index.
 *
 * A build collision is an occupied, nonmatching slot crossed while inserting
 * an entry.  A build probe includes the final empty slot.
 */
struct SymmetryLocalRepresentativeIndexStats {
  unsigned long int table_size;
  size_t table_bytes;
  uint64_t build_collisions;
  uint64_t build_probes;
  uint64_t build_max_probe;
};

/**
 * @brief Build an index over local_basis[1..local_dim].
 *
 * local_capacity is the allocated element count, including the unused
 * zero-origin element.  local_basis remains owned by the caller and must
 * outlive the index.  index_out must point to NULL and is updated only after a
 * successful build.
 */
int BuildSymmetryLocalRepresentativeIndex(
    const struct SymmetryBasisVector *local_basis,
    unsigned long int local_dim,
    unsigned long int local_capacity,
    unsigned long int local_offset,
    struct SymmetryLocalRepresentativeIndex **index_out);

/**
 * @brief Look up one representative without communication.
 *
 * A successful miss returns zero local_index/global_beta/norm.  probe_count
 * may be NULL; local_index may be NULL.  Required outputs are updated only
 * after successful validation and lookup.
 */
int SymmetryLookupLocalRepresentative(
    const struct SymmetryLocalRepresentativeIndex *index,
    unsigned long int rep_state,
    unsigned long int *local_index,
    unsigned long int *global_beta,
    double *norm,
    uint64_t *probe_count);

/**
 * @brief Copy immutable construction statistics.
 */
int GetSymmetryLocalRepresentativeIndexStats(
    const struct SymmetryLocalRepresentativeIndex *index,
    struct SymmetryLocalRepresentativeIndexStats *stats);

void FreeSymmetryLocalRepresentativeIndex(
    struct SymmetryLocalRepresentativeIndex *index);

struct SymmetryRepresentativeDirectoryInfo {
  unsigned long int dim;
  unsigned long int local_offset;
  unsigned long int local_dim;
  int rank;
  int nrank;
  int nonempty_rank_count;
  size_t splitter_bytes;
};

struct SymmetryRepresentativeBatchOptions {
  uint64_t chunk_limit;
  int force_chunked;
  /* Internal/test controls. debug_echo is boolean; -1 disables corruption. */
  int debug_echo;
  int corrupt_response_rank;
};

struct SymmetryRepresentativeBatchStats {
  uint64_t directory_batch_calls;
  uint64_t directory_request_entries_sent;
  uint64_t directory_request_entries_received;
  uint64_t directory_found_entries;
  uint64_t directory_not_found_entries;
  uint64_t directory_lookup_probe_count;
  uint64_t directory_lookup_max_probe;
  uint64_t directory_owner_peer_count_max;
  uint64_t directory_requester_peer_count_max;
  int directory_exchange_used_chunked;
  uint64_t directory_exchange_message_byte_limit;
  uint64_t directory_exchange_max_message_bytes;
  uint64_t directory_exchange_send_messages;
  uint64_t directory_exchange_recv_messages;
  size_t directory_batch_temporary_peak_bytes;
  size_t directory_batch_memory_byte_limit;
};

/**
 * @brief Collectively build a representative owner map.
 *
 * rank_offsets[0..nrank] and local_basis remain caller-owned.  The directory
 * owns a local representative index and one first-key splitter per nonempty
 * rank.  With active MPI, rank/nrank must match MPI_COMM_WORLD.  Without
 * active MPI, only rank 0 of a one-rank layout is valid.  directory_out must
 * point to NULL and is updated only after collective success.
 */
int BuildSymmetryRepresentativeDirectory(
    const struct SymmetryBasisVector *local_basis,
    unsigned long int dim,
    unsigned long int local_dim,
    unsigned long int local_capacity,
    unsigned long int local_offset,
    const unsigned long int *rank_offsets,
    int rank,
    int nrank,
    struct SymmetryRepresentativeDirectory **directory_out);

int SymmetryRepresentativeDirectoryReady(
    const struct SymmetryRepresentativeDirectory *directory);

/**
 * @brief Return the deterministic range owner for rep_state.
 *
 * A valid zero-dimensional directory returns owner -1.  Invalid calls leave
 * owner unchanged.
 */
int SymmetryRepresentativeOwner(
    const struct SymmetryRepresentativeDirectory *directory,
    unsigned long int rep_state,
    int *owner);

int SymmetryLookupDirectoryLocalRepresentative(
    const struct SymmetryRepresentativeDirectory *directory,
    unsigned long int rep_state,
    unsigned long int *local_index,
    unsigned long int *global_beta,
    double *norm,
    uint64_t *probe_count);

/**
 * @brief Resolve a sorted unique batch collectively.
 *
 * Results preserve request order. A miss is exactly {0, 0.0}. With active
 * MPI every rank must call this routine in the same collective sequence,
 * including ranks with request_count == 0. Failure leaves outputs and
 * cumulative directory statistics unchanged.
 */
int SymmetryResolveRepresentativeBatch(
    struct SymmetryRepresentativeDirectory *directory,
    const unsigned long int *request_keys,
    uint64_t request_count,
    unsigned long int *global_beta,
    double *norm);

int SymmetryResolveRepresentativeBatchWithOptions(
    struct SymmetryRepresentativeDirectory *directory,
    const unsigned long int *request_keys,
    uint64_t request_count,
    unsigned long int *global_beta,
    double *norm,
    const struct SymmetryRepresentativeBatchOptions *options);

int GetSymmetryRepresentativeDirectoryBatchStats(
    const struct SymmetryRepresentativeDirectory *directory,
    struct SymmetryRepresentativeBatchStats *stats);

int GetSymmetryRepresentativeDirectoryInfo(
    const struct SymmetryRepresentativeDirectory *directory,
    struct SymmetryRepresentativeDirectoryInfo *info);

int GetSymmetryRepresentativeDirectoryLocalIndexStats(
    const struct SymmetryRepresentativeDirectory *directory,
    struct SymmetryLocalRepresentativeIndexStats *stats);

void FreeSymmetryRepresentativeDirectory(
    struct SymmetryRepresentativeDirectory *directory);

#endif /* HPHI_SYMMETRY_DIRECTORY_H */
