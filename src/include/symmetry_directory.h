#ifndef HPHI_SYMMETRY_DIRECTORY_H
#define HPHI_SYMMETRY_DIRECTORY_H

#include <stddef.h>
#include <stdint.h>

struct SymmetryBasisVector;
struct SymmetryLocalRepresentativeIndex;

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

#endif /* HPHI_SYMMETRY_DIRECTORY_H */
