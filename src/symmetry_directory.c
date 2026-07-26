#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "symmetry_basis.h"
#include "symmetry_checked.h"
#include "symmetry_directory.h"
#include "symmetry_mpi_exchange.h"
#include "CalcTime.h"

struct SymmetryLocalRepresentativeIndex {
  const struct SymmetryBasisVector *local_basis;
  unsigned long int local_dim;
  unsigned long int local_offset;
  unsigned long int table_size;
  unsigned long int *keys;
  unsigned long int *values;
  struct SymmetryLocalRepresentativeIndexStats stats;
};

struct SymmetryRepresentativeDirectory {
  int ready;
  unsigned long int dim;
  unsigned long int local_offset;
  unsigned long int local_dim;
  int rank;
  int nrank;
  int nonempty_rank_count;
  unsigned long int *rank_first_rep;
  struct SymmetryLocalRepresentativeIndex *local_index;
  struct SymmetryRepresentativeBatchStats batch_stats;
};

static unsigned long int local_rep_state_hash(unsigned long int state)
{
#if ULONG_MAX > 0xffffffffUL
  state ^= state >> 30;
  state *= 0xbf58476d1ce4e5b9UL;
  state ^= state >> 27;
  state *= 0x94d049bb133111ebUL;
  state ^= state >> 31;
#else
  state ^= state >> 16;
  state *= 0x7feb352dUL;
  state ^= state >> 15;
  state *= 0x846ca68bUL;
  state ^= state >> 16;
#endif
  return state;
}

static unsigned long int local_next_power_of_two(unsigned long int value)
{
  unsigned long int size = 1UL;
  while (size < value) {
    if (size > ULONG_MAX / 2UL) return 0UL;
    size *= 2UL;
  }
  return size;
}

static int checked_table_shape(unsigned long int local_dim,
                               unsigned long int *table_size,
                               size_t *table_bytes)
{
  unsigned long int target_size;
  size_t one_table_bytes;
  if (table_size == NULL || table_bytes == NULL) return -1;
  if (local_dim == 0UL) {
    *table_size = 0UL;
    *table_bytes = 0U;
    return 0;
  }
  if (local_dim > (ULONG_MAX - 1UL) / 2UL) return -1;
  target_size = local_next_power_of_two(local_dim * 2UL + 1UL);
  if (target_size == 0UL) return -1;
  if (target_size < 4UL) target_size = 4UL;
#if ULONG_MAX > SIZE_MAX
  if (target_size > (unsigned long int)SIZE_MAX) return -1;
#endif
  if ((size_t)target_size > SIZE_MAX / sizeof(unsigned long int)) {
    return -1;
  }
  one_table_bytes = (size_t)target_size * sizeof(unsigned long int);
  if (one_table_bytes > SIZE_MAX - one_table_bytes) return -1;
  *table_size = target_size;
  *table_bytes = one_table_bytes + one_table_bytes;
  return 0;
}

static int checked_u64_add(uint64_t lhs, uint64_t rhs, uint64_t *result)
{
  if (result == NULL || lhs > UINT64_MAX - rhs) return -1;
  *result = lhs + rhs;
  return 0;
}

static int insert_local_representative(
    struct SymmetryLocalRepresentativeIndex *index,
    unsigned long int rep_state,
    unsigned long int local_index)
{
  unsigned long int mask;
  unsigned long int slot;
  unsigned long int probes;
  uint64_t probe_count;
  uint64_t collision_count;
  if (index == NULL || index->table_size == 0UL ||
      index->keys == NULL || index->values == NULL ||
      local_index == 0UL || local_index > index->local_dim) {
    return -1;
  }
  mask = index->table_size - 1UL;
  slot = local_rep_state_hash(rep_state) & mask;
  for (probes = 1UL; probes <= index->table_size; probes++) {
    if (index->values[slot] == 0UL) {
      probe_count = (uint64_t)probes;
      collision_count = probe_count - 1U;
      if ((unsigned long int)probe_count != probes ||
          checked_u64_add(index->stats.build_probes, probe_count,
                          &index->stats.build_probes) != 0 ||
          checked_u64_add(index->stats.build_collisions, collision_count,
                          &index->stats.build_collisions) != 0) {
        return -1;
      }
      if (probe_count > index->stats.build_max_probe) {
        index->stats.build_max_probe = probe_count;
      }
      index->keys[slot] = rep_state;
      index->values[slot] = local_index;
      return 0;
    }
    if (index->keys[slot] == rep_state) return -1;
    slot = (slot + 1UL) & mask;
  }
  return -1;
}

void FreeSymmetryLocalRepresentativeIndex(
    struct SymmetryLocalRepresentativeIndex *index)
{
  if (index == NULL) return;
  free(index->keys);
  free(index->values);
  free(index);
}

int BuildSymmetryLocalRepresentativeIndex(
    const struct SymmetryBasisVector *local_basis,
    unsigned long int local_dim,
    unsigned long int local_capacity,
    unsigned long int local_offset,
    struct SymmetryLocalRepresentativeIndex **index_out)
{
  struct SymmetryLocalRepresentativeIndex *next = NULL;
  unsigned long int table_size;
  unsigned long int local_index;
  size_t table_bytes;
  if (index_out == NULL || *index_out != NULL) return -1;
  if (local_dim == 0UL) {
    if (!((local_basis == NULL && local_capacity == 0UL) ||
          (local_basis != NULL && local_capacity >= 1UL))) {
      return -1;
    }
  } else {
    if (local_basis == NULL || local_dim == ULONG_MAX ||
        local_capacity < local_dim + 1UL ||
        local_offset > ULONG_MAX - local_dim) {
      return -1;
    }
  }
  if (checked_table_shape(local_dim, &table_size, &table_bytes) != 0) {
    return -1;
  }
  for (local_index = 1UL; local_index <= local_dim; local_index++) {
    const struct SymmetryBasisVector *entry = &local_basis[local_index];
    if (!isfinite(entry->norm) || entry->norm <= 0.0 ||
        (local_index > 1UL &&
         local_basis[local_index - 1UL].rep_state >= entry->rep_state)) {
      return -1;
    }
  }

  next = (struct SymmetryLocalRepresentativeIndex *)calloc(1U,
                                                            sizeof(*next));
  if (next == NULL) return -1;
  next->local_basis = local_basis;
  next->local_dim = local_dim;
  next->local_offset = local_offset;
  next->table_size = table_size;
  next->stats.table_size = table_size;
  next->stats.table_bytes = table_bytes;
  if (table_size > 0UL) {
    next->keys = (unsigned long int *)calloc((size_t)table_size,
                                             sizeof(*next->keys));
    next->values = (unsigned long int *)calloc((size_t)table_size,
                                               sizeof(*next->values));
    if (next->keys == NULL || next->values == NULL) {
      FreeSymmetryLocalRepresentativeIndex(next);
      return -1;
    }
  }
  for (local_index = 1UL; local_index <= local_dim; local_index++) {
    if (insert_local_representative(
            next, local_basis[local_index].rep_state, local_index) != 0) {
      FreeSymmetryLocalRepresentativeIndex(next);
      return -1;
    }
  }
  *index_out = next;
  return 0;
}

int SymmetryLookupLocalRepresentative(
    const struct SymmetryLocalRepresentativeIndex *index,
    unsigned long int rep_state,
    unsigned long int *local_index,
    unsigned long int *global_beta,
    double *norm,
    uint64_t *probe_count)
{
  unsigned long int found_local_index = 0UL;
  unsigned long int found_global_beta = 0UL;
  double found_norm = 0.0;
  uint64_t probes = 0U;
  if (index == NULL || global_beta == NULL || norm == NULL ||
      (index->local_dim > 0UL &&
       (index->local_basis == NULL || index->table_size == 0UL ||
        index->keys == NULL || index->values == NULL))) {
    return -1;
  }
  if (index->table_size > 0UL) {
    unsigned long int mask = index->table_size - 1UL;
    unsigned long int slot = local_rep_state_hash(rep_state) & mask;
    unsigned long int inspected;
    for (inspected = 1UL; inspected <= index->table_size; inspected++) {
      unsigned long int value = index->values[slot];
      if ((unsigned long int)(uint64_t)inspected != inspected) {
        return -1;
      }
      probes = (uint64_t)inspected;
      if (value == 0UL) break;
      if (index->keys[slot] == rep_state) {
        const struct SymmetryBasisVector *entry;
        if (value > index->local_dim ||
            index->local_offset > ULONG_MAX - value) {
          return -1;
        }
        entry = &index->local_basis[value];
        if (!isfinite(entry->norm) || entry->norm <= 0.0 ||
            entry->rep_state != rep_state) {
          return -1;
        }
        found_local_index = value;
        found_global_beta = index->local_offset + value;
        found_norm = entry->norm;
        break;
      }
      slot = (slot + 1UL) & mask;
    }
  }
  if (local_index != NULL) *local_index = found_local_index;
  *global_beta = found_global_beta;
  *norm = found_norm;
  if (probe_count != NULL) *probe_count = probes;
  return 0;
}

int GetSymmetryLocalRepresentativeIndexStats(
    const struct SymmetryLocalRepresentativeIndex *index,
    struct SymmetryLocalRepresentativeIndexStats *stats)
{
  struct SymmetryLocalRepresentativeIndexStats next;
  if (index == NULL || stats == NULL) return -1;
  next = index->stats;
  *stats = next;
  return 0;
}

static int directory_mpi_context(int *mpi_active,
                                 int *comm_rank,
                                 int *comm_size)
{
  if (mpi_active == NULL || comm_rank == NULL || comm_size == NULL) return -1;
  *mpi_active = 0;
  *comm_rank = 0;
  *comm_size = 1;
#ifdef MPI
  {
    int initialized = 0;
    int finalized = 0;
    if (MPI_Initialized(&initialized) != MPI_SUCCESS) return -1;
    if (initialized == 0) return 0;
    if (MPI_Finalized(&finalized) != MPI_SUCCESS) return -1;
    if (finalized != 0) return 0;
    if (MPI_Comm_rank(MPI_COMM_WORLD, comm_rank) != MPI_SUCCESS ||
        MPI_Comm_size(MPI_COMM_WORLD, comm_size) != MPI_SUCCESS) {
      return -1;
    }
    *mpi_active = 1;
  }
#endif
  return 0;
}

static int agree_directory_failure(int mpi_active, int local_error)
{
  int error_flag = local_error != 0 ? 1 : 0;
#ifdef MPI
  int global_error = error_flag;
  if (mpi_active != 0 &&
      MPI_Allreduce(&error_flag, &global_error, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
  return global_error;
#else
  (void)mpi_active;
  return error_flag;
#endif
}

static int directory_block_range(unsigned long int dim,
                                 int rank,
                                 int nrank,
                                 unsigned long int *offset,
                                 unsigned long int *count)
{
  unsigned long int quotient;
  unsigned long int remainder;
  unsigned long int rank_value;
  if (rank < 0 || nrank < 1 || rank >= nrank ||
      offset == NULL || count == NULL) {
    return -1;
  }
  rank_value = (unsigned long int)rank;
  quotient = dim / (unsigned long int)nrank;
  remainder = dim % (unsigned long int)nrank;
  *offset = quotient * rank_value +
      (rank_value < remainder ? rank_value : remainder);
  *count = quotient + (rank_value < remainder ? 1UL : 0UL);
  return 0;
}

static int validate_directory_layout(unsigned long int dim,
                                     unsigned long int local_offset,
                                     unsigned long int local_dim,
                                     const unsigned long int *rank_offsets,
                                     int rank,
                                     int nrank)
{
  int peer;
  if (rank_offsets == NULL || rank < 0 || nrank < 1 || rank >= nrank ||
      rank_offsets[0] != 0UL || rank_offsets[nrank] != dim) {
    return -1;
  }
  for (peer = 0; peer < nrank; peer++) {
    unsigned long int expected_offset;
    unsigned long int expected_count;
    if (directory_block_range(dim, peer, nrank,
                              &expected_offset, &expected_count) != 0 ||
        rank_offsets[peer] != expected_offset ||
        rank_offsets[peer + 1] != expected_offset + expected_count) {
      return -1;
    }
    if (peer == rank &&
        (local_offset != expected_offset || local_dim != expected_count)) {
      return -1;
    }
  }
  return 0;
}

static int directory_nonempty_rank_count(unsigned long int dim, int nrank)
{
  if (nrank < 1) return -1;
  if (dim < (unsigned long int)nrank) return (int)dim;
  return nrank;
}

void FreeSymmetryRepresentativeDirectory(
    struct SymmetryRepresentativeDirectory *directory)
{
  if (directory == NULL) return;
  FreeSymmetryLocalRepresentativeIndex(directory->local_index);
  free(directory->rank_first_rep);
  free(directory);
}

int BuildSymmetryRepresentativeDirectory(
    const struct SymmetryBasisVector *local_basis,
    unsigned long int dim,
    unsigned long int local_dim,
    unsigned long int local_capacity,
    unsigned long int local_offset,
    const unsigned long int *rank_offsets,
    int rank,
    int nrank,
    struct SymmetryRepresentativeDirectory **directory_out)
{
  struct SymmetryLocalRepresentativeIndex *local_index = NULL;
  struct SymmetryRepresentativeDirectory *next = NULL;
  unsigned long int local_boundary[2] = {0UL, 0UL};
  unsigned long int *all_boundaries = NULL;
  size_t boundary_elements = 0U;
  size_t boundary_bytes = 0U;
  size_t splitter_bytes = 0U;
  int mpi_active;
  int comm_rank;
  int comm_size;
  int nonempty_rank_count = 0;
  int local_error = 0;
  int global_error;
  int peer;

  if (directory_mpi_context(&mpi_active, &comm_rank, &comm_size) != 0) {
    return -1;
  }
  if (directory_out == NULL || *directory_out != NULL ||
      rank < 0 || nrank < 1 || rank >= nrank ||
      (uint64_t)HPHI_SYMMETRY_DIRECTORY_MEMORY_BYTES == 0U ||
      (uint64_t)(size_t)HPHI_SYMMETRY_DIRECTORY_MEMORY_BYTES !=
          (uint64_t)HPHI_SYMMETRY_DIRECTORY_MEMORY_BYTES ||
      (mpi_active != 0 &&
       (rank != comm_rank || nrank != comm_size)) ||
      (mpi_active == 0 && (rank != 0 || nrank != 1)) ||
      validate_directory_layout(dim, local_offset, local_dim,
                                rank_offsets, rank, nrank) != 0) {
    local_error = 1;
  }
  if (local_error == 0) {
    StartTimer(1112);
    if (BuildSymmetryLocalRepresentativeIndex(
            local_basis, local_dim, local_capacity, local_offset,
            &local_index) != 0) {
      local_error = 1;
    }
    StopTimer(1112);
  }
  global_error = agree_directory_failure(mpi_active, local_error);
  if (global_error != 0) {
    FreeSymmetryLocalRepresentativeIndex(local_index);
    return -1;
  }

  nonempty_rank_count = directory_nonempty_rank_count(dim, nrank);
  if (nonempty_rank_count < 0 ||
      (size_t)nrank > SIZE_MAX / 2U) {
    local_error = 1;
  } else {
    boundary_elements = (size_t)nrank * 2U;
    if (boundary_elements > SIZE_MAX / sizeof(*all_boundaries)) {
      local_error = 1;
    } else {
      boundary_bytes = boundary_elements * sizeof(*all_boundaries);
    }
    if ((size_t)nonempty_rank_count >
        SIZE_MAX / sizeof(unsigned long int)) {
      local_error = 1;
    } else {
      splitter_bytes =
          (size_t)nonempty_rank_count * sizeof(unsigned long int);
    }
  }
  if (local_error == 0) {
    next = (struct SymmetryRepresentativeDirectory *)calloc(1U,
                                                             sizeof(*next));
    all_boundaries =
        (unsigned long int *)malloc(boundary_bytes);
    if (next == NULL || all_boundaries == NULL) local_error = 1;
  }
  if (next != NULL) {
    next->dim = dim;
    next->local_offset = local_offset;
    next->local_dim = local_dim;
    next->rank = rank;
    next->nrank = nrank;
    next->nonempty_rank_count = nonempty_rank_count;
    next->local_index = local_index;
    next->batch_stats.directory_batch_memory_byte_limit =
        (size_t)HPHI_SYMMETRY_DIRECTORY_MEMORY_BYTES;
    local_index = NULL;
    if (nonempty_rank_count > 0) {
      next->rank_first_rep =
          (unsigned long int *)malloc(splitter_bytes);
      if (next->rank_first_rep == NULL) local_error = 1;
    }
  }
  global_error = agree_directory_failure(mpi_active, local_error);
  if (global_error != 0) {
    free(all_boundaries);
    FreeSymmetryRepresentativeDirectory(next);
    FreeSymmetryLocalRepresentativeIndex(local_index);
    return -1;
  }

  if (local_dim > 0UL) {
    local_boundary[0] = local_basis[1].rep_state;
    local_boundary[1] = local_basis[local_dim].rep_state;
  }
#ifdef MPI
  if (mpi_active != 0) {
    if (MPI_Allgather(local_boundary, 2, MPI_UNSIGNED_LONG,
                      all_boundaries, 2, MPI_UNSIGNED_LONG,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      free(all_boundaries);
      FreeSymmetryRepresentativeDirectory(next);
      return -1;
    }
  } else
#endif
  {
    all_boundaries[0] = local_boundary[0];
    all_boundaries[1] = local_boundary[1];
  }

  {
    unsigned long int previous_last = 0UL;
    int have_previous = 0;
    for (peer = 0; peer < nrank; peer++) {
      unsigned long int expected_offset;
      unsigned long int expected_count;
      unsigned long int first = all_boundaries[(size_t)peer * 2U];
      unsigned long int last = all_boundaries[(size_t)peer * 2U + 1U];
      if (directory_block_range(dim, peer, nrank,
                                &expected_offset, &expected_count) != 0) {
        local_error = 1;
        break;
      }
      (void)expected_offset;
      if (expected_count == 0UL) {
        if (first != 0UL || last != 0UL ||
            peer < nonempty_rank_count) {
          local_error = 1;
          break;
        }
        continue;
      }
      if (peer >= nonempty_rank_count || first > last ||
          (have_previous != 0 && previous_last >= first)) {
        local_error = 1;
        break;
      }
      next->rank_first_rep[peer] = first;
      previous_last = last;
      have_previous = 1;
    }
  }
  free(all_boundaries);
  all_boundaries = NULL;
  global_error = agree_directory_failure(mpi_active, local_error);
  if (global_error != 0) {
    FreeSymmetryRepresentativeDirectory(next);
    return -1;
  }
  next->ready = 1;
  *directory_out = next;
  return 0;
}

int SymmetryRepresentativeDirectoryReady(
    const struct SymmetryRepresentativeDirectory *directory)
{
  return directory != NULL && directory->ready == 1 &&
      directory->rank >= 0 && directory->nrank > 0 &&
      directory->rank < directory->nrank &&
      directory->nonempty_rank_count >= 0 &&
      directory->nonempty_rank_count <= directory->nrank &&
      directory->local_index != NULL &&
      (directory->nonempty_rank_count == 0
           ? directory->rank_first_rep == NULL
           : directory->rank_first_rep != NULL);
}

int SymmetryRepresentativeOwner(
    const struct SymmetryRepresentativeDirectory *directory,
    unsigned long int rep_state,
    int *owner)
{
  int next_owner = -1;
  int lo;
  int hi;
  if (SymmetryRepresentativeDirectoryReady(directory) == 0 ||
      owner == NULL) {
    return -1;
  }
  lo = 0;
  hi = directory->nonempty_rank_count;
  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    if (directory->rank_first_rep[mid] <= rep_state) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  if (directory->nonempty_rank_count > 0) {
    next_owner = lo == 0 ? 0 : lo - 1;
  }
  *owner = next_owner;
  return 0;
}

int SymmetryLookupDirectoryLocalRepresentative(
    const struct SymmetryRepresentativeDirectory *directory,
    unsigned long int rep_state,
    unsigned long int *local_index,
    unsigned long int *global_beta,
    double *norm,
    uint64_t *probe_count)
{
  if (SymmetryRepresentativeDirectoryReady(directory) == 0) return -1;
  return SymmetryLookupLocalRepresentative(
      directory->local_index, rep_state, local_index,
      global_beta, norm, probe_count);
}

static int directory_memory_add(size_t live_bytes,
                                size_t additional_bytes,
                                size_t byte_limit,
                                size_t *next_live)
{
  size_t next;
  if (next_live == NULL ||
      SymmetryCheckedSizeAdd(live_bytes, additional_bytes, &next) != 0 ||
      next > byte_limit) {
    return -1;
  }
  *next_live = next;
  return 0;
}

static int directory_exchange_owned_bytes(uint64_t count,
                                          size_t extent,
                                          int nrank,
                                          size_t *owned_bytes)
{
  size_t schedule_bytes;
  size_t payload_bytes;
  size_t total;
  size_t count_size;
  if (owned_bytes == NULL || nrank < 1 ||
      SymmetryCheckedU64ToSize(count, &count_size) != 0 ||
      SymmetryCheckedSizeMul((size_t)nrank, sizeof(uint64_t),
                             &schedule_bytes) != 0 ||
      SymmetryCheckedSizeMul(schedule_bytes, 2U, &schedule_bytes) != 0 ||
      SymmetryCheckedSizeMul(count_size, extent, &payload_bytes) != 0 ||
      SymmetryCheckedSizeAdd(schedule_bytes, payload_bytes, &total) != 0) {
    return -1;
  }
  *owned_bytes = total;
  return 0;
}

static int directory_exchange_workspace_bytes(
    int mpi_active,
    int nrank,
    const struct SymmetryMpiExchangeStats *stats,
    size_t *workspace_bytes)
{
#ifdef MPI
  size_t one_array;
  size_t request_count;
  size_t request_bytes;
  size_t total;
#endif
  if (stats == NULL || workspace_bytes == NULL || nrank < 1) return -1;
  *workspace_bytes = 0U;
#ifdef MPI
  if (mpi_active == 0 || nrank == 1) return 0;
  if (stats->used_chunked != 0) {
    if (SymmetryCheckedSizeMul((size_t)nrank, sizeof(uint64_t),
                               &one_array) != 0 ||
        SymmetryCheckedSizeMul((size_t)nrank, 2U, &request_count) != 0 ||
        SymmetryCheckedSizeMul(request_count, sizeof(MPI_Request),
                               &request_bytes) != 0 ||
        SymmetryCheckedSizeAdd(one_array, one_array, &total) != 0 ||
        SymmetryCheckedSizeAdd(total, request_bytes, &total) != 0) {
      return -1;
    }
  } else {
    if (SymmetryCheckedSizeMul((size_t)nrank, sizeof(int),
                               &one_array) != 0 ||
        SymmetryCheckedSizeMul(one_array, 4U, &total) != 0) {
      return -1;
    }
  }
  *workspace_bytes = total;
#else
  (void)mpi_active;
#endif
  return 0;
}

static int directory_record_exchange_peak(
    size_t caller_live,
    uint64_t result_count,
    size_t extent,
    int mpi_active,
    int nrank,
    const struct SymmetryMpiExchangeStats *exchange_stats,
    size_t *batch_peak)
{
  size_t owned_bytes;
  size_t workspace_bytes;
  size_t peak;
  if (batch_peak == NULL ||
      directory_exchange_owned_bytes(result_count, extent, nrank,
                                     &owned_bytes) != 0 ||
      directory_exchange_workspace_bytes(
          mpi_active, nrank, exchange_stats, &workspace_bytes) != 0 ||
      SymmetryCheckedSizeAdd(caller_live, owned_bytes, &peak) != 0 ||
      SymmetryCheckedSizeAdd(peak, workspace_bytes, &peak) != 0) {
    return -1;
  }
  if (peak > *batch_peak) *batch_peak = peak;
  return 0;
}

static int directory_batch_options(
    int mpi_active,
    int nrank,
    const struct SymmetryRepresentativeBatchOptions *options,
    struct SymmetryRepresentativeBatchOptions *next)
{
  uint64_t configured_entry_limit =
      (uint64_t)HPHI_SYMMETRY_EXCHANGE_MESSAGE_BYTES /
      (uint64_t)sizeof(struct SymmetryMpiLookupResponse);
  int local_error = 0;
  if (next == NULL || nrank < 1 || configured_entry_limit == 0U) return -1;
  memset(next, 0, sizeof(*next));
  next->corrupt_response_rank = -1;
  if (options != NULL) *next = *options;
  if (configured_entry_limit > (uint64_t)INT_MAX) {
    configured_entry_limit = (uint64_t)INT_MAX;
  }
  if ((next->force_chunked != 0 && next->force_chunked != 1) ||
      (next->debug_echo != 0 && next->debug_echo != 1) ||
      next->corrupt_response_rank < -1 ||
      next->corrupt_response_rank >= nrank ||
      next->chunk_limit > configured_entry_limit) {
    local_error = 1;
  }
  if (agree_directory_failure(mpi_active, local_error) != 0) return -1;
#ifdef MPI
  if (mpi_active != 0 && nrank > 1) {
    uint64_t chunk_min;
    uint64_t chunk_max;
    int force_min;
    int force_max;
    int echo_min;
    int echo_max;
    int corrupt_min;
    int corrupt_max;
    if (MPI_Allreduce(&next->chunk_limit, &chunk_min, 1, MPI_UINT64_T,
                      MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&next->chunk_limit, &chunk_max, 1, MPI_UINT64_T,
                      MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&next->force_chunked, &force_min, 1, MPI_INT,
                      MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&next->force_chunked, &force_max, 1, MPI_INT,
                      MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&next->debug_echo, &echo_min, 1, MPI_INT,
                      MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&next->debug_echo, &echo_max, 1, MPI_INT,
                      MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&next->corrupt_response_rank, &corrupt_min, 1, MPI_INT,
                      MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&next->corrupt_response_rank, &corrupt_max, 1, MPI_INT,
                      MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS) {
      return -1;
    }
    if (chunk_min != chunk_max || force_min != force_max ||
        echo_min != echo_max || corrupt_min != corrupt_max) {
      return -1;
    }
  }
#else
  (void)mpi_active;
#endif
  return 0;
}

static uint64_t directory_peer_count_max(const uint64_t *counts, int nrank)
{
  uint64_t maximum = 0U;
  int peer;
  if (counts == NULL || nrank < 1) return 0U;
  for (peer = 0; peer < nrank; peer++) {
    if (counts[peer] > maximum) maximum = counts[peer];
  }
  return maximum;
}

static int directory_include_exchange_stats(
    struct SymmetryRepresentativeBatchStats *batch,
    const struct SymmetryMpiExchangeStats *exchange)
{
  uint64_t next;
  if (batch == NULL || exchange == NULL ||
      checked_u64_add(batch->directory_exchange_send_messages,
                      exchange->send_messages, &next) != 0) {
    return -1;
  }
  batch->directory_exchange_send_messages = next;
  if (checked_u64_add(batch->directory_exchange_recv_messages,
                      exchange->recv_messages, &next) != 0) {
    return -1;
  }
  batch->directory_exchange_recv_messages = next;
  if (batch->directory_exchange_message_byte_limit == 0U ||
      exchange->message_byte_limit <
          batch->directory_exchange_message_byte_limit) {
    batch->directory_exchange_message_byte_limit =
        exchange->message_byte_limit;
  }
  if (exchange->max_message_bytes >
      batch->directory_exchange_max_message_bytes) {
    batch->directory_exchange_max_message_bytes =
        exchange->max_message_bytes;
  }
  if (exchange->used_chunked != 0) {
    batch->directory_exchange_used_chunked = 1;
  }
  return 0;
}

static int directory_publish_batch_stats(
    struct SymmetryRepresentativeDirectory *directory,
    const struct SymmetryRepresentativeBatchStats *call_stats,
    size_t batch_peak,
    struct SymmetryRepresentativeBatchStats *next_stats)
{
  uint64_t next;
  if (directory == NULL || call_stats == NULL || next_stats == NULL) {
    return -1;
  }
  *next_stats = directory->batch_stats;
#define DIRECTORY_ADD_STAT(field)                                             \
  do {                                                                        \
    if (checked_u64_add(next_stats->field, call_stats->field, &next) != 0) {  \
      return -1;                                                              \
    }                                                                         \
    next_stats->field = next;                                                 \
  } while (0)
  DIRECTORY_ADD_STAT(directory_batch_calls);
  DIRECTORY_ADD_STAT(directory_request_entries_sent);
  DIRECTORY_ADD_STAT(directory_request_entries_received);
  DIRECTORY_ADD_STAT(directory_found_entries);
  DIRECTORY_ADD_STAT(directory_not_found_entries);
  DIRECTORY_ADD_STAT(directory_lookup_probe_count);
  DIRECTORY_ADD_STAT(directory_exchange_send_messages);
  DIRECTORY_ADD_STAT(directory_exchange_recv_messages);
#undef DIRECTORY_ADD_STAT
  if (call_stats->directory_lookup_max_probe >
      next_stats->directory_lookup_max_probe) {
    next_stats->directory_lookup_max_probe =
        call_stats->directory_lookup_max_probe;
  }
  if (call_stats->directory_owner_peer_count_max >
      next_stats->directory_owner_peer_count_max) {
    next_stats->directory_owner_peer_count_max =
        call_stats->directory_owner_peer_count_max;
  }
  if (call_stats->directory_requester_peer_count_max >
      next_stats->directory_requester_peer_count_max) {
    next_stats->directory_requester_peer_count_max =
        call_stats->directory_requester_peer_count_max;
  }
  if (call_stats->directory_exchange_used_chunked != 0) {
    next_stats->directory_exchange_used_chunked = 1;
  }
  if (next_stats->directory_exchange_message_byte_limit == 0U ||
      (call_stats->directory_exchange_message_byte_limit != 0U &&
       call_stats->directory_exchange_message_byte_limit <
           next_stats->directory_exchange_message_byte_limit)) {
    next_stats->directory_exchange_message_byte_limit =
        call_stats->directory_exchange_message_byte_limit;
  }
  if (call_stats->directory_exchange_max_message_bytes >
      next_stats->directory_exchange_max_message_bytes) {
    next_stats->directory_exchange_max_message_bytes =
        call_stats->directory_exchange_max_message_bytes;
  }
  if (batch_peak > next_stats->directory_batch_temporary_peak_bytes) {
    next_stats->directory_batch_temporary_peak_bytes = batch_peak;
  }
  next_stats->directory_batch_memory_byte_limit =
      (size_t)HPHI_SYMMETRY_DIRECTORY_MEMORY_BYTES;
  return 0;
}

static int resolve_representative_batch_with_options(
    struct SymmetryRepresentativeDirectory *directory,
    const unsigned long int *request_keys,
    uint64_t request_count,
    unsigned long int *global_beta,
    double *norm,
    const struct SymmetryRepresentativeBatchOptions *options)
{
  struct SymmetryRepresentativeBatchOptions batch_options;
  struct SymmetryMpiExchangeOptions exchange_options;
  struct SymmetryMpiExchangeMemoryBudget exchange_budget;
  struct SymmetryMpiExchangeLayout request_send_layout;
  struct SymmetryMpiExchangeLayout response_send_layout;
  struct SymmetryMpiExchangeLayout response_receive_layout;
  struct SymmetryMpiUnsignedLongResult received_requests;
  struct SymmetryMpiLookupResponseResult received_responses;
  struct SymmetryMpiUnsignedLongResult received_echoes;
  struct SymmetryMpiExchangeStats request_exchange_stats;
  struct SymmetryMpiExchangeStats response_exchange_stats;
  struct SymmetryMpiExchangeStats echo_exchange_stats;
  struct SymmetryRepresentativeBatchStats call_stats;
  struct SymmetryRepresentativeBatchStats next_stats;
  struct SymmetryMpiLookupResponse *owner_responses = NULL;
  unsigned long int *packed_requests = NULL;
  uint64_t *send_counts = NULL;
  uint64_t *send_displacements = NULL;
  size_t byte_limit = (size_t)HPHI_SYMMETRY_DIRECTORY_MEMORY_BYTES;
  size_t input_request_count_size = 0U;
  size_t received_request_count_size = 0U;
  size_t schedule_bytes = 0U;
  size_t packed_bytes = 0U;
  size_t owner_response_bytes = 0U;
  size_t received_request_bytes = 0U;
  size_t received_response_bytes = 0U;
  size_t received_echo_bytes = 0U;
  size_t live_bytes = 0U;
  size_t batch_peak = 0U;
  int mpi_active;
  int comm_rank;
  int comm_size;
  int local_error = 0;
  int global_error;
  int owner = -1;
  int previous_owner = -1;
  int peer;
  uint64_t index;

  memset(&received_requests, 0, sizeof(received_requests));
  memset(&received_responses, 0, sizeof(received_responses));
  memset(&received_echoes, 0, sizeof(received_echoes));
  memset(&request_exchange_stats, 0, sizeof(request_exchange_stats));
  memset(&response_exchange_stats, 0, sizeof(response_exchange_stats));
  memset(&echo_exchange_stats, 0, sizeof(echo_exchange_stats));
  memset(&call_stats, 0, sizeof(call_stats));
  memset(&next_stats, 0, sizeof(next_stats));
  memset(&request_send_layout, 0, sizeof(request_send_layout));
  memset(&response_send_layout, 0, sizeof(response_send_layout));
  memset(&response_receive_layout, 0, sizeof(response_receive_layout));

  if (directory_mpi_context(&mpi_active, &comm_rank, &comm_size) != 0) {
    return -1;
  }
  if (SymmetryRepresentativeDirectoryReady(directory) == 0 ||
      directory->rank != comm_rank || directory->nrank != comm_size ||
      (request_count > 0U &&
       (request_keys == NULL || global_beta == NULL || norm == NULL)) ||
      SymmetryCheckedU64ToSize(
          request_count, &input_request_count_size) != 0 ||
      input_request_count_size > SIZE_MAX / sizeof(unsigned long int) ||
      input_request_count_size > SIZE_MAX / sizeof(double)) {
    local_error = 1;
  }
  if (local_error == 0) {
    for (index = 1U; index < request_count; index++) {
      if (request_keys[index - 1U] >= request_keys[index]) {
        local_error = 1;
        break;
      }
    }
  }
  global_error = agree_directory_failure(mpi_active, local_error);
  if (global_error != 0) return -1;
  if (directory_batch_options(
          mpi_active, directory->nrank, options, &batch_options) != 0) {
    return -1;
  }
  exchange_options.chunk_limit = batch_options.chunk_limit;
  exchange_options.force_chunked = batch_options.force_chunked;

  call_stats.directory_batch_calls = 1U;
  call_stats.directory_batch_memory_byte_limit = byte_limit;
  if (directory->dim == 0UL) {
    call_stats.directory_not_found_entries = request_count;
    if (directory_publish_batch_stats(
            directory, &call_stats, 0U, &next_stats) != 0) {
      local_error = 1;
    }
    global_error = agree_directory_failure(mpi_active, local_error);
    if (global_error != 0) return -1;
    for (index = 0U; index < request_count; index++) {
      global_beta[index] = 0UL;
      norm[index] = 0.0;
    }
    directory->batch_stats = next_stats;
    return 0;
  }

  if (SymmetryCheckedSizeMul((size_t)directory->nrank, sizeof(uint64_t),
                             &schedule_bytes) != 0 ||
      SymmetryCheckedSizeMul(input_request_count_size,
                             sizeof(unsigned long int),
                             &packed_bytes) != 0 ||
      SymmetryCheckedSizeMul(schedule_bytes, 2U, &live_bytes) != 0 ||
      directory_memory_add(live_bytes, packed_bytes, byte_limit,
                           &live_bytes) != 0) {
    local_error = 1;
  }
  global_error = agree_directory_failure(mpi_active, local_error);
  if (global_error != 0) goto fail;
  send_counts = (uint64_t *)calloc((size_t)directory->nrank,
                                   sizeof(*send_counts));
  send_displacements = (uint64_t *)malloc(schedule_bytes);
  if (packed_bytes > 0U) {
    packed_requests = (unsigned long int *)malloc(packed_bytes);
  }
  if (send_counts == NULL || send_displacements == NULL ||
      (packed_bytes > 0U && packed_requests == NULL)) {
    local_error = 1;
  }
  global_error = agree_directory_failure(mpi_active, local_error);
  if (global_error != 0) goto fail;
  batch_peak = live_bytes;

  for (index = 0U; index < request_count; index++) {
    if (SymmetryRepresentativeOwner(
            directory, request_keys[index], &owner) != 0 ||
        owner < 0 || owner >= directory->nrank ||
        owner < previous_owner ||
        send_counts[owner] == UINT64_MAX) {
      local_error = 1;
      break;
    }
    send_counts[owner]++;
    packed_requests[index] = request_keys[index];
    previous_owner = owner;
  }
  {
    uint64_t prefix = 0U;
    for (peer = 0; peer < directory->nrank; peer++) {
      send_displacements[peer] = prefix;
      if (checked_u64_add(prefix, send_counts[peer], &prefix) != 0) {
        local_error = 1;
        break;
      }
    }
    if (prefix != request_count) local_error = 1;
  }
  global_error = agree_directory_failure(mpi_active, local_error);
  if (global_error != 0) goto fail;

  request_send_layout.nrank = directory->nrank;
  request_send_layout.count = request_count;
  request_send_layout.counts = send_counts;
  request_send_layout.displacements = send_displacements;
  exchange_budget.live_bytes = live_bytes;
  exchange_budget.byte_limit = byte_limit;
  if (SymmetryMpiExchangeUnsignedLongsWithBudget(
          packed_requests, &request_send_layout,
          directory->rank, directory->nrank, &exchange_options,
          &exchange_budget, &received_requests,
          &request_exchange_stats) != 0) {
    goto fail;
  }
  if (directory_record_exchange_peak(
          live_bytes, received_requests.count, sizeof(unsigned long int),
          mpi_active, directory->nrank, &request_exchange_stats,
          &batch_peak) != 0 ||
      directory_exchange_owned_bytes(
          received_requests.count, sizeof(unsigned long int),
          directory->nrank, &received_request_bytes) != 0 ||
      directory_memory_add(live_bytes, received_request_bytes, byte_limit,
                           &live_bytes) != 0 ||
      SymmetryCheckedU64ToSize(received_requests.count,
                              &received_request_count_size) != 0 ||
      SymmetryCheckedSizeMul(received_request_count_size,
                             sizeof(*owner_responses),
                             &owner_response_bytes) != 0 ||
      directory_memory_add(live_bytes, owner_response_bytes, byte_limit,
                           &live_bytes) != 0) {
    local_error = 1;
  }
  global_error = agree_directory_failure(mpi_active, local_error);
  if (global_error != 0) goto fail;
  if (owner_response_bytes > 0U) {
    owner_responses =
        (struct SymmetryMpiLookupResponse *)malloc(owner_response_bytes);
    if (owner_responses == NULL) local_error = 1;
  }
  global_error = agree_directory_failure(mpi_active, local_error);
  if (global_error != 0) goto fail;
  if (live_bytes > batch_peak) batch_peak = live_bytes;

  for (index = 0U; index < received_requests.count; index++) {
    uint64_t probe_count = 0U;
    int request_owner = -1;
    if (SymmetryRepresentativeOwner(
            directory, received_requests.entries[index],
            &request_owner) != 0 ||
        request_owner != directory->rank ||
        SymmetryLookupDirectoryLocalRepresentative(
            directory, received_requests.entries[index], NULL,
            &owner_responses[index].global_beta,
            &owner_responses[index].norm, &probe_count) != 0 ||
        checked_u64_add(call_stats.directory_lookup_probe_count,
                        probe_count,
                        &call_stats.directory_lookup_probe_count) != 0) {
      local_error = 1;
      break;
    }
    if (probe_count > call_stats.directory_lookup_max_probe) {
      call_stats.directory_lookup_max_probe = probe_count;
    }
    if (owner_responses[index].global_beta == 0UL) {
      if (owner_responses[index].norm != 0.0 ||
          signbit(owner_responses[index].norm) ||
          call_stats.directory_not_found_entries == UINT64_MAX) {
        local_error = 1;
        break;
      }
      call_stats.directory_not_found_entries++;
    } else {
      if (!isfinite(owner_responses[index].norm) ||
          owner_responses[index].norm <= 0.0 ||
          call_stats.directory_found_entries == UINT64_MAX) {
        local_error = 1;
        break;
      }
      call_stats.directory_found_entries++;
    }
  }
  if (batch_options.corrupt_response_rank == directory->rank) {
    if (received_requests.count == 0U) {
      local_error = 1;
    } else {
      owner_responses[0].global_beta = ULONG_MAX;
      owner_responses[0].norm = -1.0;
    }
  }
  global_error = agree_directory_failure(mpi_active, local_error);
  if (global_error != 0) goto fail;

  response_send_layout.nrank = directory->nrank;
  response_send_layout.count = received_requests.count;
  response_send_layout.counts = received_requests.counts;
  response_send_layout.displacements = received_requests.displacements;
  response_receive_layout = request_send_layout;
  exchange_budget.live_bytes = live_bytes;
  if (SymmetryMpiExchangeLookupResponsesKnownLayoutWithBudget(
          owner_responses, &response_send_layout, &response_receive_layout,
          directory->rank, directory->nrank, &exchange_options,
          &exchange_budget, &received_responses,
          &response_exchange_stats) != 0) {
    goto fail;
  }
  if (directory_record_exchange_peak(
          live_bytes, received_responses.count,
          sizeof(struct SymmetryMpiLookupResponse),
          mpi_active, directory->nrank, &response_exchange_stats,
          &batch_peak) != 0 ||
      directory_exchange_owned_bytes(
          received_responses.count,
          sizeof(struct SymmetryMpiLookupResponse),
          directory->nrank, &received_response_bytes) != 0 ||
      directory_memory_add(live_bytes, received_response_bytes, byte_limit,
                           &live_bytes) != 0) {
    local_error = 1;
  }
  global_error = agree_directory_failure(mpi_active, local_error);
  if (global_error != 0) goto fail;
  free(owner_responses);
  owner_responses = NULL;
  live_bytes -= owner_response_bytes;
  owner_response_bytes = 0U;

  if (batch_options.debug_echo != 0) {
    exchange_budget.live_bytes = live_bytes;
    if (SymmetryMpiExchangeUnsignedLongEchoesKnownLayoutWithBudget(
            received_requests.entries, &response_send_layout,
            &response_receive_layout, directory->rank, directory->nrank,
            &exchange_options, &exchange_budget, &received_echoes,
            &echo_exchange_stats) != 0) {
      goto fail;
    }
    if (directory_record_exchange_peak(
            live_bytes, received_echoes.count, sizeof(unsigned long int),
            mpi_active, directory->nrank, &echo_exchange_stats,
            &batch_peak) != 0 ||
        directory_exchange_owned_bytes(
            received_echoes.count, sizeof(unsigned long int),
            directory->nrank, &received_echo_bytes) != 0 ||
        directory_memory_add(live_bytes, received_echo_bytes, byte_limit,
                             &live_bytes) != 0) {
      local_error = 1;
    }
    global_error = agree_directory_failure(mpi_active, local_error);
    if (global_error != 0) goto fail;
  }

  if (received_responses.count != request_count ||
      received_responses.nrank != directory->nrank ||
      (batch_options.debug_echo != 0 &&
       (received_echoes.count != request_count ||
        received_echoes.nrank != directory->nrank))) {
    local_error = 1;
  }
  for (peer = 0; local_error == 0 && peer < directory->nrank; peer++) {
    if (received_responses.counts[peer] != send_counts[peer] ||
        received_responses.displacements[peer] != send_displacements[peer] ||
        (batch_options.debug_echo != 0 &&
         (received_echoes.counts[peer] != send_counts[peer] ||
          received_echoes.displacements[peer] !=
              send_displacements[peer]))) {
      local_error = 1;
    }
  }
  for (index = 0U; local_error == 0 && index < request_count; index++) {
    unsigned long int block_offset;
    unsigned long int block_count;
    const struct SymmetryMpiLookupResponse *response =
        &received_responses.entries[index];
    if (SymmetryRepresentativeOwner(
            directory, request_keys[index], &owner) != 0 ||
        directory_block_range(directory->dim, owner, directory->nrank,
                              &block_offset, &block_count) != 0 ||
        (batch_options.debug_echo != 0 &&
         received_echoes.entries[index] != request_keys[index])) {
      local_error = 1;
      break;
    }
    if (response->global_beta == 0UL) {
      if (response->norm != 0.0 || signbit(response->norm)) {
        local_error = 1;
      }
    } else if (block_count == 0UL ||
               response->global_beta <= block_offset ||
               response->global_beta > block_offset + block_count ||
               !isfinite(response->norm) || response->norm <= 0.0) {
      local_error = 1;
    }
  }

  call_stats.directory_request_entries_sent = request_count;
  call_stats.directory_request_entries_received = received_requests.count;
  call_stats.directory_owner_peer_count_max =
      directory_peer_count_max(received_requests.counts, directory->nrank);
  call_stats.directory_requester_peer_count_max =
      directory_peer_count_max(send_counts, directory->nrank);
  if (directory_include_exchange_stats(
          &call_stats, &request_exchange_stats) != 0 ||
      directory_include_exchange_stats(
          &call_stats, &response_exchange_stats) != 0 ||
      (batch_options.debug_echo != 0 &&
       directory_include_exchange_stats(
           &call_stats, &echo_exchange_stats) != 0) ||
      directory_publish_batch_stats(
          directory, &call_stats, batch_peak, &next_stats) != 0) {
    local_error = 1;
  }
  global_error = agree_directory_failure(mpi_active, local_error);
  if (global_error != 0) goto fail;

  for (index = 0U; index < request_count; index++) {
    global_beta[index] = received_responses.entries[index].global_beta;
    norm[index] = received_responses.entries[index].norm;
  }
  directory->batch_stats = next_stats;
  FreeSymmetryMpiUnsignedLongResult(&received_echoes);
  FreeSymmetryMpiLookupResponseResult(&received_responses);
  FreeSymmetryMpiUnsignedLongResult(&received_requests);
  free(packed_requests);
  free(send_counts);
  free(send_displacements);
  return 0;

fail:
  free(owner_responses);
  FreeSymmetryMpiUnsignedLongResult(&received_echoes);
  FreeSymmetryMpiLookupResponseResult(&received_responses);
  FreeSymmetryMpiUnsignedLongResult(&received_requests);
  free(packed_requests);
  free(send_counts);
  free(send_displacements);
  return -1;
}

int SymmetryResolveRepresentativeBatchWithOptions(
    struct SymmetryRepresentativeDirectory *directory,
    const unsigned long int *request_keys,
    uint64_t request_count,
    unsigned long int *global_beta,
    double *norm,
    const struct SymmetryRepresentativeBatchOptions *options)
{
  int status;
  StartTimer(1123);
  status = resolve_representative_batch_with_options(
      directory, request_keys, request_count, global_beta, norm, options);
  StopTimer(1123);
  return status;
}

int SymmetryResolveRepresentativeBatch(
    struct SymmetryRepresentativeDirectory *directory,
    const unsigned long int *request_keys,
    uint64_t request_count,
    unsigned long int *global_beta,
    double *norm)
{
  return SymmetryResolveRepresentativeBatchWithOptions(
      directory, request_keys, request_count, global_beta, norm, NULL);
}

int GetSymmetryRepresentativeDirectoryBatchStats(
    const struct SymmetryRepresentativeDirectory *directory,
    struct SymmetryRepresentativeBatchStats *stats)
{
  struct SymmetryRepresentativeBatchStats next;
  if (SymmetryRepresentativeDirectoryReady(directory) == 0 ||
      stats == NULL) {
    return -1;
  }
  next = directory->batch_stats;
  *stats = next;
  return 0;
}

int GetSymmetryRepresentativeDirectoryInfo(
    const struct SymmetryRepresentativeDirectory *directory,
    struct SymmetryRepresentativeDirectoryInfo *info)
{
  struct SymmetryRepresentativeDirectoryInfo next;
  if (SymmetryRepresentativeDirectoryReady(directory) == 0 ||
      info == NULL ||
      (size_t)directory->nonempty_rank_count >
          SIZE_MAX / sizeof(unsigned long int)) {
    return -1;
  }
  next.dim = directory->dim;
  next.local_offset = directory->local_offset;
  next.local_dim = directory->local_dim;
  next.rank = directory->rank;
  next.nrank = directory->nrank;
  next.nonempty_rank_count = directory->nonempty_rank_count;
  next.splitter_bytes =
      (size_t)directory->nonempty_rank_count * sizeof(unsigned long int);
  *info = next;
  return 0;
}

int GetSymmetryRepresentativeDirectoryLocalIndexStats(
    const struct SymmetryRepresentativeDirectory *directory,
    struct SymmetryLocalRepresentativeIndexStats *stats)
{
  if (SymmetryRepresentativeDirectoryReady(directory) == 0) return -1;
  return GetSymmetryLocalRepresentativeIndexStats(
      directory->local_index, stats);
}
