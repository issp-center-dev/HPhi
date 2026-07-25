#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "symmetry_basis.h"
#include "symmetry_directory.h"

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
      (mpi_active != 0 &&
       (rank != comm_rank || nrank != comm_size)) ||
      (mpi_active == 0 && (rank != 0 || nrank != 1)) ||
      validate_directory_layout(dim, local_offset, local_dim,
                                rank_offsets, rank, nrank) != 0) {
    local_error = 1;
  }
  if (local_error == 0 &&
      BuildSymmetryLocalRepresentativeIndex(
          local_basis, local_dim, local_capacity, local_offset,
          &local_index) != 0) {
    local_error = 1;
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
