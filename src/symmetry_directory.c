#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

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
