#include <inttypes.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "DefCommon.h"
#include "symmetry_basis.h"
#include "symmetry_distribution.h"
#include "symmetry_mpi_exchange.h"

#define SYMMETRY_SAMPLE_COUNT_DEFAULT UINT64_C(64)
#define SYMMETRY_SAMPLE_COUNT_GLOBAL_CAP UINT64_C(131072)

struct SymmetryMergeNode {
  int source;
  uint64_t offset;
};

static int agree_distribution_failure(
    int mpi_active,
    int rank,
    const char *stage,
    const char *invariant,
    int local_error)
{
  int global_error = SymmetryMpiAgreeError(mpi_active, local_error);
  if (global_error != 0 && rank == 0) {
    if (global_error < 0) {
      fprintf(stderr,
              "HPhi symmetry distribution %s failed: "
              "MPI error agreement failed while checking %s.\n",
              stage, invariant);
    } else {
      fprintf(stderr,
              "HPhi symmetry distribution %s failed: %s.\n",
              stage, invariant);
    }
    fflush(stderr);
  }
  return global_error;
}

static int agree_distribution_limit_failure(
    int mpi_active,
    int rank,
    const char *stage,
    const char *invariant,
    int local_error,
    uint64_t observed,
    uint64_t limit)
{
  int global_error = SymmetryMpiAgreeError(mpi_active, local_error);
  uint64_t maximum_observed = observed;
  if (global_error == 0) return 0;
#ifdef MPI
  if (global_error > 0 && mpi_active != FALSE &&
      MPI_Allreduce(&observed, &maximum_observed, 1, MPI_UINT64_T,
                    MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS) {
    global_error = -1;
  }
#else
  (void)mpi_active;
#endif
  if (rank == 0) {
    if (global_error < 0) {
      fprintf(stderr,
              "HPhi symmetry distribution %s failed: "
              "MPI error agreement failed while checking %s.\n",
              stage, invariant);
    } else {
      fprintf(stderr,
              "HPhi symmetry distribution %s failed: %s "
              "(maximum observed=%" PRIu64 ", limit=%" PRIu64 ").\n",
              stage, invariant, maximum_observed, limit);
    }
    fflush(stderr);
  }
  return global_error;
}

#ifdef MPI
static void report_distribution_operation_failure(
    int rank,
    const char *stage,
    const char *operation)
{
  fprintf(stderr,
          "HPhi symmetry distribution %s failed on rank %d: %s.\n",
          stage, rank, operation);
  fflush(stderr);
}
#endif

static int checked_u64_add(uint64_t lhs, uint64_t rhs, uint64_t *result)
{
  if (result == NULL || lhs > UINT64_MAX - rhs) return -1;
  *result = lhs + rhs;
  return 0;
}

static int checked_u64_mul(uint64_t lhs, uint64_t rhs, uint64_t *result)
{
  if (result == NULL || (lhs != 0U && rhs > UINT64_MAX / lhs)) return -1;
  *result = lhs * rhs;
  return 0;
}

static int checked_size_mul(size_t lhs, size_t rhs, size_t *result)
{
  if (result == NULL || (lhs != 0U && rhs > SIZE_MAX / lhs)) return -1;
  *result = lhs * rhs;
  return 0;
}

static int checked_u64_to_size(uint64_t value, size_t *result)
{
  if (result == NULL || value > (uint64_t)SIZE_MAX) return -1;
  *result = (size_t)value;
  return 0;
}

static int checked_ulong_to_u64(unsigned long int value, uint64_t *result)
{
  uint64_t converted;
  if (result == NULL) return -1;
  converted = (uint64_t)value;
  if ((unsigned long int)converted != value) return -1;
  *result = converted;
  return 0;
}

static int symmetry_basis_run_is_valid(const struct SymmetryBasisRun *run)
{
  if (run == NULL) return FALSE;
  if (run->entries == NULL) {
    return run->count == 0UL && run->capacity == 0UL;
  }
  if (run->count == ULONG_MAX || run->capacity == 0UL) return FALSE;
  return run->capacity >= run->count + 1UL;
}

int SymmetryBlockRange(
    unsigned long int dim,
    int rank,
    int nrank,
    unsigned long int *offset,
    unsigned long int *count)
{
  unsigned long int base;
  unsigned long int remainder;
  unsigned long int unsigned_rank;
  if (rank < 0 || nrank < 1 || rank >= nrank ||
      offset == NULL || count == NULL) {
    return -1;
  }
  base = dim / (unsigned long int)nrank;
  remainder = dim % (unsigned long int)nrank;
  unsigned_rank = (unsigned long int)rank;
  *count = base + (unsigned_rank < remainder ? 1UL : 0UL);
  *offset = base * unsigned_rank +
      (unsigned_rank < remainder ? unsigned_rank : remainder);
  return 0;
}

int SymmetryCompareBasisRepState(const void *lhs, const void *rhs)
{
  const struct SymmetryBasisVector *a =
      (const struct SymmetryBasisVector *)lhs;
  const struct SymmetryBasisVector *b =
      (const struct SymmetryBasisVector *)rhs;
  if (a->rep_state < b->rep_state) return -1;
  if (a->rep_state > b->rep_state) return 1;
  return 0;
}

void FreeSymmetryBasisRun(struct SymmetryBasisRun *run)
{
  if (run == NULL) return;
  free(run->entries);
  run->entries = NULL;
  run->count = 0UL;
  run->capacity = 0UL;
}

void FreeSymmetryBasisOwnership(
    struct SymmetryBasisOwnership *ownership)
{
  if (ownership == NULL) return;
  free(ownership->rank_offsets);
  ownership->rank_offsets = NULL;
  ownership->dim = 0UL;
  ownership->local_offset = 0UL;
  ownership->local_dim = 0UL;
}

static unsigned long int symmetry_splitter_bucket(
    const unsigned long int *splitters,
    int splitter_count,
    unsigned long int key)
{
  int begin = 0;
  int end = splitter_count;
  while (begin < end) {
    int middle = begin + (end - begin) / 2;
    /*
     * upper_bound equality rule: an entry equal to a splitter belongs to
     * the bucket on the splitter's right.
     */
    if (splitters[middle] <= key) {
      begin = middle + 1;
    } else {
      end = middle;
    }
  }
  return (unsigned long int)begin;
}

static int build_regular_samples(
    const struct SymmetryBasisRun *run,
    uint64_t samples_per_nonempty_rank,
    struct SymmetryBasisVector **samples,
    uint64_t *sample_count,
    uint64_t *local_gap)
{
  uint64_t count;
  uint64_t selected;
  uint64_t base_step;
  uint64_t step_remainder;
  uint64_t offset = 0U;
  uint64_t remainder_accumulator = 0U;
  uint64_t previous_position = 0U;
  uint64_t gap_max = 0U;
  uint64_t index;
  size_t sample_bytes;

  if (samples == NULL || sample_count == NULL || local_gap == NULL ||
      checked_ulong_to_u64(run->count, &count) != 0) {
    return -1;
  }
  *samples = NULL;
  *sample_count = 0U;
  *local_gap = 0U;
  if (count == 0U) return 0;
  if (samples_per_nonempty_rank == 0U) return -1;

  selected = count < samples_per_nonempty_rank
      ? count : samples_per_nonempty_rank;
  if (checked_size_mul((size_t)selected, sizeof(**samples),
                       &sample_bytes) != 0) {
    return -1;
  }
  *samples = (struct SymmetryBasisVector *)malloc(sample_bytes);
  if (*samples == NULL) return -1;

  base_step = count / selected;
  step_remainder = count % selected;
  for (index = 0U; index < selected; index++) {
    uint64_t position;
    uint64_t gap;
    if (checked_u64_add(offset, 1U, &position) != 0 ||
        position == 0U || position > count ||
        position > (uint64_t)ULONG_MAX) {
      free(*samples);
      *samples = NULL;
      return -1;
    }
    memcpy(&(*samples)[index],
           &run->entries[(unsigned long int)position],
           sizeof(**samples));
    gap = position - previous_position;
    if (gap > gap_max) gap_max = gap;
    previous_position = position;
    if (index + 1U < selected) {
      if (checked_u64_add(offset, base_step, &offset) != 0 ||
          checked_u64_add(remainder_accumulator, step_remainder,
                          &remainder_accumulator) != 0) {
        free(*samples);
        *samples = NULL;
        return -1;
      }
      if (remainder_accumulator >= selected) {
        if (checked_u64_add(offset, 1U, &offset) != 0) {
          free(*samples);
          *samples = NULL;
          return -1;
        }
        remainder_accumulator -= selected;
      }
    }
  }
  {
    uint64_t endpoint;
    uint64_t suffix_gap;
    if (checked_u64_add(count, 1U, &endpoint) != 0) {
      free(*samples);
      *samples = NULL;
      return -1;
    }
    suffix_gap = endpoint - previous_position;
    if (suffix_gap > gap_max) gap_max = suffix_gap;
  }
  *sample_count = selected;
  *local_gap = gap_max;
  return 0;
}

static uint64_t fnv1a_bytes(uint64_t digest,
                            const void *value,
                            size_t bytes)
{
  const unsigned char *data = (const unsigned char *)value;
  size_t index;
  for (index = 0U; index < bytes; index++) {
    digest ^= (uint64_t)data[index];
    digest *= UINT64_C(1099511628211);
  }
  return digest;
}

static uint64_t digest_splitters(const unsigned long int *splitters,
                                 int splitter_count)
{
  uint64_t digest = UINT64_C(14695981039346656037);
  int index;
  for (index = 0; index < splitter_count; index++) {
    digest = fnv1a_bytes(digest, &splitters[index],
                         sizeof(splitters[index]));
  }
  return digest;
}

static uint64_t digest_basis_run(const struct SymmetryBasisVector *entries,
                                 uint64_t count)
{
  uint64_t digest = UINT64_C(14695981039346656037);
  uint64_t index;
  digest = fnv1a_bytes(digest, &count, sizeof(count));
  for (index = 0U; index < count; index++) {
    const struct SymmetryBasisVector *entry = &entries[index];
    digest = fnv1a_bytes(digest, &entry->rep_state,
                         sizeof(entry->rep_state));
    digest = fnv1a_bytes(digest, &entry->orbit_size,
                         sizeof(entry->orbit_size));
    digest = fnv1a_bytes(digest, &entry->stabilizer_size,
                         sizeof(entry->stabilizer_size));
    digest = fnv1a_bytes(digest, &entry->norm, sizeof(entry->norm));
    digest = fnv1a_bytes(digest, &entry->stabilizer_character_sum,
                         sizeof(entry->stabilizer_character_sum));
    digest = fnv1a_bytes(digest, &entry->diagonal,
                         sizeof(entry->diagonal));
  }
  return digest;
}

static int merge_node_is_less(
    const struct SymmetryMergeNode *lhs,
    const struct SymmetryMergeNode *rhs,
    const struct SymmetryMpiExchangeResult *result)
{
  const struct SymmetryBasisVector *lhs_entry =
      &result->entries[result->displacements[lhs->source] + lhs->offset];
  const struct SymmetryBasisVector *rhs_entry =
      &result->entries[result->displacements[rhs->source] + rhs->offset];
  if (lhs_entry->rep_state < rhs_entry->rep_state) return TRUE;
  if (lhs_entry->rep_state > rhs_entry->rep_state) return FALSE;
  return lhs->source < rhs->source;
}

static void merge_heap_push(
    struct SymmetryMergeNode *heap,
    size_t *heap_count,
    struct SymmetryMergeNode node,
    const struct SymmetryMpiExchangeResult *result)
{
  size_t child = *heap_count;
  heap[child] = node;
  (*heap_count)++;
  while (child > 0) {
    size_t parent = (child - 1U) / 2U;
    struct SymmetryMergeNode temporary;
    if (merge_node_is_less(&heap[parent], &heap[child], result)) break;
    temporary = heap[parent];
    heap[parent] = heap[child];
    heap[child] = temporary;
    child = parent;
  }
}

static struct SymmetryMergeNode merge_heap_pop(
    struct SymmetryMergeNode *heap,
    size_t *heap_count,
    const struct SymmetryMpiExchangeResult *result)
{
  struct SymmetryMergeNode root = heap[0];
  size_t parent = 0U;
  (*heap_count)--;
  if (*heap_count == 0U) return root;
  heap[0] = heap[*heap_count];
  for (;;) {
    size_t left = parent * 2U + 1U;
    size_t right = left + 1U;
    size_t smallest;
    struct SymmetryMergeNode temporary;
    if (left >= *heap_count) break;
    smallest = left;
    if (right < *heap_count &&
        merge_node_is_less(&heap[right], &heap[left], result)) {
      smallest = right;
    }
    if (merge_node_is_less(&heap[parent], &heap[smallest], result)) break;
    temporary = heap[parent];
    heap[parent] = heap[smallest];
    heap[smallest] = temporary;
    parent = smallest;
  }
  return root;
}

static int calculate_range_temporary_peak(
    uint64_t input_count,
    uint64_t receive_count,
    int nrank,
    uint64_t *peak_bytes)
{
  uint64_t input_elements;
  uint64_t input_bytes;
  uint64_t receive_bytes;
  uint64_t output_elements;
  uint64_t output_bytes;
  uint64_t schedule_bytes;
  uint64_t heap_bytes;
  uint64_t merge_peak;
  uint64_t exchange_overhead_per_rank;
  uint64_t exchange_overhead;
  uint64_t exchange_peak;
  uint64_t chunk_transport_overhead =
      UINT64_C(9) * (uint64_t)sizeof(uint64_t);
#ifdef MPI
  if (checked_u64_add(
          chunk_transport_overhead,
          UINT64_C(2) * (uint64_t)sizeof(MPI_Request),
          &chunk_transport_overhead) != 0) {
    return -1;
  }
#endif
  exchange_overhead_per_rank =
      UINT64_C(7) * (uint64_t)sizeof(uint64_t) +
      UINT64_C(4) * (uint64_t)sizeof(int);
  if (chunk_transport_overhead > exchange_overhead_per_rank) {
    exchange_overhead_per_rank = chunk_transport_overhead;
  }
  if (nrank < 1 ||
      checked_u64_add(input_count, 1U, &input_elements) != 0 ||
      checked_u64_mul(input_elements,
                      (uint64_t)sizeof(struct SymmetryBasisVector),
                      &input_bytes) != 0 ||
      checked_u64_mul(receive_count,
                      (uint64_t)sizeof(struct SymmetryBasisVector),
                      &receive_bytes) != 0 ||
      checked_u64_add(receive_count, 1U, &output_elements) != 0 ||
      checked_u64_mul(output_elements,
                      (uint64_t)sizeof(struct SymmetryBasisVector),
                      &output_bytes) != 0 ||
      checked_u64_mul((uint64_t)nrank,
                      UINT64_C(5) * (uint64_t)sizeof(uint64_t),
                      &schedule_bytes) != 0 ||
      checked_u64_mul((uint64_t)nrank,
                      (uint64_t)sizeof(struct SymmetryMergeNode),
                      &heap_bytes) != 0 ||
      checked_u64_add(input_bytes, receive_bytes, &merge_peak) != 0 ||
      checked_u64_add(merge_peak, output_bytes, &merge_peak) != 0 ||
      checked_u64_add(merge_peak, schedule_bytes, &merge_peak) != 0 ||
      checked_u64_add(merge_peak, heap_bytes, &merge_peak) != 0 ||
      checked_u64_mul((uint64_t)nrank, exchange_overhead_per_rank,
                      &exchange_overhead) != 0 ||
      checked_u64_add(input_bytes, receive_bytes,
                      &exchange_peak) != 0) {
    return -1;
  }
  if (checked_u64_add(exchange_peak, exchange_overhead,
                      &exchange_peak) != 0) {
    return -1;
  }
  *peak_bytes = merge_peak > exchange_peak ? merge_peak : exchange_peak;
  return 0;
}

int SymmetrySampleSortBasisRun(
    struct SymmetryBasisRun *run,
    int rank,
    int nrank,
    struct SymmetryBasisDistributionStats *stats)
{
  struct SymmetryBasisDistributionStats next_stats;
  struct SymmetryBasisVector *local_samples = NULL;
  struct SymmetryBasisVector *merged_entries = NULL;
  struct SymmetryMergeNode *merge_heap = NULL;
  struct SymmetryMpiExchangeResult sample_result;
  struct SymmetryMpiExchangeResult range_result;
  struct SymmetryMpiExchangeLayout sample_layout;
  struct SymmetryMpiExchangeLayout range_layout;
  uint64_t *all_counts = NULL;
  uint64_t *all_gaps = NULL;
  uint64_t *all_last_keys = NULL;
  uint64_t *sample_send_counts = NULL;
  uint64_t *sample_send_displacements = NULL;
  uint64_t *bucket_sample_counts = NULL;
  uint64_t *send_counts = NULL;
  uint64_t *send_displacements = NULL;
  uint64_t *preflight_recv_counts = NULL;
  unsigned long int *splitters = NULL;
  uint64_t local_count = 0U;
  uint64_t global_count = 0U;
  uint64_t local_sample_count = 0U;
  uint64_t global_sample_count = 0U;
  uint64_t local_gap = 0U;
  uint64_t global_gap_max = 0U;
  uint64_t global_gap_sum = 0U;
  uint64_t sample_limit = 0U;
  uint64_t nonempty_ranks = 0U;
  uint64_t receive_count = 0U;
  uint64_t bucket_upper_bound = 0U;
  uint64_t sample_peak = 0U;
  uint64_t temporary_peak = 0U;
  uint64_t merged_count = 0U;
  uint64_t range_global_count = 0U;
  uint64_t range_max_count = 0U;
  size_t rank_array_bytes = 0U;
  size_t splitter_bytes = 0U;
  size_t merged_bytes = 0U;
  size_t heap_bytes = 0U;
  size_t sample_peak_size = 0U;
  int splitter_count = 0;
  int mpi_active;
  int local_error = 0;
  int global_error;
  int peer;
  size_t merge_heap_count = 0U;
#ifdef MPI
  int comm_rank = 0;
  int comm_size = 1;
#endif

  memset(&next_stats, 0, sizeof(next_stats));
  memset(&sample_result, 0, sizeof(sample_result));
  memset(&range_result, 0, sizeof(range_result));
  if (stats != NULL) memset(stats, 0, sizeof(*stats));
  mpi_active = SymmetryMpiCollectivesActive();

  if (rank < 0 || nrank < 1 || rank >= nrank) {
    local_error = 1;
  }
#ifdef MPI
  if (mpi_active != FALSE) {
    if (MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank) != MPI_SUCCESS ||
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size) != MPI_SUCCESS ||
        comm_rank != rank || comm_size != nrank) {
      local_error = 1;
    }
  } else if (rank != 0 || nrank != 1) {
    local_error = 1;
  }
#else
  if (rank != 0 || nrank != 1) local_error = 1;
#endif
  if (nrank > 0 &&
      checked_size_mul((size_t)nrank, sizeof(uint64_t),
                       &rank_array_bytes) != 0) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort", "invalid rank or communicator",
      local_error);
  if (global_error != 0) return -1;
  if (symmetry_basis_run_is_valid(run) != TRUE ||
      checked_ulong_to_u64(run != NULL ? run->count : 0UL,
                           &local_count) != 0 ||
      local_count > (uint64_t)(SIZE_MAX /
          sizeof(struct SymmetryBasisVector)) - 1U ||
      HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES == 0U) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort",
      "invalid input run or distribution memory configuration",
      local_error);
  if (global_error != 0) return -1;

  if (local_count > 1U) {
    qsort(run->entries + 1, (size_t)local_count,
          sizeof(struct SymmetryBasisVector),
          SymmetryCompareBasisRepState);
  }
  for (merged_count = 1U; merged_count < local_count; merged_count++) {
    if (run->entries[merged_count].rep_state >=
        run->entries[merged_count + 1U].rep_state) {
      local_error = 1;
      break;
    }
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort",
      "local input keys are not strictly increasing", local_error);
  if (global_error != 0) return -1;

  all_counts = (uint64_t *)malloc(rank_array_bytes);
  all_gaps = (uint64_t *)malloc(rank_array_bytes);
  if (all_counts == NULL || all_gaps == NULL) local_error = 1;
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort", "rank metadata allocation",
      local_error);
  if (global_error != 0) goto fail;

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    if (MPI_Allgather(&local_count, 1, MPI_UINT64_T,
                      all_counts, 1, MPI_UINT64_T,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      report_distribution_operation_failure(
          rank, "sample sort", "MPI_Allgather of local entry counts");
      goto fail;
    }
  } else
#endif
  {
    all_counts[0] = local_count;
  }
  for (peer = 0; peer < nrank; peer++) {
    if (all_counts[peer] > 0U) nonempty_ranks++;
    if (checked_u64_add(global_count, all_counts[peer],
                        &global_count) != 0) {
      local_error = 1;
      break;
    }
  }
  if (nonempty_ranks > SYMMETRY_SAMPLE_COUNT_GLOBAL_CAP) {
    local_error = 1;
  } else if (nonempty_ranks > 0U) {
    sample_limit = SYMMETRY_SAMPLE_COUNT_DEFAULT;
    if (nonempty_ranks >
        SYMMETRY_SAMPLE_COUNT_GLOBAL_CAP / sample_limit) {
      sample_limit =
          SYMMETRY_SAMPLE_COUNT_GLOBAL_CAP / nonempty_ranks;
      if (sample_limit == 0U) sample_limit = 1U;
    }
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort",
      "global entry count or sampling limit", local_error);
  if (global_error != 0) goto fail;

  if (build_regular_samples(run, sample_limit, &local_samples,
                            &local_sample_count, &local_gap) != 0) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort", "regular sample construction",
      local_error);
  if (global_error != 0) goto fail;

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    if (MPI_Allgather(&local_gap, 1, MPI_UINT64_T,
                      all_gaps, 1, MPI_UINT64_T,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      report_distribution_operation_failure(
          rank, "sample sort", "MPI_Allgather of sample gaps");
      goto fail;
    }
  } else
#endif
  {
    all_gaps[0] = local_gap;
  }
  for (peer = 0; peer < nrank; peer++) {
    uint64_t peer_samples = all_counts[peer] < sample_limit
        ? all_counts[peer] : sample_limit;
    if (all_gaps[peer] > global_gap_max) global_gap_max = all_gaps[peer];
    if (checked_u64_add(global_gap_sum, all_gaps[peer],
                        &global_gap_sum) != 0 ||
        checked_u64_add(global_sample_count, peer_samples,
                        &global_sample_count) != 0) {
      local_error = 1;
      break;
    }
  }
  if (global_sample_count > SYMMETRY_SAMPLE_COUNT_GLOBAL_CAP) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort",
      "global sample count or gap arithmetic", local_error);
  if (global_error != 0) goto fail;

  sample_send_counts = (uint64_t *)calloc(
      (size_t)nrank, sizeof(*sample_send_counts));
  sample_send_displacements = (uint64_t *)calloc(
      (size_t)nrank, sizeof(*sample_send_displacements));
  if (sample_send_counts == NULL || sample_send_displacements == NULL) {
    local_error = 1;
  } else {
    uint64_t displacement = 0U;
    sample_send_counts[0] = local_sample_count;
    for (peer = 0; peer < nrank; peer++) {
      sample_send_displacements[peer] = displacement;
      if (checked_u64_add(displacement, sample_send_counts[peer],
                          &displacement) != 0) {
        local_error = 1;
        break;
      }
    }
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort", "sample gather schedule",
      local_error);
  if (global_error != 0) goto fail;

  sample_layout.nrank = nrank;
  sample_layout.count = local_sample_count;
  sample_layout.counts = sample_send_counts;
  sample_layout.displacements = sample_send_displacements;
  if (SymmetryMpiExchangeBasisVectors(
          local_samples, &sample_layout, rank, nrank, NULL,
          &sample_result, NULL) != 0) {
    (void)agree_distribution_failure(
        mpi_active, rank, "sample sort", "sample gather exchange", 1);
    goto fail;
  }
  if ((rank == 0 && sample_result.count != global_sample_count) ||
      (rank != 0 && sample_result.count != 0U)) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort", "sample gather receive count",
      local_error);
  if (global_error != 0) goto fail;

  splitter_count = nrank - 1;
  if (splitter_count > 0 &&
      checked_size_mul((size_t)splitter_count, sizeof(*splitters),
                       &splitter_bytes) != 0) {
    local_error = 1;
  }
  splitters = splitter_count > 0
      ? (unsigned long int *)malloc(splitter_bytes) : NULL;
  bucket_sample_counts = (uint64_t *)calloc(
      (size_t)nrank, sizeof(*bucket_sample_counts));
  if ((splitter_count > 0 && splitters == NULL) ||
      bucket_sample_counts == NULL) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort", "splitter allocation",
      local_error);
  if (global_error != 0) goto fail;

  if (rank == 0) {
    uint64_t index;
    if (sample_result.count > 1U) {
      qsort(sample_result.entries, (size_t)sample_result.count,
            sizeof(*sample_result.entries),
            SymmetryCompareBasisRepState);
    }
    if (sample_result.count == 0U) {
      for (peer = 0; peer < splitter_count; peer++) {
        splitters[peer] = ULONG_MAX;
      }
    } else {
      for (peer = 1; peer < nrank; peer++) {
        uint64_t sample_index =
            ((uint64_t)peer * sample_result.count) /
            (uint64_t)nrank;
        if (sample_index >= sample_result.count) {
          sample_index = sample_result.count - 1U;
        }
        splitters[peer - 1] =
            sample_result.entries[sample_index].rep_state;
      }
    }
    for (index = 0U; index < sample_result.count; index++) {
      unsigned long int bucket = symmetry_splitter_bucket(
          splitters, splitter_count,
          sample_result.entries[index].rep_state);
      bucket_sample_counts[bucket]++;
    }
  }

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    if ((splitter_count > 0 &&
         MPI_Bcast(splitters, splitter_count, MPI_UNSIGNED_LONG,
                   0, MPI_COMM_WORLD) != MPI_SUCCESS) ||
        MPI_Bcast(bucket_sample_counts, nrank, MPI_UINT64_T,
                  0, MPI_COMM_WORLD) != MPI_SUCCESS) {
      report_distribution_operation_failure(
          rank, "sample sort", "MPI_Bcast of splitters or bucket counts");
      goto fail;
    }
  }
#endif
  {
    uint64_t sample_count_check = 0U;
    for (peer = 0; peer < nrank; peer++) {
      if (checked_u64_add(sample_count_check,
                          bucket_sample_counts[peer],
                          &sample_count_check) != 0) {
        local_error = 1;
        break;
      }
    }
    if (sample_count_check != global_sample_count) local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort",
      "splitter bucket sample count", local_error);
  if (global_error != 0) goto fail;

  {
    uint64_t sample_bytes = 0U;
    uint64_t sample_receive_bytes = 0U;
    uint64_t input_elements = 0U;
    uint64_t input_bytes = 0U;
    uint64_t schedule_bytes = 0U;
    uint64_t splitter_and_bucket_bytes = 0U;
    uint64_t overhead_per_rank =
        UINT64_C(8) * (uint64_t)sizeof(uint64_t);
#ifdef MPI
    if (checked_u64_add(
            overhead_per_rank,
            UINT64_C(2) * (uint64_t)sizeof(MPI_Request),
            &overhead_per_rank) != 0) {
      local_error = 1;
    }
#endif
    if (checked_u64_add(local_count, 1U, &input_elements) != 0 ||
        checked_u64_mul(input_elements,
                        (uint64_t)sizeof(*run->entries),
                        &input_bytes) != 0 ||
        checked_u64_mul(local_sample_count,
                        (uint64_t)sizeof(*local_samples),
                        &sample_bytes) != 0 ||
        checked_u64_mul(sample_result.count,
                        (uint64_t)sizeof(*sample_result.entries),
                        &sample_receive_bytes) != 0 ||
        checked_u64_mul((uint64_t)nrank, overhead_per_rank,
                        &schedule_bytes) != 0 ||
        checked_u64_add((uint64_t)splitter_bytes,
                        (uint64_t)rank_array_bytes,
                        &splitter_and_bucket_bytes) != 0 ||
        checked_u64_add(input_bytes, sample_bytes,
                        &sample_peak) != 0 ||
        checked_u64_add(sample_peak, sample_receive_bytes,
                        &sample_peak) != 0 ||
        checked_u64_add(sample_peak, schedule_bytes,
                        &sample_peak) != 0 ||
        checked_u64_add(sample_peak, splitter_and_bucket_bytes,
                        &sample_peak) != 0 ||
        checked_u64_to_size(sample_peak, &sample_peak_size) != 0) {
      local_error = 1;
    }
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort",
      "sample-stage temporary memory arithmetic", local_error);
  if (global_error != 0) goto fail;
  local_error =
      sample_peak >
      (uint64_t)HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES;
  global_error = agree_distribution_limit_failure(
      mpi_active, rank, "sample sort",
      "sample-stage temporary memory limit exceeded; increase MPI ranks "
      "or HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES", local_error,
      sample_peak,
      (uint64_t)HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES);
  if (global_error != 0) goto fail;

  send_counts = (uint64_t *)calloc((size_t)nrank, sizeof(*send_counts));
  send_displacements = (uint64_t *)calloc(
      (size_t)nrank, sizeof(*send_displacements));
  preflight_recv_counts = (uint64_t *)calloc(
      (size_t)nrank, sizeof(*preflight_recv_counts));
  if (send_counts == NULL || send_displacements == NULL ||
      preflight_recv_counts == NULL) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort", "range schedule allocation",
      local_error);
  if (global_error != 0) goto fail;

  for (merged_count = 0U; merged_count < local_count; merged_count++) {
    unsigned long int bucket = symmetry_splitter_bucket(
        splitters, splitter_count,
        run->entries[merged_count + 1U].rep_state);
    send_counts[bucket]++;
  }
  receive_count = 0U;
  for (peer = 0; peer < nrank; peer++) {
    send_displacements[peer] = receive_count;
    if (checked_u64_add(receive_count, send_counts[peer],
                        &receive_count) != 0) {
      local_error = 1;
      break;
    }
  }
  if (receive_count != local_count) local_error = 1;
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort", "range send schedule arithmetic",
      local_error);
  if (global_error != 0) goto fail;
  receive_count = 0U;

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    if (MPI_Alltoall(send_counts, 1, MPI_UINT64_T,
                     preflight_recv_counts, 1, MPI_UINT64_T,
                     MPI_COMM_WORLD) != MPI_SUCCESS) {
      report_distribution_operation_failure(
          rank, "sample sort", "MPI_Alltoall range preflight");
      goto fail;
    }
  } else
#endif
  {
    preflight_recv_counts[0] = send_counts[0];
  }
  for (peer = 0; peer < nrank; peer++) {
    if (checked_u64_add(receive_count, preflight_recv_counts[peer],
                        &receive_count) != 0) {
      local_error = 1;
      break;
    }
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort",
      "range receive count arithmetic", local_error);
  if (global_error != 0) goto fail;
  {
    uint64_t sampled_gap_bound;
    if (checked_u64_mul(bucket_sample_counts[rank],
                        global_gap_max,
                        &sampled_gap_bound) != 0 ||
        checked_u64_add(sampled_gap_bound, global_gap_sum,
                        &bucket_upper_bound) != 0) {
      local_error = 1;
    }
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort",
      "sample-derived bucket upper-bound arithmetic", local_error);
  if (global_error != 0) goto fail;
  local_error = receive_count > bucket_upper_bound;
  global_error = agree_distribution_limit_failure(
      mpi_active, rank, "sample sort",
      "range receive count exceeds the sample-derived bucket bound",
      local_error, receive_count, bucket_upper_bound);
  if (global_error != 0) goto fail;
  if (calculate_range_temporary_peak(
          local_count, receive_count, nrank, &temporary_peak) != 0 ||
      temporary_peak > (uint64_t)SIZE_MAX ||
      receive_count >= (uint64_t)ULONG_MAX) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort",
      "range-stage size or temporary memory arithmetic", local_error);
  if (global_error != 0) goto fail;
  local_error =
      temporary_peak >
      (uint64_t)HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES;
  global_error = agree_distribution_limit_failure(
      mpi_active, rank, "sample sort",
      "range-stage temporary memory limit exceeded; increase MPI ranks "
      "or HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES", local_error,
      temporary_peak,
      (uint64_t)HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES);
  if (global_error != 0) goto fail;

  next_stats.local_survivor_entries = local_count;
  next_stats.nonempty_rank_count = nonempty_ranks;
  next_stats.samples_per_nonempty_rank = sample_limit;
  next_stats.local_sample_entries = local_sample_count;
  next_stats.global_sample_entries = global_sample_count;
  next_stats.local_sample_gap = local_gap;
  next_stats.global_sample_gap_max = global_gap_max;
  next_stats.global_sample_gap_sum = global_gap_sum;
  next_stats.bucket_sample_entries = bucket_sample_counts[rank];
  next_stats.bucket_entry_upper_bound = bucket_upper_bound;
  next_stats.global_entries = global_count;
  next_stats.distribution_memory_byte_limit =
      (uint64_t)HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES;
  next_stats.splitter_digest =
      digest_splitters(splitters, splitter_count);
  next_stats.sample_temporary_peak_bytes = sample_peak_size;

  FreeSymmetryMpiExchangeResult(&sample_result);
  free(local_samples);
  local_samples = NULL;
  free(sample_send_counts);
  sample_send_counts = NULL;
  free(sample_send_displacements);
  sample_send_displacements = NULL;
  free(splitters);
  splitters = NULL;
  free(bucket_sample_counts);
  bucket_sample_counts = NULL;

  range_layout.nrank = nrank;
  range_layout.count = local_count;
  range_layout.counts = send_counts;
  range_layout.displacements = send_displacements;
  if (SymmetryMpiExchangeBasisVectors(
          local_count > 0U ? run->entries + 1 : run->entries,
          &range_layout, rank, nrank, NULL,
          &range_result, NULL) != 0) {
    (void)agree_distribution_failure(
        mpi_active, rank, "sample sort", "range exchange", 1);
    goto fail;
  }
  if (range_result.count != receive_count) local_error = 1;
  for (peer = 0; peer < nrank; peer++) {
    if (range_result.counts[peer] != preflight_recv_counts[peer]) {
      local_error = 1;
      break;
    }
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort",
      "range exchange receive schedule", local_error);
  if (global_error != 0) goto fail;

  free(send_counts);
  send_counts = NULL;
  free(send_displacements);
  send_displacements = NULL;
  free(preflight_recv_counts);
  preflight_recv_counts = NULL;

  {
    uint64_t merged_elements;
    uint64_t merged_byte_count;
    if (checked_u64_add(receive_count, 1U, &merged_elements) != 0 ||
        checked_u64_mul(merged_elements,
                        (uint64_t)sizeof(*merged_entries),
                        &merged_byte_count) != 0 ||
        checked_u64_to_size(merged_byte_count, &merged_bytes) != 0) {
      local_error = 1;
    }
  }
  if (local_error != 0 ||
      checked_size_mul((size_t)nrank, sizeof(*merge_heap),
                       &heap_bytes) != 0) {
    local_error = 1;
  } else {
    merged_entries =
        (struct SymmetryBasisVector *)calloc(1U, merged_bytes);
    merge_heap =
        (struct SymmetryMergeNode *)malloc(heap_bytes);
    all_last_keys = (uint64_t *)malloc(rank_array_bytes);
    if (merged_entries == NULL || merge_heap == NULL ||
        all_last_keys == NULL) {
      local_error = 1;
    }
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort", "range merge allocation",
      local_error);
  if (global_error != 0) goto fail;

  for (peer = 0; peer < nrank; peer++) {
    if (range_result.counts[peer] > 0U) {
      struct SymmetryMergeNode node;
      node.source = peer;
      node.offset = 0U;
      merge_heap_push(merge_heap, &merge_heap_count, node, &range_result);
    }
  }
  merged_count = 0U;
  while (merge_heap_count > 0) {
    struct SymmetryMergeNode node =
        merge_heap_pop(merge_heap, &merge_heap_count, &range_result);
    const struct SymmetryBasisVector *entry =
        &range_result.entries[
            range_result.displacements[node.source] + node.offset];
    if (merged_count > 0U &&
        merged_entries[merged_count].rep_state >= entry->rep_state) {
      local_error = 1;
      break;
    }
    merged_count++;
    memcpy(&merged_entries[merged_count], entry, sizeof(*merged_entries));
    node.offset++;
    if (node.offset < range_result.counts[node.source]) {
      merge_heap_push(merge_heap, &merge_heap_count, node, &range_result);
    }
  }
  if (merged_count != receive_count) local_error = 1;
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort",
      "merged keys are not globally strictly increasing", local_error);
  if (global_error != 0) goto fail;

  {
    unsigned long int first_key_ulong = receive_count > 0U
        ? merged_entries[1].rep_state : 0UL;
    unsigned long int last_key_ulong = receive_count > 0U
        ? merged_entries[receive_count].rep_state : 0UL;
    uint64_t first_key = (uint64_t)first_key_ulong;
    uint64_t last_key = (uint64_t)last_key_ulong;
    if ((unsigned long int)first_key != first_key_ulong ||
        (unsigned long int)last_key != last_key_ulong) {
      local_error = 1;
    }
    global_error = agree_distribution_failure(
        mpi_active, rank, "sample sort",
        "range boundary key representation", local_error);
    if (global_error != 0) goto fail;
#ifdef MPI
    if (mpi_active != FALSE && nrank > 1) {
      if (MPI_Allgather(&receive_count, 1, MPI_UINT64_T,
                        all_counts, 1, MPI_UINT64_T,
                        MPI_COMM_WORLD) != MPI_SUCCESS ||
          MPI_Allgather(&first_key, 1, MPI_UINT64_T,
                        all_gaps, 1, MPI_UINT64_T,
                        MPI_COMM_WORLD) != MPI_SUCCESS ||
          MPI_Allgather(&last_key, 1, MPI_UINT64_T,
                        all_last_keys, 1, MPI_UINT64_T,
                        MPI_COMM_WORLD) != MPI_SUCCESS) {
        report_distribution_operation_failure(
            rank, "sample sort", "MPI_Allgather of range boundaries");
        goto fail;
      }
    } else
#endif
    {
      all_counts[0] = receive_count;
      all_gaps[0] = first_key;
      all_last_keys[0] = last_key;
    }
  }
  {
    int have_previous = FALSE;
    uint64_t previous_last = 0U;
    for (peer = 0; peer < nrank; peer++) {
      if (checked_u64_add(range_global_count, all_counts[peer],
                          &range_global_count) != 0) {
        local_error = 1;
        break;
      }
      if (all_counts[peer] > range_max_count) {
        range_max_count = all_counts[peer];
      }
      if (all_counts[peer] == 0U) continue;
      if (have_previous != FALSE &&
          previous_last >= all_gaps[peer]) {
        local_error = 1;
        break;
      }
      previous_last = all_last_keys[peer];
      have_previous = TRUE;
    }
    if (range_global_count != global_count) local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "sample sort",
      "global range count or adjacent-rank boundary ordering",
      local_error);
  if (global_error != 0) goto fail;

  next_stats.range_entries = receive_count;
  next_stats.range_max_entries = range_max_count;
  next_stats.sample_send_entries = local_count;
  next_stats.sample_recv_entries = receive_count;
  next_stats.range_max_over_mean = global_count > 0U
      ? ((double)range_max_count * (double)nrank) /
          (double)global_count
      : 0.0;
  next_stats.range_digest =
      digest_basis_run(merged_entries + 1, receive_count);
  next_stats.sample_temporary_peak_bytes = sample_peak_size >
      (size_t)temporary_peak
      ? sample_peak_size : (size_t)temporary_peak;

  free(run->entries);
  run->entries = merged_entries;
  run->count = (unsigned long int)receive_count;
  run->capacity = run->count + 1UL;
  merged_entries = NULL;
  if (stats != NULL) *stats = next_stats;

  FreeSymmetryMpiExchangeResult(&range_result);
  free(merge_heap);
  free(all_counts);
  free(all_gaps);
  free(all_last_keys);
  return 0;

fail:
  FreeSymmetryMpiExchangeResult(&sample_result);
  FreeSymmetryMpiExchangeResult(&range_result);
  free(local_samples);
  free(merged_entries);
  free(merge_heap);
  free(all_counts);
  free(all_gaps);
  free(all_last_keys);
  free(sample_send_counts);
  free(sample_send_displacements);
  free(bucket_sample_counts);
  free(send_counts);
  free(send_displacements);
  free(preflight_recv_counts);
  free(splitters);
  return -1;
}

static int calculate_rebalance_temporary_peak(
    uint64_t input_count,
    uint64_t receive_count,
    int nrank,
    size_t *peak_bytes)
{
  uint64_t input_elements;
  uint64_t input_bytes;
  uint64_t receive_bytes;
  uint64_t output_elements;
  uint64_t output_bytes;
  uint64_t rank_elements;
  uint64_t schedule_bytes;
  uint64_t request_bytes = 0U;
  uint64_t peak;
  if (peak_bytes == NULL || nrank < 1 ||
      checked_u64_add(input_count, 1U, &input_elements) != 0 ||
      checked_u64_mul(input_elements,
                      (uint64_t)sizeof(struct SymmetryBasisVector),
                      &input_bytes) != 0 ||
      checked_u64_mul(receive_count,
                      (uint64_t)sizeof(struct SymmetryBasisVector),
                      &receive_bytes) != 0 ||
      checked_u64_add(receive_count, 1U, &output_elements) != 0 ||
      checked_u64_mul(output_elements,
                      (uint64_t)sizeof(struct SymmetryBasisVector),
                      &output_bytes) != 0 ||
      checked_u64_mul((uint64_t)nrank, UINT64_C(10),
                      &rank_elements) != 0 ||
      checked_u64_add(rank_elements, 4U, &rank_elements) != 0 ||
      checked_u64_mul(rank_elements, (uint64_t)sizeof(uint64_t),
                      &schedule_bytes) != 0) {
    return -1;
  }
#ifdef MPI
  if (checked_u64_mul(
          (uint64_t)nrank,
          UINT64_C(2) * (uint64_t)sizeof(MPI_Request),
          &request_bytes) != 0) {
    return -1;
  }
#endif
  if (checked_u64_add(input_bytes, receive_bytes, &peak) != 0 ||
      checked_u64_add(peak, output_bytes, &peak) != 0 ||
      checked_u64_add(peak, schedule_bytes, &peak) != 0 ||
      checked_u64_add(peak, request_bytes, &peak) != 0 ||
      checked_u64_to_size(peak, peak_bytes) != 0) {
    return -1;
  }
  return 0;
}

int SymmetryExactRebalanceBasisRun(
    struct SymmetryBasisRun *run,
    int rank,
    int nrank,
    struct SymmetryBasisOwnership *ownership,
    struct SymmetryBasisDistributionStats *stats)
{
  struct SymmetryBasisOwnership next_ownership;
  struct SymmetryMpiExchangeLayout layout;
  struct SymmetryMpiExchangeResult result;
  struct SymmetryMpiExchangeStats exchange_stats;
  struct SymmetryBasisVector *next_entries = NULL;
  uint64_t *all_counts = NULL;
  uint64_t *source_offsets = NULL;
  uint64_t *send_counts = NULL;
  uint64_t *send_displacements = NULL;
  uint64_t *preflight_recv_counts = NULL;
  unsigned long int *all_first_keys = NULL;
  unsigned long int *all_last_keys = NULL;
  uint64_t local_count = 0U;
  uint64_t global_count = 0U;
  uint64_t source_begin;
  uint64_t source_end;
  uint64_t target_begin;
  uint64_t target_end;
  uint64_t send_total = 0U;
  uint64_t receive_total = 0U;
  unsigned long int global_dim = 0UL;
  unsigned long int local_offset = 0UL;
  unsigned long int local_dim = 0UL;
  unsigned long int first_key = 0UL;
  unsigned long int last_key = 0UL;
  size_t rank_bytes = 0U;
  size_t offset_bytes = 0U;
  size_t key_rank_bytes = 0U;
  size_t ownership_bytes = 0U;
  size_t output_elements = 0U;
  size_t output_bytes = 0U;
  size_t temporary_peak = 0U;
  uint64_t entry_index;
  int mpi_active;
  int local_error = 0;
  int global_error;
  int peer;
#ifdef MPI
  int comm_rank = 0;
  int comm_size = 1;
#endif

  memset(&next_ownership, 0, sizeof(next_ownership));
  memset(&result, 0, sizeof(result));
  memset(&exchange_stats, 0, sizeof(exchange_stats));
  mpi_active = SymmetryMpiCollectivesActive();

  if (rank < 0 || nrank < 1 || rank >= nrank) {
    local_error = 1;
  }
#ifdef MPI
  if (mpi_active != FALSE) {
    if (MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank) != MPI_SUCCESS ||
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size) != MPI_SUCCESS ||
        comm_rank != rank || comm_size != nrank) {
      local_error = 1;
    }
  } else if (rank != 0 || nrank != 1) {
    local_error = 1;
  }
#else
  if (rank != 0 || nrank != 1) local_error = 1;
#endif
  if (nrank > 0 &&
      (checked_size_mul((size_t)nrank, sizeof(uint64_t),
                        &rank_bytes) != 0 ||
       checked_size_mul((size_t)nrank, sizeof(unsigned long int),
                        &key_rank_bytes) != 0 ||
       (size_t)nrank == SIZE_MAX ||
       checked_size_mul((size_t)nrank + 1U, sizeof(uint64_t),
                        &offset_bytes) != 0 ||
       checked_size_mul((size_t)nrank + 1U,
                        sizeof(unsigned long int),
                        &ownership_bytes) != 0)) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance",
      "invalid rank or communicator", local_error);
  if (global_error != 0) return -1;
  if (symmetry_basis_run_is_valid(run) != TRUE ||
      checked_ulong_to_u64(run != NULL ? run->count : 0UL,
                           &local_count) != 0 ||
      HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES == 0U) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance",
      "invalid input run or distribution memory configuration",
      local_error);
  if (global_error != 0) return -1;
  if (ownership == NULL || ownership->dim != 0UL ||
      ownership->local_offset != 0UL || ownership->local_dim != 0UL ||
      ownership->rank_offsets != NULL) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance",
      "ownership must be empty on entry", local_error);
  if (global_error != 0) return -1;

  for (entry_index = 1U; entry_index < local_count; entry_index++) {
    if (run->entries[entry_index].rep_state >=
        run->entries[entry_index + 1U].rep_state) {
      local_error = 1;
      break;
    }
  }
  first_key = local_count > 0U ? run->entries[1].rep_state : 0UL;
  last_key = local_count > 0U
      ? run->entries[(unsigned long int)local_count].rep_state : 0UL;
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance",
      "local input keys are not strictly increasing", local_error);
  if (global_error != 0) return -1;

  all_counts = (uint64_t *)malloc(rank_bytes);
  source_offsets = (uint64_t *)malloc(offset_bytes);
  all_first_keys = (unsigned long int *)malloc(key_rank_bytes);
  all_last_keys = (unsigned long int *)malloc(key_rank_bytes);
  if (all_counts == NULL || source_offsets == NULL ||
      all_first_keys == NULL || all_last_keys == NULL) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance", "rank metadata allocation",
      local_error);
  if (global_error != 0) goto fail_rebalance;

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    if (MPI_Allgather(&local_count, 1, MPI_UINT64_T,
                      all_counts, 1, MPI_UINT64_T,
                      MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allgather(&first_key, 1, MPI_UNSIGNED_LONG,
                      all_first_keys, 1, MPI_UNSIGNED_LONG,
                      MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allgather(&last_key, 1, MPI_UNSIGNED_LONG,
                      all_last_keys, 1, MPI_UNSIGNED_LONG,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      report_distribution_operation_failure(
          rank, "exact rebalance",
          "MPI_Allgather of input counts or boundaries");
      goto fail_rebalance;
    }
  } else
#endif
  {
    all_counts[0] = local_count;
    all_first_keys[0] = first_key;
    all_last_keys[0] = last_key;
  }

  source_offsets[0] = 0U;
  {
    int have_previous = FALSE;
    unsigned long int previous_last = 0UL;
    for (peer = 0; peer < nrank; peer++) {
      if (checked_u64_add(source_offsets[peer], all_counts[peer],
                          &source_offsets[peer + 1]) != 0) {
        local_error = 1;
        break;
      }
      if (all_counts[peer] == 0U) continue;
      if (have_previous != FALSE &&
          previous_last >= all_first_keys[peer]) {
        local_error = 1;
        break;
      }
      previous_last = all_last_keys[peer];
      have_previous = TRUE;
    }
  }
  if (local_error == 0) global_count = source_offsets[nrank];
  if (local_error == 0 &&
      (global_count > (uint64_t)ULONG_MAX ||
       (uint64_t)(unsigned long int)global_count != global_count)) {
    local_error = 1;
  } else if (local_error == 0) {
    global_dim = (unsigned long int)global_count;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance",
      "global input count or adjacent-rank boundary ordering",
      local_error);
  if (global_error != 0) goto fail_rebalance;

  next_ownership.rank_offsets =
      (unsigned long int *)malloc(ownership_bytes);
  if (next_ownership.rank_offsets == NULL) local_error = 1;
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance", "ownership offset allocation",
      local_error);
  if (global_error != 0) goto fail_rebalance;

  for (peer = 0; peer < nrank; peer++) {
    unsigned long int offset;
    unsigned long int count;
    if (SymmetryBlockRange(global_dim, peer, nrank,
                           &offset, &count) != 0 ||
        offset > global_dim || count > global_dim - offset) {
      local_error = 1;
      break;
    }
    next_ownership.rank_offsets[peer] = offset;
    if (peer + 1 == nrank) {
      next_ownership.rank_offsets[peer + 1] = offset + count;
    }
  }
  if (local_error == 0 &&
      (next_ownership.rank_offsets[0] != 0UL ||
       next_ownership.rank_offsets[nrank] != global_dim ||
       SymmetryBlockRange(global_dim, rank, nrank,
                          &local_offset, &local_dim) != 0 ||
       local_dim == ULONG_MAX)) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance", "exact block range construction",
      local_error);
  if (global_error != 0) goto fail_rebalance;

  send_counts = (uint64_t *)calloc((size_t)nrank, sizeof(*send_counts));
  send_displacements = (uint64_t *)calloc(
      (size_t)nrank, sizeof(*send_displacements));
  preflight_recv_counts = (uint64_t *)calloc(
      (size_t)nrank, sizeof(*preflight_recv_counts));
  if (send_counts == NULL || send_displacements == NULL ||
      preflight_recv_counts == NULL) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance", "exchange schedule allocation",
      local_error);
  if (global_error != 0) goto fail_rebalance;

  source_begin = source_offsets[rank];
  source_end = source_offsets[rank + 1];
  for (peer = 0; peer < nrank; peer++) {
    uint64_t block_begin =
        (uint64_t)next_ownership.rank_offsets[peer];
    uint64_t block_end =
        (uint64_t)next_ownership.rank_offsets[peer + 1];
    uint64_t intersection_begin =
        source_begin > block_begin ? source_begin : block_begin;
    uint64_t intersection_end =
        source_end < block_end ? source_end : block_end;
    uint64_t intersection_count = intersection_end > intersection_begin
        ? intersection_end - intersection_begin : 0U;
    send_displacements[peer] = send_total;
    send_counts[peer] = intersection_count;
    if ((intersection_count > 0U &&
         intersection_begin - source_begin != send_total) ||
        checked_u64_add(send_total, intersection_count,
                        &send_total) != 0) {
      local_error = 1;
      break;
    }
  }
  if (send_total != local_count) local_error = 1;
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance",
      "destination send slices do not form a contiguous source range",
      local_error);
  if (global_error != 0) goto fail_rebalance;

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    if (MPI_Alltoall(send_counts, 1, MPI_UINT64_T,
                     preflight_recv_counts, 1, MPI_UINT64_T,
                     MPI_COMM_WORLD) != MPI_SUCCESS) {
      report_distribution_operation_failure(
          rank, "exact rebalance", "MPI_Alltoall exchange preflight");
      goto fail_rebalance;
    }
  } else
#endif
  {
    preflight_recv_counts[0] = send_counts[0];
  }

  target_begin = (uint64_t)local_offset;
  target_end = target_begin + (uint64_t)local_dim;
  receive_total = 0U;
  for (peer = 0; peer < nrank; peer++) {
    uint64_t peer_begin = source_offsets[peer];
    uint64_t peer_end = source_offsets[peer + 1];
    uint64_t intersection_begin =
        peer_begin > target_begin ? peer_begin : target_begin;
    uint64_t intersection_end =
        peer_end < target_end ? peer_end : target_end;
    uint64_t expected_count = intersection_end > intersection_begin
        ? intersection_end - intersection_begin : 0U;
    if (preflight_recv_counts[peer] != expected_count ||
        checked_u64_add(receive_total, preflight_recv_counts[peer],
                        &receive_total) != 0) {
      local_error = 1;
      break;
    }
  }
  if (receive_total != (uint64_t)local_dim) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance",
      "source receive slices do not match the exact local block",
      local_error);
  if (global_error != 0) goto fail_rebalance;
  if (calculate_rebalance_temporary_peak(
          local_count, receive_total, nrank, &temporary_peak) != 0) {
    local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance",
      "temporary memory arithmetic", local_error);
  if (global_error != 0) goto fail_rebalance;
  local_error =
      (uint64_t)temporary_peak >
      (uint64_t)HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES;
  global_error = agree_distribution_limit_failure(
      mpi_active, rank, "exact rebalance",
      "temporary memory limit exceeded; increase MPI ranks "
      "or HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES", local_error,
      (uint64_t)temporary_peak,
      (uint64_t)HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES);
  if (global_error != 0) goto fail_rebalance;

  layout.nrank = nrank;
  layout.count = local_count;
  layout.counts = send_counts;
  layout.displacements = send_displacements;
  if (SymmetryMpiExchangeBasisVectors(
          local_count > 0U ? run->entries + 1 : run->entries,
          &layout, rank, nrank, NULL, &result, &exchange_stats) != 0) {
    (void)agree_distribution_failure(
        mpi_active, rank, "exact rebalance", "basis exchange", 1);
    goto fail_rebalance;
  }
  if (result.count != receive_total ||
      exchange_stats.send_entries != local_count ||
      exchange_stats.recv_entries != receive_total) {
    local_error = 1;
  }
  send_total = 0U;
  for (peer = 0; peer < nrank; peer++) {
    if (result.counts[peer] != preflight_recv_counts[peer] ||
        result.displacements[peer] != send_total ||
        checked_u64_add(send_total, result.counts[peer],
                        &send_total) != 0) {
      local_error = 1;
      break;
    }
  }
  if (send_total != receive_total) local_error = 1;
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance",
      "receive count or source-rank segment layout", local_error);
  if (global_error != 0) goto fail_rebalance;

  {
    uint64_t output_element_count;
    if (checked_u64_add(receive_total, 1U, &output_element_count) != 0 ||
        checked_u64_to_size(output_element_count, &output_elements) != 0) {
      local_error = 1;
    }
  }
  if (local_error != 0 ||
      checked_size_mul(output_elements, sizeof(*next_entries),
                       &output_bytes) != 0) {
    local_error = 1;
  } else {
    next_entries = (struct SymmetryBasisVector *)calloc(1U, output_bytes);
    if (next_entries == NULL) local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance", "output allocation",
      local_error);
  if (global_error != 0) goto fail_rebalance;

  if (receive_total > 0U) {
    memcpy(next_entries + 1, result.entries,
           (size_t)receive_total * sizeof(*next_entries));
  }
  for (send_total = 1U; send_total < receive_total; send_total++) {
    if (next_entries[send_total].rep_state >=
        next_entries[send_total + 1U].rep_state) {
      local_error = 1;
      break;
    }
  }
  first_key = receive_total > 0U ? next_entries[1].rep_state : 0UL;
  last_key = receive_total > 0U
      ? next_entries[(unsigned long int)receive_total].rep_state : 0UL;
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance",
      "received keys are not strictly increasing", local_error);
  if (global_error != 0) goto fail_rebalance;

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    if (MPI_Allgather(&receive_total, 1, MPI_UINT64_T,
                      all_counts, 1, MPI_UINT64_T,
                      MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allgather(&first_key, 1, MPI_UNSIGNED_LONG,
                      all_first_keys, 1, MPI_UNSIGNED_LONG,
                      MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allgather(&last_key, 1, MPI_UNSIGNED_LONG,
                      all_last_keys, 1, MPI_UNSIGNED_LONG,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      report_distribution_operation_failure(
          rank, "exact rebalance",
          "MPI_Allgather of output counts or boundaries");
      goto fail_rebalance;
    }
  } else
#endif
  {
    all_counts[0] = receive_total;
    all_first_keys[0] = first_key;
    all_last_keys[0] = last_key;
  }
  {
    uint64_t verified_global_count = 0U;
    int have_previous = FALSE;
    unsigned long int previous_last = 0UL;
    for (peer = 0; peer < nrank; peer++) {
      uint64_t expected_count =
          (uint64_t)next_ownership.rank_offsets[peer + 1] -
          (uint64_t)next_ownership.rank_offsets[peer];
      if (all_counts[peer] != expected_count ||
          checked_u64_add(verified_global_count, all_counts[peer],
                          &verified_global_count) != 0) {
        local_error = 1;
        break;
      }
      if (all_counts[peer] == 0U) continue;
      if (have_previous != FALSE &&
          previous_last >= all_first_keys[peer]) {
        local_error = 1;
        break;
      }
      previous_last = all_last_keys[peer];
      have_previous = TRUE;
    }
    if (verified_global_count != global_count) local_error = 1;
  }
  global_error = agree_distribution_failure(
      mpi_active, rank, "exact rebalance",
      "final exact counts or adjacent-rank boundary ordering",
      local_error);
  if (global_error != 0) goto fail_rebalance;

  next_ownership.dim = global_dim;
  next_ownership.local_offset = local_offset;
  next_ownership.local_dim = local_dim;
  free(run->entries);
  run->entries = next_entries;
  run->count = local_dim;
  run->capacity = local_dim + 1UL;
  next_entries = NULL;
  *ownership = next_ownership;
  next_ownership.rank_offsets = NULL;
  if (stats != NULL) {
    stats->rebalance_send_entries = local_count;
    stats->rebalance_recv_entries = receive_total;
    stats->rebalance_temporary_peak_bytes = temporary_peak;
  }

  FreeSymmetryMpiExchangeResult(&result);
  free(all_counts);
  free(source_offsets);
  free(send_counts);
  free(send_displacements);
  free(preflight_recv_counts);
  free(all_first_keys);
  free(all_last_keys);
  return 0;

fail_rebalance:
  FreeSymmetryMpiExchangeResult(&result);
  free(next_entries);
  free(next_ownership.rank_offsets);
  free(all_counts);
  free(source_offsets);
  free(send_counts);
  free(send_displacements);
  free(preflight_recv_counts);
  free(all_first_keys);
  free(all_last_keys);
  return -1;
}
