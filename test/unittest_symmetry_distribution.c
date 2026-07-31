#include <complex.h>
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

static int test_rank = 0;
static int test_nrank = 1;
static uint64_t test_message_entry_limit = 0U;
static uint64_t test_request_message_entry_limit = 0U;
static uint64_t test_response_message_entry_limit = 0U;
#ifdef MPI
static int test_mpi_active = FALSE;
#endif

static void fail_test(const char *label)
{
  fprintf(stderr, "rank %d: %s\n", test_rank, label);
#ifdef MPI
  if (test_mpi_active != FALSE) MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  exit(1);
}

static void require_true(int condition, const char *label)
{
  if (condition == 0) fail_test(label);
}

static int field_bits_equal(const struct SymmetryBasisVector *lhs,
                            const struct SymmetryBasisVector *rhs)
{
  return memcmp(&lhs->rep_state, &rhs->rep_state,
                sizeof(lhs->rep_state)) == 0 &&
         memcmp(&lhs->orbit_size, &rhs->orbit_size,
                sizeof(lhs->orbit_size)) == 0 &&
         memcmp(&lhs->stabilizer_size, &rhs->stabilizer_size,
                sizeof(lhs->stabilizer_size)) == 0 &&
         memcmp(&lhs->norm, &rhs->norm, sizeof(lhs->norm)) == 0 &&
         memcmp(&lhs->stabilizer_character_sum,
                &rhs->stabilizer_character_sum,
                sizeof(lhs->stabilizer_character_sum)) == 0 &&
         memcmp(&lhs->diagonal, &rhs->diagonal,
                sizeof(lhs->diagonal)) == 0;
}

static struct SymmetryBasisVector make_entry(int source,
                                             int destination,
                                             uint64_t index)
{
  struct SymmetryBasisVector entry;
  unsigned long int identifier =
      (unsigned long int)(1U + (uint64_t)source * 100000U +
                          (uint64_t)destination * 1000U + index);
  memset(&entry, 0, sizeof(entry));
  entry.rep_state = identifier;
  entry.orbit_size = (unsigned int)(2U + (identifier % 17U));
  entry.stabilizer_size = (unsigned int)(1U + (identifier % 7U));
  entry.norm = 0.125 + (double)identifier * 0.5;
  entry.stabilizer_character_sum =
      (0.25 + (double)identifier) +
      (0.75 - (double)identifier * 0.125) * I;
  entry.diagonal = -0.5 - (double)identifier * 0.25;
  return entry;
}

static uint64_t asymmetric_count(int source, int destination)
{
  return (uint64_t)((source * 3 + destination * 2 + 1) % 4);
}

static void build_layout(
    uint64_t (*count_function)(int, int),
    struct SymmetryBasisVector **entries,
    uint64_t **counts,
    uint64_t **displacements,
    struct SymmetryMpiExchangeLayout *layout)
{
  uint64_t total = 0U;
  int destination;
  *counts = (uint64_t *)calloc((size_t)test_nrank, sizeof(**counts));
  *displacements =
      (uint64_t *)calloc((size_t)test_nrank, sizeof(**displacements));
  require_true(*counts != NULL && *displacements != NULL,
               "layout allocation failed");
  for (destination = 0; destination < test_nrank; destination++) {
    (*counts)[destination] =
        count_function(test_rank, destination);
    (*displacements)[destination] = total;
    total += (*counts)[destination];
  }
  *entries = total > 0U
      ? (struct SymmetryBasisVector *)malloc(
            (size_t)total * sizeof(**entries))
      : NULL;
  require_true(total == 0U || *entries != NULL,
               "send entry allocation failed");
  for (destination = 0; destination < test_nrank; destination++) {
    uint64_t index;
    for (index = 0U; index < (*counts)[destination]; index++) {
      (*entries)[(*displacements)[destination] + index] =
          make_entry(test_rank, destination, index);
    }
  }
  layout->nrank = test_nrank;
  layout->count = total;
  layout->counts = *counts;
  layout->displacements = *displacements;
}

static void free_layout(struct SymmetryBasisVector *entries,
                        uint64_t *counts,
                        uint64_t *displacements)
{
  free(entries);
  free(counts);
  free(displacements);
}

static void validate_result(
    const struct SymmetryMpiExchangeResult *result,
    uint64_t (*count_function)(int, int))
{
  uint64_t total = 0U;
  int source;
  require_true(result->nrank == test_nrank,
               "receive rank count mismatch");
  for (source = 0; source < test_nrank; source++) {
    uint64_t expected_count = count_function(source, test_rank);
    uint64_t index;
    require_true(result->counts[source] == expected_count,
                 "receive peer count mismatch");
    require_true(result->displacements[source] == total,
                 "receive source-rank segment order mismatch");
    for (index = 0U; index < expected_count; index++) {
      struct SymmetryBasisVector expected =
          make_entry(source, test_rank, index);
      require_true(
          field_bits_equal(
              &result->entries[result->displacements[source] + index],
              &expected),
          "received payload field mismatch");
    }
    total += expected_count;
  }
  require_true(result->count == total, "receive total mismatch");
}

static void validate_asymmetric_result(
    const struct SymmetryMpiExchangeResult *result)
{
  validate_result(result, asymmetric_count);
}

static void validate_results_equal(
    const struct SymmetryMpiExchangeResult *lhs,
    const struct SymmetryMpiExchangeResult *rhs)
{
  uint64_t index;
  int peer;
  require_true(lhs->nrank == rhs->nrank && lhs->count == rhs->count,
               "fast/chunk result shape mismatch");
  for (peer = 0; peer < lhs->nrank; peer++) {
    require_true(lhs->counts[peer] == rhs->counts[peer] &&
                     lhs->displacements[peer] == rhs->displacements[peer],
                 "fast/chunk source segment mismatch");
  }
  for (index = 0U; index < lhs->count; index++) {
    require_true(field_bits_equal(&lhs->entries[index],
                                  &rhs->entries[index]),
                 "fast/chunk payload mismatch");
  }
}

static void assert_asymmetric_fast_and_chunked(void)
{
  struct SymmetryBasisVector *entries = NULL;
  uint64_t *counts = NULL;
  uint64_t *displacements = NULL;
  struct SymmetryMpiExchangeLayout layout;
  struct SymmetryMpiExchangeOptions fast_options;
  struct SymmetryMpiExchangeOptions chunk_options;
  struct SymmetryMpiExchangeResult fast_result;
  struct SymmetryMpiExchangeResult chunk_result;
  struct SymmetryMpiExchangeStats fast_stats;
  struct SymmetryMpiExchangeStats chunk_stats;
  uint64_t expected_send_messages = 0U;
  uint64_t expected_recv_messages = 0U;
  int peer;

  memset(&fast_result, 0, sizeof(fast_result));
  memset(&chunk_result, 0, sizeof(chunk_result));
  fast_options.chunk_limit = 0U;
  fast_options.force_chunked = FALSE;
  chunk_options.chunk_limit = 2U;
  chunk_options.force_chunked = TRUE;
  build_layout(asymmetric_count, &entries, &counts, &displacements, &layout);

  require_true(
      SymmetryMpiExchangeBasisVectors(
          entries, &layout, test_rank, test_nrank, &fast_options,
          &fast_result, &fast_stats) == 0,
      "fast exchange failed");
  validate_asymmetric_result(&fast_result);
  require_true(fast_stats.send_entries == layout.count &&
                   fast_stats.recv_entries == fast_result.count,
               "fast exchange stats mismatch");
  require_true(fast_stats.used_chunked == FALSE,
               "fast exchange unexpectedly used chunking");

  require_true(
      SymmetryMpiExchangeBasisVectors(
          entries, &layout, test_rank, test_nrank, &chunk_options,
          &chunk_result, &chunk_stats) == 0,
      "forced chunk exchange failed");
  validate_asymmetric_result(&chunk_result);
  validate_results_equal(&fast_result, &chunk_result);
  if (test_nrank > 1) {
    for (peer = 0; peer < test_nrank; peer++) {
      if (peer == test_rank) continue;
      expected_send_messages += (counts[peer] + 1U) / 2U;
      expected_recv_messages +=
          (asymmetric_count(peer, test_rank) + 1U) / 2U;
    }
    require_true(chunk_stats.used_chunked == TRUE,
                 "forced chunk exchange did not use chunking");
    require_true(chunk_stats.send_messages == expected_send_messages &&
                     chunk_stats.recv_messages == expected_recv_messages,
                 "forced chunk message count mismatch");
  }

  FreeSymmetryMpiExchangeResult(&fast_result);
  FreeSymmetryMpiExchangeResult(&chunk_result);
  free_layout(entries, counts, displacements);
}

static uint64_t byte_cap_boundary_count(int source, int destination)
{
  return source == destination ? 0U : test_message_entry_limit;
}

static uint64_t deep_chunk_count(int source, int destination)
{
  uint64_t extra_round;
  if (source == destination) return 0U;
  extra_round = source % 2 == 0 ? 0U : test_message_entry_limit;
  return 3U * test_message_entry_limit + 1U + extra_round;
}

static void assert_byte_cap_boundary_and_deep_chunks(void)
{
  struct SymmetryBasisVector *entries = NULL;
  uint64_t *counts = NULL;
  uint64_t *displacements = NULL;
  struct SymmetryMpiExchangeLayout layout;
  struct SymmetryMpiExchangeResult result;
  struct SymmetryMpiExchangeStats stats;
  uint64_t expected_send_messages = 0U;
  uint64_t expected_recv_messages = 0U;
  int peer;

  memset(&result, 0, sizeof(result));
  build_layout(byte_cap_boundary_count, &entries, &counts, &displacements,
               &layout);
  require_true(
      SymmetryMpiExchangeBasisVectors(
          entries, &layout, test_rank, test_nrank, NULL, &result,
          &stats) == 0,
      "message byte cap boundary exchange failed");
  validate_result(&result, byte_cap_boundary_count);
  require_true(stats.message_entry_limit == test_message_entry_limit &&
                   stats.message_byte_limit ==
                       test_message_entry_limit *
                           sizeof(struct SymmetryBasisVector) &&
                   stats.message_byte_limit <=
                       HPHI_SYMMETRY_EXCHANGE_MESSAGE_BYTES,
               "message byte cap stats mismatch");
  require_true(stats.used_chunked == FALSE,
               "message byte cap boundary unexpectedly used chunking");
  FreeSymmetryMpiExchangeResult(&result);
  free_layout(entries, counts, displacements);

  entries = NULL;
  counts = NULL;
  displacements = NULL;
  memset(&result, 0, sizeof(result));
  build_layout(deep_chunk_count, &entries, &counts, &displacements, &layout);
  require_true(
      SymmetryMpiExchangeBasisVectors(
          entries, &layout, test_rank, test_nrank, NULL, &result,
          &stats) == 0,
      "automatic message byte chunking failed");
  validate_result(&result, deep_chunk_count);
  if (test_nrank > 1) {
    for (peer = 0; peer < test_nrank; peer++) {
      if (peer == test_rank) continue;
      expected_send_messages +=
          (counts[peer] + test_message_entry_limit - 1U) /
          test_message_entry_limit;
      expected_recv_messages +=
          (deep_chunk_count(peer, test_rank) +
           test_message_entry_limit - 1U) /
          test_message_entry_limit;
    }
    require_true(stats.used_chunked == TRUE,
                 "byte cap overflow did not select chunking");
    require_true(stats.send_messages == expected_send_messages &&
                     stats.recv_messages == expected_recv_messages,
                 "deep chunk message count mismatch");
    require_true(stats.max_message_entries == test_message_entry_limit &&
                     stats.max_message_bytes ==
                         test_message_entry_limit *
                             sizeof(struct SymmetryBasisVector) &&
                     stats.max_message_bytes <= stats.message_byte_limit,
                 "deep chunk maximum message size mismatch");
  }
  FreeSymmetryMpiExchangeResult(&result);
  free_layout(entries, counts, displacements);
}

static uint64_t empty_count(int source, int destination)
{
  (void)source;
  (void)destination;
  return 0U;
}

static uint64_t self_count(int source, int destination)
{
  return source == destination ? 5U : 0U;
}

static void assert_empty_and_self_only(void)
{
  uint64_t (*fixtures[2])(int, int) = {empty_count, self_count};
  int fixture;
  for (fixture = 0; fixture < 2; fixture++) {
    struct SymmetryBasisVector *entries = NULL;
    uint64_t *counts = NULL;
    uint64_t *displacements = NULL;
    struct SymmetryMpiExchangeLayout layout;
    struct SymmetryMpiExchangeOptions options;
    struct SymmetryMpiExchangeResult result;
    struct SymmetryMpiExchangeStats stats;
    uint64_t index;
    memset(&result, 0, sizeof(result));
    options.chunk_limit = 2U;
    options.force_chunked = TRUE;
    build_layout(fixtures[fixture], &entries, &counts, &displacements,
                 &layout);
    require_true(
        SymmetryMpiExchangeBasisVectors(
            entries, &layout, test_rank, test_nrank, &options,
            &result, &stats) == 0,
        "empty/self-only exchange failed");
    if (fixture == 0) {
      require_true(result.count == 0U && result.entries == NULL,
                   "all-empty exchange result mismatch");
    } else {
      require_true(result.count == 5U &&
                       result.counts[test_rank] == 5U,
                   "self-only exchange count mismatch");
      for (index = 0U; index < 5U; index++) {
        struct SymmetryBasisVector expected =
            make_entry(test_rank, test_rank, index);
        require_true(field_bits_equal(&result.entries[index], &expected),
                     "self-only payload mismatch");
      }
    }
    require_true(stats.send_messages == 0U && stats.recv_messages == 0U,
                 "empty/self-only exchange posted remote messages");
    FreeSymmetryMpiExchangeResult(&result);
    free_layout(entries, counts, displacements);
  }
}

static void require_all_ranks_failed(int status, const char *label)
{
  int failed = status != 0 ? 1 : 0;
#ifdef MPI
  if (test_mpi_active != FALSE) {
    int failed_min = 0;
    int failed_max = 0;
    require_true(
        MPI_Allreduce(&failed, &failed_min, 1, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD) == MPI_SUCCESS &&
            MPI_Allreduce(&failed, &failed_max, 1, MPI_INT, MPI_MAX,
                          MPI_COMM_WORLD) == MPI_SUCCESS,
        "failure agreement verification failed");
    require_true(failed_min == 1 && failed_max == 1, label);
    return;
  }
#endif
  require_true(failed == 1, label);
}

static void assert_collective_validation_failure(void)
{
  struct SymmetryBasisVector *entries = NULL;
  uint64_t *counts = NULL;
  uint64_t *displacements = NULL;
  struct SymmetryMpiExchangeLayout layout;
  struct SymmetryMpiExchangeResult result;
  uint64_t saved_displacement;
  int failing_rank = test_nrank > 1 ? 1 : 0;
  int status;

  memset(&result, 0, sizeof(result));
  build_layout(asymmetric_count, &entries, &counts, &displacements, &layout);
  saved_displacement = displacements[0];
  if (test_rank == failing_rank) displacements[0] = 1U;
  status = SymmetryMpiExchangeBasisVectors(
      entries, &layout, test_rank, test_nrank, NULL, &result, NULL);
  require_all_ranks_failed(status,
                           "one-rank validation failure was not agreed");
  require_true(result.entries == NULL && result.counts == NULL &&
                   result.displacements == NULL,
               "failed exchange published partial result");
  displacements[0] = saved_displacement;

  require_true(
      SymmetryMpiExchangeBasisVectors(
          entries, &layout, test_rank, test_nrank, NULL, &result, NULL) == 0,
      "valid exchange after agreed failure failed");
  validate_asymmetric_result(&result);
  FreeSymmetryMpiExchangeResult(&result);
  free_layout(entries, counts, displacements);
}

static void assert_chunk_limit_agreement(void)
{
  struct SymmetryBasisVector *entries = NULL;
  uint64_t *counts = NULL;
  uint64_t *displacements = NULL;
  struct SymmetryMpiExchangeLayout layout;
  struct SymmetryMpiExchangeOptions options;
  struct SymmetryMpiExchangeResult result;
  int status;

  if (test_nrank < 2) return;
  memset(&result, 0, sizeof(result));
  build_layout(asymmetric_count, &entries, &counts, &displacements, &layout);
  options.chunk_limit = test_rank == 0 ? 2U : 3U;
  options.force_chunked = TRUE;
  status = SymmetryMpiExchangeBasisVectors(
      entries, &layout, test_rank, test_nrank, &options, &result, NULL);
  require_all_ranks_failed(status,
                           "chunk limit mismatch was not rejected");
  require_true(result.entries == NULL && result.counts == NULL &&
                   result.displacements == NULL,
               "chunk limit mismatch published partial result");

  options.chunk_limit = 2U;
  require_true(
      SymmetryMpiExchangeBasisVectors(
          entries, &layout, test_rank, test_nrank, &options, &result,
          NULL) == 0,
      "valid exchange after chunk limit mismatch failed");
  validate_asymmetric_result(&result);
  FreeSymmetryMpiExchangeResult(&result);
  free_layout(entries, counts, displacements);
}

static void assert_overflow_and_option_validation(void)
{
  struct SymmetryBasisVector dummy;
  uint64_t *counts =
      (uint64_t *)calloc((size_t)test_nrank, sizeof(*counts));
  uint64_t *displacements =
      (uint64_t *)calloc((size_t)test_nrank, sizeof(*displacements));
  struct SymmetryMpiExchangeLayout layout;
  struct SymmetryMpiExchangeOptions options;
  struct SymmetryMpiExchangeResult result;
  int peer;
  int status;

  require_true(counts != NULL && displacements != NULL,
               "overflow fixture allocation failed");
  memset(&dummy, 0, sizeof(dummy));
  memset(&result, 0, sizeof(result));
  counts[0] = UINT64_MAX;
  for (peer = 1; peer < test_nrank; peer++) {
    displacements[peer] = UINT64_MAX;
  }
  layout.nrank = test_nrank;
  layout.count = UINT64_MAX;
  layout.counts = counts;
  layout.displacements = displacements;
  status = SymmetryMpiExchangeBasisVectors(
      &dummy, &layout, test_rank, test_nrank, NULL, &result, NULL);
  require_all_ranks_failed(status, "byte-extent overflow did not fail");

  if (test_nrank > 1) {
    memset(counts, 0, (size_t)test_nrank * sizeof(*counts));
    memset(displacements, 0,
           (size_t)test_nrank * sizeof(*displacements));
    counts[0] = UINT64_MAX;
    counts[1] = 1U;
    displacements[1] = UINT64_MAX;
    layout.count = UINT64_MAX;
    status = SymmetryMpiExchangeBasisVectors(
        &dummy, &layout, test_rank, test_nrank, NULL, &result, NULL);
    require_all_ranks_failed(status, "prefix overflow did not fail");
  }

  memset(counts, 0, (size_t)test_nrank * sizeof(*counts));
  memset(displacements, 0, (size_t)test_nrank * sizeof(*displacements));
  layout.count = 0U;
  options.chunk_limit = (uint64_t)INT_MAX + 1U;
  options.force_chunked = TRUE;
  status = SymmetryMpiExchangeBasisVectors(
      NULL, &layout, test_rank, test_nrank, &options, &result, NULL);
  require_all_ranks_failed(status, "invalid chunk limit did not fail");

  options.chunk_limit = test_message_entry_limit + 1U;
  options.force_chunked = FALSE;
  status = SymmetryMpiExchangeBasisVectors(
      NULL, &layout, test_rank, test_nrank, &options, &result, NULL);
  require_all_ranks_failed(
      status, "chunk limit above compiled message cap did not fail");

  free(counts);
  free(displacements);
}

static unsigned long int make_request_entry(int source,
                                            int destination,
                                            uint64_t index,
                                            unsigned int generation)
{
  return 1UL +
      (unsigned long int)generation * 100000000UL +
      (unsigned long int)source * 100000UL +
      (unsigned long int)destination * 1000UL +
      (unsigned long int)index;
}

static struct SymmetryMpiLookupResponse make_lookup_response(
    unsigned long int request)
{
  struct SymmetryMpiLookupResponse response;
  memset(&response, 0, sizeof(response));
  response.global_beta = request + 17UL;
  response.norm = (double)request * 0.125 + 0.375;
  return response;
}

static int lookup_response_bits_equal(
    const struct SymmetryMpiLookupResponse *lhs,
    const struct SymmetryMpiLookupResponse *rhs)
{
  return memcmp(&lhs->global_beta, &rhs->global_beta,
                sizeof(lhs->global_beta)) == 0 &&
         memcmp(&lhs->norm, &rhs->norm, sizeof(lhs->norm)) == 0;
}

static void build_unsigned_long_layout(
    uint64_t (*count_function)(int, int),
    unsigned int generation,
    unsigned long int **entries,
    uint64_t **counts,
    uint64_t **displacements,
    struct SymmetryMpiExchangeLayout *layout)
{
  uint64_t total = 0U;
  int destination;
  *counts = (uint64_t *)calloc((size_t)test_nrank, sizeof(**counts));
  *displacements =
      (uint64_t *)calloc((size_t)test_nrank, sizeof(**displacements));
  require_true(*counts != NULL && *displacements != NULL,
               "typed request layout allocation failed");
  for (destination = 0; destination < test_nrank; destination++) {
    (*counts)[destination] =
        count_function(test_rank, destination);
    (*displacements)[destination] = total;
    total += (*counts)[destination];
  }
  *entries = total > 0U
      ? (unsigned long int *)malloc((size_t)total * sizeof(**entries))
      : NULL;
  require_true(total == 0U || *entries != NULL,
               "typed request entry allocation failed");
  for (destination = 0; destination < test_nrank; destination++) {
    uint64_t index;
    for (index = 0U; index < (*counts)[destination]; index++) {
      (*entries)[(*displacements)[destination] + index] =
          make_request_entry(test_rank, destination, index, generation);
    }
  }
  layout->nrank = test_nrank;
  layout->count = total;
  layout->counts = *counts;
  layout->displacements = *displacements;
}

static void validate_unsigned_long_result(
    const struct SymmetryMpiUnsignedLongResult *result,
    uint64_t (*count_function)(int, int),
    unsigned int generation)
{
  uint64_t total = 0U;
  int source;
  require_true(result->nrank == test_nrank,
               "typed request receive rank count mismatch");
  for (source = 0; source < test_nrank; source++) {
    uint64_t expected_count = count_function(source, test_rank);
    uint64_t index;
    require_true(result->counts[source] == expected_count,
                 "typed request receive peer count mismatch");
    require_true(result->displacements[source] == total,
                 "typed request source-rank segment order mismatch");
    for (index = 0U; index < expected_count; index++) {
      unsigned long int expected =
          make_request_entry(source, test_rank, index, generation);
      require_true(
          result->entries[result->displacements[source] + index] == expected,
          "typed request payload mismatch");
    }
    total += expected_count;
  }
  require_true(result->count == total,
               "typed request receive total mismatch");
}

static struct SymmetryMpiLookupResponse *build_lookup_responses(
    const struct SymmetryMpiUnsignedLongResult *requests)
{
  struct SymmetryMpiLookupResponse *responses =
      requests->count > 0U
          ? (struct SymmetryMpiLookupResponse *)malloc(
                (size_t)requests->count * sizeof(*responses))
          : NULL;
  uint64_t index;
  require_true(requests->count == 0U || responses != NULL,
               "typed response allocation failed");
  for (index = 0U; index < requests->count; index++) {
    responses[index] = make_lookup_response(requests->entries[index]);
  }
  return responses;
}

static void validate_lookup_response_result(
    const struct SymmetryMpiLookupResponseResult *result,
    const unsigned long int *original_requests,
    const struct SymmetryMpiExchangeLayout *original_layout)
{
  uint64_t total = 0U;
  int source;
  require_true(result->nrank == test_nrank &&
                   result->count == original_layout->count,
               "typed response result shape mismatch");
  for (source = 0; source < test_nrank; source++) {
    uint64_t index;
    require_true(
        result->counts[source] == original_layout->counts[source] &&
            result->displacements[source] ==
                original_layout->displacements[source] &&
            result->displacements[source] == total,
        "typed response known source-rank layout mismatch");
    for (index = 0U; index < result->counts[source]; index++) {
      uint64_t offset = result->displacements[source] + index;
      struct SymmetryMpiLookupResponse expected =
          make_lookup_response(original_requests[offset]);
      require_true(
          lookup_response_bits_equal(&result->entries[offset], &expected),
          "typed response payload mismatch");
    }
    total += result->counts[source];
  }
  require_true(total == result->count,
               "typed response total mismatch");
}

static void validate_unsigned_long_results_equal(
    const struct SymmetryMpiUnsignedLongResult *lhs,
    const struct SymmetryMpiUnsignedLongResult *rhs)
{
  uint64_t index;
  int peer;
  require_true(lhs->nrank == rhs->nrank && lhs->count == rhs->count,
               "typed request fast/chunk shape mismatch");
  for (peer = 0; peer < lhs->nrank; peer++) {
    require_true(lhs->counts[peer] == rhs->counts[peer] &&
                     lhs->displacements[peer] ==
                         rhs->displacements[peer],
                 "typed request fast/chunk schedule mismatch");
  }
  for (index = 0U; index < lhs->count; index++) {
    require_true(lhs->entries[index] == rhs->entries[index],
                 "typed request fast/chunk payload mismatch");
  }
}

static void validate_lookup_response_results_equal(
    const struct SymmetryMpiLookupResponseResult *lhs,
    const struct SymmetryMpiLookupResponseResult *rhs)
{
  uint64_t index;
  int peer;
  require_true(lhs->nrank == rhs->nrank && lhs->count == rhs->count,
               "typed response fast/chunk shape mismatch");
  for (peer = 0; peer < lhs->nrank; peer++) {
    require_true(lhs->counts[peer] == rhs->counts[peer] &&
                     lhs->displacements[peer] ==
                         rhs->displacements[peer],
                 "typed response fast/chunk schedule mismatch");
  }
  for (index = 0U; index < lhs->count; index++) {
    require_true(
        lookup_response_bits_equal(&lhs->entries[index],
                                   &rhs->entries[index]),
        "typed response fast/chunk payload mismatch");
  }
}

static void execute_typed_roundtrip(
    uint64_t (*count_function)(int, int),
    unsigned int generation,
    const struct SymmetryMpiExchangeOptions *request_options,
    const struct SymmetryMpiExchangeOptions *response_options,
    struct SymmetryMpiUnsignedLongResult *request_result,
    struct SymmetryMpiLookupResponseResult *response_result,
    struct SymmetryMpiExchangeStats *request_stats,
    struct SymmetryMpiExchangeStats *response_stats)
{
  unsigned long int *request_entries = NULL;
  struct SymmetryMpiLookupResponse *response_entries = NULL;
  uint64_t *counts = NULL;
  uint64_t *displacements = NULL;
  struct SymmetryMpiExchangeLayout request_layout;
  struct SymmetryMpiExchangeLayout response_send_layout;

  build_unsigned_long_layout(
      count_function, generation, &request_entries,
      &counts, &displacements, &request_layout);
  require_true(
      SymmetryMpiExchangeUnsignedLongs(
          request_entries, &request_layout, test_rank, test_nrank,
          request_options, request_result, request_stats) == 0,
      "typed request exchange failed");
  validate_unsigned_long_result(
      request_result, count_function, generation);

  response_entries = build_lookup_responses(request_result);
  response_send_layout.nrank = request_result->nrank;
  response_send_layout.count = request_result->count;
  response_send_layout.counts = request_result->counts;
  response_send_layout.displacements = request_result->displacements;
  require_true(
      SymmetryMpiExchangeLookupResponsesKnownLayout(
          response_entries, &response_send_layout, &request_layout,
          test_rank, test_nrank, response_options,
          response_result, response_stats) == 0,
      "typed known-layout response exchange failed");
  validate_lookup_response_result(
      response_result, request_entries, &request_layout);

  free(response_entries);
  free(request_entries);
  free(counts);
  free(displacements);
}

static void assert_typed_fast_and_chunked(void)
{
  struct SymmetryMpiUnsignedLongResult fast_requests;
  struct SymmetryMpiUnsignedLongResult chunk_requests;
  struct SymmetryMpiLookupResponseResult fast_responses;
  struct SymmetryMpiLookupResponseResult chunk_responses;
  struct SymmetryMpiExchangeOptions fast_options;
  struct SymmetryMpiExchangeOptions chunk_options;
  struct SymmetryMpiExchangeStats fast_request_stats;
  struct SymmetryMpiExchangeStats chunk_request_stats;
  struct SymmetryMpiExchangeStats fast_response_stats;
  struct SymmetryMpiExchangeStats chunk_response_stats;

  memset(&fast_requests, 0, sizeof(fast_requests));
  memset(&chunk_requests, 0, sizeof(chunk_requests));
  memset(&fast_responses, 0, sizeof(fast_responses));
  memset(&chunk_responses, 0, sizeof(chunk_responses));
  fast_options.chunk_limit = 0U;
  fast_options.force_chunked = FALSE;
  chunk_options.chunk_limit = 2U;
  chunk_options.force_chunked = TRUE;

  execute_typed_roundtrip(
      asymmetric_count, 0U, &fast_options, &fast_options,
      &fast_requests, &fast_responses,
      &fast_request_stats, &fast_response_stats);
  execute_typed_roundtrip(
      asymmetric_count, 0U, &chunk_options, &chunk_options,
      &chunk_requests, &chunk_responses,
      &chunk_request_stats, &chunk_response_stats);
  validate_unsigned_long_results_equal(&fast_requests, &chunk_requests);
  validate_lookup_response_results_equal(&fast_responses, &chunk_responses);
  require_true(
      fast_request_stats.message_entry_limit ==
          test_request_message_entry_limit &&
          fast_response_stats.message_entry_limit ==
              test_response_message_entry_limit &&
          fast_request_stats.message_byte_limit ==
              test_request_message_entry_limit *
                  sizeof(unsigned long int) &&
          fast_response_stats.message_byte_limit ==
              test_response_message_entry_limit *
                  sizeof(struct SymmetryMpiLookupResponse),
      "typed message byte cap stats mismatch");
  if (test_nrank > 1) {
    require_true(
        fast_request_stats.used_chunked == FALSE &&
            fast_response_stats.used_chunked == FALSE &&
            chunk_request_stats.used_chunked == TRUE &&
            chunk_response_stats.used_chunked == TRUE &&
            chunk_request_stats.max_message_bytes <=
                chunk_request_stats.message_byte_limit &&
            chunk_response_stats.max_message_bytes <=
                chunk_response_stats.message_byte_limit,
        "typed fast/chunk stats mismatch");
  }

  FreeSymmetryMpiUnsignedLongResult(&fast_requests);
  FreeSymmetryMpiUnsignedLongResult(&chunk_requests);
  FreeSymmetryMpiLookupResponseResult(&fast_responses);
  FreeSymmetryMpiLookupResponseResult(&chunk_responses);
}

static uint64_t typed_deep_chunk_count(int source, int destination)
{
  uint64_t extra_round;
  if (source == destination) return 0U;
  extra_round =
      source % 2 == 0 ? 0U : test_response_message_entry_limit;
  return 3U * test_response_message_entry_limit + 1U + extra_round;
}

static void assert_typed_deep_chunks(void)
{
  struct SymmetryMpiUnsignedLongResult requests;
  struct SymmetryMpiLookupResponseResult responses;
  struct SymmetryMpiExchangeStats request_stats;
  struct SymmetryMpiExchangeStats response_stats;

  memset(&requests, 0, sizeof(requests));
  memset(&responses, 0, sizeof(responses));
  execute_typed_roundtrip(
      typed_deep_chunk_count, 1U, NULL, NULL,
      &requests, &responses, &request_stats, &response_stats);
  if (test_nrank > 1) {
    require_true(
        request_stats.used_chunked == TRUE &&
            response_stats.used_chunked == TRUE &&
            request_stats.max_message_entries ==
                test_request_message_entry_limit &&
            response_stats.max_message_entries ==
                test_response_message_entry_limit &&
            request_stats.max_message_bytes <=
                HPHI_SYMMETRY_EXCHANGE_MESSAGE_BYTES &&
            response_stats.max_message_bytes <=
                HPHI_SYMMETRY_EXCHANGE_MESSAGE_BYTES,
        "typed deep chunk stats mismatch");
  }
  FreeSymmetryMpiUnsignedLongResult(&requests);
  FreeSymmetryMpiLookupResponseResult(&responses);
}

static void assert_typed_empty_and_self_only(void)
{
  uint64_t (*fixtures[2])(int, int) = {empty_count, self_count};
  int fixture;
  for (fixture = 0; fixture < 2; fixture++) {
    struct SymmetryMpiUnsignedLongResult requests;
    struct SymmetryMpiLookupResponseResult responses;
    struct SymmetryMpiExchangeOptions options;
    struct SymmetryMpiExchangeStats request_stats;
    struct SymmetryMpiExchangeStats response_stats;
    memset(&requests, 0, sizeof(requests));
    memset(&responses, 0, sizeof(responses));
    options.chunk_limit = 1U;
    options.force_chunked = TRUE;
    execute_typed_roundtrip(
        fixtures[fixture], (unsigned int)(fixture + 2),
        &options, &options, &requests, &responses,
        &request_stats, &response_stats);
    require_true(
        request_stats.send_messages == 0U &&
            request_stats.recv_messages == 0U &&
            response_stats.send_messages == 0U &&
            response_stats.recv_messages == 0U,
        "typed empty/self exchange posted remote messages");
    if (fixture == 0) {
      require_true(
          requests.entries == NULL && requests.count == 0U &&
              responses.entries == NULL && responses.count == 0U,
          "typed all-empty result mismatch");
    } else {
      require_true(
          requests.count == 5U && responses.count == 5U,
          "typed self-only result count mismatch");
    }
    FreeSymmetryMpiUnsignedLongResult(&requests);
    FreeSymmetryMpiLookupResponseResult(&responses);
  }
}

static void assert_typed_failure_recovery(void)
{
  unsigned long int *request_entries = NULL;
  struct SymmetryMpiLookupResponse *response_entries = NULL;
  uint64_t *counts = NULL;
  uint64_t *displacements = NULL;
  uint64_t *known_counts = NULL;
  uint64_t *known_displacements = NULL;
  struct SymmetryMpiExchangeLayout request_layout;
  struct SymmetryMpiExchangeLayout response_send_layout;
  struct SymmetryMpiExchangeLayout known_receive_layout;
  struct SymmetryMpiUnsignedLongResult requests;
  struct SymmetryMpiLookupResponseResult responses;
  struct SymmetryMpiExchangeOptions options;
  uint64_t saved_displacement;
  int failing_rank = test_nrank > 1 ? 1 : 0;
  int variant;
  int status;

  memset(&requests, 0, sizeof(requests));
  memset(&responses, 0, sizeof(responses));
  build_unsigned_long_layout(
      asymmetric_count, 4U, &request_entries,
      &counts, &displacements, &request_layout);
  saved_displacement = displacements[0];
  if (test_rank == failing_rank) displacements[0] = 1U;
  status = SymmetryMpiExchangeUnsignedLongs(
      request_entries, &request_layout, test_rank, test_nrank,
      NULL, &requests, NULL);
  require_all_ranks_failed(
      status, "typed request one-rank layout failure was not agreed");
  require_true(
      requests.entries == NULL && requests.counts == NULL &&
          requests.displacements == NULL,
      "typed request failure published partial result");
  displacements[0] = saved_displacement;

  if (test_nrank > 1) {
    options.chunk_limit = test_rank == 0 ? 1U : 2U;
    options.force_chunked = TRUE;
    status = SymmetryMpiExchangeUnsignedLongs(
        request_entries, &request_layout, test_rank, test_nrank,
        &options, &requests, NULL);
    require_all_ranks_failed(
        status, "typed request option mismatch was not rejected");
    require_true(
        requests.entries == NULL && requests.counts == NULL &&
            requests.displacements == NULL,
        "typed request option failure published partial result");
  }

  options.chunk_limit = 1U;
  options.force_chunked = TRUE;
  require_true(
      SymmetryMpiExchangeUnsignedLongs(
          request_entries, &request_layout, test_rank, test_nrank,
          &options, &requests, NULL) == 0,
      "valid typed request after collective failure failed");
  validate_unsigned_long_result(&requests, asymmetric_count, 4U);
  response_entries = build_lookup_responses(&requests);
  response_send_layout.nrank = requests.nrank;
  response_send_layout.count = requests.count;
  response_send_layout.counts = requests.counts;
  response_send_layout.displacements = requests.displacements;

  known_counts =
      (uint64_t *)malloc((size_t)test_nrank * sizeof(*known_counts));
  known_displacements = (uint64_t *)malloc(
      (size_t)test_nrank * sizeof(*known_displacements));
  require_true(known_counts != NULL && known_displacements != NULL,
               "known response layout copy allocation failed");
  memcpy(known_counts, counts,
         (size_t)test_nrank * sizeof(*known_counts));
  memcpy(known_displacements, displacements,
         (size_t)test_nrank * sizeof(*known_displacements));
  known_receive_layout.nrank = test_nrank;
  known_receive_layout.count = request_layout.count;
  known_receive_layout.counts = known_counts;
  known_receive_layout.displacements = known_displacements;

  for (variant = 0; variant < 4; variant++) {
    if (test_rank == failing_rank) {
      if (variant == 0) {
        int peer;
        known_counts[test_rank]++;
        for (peer = test_rank + 1; peer < test_nrank; peer++) {
          known_displacements[peer]++;
        }
        known_receive_layout.count++;
      } else if (variant == 1) {
        known_displacements[0] = 1U;
      } else if (variant == 2) {
        known_receive_layout.count++;
      } else {
        known_receive_layout.nrank = test_nrank + 1;
      }
    }
    status = SymmetryMpiExchangeLookupResponsesKnownLayout(
        response_entries, &response_send_layout, &known_receive_layout,
        test_rank, test_nrank, &options, &responses, NULL);
    require_all_ranks_failed(
        status, "invalid known response layout was not rejected");
    require_true(
        responses.entries == NULL && responses.counts == NULL &&
            responses.displacements == NULL,
        "known response layout failure published partial result");
    known_receive_layout.nrank = test_nrank;
    known_receive_layout.count = request_layout.count;
    memcpy(known_counts, counts,
           (size_t)test_nrank * sizeof(*known_counts));
    memcpy(known_displacements, displacements,
           (size_t)test_nrank * sizeof(*known_displacements));
  }

  require_true(
      SymmetryMpiExchangeLookupResponsesKnownLayout(
          response_entries, &response_send_layout, &known_receive_layout,
          test_rank, test_nrank, &options, &responses, NULL) == 0,
      "valid known response exchange after failure failed");
  validate_lookup_response_result(
      &responses, request_entries, &request_layout);

  FreeSymmetryMpiLookupResponseResult(&responses);
  FreeSymmetryMpiUnsignedLongResult(&requests);
  free(response_entries);
  free(known_counts);
  free(known_displacements);
  free(request_entries);
  free(counts);
  free(displacements);
}

static void assert_typed_overflow_failure(void)
{
  unsigned long int dummy = 0UL;
  uint64_t *counts =
      (uint64_t *)calloc((size_t)test_nrank, sizeof(*counts));
  uint64_t *displacements =
      (uint64_t *)calloc((size_t)test_nrank, sizeof(*displacements));
  struct SymmetryMpiExchangeLayout layout;
  struct SymmetryMpiUnsignedLongResult result;
  int peer;
  int status;

  require_true(counts != NULL && displacements != NULL,
               "typed overflow fixture allocation failed");
  memset(&result, 0, sizeof(result));
  counts[0] = UINT64_MAX;
  for (peer = 1; peer < test_nrank; peer++) {
    displacements[peer] = UINT64_MAX;
  }
  layout.nrank = test_nrank;
  layout.count = UINT64_MAX;
  layout.counts = counts;
  layout.displacements = displacements;
  status = SymmetryMpiExchangeUnsignedLongs(
      &dummy, &layout, test_rank, test_nrank, NULL, &result, NULL);
  require_all_ranks_failed(status,
                           "typed request byte overflow did not fail");
  require_true(
      result.entries == NULL && result.counts == NULL &&
          result.displacements == NULL,
      "typed overflow failure published partial result");
  free(counts);
  free(displacements);
}

enum SampleSortFixture {
  SAMPLE_SORT_RANDOM = 0,
  SAMPLE_SORT_SKEW = 1,
  SAMPLE_SORT_ALREADY_SORTED = 2,
  SAMPLE_SORT_REVERSE_SORTED = 3,
  SAMPLE_SORT_ONE_RANK_ONLY = 4,
  SAMPLE_SORT_LOCAL_EMPTY = 5,
  SAMPLE_SORT_ALL_EMPTY = 6
};

static uint64_t sample_fixture_count(enum SampleSortFixture fixture)
{
  if (fixture == SAMPLE_SORT_ALL_EMPTY) return 0U;
  if (fixture == SAMPLE_SORT_ONE_RANK_ONLY) return 96U;
  return (uint64_t)test_nrank * 3U + 11U;
}

static void assert_sample_count_warning_policy(void)
{
  const uint64_t default_count =
      (uint64_t)HPHI_SYMMETRY_SAMPLE_COUNT_DEFAULT;
  const uint64_t warning_threshold =
      (uint64_t)HPHI_SYMMETRY_SAMPLE_COUNT_GLOBAL_WARNING_THRESHOLD;
  uint64_t expected_single =
      default_count < warning_threshold ? default_count : warning_threshold;
  uint64_t selected = UINT64_MAX;
  int warning = -1;

  require_true(
      SymmetrySelectSampleCountPolicy(0U, &selected, &warning) == 0 &&
          selected == 0U && warning == FALSE,
      "zero-rank sample policy mismatch");
  require_true(
      SymmetrySelectSampleCountPolicy(1U, &selected, &warning) == 0 &&
          selected == expected_single && warning == FALSE,
      "single-rank sample policy mismatch");
  if (warning_threshold >= default_count) {
    require_true(
        SymmetrySelectSampleCountPolicy(
            warning_threshold / default_count,
            &selected, &warning) == 0 &&
            selected == default_count && warning == FALSE,
        "default sample policy boundary mismatch");
  }
  require_true(
      SymmetrySelectSampleCountPolicy(
          warning_threshold, &selected, &warning) == 0 &&
          selected == 1U && warning == FALSE,
      "global sample warning threshold mismatch");
  require_true(
      warning_threshold < UINT64_MAX &&
          SymmetrySelectSampleCountPolicy(
              warning_threshold + 1U,
              &selected, &warning) == 0 &&
          selected == 1U && warning == TRUE,
      "global sample warning continuation mismatch");
  require_true(
      SymmetrySelectSampleCountPolicy(
          UINT64_MAX, &selected, &warning) == 0 &&
          selected == 1U && warning == TRUE,
      "maximum-rank sample warning continuation mismatch");
  require_true(
      SymmetrySelectSampleCountPolicy(1U, NULL, &warning) == -1 &&
          SymmetrySelectSampleCountPolicy(1U, &selected, NULL) == -1,
      "invalid sample policy output was not rejected");
}

static int sample_fixture_owner(enum SampleSortFixture fixture,
                                uint64_t ordinal,
                                uint64_t global_count)
{
  int even_rank_count;
  if (test_nrank == 1) return 0;
  switch (fixture) {
  case SAMPLE_SORT_RANDOM:
    return (int)((ordinal * 5U + 3U) % (uint64_t)test_nrank);
  case SAMPLE_SORT_SKEW:
    if (ordinal + (uint64_t)(test_nrank - 1) < global_count) return 0;
    return 1 + (int)(
        ordinal - (global_count - (uint64_t)(test_nrank - 1)));
  case SAMPLE_SORT_ALREADY_SORTED:
    return (int)((ordinal * (uint64_t)test_nrank) / global_count);
  case SAMPLE_SORT_REVERSE_SORTED:
    return (int)((global_count - 1U - ordinal) %
                 (uint64_t)test_nrank);
  case SAMPLE_SORT_ONE_RANK_ONLY:
    return test_nrank - 1;
  case SAMPLE_SORT_LOCAL_EMPTY:
    even_rank_count = (test_nrank + 1) / 2;
    return 2 * (int)(ordinal % (uint64_t)even_rank_count);
  case SAMPLE_SORT_ALL_EMPTY:
    return 0;
  }
  return 0;
}

static struct SymmetryBasisVector make_sample_sort_entry(uint64_t ordinal)
{
  struct SymmetryBasisVector entry;
  unsigned long int key =
      (unsigned long int)(UINT64_C(1009) + ordinal * UINT64_C(17));
  memset(&entry, 0, sizeof(entry));
  entry.rep_state = key;
  entry.orbit_size = (unsigned int)(1U + ordinal % 29U);
  entry.stabilizer_size = (unsigned int)(1U + ordinal % 11U);
  entry.norm = 0.375 + (double)ordinal * 0.125;
  entry.stabilizer_character_sum =
      (0.625 + (double)ordinal * 0.25) +
      (-0.875 + (double)ordinal * 0.0625) * I;
  entry.diagonal = -1.25 + (double)ordinal * 0.5;
  return entry;
}

static void build_sample_sort_run(enum SampleSortFixture fixture,
                                  struct SymmetryBasisRun *run)
{
  uint64_t global_count = sample_fixture_count(fixture);
  uint64_t local_count = 0U;
  uint64_t ordinal;
  uint64_t local_index = 0U;
  int reverse =
      fixture == SAMPLE_SORT_RANDOM ||
      fixture == SAMPLE_SORT_SKEW ||
      fixture == SAMPLE_SORT_REVERSE_SORTED ||
      fixture == SAMPLE_SORT_ONE_RANK_ONLY ||
      fixture == SAMPLE_SORT_LOCAL_EMPTY;

  memset(run, 0, sizeof(*run));
  for (ordinal = 0U; ordinal < global_count; ordinal++) {
    if (sample_fixture_owner(fixture, ordinal, global_count) == test_rank) {
      local_count++;
    }
  }
  require_true(local_count < (uint64_t)ULONG_MAX,
               "sample-sort fixture count overflow");
  if (fixture == SAMPLE_SORT_ALL_EMPTY && test_rank % 2 == 0) {
    return;
  }
  run->entries = (struct SymmetryBasisVector *)calloc(
      (size_t)local_count + 1U, sizeof(*run->entries));
  require_true(run->entries != NULL,
               "sample-sort fixture allocation failed");
  run->count = (unsigned long int)local_count;
  run->capacity = run->count + 1UL;

  if (reverse != FALSE) {
    ordinal = global_count;
    while (ordinal > 0U) {
      ordinal--;
      if (sample_fixture_owner(fixture, ordinal, global_count) != test_rank) {
        continue;
      }
      local_index++;
      run->entries[local_index] = make_sample_sort_entry(ordinal);
    }
  } else {
    for (ordinal = 0U; ordinal < global_count; ordinal++) {
      if (sample_fixture_owner(fixture, ordinal, global_count) != test_rank) {
        continue;
      }
      local_index++;
      run->entries[local_index] = make_sample_sort_entry(ordinal);
    }
  }
  require_true(local_index == local_count,
               "sample-sort fixture fill count mismatch");
}

static void require_all_ranks_u64_equal(uint64_t value, const char *label)
{
#ifdef MPI
  if (test_mpi_active != FALSE) {
    uint64_t minimum = 0U;
    uint64_t maximum = 0U;
    require_true(
        MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T,
                      MPI_MIN, MPI_COMM_WORLD) == MPI_SUCCESS &&
            MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T,
                          MPI_MAX, MPI_COMM_WORLD) == MPI_SUCCESS,
        "u64 rank agreement reduction failed");
    require_true(minimum == maximum, label);
    return;
  }
#endif
  (void)value;
  (void)label;
}

static void assert_run_fields_equal(const struct SymmetryBasisRun *lhs,
                                    const struct SymmetryBasisRun *rhs,
                                    const char *label)
{
  unsigned long int index;
  require_true(lhs->count == rhs->count &&
                   lhs->capacity == rhs->capacity,
               label);
  for (index = 0UL; index <= lhs->count; index++) {
    require_true(field_bits_equal(&lhs->entries[index],
                                  &rhs->entries[index]),
                 label);
  }
}

static void validate_sample_sort_global_result(
    const struct SymmetryBasisRun *run,
    enum SampleSortFixture fixture,
    const char *label)
{
  struct SymmetryMpiExchangeLayout layout;
  struct SymmetryMpiExchangeResult result;
  uint64_t *counts = (uint64_t *)calloc(
      (size_t)test_nrank, sizeof(*counts));
  uint64_t *displacements = (uint64_t *)calloc(
      (size_t)test_nrank, sizeof(*displacements));
  uint64_t global_count = sample_fixture_count(fixture);
  uint64_t displacement = 0U;
  uint64_t ordinal;
  int destination;

  require_true(counts != NULL && displacements != NULL,
               "sample-sort gather schedule allocation failed");
  counts[0] = (uint64_t)run->count;
  for (destination = 0; destination < test_nrank; destination++) {
    displacements[destination] = displacement;
    displacement += counts[destination];
  }
  memset(&result, 0, sizeof(result));
  layout.nrank = test_nrank;
  layout.count = (uint64_t)run->count;
  layout.counts = counts;
  layout.displacements = displacements;
  require_true(
      SymmetryMpiExchangeBasisVectors(
          run->count > 0UL ? run->entries + 1 : run->entries,
          &layout, test_rank, test_nrank, NULL, &result, NULL) == 0,
      "sample-sort result gather failed");
  if (test_rank == 0) {
    require_true(result.count == global_count, label);
    for (ordinal = 0U; ordinal < global_count; ordinal++) {
      struct SymmetryBasisVector expected =
          make_sample_sort_entry(ordinal);
      require_true(field_bits_equal(&result.entries[ordinal], &expected),
                   label);
    }
  } else {
    require_true(result.count == 0U, label);
  }
  FreeSymmetryMpiExchangeResult(&result);
  free(counts);
  free(displacements);
}

static void assert_sample_sort_fixture(enum SampleSortFixture fixture,
                                       const char *label)
{
  struct SymmetryBasisRun first;
  struct SymmetryBasisRun second;
  struct SymmetryBasisDistributionStats first_stats;
  struct SymmetryBasisDistributionStats second_stats;
  uint64_t expected_global_count = sample_fixture_count(fixture);
  unsigned long int index;

  build_sample_sort_run(fixture, &first);
  build_sample_sort_run(fixture, &second);
  require_true(
      SymmetrySampleSortBasisRun(
          &first, test_rank, test_nrank, &first_stats) == 0,
      label);
  require_true(first.entries != NULL &&
                   first.capacity == first.count + 1UL,
               "sample-sort output ownership mismatch");
  for (index = 1UL; index < first.count; index++) {
    require_true(first.entries[index].rep_state <
                     first.entries[index + 1UL].rep_state,
                 "sample-sort local ordering mismatch");
  }
  require_true(first_stats.global_entries == expected_global_count &&
                   first_stats.range_entries == (uint64_t)first.count &&
                   first_stats.range_entries <=
                       first_stats.bucket_entry_upper_bound &&
                   first_stats.bucket_entry_upper_bound ==
                       first_stats.bucket_sample_entries *
                           first_stats.global_sample_gap_max +
                       first_stats.global_sample_gap_sum &&
                   first_stats
                           .distribution_memory_warning_byte_threshold ==
                       (uint64_t)HPHI_SYMMETRY_MEMORY_WARN_BYTES &&
                   first_stats.distribution_memory_byte_limit ==
                       HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES &&
                   first_stats.sample_gather_message_byte_limit ==
                       test_message_entry_limit *
                           sizeof(struct SymmetryBasisVector) &&
                   first_stats.sample_gather_max_message_bytes <=
                       first_stats.sample_gather_message_byte_limit &&
                   first_stats.range_send_entries ==
                       first_stats.local_survivor_entries &&
                   first_stats.range_recv_entries ==
                       first_stats.range_entries &&
                   first_stats.range_exchange_message_byte_limit ==
                       test_message_entry_limit *
                           sizeof(struct SymmetryBasisVector) &&
                   first_stats.range_exchange_max_message_bytes <=
                       first_stats.range_exchange_message_byte_limit &&
                   first_stats.sort_temporary_peak_bytes <=
                       HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES &&
                   first_stats.rebalance_send_entries == 0U &&
                   first_stats.rebalance_recv_entries == 0U &&
                   first_stats.rebalance_temporary_peak_bytes == 0U,
               "sample-sort stats contract mismatch");
  if (fixture == SAMPLE_SORT_ONE_RANK_ONLY) {
    uint64_t expected_samples = 0U;
    int expected_warning = FALSE;
    require_true(
        SymmetrySelectSampleCountPolicy(
            first_stats.nonempty_rank_count,
            &expected_samples, &expected_warning) == 0 &&
            expected_warning == FALSE &&
            first_stats.samples_per_nonempty_rank == expected_samples &&
            first_stats.global_sample_entries == expected_samples &&
            first_stats.global_sample_gap_max > 1U &&
            (test_nrank == 1 ||
             (first_stats.sample_gather_used_chunked ==
                  (expected_samples > test_message_entry_limit
                       ? TRUE : FALSE) &&
              first_stats.range_exchange_used_chunked == TRUE)),
                 "regular sample gap fixture did not downsample");
  }
  require_all_ranks_u64_equal(
      first_stats.splitter_digest,
      "sample-sort splitter digest differs across ranks");
  validate_sample_sort_global_result(&first, fixture, label);

  require_true(
      SymmetrySampleSortBasisRun(
          &second, test_rank, test_nrank, &second_stats) == 0,
      "repeat sample-sort execution failed");
  assert_run_fields_equal(&first, &second,
                          "repeat sample-sort range mismatch");
  require_true(
      first_stats.nonempty_rank_count == second_stats.nonempty_rank_count &&
          first_stats.samples_per_nonempty_rank ==
              second_stats.samples_per_nonempty_rank &&
          first_stats.global_sample_entries ==
              second_stats.global_sample_entries &&
          first_stats.sample_gather_used_chunked ==
              second_stats.sample_gather_used_chunked &&
          first_stats.sample_gather_message_byte_limit ==
              second_stats.sample_gather_message_byte_limit &&
          first_stats.sample_gather_max_message_bytes ==
              second_stats.sample_gather_max_message_bytes &&
          first_stats.local_sample_gap == second_stats.local_sample_gap &&
          first_stats.bucket_sample_entries ==
              second_stats.bucket_sample_entries &&
          first_stats.bucket_entry_upper_bound ==
              second_stats.bucket_entry_upper_bound &&
          first_stats.range_send_entries ==
              second_stats.range_send_entries &&
          first_stats.range_recv_entries ==
              second_stats.range_recv_entries &&
          first_stats.range_exchange_used_chunked ==
              second_stats.range_exchange_used_chunked &&
          first_stats.range_exchange_message_byte_limit ==
              second_stats.range_exchange_message_byte_limit &&
          first_stats.range_exchange_max_message_bytes ==
              second_stats.range_exchange_max_message_bytes &&
          first_stats.splitter_digest == second_stats.splitter_digest &&
          first_stats.range_digest == second_stats.range_digest,
      "sample-sort determinism stats mismatch");
  validate_sample_sort_global_result(&second, fixture, label);

  FreeSymmetryBasisRun(&first);
  FreeSymmetryBasisRun(&second);
}

static void assert_sample_sort_duplicate_failure(void)
{
  struct SymmetryBasisRun run;
  struct SymmetryBasisVector *original_entries;
  unsigned long int original_count;
  int status;

  memset(&run, 0, sizeof(run));
  run.count = test_nrank > 1
      ? (test_rank < 2 ? 1UL : 0UL)
      : 2UL;
  run.capacity = run.count + 1UL;
  run.entries = (struct SymmetryBasisVector *)calloc(
      (size_t)run.capacity, sizeof(*run.entries));
  require_true(run.entries != NULL,
               "duplicate sample-sort fixture allocation failed");
  if (run.count > 0UL) {
    run.entries[1] = make_sample_sort_entry(7U);
  }
  if (run.count > 1UL) {
    run.entries[2] = make_sample_sort_entry(7U);
  }
  original_entries = run.entries;
  original_count = run.count;
  status = SymmetrySampleSortBasisRun(
      &run, test_rank, test_nrank, NULL);
  require_all_ranks_failed(
      status, "sample-sort duplicate was not rejected collectively");
  require_true(run.entries == original_entries &&
                   run.count == original_count &&
                   run.capacity == original_count + 1UL,
               "duplicate failure consumed run ownership");
  FreeSymmetryBasisRun(&run);

  build_sample_sort_run(SAMPLE_SORT_ALREADY_SORTED, &run);
  require_true(
      SymmetrySampleSortBasisRun(
          &run, test_rank, test_nrank, NULL) == 0,
      "valid sample sort after duplicate failure failed");
  FreeSymmetryBasisRun(&run);
}

static void assert_sample_sort_memory_cap_failure(void)
{
  struct SymmetryBasisRun run;
  struct SymmetryBasisVector *original_entries;
  unsigned long int original_count;
  const uint64_t global_count = UINT64_C(5000);
  uint64_t ordinal;
  int status;

  memset(&run, 0, sizeof(run));
  run.count = test_rank == 0
      ? (unsigned long int)global_count : 0UL;
  run.capacity = run.count + 1UL;
  run.entries = (struct SymmetryBasisVector *)calloc(
      (size_t)run.capacity, sizeof(*run.entries));
  require_true(run.entries != NULL,
               "memory-cap sample-sort fixture allocation failed");
  if (test_rank == 0) {
    for (ordinal = 0U; ordinal < global_count; ordinal++) {
      run.entries[ordinal + 1U] =
          make_sample_sort_entry(global_count - 1U - ordinal);
    }
  }
  original_entries = run.entries;
  original_count = run.count;
  status = SymmetrySampleSortBasisRun(
      &run, test_rank, test_nrank, NULL);
  require_all_ranks_failed(
      status, "sample-sort memory cap was not rejected collectively");
  require_true(run.entries == original_entries &&
                   run.count == original_count &&
                   run.capacity == original_count + 1UL,
               "memory-cap failure consumed run ownership");
  FreeSymmetryBasisRun(&run);

  build_sample_sort_run(SAMPLE_SORT_RANDOM, &run);
  require_true(
      SymmetrySampleSortBasisRun(
          &run, test_rank, test_nrank, NULL) == 0,
      "valid sample sort after memory-cap failure failed");
  FreeSymmetryBasisRun(&run);
}

static void exact_source_range(uint64_t global_count,
                               int one_rank_only,
                               int rank,
                               uint64_t *offset,
                               uint64_t *count)
{
  if (one_rank_only != FALSE) {
    *offset = 0U;
    *count = rank == test_nrank - 1 ? global_count : 0U;
    return;
  }
  if (test_nrank == 1) {
    *offset = 0U;
    *count = global_count;
    return;
  }
  {
    uint64_t head = (global_count + 1U) / 2U;
    unsigned long int tail_offset = 0UL;
    unsigned long int tail_count = 0UL;
    if (rank == 0) {
      *offset = 0U;
      *count = head;
      return;
    }
    require_true(
        global_count - head <= (uint64_t)ULONG_MAX &&
            SymmetryBlockRange(
                (unsigned long int)(global_count - head),
                rank - 1, test_nrank - 1,
                &tail_offset, &tail_count) == 0,
        "exact source range calculation failed");
    *offset = head + (uint64_t)tail_offset;
    *count = (uint64_t)tail_count;
  }
}

static void build_exact_range_run(uint64_t global_count,
                                  int one_rank_only,
                                  struct SymmetryBasisRun *run)
{
  uint64_t offset;
  uint64_t count;
  uint64_t index;

  memset(run, 0, sizeof(*run));
  exact_source_range(global_count, one_rank_only, test_rank,
                     &offset, &count);
  require_true(count < (uint64_t)ULONG_MAX &&
                   count <= (uint64_t)SIZE_MAX /
                       sizeof(*run->entries) - 1U,
               "exact range fixture count overflow");
  if (global_count == 0U && test_rank % 2 == 0) return;
  run->entries = (struct SymmetryBasisVector *)calloc(
      (size_t)count + 1U, sizeof(*run->entries));
  require_true(run->entries != NULL,
               "exact range fixture allocation failed");
  run->count = (unsigned long int)count;
  run->capacity = run->count + 1UL;
  for (index = 0U; index < count; index++) {
    run->entries[index + 1U] =
        make_sample_sort_entry(offset + index);
  }
}

static void validate_exact_global_result(
    const struct SymmetryBasisRun *run,
    uint64_t global_count,
    const char *label)
{
  struct SymmetryMpiExchangeLayout layout;
  struct SymmetryMpiExchangeResult result;
  uint64_t *counts = (uint64_t *)calloc(
      (size_t)test_nrank, sizeof(*counts));
  uint64_t *displacements = (uint64_t *)calloc(
      (size_t)test_nrank, sizeof(*displacements));
  uint64_t displacement = 0U;
  uint64_t ordinal;
  int destination;

  require_true(counts != NULL && displacements != NULL,
               "exact result gather schedule allocation failed");
  counts[0] = (uint64_t)run->count;
  for (destination = 0; destination < test_nrank; destination++) {
    displacements[destination] = displacement;
    displacement += counts[destination];
  }
  memset(&result, 0, sizeof(result));
  layout.nrank = test_nrank;
  layout.count = (uint64_t)run->count;
  layout.counts = counts;
  layout.displacements = displacements;
  require_true(
      SymmetryMpiExchangeBasisVectors(
          run->count > 0UL ? run->entries + 1 : run->entries,
          &layout, test_rank, test_nrank, NULL, &result, NULL) == 0,
      "exact result gather failed");
  if (test_rank == 0) {
    require_true(result.count == global_count, label);
    for (ordinal = 0U; ordinal < global_count; ordinal++) {
      struct SymmetryBasisVector expected =
          make_sample_sort_entry(ordinal);
      require_true(field_bits_equal(&result.entries[ordinal], &expected),
                   label);
    }
  } else {
    require_true(result.count == 0U, label);
  }
  FreeSymmetryMpiExchangeResult(&result);
  free(counts);
  free(displacements);
}

static void assert_exact_rebalance_fixture(uint64_t global_count,
                                           int one_rank_only,
                                           const char *label)
{
  struct SymmetryBasisRun run;
  struct SymmetryBasisOwnership ownership;
  struct SymmetryBasisDistributionStats stats;
  struct SymmetryBasisVector sentinel;
  uint64_t input_count;
  unsigned long int expected_offset;
  unsigned long int expected_count;
  unsigned long int index;
  int peer;

  require_true(global_count <= (uint64_t)ULONG_MAX,
               "exact fixture global count overflow");
  build_exact_range_run(global_count, one_rank_only, &run);
  input_count = (uint64_t)run.count;
  memset(&ownership, 0, sizeof(ownership));
  memset(&stats, 0, sizeof(stats));
  stats.global_entries = global_count;
  stats.range_send_entries = UINT64_C(101);
  stats.range_recv_entries = UINT64_C(103);
  stats.range_exchange_used_chunked = TRUE;
  stats.range_exchange_message_byte_limit = UINT64_C(105);
  stats.range_exchange_max_message_bytes = UINT64_C(106);
  stats.sort_temporary_peak_bytes = 107U;
  stats.distribution_memory_warning_byte_threshold = UINT64_C(108);
  stats.distribution_memory_byte_limit = UINT64_C(109);
  stats.splitter_digest = UINT64_C(110);
  stats.range_digest = UINT64_C(113);

  require_true(
      SymmetryExactRebalanceBasisRun(
          &run, test_rank, test_nrank, &ownership, &stats) == 0,
      label);
  require_true(
      SymmetryBlockRange((unsigned long int)global_count,
                         test_rank, test_nrank,
                         &expected_offset, &expected_count) == 0,
      "exact expected block calculation failed");
  require_true(ownership.dim == (unsigned long int)global_count &&
                   ownership.local_offset == expected_offset &&
                   ownership.local_dim == expected_count &&
                   ownership.rank_offsets != NULL &&
                   run.count == expected_count &&
                   run.capacity == expected_count + 1UL &&
                   run.entries != NULL,
               "exact ownership contract mismatch");
  for (peer = 0; peer < test_nrank; peer++) {
    unsigned long int peer_offset;
    unsigned long int peer_count;
    require_true(
        SymmetryBlockRange((unsigned long int)global_count,
                           peer, test_nrank,
                           &peer_offset, &peer_count) == 0 &&
            ownership.rank_offsets[peer] == peer_offset &&
            ownership.rank_offsets[peer + 1] ==
                peer_offset + peer_count,
        "exact rank offset mismatch");
  }
  memset(&sentinel, 0, sizeof(sentinel));
  require_true(field_bits_equal(&run.entries[0], &sentinel),
               "exact output sentinel mismatch");
  for (index = 1UL; index <= run.count; index++) {
    struct SymmetryBasisVector expected = make_sample_sort_entry(
        (uint64_t)expected_offset + (uint64_t)index - 1U);
    require_true(field_bits_equal(&run.entries[index], &expected),
                 label);
  }
  require_true(stats.global_entries == global_count &&
                   stats.range_send_entries == UINT64_C(101) &&
                   stats.range_recv_entries == UINT64_C(103) &&
                   stats.range_exchange_used_chunked == TRUE &&
                   stats.range_exchange_message_byte_limit ==
                       UINT64_C(105) &&
                   stats.range_exchange_max_message_bytes ==
                       UINT64_C(106) &&
                   stats.sort_temporary_peak_bytes == 107U &&
                   stats.distribution_memory_warning_byte_threshold ==
                       (uint64_t)HPHI_SYMMETRY_MEMORY_WARN_BYTES &&
                   stats.distribution_memory_byte_limit ==
                       HPHI_SYMMETRY_DISTRIBUTION_MEMORY_BYTES &&
                   stats.splitter_digest == UINT64_C(110) &&
                   stats.range_digest == UINT64_C(113) &&
                   stats.rebalance_send_entries == input_count &&
                   stats.rebalance_recv_entries ==
                       (uint64_t)expected_count &&
                   stats.rebalance_exchange_message_byte_limit ==
                       test_message_entry_limit *
                           sizeof(struct SymmetryBasisVector) &&
                   stats.rebalance_exchange_max_message_bytes <=
                       stats.rebalance_exchange_message_byte_limit &&
                   stats.rebalance_temporary_peak_bytes >=
                       ((size_t)input_count +
                        (size_t)expected_count * 2U + 2U) *
                           sizeof(struct SymmetryBasisVector),
               "exact rebalance stats mismatch");
  if (one_rank_only != FALSE && test_nrank > 1 &&
      global_count > (uint64_t)test_nrank *
          test_message_entry_limit) {
    require_true(
        (uint64_t)expected_count > test_message_entry_limit &&
            stats.rebalance_exchange_used_chunked == TRUE,
        "exact deep-chunk fixture did not exceed the message cap");
  }
  validate_exact_global_result(&run, global_count, label);

  FreeSymmetryBasisOwnership(&ownership);
  require_true(ownership.rank_offsets == NULL &&
                   ownership.dim == 0UL &&
                   ownership.local_offset == 0UL &&
                   ownership.local_dim == 0UL,
               "exact ownership cleanup mismatch");
  FreeSymmetryBasisRun(&run);
}

static void assert_block_range_contract(void)
{
  static const unsigned long int dimensions[] = {0UL, 1UL, 7UL, 23UL};
  size_t dimension_index;
  for (dimension_index = 0U;
       dimension_index < sizeof(dimensions) / sizeof(dimensions[0]);
       dimension_index++) {
    unsigned long int previous_end = 0UL;
    int rank;
    for (rank = 0; rank < test_nrank; rank++) {
      unsigned long int offset;
      unsigned long int count;
      require_true(
          SymmetryBlockRange(dimensions[dimension_index],
                             rank, test_nrank,
                             &offset, &count) == 0 &&
              offset == previous_end,
          "block range prefix mismatch");
      if (rank > 0) {
        unsigned long int previous_offset;
        unsigned long int previous_count;
        require_true(
            SymmetryBlockRange(dimensions[dimension_index],
                               rank - 1, test_nrank,
                               &previous_offset, &previous_count) == 0 &&
                previous_count >= count,
            "block range balance mismatch");
      }
      previous_end = offset + count;
    }
    require_true(previous_end == dimensions[dimension_index],
                 "block range final dimension mismatch");
  }
  {
    unsigned long int offset = 0UL;
    unsigned long int count = 0UL;
    require_true(
        SymmetryBlockRange(1UL, -1, test_nrank,
                           &offset, &count) != 0 &&
            SymmetryBlockRange(1UL, test_nrank, test_nrank,
                               &offset, &count) != 0 &&
            SymmetryBlockRange(1UL, 0, 0,
                               &offset, &count) != 0 &&
            SymmetryBlockRange(1UL, 0, test_nrank,
                               NULL, &count) != 0 &&
            SymmetryBlockRange(1UL, 0, test_nrank,
                               &offset, NULL) != 0,
        "invalid block range input was accepted");
  }
}

static void assert_exact_rebalance_failure_recovery(void)
{
  struct SymmetryBasisRun run;
  struct SymmetryBasisOwnership ownership;
  struct SymmetryBasisDistributionStats stats;
  struct SymmetryBasisDistributionStats saved_stats;
  struct SymmetryBasisVector *saved_entries;
  unsigned long int saved_count;
  unsigned long int saved_capacity;
  uint64_t global_count = (uint64_t)test_nrank * 3U + 1U;
  int failing_rank = test_nrank > 1 ? 1 : 0;
  int status;

  build_exact_range_run(global_count, FALSE, &run);
  saved_entries = run.entries;
  saved_count = run.count;
  saved_capacity = run.capacity;
  memset(&ownership, 0, sizeof(ownership));
  if (test_rank == failing_rank) {
    ownership.rank_offsets =
        (unsigned long int *)calloc(1U, sizeof(*ownership.rank_offsets));
    require_true(ownership.rank_offsets != NULL,
                 "invalid ownership fixture allocation failed");
  }
  memset(&stats, 0x5a, sizeof(stats));
  memcpy(&saved_stats, &stats, sizeof(saved_stats));
  status = SymmetryExactRebalanceBasisRun(
      &run, test_rank, test_nrank, &ownership, &stats);
  require_all_ranks_failed(
      status, "one-rank nonempty ownership was not rejected");
  require_true(run.entries == saved_entries &&
                   run.count == saved_count &&
                   run.capacity == saved_capacity &&
                   memcmp(&stats, &saved_stats, sizeof(stats)) == 0,
               "failed exact rebalance changed run or stats");
  if (test_rank == failing_rank) {
    require_true(ownership.rank_offsets != NULL,
                 "failed exact rebalance consumed ownership");
  } else {
    require_true(ownership.rank_offsets == NULL,
                 "failed exact rebalance published ownership");
  }
  FreeSymmetryBasisOwnership(&ownership);

  memset(&stats, 0, sizeof(stats));
  require_true(
      SymmetryExactRebalanceBasisRun(
          &run, test_rank, test_nrank, &ownership, &stats) == 0,
      "valid exact rebalance after ownership failure failed");
  FreeSymmetryBasisOwnership(&ownership);
  FreeSymmetryBasisRun(&run);
}

static void assert_exact_rebalance_order_failure_recovery(void)
{
  struct SymmetryBasisRun run;
  struct SymmetryBasisOwnership ownership;
  struct SymmetryBasisDistributionStats stats;
  struct SymmetryBasisDistributionStats saved_stats;
  struct SymmetryBasisVector saved_entry;
  struct SymmetryBasisVector *saved_entries;
  uint64_t global_count = (uint64_t)test_nrank * 4U + 3U;
  uint64_t source_offset;
  uint64_t source_count;
  int failing_rank = test_nrank > 1 ? 1 : 0;
  int status;

  build_exact_range_run(global_count, FALSE, &run);
  exact_source_range(global_count, FALSE, test_rank,
                     &source_offset, &source_count);
  saved_entries = run.entries;
  memset(&ownership, 0, sizeof(ownership));
  memset(&stats, 0xa5, sizeof(stats));
  memcpy(&saved_stats, &stats, sizeof(saved_stats));
  if (test_rank == failing_rank) {
    require_true(source_count >= 2U,
                 "invalid exact order fixture is too short");
    if (test_nrank > 1) {
      require_true(source_offset > 0U,
                   "invalid exact boundary fixture lacks a predecessor");
      saved_entry = run.entries[1];
      run.entries[1] =
          make_sample_sort_entry(source_offset - 1U);
    } else {
      saved_entry = run.entries[2];
      run.entries[2] = run.entries[1];
    }
  }
  status = SymmetryExactRebalanceBasisRun(
      &run, test_rank, test_nrank, &ownership, &stats);
  require_all_ranks_failed(
      status, "globally non-strict exact input was not rejected");
  require_true(run.entries == saved_entries &&
                   ownership.rank_offsets == NULL &&
                   ownership.dim == 0UL &&
                   ownership.local_offset == 0UL &&
                   ownership.local_dim == 0UL &&
                   memcmp(&stats, &saved_stats, sizeof(stats)) == 0,
               "order failure changed exact ownership or stats");
  if (test_rank == failing_rank) {
    if (test_nrank > 1) {
      run.entries[1] = saved_entry;
    } else {
      run.entries[2] = saved_entry;
    }
  }

  memset(&stats, 0, sizeof(stats));
  require_true(
      SymmetryExactRebalanceBasisRun(
          &run, test_rank, test_nrank, &ownership, &stats) == 0,
      "valid exact rebalance after order failure failed");
  FreeSymmetryBasisOwnership(&ownership);
  FreeSymmetryBasisRun(&run);
}

static void assert_exact_rebalance_memory_cap_failure(void)
{
  struct SymmetryBasisRun run;
  struct SymmetryBasisOwnership ownership;
  struct SymmetryBasisDistributionStats stats;
  struct SymmetryBasisDistributionStats saved_stats;
  struct SymmetryBasisVector *saved_entries;
  unsigned long int saved_count;
  unsigned long int saved_capacity;
  uint64_t global_count = (uint64_t)test_nrank * UINT64_C(500);
  int status;

  build_exact_range_run(global_count, FALSE, &run);
  saved_entries = run.entries;
  saved_count = run.count;
  saved_capacity = run.capacity;
  memset(&ownership, 0, sizeof(ownership));
  memset(&stats, 0x3c, sizeof(stats));
  memcpy(&saved_stats, &stats, sizeof(saved_stats));

  status = SymmetryExactRebalanceBasisRun(
      &run, test_rank, test_nrank, &ownership, &stats);
  require_all_ranks_failed(
      status, "exact rebalance memory cap was not rejected collectively");
  require_true(run.entries == saved_entries &&
                   run.count == saved_count &&
                   run.capacity == saved_capacity &&
                   ownership.rank_offsets == NULL &&
                   ownership.dim == 0UL &&
                   ownership.local_offset == 0UL &&
                   ownership.local_dim == 0UL &&
                   memcmp(&stats, &saved_stats, sizeof(stats)) == 0,
               "memory-cap failure changed exact ownership or stats");
  FreeSymmetryBasisRun(&run);

  build_exact_range_run((uint64_t)test_nrank * 3U + 1U, FALSE, &run);
  memset(&stats, 0, sizeof(stats));
  require_true(
      SymmetryExactRebalanceBasisRun(
          &run, test_rank, test_nrank, &ownership, &stats) == 0,
      "valid exact rebalance after memory-cap failure failed");
  FreeSymmetryBasisOwnership(&ownership);
  FreeSymmetryBasisRun(&run);
}

int main(int argc, char **argv)
{
#ifdef MPI
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS ||
      MPI_Comm_rank(MPI_COMM_WORLD, &test_rank) != MPI_SUCCESS ||
      MPI_Comm_size(MPI_COMM_WORLD, &test_nrank) != MPI_SUCCESS) {
    fprintf(stderr, "MPI initialization failed\n");
    return 1;
  }
  test_mpi_active = TRUE;
#else
  if (argc == 2 && strcmp(argv[1], "--mpi") == 0) {
    fprintf(stderr, "MPI exchange test requested from an MPI-disabled build\n");
    return 1;
  }
#endif

  test_message_entry_limit =
      HPHI_SYMMETRY_EXCHANGE_MESSAGE_BYTES /
      sizeof(struct SymmetryBasisVector);
  test_request_message_entry_limit =
      HPHI_SYMMETRY_EXCHANGE_MESSAGE_BYTES /
      sizeof(unsigned long int);
  test_response_message_entry_limit =
      HPHI_SYMMETRY_EXCHANGE_MESSAGE_BYTES /
      sizeof(struct SymmetryMpiLookupResponse);
  require_true(test_message_entry_limit > 0U &&
                   test_message_entry_limit < (uint64_t)INT_MAX &&
                   test_request_message_entry_limit > 0U &&
                   test_request_message_entry_limit < (uint64_t)INT_MAX &&
                   test_response_message_entry_limit > 0U &&
                   test_response_message_entry_limit < (uint64_t)INT_MAX,
               "invalid test message byte cap");
  assert_sample_count_warning_policy();
  assert_asymmetric_fast_and_chunked();
  assert_byte_cap_boundary_and_deep_chunks();
  assert_empty_and_self_only();
  assert_collective_validation_failure();
  assert_chunk_limit_agreement();
  assert_overflow_and_option_validation();
  assert_typed_fast_and_chunked();
  assert_typed_deep_chunks();
  assert_typed_empty_and_self_only();
  assert_typed_failure_recovery();
  assert_typed_overflow_failure();
  assert_sample_sort_fixture(
      SAMPLE_SORT_RANDOM, "unique random sample sort mismatch");
  assert_sample_sort_fixture(
      SAMPLE_SORT_SKEW, "skew sample sort mismatch");
  assert_sample_sort_fixture(
      SAMPLE_SORT_ALREADY_SORTED, "already-sorted sample sort mismatch");
  assert_sample_sort_fixture(
      SAMPLE_SORT_REVERSE_SORTED, "reverse-sorted sample sort mismatch");
  assert_sample_sort_fixture(
      SAMPLE_SORT_ONE_RANK_ONLY, "one-rank-only sample sort mismatch");
  assert_sample_sort_fixture(
      SAMPLE_SORT_LOCAL_EMPTY, "local-empty sample sort mismatch");
  assert_sample_sort_fixture(
      SAMPLE_SORT_ALL_EMPTY, "all-empty sample sort mismatch");
  assert_sample_sort_duplicate_failure();
  assert_sample_sort_memory_cap_failure();
  assert_block_range_contract();
  assert_exact_rebalance_fixture(
      0U, FALSE, "zero-dimension exact rebalance mismatch");
  assert_exact_rebalance_fixture(
      1U, TRUE, "one-entry exact rebalance mismatch");
  if (test_nrank > 1) {
    assert_exact_rebalance_fixture(
        (uint64_t)test_nrank - 1U, FALSE,
        "dimension-below-rank-count exact rebalance mismatch");
  }
  assert_exact_rebalance_fixture(
      (uint64_t)test_nrank * 3U, FALSE,
      "divisible exact rebalance mismatch");
  assert_exact_rebalance_fixture(
      (uint64_t)test_nrank * 3U + 1U, FALSE,
      "remainder-one exact rebalance mismatch");
  if (test_nrank > 2) {
    assert_exact_rebalance_fixture(
        (uint64_t)test_nrank * 3U +
            (uint64_t)test_nrank / 2U,
        FALSE, "intermediate-remainder exact rebalance mismatch");
  }
  assert_exact_rebalance_fixture(
      (uint64_t)test_nrank * 8U + 3U, TRUE,
      "one-range-rank exact rebalance mismatch");
  assert_exact_rebalance_failure_recovery();
  assert_exact_rebalance_order_failure_recovery();
  assert_exact_rebalance_memory_cap_failure();

  if (test_rank == 0) {
    fprintf(stdout,
            "checked symmetry basis exchange gate: PASS (%d rank%s)\n",
            test_nrank, test_nrank == 1 ? "" : "s");
  }
#ifdef MPI
  if (MPI_Finalize() != MPI_SUCCESS) return 1;
#endif
  return 0;
}
