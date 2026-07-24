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
#include "symmetry_mpi_exchange.h"

static int test_rank = 0;
static int test_nrank = 1;
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

static void validate_asymmetric_result(
    const struct SymmetryMpiExchangeResult *result)
{
  uint64_t total = 0U;
  int source;
  require_true(result->nrank == test_nrank,
               "receive rank count mismatch");
  for (source = 0; source < test_nrank; source++) {
    uint64_t expected_count = asymmetric_count(source, test_rank);
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
  fast_options.chunk_limit = 2U;
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

  free(counts);
  free(displacements);
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

  assert_asymmetric_fast_and_chunked();
  assert_empty_and_self_only();
  assert_collective_validation_failure();
  assert_overflow_and_option_validation();

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
