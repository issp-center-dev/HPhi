#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "symmetry_basis.h"
#include "symmetry_directory.h"

#ifdef MPI
#include <mpi.h>
#endif

static int test_rank = 0;
static int test_nrank = 1;
#ifdef MPI
static int test_mpi_active = 0;
#endif

static void fail_test(const char *label)
{
  fprintf(stderr, "rank %d: %s\n", test_rank, label);
#ifdef MPI
  if (test_mpi_active != 0) MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  exit(1);
}

static void require_true(int condition, const char *label)
{
  if (condition == 0) fail_test(label);
}

static struct SymmetryBasisVector *make_basis(unsigned long int count)
{
  struct SymmetryBasisVector *basis;
  size_t elements;
  if (count == ULONG_MAX ||
      count + 1UL > (unsigned long int)(SIZE_MAX /
                                       sizeof(struct SymmetryBasisVector))) {
    return NULL;
  }
  elements = (size_t)(count + 1UL);
  basis = (struct SymmetryBasisVector *)calloc(elements, sizeof(*basis));
  return basis;
}

static void set_entry(struct SymmetryBasisVector *entry,
                      unsigned long int rep_state,
                      double norm)
{
  memset(entry, 0, sizeof(*entry));
  entry->rep_state = rep_state;
  entry->orbit_size = 1U;
  entry->stabilizer_size = 1U;
  entry->norm = norm;
  entry->stabilizer_character_sum = 1.0;
}

static void assert_basic_lookup(void)
{
  const unsigned long int keys[] = {0UL, 7UL, 19UL, 50UL, 91UL};
  const double norms[] = {0.5, 1.25, 2.5, 4.0, 8.5};
  struct SymmetryBasisVector *basis = make_basis(5UL);
  struct SymmetryLocalRepresentativeIndex *index = NULL;
  struct SymmetryLocalRepresentativeIndexStats stats;
  unsigned long int entry;
  require_true(basis != NULL, "basic basis allocation failed");
  for (entry = 1UL; entry <= 5UL; entry++) {
    set_entry(&basis[entry], keys[entry - 1UL], norms[entry - 1UL]);
  }
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 5UL, 6UL, 100UL, &index) == 0 &&
                   index != NULL,
               "basic local index build failed");
  require_true(GetSymmetryLocalRepresentativeIndexStats(index, &stats) == 0,
               "basic stats query failed");
  require_true(stats.table_size >= 10UL &&
                   (stats.table_size & (stats.table_size - 1UL)) == 0UL &&
                   stats.table_bytes ==
                       (size_t)stats.table_size *
                           2U * sizeof(unsigned long int) &&
                   stats.build_probes ==
                       UINT64_C(5) + stats.build_collisions &&
                   stats.build_max_probe >= UINT64_C(1),
               "basic construction stats mismatch");
  for (entry = 1UL; entry <= 5UL; entry++) {
    unsigned long int local_index = ULONG_MAX;
    unsigned long int global_beta = ULONG_MAX;
    double norm = -1.0;
    uint64_t probes = UINT64_MAX;
    require_true(SymmetryLookupLocalRepresentative(
                     index, keys[entry - 1UL], &local_index,
                     &global_beta, &norm, &probes) == 0,
                 "basic found lookup failed");
    require_true(local_index == entry &&
                     global_beta == 100UL + entry &&
                     memcmp(&norm, &norms[entry - 1UL],
                            sizeof(norm)) == 0 &&
                     probes >= UINT64_C(1),
                 "basic found lookup mismatch");
  }
  {
    unsigned long int local_index = ULONG_MAX;
    unsigned long int global_beta = ULONG_MAX;
    double norm = -1.0;
    uint64_t probes = UINT64_MAX;
    require_true(SymmetryLookupLocalRepresentative(
                     index, 49UL, &local_index, &global_beta,
                     &norm, &probes) == 0,
                 "basic miss lookup failed");
    require_true(local_index == 0UL && global_beta == 0UL &&
                     norm == 0.0 && probes >= UINT64_C(1),
                 "basic miss result mismatch");
  }
  FreeSymmetryLocalRepresentativeIndex(index);
  free(basis);
}

static void assert_empty_index(void)
{
  struct SymmetryBasisVector *allocated_empty = NULL;
  struct SymmetryLocalRepresentativeIndex *index = NULL;
  struct SymmetryLocalRepresentativeIndexStats stats;
  unsigned long int local_index = ULONG_MAX;
  unsigned long int global_beta = ULONG_MAX;
  double norm = -1.0;
  uint64_t probes = UINT64_MAX;
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   NULL, 0UL, 0UL, ULONG_MAX, &index) == 0 &&
                   index != NULL,
               "empty local index build failed");
  require_true(GetSymmetryLocalRepresentativeIndexStats(index, &stats) == 0 &&
                   stats.table_size == 0UL &&
                   stats.table_bytes == 0U &&
                   stats.build_collisions == 0U &&
                   stats.build_probes == 0U &&
                   stats.build_max_probe == 0U,
               "empty construction stats mismatch");
  require_true(SymmetryLookupLocalRepresentative(
                   index, 0UL, &local_index, &global_beta,
                   &norm, &probes) == 0 &&
                   local_index == 0UL && global_beta == 0UL &&
                   norm == 0.0 && probes == 0U,
               "empty lookup mismatch");
  FreeSymmetryLocalRepresentativeIndex(index);
  index = NULL;

  allocated_empty = make_basis(0UL);
  require_true(allocated_empty != NULL,
               "allocated empty basis allocation failed");
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   allocated_empty, 0UL, 1UL, 0UL, &index) == 0 &&
                   index != NULL,
               "allocated empty local index build failed");
  FreeSymmetryLocalRepresentativeIndex(index);
  free(allocated_empty);
}

static void assert_collision_accounting(void)
{
  struct SymmetryBasisVector *basis = make_basis(2UL);
  unsigned long int first;
  int found_collision = 0;
  require_true(basis != NULL, "collision basis allocation failed");
  for (first = 0UL; first < 64UL && found_collision == 0; first++) {
    unsigned long int second;
    for (second = first + 1UL; second < 64UL; second++) {
      struct SymmetryLocalRepresentativeIndex *index = NULL;
      struct SymmetryLocalRepresentativeIndexStats stats;
      set_entry(&basis[1], first, 1.0);
      set_entry(&basis[2], second, 2.0);
      require_true(BuildSymmetryLocalRepresentativeIndex(
                       basis, 2UL, 3UL, 0UL, &index) == 0,
                   "collision candidate build failed");
      require_true(GetSymmetryLocalRepresentativeIndexStats(
                       index, &stats) == 0,
                   "collision candidate stats failed");
      if (stats.build_collisions > 0U) {
        unsigned long int local_index = 0UL;
        unsigned long int global_beta = 0UL;
        double norm = 0.0;
        uint64_t probes = 0U;
        require_true(stats.build_probes ==
                         UINT64_C(2) + stats.build_collisions &&
                         stats.build_max_probe > UINT64_C(1),
                     "collision construction accounting mismatch");
        require_true(SymmetryLookupLocalRepresentative(
                         index, second, &local_index, &global_beta,
                         &norm, &probes) == 0 &&
                         local_index == 2UL && global_beta == 2UL &&
                         norm == 2.0 && probes > UINT64_C(1),
                     "collision lookup mismatch");
        found_collision = 1;
      }
      FreeSymmetryLocalRepresentativeIndex(index);
    }
  }
  require_true(found_collision != 0,
               "failed to construct a collision fixture");
  free(basis);
}

static void assert_invalid_builds_and_recovery(void)
{
  struct SymmetryBasisVector *basis = make_basis(3UL);
  struct SymmetryLocalRepresentativeIndex *index = NULL;
  unsigned long int huge_dim;
  require_true(basis != NULL, "invalid basis allocation failed");
  set_entry(&basis[1], 1UL, 1.0);
  set_entry(&basis[2], 2UL, 2.0);
  set_entry(&basis[3], 3UL, 3.0);

  require_true(BuildSymmetryLocalRepresentativeIndex(
                   NULL, 1UL, 2UL, 0UL, &index) != 0 &&
                   index == NULL,
               "null nonempty basis accepted");
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 3UL, 3UL, 0UL, &index) != 0 &&
                   index == NULL,
               "short capacity accepted");
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 3UL, 4UL, ULONG_MAX - 2UL, &index) != 0 &&
                   index == NULL,
               "global beta overflow accepted");
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   NULL, 0UL, 1UL, 0UL, &index) != 0 &&
                   index == NULL,
               "inconsistent empty capacity accepted");
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 0UL, 0UL, 0UL, &index) != 0 &&
                   index == NULL,
               "inconsistent empty storage accepted");

  set_entry(&basis[2], 1UL, 2.0);
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 3UL, 4UL, 0UL, &index) != 0 &&
                   index == NULL,
               "duplicate key accepted");
  set_entry(&basis[2], 0UL, 2.0);
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 3UL, 4UL, 0UL, &index) != 0 &&
                   index == NULL,
               "descending key accepted");
  set_entry(&basis[2], 2UL, NAN);
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 3UL, 4UL, 0UL, &index) != 0 &&
                   index == NULL,
               "NaN norm accepted");
  set_entry(&basis[2], 2UL, INFINITY);
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 3UL, 4UL, 0UL, &index) != 0 &&
                   index == NULL,
               "infinite norm accepted");
  set_entry(&basis[2], 2UL, 0.0);
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 3UL, 4UL, 0UL, &index) != 0 &&
                   index == NULL,
               "zero norm accepted");
  set_entry(&basis[2], 2UL, -1.0);
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 3UL, 4UL, 0UL, &index) != 0 &&
                   index == NULL,
               "negative norm accepted");
  set_entry(&basis[2], 2UL, 2.0);

  huge_dim = (ULONG_MAX - 1UL) / 2UL + 1UL;
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, huge_dim, huge_dim + 1UL, 0UL, &index) != 0 &&
                   index == NULL,
               "table shape overflow accepted");
  huge_dim =
      (unsigned long int)(SIZE_MAX / sizeof(unsigned long int) / 2U);
  if (huge_dim > 0UL && huge_dim < (ULONG_MAX - 1UL) / 2UL) {
    require_true(BuildSymmetryLocalRepresentativeIndex(
                     basis, huge_dim, huge_dim + 1UL, 0UL, &index) != 0 &&
                     index == NULL,
                 "table byte overflow accepted");
  }

  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 3UL, 4UL, 9UL, &index) == 0 &&
                   index != NULL,
               "valid build after failures failed");
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 3UL, 4UL, 9UL, &index) != 0 &&
                   index != NULL,
               "nonempty output object accepted");
  FreeSymmetryLocalRepresentativeIndex(index);
  free(basis);
}

static void assert_lookup_failure_atomicity(void)
{
  struct SymmetryBasisVector *basis = make_basis(1UL);
  struct SymmetryLocalRepresentativeIndex *index = NULL;
  struct SymmetryLocalRepresentativeIndexStats stats;
  unsigned long int local_index = 71UL;
  unsigned long int global_beta = 72UL;
  double norm = -73.0;
  uint64_t probes = UINT64_C(74);
  require_true(basis != NULL, "atomicity basis allocation failed");
  set_entry(&basis[1], 0UL, 1.5);
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 1UL, 2UL, 4UL, &index) == 0,
               "atomicity index build failed");
  require_true(SymmetryLookupLocalRepresentative(
                   NULL, 0UL, &local_index, &global_beta,
                   &norm, &probes) != 0 &&
                   local_index == 71UL && global_beta == 72UL &&
                   norm == -73.0 && probes == UINT64_C(74),
               "null lookup changed outputs");
  require_true(SymmetryLookupLocalRepresentative(
                   index, 0UL, &local_index, NULL,
                   &norm, &probes) != 0 &&
                   local_index == 71UL && norm == -73.0 &&
                   probes == UINT64_C(74),
               "invalid output lookup changed outputs");
  memset(&stats, 0xa5, sizeof(stats));
  {
    struct SymmetryLocalRepresentativeIndexStats saved = stats;
    require_true(GetSymmetryLocalRepresentativeIndexStats(
                     NULL, &stats) != 0 &&
                     memcmp(&stats, &saved, sizeof(stats)) == 0,
                 "invalid stats query changed output");
  }
  require_true(SymmetryLookupLocalRepresentative(
                   index, 0UL, NULL, &global_beta,
                   &norm, NULL) == 0 &&
                   global_beta == 5UL && norm == 1.5,
               "optional lookup outputs failed");
  FreeSymmetryLocalRepresentativeIndex(index);
  free(basis);
}

static void assert_parallel_read_lookup(void)
{
  struct SymmetryBasisVector *basis = make_basis(32UL);
  struct SymmetryLocalRepresentativeIndex *index = NULL;
  unsigned long int entry;
  int failed = 0;
  require_true(basis != NULL, "parallel basis allocation failed");
  for (entry = 1UL; entry <= 32UL; entry++) {
    set_entry(&basis[entry], entry * 3UL, 0.25 * (double)entry);
  }
  require_true(BuildSymmetryLocalRepresentativeIndex(
                   basis, 32UL, 33UL, 500UL, &index) == 0,
               "parallel index build failed");
#ifdef _OPENMP
#pragma omp parallel for reduction(| : failed)
#endif
  for (entry = 1UL; entry <= 4096UL; entry++) {
    unsigned long int expected = (entry - 1UL) % 32UL + 1UL;
    unsigned long int local_index = 0UL;
    unsigned long int global_beta = 0UL;
    double norm = 0.0;
    uint64_t probes = 0U;
    if (SymmetryLookupLocalRepresentative(
            index, expected * 3UL, &local_index, &global_beta,
            &norm, &probes) != 0 ||
        local_index != expected ||
        global_beta != 500UL + expected ||
        norm != 0.25 * (double)expected ||
        probes == 0U) {
      failed = 1;
    }
  }
  require_true(failed == 0, "parallel read-only lookup mismatch");
  FreeSymmetryLocalRepresentativeIndex(index);
  free(basis);
}

struct DirectoryFixture {
  struct SymmetryBasisVector *basis;
  unsigned long int *rank_offsets;
  unsigned long int dim;
  unsigned long int local_offset;
  unsigned long int local_dim;
  unsigned long int local_capacity;
  unsigned long int key_base;
};

static int test_block_range(unsigned long int dim,
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

static unsigned long int fixture_key(unsigned long int key_base,
                                     unsigned long int ordinal)
{
  require_true(ordinal <= (ULONG_MAX - key_base) / 10UL,
               "fixture key overflow");
  return key_base + ordinal * 10UL;
}

static void init_directory_fixture(unsigned long int dim,
                                   unsigned long int key_base,
                                   struct DirectoryFixture *fixture)
{
  int peer;
  memset(fixture, 0, sizeof(*fixture));
  fixture->dim = dim;
  fixture->key_base = key_base;
  fixture->rank_offsets = (unsigned long int *)calloc(
      (size_t)test_nrank + 1U, sizeof(*fixture->rank_offsets));
  require_true(fixture->rank_offsets != NULL,
               "directory rank offset allocation failed");
  for (peer = 0; peer < test_nrank; peer++) {
    unsigned long int offset;
    unsigned long int count;
    require_true(test_block_range(dim, peer, test_nrank,
                                  &offset, &count) == 0,
                 "directory fixture block range failed");
    fixture->rank_offsets[peer] = offset;
    fixture->rank_offsets[peer + 1] = offset + count;
    if (peer == test_rank) {
      fixture->local_offset = offset;
      fixture->local_dim = count;
    }
  }
  if (fixture->local_dim == 0UL && (test_rank % 2) == 0) {
    fixture->basis = NULL;
    fixture->local_capacity = 0UL;
  } else {
    unsigned long int local_index;
    fixture->basis = make_basis(fixture->local_dim);
    require_true(fixture->basis != NULL,
                 "directory local basis allocation failed");
    fixture->local_capacity = fixture->local_dim + 1UL;
    for (local_index = 1UL;
         local_index <= fixture->local_dim;
         local_index++) {
      unsigned long int ordinal =
          fixture->local_offset + local_index - 1UL;
      set_entry(&fixture->basis[local_index],
                fixture_key(key_base, ordinal),
                0.5 + (double)ordinal);
    }
  }
}

static void free_directory_fixture(struct DirectoryFixture *fixture)
{
  if (fixture == NULL) return;
  free(fixture->basis);
  free(fixture->rank_offsets);
  memset(fixture, 0, sizeof(*fixture));
}

static int expected_owner_for_ordinal(
    const struct DirectoryFixture *fixture,
    unsigned long int ordinal)
{
  int peer;
  for (peer = 0; peer < test_nrank; peer++) {
    if (ordinal < fixture->rank_offsets[peer + 1]) return peer;
  }
  return -1;
}

static int expected_nonempty_rank_count(unsigned long int dim)
{
  if (dim < (unsigned long int)test_nrank) return (int)dim;
  return test_nrank;
}

static void assert_directory_fixture(unsigned long int dim,
                                     unsigned long int key_base,
                                     const char *label)
{
  struct DirectoryFixture fixture;
  struct SymmetryRepresentativeDirectory *directory = NULL;
  struct SymmetryRepresentativeDirectoryInfo info;
  struct SymmetryLocalRepresentativeIndexStats index_stats;
  unsigned long int ordinal;
  int nonempty_rank_count = expected_nonempty_rank_count(dim);
  int owner = 77;
  init_directory_fixture(dim, key_base, &fixture);
  require_true(BuildSymmetryRepresentativeDirectory(
                   fixture.basis, fixture.dim, fixture.local_dim,
                   fixture.local_capacity, fixture.local_offset,
                   fixture.rank_offsets, test_rank, test_nrank,
                   &directory) == 0 &&
                   SymmetryRepresentativeDirectoryReady(directory) != 0,
               label);
  require_true(GetSymmetryRepresentativeDirectoryInfo(
                   directory, &info) == 0 &&
                   info.dim == fixture.dim &&
                   info.local_offset == fixture.local_offset &&
                   info.local_dim == fixture.local_dim &&
                   info.rank == test_rank &&
                   info.nrank == test_nrank &&
                   info.nonempty_rank_count == nonempty_rank_count &&
                   info.splitter_bytes ==
                       (size_t)nonempty_rank_count *
                           sizeof(unsigned long int),
               "directory info mismatch");
  require_true(GetSymmetryRepresentativeDirectoryLocalIndexStats(
                   directory, &index_stats) == 0 &&
                   (fixture.local_dim == 0UL
                        ? index_stats.table_size == 0UL
                        : index_stats.table_size >=
                              fixture.local_dim * 2UL),
               "directory local index stats mismatch");

  if (dim == 0UL) {
    unsigned long int local_index = ULONG_MAX;
    unsigned long int global_beta = ULONG_MAX;
    double norm = -1.0;
    uint64_t probes = UINT64_MAX;
    require_true(SymmetryRepresentativeOwner(
                     directory, 0UL, &owner) == 0 &&
                     owner == -1,
                 "zero-dimensional owner mismatch");
    require_true(SymmetryLookupDirectoryLocalRepresentative(
                     directory, 0UL, &local_index, &global_beta,
                     &norm, &probes) == 0 &&
                     local_index == 0UL && global_beta == 0UL &&
                     norm == 0.0 && probes == 0U,
                 "zero-dimensional local lookup mismatch");
  } else {
    for (ordinal = 0UL; ordinal < dim; ordinal++) {
      unsigned long int local_index = ULONG_MAX;
      unsigned long int global_beta = ULONG_MAX;
      double norm = -1.0;
      uint64_t probes = UINT64_MAX;
      unsigned long int key = fixture_key(key_base, ordinal);
      int expected_owner =
          expected_owner_for_ordinal(&fixture, ordinal);
      require_true(SymmetryRepresentativeOwner(
                       directory, key, &owner) == 0 &&
                       owner == expected_owner,
                   "exact representative owner mismatch");
      require_true(SymmetryLookupDirectoryLocalRepresentative(
                       directory, key, &local_index, &global_beta,
                       &norm, &probes) == 0,
                   "directory local lookup failed");
      if (expected_owner == test_rank) {
        unsigned long int expected_local =
            ordinal - fixture.local_offset + 1UL;
        double expected_norm = 0.5 + (double)ordinal;
        require_true(local_index == expected_local &&
                         global_beta == ordinal + 1UL &&
                         memcmp(&norm, &expected_norm,
                                sizeof(norm)) == 0 &&
                         probes >= UINT64_C(1),
                     "directory local found mismatch");
      } else {
        require_true(local_index == 0UL && global_beta == 0UL &&
                         norm == 0.0 &&
                         (fixture.local_dim == 0UL
                              ? probes == 0U
                              : probes >= UINT64_C(1)),
                     "directory remote key was found locally");
      }
    }
    if (key_base > 0UL) {
      require_true(SymmetryRepresentativeOwner(
                       directory, key_base - 1UL, &owner) == 0 &&
                       owner == 0,
                   "owner below global first mismatch");
    }
    {
      int peer;
      for (peer = 0; peer + 1 < nonempty_rank_count; peer++) {
        unsigned long int previous_ordinal =
            fixture.rank_offsets[peer + 1] - 1UL;
        unsigned long int gap_key =
            fixture_key(key_base, previous_ordinal) + 1UL;
        require_true(SymmetryRepresentativeOwner(
                         directory, gap_key, &owner) == 0 &&
                         owner == peer,
                     "rank-boundary gap owner mismatch");
      }
    }
    require_true(SymmetryRepresentativeOwner(
                     directory,
                     fixture_key(key_base, dim - 1UL) + 1UL,
                     &owner) == 0 &&
                     owner == nonempty_rank_count - 1,
                 "owner above global last mismatch");
    require_true(SymmetryRepresentativeOwner(
                     directory, ULONG_MAX, &owner) == 0 &&
                     owner == nonempty_rank_count - 1,
                 "ULONG_MAX owner mismatch");
  }

  owner = 71;
  require_true(SymmetryRepresentativeOwner(
                   NULL, 0UL, &owner) != 0 &&
                   owner == 71,
               "invalid owner query changed output");
  require_true(SymmetryRepresentativeOwner(
                   directory, 0UL, NULL) != 0,
               "null owner output accepted");
  FreeSymmetryRepresentativeDirectory(directory);
  free_directory_fixture(&fixture);
}

enum DirectoryFailureKind {
  DIRECTORY_FAILURE_DUPLICATE = 0,
  DIRECTORY_FAILURE_DESCENDING = 1,
  DIRECTORY_FAILURE_NORM = 2,
  DIRECTORY_FAILURE_CAPACITY = 3,
  DIRECTORY_FAILURE_OFFSET = 4,
  DIRECTORY_FAILURE_COUNT = 5,
  DIRECTORY_FAILURE_RANK_OFFSETS = 6,
  DIRECTORY_FAILURE_BOUNDARY = 7
};

static void assert_directory_valid_recovery(unsigned long int dim)
{
  struct DirectoryFixture fixture;
  struct SymmetryRepresentativeDirectory *directory = NULL;
  init_directory_fixture(dim, 100UL, &fixture);
  require_true(BuildSymmetryRepresentativeDirectory(
                   fixture.basis, fixture.dim, fixture.local_dim,
                   fixture.local_capacity, fixture.local_offset,
                   fixture.rank_offsets, test_rank, test_nrank,
                   &directory) == 0 &&
                   SymmetryRepresentativeDirectoryReady(directory) != 0,
               "valid directory recovery failed");
  FreeSymmetryRepresentativeDirectory(directory);
  free_directory_fixture(&fixture);
}

static void assert_directory_failure_kind(
    enum DirectoryFailureKind kind,
    unsigned long int dim)
{
  struct DirectoryFixture fixture;
  struct SymmetryRepresentativeDirectory *directory = NULL;
  unsigned long int local_dim;
  unsigned long int local_capacity;
  unsigned long int local_offset;
  int failing_rank = test_nrank > 1 ? 1 : 0;
  int status;
  init_directory_fixture(dim, 100UL, &fixture);
  local_dim = fixture.local_dim;
  local_capacity = fixture.local_capacity;
  local_offset = fixture.local_offset;
  if (test_rank == failing_rank) {
    switch (kind) {
    case DIRECTORY_FAILURE_DUPLICATE:
      fixture.basis[2].rep_state = fixture.basis[1].rep_state;
      break;
    case DIRECTORY_FAILURE_DESCENDING:
      fixture.basis[2].rep_state =
          fixture.basis[1].rep_state - 1UL;
      break;
    case DIRECTORY_FAILURE_NORM:
      fixture.basis[2].norm = 0.0;
      break;
    case DIRECTORY_FAILURE_CAPACITY:
      local_capacity = fixture.local_dim;
      break;
    case DIRECTORY_FAILURE_OFFSET:
      local_offset++;
      break;
    case DIRECTORY_FAILURE_COUNT:
      local_dim--;
      break;
    case DIRECTORY_FAILURE_RANK_OFFSETS:
      fixture.rank_offsets[0] = 1UL;
      break;
    case DIRECTORY_FAILURE_BOUNDARY:
      {
        unsigned long int local_index;
        unsigned long int previous_last_ordinal =
            fixture.rank_offsets[failing_rank] - 1UL;
        unsigned long int previous_last =
            fixture_key(fixture.key_base, previous_last_ordinal);
        for (local_index = 1UL;
             local_index <= fixture.local_dim;
             local_index++) {
          fixture.basis[local_index].rep_state =
              previous_last + local_index - 1UL;
        }
      }
      break;
    }
  }
  status = BuildSymmetryRepresentativeDirectory(
      fixture.basis, fixture.dim, local_dim, local_capacity,
      local_offset, fixture.rank_offsets, test_rank, test_nrank,
      &directory);
  require_true(status != 0 && directory == NULL,
               "invalid directory input was not rejected atomically");
  free_directory_fixture(&fixture);
  assert_directory_valid_recovery(dim);
}

static void assert_directory_failure_recovery(void)
{
  unsigned long int dim = (unsigned long int)test_nrank * 4UL;
  enum DirectoryFailureKind kind;
  for (kind = DIRECTORY_FAILURE_DUPLICATE;
       kind <= DIRECTORY_FAILURE_RANK_OFFSETS;
       kind = (enum DirectoryFailureKind)((int)kind + 1)) {
    assert_directory_failure_kind(kind, dim);
  }
  if (test_nrank > 1) {
    assert_directory_failure_kind(DIRECTORY_FAILURE_BOUNDARY, dim);
  }
  {
    struct DirectoryFixture fixture;
    struct SymmetryRepresentativeDirectory *directory = NULL;
    struct SymmetryRepresentativeDirectory *existing = NULL;
    int failing_rank = test_nrank > 1 ? 1 : 0;
    int status;
    init_directory_fixture(dim, 100UL, &fixture);
    require_true(BuildSymmetryRepresentativeDirectory(
                     fixture.basis, fixture.dim, fixture.local_dim,
                     fixture.local_capacity, fixture.local_offset,
                     fixture.rank_offsets, test_rank, test_nrank,
                     &directory) == 0,
                 "existing-directory fixture build failed");
    if (test_rank == failing_rank) {
      existing = directory;
    } else {
      FreeSymmetryRepresentativeDirectory(directory);
    }
    directory = existing;
    status = BuildSymmetryRepresentativeDirectory(
        fixture.basis, fixture.dim, fixture.local_dim,
        fixture.local_capacity, fixture.local_offset,
        fixture.rank_offsets, test_rank, test_nrank,
        &directory);
    require_true(status != 0 &&
                     (test_rank == failing_rank
                          ? directory == existing
                          : directory == NULL),
                 "nonempty directory output was not rejected atomically");
    if (test_rank == failing_rank) {
      FreeSymmetryRepresentativeDirectory(existing);
    }
    free_directory_fixture(&fixture);
    assert_directory_valid_recovery(dim);
  }
}

int main(int argc, char **argv)
{
  int mpi_requested =
      argc == 2 && strcmp(argv[1], "--mpi") == 0;
  if (argc > 2 || (argc == 2 && mpi_requested == 0)) {
    fprintf(stderr, "Usage: %s [--mpi]\n", argv[0]);
    return 1;
  }
#ifdef MPI
  if (mpi_requested != 0) {
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS ||
        MPI_Comm_rank(MPI_COMM_WORLD, &test_rank) != MPI_SUCCESS ||
        MPI_Comm_size(MPI_COMM_WORLD, &test_nrank) != MPI_SUCCESS) {
      fprintf(stderr, "MPI initialization failed\n");
      return 1;
    }
    test_mpi_active = 1;
  }
#else
  if (mpi_requested != 0) {
    fprintf(stderr, "MPI directory test requested from an MPI-disabled build\n");
    return 1;
  }
#endif
  assert_basic_lookup();
  assert_empty_index();
  assert_collision_accounting();
  assert_invalid_builds_and_recovery();
  assert_lookup_failure_atomicity();
  assert_parallel_read_lookup();
  assert_directory_fixture(
      0UL, 0UL, "zero-dimensional directory build failed");
  assert_directory_fixture(
      1UL, 0UL, "one-entry key-zero directory build failed");
  if (test_nrank > 1) {
    assert_directory_fixture(
        (unsigned long int)test_nrank - 1UL, 10UL,
        "dimension-below-rank-count directory build failed");
  }
  assert_directory_fixture(
      (unsigned long int)test_nrank * 2UL, 10UL,
      "divisible directory build failed");
  assert_directory_fixture(
      (unsigned long int)test_nrank * 2UL + 1UL, 10UL,
      "remainder-one directory build failed");
  if (test_nrank > 2) {
    assert_directory_fixture(
        (unsigned long int)test_nrank * 2UL +
            (unsigned long int)test_nrank - 1UL,
        10UL, "intermediate-remainder directory build failed");
  }
  assert_directory_failure_recovery();
  if (test_rank == 0) {
    printf("symmetry directory metadata gate: PASS (%d rank%s)\n",
           test_nrank, test_nrank == 1 ? "" : "s");
  }
#ifdef MPI
  if (test_mpi_active != 0 &&
      MPI_Finalize() != MPI_SUCCESS) {
    return 1;
  }
#endif
  return 0;
}
