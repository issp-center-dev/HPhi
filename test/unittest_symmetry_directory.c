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

static void init_explicit_directory_fixture(
    const unsigned long int *global_keys,
    unsigned long int dim,
    struct DirectoryFixture *fixture)
{
  int peer;
  memset(fixture, 0, sizeof(*fixture));
  fixture->dim = dim;
  fixture->rank_offsets = (unsigned long int *)calloc(
      (size_t)test_nrank + 1U, sizeof(*fixture->rank_offsets));
  require_true(global_keys != NULL && dim > 0UL &&
                   fixture->rank_offsets != NULL,
               "explicit directory fixture allocation failed");
  for (peer = 0; peer < test_nrank; peer++) {
    unsigned long int offset;
    unsigned long int count;
    require_true(test_block_range(dim, peer, test_nrank,
                                  &offset, &count) == 0,
                 "explicit directory block range failed");
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
                 "explicit local basis allocation failed");
    fixture->local_capacity = fixture->local_dim + 1UL;
    for (local_index = 1UL;
         local_index <= fixture->local_dim;
         local_index++) {
      unsigned long int ordinal =
          fixture->local_offset + local_index - 1UL;
      set_entry(&fixture->basis[local_index], global_keys[ordinal],
                1.0 + 0.125 * (double)ordinal);
    }
  }
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

static struct SymmetryRepresentativeDirectory *build_batch_directory(
    struct DirectoryFixture *fixture,
    unsigned long int dim,
    unsigned long int key_base,
    const char *label)
{
  struct SymmetryRepresentativeDirectory *directory = NULL;
  init_directory_fixture(dim, key_base, fixture);
  require_true(BuildSymmetryRepresentativeDirectory(
                   fixture->basis, fixture->dim, fixture->local_dim,
                   fixture->local_capacity, fixture->local_offset,
                   fixture->rank_offsets, test_rank, test_nrank,
                   &directory) == 0 &&
                   SymmetryRepresentativeDirectoryReady(directory) != 0,
               label);
  return directory;
}

static void assert_batch_result(
    const struct DirectoryFixture *fixture,
    const unsigned long int *keys,
    uint64_t count,
    const unsigned long int *global_beta,
    const double *norm,
    const char *label)
{
  uint64_t index;
  for (index = 0U; index < count; index++) {
    unsigned long int expected_beta = 0UL;
    double expected_norm = 0.0;
    if (keys[index] >= fixture->key_base) {
      unsigned long int delta = keys[index] - fixture->key_base;
      if (delta % 10UL == 0UL &&
          delta / 10UL < fixture->dim) {
        unsigned long int ordinal = delta / 10UL;
        expected_beta = ordinal + 1UL;
        expected_norm = 0.5 + (double)ordinal;
      }
    }
    if (global_beta[index] != expected_beta ||
        memcmp(&norm[index], &expected_norm, sizeof(expected_norm)) != 0 ||
        (expected_beta == 0UL &&
         (norm[index] != 0.0 || signbit(norm[index])))) {
      fail_test(label);
    }
  }
}

static void assert_batch_resolution(void)
{
  struct DirectoryFixture fixture;
  struct SymmetryRepresentativeDirectory *directory;
  struct SymmetryRepresentativeBatchOptions options;
  struct SymmetryRepresentativeBatchStats first_stats;
  struct SymmetryRepresentativeBatchStats second_stats;
  struct SymmetryRepresentativeBatchStats third_stats;
  unsigned long int dim = (unsigned long int)test_nrank * 4UL + 1UL;
  unsigned long int key_base = 100UL;
  uint64_t count = (uint64_t)dim * UINT64_C(2) + UINT64_C(1);
  unsigned long int *keys =
      (unsigned long int *)malloc((size_t)count * sizeof(*keys));
  unsigned long int *global_beta =
      (unsigned long int *)malloc((size_t)count * sizeof(*global_beta));
  double *norm = (double *)malloc((size_t)count * sizeof(*norm));
  uint64_t index = 0U;
  unsigned long int ordinal;
  require_true(keys != NULL && global_beta != NULL && norm != NULL,
               "batch fixture allocation failed");
  directory = build_batch_directory(
      &fixture, dim, key_base, "batch directory build failed");

  keys[index++] = key_base - 1UL;
  for (ordinal = 0UL; ordinal < dim; ordinal++) {
    keys[index++] = fixture_key(key_base, ordinal);
    if (ordinal + 1UL < dim) {
      keys[index++] = fixture_key(key_base, ordinal) + 1UL;
    }
  }
  keys[index++] = fixture_key(key_base, dim - 1UL) + 1UL;
  require_true(index == count, "batch request fixture count mismatch");
  for (index = 0U; index < count; index++) {
    global_beta[index] = ULONG_MAX;
    norm[index] = -1.0;
  }
  require_true(SymmetryResolveRepresentativeBatch(
                   directory, keys, count, global_beta, norm) == 0,
               "batch fast resolution failed");
  assert_batch_result(
      &fixture, keys, count, global_beta, norm,
      "batch fast result mismatch");
  require_true(GetSymmetryRepresentativeDirectoryBatchStats(
                   directory, &first_stats) == 0 &&
                   first_stats.directory_batch_calls == UINT64_C(1) &&
                   first_stats.directory_request_entries_sent == count &&
                   first_stats.directory_request_entries_received > 0U &&
                   first_stats.directory_found_entries > 0U &&
                   first_stats.directory_not_found_entries > 0U &&
                   first_stats.directory_lookup_probe_count > 0U &&
                   first_stats.directory_lookup_max_probe > 0U &&
                   first_stats.directory_owner_peer_count_max > 0U &&
                   first_stats.directory_requester_peer_count_max > 0U &&
                   first_stats.directory_exchange_message_byte_limit > 0U &&
                   first_stats.directory_batch_temporary_peak_bytes > 0U &&
                   first_stats.directory_batch_temporary_peak_bytes <=
                       first_stats.directory_batch_memory_byte_limit &&
                   first_stats.directory_batch_memory_byte_limit ==
                       (size_t)HPHI_SYMMETRY_DIRECTORY_MEMORY_BYTES,
               "batch fast stats mismatch");

  memset(&options, 0, sizeof(options));
  options.force_chunked = 1;
  options.debug_echo = 1;
  options.corrupt_response_rank = -1;
  for (index = 0U; index < count; index++) {
    global_beta[index] = ULONG_MAX;
    norm[index] = -1.0;
  }
  require_true(SymmetryResolveRepresentativeBatchWithOptions(
                   directory, keys, count, global_beta, norm,
                   &options) == 0,
               "batch forced-chunk resolution failed");
  assert_batch_result(
      &fixture, keys, count, global_beta, norm,
      "batch forced-chunk result mismatch");
  require_true(GetSymmetryRepresentativeDirectoryBatchStats(
                   directory, &second_stats) == 0 &&
                   second_stats.directory_batch_calls == UINT64_C(2) &&
                   second_stats.directory_request_entries_sent ==
                       count * UINT64_C(2) &&
                   second_stats.directory_request_entries_received ==
                       first_stats.directory_request_entries_received *
                           UINT64_C(2) &&
                   second_stats.directory_found_entries ==
                       first_stats.directory_found_entries * UINT64_C(2) &&
                   second_stats.directory_not_found_entries ==
                       first_stats.directory_not_found_entries * UINT64_C(2) &&
                   second_stats.directory_batch_temporary_peak_bytes >=
                       first_stats.directory_batch_temporary_peak_bytes &&
                   (test_nrank == 1 ||
                    second_stats.directory_exchange_used_chunked == 1),
               "batch forced-chunk stats mismatch");

  require_true(SymmetryResolveRepresentativeBatch(
                   directory, NULL, 0U, NULL, NULL) == 0,
               "zero-request batch failed");
  require_true(GetSymmetryRepresentativeDirectoryBatchStats(
                   directory, &third_stats) == 0 &&
                   third_stats.directory_batch_calls == UINT64_C(3) &&
                   third_stats.directory_batch_temporary_peak_bytes ==
                       second_stats.directory_batch_temporary_peak_bytes,
               "zero-request batch stats mismatch");

  FreeSymmetryRepresentativeDirectory(directory);
  free_directory_fixture(&fixture);
  free(keys);
  free(global_beta);
  free(norm);
}

static void assert_zero_dimension_batch(void)
{
  const unsigned long int keys[3] = {0UL, 7UL, ULONG_MAX};
  unsigned long int global_beta[3] = {11UL, 12UL, 13UL};
  double norm[3] = {-1.0, -2.0, -3.0};
  struct DirectoryFixture fixture;
  struct SymmetryRepresentativeDirectory *directory =
      build_batch_directory(
          &fixture, 0UL, 0UL, "zero-dimension batch directory failed");
  struct SymmetryRepresentativeBatchStats stats;
  require_true(SymmetryResolveRepresentativeBatch(
                   directory, keys, UINT64_C(3),
                   global_beta, norm) == 0,
               "zero-dimension batch resolution failed");
  require_true(global_beta[0] == 0UL && global_beta[1] == 0UL &&
                   global_beta[2] == 0UL &&
                   norm[0] == 0.0 && !signbit(norm[0]) &&
                   norm[1] == 0.0 && !signbit(norm[1]) &&
                   norm[2] == 0.0 && !signbit(norm[2]) &&
                   GetSymmetryRepresentativeDirectoryBatchStats(
                       directory, &stats) == 0 &&
                   stats.directory_batch_calls == UINT64_C(1) &&
                   stats.directory_not_found_entries == UINT64_C(3) &&
                   stats.directory_request_entries_sent == 0U &&
                   stats.directory_request_entries_received == 0U &&
                   stats.directory_batch_temporary_peak_bytes == 0U,
               "zero-dimension batch result or stats mismatch");
  FreeSymmetryRepresentativeDirectory(directory);
  free_directory_fixture(&fixture);
}

static void assert_mixed_zero_request_batch(void)
{
  struct DirectoryFixture fixture;
  struct SymmetryRepresentativeDirectory *directory;
  struct SymmetryRepresentativeBatchStats stats;
  unsigned long int key = 700UL;
  unsigned long int global_beta = ULONG_MAX;
  double norm = -1.0;
  uint64_t count = (test_nrank == 1 || (test_rank % 2) == 0)
                       ? UINT64_C(1)
                       : UINT64_C(0);
  directory = build_batch_directory(
      &fixture, (unsigned long int)test_nrank * 2UL, key,
      "mixed-zero batch directory failed");
  require_true(SymmetryResolveRepresentativeBatch(
                   directory, count > 0U ? &key : NULL, count,
                   count > 0U ? &global_beta : NULL,
                   count > 0U ? &norm : NULL) == 0,
               "mixed-zero batch resolution failed");
  if (count > 0U) {
    require_true(global_beta == 1UL && norm == 0.5,
                 "mixed-zero batch result mismatch");
  }
  require_true(GetSymmetryRepresentativeDirectoryBatchStats(
                   directory, &stats) == 0 &&
                   stats.directory_batch_calls == UINT64_C(1) &&
                   stats.directory_request_entries_sent == count,
               "mixed-zero batch stats mismatch");
  FreeSymmetryRepresentativeDirectory(directory);
  free_directory_fixture(&fixture);
}

static void assert_batch_failure_atomicity(void)
{
  struct DirectoryFixture fixture;
  struct SymmetryRepresentativeDirectory *directory;
  struct SymmetryRepresentativeBatchStats before;
  struct SymmetryRepresentativeBatchStats after;
  struct SymmetryRepresentativeBatchOptions options;
  unsigned long int dim = (unsigned long int)test_nrank * 4UL;
  unsigned long int first_key = 1000UL;
  unsigned long int invalid_keys[2];
  unsigned long int global_beta[2] = {71UL, 72UL};
  double norm[2] = {-73.0, -74.0};
  int failing_rank = test_nrank > 1 ? 1 : 0;
  directory = build_batch_directory(
      &fixture, dim, first_key, "atomic batch directory build failed");
  require_true(GetSymmetryRepresentativeDirectoryBatchStats(
                   directory, &before) == 0,
               "atomic batch initial stats failed");

  if (test_rank == failing_rank) {
    invalid_keys[0] = fixture_key(first_key, 1UL);
    invalid_keys[1] = fixture_key(first_key, 0UL);
  } else {
    invalid_keys[0] = fixture_key(first_key, 0UL);
    invalid_keys[1] = fixture_key(first_key, 1UL);
  }
  require_true(SymmetryResolveRepresentativeBatch(
                   directory, invalid_keys, UINT64_C(2),
                   global_beta, norm) != 0 &&
                   global_beta[0] == 71UL && global_beta[1] == 72UL &&
                   norm[0] == -73.0 && norm[1] == -74.0 &&
                   GetSymmetryRepresentativeDirectoryBatchStats(
                       directory, &after) == 0 &&
                   memcmp(&before, &after, sizeof(before)) == 0,
               "unsorted batch failure was not atomic");

  if (test_rank == failing_rank) {
    invalid_keys[0] = fixture_key(first_key, 0UL);
    invalid_keys[1] = fixture_key(first_key, 0UL);
  }
  require_true(SymmetryResolveRepresentativeBatch(
                   directory, invalid_keys, UINT64_C(2),
                   global_beta, norm) != 0 &&
                   global_beta[0] == 71UL && global_beta[1] == 72UL &&
                   norm[0] == -73.0 && norm[1] == -74.0 &&
                   GetSymmetryRepresentativeDirectoryBatchStats(
                       directory, &after) == 0 &&
                   memcmp(&before, &after, sizeof(before)) == 0,
               "duplicate batch failure was not atomic");

  invalid_keys[0] = first_key;
  memset(&options, 0, sizeof(options));
  options.corrupt_response_rank = 0;
  require_true(SymmetryResolveRepresentativeBatchWithOptions(
                   directory, invalid_keys, UINT64_C(1),
                   global_beta, norm, &options) != 0 &&
                   global_beta[0] == 71UL && norm[0] == -73.0 &&
                   GetSymmetryRepresentativeDirectoryBatchStats(
                       directory, &after) == 0 &&
                   memcmp(&before, &after, sizeof(before)) == 0,
               "corrupt response failure was not atomic");

  require_true(SymmetryResolveRepresentativeBatch(
                   directory, invalid_keys, UINT64_C(1),
                   global_beta, norm) == 0 &&
                   global_beta[0] == 1UL && norm[0] == 0.5,
               "valid batch recovery after failures failed");
  FreeSymmetryRepresentativeDirectory(directory);
  free_directory_fixture(&fixture);
}

static void assert_batch_memory_cap_recovery(void)
{
  struct DirectoryFixture fixture;
  struct SymmetryRepresentativeDirectory *directory;
  struct SymmetryRepresentativeBatchStats before;
  struct SymmetryRepresentativeBatchStats after;
  uint64_t count =
      test_nrank == 1 ? UINT64_C(2050) : UINT64_C(300);
  unsigned long int *keys =
      (unsigned long int *)malloc((size_t)count * sizeof(*keys));
  unsigned long int *global_beta =
      (unsigned long int *)malloc((size_t)count * sizeof(*global_beta));
  double *norm = (double *)malloc((size_t)count * sizeof(*norm));
  uint64_t index;
  require_true(keys != NULL && global_beta != NULL && norm != NULL,
               "cap fixture allocation failed");
  directory = build_batch_directory(
      &fixture, (unsigned long int)test_nrank * 4UL, 100000UL,
      "cap directory build failed");
  for (index = 0U; index < count; index++) {
    keys[index] = (unsigned long int)index;
    global_beta[index] = ULONG_MAX;
    norm[index] = -1.0;
  }
  require_true(GetSymmetryRepresentativeDirectoryBatchStats(
                   directory, &before) == 0 &&
                   SymmetryResolveRepresentativeBatch(
                       directory, keys, count, global_beta, norm) != 0 &&
                   global_beta[0] == ULONG_MAX &&
                   global_beta[count - 1U] == ULONG_MAX &&
                   norm[0] == -1.0 && norm[count - 1U] == -1.0 &&
                   GetSymmetryRepresentativeDirectoryBatchStats(
                       directory, &after) == 0 &&
                   memcmp(&before, &after, sizeof(before)) == 0,
               "batch memory cap failure was not atomic");
  keys[0] = 100000UL;
  require_true(SymmetryResolveRepresentativeBatch(
                   directory, keys, UINT64_C(1),
                   global_beta, norm) == 0 &&
                   global_beta[0] == 1UL && norm[0] == 0.5,
               "batch recovery after cap failure failed");
  FreeSymmetryRepresentativeDirectory(directory);
  free_directory_fixture(&fixture);
  free(keys);
  free(global_beta);
  free(norm);
}

static unsigned long int apply_hamiltonian_hop(
    unsigned long int state,
    unsigned int from,
    unsigned int to,
    const char *label)
{
  unsigned long int from_mask;
  unsigned long int to_mask;
  require_true(from < sizeof(state) * CHAR_BIT &&
                   to < sizeof(state) * CHAR_BIT && from != to,
               label);
  from_mask = 1UL << from;
  to_mask = 1UL << to;
  require_true((state & from_mask) != 0UL &&
                   (state & to_mask) == 0UL,
               label);
  return state ^ from_mask ^ to_mask;
}

static void replicated_canonical_lookup(
    const unsigned long int *basis_keys,
    unsigned long int dim,
    unsigned long int target,
    unsigned long int *global_beta,
    double *norm)
{
  unsigned long int lo = 0UL;
  unsigned long int hi = dim;
  *global_beta = 0UL;
  *norm = 0.0;
  while (lo < hi) {
    unsigned long int mid = lo + (hi - lo) / 2UL;
    if (basis_keys[mid] < target) {
      lo = mid + 1UL;
    } else {
      hi = mid;
    }
  }
  if (lo < dim && basis_keys[lo] == target) {
    *global_beta = lo + 1UL;
    *norm = 1.0 + 0.125 * (double)lo;
  }
}

static void assert_transition_model_batch(
    const unsigned long int *basis_keys,
    unsigned long int dim,
    const unsigned long int targets[3],
    const char *label)
{
  struct DirectoryFixture fixture;
  struct SymmetryRepresentativeDirectory *directory = NULL;
  unsigned long int global_beta[3] = {ULONG_MAX, ULONG_MAX, ULONG_MAX};
  double norm[3] = {-1.0, -1.0, -1.0};
  int category_count[3] = {0, 0, 0};
  uint64_t index;
  init_explicit_directory_fixture(basis_keys, dim, &fixture);
  require_true(BuildSymmetryRepresentativeDirectory(
                   fixture.basis, fixture.dim, fixture.local_dim,
                   fixture.local_capacity, fixture.local_offset,
                   fixture.rank_offsets, test_rank, test_nrank,
                   &directory) == 0,
               label);
  require_true(targets[0] < targets[1] && targets[1] < targets[2] &&
                   SymmetryResolveRepresentativeBatch(
                       directory, targets, UINT64_C(3),
                       global_beta, norm) == 0,
               label);
  for (index = 0U; index < UINT64_C(3); index++) {
    unsigned long int expected_beta;
    double expected_norm;
    replicated_canonical_lookup(
        basis_keys, dim, targets[index], &expected_beta, &expected_norm);
    require_true(global_beta[index] == expected_beta &&
                     memcmp(&norm[index], &expected_norm,
                            sizeof(expected_norm)) == 0,
                 label);
    if (expected_beta != 0UL) {
      category_count[0]++;
    } else if (targets[index] > basis_keys[0] &&
               targets[index] < basis_keys[dim - 1UL]) {
      category_count[1]++;
    } else {
      category_count[2]++;
    }
  }
  require_true(category_count[0] > 0 &&
                   category_count[1] > 0 &&
                   category_count[2] > 0,
               label);
  FreeSymmetryRepresentativeDirectory(directory);
  free_directory_fixture(&fixture);
}

static void assert_hamiltonian_transition_batches(void)
{
  static const unsigned long int spin_basis[] = {
      5UL, 9UL, 10UL, 12UL};
  static const unsigned long int spinless_basis[] = {
      6UL, 9UL, 12UL, 17UL, 18UL, 20UL, 24UL};
  static const unsigned long int hubbard_basis[] = {
      10UL, 12UL, 17UL, 18UL, 24UL, 33UL, 34UL, 40UL, 48UL};
  unsigned long int spin_targets[3];
  unsigned long int spinless_targets[3];
  unsigned long int hubbard_targets[3];

  /*
   * Identity canonicalization makes each Hamiltonian target its own
   * representative. The sparse replicated lists intentionally contain one
   * found target, one sector/stabilizer-incompatible interior target, and one
   * target outside the represented range for each model.
   */
  spin_targets[0] = apply_hamiltonian_hop(
      10UL, 3U, 0U, "spin range transition failed");
  spin_targets[1] = apply_hamiltonian_hop(
      9UL, 3U, 2U, "spin found transition failed");
  spin_targets[2] = apply_hamiltonian_hop(
      5UL, 0U, 1U, "spin interior transition failed");

  spinless_targets[0] = apply_hamiltonian_hop(
      17UL, 4U, 1U, "spinless range transition failed");
  spinless_targets[1] = apply_hamiltonian_hop(
      18UL, 4U, 3U, "spinless interior transition failed");
  spinless_targets[2] = apply_hamiltonian_hop(
      6UL, 1U, 3U, "spinless found transition failed");

  hubbard_targets[0] = apply_hamiltonian_hop(
      17UL, 4U, 3U, "Hubbard range transition failed");
  hubbard_targets[1] = apply_hamiltonian_hop(
      10UL, 1U, 2U, "Hubbard found transition failed");
  hubbard_targets[2] = apply_hamiltonian_hop(
      18UL, 1U, 2U, "Hubbard interior transition failed");

  assert_transition_model_batch(
      spin_basis, sizeof(spin_basis) / sizeof(spin_basis[0]),
      spin_targets, "spin transition batch mismatch");
  assert_transition_model_batch(
      spinless_basis,
      sizeof(spinless_basis) / sizeof(spinless_basis[0]),
      spinless_targets, "spinless transition batch mismatch");
  assert_transition_model_batch(
      hubbard_basis, sizeof(hubbard_basis) / sizeof(hubbard_basis[0]),
      hubbard_targets, "Hubbard transition batch mismatch");
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
  assert_zero_dimension_batch();
  assert_mixed_zero_request_batch();
  assert_batch_resolution();
  assert_batch_failure_atomicity();
  assert_batch_memory_cap_recovery();
  assert_hamiltonian_transition_batches();
  if (test_rank == 0) {
    printf("symmetry directory batch gate: PASS (%d rank%s)\n",
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
