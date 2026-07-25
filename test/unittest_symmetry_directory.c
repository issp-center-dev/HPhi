#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "symmetry_basis.h"
#include "symmetry_directory.h"

static void fail_test(const char *label)
{
  fprintf(stderr, "%s\n", label);
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

int main(void)
{
  assert_basic_lookup();
  assert_empty_index();
  assert_collision_accounting();
  assert_invalid_builds_and_recovery();
  assert_lookup_failure_atomicity();
  assert_parallel_read_lookup();
  printf("symmetry directory local index tests passed\n");
  return 0;
}
