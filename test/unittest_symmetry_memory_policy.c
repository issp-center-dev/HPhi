#define _POSIX_C_SOURCE 200809L

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "symmetry_memory_policy.h"

#ifdef MPI
static int test_rank = 0;
#endif

static void require_true(int condition, const char *message)
{
  if (condition) return;
  fprintf(stderr, "symmetry memory policy test failed: %s\n", message);
#ifdef MPI
  {
    int initialized = 0;
    if (MPI_Initialized(&initialized) == MPI_SUCCESS && initialized != 0) {
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
#endif
  exit(1);
}

static void set_memory_environment(const char *warning, const char *limit)
{
  if (warning == NULL) {
    require_true(
        unsetenv(HPHI_SYMMETRY_MEMORY_WARN_ENV) == 0,
        "failed to unset warning threshold");
  } else {
    require_true(
        setenv(HPHI_SYMMETRY_MEMORY_WARN_ENV, warning, 1) == 0,
        "failed to set warning threshold");
  }
  if (limit == NULL) {
    require_true(
        unsetenv(HPHI_SYMMETRY_MEMORY_LIMIT_ENV) == 0,
        "failed to unset hard limit");
  } else {
    require_true(
        setenv(HPHI_SYMMETRY_MEMORY_LIMIT_ENV, limit, 1) == 0,
        "failed to set hard limit");
  }
}

static void test_serial_policy(void)
{
  struct SymmetryMemoryPolicy policy;
  set_memory_environment(NULL, NULL);
  require_true(
      SymmetryLoadMemoryPolicy(
          0U, 0, 0, "test", &policy) == 0 &&
          policy.warning_byte_threshold ==
              (uint64_t)HPHI_SYMMETRY_MEMORY_WARN_BYTES &&
          policy.hard_byte_limit == 0U,
      "default policy mismatch");
  require_true(
      SymmetryCheckMemoryPolicy(
          &policy, policy.warning_byte_threshold + 1U,
          0, 0, "test", "default-warning", 1) == 1,
      "default threshold did not warn");
  require_true(
      SymmetryCheckMemoryPolicy(
          &policy, policy.warning_byte_threshold + 1U,
          0, 0, "test", "suppressed-warning", 0) == 0,
      "warning suppression failed");

  set_memory_environment("1024", "2048");
  require_true(
      SymmetryLoadMemoryPolicy(
          4096U, 0, 0, "test", &policy) == 0 &&
          policy.warning_byte_threshold == 1024U &&
          policy.hard_byte_limit == 2048U,
      "runtime policy did not override compile-time fallback");
  require_true(
      SymmetryCheckMemoryPolicy(
          &policy, 1536U, 0, 0, "test", "runtime-warning", 1) == 1,
      "runtime warning threshold failed");
  require_true(
      SymmetryCheckMemoryPolicy(
          &policy, 2049U, 0, 0, "test", "runtime-limit", 1) == -1,
      "runtime hard limit failed");

  set_memory_environment("0", NULL);
  require_true(
      SymmetryLoadMemoryPolicy(
          4096U, 0, 0, "test", &policy) == 0 &&
          policy.warning_byte_threshold == 0U &&
          policy.hard_byte_limit == 4096U,
      "compile-time hard-limit fallback mismatch");
  require_true(
      SymmetryCheckMemoryPolicy(
          &policy, 4097U, 0, 0, "test", "compile-limit", 1) == -1,
      "compile-time hard limit was not retained");

  set_memory_environment("0", "0");
  require_true(
      SymmetryLoadMemoryPolicy(
          4096U, 0, 0, "test", &policy) == 0 &&
          policy.warning_byte_threshold == 0U &&
          policy.hard_byte_limit == 0U &&
          SymmetryCheckMemoryPolicy(
              &policy, UINT64_MAX, 0, 0,
              "test", "unlimited", 1) == 0,
      "explicit unlimited policy failed");

  set_memory_environment("1GiB", "0");
  require_true(
      SymmetryLoadMemoryPolicy(
          0U, 0, 0, "test", &policy) == -1,
      "invalid warning value was accepted");
  set_memory_environment(NULL, NULL);
}

#ifdef MPI
static void test_mpi_policy(int rank, int nrank)
{
  struct SymmetryMemoryPolicy policy;
  uint64_t observed;
  require_true(nrank >= 2, "MPI policy test requires at least two ranks");
  if (rank == 0) {
    set_memory_environment("1024", "0");
  } else {
    set_memory_environment("2048", "4096");
  }
  require_true(
      SymmetryLoadMemoryPolicy(
          0U, 1, rank, "test-mpi", &policy) == 0 &&
          policy.warning_byte_threshold == 1024U &&
          policy.hard_byte_limit == 0U,
      "rank-zero policy was not broadcast");
  observed = rank == nrank - 1 ? 2048U : 512U;
  require_true(
      SymmetryCheckMemoryPolicy(
          &policy, observed, 1, rank,
          "test-mpi", "collective-warning", 1) == 1,
      "collective maximum did not warn on every rank");

  if (rank == 0) {
    set_memory_environment("1024", "1536");
  } else {
    set_memory_environment("0", "0");
  }
  require_true(
      SymmetryLoadMemoryPolicy(
          0U, 1, rank, "test-mpi", &policy) == 0 &&
          policy.warning_byte_threshold == 1024U &&
          policy.hard_byte_limit == 1536U,
      "collective hard-limit policy mismatch");
  require_true(
      SymmetryCheckMemoryPolicy(
          &policy, observed, 1, rank,
          "test-mpi", "collective-limit", 1) == -1,
      "collective hard limit did not fail on every rank");
  set_memory_environment(NULL, NULL);
}
#endif

int main(int argc, char **argv)
{
  int mpi_mode = argc == 2 && strcmp(argv[1], "--mpi") == 0;
  if (argc > 2 || (argc == 2 && !mpi_mode)) {
    fprintf(stderr, "usage: %s [--mpi]\n", argv[0]);
    return 2;
  }
#ifdef MPI
  if (mpi_mode) {
    int nrank = 1;
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS ||
        MPI_Comm_rank(MPI_COMM_WORLD, &test_rank) != MPI_SUCCESS ||
        MPI_Comm_size(MPI_COMM_WORLD, &nrank) != MPI_SUCCESS) {
      return 1;
    }
    test_mpi_policy(test_rank, nrank);
    if (MPI_Finalize() != MPI_SUCCESS) return 1;
    return 0;
  }
#else
  if (mpi_mode) {
    fprintf(stderr, "MPI mode requested from a non-MPI build\n");
    return 2;
  }
#endif
  test_serial_policy();
  return 0;
}
