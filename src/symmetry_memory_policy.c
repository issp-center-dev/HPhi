#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "symmetry_memory_policy.h"

static int parse_memory_bytes(
    const char *name,
    uint64_t fallback,
    uint64_t *value)
{
  const char *text;
  char *end = NULL;
  unsigned long long parsed;
  if (name == NULL || value == NULL) return -1;
  text = getenv(name);
  if (text == NULL) {
    *value = fallback;
    return 0;
  }
  if (*text == '\0' || *text == '+' || *text == '-') return -1;
  errno = 0;
  parsed = strtoull(text, &end, 10);
  if (errno == ERANGE || end == text || *end != '\0' ||
      parsed > (unsigned long long)UINT64_MAX) {
    return -1;
  }
  *value = (uint64_t)parsed;
  return 0;
}
int SymmetryLoadMemoryPolicy(
    uint64_t compile_time_hard_limit,
    int mpi_active,
    int rank,
    const char *component,
    struct SymmetryMemoryPolicy *policy)
{
  struct SymmetryMemoryPolicy next;
  int status = 0;
  if (policy == NULL || component == NULL || rank < 0) return -1;
  next.warning_byte_threshold =
      (uint64_t)HPHI_SYMMETRY_MEMORY_WARN_BYTES;
  next.hard_byte_limit = compile_time_hard_limit;
  if (rank == 0 &&
      (parse_memory_bytes(
           HPHI_SYMMETRY_MEMORY_WARN_ENV,
           next.warning_byte_threshold,
           &next.warning_byte_threshold) != 0 ||
       parse_memory_bytes(
           HPHI_SYMMETRY_MEMORY_LIMIT_ENV,
           next.hard_byte_limit,
           &next.hard_byte_limit) != 0)) {
    fprintf(stderr,
            "Error: HPhi symmetry %s memory policy requires unsigned "
            "decimal byte values for %s and %s; zero disables the "
            "corresponding threshold.\n",
            component, HPHI_SYMMETRY_MEMORY_WARN_ENV,
            HPHI_SYMMETRY_MEMORY_LIMIT_ENV);
    fflush(stderr);
    status = -1;
  }
#ifdef MPI
  if (mpi_active != 0) {
    if (MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
      return -1;
    }
    if (status == 0 &&
        (MPI_Bcast(
             &next.warning_byte_threshold, 1, MPI_UINT64_T,
             0, MPI_COMM_WORLD) != MPI_SUCCESS ||
         MPI_Bcast(
             &next.hard_byte_limit, 1, MPI_UINT64_T,
             0, MPI_COMM_WORLD) != MPI_SUCCESS)) {
      return -1;
    }
  }
#else
  (void)mpi_active;
#endif
  if (status != 0) return -1;
  *policy = next;
  return 0;
}

int SymmetryCheckMemoryPolicy(
    const struct SymmetryMemoryPolicy *policy,
    uint64_t observed_bytes,
    int mpi_active,
    int rank,
    const char *component,
    const char *stage,
    int emit_warning)
{
  uint64_t maximum_observed = observed_bytes;
  int all_emit_warning = emit_warning != 0 ? 1 : 0;
  if (policy == NULL || component == NULL || stage == NULL || rank < 0) {
    return -1;
  }
#ifdef MPI
  if (mpi_active != 0) {
    if (MPI_Allreduce(
            &observed_bytes, &maximum_observed, 1, MPI_UINT64_T,
            MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(
            MPI_IN_PLACE, &all_emit_warning, 1, MPI_INT,
            MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS) {
      if (rank == 0) {
        fprintf(stderr,
                "Error: HPhi symmetry %s %s memory-policy agreement "
                "failed.\n",
                component, stage);
        fflush(stderr);
      }
      return -1;
    }
  }
#else
  (void)mpi_active;
#endif
  if (policy->hard_byte_limit != 0U &&
      maximum_observed > policy->hard_byte_limit) {
    if (rank == 0) {
      fprintf(stderr,
              "Error: HPhi symmetry %s %s temporary memory may reach "
              "%" PRIu64 " bytes/rank, exceeding hard limit=%" PRIu64
              ". Increase MPI ranks or set %s to a larger value "
              "(zero means unlimited).\n",
              component, stage, maximum_observed,
              policy->hard_byte_limit,
              HPHI_SYMMETRY_MEMORY_LIMIT_ENV);
      fflush(stderr);
    }
    return -1;
  }
  if (all_emit_warning != 0 &&
      policy->warning_byte_threshold != 0U &&
      maximum_observed > policy->warning_byte_threshold) {
    if (rank == 0) {
      if (policy->hard_byte_limit == 0U) {
        fprintf(stderr,
                "Warning: HPhi symmetry %s %s temporary memory may "
                "reach %" PRIu64 " bytes/rank; warning threshold=%"
                PRIu64 ", hard limit=unlimited. Set %s to change or "
                "disable the warning, and %s to set an optional hard "
                "limit.\n",
                component, stage, maximum_observed,
                policy->warning_byte_threshold,
                HPHI_SYMMETRY_MEMORY_WARN_ENV,
                HPHI_SYMMETRY_MEMORY_LIMIT_ENV);
      } else {
        fprintf(stderr,
                "Warning: HPhi symmetry %s %s temporary memory may "
                "reach %" PRIu64 " bytes/rank; warning threshold=%"
                PRIu64 ", hard limit=%" PRIu64 ".\n",
                component, stage, maximum_observed,
                policy->warning_byte_threshold,
                policy->hard_byte_limit);
      }
      fflush(stderr);
    }
    return 1;
  }
  return 0;
}
