#ifndef HPHI_SYMMETRY_MEMORY_POLICY_H
#define HPHI_SYMMETRY_MEMORY_POLICY_H

#include <stdint.h>

#ifndef HPHI_SYMMETRY_MEMORY_WARN_BYTES
#define HPHI_SYMMETRY_MEMORY_WARN_BYTES UINT64_C(1073741824)
#endif

#ifndef HPHI_SYMMETRY_MEMORY_LIMIT_BYTES
#define HPHI_SYMMETRY_MEMORY_LIMIT_BYTES UINT64_C(0)
#endif

#define HPHI_SYMMETRY_MEMORY_WARN_ENV \
  "HPHI_SYMMETRY_MEMORY_WARN_BYTES"
#define HPHI_SYMMETRY_MEMORY_LIMIT_ENV \
  "HPHI_SYMMETRY_MEMORY_LIMIT_BYTES"

struct SymmetryMemoryPolicy {
  uint64_t warning_byte_threshold;
  uint64_t hard_byte_limit;
};

/*
 * Rank zero reads the runtime environment and broadcasts one policy.
 *
 * HPHI_SYMMETRY_MEMORY_WARN_BYTES defaults to 1 GiB.  Zero disables
 * warnings.  HPHI_SYMMETRY_MEMORY_LIMIT_BYTES defaults to
 * compile_time_hard_limit; zero means unlimited.  The component-specific
 * compile-time limits are retained as backward-compatible hard limits.
 */
int SymmetryLoadMemoryPolicy(
    uint64_t compile_time_hard_limit,
    int mpi_active,
    int rank,
    const char *component,
    struct SymmetryMemoryPolicy *policy);

/*
 * Collectively compare the maximum observed per-rank temporary peak with the
 * policy.  Returns -1 on MPI error or hard-limit violation, 1 when a warning
 * was emitted, and 0 otherwise.  emit_warning lets a caller suppress repeated
 * warnings while preserving hard-limit enforcement.
 */
int SymmetryCheckMemoryPolicy(
    const struct SymmetryMemoryPolicy *policy,
    uint64_t observed_bytes,
    int mpi_active,
    int rank,
    const char *component,
    const char *stage,
    int emit_warning);

#endif /* HPHI_SYMMETRY_MEMORY_POLICY_H */
