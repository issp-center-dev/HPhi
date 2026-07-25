#ifndef HPHI_SYMMETRY_MPI_EXCHANGE_H
#define HPHI_SYMMETRY_MPI_EXCHANGE_H

#include <stdint.h>

#ifndef HPHI_SYMMETRY_EXCHANGE_MESSAGE_BYTES
#define HPHI_SYMMETRY_EXCHANGE_MESSAGE_BYTES UINT64_C(1073741824)
#endif

#ifdef MPI
#include <mpi.h>
#endif

struct SymmetryBasisVector;

struct SymmetryMpiExchangeLayout {
  int nrank;
  uint64_t count;
  const uint64_t *counts;
  /* Strict source/destination-rank prefix layout, starting at zero. */
  const uint64_t *displacements;
};

struct SymmetryMpiExchangeOptions {
  /*
   * Zero selects the entry limit derived from the compile-time message byte
   * cap. A nonzero value may lower, but not raise, that limit.
   * The effective value must be identical on every active MPI rank;
   * the exchange checks this collectively and rejects a mismatch.
   */
  uint64_t chunk_limit;
  /* Internal/test control. Production callers normally pass FALSE. */
  int force_chunked;
};

struct SymmetryMpiExchangeStats {
  int used_chunked;
  uint64_t send_entries;
  uint64_t recv_entries;
  uint64_t send_messages;
  uint64_t recv_messages;
  uint64_t message_entry_limit;
  uint64_t message_byte_limit;
  uint64_t max_message_entries;
  uint64_t max_message_bytes;
};

struct SymmetryMpiExchangeResult {
  /* Zero-origin payload segmented by ascending source rank. */
  struct SymmetryBasisVector *entries;
  uint64_t count;
  uint64_t *counts;
  uint64_t *displacements;
  int nrank;
};

/* Shared collective preflight helpers for symmetry MPI modules. */
int SymmetryMpiCollectivesActive(void);
int SymmetryMpiAgreeError(int mpi_active, int local_error);

/*
 * The send layout must be a zero-based prefix partition by destination rank.
 * On success result owns entries/counts/displacements and must be freed with
 * FreeSymmetryMpiExchangeResult(). The caller must pass an empty result.
 * On failure result is unchanged, and collective preflight failures are
 * reported consistently by all active MPI ranks.
 */
int SymmetryMpiExchangeBasisVectors(
    const struct SymmetryBasisVector *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    int rank,
    int nrank,
    const struct SymmetryMpiExchangeOptions *options,
    struct SymmetryMpiExchangeResult *result,
    struct SymmetryMpiExchangeStats *stats);

void FreeSymmetryMpiExchangeResult(
    struct SymmetryMpiExchangeResult *result);

#ifdef MPI
int SymmetryMpiCreateBasisVectorType(MPI_Datatype *vector_type);
#endif

#endif /* HPHI_SYMMETRY_MPI_EXCHANGE_H */
