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

struct SymmetryMpiUnsignedLongResult {
  /* Zero-origin payload segmented by ascending source rank. */
  unsigned long int *entries;
  uint64_t count;
  uint64_t *counts;
  uint64_t *displacements;
  int nrank;
};

struct SymmetryMpiLookupResponse {
  /* Global basis indices are one-origin. Zero denotes not found. */
  unsigned long int global_beta;
  double norm;
};

struct SymmetryMpiLookupResponseResult {
  /* Zero-origin payload segmented by ascending source rank. */
  struct SymmetryMpiLookupResponse *entries;
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

/*
 * Exchange unsigned-long payloads with the same checked schedule, message
 * cap, source-rank ordering, and failure contract as the basis-vector API.
 * Receive counts are discovered collectively from send_layout.
 */
int SymmetryMpiExchangeUnsignedLongs(
    const unsigned long int *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    int rank,
    int nrank,
    const struct SymmetryMpiExchangeOptions *options,
    struct SymmetryMpiUnsignedLongResult *result,
    struct SymmetryMpiExchangeStats *stats);

void FreeSymmetryMpiUnsignedLongResult(
    struct SymmetryMpiUnsignedLongResult *result);

/*
 * Exchange lookup responses using the exact reverse schedule already known
 * from a preceding request exchange. known_receive_layout is copied into
 * result and must be the zero-based source-rank prefix layout expected by
 * this rank, including the peerwise transpose of send_layout. No counts
 * exchange is performed on this path. The caller must pass an empty result;
 * failure leaves it unchanged.
 */
int SymmetryMpiExchangeLookupResponsesKnownLayout(
    const struct SymmetryMpiLookupResponse *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    const struct SymmetryMpiExchangeLayout *known_receive_layout,
    int rank,
    int nrank,
    const struct SymmetryMpiExchangeOptions *options,
    struct SymmetryMpiLookupResponseResult *result,
    struct SymmetryMpiExchangeStats *stats);

void FreeSymmetryMpiLookupResponseResult(
    struct SymmetryMpiLookupResponseResult *result);

#ifdef MPI
int SymmetryMpiCreateBasisVectorType(MPI_Datatype *vector_type);
int SymmetryMpiCreateLookupResponseType(MPI_Datatype *response_type);
#endif

#endif /* HPHI_SYMMETRY_MPI_EXCHANGE_H */
