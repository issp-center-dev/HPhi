#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "DefCommon.h"
#include "symmetry_basis.h"
#include "symmetry_checked.h"
#include "symmetry_mpi_exchange.h"

#define SYMMETRY_BASIS_EXCHANGE_TAG 23173
#define SYMMETRY_UNSIGNED_LONG_EXCHANGE_TAG 23180
#define SYMMETRY_LOOKUP_RESPONSE_EXCHANGE_TAG 23181
#define SYMMETRY_UNSIGNED_LONG_ECHO_EXCHANGE_TAG 23182

struct SymmetryMpiRawResult {
  void *entries;
  uint64_t count;
  uint64_t *counts;
  uint64_t *displacements;
  int nrank;
};

struct SymmetryMpiPayloadDescriptor {
  size_t extent;
  int tag;
#ifdef MPI
  MPI_Datatype builtin_type;
  int (*create_type)(MPI_Datatype *);
#endif
};

static int checked_budget_peak(
    const struct SymmetryMpiExchangeMemoryBudget *budget,
    size_t result_bytes,
    size_t workspace_bytes)
{
  size_t peak;
  if (budget == NULL) return 0;
  if (budget->byte_limit == 0U ||
      budget->live_bytes > budget->byte_limit ||
      SymmetryCheckedSizeAdd(budget->live_bytes, result_bytes, &peak) != 0 ||
      SymmetryCheckedSizeAdd(peak, workspace_bytes, &peak) != 0 ||
      peak > budget->byte_limit) {
    return -1;
  }
  return 0;
}

static int exchange_workspace_bytes(int nrank,
                                    int mpi_active,
                                    int use_chunked,
                                    size_t *workspace_bytes)
{
#ifdef MPI
  size_t one_array;
  size_t request_count;
  size_t request_bytes;
  size_t total;
#endif
  if (workspace_bytes == NULL || nrank < 1) return -1;
  *workspace_bytes = 0U;
#ifdef MPI
  if (mpi_active == FALSE || nrank == 1) return 0;
  if (use_chunked != FALSE) {
    if (SymmetryCheckedSizeMul((size_t)nrank, sizeof(uint64_t),
                               &one_array) != 0 ||
        SymmetryCheckedSizeMul((size_t)nrank, 2U, &request_count) != 0 ||
        SymmetryCheckedSizeMul(request_count, sizeof(MPI_Request),
                               &request_bytes) != 0 ||
        SymmetryCheckedSizeAdd(one_array, one_array, &total) != 0 ||
        SymmetryCheckedSizeAdd(total, request_bytes, &total) != 0) {
      return -1;
    }
  } else {
    if (SymmetryCheckedSizeMul((size_t)nrank, sizeof(int),
                               &one_array) != 0 ||
        SymmetryCheckedSizeMul(one_array, 4U, &total) != 0) {
      return -1;
    }
  }
  *workspace_bytes = total;
  return 0;
#else
  (void)mpi_active;
  (void)use_chunked;
  return 0;
#endif
}

int SymmetryMpiCollectivesActive(void)
{
#ifdef MPI
  int initialized = 0;
  int finalized = 0;
  if (MPI_Initialized(&initialized) != MPI_SUCCESS || initialized == 0) {
    return FALSE;
  }
  if (MPI_Finalized(&finalized) != MPI_SUCCESS || finalized != 0) {
    return FALSE;
  }
  return TRUE;
#else
  return FALSE;
#endif
}

int SymmetryMpiAgreeError(int mpi_active, int local_error)
{
  int error_flag = local_error != 0 ? 1 : 0;
#ifdef MPI
  int global_error = error_flag;
  if (mpi_active != FALSE &&
      MPI_Allreduce(&error_flag, &global_error, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
  return global_error;
#else
  (void)mpi_active;
  return error_flag;
#endif
}

static int get_message_limits(size_t extent,
                              uint64_t *entry_limit,
                              uint64_t *byte_limit)
{
  const uint64_t configured_bytes =
      (uint64_t)HPHI_SYMMETRY_EXCHANGE_MESSAGE_BYTES;
  uint64_t entries;
  if (extent == 0U || entry_limit == NULL || byte_limit == NULL ||
      configured_bytes < (uint64_t)extent) {
    return -1;
  }
  entries = configured_bytes / (uint64_t)extent;
  if (entries > (uint64_t)INT_MAX) entries = (uint64_t)INT_MAX;
  if (entries == 0U ||
      entries > UINT64_MAX / (uint64_t)extent) {
    return -1;
  }
  *entry_limit = entries;
  *byte_limit = entries * (uint64_t)extent;
  return 0;
}

static int validate_layout(const struct SymmetryMpiExchangeLayout *layout,
                           const void *entries,
                           int nrank,
                           size_t extent,
                           int require_entries)
{
  uint64_t total = 0U;
  int peer;
  if (layout == NULL || layout->nrank != nrank || nrank < 1 ||
      extent == 0U || layout->counts == NULL ||
      layout->displacements == NULL) {
    return -1;
  }
  for (peer = 0; peer < nrank; peer++) {
    if (layout->displacements[peer] != total ||
        SymmetryCheckedU64Add(total, layout->counts[peer], &total) != 0) {
      return -1;
    }
  }
  if (total != layout->count ||
      total > (uint64_t)(SIZE_MAX / extent) ||
      (require_entries != FALSE && total > 0U && entries == NULL)) {
    return -1;
  }
  return 0;
}

static void reset_raw_result(struct SymmetryMpiRawResult *result)
{
  if (result == NULL) return;
  result->entries = NULL;
  result->count = 0U;
  result->counts = NULL;
  result->displacements = NULL;
  result->nrank = 0;
}

static void free_raw_result(struct SymmetryMpiRawResult *result)
{
  if (result == NULL) return;
  free(result->entries);
  free(result->counts);
  free(result->displacements);
  reset_raw_result(result);
}

static void reset_stats(struct SymmetryMpiExchangeStats *stats)
{
  if (stats == NULL) return;
  memset(stats, 0, sizeof(*stats));
}

static void reset_basis_result(struct SymmetryMpiExchangeResult *result)
{
  if (result == NULL) return;
  result->entries = NULL;
  result->count = 0U;
  result->counts = NULL;
  result->displacements = NULL;
  result->nrank = 0;
}

static void reset_unsigned_long_result(
    struct SymmetryMpiUnsignedLongResult *result)
{
  if (result == NULL) return;
  result->entries = NULL;
  result->count = 0U;
  result->counts = NULL;
  result->displacements = NULL;
  result->nrank = 0;
}

static void reset_lookup_response_result(
    struct SymmetryMpiLookupResponseResult *result)
{
  if (result == NULL) return;
  result->entries = NULL;
  result->count = 0U;
  result->counts = NULL;
  result->displacements = NULL;
  result->nrank = 0;
}

void FreeSymmetryMpiExchangeResult(
    struct SymmetryMpiExchangeResult *result)
{
  if (result == NULL) return;
  free(result->entries);
  free(result->counts);
  free(result->displacements);
  reset_basis_result(result);
}

void FreeSymmetryMpiUnsignedLongResult(
    struct SymmetryMpiUnsignedLongResult *result)
{
  if (result == NULL) return;
  free(result->entries);
  free(result->counts);
  free(result->displacements);
  reset_unsigned_long_result(result);
}

void FreeSymmetryMpiLookupResponseResult(
    struct SymmetryMpiLookupResponseResult *result)
{
  if (result == NULL) return;
  free(result->entries);
  free(result->counts);
  free(result->displacements);
  reset_lookup_response_result(result);
}

#ifdef MPI
int SymmetryMpiCreateBasisVectorType(MPI_Datatype *vector_type)
{
  struct SymmetryBasisVector sample;
  int block_lengths[6] = {1, 1, 1, 1, 1, 1};
  MPI_Aint base;
  MPI_Aint displacements[6];
  MPI_Datatype member_types[6] = {
    MPI_UNSIGNED_LONG, MPI_UNSIGNED, MPI_UNSIGNED,
    MPI_DOUBLE, MPI_DOUBLE_COMPLEX, MPI_DOUBLE
  };
  MPI_Datatype packed_type = MPI_DATATYPE_NULL;
  int member;
  int ierr;

  if (vector_type == NULL) return -1;
  *vector_type = MPI_DATATYPE_NULL;
  if (MPI_Get_address(&sample, &base) != MPI_SUCCESS ||
      MPI_Get_address(&sample.rep_state, &displacements[0]) != MPI_SUCCESS ||
      MPI_Get_address(&sample.orbit_size, &displacements[1]) != MPI_SUCCESS ||
      MPI_Get_address(&sample.stabilizer_size, &displacements[2]) !=
          MPI_SUCCESS ||
      MPI_Get_address(&sample.norm, &displacements[3]) != MPI_SUCCESS ||
      MPI_Get_address(&sample.stabilizer_character_sum,
                      &displacements[4]) != MPI_SUCCESS ||
      MPI_Get_address(&sample.diagonal, &displacements[5]) != MPI_SUCCESS) {
    return -1;
  }
  for (member = 0; member < 6; member++) {
    displacements[member] -= base;
  }

  ierr = MPI_Type_create_struct(6, block_lengths, displacements, member_types,
                                &packed_type);
  if (ierr != MPI_SUCCESS) return -1;
  ierr = MPI_Type_create_resized(packed_type, 0,
                                 (MPI_Aint)sizeof(struct SymmetryBasisVector),
                                 vector_type);
  MPI_Type_free(&packed_type);
  if (ierr != MPI_SUCCESS) return -1;
  ierr = MPI_Type_commit(vector_type);
  if (ierr != MPI_SUCCESS) {
    MPI_Type_free(vector_type);
    *vector_type = MPI_DATATYPE_NULL;
    return -1;
  }
  return 0;
}

int SymmetryMpiCreateLookupResponseType(MPI_Datatype *response_type)
{
  struct SymmetryMpiLookupResponse sample;
  int block_lengths[2] = {1, 1};
  MPI_Aint base;
  MPI_Aint displacements[2];
  MPI_Datatype member_types[2] = {MPI_UNSIGNED_LONG, MPI_DOUBLE};
  MPI_Datatype packed_type = MPI_DATATYPE_NULL;
  int ierr;

  if (response_type == NULL) return -1;
  *response_type = MPI_DATATYPE_NULL;
  if (MPI_Get_address(&sample, &base) != MPI_SUCCESS ||
      MPI_Get_address(&sample.global_beta, &displacements[0]) != MPI_SUCCESS ||
      MPI_Get_address(&sample.norm, &displacements[1]) != MPI_SUCCESS) {
    return -1;
  }
  displacements[0] -= base;
  displacements[1] -= base;
  ierr = MPI_Type_create_struct(2, block_lengths, displacements, member_types,
                                &packed_type);
  if (ierr != MPI_SUCCESS) return -1;
  ierr = MPI_Type_create_resized(
      packed_type, 0, (MPI_Aint)sizeof(struct SymmetryMpiLookupResponse),
      response_type);
  MPI_Type_free(&packed_type);
  if (ierr != MPI_SUCCESS) return -1;
  ierr = MPI_Type_commit(response_type);
  if (ierr != MPI_SUCCESS) {
    MPI_Type_free(response_type);
    *response_type = MPI_DATATYPE_NULL;
    return -1;
  }
  return 0;
}

static int layout_requires_chunking(
    const struct SymmetryMpiExchangeLayout *send_layout,
    const struct SymmetryMpiRawResult *recv_result,
    uint64_t chunk_limit,
    int force_chunked)
{
  int peer;
  if (force_chunked != FALSE) return TRUE;
  for (peer = 0; peer < send_layout->nrank; peer++) {
    if (send_layout->counts[peer] > chunk_limit ||
        recv_result->counts[peer] > chunk_limit ||
        send_layout->displacements[peer] > (uint64_t)INT_MAX ||
        recv_result->displacements[peer] > (uint64_t)INT_MAX) {
      return TRUE;
    }
  }
  return FALSE;
}

static const void *const_entry_offset(const void *entries,
                                      uint64_t offset,
                                      size_t extent)
{
  return (const unsigned char *)entries + (size_t)offset * extent;
}

static void *entry_offset(void *entries, uint64_t offset, size_t extent)
{
  return (unsigned char *)entries + (size_t)offset * extent;
}

static int exchange_alltoallv(
    const void *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    struct SymmetryMpiRawResult *recv_result,
    MPI_Datatype payload_type)
{
  unsigned char send_dummy = 0U;
  unsigned char recv_dummy = 0U;
  const void *send_buffer =
      send_layout->count > 0U ? send_entries : &send_dummy;
  void *recv_buffer =
      recv_result->count > 0U ? recv_result->entries : &recv_dummy;
  int *send_counts = NULL;
  int *send_displacements = NULL;
  int *recv_counts = NULL;
  int *recv_displacements = NULL;
  size_t array_bytes;
  int peer;
  int ierr;

  if (SymmetryCheckedSizeMul((size_t)send_layout->nrank,
                             sizeof(*send_counts), &array_bytes) != 0) {
    return -1;
  }
  send_counts = (int *)malloc(array_bytes);
  send_displacements = (int *)malloc(array_bytes);
  recv_counts = (int *)malloc(array_bytes);
  recv_displacements = (int *)malloc(array_bytes);
  {
    int allocation_error =
        send_counts == NULL || send_displacements == NULL ||
        recv_counts == NULL || recv_displacements == NULL;
    if (SymmetryMpiAgreeError(TRUE, allocation_error) != 0) {
      free(send_counts);
      free(send_displacements);
      free(recv_counts);
      free(recv_displacements);
      return -1;
    }
  }
  for (peer = 0; peer < send_layout->nrank; peer++) {
    send_counts[peer] = (int)send_layout->counts[peer];
    send_displacements[peer] = (int)send_layout->displacements[peer];
    recv_counts[peer] = (int)recv_result->counts[peer];
    recv_displacements[peer] = (int)recv_result->displacements[peer];
  }

  ierr = MPI_Alltoallv(send_buffer, send_counts, send_displacements,
                       payload_type,
                       recv_buffer, recv_counts, recv_displacements,
                       payload_type, MPI_COMM_WORLD);
  free(send_counts);
  free(send_displacements);
  free(recv_counts);
  free(recv_displacements);
  return ierr == MPI_SUCCESS ? 0 : -1;
}

static int exchange_chunked(
    const void *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    int rank,
    uint64_t chunk_limit,
    struct SymmetryMpiRawResult *recv_result,
    struct SymmetryMpiExchangeStats *stats,
    const struct SymmetryMpiPayloadDescriptor *descriptor,
    MPI_Datatype payload_type)
{
  uint64_t *send_progress = NULL;
  uint64_t *recv_progress = NULL;
  MPI_Request *requests = NULL;
  size_t progress_bytes;
  size_t request_count_max;
  size_t request_bytes;
  int peer;
  int request_count;
  int ierr;

  if (send_layout->nrank > INT_MAX / 2 ||
      SymmetryCheckedSizeMul((size_t)send_layout->nrank,
                             sizeof(*send_progress), &progress_bytes) != 0 ||
      SymmetryCheckedSizeMul((size_t)send_layout->nrank, 2U,
                             &request_count_max) != 0 ||
      SymmetryCheckedSizeMul(request_count_max, sizeof(*requests),
                             &request_bytes) != 0) {
    return -1;
  }
  send_progress = (uint64_t *)malloc(progress_bytes);
  recv_progress = (uint64_t *)malloc(progress_bytes);
  requests = (MPI_Request *)malloc(request_bytes);
  {
    int allocation_error =
        send_progress == NULL || recv_progress == NULL || requests == NULL ||
        send_layout->counts[rank] != recv_result->counts[rank];
    if (SymmetryMpiAgreeError(TRUE, allocation_error) != 0) {
      free(send_progress);
      free(recv_progress);
      free(requests);
      return -1;
    }
  }
  memset(send_progress, 0, progress_bytes);
  memset(recv_progress, 0, progress_bytes);
  if (send_layout->counts[rank] > 0U) {
    memcpy(entry_offset(recv_result->entries,
                        recv_result->displacements[rank],
                        descriptor->extent),
           const_entry_offset(send_entries,
                              send_layout->displacements[rank],
                              descriptor->extent),
           (size_t)send_layout->counts[rank] * descriptor->extent);
    send_progress[rank] = send_layout->counts[rank];
    recv_progress[rank] = recv_result->counts[rank];
  }

  for (;;) {
    request_count = 0;
    for (peer = 0; peer < send_layout->nrank; peer++) {
      uint64_t remaining;
      uint64_t chunk;
      if (peer == rank) continue;
      remaining = recv_result->counts[peer] - recv_progress[peer];
      if (remaining == 0U) continue;
      chunk = remaining < chunk_limit ? remaining : chunk_limit;
      ierr = MPI_Irecv(
          entry_offset(recv_result->entries,
                       recv_result->displacements[peer] +
                           recv_progress[peer],
                       descriptor->extent),
          (int)chunk, payload_type, peer, descriptor->tag,
          MPI_COMM_WORLD, &requests[request_count]);
      if (ierr != MPI_SUCCESS) goto fail;
      request_count++;
      recv_progress[peer] += chunk;
      if (stats != NULL) stats->recv_messages++;
    }
    for (peer = 0; peer < send_layout->nrank; peer++) {
      uint64_t remaining;
      uint64_t chunk;
      if (peer == rank) continue;
      remaining = send_layout->counts[peer] - send_progress[peer];
      if (remaining == 0U) continue;
      chunk = remaining < chunk_limit ? remaining : chunk_limit;
      ierr = MPI_Isend(
          const_entry_offset(send_entries,
                             send_layout->displacements[peer] +
                                 send_progress[peer],
                             descriptor->extent),
          (int)chunk, payload_type, peer, descriptor->tag,
          MPI_COMM_WORLD, &requests[request_count]);
      if (ierr != MPI_SUCCESS) goto fail;
      request_count++;
      send_progress[peer] += chunk;
      if (stats != NULL) stats->send_messages++;
    }
    if (request_count == 0) break;
    ierr = MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
    if (ierr != MPI_SUCCESS) goto fail;
  }

  free(send_progress);
  free(recv_progress);
  free(requests);
  return 0;

fail:
  if (request_count > 0) {
    MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
  }
  free(send_progress);
  free(recv_progress);
  free(requests);
  return -1;
}
#endif

static int exchange_checked_payload(
    const void *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    const struct SymmetryMpiExchangeLayout *known_receive_layout,
    int rank,
    int nrank,
    const struct SymmetryMpiExchangeOptions *options,
    const struct SymmetryMpiExchangeMemoryBudget *budget,
    int output_is_empty,
    const struct SymmetryMpiPayloadDescriptor *descriptor,
    struct SymmetryMpiRawResult *result,
    struct SymmetryMpiExchangeStats *stats)
{
  struct SymmetryMpiRawResult next_result;
  struct SymmetryMpiExchangeStats next_stats;
  uint64_t recv_total = 0U;
  uint64_t chunk_limit = 0U;
  uint64_t message_byte_limit = 0U;
  size_t schedule_bytes = 0U;
  size_t result_bytes = 0U;
  size_t payload_bytes = 0U;
  size_t workspace_bytes = 0U;
  int force_chunked = FALSE;
  int mpi_active;
  int use_chunked = FALSE;
  int local_error = 0;
  int peer;
#ifdef MPI
  int comm_rank = 0;
  int comm_size = 1;
  uint64_t chunk_limit_min = 0U;
  uint64_t chunk_limit_max = 0U;
  int global_use_chunked = FALSE;
  int owns_payload_type = FALSE;
  MPI_Aint type_lower_bound = 0;
  MPI_Aint type_extent = 0;
  MPI_Datatype payload_type = MPI_DATATYPE_NULL;
#endif

  reset_raw_result(&next_result);
  reset_stats(&next_stats);
  reset_stats(stats);
  mpi_active = SymmetryMpiCollectivesActive();
  if (descriptor == NULL ||
      get_message_limits(
          descriptor != NULL ? descriptor->extent : 0U,
          &chunk_limit, &message_byte_limit) != 0 ||
      descriptor->tag < 0 || descriptor->tag > 32767) {
    local_error = 1;
  }
  if (result == NULL || output_is_empty == FALSE ||
      rank < 0 || rank >= nrank ||
      validate_layout(send_layout, send_entries, nrank,
                      descriptor != NULL ? descriptor->extent : 0U,
                      TRUE) != 0 ||
      (known_receive_layout != NULL &&
       validate_layout(known_receive_layout, NULL, nrank,
                       descriptor != NULL ? descriptor->extent : 0U,
                       FALSE) != 0)) {
    local_error = 1;
  }
  if (known_receive_layout != NULL &&
      send_layout != NULL &&
      rank >= 0 && rank < nrank &&
      send_layout->nrank == nrank &&
      known_receive_layout->nrank == nrank &&
      send_layout->counts != NULL &&
      known_receive_layout->counts != NULL &&
      send_layout->counts[rank] !=
          known_receive_layout->counts[rank]) {
    local_error = 1;
  }
  if (options != NULL) {
    if ((options->force_chunked != FALSE &&
         options->force_chunked != TRUE) ||
        options->chunk_limit > chunk_limit) {
      local_error = 1;
    } else {
      force_chunked = options->force_chunked;
      if (options->chunk_limit > 0U) chunk_limit = options->chunk_limit;
    }
  }
  if (budget != NULL &&
      (budget->byte_limit == 0U ||
       budget->live_bytes > budget->byte_limit)) {
    local_error = 1;
  }
  if (descriptor != NULL && chunk_limit > 0U) {
    message_byte_limit =
        chunk_limit * (uint64_t)descriptor->extent;
  }
#ifdef MPI
  if (mpi_active != FALSE) {
    if (MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank) != MPI_SUCCESS ||
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size) != MPI_SUCCESS ||
        comm_rank != rank || comm_size != nrank) {
      local_error = 1;
    }
  } else if (nrank != 1 || rank != 0) {
    local_error = 1;
  }
#else
  if (nrank != 1 || rank != 0) local_error = 1;
  (void)force_chunked;
#endif
  if (nrank < 1 ||
      SymmetryCheckedSizeMul((size_t)nrank, sizeof(*next_result.counts),
                             &schedule_bytes) != 0 ||
      SymmetryCheckedSizeAdd(schedule_bytes, schedule_bytes,
                             &result_bytes) != 0 ||
      checked_budget_peak(budget, result_bytes, 0U) != 0) {
    local_error = 1;
  }
  if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) return -1;

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    if (MPI_Allreduce(&chunk_limit, &chunk_limit_min, 1, MPI_UINT64_T,
                      MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&chunk_limit, &chunk_limit_max, 1, MPI_UINT64_T,
                      MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS) {
      return -1;
    }
    if (chunk_limit_min != chunk_limit_max) return -1;
  }
#endif

  next_result.counts = (uint64_t *)malloc(schedule_bytes);
  next_result.displacements = (uint64_t *)malloc(schedule_bytes);
  if (next_result.counts == NULL || next_result.displacements == NULL) {
    local_error = 1;
  }
  if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;

  if (known_receive_layout != NULL) {
    memcpy(next_result.counts, known_receive_layout->counts, schedule_bytes);
    memcpy(next_result.displacements,
           known_receive_layout->displacements, schedule_bytes);
  }
#ifdef MPI
  else if (mpi_active != FALSE && nrank > 1) {
    if (MPI_Alltoall(send_layout->counts, 1, MPI_UINT64_T,
                     next_result.counts, 1, MPI_UINT64_T,
                     MPI_COMM_WORLD) != MPI_SUCCESS) {
      goto fail;
    }
  }
#endif
  else {
    next_result.counts[0] = send_layout->counts[0];
    next_result.displacements[0] = 0U;
  }

  for (peer = 0; peer < nrank; peer++) {
    if (known_receive_layout != NULL &&
        next_result.displacements[peer] != recv_total) {
      local_error = 1;
      break;
    }
    if (known_receive_layout == NULL) {
      next_result.displacements[peer] = recv_total;
    }
    if (SymmetryCheckedU64Add(recv_total, next_result.counts[peer],
                              &recv_total) != 0) {
      local_error = 1;
      break;
    }
  }
  if (known_receive_layout != NULL &&
      recv_total != known_receive_layout->count) {
    local_error = 1;
  }
  if (descriptor == NULL ||
      recv_total > (uint64_t)(SIZE_MAX / descriptor->extent) ||
      SymmetryCheckedSizeMul((size_t)recv_total,
                             descriptor != NULL ? descriptor->extent : 0U,
                             &payload_bytes) != 0) {
    local_error = 1;
  }
  next_result.count = recv_total;
  next_result.nrank = nrank;

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    use_chunked =
        local_error == 0
            ? layout_requires_chunking(send_layout, &next_result, chunk_limit,
                                       force_chunked)
            : TRUE;
    if (MPI_Allreduce(&use_chunked, &global_use_chunked, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      goto fail;
    }
    use_chunked = global_use_chunked;
  }
#endif
  if (SymmetryCheckedSizeAdd(result_bytes, payload_bytes,
                             &result_bytes) != 0 ||
      exchange_workspace_bytes(nrank, mpi_active, use_chunked,
                               &workspace_bytes) != 0 ||
      checked_budget_peak(budget, result_bytes, workspace_bytes) != 0) {
    local_error = 1;
  }
  if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;
  if (payload_bytes > 0U) {
    next_result.entries = malloc(payload_bytes);
    if (next_result.entries == NULL) local_error = 1;
  }
  if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;

  next_stats.send_entries = send_layout->count;
  next_stats.recv_entries = recv_total;
  next_stats.message_entry_limit = chunk_limit;
  next_stats.message_byte_limit = message_byte_limit;

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    {
      uint64_t max_message_entries = 0U;
      for (peer = 0; peer < nrank; peer++) {
        uint64_t send_count;
        uint64_t recv_count;
        if (peer == rank) continue;
        send_count = send_layout->counts[peer];
        recv_count = next_result.counts[peer];
        if (use_chunked != FALSE) {
          if (send_count > chunk_limit) send_count = chunk_limit;
          if (recv_count > chunk_limit) recv_count = chunk_limit;
        }
        if (send_count > max_message_entries) {
          max_message_entries = send_count;
        }
        if (recv_count > max_message_entries) {
          max_message_entries = recv_count;
        }
      }
      next_stats.max_message_entries = max_message_entries;
      next_stats.max_message_bytes =
          max_message_entries * (uint64_t)descriptor->extent;
    }
    if (descriptor->create_type != NULL) {
      if (descriptor->create_type(&payload_type) != 0) {
        local_error = 1;
      } else {
        owns_payload_type = TRUE;
      }
    } else {
      payload_type = descriptor->builtin_type;
      if (payload_type == MPI_DATATYPE_NULL) local_error = 1;
    }
    if (local_error == 0 &&
        (MPI_Type_get_extent(payload_type, &type_lower_bound, &type_extent) !=
             MPI_SUCCESS ||
         type_lower_bound != 0 ||
         type_extent != (MPI_Aint)descriptor->extent)) {
      local_error = 1;
    }
    if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;
    if (use_chunked != FALSE) {
      if (exchange_chunked(send_entries, send_layout, rank, chunk_limit,
                           &next_result, &next_stats, descriptor,
                           payload_type) != 0) {
        goto fail;
      }
      next_stats.used_chunked = TRUE;
    } else {
      if (exchange_alltoallv(send_entries, send_layout, &next_result,
                             payload_type) != 0) {
        goto fail;
      }
    }
    if (owns_payload_type != FALSE) {
      MPI_Type_free(&payload_type);
      owns_payload_type = FALSE;
    }
  } else
#endif
  {
    if (send_layout->count != recv_total) goto fail;
    if (payload_bytes > 0U) {
      memcpy(next_result.entries, send_entries, payload_bytes);
    }
  }

  *result = next_result;
  if (stats != NULL) *stats = next_stats;
  return 0;

fail:
#ifdef MPI
  if (owns_payload_type != FALSE &&
      payload_type != MPI_DATATYPE_NULL) {
    MPI_Type_free(&payload_type);
  }
#endif
  free_raw_result(&next_result);
  return -1;
}

static int basis_result_is_empty(
    const struct SymmetryMpiExchangeResult *result)
{
  return result != NULL &&
         result->entries == NULL && result->count == 0U &&
         result->counts == NULL && result->displacements == NULL &&
         result->nrank == 0;
}

static int unsigned_long_result_is_empty(
    const struct SymmetryMpiUnsignedLongResult *result)
{
  return result != NULL &&
         result->entries == NULL && result->count == 0U &&
         result->counts == NULL && result->displacements == NULL &&
         result->nrank == 0;
}

static int lookup_response_result_is_empty(
    const struct SymmetryMpiLookupResponseResult *result)
{
  return result != NULL &&
         result->entries == NULL && result->count == 0U &&
         result->counts == NULL && result->displacements == NULL &&
         result->nrank == 0;
}

int SymmetryMpiExchangeBasisVectors(
    const struct SymmetryBasisVector *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    int rank,
    int nrank,
    const struct SymmetryMpiExchangeOptions *options,
    struct SymmetryMpiExchangeResult *result,
    struct SymmetryMpiExchangeStats *stats)
{
  struct SymmetryMpiPayloadDescriptor descriptor;
  struct SymmetryMpiRawResult raw_result;
  int status;
  memset(&descriptor, 0, sizeof(descriptor));
  descriptor.extent = sizeof(*send_entries);
  descriptor.tag = SYMMETRY_BASIS_EXCHANGE_TAG;
#ifdef MPI
  descriptor.builtin_type = MPI_DATATYPE_NULL;
  descriptor.create_type = SymmetryMpiCreateBasisVectorType;
#endif
  reset_raw_result(&raw_result);
  status = exchange_checked_payload(
      send_entries, send_layout, NULL, rank, nrank, options, NULL,
      basis_result_is_empty(result), &descriptor, &raw_result, stats);
  if (status != 0) return -1;
  result->entries = (struct SymmetryBasisVector *)raw_result.entries;
  result->count = raw_result.count;
  result->counts = raw_result.counts;
  result->displacements = raw_result.displacements;
  result->nrank = raw_result.nrank;
  return 0;
}

int SymmetryMpiExchangeUnsignedLongs(
    const unsigned long int *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    int rank,
    int nrank,
    const struct SymmetryMpiExchangeOptions *options,
    struct SymmetryMpiUnsignedLongResult *result,
    struct SymmetryMpiExchangeStats *stats)
{
  return SymmetryMpiExchangeUnsignedLongsWithBudget(
      send_entries, send_layout, rank, nrank, options, NULL, result, stats);
}

int SymmetryMpiExchangeUnsignedLongsWithBudget(
    const unsigned long int *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    int rank,
    int nrank,
    const struct SymmetryMpiExchangeOptions *options,
    const struct SymmetryMpiExchangeMemoryBudget *budget,
    struct SymmetryMpiUnsignedLongResult *result,
    struct SymmetryMpiExchangeStats *stats)
{
  struct SymmetryMpiPayloadDescriptor descriptor;
  struct SymmetryMpiRawResult raw_result;
  int status;
  memset(&descriptor, 0, sizeof(descriptor));
  descriptor.extent = sizeof(*send_entries);
  descriptor.tag = SYMMETRY_UNSIGNED_LONG_EXCHANGE_TAG;
#ifdef MPI
  descriptor.builtin_type = MPI_UNSIGNED_LONG;
  descriptor.create_type = NULL;
#endif
  reset_raw_result(&raw_result);
  status = exchange_checked_payload(
      send_entries, send_layout, NULL, rank, nrank, options, budget,
      unsigned_long_result_is_empty(result), &descriptor, &raw_result, stats);
  if (status != 0) return -1;
  result->entries = (unsigned long int *)raw_result.entries;
  result->count = raw_result.count;
  result->counts = raw_result.counts;
  result->displacements = raw_result.displacements;
  result->nrank = raw_result.nrank;
  return 0;
}

int SymmetryMpiExchangeLookupResponsesKnownLayout(
    const struct SymmetryMpiLookupResponse *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    const struct SymmetryMpiExchangeLayout *known_receive_layout,
    int rank,
    int nrank,
    const struct SymmetryMpiExchangeOptions *options,
    struct SymmetryMpiLookupResponseResult *result,
    struct SymmetryMpiExchangeStats *stats)
{
  return SymmetryMpiExchangeLookupResponsesKnownLayoutWithBudget(
      send_entries, send_layout, known_receive_layout,
      rank, nrank, options, NULL, result, stats);
}

int SymmetryMpiExchangeLookupResponsesKnownLayoutWithBudget(
    const struct SymmetryMpiLookupResponse *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    const struct SymmetryMpiExchangeLayout *known_receive_layout,
    int rank,
    int nrank,
    const struct SymmetryMpiExchangeOptions *options,
    const struct SymmetryMpiExchangeMemoryBudget *budget,
    struct SymmetryMpiLookupResponseResult *result,
    struct SymmetryMpiExchangeStats *stats)
{
  struct SymmetryMpiPayloadDescriptor descriptor;
  struct SymmetryMpiRawResult raw_result;
  int status;
  memset(&descriptor, 0, sizeof(descriptor));
  descriptor.extent = sizeof(*send_entries);
  descriptor.tag = SYMMETRY_LOOKUP_RESPONSE_EXCHANGE_TAG;
#ifdef MPI
  descriptor.builtin_type = MPI_DATATYPE_NULL;
  descriptor.create_type = SymmetryMpiCreateLookupResponseType;
#endif
  reset_raw_result(&raw_result);
  status = exchange_checked_payload(
      send_entries, send_layout, known_receive_layout,
      rank, nrank, options, budget, lookup_response_result_is_empty(result),
      &descriptor, &raw_result, stats);
  if (status != 0) return -1;
  result->entries =
      (struct SymmetryMpiLookupResponse *)raw_result.entries;
  result->count = raw_result.count;
  result->counts = raw_result.counts;
  result->displacements = raw_result.displacements;
  result->nrank = raw_result.nrank;
  return 0;
}

int SymmetryMpiExchangeUnsignedLongEchoesKnownLayoutWithBudget(
    const unsigned long int *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    const struct SymmetryMpiExchangeLayout *known_receive_layout,
    int rank,
    int nrank,
    const struct SymmetryMpiExchangeOptions *options,
    const struct SymmetryMpiExchangeMemoryBudget *budget,
    struct SymmetryMpiUnsignedLongResult *result,
    struct SymmetryMpiExchangeStats *stats)
{
  struct SymmetryMpiPayloadDescriptor descriptor;
  struct SymmetryMpiRawResult raw_result;
  int status;
  memset(&descriptor, 0, sizeof(descriptor));
  descriptor.extent = sizeof(*send_entries);
  descriptor.tag = SYMMETRY_UNSIGNED_LONG_ECHO_EXCHANGE_TAG;
#ifdef MPI
  descriptor.builtin_type = MPI_UNSIGNED_LONG;
  descriptor.create_type = NULL;
#endif
  reset_raw_result(&raw_result);
  status = exchange_checked_payload(
      send_entries, send_layout, known_receive_layout,
      rank, nrank, options, budget, unsigned_long_result_is_empty(result),
      &descriptor, &raw_result, stats);
  if (status != 0) return -1;
  result->entries = (unsigned long int *)raw_result.entries;
  result->count = raw_result.count;
  result->counts = raw_result.counts;
  result->displacements = raw_result.displacements;
  result->nrank = raw_result.nrank;
  return 0;
}
