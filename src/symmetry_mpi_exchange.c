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

static int get_message_limits(uint64_t *entry_limit, uint64_t *byte_limit)
{
  const uint64_t configured_bytes =
      (uint64_t)HPHI_SYMMETRY_EXCHANGE_MESSAGE_BYTES;
  uint64_t entries;
  if (entry_limit == NULL || byte_limit == NULL ||
      configured_bytes < (uint64_t)sizeof(struct SymmetryBasisVector)) {
    return -1;
  }
  entries = configured_bytes /
      (uint64_t)sizeof(struct SymmetryBasisVector);
  if (entries > (uint64_t)INT_MAX) entries = (uint64_t)INT_MAX;
  if (entries == 0U ||
      entries > UINT64_MAX /
          (uint64_t)sizeof(struct SymmetryBasisVector)) {
    return -1;
  }
  *entry_limit = entries;
  *byte_limit =
      entries * (uint64_t)sizeof(struct SymmetryBasisVector);
  return 0;
}

static int validate_layout(const struct SymmetryMpiExchangeLayout *layout,
                           const struct SymmetryBasisVector *entries,
                           int nrank)
{
  uint64_t total = 0U;
  int peer;
  if (layout == NULL || layout->nrank != nrank || nrank < 1 ||
      layout->counts == NULL || layout->displacements == NULL) {
    return -1;
  }
  for (peer = 0; peer < nrank; peer++) {
    if (layout->displacements[peer] != total ||
        SymmetryCheckedU64Add(total, layout->counts[peer], &total) != 0) {
      return -1;
    }
  }
  if (total != layout->count ||
      total > (uint64_t)(SIZE_MAX / sizeof(*entries)) ||
      (total > 0U && entries == NULL)) {
    return -1;
  }
  return 0;
}

static void reset_result(struct SymmetryMpiExchangeResult *result)
{
  if (result == NULL) return;
  result->entries = NULL;
  result->count = 0U;
  result->counts = NULL;
  result->displacements = NULL;
  result->nrank = 0;
}

static void reset_stats(struct SymmetryMpiExchangeStats *stats)
{
  if (stats == NULL) return;
  memset(stats, 0, sizeof(*stats));
}

void FreeSymmetryMpiExchangeResult(
    struct SymmetryMpiExchangeResult *result)
{
  if (result == NULL) return;
  free(result->entries);
  free(result->counts);
  free(result->displacements);
  reset_result(result);
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
      MPI_Get_address(&sample.stabilizer_size, &displacements[2]) != MPI_SUCCESS ||
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

static int layout_requires_chunking(
    const struct SymmetryMpiExchangeLayout *send_layout,
    const struct SymmetryMpiExchangeResult *recv_result,
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

static int exchange_alltoallv(
    const struct SymmetryBasisVector *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    struct SymmetryMpiExchangeResult *recv_result,
    MPI_Datatype vector_type)
{
  struct SymmetryBasisVector send_dummy;
  struct SymmetryBasisVector recv_dummy;
  const struct SymmetryBasisVector *send_buffer =
      send_layout->count > 0U ? send_entries : &send_dummy;
  struct SymmetryBasisVector *recv_buffer =
      recv_result->count > 0U ? recv_result->entries : &recv_dummy;
  int *send_counts = NULL;
  int *send_displacements = NULL;
  int *recv_counts = NULL;
  int *recv_displacements = NULL;
  size_t array_bytes;
  int peer;
  int ierr;

  if (SymmetryCheckedSizeMul((size_t)send_layout->nrank, sizeof(*send_counts),
                       &array_bytes) != 0) {
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
                       vector_type,
                       recv_buffer, recv_counts, recv_displacements,
                       vector_type, MPI_COMM_WORLD);
  free(send_counts);
  free(send_displacements);
  free(recv_counts);
  free(recv_displacements);
  return ierr == MPI_SUCCESS ? 0 : -1;
}

static int exchange_chunked(
    const struct SymmetryBasisVector *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    int rank,
    uint64_t chunk_limit,
    struct SymmetryMpiExchangeResult *recv_result,
    struct SymmetryMpiExchangeStats *stats,
    MPI_Datatype vector_type)
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
      SymmetryCheckedSizeMul((size_t)send_layout->nrank, sizeof(*send_progress),
                       &progress_bytes) != 0 ||
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
    memcpy(recv_result->entries + recv_result->displacements[rank],
           send_entries + send_layout->displacements[rank],
           (size_t)send_layout->counts[rank] * sizeof(*send_entries));
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
          recv_result->entries + recv_result->displacements[peer] +
              recv_progress[peer],
          (int)chunk, vector_type, peer, SYMMETRY_BASIS_EXCHANGE_TAG,
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
          send_entries + send_layout->displacements[peer] +
              send_progress[peer],
          (int)chunk, vector_type, peer, SYMMETRY_BASIS_EXCHANGE_TAG,
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
  if (request_count > 0) MPI_Waitall(request_count, requests,
                                     MPI_STATUSES_IGNORE);
  free(send_progress);
  free(recv_progress);
  free(requests);
  return -1;
}
#endif

int SymmetryMpiExchangeBasisVectors(
    const struct SymmetryBasisVector *send_entries,
    const struct SymmetryMpiExchangeLayout *send_layout,
    int rank,
    int nrank,
    const struct SymmetryMpiExchangeOptions *options,
    struct SymmetryMpiExchangeResult *result,
    struct SymmetryMpiExchangeStats *stats)
{
  struct SymmetryMpiExchangeResult next_result;
  uint64_t recv_total = 0U;
  uint64_t chunk_limit = 0U;
  uint64_t message_byte_limit = 0U;
  size_t schedule_bytes = 0U;
  int force_chunked = FALSE;
  int mpi_active;
  int local_error = 0;
  int peer;
#ifdef MPI
  int comm_rank = 0;
  int comm_size = 1;
  uint64_t chunk_limit_min = 0U;
  uint64_t chunk_limit_max = 0U;
  int use_chunked = FALSE;
  int global_use_chunked = FALSE;
  MPI_Datatype vector_type = MPI_DATATYPE_NULL;
#endif

  reset_result(&next_result);
  reset_stats(stats);
  mpi_active = SymmetryMpiCollectivesActive();
  if (get_message_limits(&chunk_limit, &message_byte_limit) != 0) {
    local_error = 1;
  }
  if (result == NULL ||
      (result != NULL &&
       (result->entries != NULL || result->count != 0U ||
        result->counts != NULL || result->displacements != NULL ||
        result->nrank != 0)) ||
      rank < 0 || rank >= nrank ||
      validate_layout(send_layout, send_entries, nrank) != 0) {
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
  if (chunk_limit > 0U) {
    message_byte_limit =
        chunk_limit * (uint64_t)sizeof(struct SymmetryBasisVector);
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
  if (SymmetryCheckedSizeMul((size_t)nrank, sizeof(*next_result.counts),
                       &schedule_bytes) != 0) {
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

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    if (MPI_Alltoall(send_layout->counts, 1, MPI_UINT64_T,
                     next_result.counts, 1, MPI_UINT64_T,
                     MPI_COMM_WORLD) != MPI_SUCCESS) {
      goto fail;
    }
  } else
#endif
  {
    next_result.counts[0] = send_layout->counts[0];
  }

  for (peer = 0; peer < nrank; peer++) {
    next_result.displacements[peer] = recv_total;
    if (SymmetryCheckedU64Add(recv_total, next_result.counts[peer],
                        &recv_total) != 0) {
      local_error = 1;
      break;
    }
  }
  if (recv_total > (uint64_t)(SIZE_MAX / sizeof(*next_result.entries))) {
    local_error = 1;
  } else if (recv_total > 0U) {
    next_result.entries = (struct SymmetryBasisVector *)malloc(
        (size_t)recv_total * sizeof(*next_result.entries));
    if (next_result.entries == NULL) local_error = 1;
  }
  if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;
  next_result.count = recv_total;
  next_result.nrank = nrank;

  if (stats != NULL) {
    stats->send_entries = send_layout->count;
    stats->recv_entries = recv_total;
    stats->message_entry_limit = chunk_limit;
    stats->message_byte_limit = message_byte_limit;
  }

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    use_chunked =
        layout_requires_chunking(send_layout, &next_result, chunk_limit,
                                 force_chunked);
    if (MPI_Allreduce(&use_chunked, &global_use_chunked, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      goto fail;
    }
    use_chunked = global_use_chunked;
    if (stats != NULL) {
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
      stats->max_message_entries = max_message_entries;
      stats->max_message_bytes =
          max_message_entries *
          (uint64_t)sizeof(struct SymmetryBasisVector);
    }
    if (SymmetryMpiCreateBasisVectorType(&vector_type) != 0) local_error = 1;
    if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;
    if (use_chunked != FALSE) {
      if (exchange_chunked(send_entries, send_layout, rank, chunk_limit,
                           &next_result, stats, vector_type) != 0) {
        goto fail;
      }
      if (stats != NULL) stats->used_chunked = TRUE;
    } else {
      if (exchange_alltoallv(send_entries, send_layout, &next_result,
                             vector_type) != 0) {
        goto fail;
      }
    }
    MPI_Type_free(&vector_type);
  } else
#endif
  {
    if (send_layout->count != recv_total) goto fail;
    if (recv_total > 0U) {
      memcpy(next_result.entries, send_entries,
             (size_t)recv_total * sizeof(*send_entries));
    }
  }

  *result = next_result;
  return 0;

fail:
#ifdef MPI
  if (vector_type != MPI_DATATYPE_NULL) MPI_Type_free(&vector_type);
#endif
  FreeSymmetryMpiExchangeResult(&next_result);
  return -1;
}
