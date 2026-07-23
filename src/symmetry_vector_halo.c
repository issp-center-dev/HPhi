#include <limits.h>
#include <stdint.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "DefCommon.h"
struct BindStruct;
#include "CalcTime.h"
#include "symmetry_vector_halo.h"

#ifndef HPHI_SYMMETRY_HALO_BITSET_CAP_BYTES
#define HPHI_SYMMETRY_HALO_BITSET_CAP_BYTES (64U * 1024U * 1024U)
#endif

static int checked_size_add(size_t lhs, size_t rhs, size_t *result)
{
  if (result == NULL || lhs > SIZE_MAX - rhs) return -1;
  *result = lhs + rhs;
  return 0;
}

static int checked_size_mul(size_t lhs, size_t rhs, size_t *result)
{
  if (result == NULL || (lhs != 0U && rhs > SIZE_MAX / lhs)) return -1;
  *result = lhs * rhs;
  return 0;
}

static int halo_mpi_collectives_active(void)
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

static int agree_halo_error(int mpi_active, int local_error)
{
#ifdef MPI
  int global_error = local_error;
  if (mpi_active != FALSE &&
      MPI_Allreduce(&local_error, &global_error, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
  return global_error;
#else
  (void)mpi_active;
  return local_error;
#endif
}

static unsigned long long hash_bytes(unsigned long long hash,
                                     const void *data,
                                     size_t size)
{
  const unsigned char *bytes = (const unsigned char *)data;
  size_t index;
  for (index = 0U; index < size; index++) {
    hash ^= (unsigned long long)bytes[index];
    hash *= 1099511628211ULL;
  }
  return hash;
}

static unsigned long long halo_schedule_checksum(
    const struct SymmetryVectorHaloPlan *halo)
{
  unsigned long long hash = 14695981039346656037ULL;
  hash = hash_bytes(hash, &halo->dim, sizeof(halo->dim));
  hash = hash_bytes(hash, &halo->local_offset, sizeof(halo->local_offset));
  hash = hash_bytes(hash, &halo->local_dim, sizeof(halo->local_dim));
  hash = hash_bytes(hash, &halo->nrank, sizeof(halo->nrank));
  hash = hash_bytes(hash, &halo->rank, sizeof(halo->rank));
  hash = hash_bytes(hash, &halo->ghost_count, sizeof(halo->ghost_count));
  hash = hash_bytes(hash, &halo->send_value_count,
                    sizeof(halo->send_value_count));
  if (halo->nrank > 0) {
    size_t count_bytes = (size_t)halo->nrank * sizeof(*halo->recv_counts);
    hash = hash_bytes(hash, halo->recv_counts, count_bytes);
    hash = hash_bytes(hash, halo->recv_displs, count_bytes);
    hash = hash_bytes(hash, halo->send_counts, count_bytes);
    hash = hash_bytes(hash, halo->send_displs, count_bytes);
  }
  if (halo->ghost_count > 0U) {
    hash = hash_bytes(hash, halo->ghost_global_index,
                      halo->ghost_count * sizeof(*halo->ghost_global_index));
  }
  if (halo->send_value_count > 0U) {
    hash = hash_bytes(hash, halo->send_local_index,
                      halo->send_value_count *
                          sizeof(*halo->send_local_index));
  }
  return hash;
}

int SymmetryVectorOwnerOfGlobalIndex(unsigned long int dim,
                                     int nrank,
                                     unsigned long int global_index)
{
  unsigned long int base;
  unsigned long int remainder;
  unsigned long int zero_index;
  unsigned long int large_block_end;
  if (nrank < 1 || global_index == 0UL || global_index > dim) return -1;
  base = dim / (unsigned long int)nrank;
  remainder = dim % (unsigned long int)nrank;
  zero_index = global_index - 1UL;
  large_block_end =
      remainder == 0UL ? 0UL : (base + 1UL) * remainder;
  if (zero_index < large_block_end) {
    return (int)(zero_index / (base + 1UL));
  }
  if (base == 0UL) return -1;
  return (int)(remainder + (zero_index - large_block_end) / base);
}

void FreeSymmetryVectorHaloPlan(struct SymmetryVectorHaloPlan *halo)
{
  if (halo == NULL) return;
  free(halo->send_counts);
  free(halo->send_displs);
  free(halo->recv_counts);
  free(halo->recv_displs);
  free(halo->ghost_global_index);
  free(halo->send_local_index);
  free(halo->send_values);
  free(halo->ghost_values);
  memset(halo, 0, sizeof(*halo));
}

int BuildSymmetryVectorHaloPlan(struct SymmetryVectorHaloPlan *halo,
                                unsigned long int dim,
                                unsigned long int local_offset,
                                unsigned long int local_dim,
                                const unsigned long int *global_columns,
                                size_t column_count,
                                int nrank,
                                int rank,
                                size_t *local_column_count,
                                size_t *remote_column_count)
{
  unsigned char *remote_bits = NULL;
  unsigned long long *request_counts = NULL;
  unsigned long long *serve_counts = NULL;
  unsigned long int *requested_global_index = NULL;
#ifdef MPI
  unsigned long int request_dummy = 0UL;
#endif
  unsigned long int first_local;
  unsigned long int last_local;
  unsigned long int expected_local_offset;
  unsigned long int expected_local_dim;
  unsigned long int block_base;
  unsigned long int block_remainder;
  unsigned long int unsigned_rank;
  unsigned long int window_capacity;
  unsigned long int window_zero;
  unsigned long int window_dim;
  size_t bitset_bytes = 0U;
  size_t window_bytes = 0U;
  size_t count64_bytes = 0U;
  size_t count_int_bytes = 0U;
  size_t request_bytes = 0U;
  size_t index_bytes = 0U;
  size_t schedule_scratch_bytes = 0U;
  size_t total_count;
  size_t column;
  size_t ghost_position = 0U;
  size_t byte_index;
  int mpi_active;
  int local_error = 0;
  int owner;
  int peer;

  mpi_active = halo_mpi_collectives_active();
  if (halo == NULL || local_column_count == NULL ||
      remote_column_count == NULL || nrank < 1 || rank < 0 || rank >= nrank ||
      local_offset > dim || local_dim > dim - local_offset ||
      (column_count > 0U && global_columns == NULL)) {
    local_error = 1;
  }
#ifdef MPI
  if (mpi_active != FALSE) {
    int comm_size = 0;
    if (MPI_Comm_size(MPI_COMM_WORLD, &comm_size) != MPI_SUCCESS ||
        comm_size != nrank) {
      local_error = 1;
    }
  }
#endif
  if (agree_halo_error(mpi_active, local_error) != 0) return -1;

  block_base = dim / (unsigned long int)nrank;
  block_remainder = dim % (unsigned long int)nrank;
  unsigned_rank = (unsigned long int)rank;
  expected_local_dim =
      block_base + (unsigned_rank < block_remainder ? 1UL : 0UL);
  expected_local_offset =
      block_base * unsigned_rank +
      (unsigned_rank < block_remainder ? unsigned_rank : block_remainder);
  if (local_offset != expected_local_offset ||
      local_dim != expected_local_dim) {
    local_error = 1;
  }
  if (agree_halo_error(mpi_active, local_error) != 0) return -1;

  memset(halo, 0, sizeof(*halo));
  *local_column_count = 0U;
  *remote_column_count = 0U;
  halo->nrank = nrank;
  halo->rank = rank;
  halo->dim = dim;
  halo->local_offset = local_offset;
  halo->local_dim = local_dim;
  first_local = local_offset + 1UL;
  last_local = local_offset + local_dim;

  StartTimer(1130);
  for (column = 0U; column < column_count; column++) {
    unsigned long int global_index = global_columns[column];
    if (global_index == 0UL || global_index > dim) {
      local_error = 1;
      break;
    }
    if (local_dim > 0UL &&
        global_index >= first_local && global_index <= last_local) {
      if (*local_column_count == SIZE_MAX) {
        local_error = 1;
        break;
      }
      (*local_column_count)++;
    } else {
      if (*remote_column_count == SIZE_MAX) {
        local_error = 1;
        break;
      }
      (*remote_column_count)++;
    }
  }
  if (agree_halo_error(mpi_active, local_error) != 0) goto fail_topology;
  if (*remote_column_count > 0U) {
    size_t bitset_cap = (size_t)HPHI_SYMMETRY_HALO_BITSET_CAP_BYTES;
    if (bitset_cap == 0U) {
      local_error = 1;
    } else if (bitset_cap > ULONG_MAX / 8UL) {
      window_capacity = ULONG_MAX;
    } else {
      window_capacity = (unsigned long int)bitset_cap * 8UL;
    }
    if (local_error == 0) {
      window_dim = dim < window_capacity ? dim : window_capacity;
      bitset_bytes = (size_t)(window_dim / 8UL);
      if (window_dim % 8UL != 0UL) {
        bitset_bytes++;
      }
    }
    if (local_error == 0) {
      remote_bits = (unsigned char *)calloc(bitset_bytes,
                                             sizeof(*remote_bits));
      if (remote_bits == NULL) local_error = 1;
    }
  }
  if ((size_t)nrank > SIZE_MAX / sizeof(*request_counts)) {
    local_error = 1;
  } else {
    request_counts = (unsigned long long *)calloc(
        (size_t)nrank, sizeof(*request_counts));
    serve_counts = (unsigned long long *)calloc(
        (size_t)nrank, sizeof(*serve_counts));
    halo->recv_counts = (int *)calloc((size_t)nrank,
                                      sizeof(*halo->recv_counts));
    halo->recv_displs = (int *)calloc((size_t)nrank,
                                      sizeof(*halo->recv_displs));
    halo->send_counts = (int *)calloc((size_t)nrank,
                                      sizeof(*halo->send_counts));
    halo->send_displs = (int *)calloc((size_t)nrank,
                                      sizeof(*halo->send_displs));
    if (request_counts == NULL || serve_counts == NULL ||
        halo->recv_counts == NULL || halo->recv_displs == NULL ||
        halo->send_counts == NULL || halo->send_displs == NULL) {
      local_error = 1;
    }
  }
  if (checked_size_mul(2U * sizeof(*request_counts), (size_t)nrank,
                       &count64_bytes) != 0 ||
      checked_size_mul(4U * sizeof(*halo->recv_counts), (size_t)nrank,
                       &count_int_bytes) != 0) {
    local_error = 1;
  }
  if (agree_halo_error(mpi_active, local_error) != 0) goto fail_topology;

  window_zero = 0UL;
  while (window_zero < dim && *remote_column_count > 0U) {
    unsigned long int remaining = dim - window_zero;
    window_dim =
        remaining < window_capacity ? remaining : window_capacity;
    window_bytes = (size_t)(window_dim / 8UL);
    if (window_dim % 8UL != 0UL) window_bytes++;
    memset(remote_bits, 0, bitset_bytes);
    for (column = 0U; column < column_count; column++) {
      unsigned long int global_index = global_columns[column];
      if (local_dim == 0UL ||
          global_index < first_local || global_index > last_local) {
        unsigned long int zero_index = global_index - 1UL;
        if (zero_index >= window_zero &&
            zero_index - window_zero < window_dim) {
          unsigned long int relative_index = zero_index - window_zero;
          remote_bits[(size_t)(relative_index / 8UL)] |=
              (unsigned char)(1U <<
                              (unsigned int)(relative_index % 8UL));
        }
      }
    }
    for (byte_index = 0U; byte_index < window_bytes; byte_index++) {
      unsigned char bits = remote_bits[byte_index];
      unsigned int bit;
      for (bit = 0U; bits != 0U && bit < 8U; bit++) {
        if ((bits & (unsigned char)(1U << bit)) == 0U) continue;
        {
          unsigned long int global_index =
              window_zero + (unsigned long int)byte_index * 8UL +
              (unsigned long int)bit + 1UL;
          owner = SymmetryVectorOwnerOfGlobalIndex(
              dim, nrank, global_index);
          if (owner < 0 || owner == rank ||
              request_counts[owner] == ULLONG_MAX ||
              halo->ghost_count == SIZE_MAX) {
            local_error = 1;
            break;
          }
          request_counts[owner]++;
          halo->ghost_count++;
        }
      }
      if (local_error != 0) break;
    }
    if (local_error != 0) break;
    window_zero += window_dim;
  }
  if (agree_halo_error(mpi_active, local_error) != 0) goto fail_topology;

  if (halo->ghost_count > 0U) {
    if (halo->ghost_count >
        SIZE_MAX / sizeof(*halo->ghost_global_index)) {
      local_error = 1;
    } else {
      halo->ghost_global_index = (unsigned long int *)malloc(
          halo->ghost_count * sizeof(*halo->ghost_global_index));
      if (halo->ghost_global_index == NULL) local_error = 1;
    }
  }
  if (agree_halo_error(mpi_active, local_error) != 0) goto fail_topology;

  window_zero = 0UL;
  while (window_zero < dim && *remote_column_count > 0U) {
    unsigned long int remaining = dim - window_zero;
    window_dim =
        remaining < window_capacity ? remaining : window_capacity;
    window_bytes = (size_t)(window_dim / 8UL);
    if (window_dim % 8UL != 0UL) window_bytes++;
    memset(remote_bits, 0, bitset_bytes);
    for (column = 0U; column < column_count; column++) {
      unsigned long int global_index = global_columns[column];
      if (local_dim == 0UL ||
          global_index < first_local || global_index > last_local) {
        unsigned long int zero_index = global_index - 1UL;
        if (zero_index >= window_zero &&
            zero_index - window_zero < window_dim) {
          unsigned long int relative_index = zero_index - window_zero;
          remote_bits[(size_t)(relative_index / 8UL)] |=
              (unsigned char)(1U <<
                              (unsigned int)(relative_index % 8UL));
        }
      }
    }
    for (byte_index = 0U; byte_index < window_bytes; byte_index++) {
      unsigned char bits = remote_bits[byte_index];
      unsigned int bit;
      for (bit = 0U; bits != 0U && bit < 8U; bit++) {
        if ((bits & (unsigned char)(1U << bit)) != 0U) {
          halo->ghost_global_index[ghost_position++] =
              window_zero + (unsigned long int)byte_index * 8UL +
              (unsigned long int)bit + 1UL;
        }
      }
    }
    window_zero += window_dim;
  }
  if (ghost_position != halo->ghost_count) local_error = 1;

  total_count = 0U;
  for (peer = 0; peer < nrank; peer++) {
    if (request_counts[peer] > (unsigned long long)INT_MAX ||
        total_count > (size_t)INT_MAX -
                          (size_t)request_counts[peer]) {
      local_error = 1;
      break;
    }
    halo->recv_counts[peer] = (int)request_counts[peer];
    halo->recv_displs[peer] = (int)total_count;
    total_count += (size_t)request_counts[peer];
    if (request_counts[peer] > 0ULL) {
      halo->incoming_peer_count++;
      if ((size_t)request_counts[peer] > halo->max_recv_from_peer) {
        halo->max_recv_from_peer = (size_t)request_counts[peer];
      }
    }
  }
  if (total_count != halo->ghost_count) local_error = 1;
  if (checked_size_add(bitset_bytes, count64_bytes,
                       &halo->topology_scratch_bytes) != 0) {
    local_error = 1;
  }
  if (agree_halo_error(mpi_active, local_error) != 0) goto fail_topology;
  halo->request_layout_ready = TRUE;
  free(remote_bits);
  remote_bits = NULL;
  StopTimer(1130);

  StartTimer(1131);
  if (nrank > 1 && mpi_active == FALSE) {
    halo->schedule_checksum = halo_schedule_checksum(halo);
    free(request_counts);
    free(serve_counts);
    StopTimer(1131);
    return 0;
  }
#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    local_error =
        MPI_Alltoall(request_counts, 1, MPI_UNSIGNED_LONG_LONG,
                     serve_counts, 1, MPI_UNSIGNED_LONG_LONG,
                     MPI_COMM_WORLD) == MPI_SUCCESS ? 0 : 1;
  }
#endif
  if (agree_halo_error(mpi_active, local_error) != 0) goto fail_schedule;

  total_count = 0U;
  for (peer = 0; peer < nrank; peer++) {
    if (serve_counts[peer] > (unsigned long long)INT_MAX ||
        total_count > (size_t)INT_MAX - (size_t)serve_counts[peer]) {
      local_error = 1;
      break;
    }
    halo->send_counts[peer] = (int)serve_counts[peer];
    halo->send_displs[peer] = (int)total_count;
    total_count += (size_t)serve_counts[peer];
    if (serve_counts[peer] > 0ULL) {
      halo->outgoing_peer_count++;
      if ((size_t)serve_counts[peer] > halo->max_send_to_peer) {
        halo->max_send_to_peer = (size_t)serve_counts[peer];
      }
    }
  }
  halo->send_value_count = total_count;
  if (agree_halo_error(mpi_active, local_error) != 0) goto fail_schedule;

  if (halo->send_value_count > 0U) {
    if (halo->send_value_count >
        SIZE_MAX / sizeof(*requested_global_index)) {
      local_error = 1;
    } else {
      requested_global_index = (unsigned long int *)malloc(
          halo->send_value_count * sizeof(*requested_global_index));
      halo->send_local_index = (unsigned long int *)malloc(
          halo->send_value_count * sizeof(*halo->send_local_index));
      if (requested_global_index == NULL ||
          halo->send_local_index == NULL) {
        local_error = 1;
      }
    }
  }
  if (agree_halo_error(mpi_active, local_error) != 0) goto fail_schedule;

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    unsigned long int *request_send =
        halo->ghost_count > 0U ? halo->ghost_global_index : &request_dummy;
    unsigned long int *request_recv =
        halo->send_value_count > 0U ? requested_global_index : &request_dummy;
    local_error =
        MPI_Alltoallv(request_send, halo->recv_counts, halo->recv_displs,
                      MPI_UNSIGNED_LONG,
                      request_recv, halo->send_counts, halo->send_displs,
                      MPI_UNSIGNED_LONG, MPI_COMM_WORLD) == MPI_SUCCESS
            ? 0
            : 1;
  }
#endif
  if (agree_halo_error(mpi_active, local_error) != 0) goto fail_schedule;

  for (column = 0U; column < halo->send_value_count; column++) {
    unsigned long int global_index = requested_global_index[column];
    owner = SymmetryVectorOwnerOfGlobalIndex(dim, nrank, global_index);
    if (owner != rank || global_index <= local_offset ||
        global_index > local_offset + local_dim) {
      local_error = 1;
      break;
    }
    halo->send_local_index[column] = global_index - local_offset;
  }
  if (agree_halo_error(mpi_active, local_error) != 0) goto fail_schedule;

#ifdef MPI
  if (mpi_active != FALSE && nrank > 1) {
    unsigned long long local_ghost_count =
        (unsigned long long)halo->ghost_count;
    unsigned long long local_send_count =
        (unsigned long long)halo->send_value_count;
    unsigned long long ghost_sum = 0ULL;
    unsigned long long send_sum = 0ULL;
    if (MPI_Allreduce(&local_ghost_count, &ghost_sum, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&local_send_count, &send_sum, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS ||
        ghost_sum != send_sum) {
      local_error = 1;
    }
  }
#endif
  if (agree_halo_error(mpi_active, local_error) != 0) goto fail_schedule;

  if (halo->send_value_count > 0U) {
    halo->send_values = (double complex *)malloc(
        halo->send_value_count * sizeof(*halo->send_values));
    if (halo->send_values == NULL) local_error = 1;
  }
  if (halo->ghost_count > 0U) {
    halo->ghost_values = (double complex *)malloc(
        halo->ghost_count * sizeof(*halo->ghost_values));
    if (halo->ghost_values == NULL) local_error = 1;
  }
  if (checked_size_add(halo->ghost_count, halo->send_value_count,
                       &total_count) != 0 ||
      checked_size_mul(total_count, sizeof(unsigned long int),
                       &index_bytes) != 0 ||
      checked_size_add(count_int_bytes, index_bytes,
                       &halo->schedule_bytes) != 0 ||
      checked_size_mul(total_count, sizeof(double complex),
                       &halo->runtime_buffer_bytes) != 0 ||
      checked_size_mul(halo->send_value_count,
                       sizeof(*requested_global_index),
                       &request_bytes) != 0 ||
      checked_size_add(count64_bytes, request_bytes,
                       &schedule_scratch_bytes) != 0) {
    local_error = 1;
  }
  if (schedule_scratch_bytes > halo->topology_scratch_bytes) {
    halo->topology_scratch_bytes = schedule_scratch_bytes;
  }
  if (agree_halo_error(mpi_active, local_error) != 0) goto fail_schedule;
  halo->ready = TRUE;
  halo->schedule_checksum = halo_schedule_checksum(halo);
  free(requested_global_index);
  free(request_counts);
  free(serve_counts);
  StopTimer(1131);
  return 0;

fail_schedule:
  free(requested_global_index);
  free(request_counts);
  free(serve_counts);
  StopTimer(1131);
  FreeSymmetryVectorHaloPlan(halo);
  return -1;

fail_topology:
  free(remote_bits);
  free(request_counts);
  free(serve_counts);
  StopTimer(1130);
  FreeSymmetryVectorHaloPlan(halo);
  return -1;
}

int ExchangeSymmetryVectorHaloReference(
    struct SymmetryVectorHaloPlan *halo,
    const double complex *local_vector,
    const double complex *full_vector)
{
  double complex send_dummy = 0.0;
  double complex recv_dummy = 0.0;
  size_t index;
  int mpi_active;
  int local_error = 0;

  mpi_active = halo_mpi_collectives_active();
  if (halo == NULL || halo->ready != TRUE || local_vector == NULL ||
      full_vector == NULL ||
      (halo != NULL && halo->nrank > 1 && mpi_active == FALSE)) {
    local_error = 1;
  }
#ifdef MPI
  if (mpi_active != FALSE && halo != NULL) {
    int comm_rank = -1;
    int comm_size = 0;
    if (MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank) != MPI_SUCCESS ||
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size) != MPI_SUCCESS ||
        comm_rank != halo->rank || comm_size != halo->nrank) {
      local_error = 1;
    }
  }
#endif
  if (agree_halo_error(mpi_active, local_error) != 0) return -1;

  StartTimer(1510);
  for (index = 0U; index < halo->send_value_count; index++) {
    unsigned long int local_index = halo->send_local_index[index];
    if (local_index == 0UL || local_index > halo->local_dim) {
      local_error = 1;
      break;
    }
    halo->send_values[index] = local_vector[local_index];
  }
  StopTimer(1510);
  if (agree_halo_error(mpi_active, local_error) != 0) return -1;

  StartTimer(1511);
#ifdef MPI
  if (mpi_active != FALSE && halo->nrank > 1) {
    double complex *send_buffer =
        halo->send_value_count > 0U ? halo->send_values : &send_dummy;
    double complex *recv_buffer =
        halo->ghost_count > 0U ? halo->ghost_values : &recv_dummy;
    local_error =
        MPI_Alltoallv(send_buffer, halo->send_counts, halo->send_displs,
                      MPI_DOUBLE_COMPLEX,
                      recv_buffer, halo->recv_counts, halo->recv_displs,
                      MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD) == MPI_SUCCESS
            ? 0
            : 1;
  }
#else
  (void)send_dummy;
  (void)recv_dummy;
#endif
  StopTimer(1511);
  if (agree_halo_error(mpi_active, local_error) != 0) return -1;

  StartTimer(1512);
  for (index = 0U; index < halo->ghost_count; index++) {
    unsigned long int global_index = halo->ghost_global_index[index];
    if (memcmp(&halo->ghost_values[index], &full_vector[global_index],
               sizeof(halo->ghost_values[index])) != 0) {
      local_error = 1;
      break;
    }
  }
  StopTimer(1512);
  if (agree_halo_error(mpi_active, local_error) != 0) return -1;
  halo->reference_exchange_calls++;
  return 0;
}
