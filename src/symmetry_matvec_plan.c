#include <limits.h>
#include <stdint.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "DefCommon.h"
#include "bitcalc.h"
#include "global.h"
#include "struct.h"
#include "CalcTime.h"
#include "symmetry_basis.h"
#include "symmetry_matvec_plan.h"
#include "wrapperMPI.h"

static int apply_exchange_halfspin(unsigned long int state,
                                   int site0,
                                   int site1,
                                   unsigned long int *out_state)
{
  unsigned long int b0 = (state >> (unsigned int)site0) & 1UL;
  unsigned long int b1 = (state >> (unsigned int)site1) & 1UL;
  if (b0 == b1) return FALSE;
  *out_state = state ^ (1UL << (unsigned int)site0) ^ (1UL << (unsigned int)site1);
  return TRUE;
}

static unsigned long int mask_between_sites(unsigned int site0, unsigned int site1)
{
  unsigned int lo = site0 < site1 ? site0 : site1;
  unsigned int hi = site0 < site1 ? site1 : site0;
  if (hi <= lo + 1U) return 0UL;
  return (1UL << hi) - (1UL << (lo + 1U));
}

static int apply_spinless_hopping_hermite(unsigned long int state,
                                          unsigned int site1,
                                          unsigned int site2,
                                          double complex trans,
                                          unsigned long int *out_state,
                                          double complex *hval)
{
  unsigned long int mask1 = 1UL << site1;
  unsigned long int mask2 = 1UL << site2;
  unsigned long int occupied1 = state & mask1;
  unsigned long int occupied2 = state & mask2;
  int sgn = 1;
  if (site1 == site2) return FALSE;
  if ((occupied1 == 0UL && occupied2 == 0UL) ||
      (occupied1 != 0UL && occupied2 != 0UL)) {
    return FALSE;
  }
  SgnBit(state & mask_between_sites(site1, site2), &sgn);
  *out_state = state ^ mask1 ^ mask2;
  if (occupied1 != 0UL && occupied2 == 0UL) {
    *hval = (double)sgn * conj(trans);
  } else {
    *hval = (double)sgn * trans;
  }
  return TRUE;
}

static int apply_hubbard_hopping_hermite(unsigned long int state,
                                         unsigned int site1,
                                         unsigned int spin1,
                                         unsigned int site2,
                                         unsigned int spin2,
                                         double complex trans,
                                         unsigned long int *out_state,
                                         double complex *hval)
{
  const unsigned int max_bits = (unsigned int)(sizeof(unsigned long int) * CHAR_BIT);
  unsigned int orbital1;
  unsigned int orbital2;
  if (spin1 > 1U || spin2 > 1U) return FALSE;
  if (site1 > (max_bits - spin1) / 2U ||
      site2 > (max_bits - spin2) / 2U) {
    return FALSE;
  }
  orbital1 = 2U * site1 + spin1;
  orbital2 = 2U * site2 + spin2;
  if (orbital1 >= max_bits || orbital2 >= max_bits) return FALSE;
  return apply_spinless_hopping_hermite(state, orbital1, orbital2, trans,
                                        out_state, hval);
}

static int emit_canonicalized_transition(const struct BindStruct *X,
                                         unsigned long int beta,
                                         unsigned long int to_state,
                                         double complex hval,
                                         SymmetryEntryCallback callback,
                                         void *context)
{
  double norm_factor;
  double complex coefficient;
  struct SymmetryCanonicalResult result;
  if (SymmetryCanonicalizeState(X, to_state, &result) != 0) return -1;
  if (result.found != TRUE) return 0;
  if (X->Sym->basis[beta].norm == 0.0) return -1;
  norm_factor = X->Sym->basis[result.basis_index].norm /
                X->Sym->basis[beta].norm;
  coefficient = hval * result.phase * norm_factor;
  return callback(result.basis_index, coefficient, context);
}

int SymmetryEnumerateColumn(const struct BindStruct *X,
                            unsigned long int beta,
                            SymmetryEntryCallback callback,
                            void *context)
{
  unsigned int p;
  if (X == NULL || X->Sym == NULL || callback == NULL || beta == 0UL ||
      beta > X->Sym->dim || X->Sym->sym_diagonal == NULL) {
    return -1;
  }

  if (callback(beta, X->Sym->sym_diagonal[beta], context) != 0) return -1;

  if (X->Def.iCalcModel == Spin) {
    for (p = 0; p < X->Def.NExchangeCoupling; p++) {
      unsigned long int out_state;
      if (apply_exchange_halfspin(X->Sym->basis[beta].rep_state,
                                  X->Def.ExchangeCoupling[p][0],
                                  X->Def.ExchangeCoupling[p][1],
                                  &out_state) == TRUE &&
          emit_canonicalized_transition(X, beta, out_state,
                                        X->Def.ParaExchangeCoupling[p],
                                        callback, context) != 0) {
        return -1;
      }
    }
    return 0;
  }

  if (X->Def.iCalcModel == SpinlessFermion) {
    for (p = 0; p < X->Def.EDNTransfer; p += 2U) {
      unsigned long int out_state;
      double complex hval;
      double complex trans = -X->Def.EDParaGeneralTransfer[p];
      unsigned int site1 = (unsigned int)X->Def.EDGeneralTransfer[p][0];
      unsigned int site2 = (unsigned int)X->Def.EDGeneralTransfer[p][2];
      if (apply_spinless_hopping_hermite(X->Sym->basis[beta].rep_state,
                                         site1, site2, trans,
                                         &out_state, &hval) == TRUE &&
          emit_canonicalized_transition(X, beta, out_state, hval,
                                        callback, context) != 0) {
        return -1;
      }
    }
    return 0;
  }

  if (X->Def.iCalcModel == Hubbard) {
    for (p = 0; p < X->Def.EDNTransfer; p += 2U) {
      unsigned long int out_state;
      double complex hval;
      double complex trans = -X->Def.EDParaGeneralTransfer[p];
      unsigned int site1 = (unsigned int)X->Def.EDGeneralTransfer[p][0];
      unsigned int spin1 = (unsigned int)X->Def.EDGeneralTransfer[p][1];
      unsigned int site2 = (unsigned int)X->Def.EDGeneralTransfer[p][2];
      unsigned int spin2 = (unsigned int)X->Def.EDGeneralTransfer[p][3];
      if (apply_hubbard_hopping_hermite(X->Sym->basis[beta].rep_state,
                                        site1, spin1, site2, spin2, trans,
                                        &out_state, &hval) == TRUE &&
          emit_canonicalized_transition(X, beta, out_state, hval,
                                        callback, context) != 0) {
        return -1;
      }
    }
    return 0;
  }

  return -1;
}

struct CountEntriesContext {
  size_t count;
};

static int count_entry(unsigned long int out_index,
                       double complex coefficient,
                       void *context)
{
  struct CountEntriesContext *count = (struct CountEntriesContext *)context;
  (void)out_index;
  (void)coefficient;
  if (count->count == SIZE_MAX) return -1;
  count->count++;
  return 0;
}

struct FillEntriesContext {
  struct SymmetryMatvecPlan *plan;
  size_t next;
  size_t end;
};

static int fill_transposed_entry(unsigned long int out_index,
                                 double complex coefficient,
                                 void *context)
{
  struct FillEntriesContext *fill = (struct FillEntriesContext *)context;
  if (fill->next >= fill->end || out_index == 0UL ||
      out_index > fill->plan->dim) {
    return -1;
  }
  fill->plan->col_index[fill->next] = out_index;
  fill->plan->values[fill->next] = conj(coefficient);
  fill->next++;
  return 0;
}

static int parse_matvec_mode(void)
{
  const char *value = getenv("HPHI_SYMMETRY_MATVEC");
  if (value == NULL || strcmp(value, "plan") == 0) {
    return SYMMETRY_MATVEC_MODE_PLAN;
  }
  if (strcmp(value, "legacy") == 0) {
    return SYMMETRY_MATVEC_MODE_LEGACY;
  }
  fprintf(stdoutMPI,
          "Error: HPHI_SYMMETRY_MATVEC must be 'plan' or 'legacy', got '%s'.\n",
          value);
  return -1;
}

static int select_matvec_mode(void)
{
  int mode = SYMMETRY_MATVEC_MODE_PLAN;
  if (myrank == 0) mode = parse_matvec_mode();
  return BcastMPI_i(0, mode);
}

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

static int owner_of_global_index(unsigned long int dim,
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
  large_block_end = (base + 1UL) * remainder;
  if (zero_index < large_block_end) {
    return (int)(zero_index / (base + 1UL));
  }
  if (base == 0UL) return -1;
  return (int)(remainder + (zero_index - large_block_end) / base);
}

static int mpi_collectives_active(void)
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

static int agree_plan_error(int mpi_active, int local_error)
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

static int BuildSymmetryMatvecTopology(struct SymmetryMatvecPlan *plan)
{
  unsigned char *remote_bits = NULL;
  unsigned long long *recv_counts = NULL;
  unsigned long long *send_counts = NULL;
  unsigned long int first_local;
  unsigned long int last_local;
  size_t bitset_bytes = 0U;
  size_t peer_bytes = 0U;
  size_t schedule_count_bytes = 0U;
  size_t p;
  size_t byte_index;
  int mpi_active;
  int local_error = 0;
  int rank;

  if (plan == NULL || nproc < 1 || myrank < 0 || myrank >= nproc ||
      plan->local_offset > plan->dim ||
      plan->local_dim > plan->dim - plan->local_offset) {
    return -1;
  }
  StartTimer(1130);
  first_local = plan->local_offset + 1UL;
  last_local = plan->local_offset + plan->local_dim;
  if (plan->dim - plan->local_dim > (unsigned long int)SIZE_MAX) {
    local_error = 1;
  } else {
    plan->allgather_nonlocal_values_per_call =
        (size_t)(plan->dim - plan->local_dim);
    if (checked_size_mul(plan->allgather_nonlocal_values_per_call,
                         sizeof(double complex),
                         &plan->allgather_payload_bytes_per_call) != 0) {
      local_error = 1;
    }
  }

  for (p = 0U; p < plan->nnz; p++) {
    unsigned long int column = plan->col_index[p];
    if (column == 0UL || column > plan->dim) {
      local_error = 1;
      break;
    }
    if (plan->local_dim > 0UL &&
        column >= first_local && column <= last_local) {
      plan->local_column_nnz++;
    } else {
      plan->remote_column_nnz++;
    }
  }
  mpi_active = mpi_collectives_active();
  if (agree_plan_error(mpi_active, local_error) != 0) goto fail;
  if (mpi_active != FALSE) {
#ifdef MPI
    int comm_size = 0;
    if (MPI_Comm_size(MPI_COMM_WORLD, &comm_size) != MPI_SUCCESS ||
        comm_size != nproc) {
      local_error = 1;
    }
#endif
  }
  if (agree_plan_error(mpi_active, local_error) != 0) goto fail;

  if (plan->remote_column_nnz > 0U) {
    if (plan->dim / 8UL > (unsigned long int)SIZE_MAX) {
      local_error = 1;
    } else {
      bitset_bytes = (size_t)(plan->dim / 8UL);
    }
    if (local_error == 0 && plan->dim % 8UL != 0UL) {
      if (bitset_bytes == SIZE_MAX) {
        local_error = 1;
      } else {
        bitset_bytes++;
      }
    }
    if (local_error == 0) {
      remote_bits = (unsigned char *)calloc(bitset_bytes, sizeof(*remote_bits));
      if (remote_bits == NULL) local_error = 1;
    }
  }
  if (nproc > 1) {
    if ((size_t)nproc > SIZE_MAX / sizeof(*recv_counts)) {
      local_error = 1;
    } else {
      recv_counts = (unsigned long long *)calloc((size_t)nproc,
                                                 sizeof(*recv_counts));
      send_counts = (unsigned long long *)calloc((size_t)nproc,
                                                 sizeof(*send_counts));
      if (recv_counts == NULL || send_counts == NULL) local_error = 1;
      if (checked_size_mul(2U * sizeof(*recv_counts), (size_t)nproc,
                           &peer_bytes) != 0) {
        local_error = 1;
      }
    }
  }
  if (agree_plan_error(mpi_active, local_error) != 0) goto fail;

  for (p = 0U; p < plan->nnz; p++) {
    unsigned long int column = plan->col_index[p];
    if (plan->local_dim == 0UL ||
        column < first_local || column > last_local) {
      unsigned long int zero_index = column - 1UL;
      remote_bits[(size_t)(zero_index / 8UL)] |=
          (unsigned char)(1U << (unsigned int)(zero_index % 8UL));
    }
  }
  for (byte_index = 0U; byte_index < bitset_bytes; byte_index++) {
    unsigned char bits = remote_bits[byte_index];
    unsigned int bit;
    for (bit = 0U; bits != 0U && bit < 8U; bit++) {
      if ((bits & (unsigned char)(1U << bit)) != 0U) {
        unsigned long int global_index =
            (unsigned long int)(byte_index * 8U + (size_t)bit + 1U);
        int owner = owner_of_global_index(plan->dim, nproc, global_index);
        if (owner < 0 || owner == myrank) {
          local_error = 1;
          break;
        }
        plan->ghost_count++;
        recv_counts[owner]++;
      }
    }
    if (local_error != 0) break;
  }
  if (agree_plan_error(mpi_active, local_error) != 0) goto fail;

  for (rank = 0; rank < nproc; rank++) {
    if (recv_counts != NULL && recv_counts[rank] > 0ULL) {
      plan->incoming_peer_count++;
      if (recv_counts[rank] > (unsigned long long)SIZE_MAX) {
        local_error = 1;
      } else if ((size_t)recv_counts[rank] > plan->max_recv_from_peer) {
        plan->max_recv_from_peer = (size_t)recv_counts[rank];
      }
    }
  }
  if (agree_plan_error(mpi_active, local_error) != 0) goto fail;

#ifdef MPI
  if (mpi_active != FALSE && nproc > 1) {
    local_error =
        MPI_Alltoall(recv_counts, 1, MPI_UNSIGNED_LONG_LONG,
                     send_counts, 1, MPI_UNSIGNED_LONG_LONG,
                     MPI_COMM_WORLD) == MPI_SUCCESS ? 0 : 1;
  }
#endif
  if (agree_plan_error(mpi_active, local_error) != 0) goto fail;
  if (mpi_active != FALSE && nproc > 1) {
    for (rank = 0; rank < nproc; rank++) {
      if (send_counts[rank] > 0ULL) {
        if (send_counts[rank] > (unsigned long long)SIZE_MAX ||
            plan->send_value_count >
                SIZE_MAX - (size_t)send_counts[rank]) {
          local_error = 1;
          break;
        }
        plan->outgoing_peer_count++;
        plan->send_value_count += (size_t)send_counts[rank];
        if ((size_t)send_counts[rank] > plan->max_send_to_peer) {
          plan->max_send_to_peer = (size_t)send_counts[rank];
        }
      }
    }
  }
  if (agree_plan_error(mpi_active, local_error) != 0) goto fail;

  if (checked_size_mul(4U * sizeof(unsigned long long), (size_t)nproc,
                       &schedule_count_bytes) != 0 ||
      checked_size_add(bitset_bytes, peer_bytes,
                       &plan->topology_scratch_bytes) != 0 ||
      checked_size_add(plan->ghost_count, plan->send_value_count, &p) != 0 ||
      checked_size_mul(p, sizeof(unsigned long int),
                       &plan->halo_schedule_bytes_estimate) != 0 ||
      checked_size_add(plan->halo_schedule_bytes_estimate,
                       schedule_count_bytes,
                       &plan->halo_schedule_bytes_estimate) != 0 ||
      checked_size_mul(p, sizeof(double complex),
                       &plan->halo_runtime_buffer_bytes_estimate) != 0 ||
      checked_size_add((size_t)plan->local_dim, plan->ghost_count, &p) != 0) {
    local_error = 1;
  }
  if (agree_plan_error(mpi_active, local_error) != 0) goto fail;
  plan->column_slot_width = p <= (size_t)UINT32_MAX ? 32U : 64U;
  free(remote_bits);
  free(recv_counts);
  free(send_counts);
  StopTimer(1130);
  return 0;

fail:
  free(remote_bits);
  free(recv_counts);
  free(send_counts);
  StopTimer(1130);
  return -1;
}

void FreeSymmetryMatvecPlan(struct SymmetryMatvecPlan *plan)
{
  if (plan == NULL) return;
  free(plan->row_ptr);
  free(plan->col_index);
  free(plan->values);
  free(plan);
}

int BuildSymmetryMatvecPlan(struct BindStruct *X)
{
  unsigned long int local_row;
  size_t row_ptr_bytes = 0U;
  size_t col_bytes = 0U;
  size_t value_bytes = 0U;
  size_t plan_bytes = 0U;
  size_t min_row_nnz = SIZE_MAX;
  size_t max_row_nnz = 0U;
  size_t *row_counts = NULL;
  struct SymmetryMatvecPlan *plan = NULL;
  int count_error = 0;
  int fill_error = 0;
  int local_error = 0;
  int mpi_active;
  int mode;

  if (X == NULL || X->Sym == NULL || X->Sym->enabled != TRUE) return -1;
  FreeSymmetryMatvecPlan(X->Sym->matvec_plan);
  X->Sym->matvec_plan = NULL;

  mode = select_matvec_mode();
  if (mode < 0) return -1;
  X->Sym->matvec_mode = mode;
  if (mode == SYMMETRY_MATVEC_MODE_LEGACY) {
    fprintf(stdoutMPI,
            "Symmetry matvec: mode=legacy (replicated beta scan).\n");
    return 0;
  }

  mpi_active = mpi_collectives_active();
  plan = (struct SymmetryMatvecPlan *)calloc(1, sizeof(*plan));
  local_error = plan == NULL ? 1 : 0;
  if (agree_plan_error(mpi_active, local_error) != 0) goto fail;
  plan->dim = X->Sym->dim;
  plan->local_offset = X->Sym->local_offset;
  plan->local_dim = X->Sym->local_dim;

  if (plan->local_offset > plan->dim ||
      plan->local_dim > plan->dim - plan->local_offset ||
      plan->local_dim > (unsigned long int)(SIZE_MAX - 1U) ||
      (size_t)plan->local_dim + 1U > SIZE_MAX / sizeof(*plan->row_ptr)) {
    local_error = 1;
  } else {
    row_ptr_bytes = ((size_t)plan->local_dim + 1U) * sizeof(*plan->row_ptr);
    plan->row_ptr = (size_t *)calloc((size_t)plan->local_dim + 1U,
                                    sizeof(*plan->row_ptr));
    if (plan->row_ptr == NULL) local_error = 1;
  }
  if (local_error == 0 && plan->local_dim > 0UL) {
    if ((size_t)plan->local_dim > SIZE_MAX / sizeof(*row_counts)) {
      local_error = 1;
    }
  }
  if (local_error == 0 && plan->local_dim > 0UL) {
    row_counts = (size_t *)malloc((size_t)plan->local_dim *
                                 sizeof(*row_counts));
    if (row_counts == NULL) local_error = 1;
  }
  if (agree_plan_error(mpi_active, local_error) != 0) goto fail;

  StartTimer(1120);
#pragma omp parallel for default(none) schedule(static) reduction(|:count_error) \
  shared(X, plan, row_counts)
  for (local_row = 0UL; local_row < plan->local_dim; local_row++) {
    unsigned long int alpha = plan->local_offset + local_row + 1UL;
    struct CountEntriesContext count = {0U};
    if (SymmetryEnumerateColumn(X, alpha, count_entry, &count) != 0) {
      count_error = 1;
      row_counts[local_row] = 0U;
    } else {
      row_counts[local_row] = count.count;
    }
  }
  if (count_error == 0) {
    for (local_row = 0UL; local_row < plan->local_dim; local_row++) {
      size_t row_count = row_counts[local_row];
      if (plan->row_ptr[local_row] > SIZE_MAX - row_count) {
        count_error = 1;
        break;
      }
      plan->row_ptr[local_row + 1UL] = plan->row_ptr[local_row] + row_count;
      if (row_count < min_row_nnz) min_row_nnz = row_count;
      if (row_count > max_row_nnz) max_row_nnz = row_count;
    }
  }
  StopTimer(1120);
  if (agree_plan_error(mpi_active, count_error) != 0) goto fail;
  free(row_counts);
  row_counts = NULL;
  plan->nnz = plan->row_ptr[plan->local_dim];
  plan->row_nnz_max = max_row_nnz;

  StartTimer(1121);
  if (plan->nnz > SIZE_MAX / sizeof(*plan->col_index) ||
      plan->nnz > SIZE_MAX / sizeof(*plan->values)) {
    local_error = 1;
  } else {
    col_bytes = plan->nnz * sizeof(*plan->col_index);
    value_bytes = plan->nnz * sizeof(*plan->values);
  }
  if (local_error == 0 && plan->nnz > 0U) {
    plan->col_index = (unsigned long int *)malloc(col_bytes);
    plan->values = (double complex *)malloc(value_bytes);
    if (plan->col_index == NULL || plan->values == NULL) {
      local_error = 1;
    }
  }
  StopTimer(1121);
  if (agree_plan_error(mpi_active, local_error) != 0) goto fail;

  StartTimer(1122);
#pragma omp parallel for default(none) schedule(static) reduction(|:fill_error) \
  shared(X, plan)
  for (local_row = 0UL; local_row < plan->local_dim; local_row++) {
    unsigned long int alpha = plan->local_offset + local_row + 1UL;
    struct FillEntriesContext fill;
    fill.plan = plan;
    fill.next = plan->row_ptr[local_row];
    fill.end = plan->row_ptr[local_row + 1UL];
    if (SymmetryEnumerateColumn(X, alpha, fill_transposed_entry, &fill) != 0 ||
        fill.next != fill.end) {
      fill_error = 1;
    }
  }
  StopTimer(1122);
  if (agree_plan_error(mpi_active, fill_error) != 0) goto fail;

  if (BuildSymmetryMatvecTopology(plan) != 0) goto fail;

  if (row_ptr_bytes > SIZE_MAX - col_bytes ||
      row_ptr_bytes + col_bytes > SIZE_MAX - value_bytes) {
    local_error = 1;
  } else {
    plan_bytes = row_ptr_bytes + col_bytes + value_bytes;
  }
  if (agree_plan_error(mpi_active, local_error) != 0) goto fail;
  if (plan->local_dim == 0UL) min_row_nnz = 0U;
  plan->ready = TRUE;
  X->Sym->matvec_plan = plan;
  fprintf(stdoutMPI,
          "Symmetry matvec: mode=plan local_rows=%lu local_nnz=%zu "
          "row_nnz_min=%zu row_nnz_max=%zu row_nnz_mean=%.3f bytes=%zu "
          "local_column_nnz=%zu remote_column_nnz=%zu ghost_count=%zu "
          "incoming_peers=%zu outgoing_peers=%zu.\n",
          plan->local_dim, plan->nnz, min_row_nnz, max_row_nnz,
          plan->local_dim > 0UL ? (double)plan->nnz / (double)plan->local_dim : 0.0,
          plan_bytes, plan->local_column_nnz, plan->remote_column_nnz,
          plan->ghost_count, plan->incoming_peer_count,
          plan->outgoing_peer_count);
  return 0;

fail:
  free(row_counts);
  FreeSymmetryMatvecPlan(plan);
  fprintf(stdoutMPI, "Error: failed to build symmetry matvec plan.\n");
  return -1;
}

int ApplySymmetryMatvecPlan(const struct BindStruct *X,
                            double complex *tmp_v0,
                            const double complex *full_v1,
                            double complex *local_prdct)
{
  unsigned long int local_row;
  double complex prdct = 0.0;
  const struct SymmetryMatvecPlan *plan;
  if (X == NULL || X->Sym == NULL || tmp_v0 == NULL || full_v1 == NULL ||
      local_prdct == NULL) {
    return -1;
  }
  plan = X->Sym->matvec_plan;
  if (plan == NULL || plan->ready != TRUE || plan->dim != X->Sym->dim ||
      plan->local_offset != X->Sym->local_offset ||
      plan->local_dim != X->Sym->local_dim) {
    fprintf(stdoutMPI, "Error: symmetry matvec plan does not match the active sector.\n");
    return -1;
  }

  /* unsigned long canonical loops require OpenMP 3.0 or later. */
#pragma omp parallel for default(none) schedule(static) reduction(+:prdct) \
  shared(plan, tmp_v0, full_v1)
  for (local_row = 0UL; local_row < plan->local_dim; local_row++) {
    size_t p;
    double complex sum = 0.0;
    unsigned long int global_alpha = plan->local_offset + local_row + 1UL;
    for (p = plan->row_ptr[local_row]; p < plan->row_ptr[local_row + 1UL]; p++) {
      sum += plan->values[p] * full_v1[plan->col_index[p]];
    }
    tmp_v0[local_row + 1UL] += sum;
    prdct += conj(full_v1[global_alpha]) * sum;
  }
  *local_prdct = prdct;
  return 0;
}
