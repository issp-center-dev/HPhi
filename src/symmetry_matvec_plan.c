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
#include "symmetry_distribution.h"
#include "symmetry_matvec_plan.h"
#include "symmetry_mpi_exchange.h"
#include "symmetry_vector_halo.h"
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
  const struct SymmetryBasisVector *source;
  const struct SymmetryBasisVector *target;
  double norm_factor;
  double complex coefficient;
  struct SymmetryCanonicalResult result;
  if (SymmetryCanonicalizeState(X, to_state, &result) != 0) return -1;
  if (result.found != TRUE) return 0;
  source = SymmetryBasisReplicatedGlobalEntry(X->Sym, beta);
  target =
      SymmetryBasisReplicatedGlobalEntry(X->Sym, result.basis_index);
  if (source == NULL || target == NULL || source->norm == 0.0) return -1;
  norm_factor = target->norm / source->norm;
  coefficient = hval * result.phase * norm_factor;
  return callback(result.basis_index, coefficient, context);
}

int SymmetryEnumerateColumn(const struct BindStruct *X,
                            unsigned long int beta,
                            SymmetryEntryCallback callback,
                            void *context)
{
  const struct SymmetryBasisVector *source;
  unsigned int p;
  if (X == NULL || X->Sym == NULL || callback == NULL || beta == 0UL ||
      beta > X->Sym->dim) {
    return -1;
  }
  if (X->Sym->basis_layout != SYMMETRY_BASIS_REPLICATED) {
    fprintf(stdoutMPI,
            "Error: distributed symmetry basis column enumeration is "
            "staged for the B4 block plan.\n");
    return -1;
  }
  source = SymmetryBasisReplicatedGlobalEntry(X->Sym, beta);
  if (source == NULL) return -1;

  if (callback(beta, source->diagonal, context) != 0) return -1;

  if (X->Def.iCalcModel == Spin) {
    for (p = 0; p < X->Def.NExchangeCoupling; p++) {
      unsigned long int out_state;
      if (apply_exchange_halfspin(source->rep_state,
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
      if (apply_spinless_hopping_hermite(source->rep_state,
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
      if (apply_hubbard_hopping_hermite(source->rep_state,
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
  struct SymmetryMatvecBlock *block;
  unsigned long int dim;
  size_t next;
  size_t end;
};

static int fill_transposed_entry(unsigned long int out_index,
                                 double complex coefficient,
                                 void *context)
{
  struct FillEntriesContext *fill = (struct FillEntriesContext *)context;
  if (fill->next >= fill->end || out_index == 0UL ||
      out_index > fill->dim) {
    return -1;
  }
  fill->block->global_columns[fill->next] = out_index;
  fill->block->values[fill->next] = conj(coefficient);
  fill->next++;
  return 0;
}

static void sync_single_block_aliases(struct SymmetryMatvecPlan *plan)
{
  struct SymmetryMatvecBlock *block;
  if (plan == NULL) return;
  plan->row_ptr = NULL;
  plan->col_index = NULL;
  plan->column_slot32 = NULL;
  plan->column_slot64 = NULL;
  plan->values = NULL;
  if (plan->block_count != 1U || plan->blocks == NULL) return;
  block = &plan->blocks[0];
  plan->row_ptr = block->row_ptr;
  plan->col_index = block->global_columns;
  plan->column_slot32 = block->column_slot32;
  plan->column_slot64 = block->column_slot64;
  plan->values = block->values;
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

static int parse_vector_exchange_mode(int matvec_mode)
{
  const char *value = getenv("HPHI_SYMMETRY_VECTOR_EXCHANGE");
  if (value == NULL) {
    /* Keep the one-switch legacy rollback usable while plan uses halo. */
    return matvec_mode == SYMMETRY_MATVEC_MODE_LEGACY
               ? SYMMETRY_VECTOR_EXCHANGE_ALLGATHER
               : SYMMETRY_VECTOR_EXCHANGE_HALO;
  }
  if (strcmp(value, "allgather") == 0) {
    return SYMMETRY_VECTOR_EXCHANGE_ALLGATHER;
  }
  if (strcmp(value, "halo") == 0) {
    return SYMMETRY_VECTOR_EXCHANGE_HALO;
  }
  fprintf(stdoutMPI,
          "Error: HPHI_SYMMETRY_VECTOR_EXCHANGE must be "
          "'allgather' or 'halo', got '%s'.\n",
          value);
  return -1;
}

static int select_vector_exchange_mode(int matvec_mode)
{
  int mode = SYMMETRY_VECTOR_EXCHANGE_ALLGATHER;
  if (myrank == 0) mode = parse_vector_exchange_mode(matvec_mode);
  return BcastMPI_i(0, mode);
}

static int parse_halo_reference_mode(void)
{
  const char *value = getenv("HPHI_SYMMETRY_HALO_REFERENCE");
  if (value == NULL || strcmp(value, "0") == 0 ||
      strcmp(value, "off") == 0) {
    return FALSE;
  }
  if (strcmp(value, "1") == 0 || strcmp(value, "on") == 0) {
    return TRUE;
  }
  fprintf(stdoutMPI,
          "Error: HPHI_SYMMETRY_HALO_REFERENCE must be "
          "'0'/'off' or '1'/'on', got '%s'.\n",
          value);
  return -1;
}

static int select_halo_reference_mode(void)
{
  int enabled = FALSE;
  if (myrank == 0) enabled = parse_halo_reference_mode();
  return BcastMPI_i(0, enabled);
}

static int configure_full_input_vector(struct BindStruct *X,
                                       int required,
                                       int mpi_active)
{
#ifdef MPI
  int local_error = 0;
  int rank;
  free(X->Sym->mpi_recvcounts);
  free(X->Sym->mpi_displs);
  free(X->Sym->mpi_full_v1);
  X->Sym->mpi_recvcounts = NULL;
  X->Sym->mpi_displs = NULL;
  X->Sym->mpi_full_v1 = NULL;
  if (required == FALSE || nproc <= 1) return 0;
  if (X->Sym->dim > (unsigned long int)INT_MAX ||
      X->Sym->dim > (unsigned long int)(SIZE_MAX / sizeof(double complex) -
                                        1U)) {
    local_error = 1;
  } else {
    X->Sym->mpi_recvcounts =
        (int *)calloc((size_t)nproc, sizeof(*X->Sym->mpi_recvcounts));
    X->Sym->mpi_displs =
        (int *)calloc((size_t)nproc, sizeof(*X->Sym->mpi_displs));
    X->Sym->mpi_full_v1 =
        (double complex *)calloc((size_t)X->Sym->dim + 1U,
                                 sizeof(*X->Sym->mpi_full_v1));
    if (X->Sym->mpi_recvcounts == NULL || X->Sym->mpi_displs == NULL ||
        X->Sym->mpi_full_v1 == NULL) {
      local_error = 1;
    }
  }
  if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) {
    free(X->Sym->mpi_recvcounts);
    free(X->Sym->mpi_displs);
    free(X->Sym->mpi_full_v1);
    X->Sym->mpi_recvcounts = NULL;
    X->Sym->mpi_displs = NULL;
    X->Sym->mpi_full_v1 = NULL;
    fprintf(stdoutMPI,
            "Error: failed to allocate the symmetry full input vector.\n");
    return -1;
  }
  for (rank = 0; rank < nproc; rank++) {
    unsigned long int offset;
    unsigned long int count;
    if (SymmetryBlockRange(X->Sym->dim, rank, nproc,
                           &offset, &count) != 0) {
      return -1;
    }
    X->Sym->mpi_recvcounts[rank] = (int)count;
    X->Sym->mpi_displs[rank] = (int)offset;
  }
#else
  (void)X;
  (void)required;
  (void)mpi_active;
#endif
  return 0;
}

static int find_ghost_position(const struct SymmetryVectorHaloPlan *halo,
                               unsigned long int global_index,
                               size_t *position)
{
  size_t left = 0U;
  size_t right = halo->ghost_count;
  while (left < right) {
    size_t middle = left + (right - left) / 2U;
    unsigned long int candidate = halo->ghost_global_index[middle];
    if (candidate < global_index) {
      left = middle + 1U;
    } else {
      right = middle;
    }
  }
  if (left >= halo->ghost_count ||
      halo->ghost_global_index[left] != global_index) {
    return -1;
  }
  *position = left;
  return 0;
}

static int build_symmetry_matvec_plan_halo(
    struct SymmetryMatvecPlan *plan,
    int nrank,
    int rank)
{
  struct SymmetryGlobalColumnSpan *spans;
  size_t block_index;
  int status;
  if (plan == NULL || plan->block_count == 0U || plan->blocks == NULL ||
      plan->block_count > SIZE_MAX / sizeof(*spans)) {
    return -1;
  }
  spans = (struct SymmetryGlobalColumnSpan *)calloc(
      plan->block_count, sizeof(*spans));
  if (spans == NULL) return -1;
  for (block_index = 0U; block_index < plan->block_count; block_index++) {
    const struct SymmetryMatvecBlock *block = &plan->blocks[block_index];
    if (block->nnz > 0U && block->global_columns == NULL) {
      free(spans);
      return -1;
    }
    spans[block_index].columns = block->global_columns;
    spans[block_index].count = block->nnz;
  }
  status = BuildSymmetryVectorHaloPlan(
      &plan->halo, plan->dim, plan->local_offset, plan->local_dim,
      spans, plan->block_count, nrank, rank,
      &plan->local_column_nnz, &plan->remote_column_nnz);
  free(spans);
  return status;
}

static int remap_symmetry_matvec_block_columns(
    const struct SymmetryMatvecPlan *plan,
    const struct SymmetryMatvecBlockView *view,
    size_t slot_count,
    enum SymmetryColumnWidth column_width,
    uint32_t *column_slot32,
    uint64_t *column_slot64)
{
  size_t column;
  int remap_error = 0;
#pragma omp parallel for default(none) schedule(static) reduction(|:remap_error) \
  shared(plan, view, slot_count, column_width, column_slot32, column_slot64)
  for (column = 0U; column < view->nnz; column++) {
    unsigned long int global_index = view->global_columns[column];
    size_t slot = 0U;
    if (global_index == 0UL || global_index > plan->dim) {
      remap_error = 1;
      continue;
    }
    if (global_index > plan->local_offset &&
        global_index <= plan->local_offset + plan->local_dim) {
      slot = (size_t)(global_index - plan->local_offset - 1UL);
    } else {
      size_t ghost_position = 0U;
      if (find_ghost_position(&plan->halo, global_index,
                              &ghost_position) != 0 ||
          (size_t)plan->local_dim > SIZE_MAX - ghost_position) {
        remap_error = 1;
        continue;
      }
      slot = (size_t)plan->local_dim + ghost_position;
    }
    if (slot >= slot_count) {
      remap_error = 1;
      continue;
    }
    if (column_width == SYMMETRY_COLUMN_U32) {
      if (slot > (size_t)UINT32_MAX) {
        remap_error = 1;
        continue;
      }
      column_slot32[column] = (uint32_t)slot;
    } else {
      column_slot64[column] = (uint64_t)slot;
    }
  }
  return remap_error == 0 ? 0 : -1;
}

int RemapSymmetryMatvecPlanColumns(struct SymmetryMatvecPlan *plan)
{
  size_t slot_count;
  size_t block_index;
  size_t total_nnz = 0U;
  unsigned long int covered_rows = 0UL;
  enum SymmetryColumnWidth column_width;
  int remap_error = 0;
  if (plan == NULL || plan->columns_remapped == TRUE ||
      plan->local_offset > plan->dim ||
      plan->local_dim > plan->dim - plan->local_offset ||
      plan->block_count == 0U || plan->blocks == NULL ||
      (size_t)plan->local_dim > SIZE_MAX - plan->halo.ghost_count) {
    return -1;
  }
  for (block_index = 0U; block_index < plan->block_count; block_index++) {
    const struct SymmetryMatvecBlock *block = &plan->blocks[block_index];
    if (block->local_row_begin != covered_rows ||
        block->local_row_count > plan->local_dim - covered_rows ||
        (block->nnz > 0U && block->global_columns == NULL) ||
        block->column_slot32 != NULL || block->column_slot64 != NULL ||
        total_nnz > SIZE_MAX - block->nnz) {
      return -1;
    }
    covered_rows += block->local_row_count;
    total_nnz += block->nnz;
  }
  if (covered_rows != plan->local_dim || total_nnz != plan->nnz) return -1;
  slot_count = (size_t)plan->local_dim + plan->halo.ghost_count;
  column_width =
      slot_count <= (size_t)UINT32_MAX ? SYMMETRY_COLUMN_U32
                                       : SYMMETRY_COLUMN_U64;
  StartTimer(1132);
  for (block_index = 0U; block_index < plan->block_count; block_index++) {
    struct SymmetryMatvecBlock *block = &plan->blocks[block_index];
    struct SymmetryMatvecBlockView view;
    memset(&view, 0, sizeof(view));
    view.local_row_begin = block->local_row_begin;
    view.local_row_count = block->local_row_count;
    view.nnz = block->nnz;
    view.row_ptr = block->row_ptr;
    view.global_columns = block->global_columns;
    view.values = block->values;
    if (block->nnz == 0U) continue;
    if (column_width == SYMMETRY_COLUMN_U32) {
      if (block->nnz > SIZE_MAX / sizeof(*block->column_slot32)) {
        remap_error = 1;
      } else {
        block->column_slot32 =
            (uint32_t *)malloc(block->nnz * sizeof(*block->column_slot32));
        if (block->column_slot32 == NULL) remap_error = 1;
      }
    } else {
      if (block->nnz > SIZE_MAX / sizeof(*block->column_slot64)) {
        remap_error = 1;
      } else {
        block->column_slot64 =
            (uint64_t *)malloc(block->nnz * sizeof(*block->column_slot64));
        if (block->column_slot64 == NULL) remap_error = 1;
      }
    }
    if (remap_error == 0 &&
        remap_symmetry_matvec_block_columns(
            plan, &view, slot_count, column_width,
            block->column_slot32, block->column_slot64) != 0) {
      remap_error = 1;
    }
    if (remap_error != 0) break;
  }
  StopTimer(1132);
  if (remap_error != 0) {
    for (block_index = 0U; block_index < plan->block_count; block_index++) {
      free(plan->blocks[block_index].column_slot32);
      free(plan->blocks[block_index].column_slot64);
      plan->blocks[block_index].column_slot32 = NULL;
      plan->blocks[block_index].column_slot64 = NULL;
    }
    return -1;
  }
  for (block_index = 0U; block_index < plan->block_count; block_index++) {
    free(plan->blocks[block_index].global_columns);
    plan->blocks[block_index].global_columns = NULL;
  }
  plan->column_slot_width = column_width;
  plan->columns_remapped = TRUE;
  sync_single_block_aliases(plan);
  return 0;
}

size_t SymmetryMatvecPlanBlockCount(
    const struct SymmetryMatvecPlan *plan)
{
  if (plan == NULL || plan->ready != TRUE) return 0U;
  return plan->block_count;
}

int SymmetryMatvecPlanGetBlockView(
    const struct SymmetryMatvecPlan *plan,
    size_t block_index,
    struct SymmetryMatvecBlockView *view)
{
  const struct SymmetryMatvecBlock *block;
  if (view == NULL) return -1;
  memset(view, 0, sizeof(*view));
  if (plan == NULL || plan->ready != TRUE || plan->blocks == NULL ||
      block_index >= plan->block_count) {
    return -1;
  }
  block = &plan->blocks[block_index];
  if (block->local_row_begin > plan->local_dim ||
      block->local_row_count >
          plan->local_dim - block->local_row_begin ||
      block->row_ptr == NULL || block->row_ptr[0] != 0U ||
      block->row_ptr[block->local_row_count] != block->nnz ||
      (block->nnz > 0U && block->values == NULL)) {
    return -1;
  }
  if (plan->columns_remapped == TRUE) {
    if (block->global_columns != NULL ||
        (plan->column_slot_width == SYMMETRY_COLUMN_U32 &&
         (block->column_slot64 != NULL ||
          (block->nnz > 0U && block->column_slot32 == NULL))) ||
        (plan->column_slot_width == SYMMETRY_COLUMN_U64 &&
         (block->column_slot32 != NULL ||
          (block->nnz > 0U && block->column_slot64 == NULL))) ||
        (plan->column_slot_width != SYMMETRY_COLUMN_U32 &&
         plan->column_slot_width != SYMMETRY_COLUMN_U64)) {
      return -1;
    }
  } else if ((block->nnz > 0U && block->global_columns == NULL) ||
             block->column_slot32 != NULL ||
             block->column_slot64 != NULL) {
    return -1;
  }
  view->local_row_begin = block->local_row_begin;
  view->local_row_count = block->local_row_count;
  view->nnz = block->nnz;
  view->row_ptr = block->row_ptr;
  view->global_columns = block->global_columns;
  view->column_slot32 = block->column_slot32;
  view->column_slot64 = block->column_slot64;
  view->values = block->values;
  return 0;
}

void FreeSymmetryMatvecPlan(struct SymmetryMatvecPlan *plan)
{
  size_t block_index;
  if (plan == NULL) return;
  FreeSymmetryVectorHaloPlan(&plan->halo);
  if (plan->blocks != NULL) {
    for (block_index = 0U; block_index < plan->block_count; block_index++) {
      free(plan->blocks[block_index].row_ptr);
      free(plan->blocks[block_index].global_columns);
      free(plan->blocks[block_index].column_slot32);
      free(plan->blocks[block_index].column_slot64);
      free(plan->blocks[block_index].values);
    }
  }
  free(plan->blocks);
  free(plan);
}

int BuildSymmetryMatvecPlan(struct BindStruct *X)
{
  unsigned long int local_row;
  size_t row_ptr_bytes = 0U;
  size_t col_bytes = 0U;
  size_t ghost_index_bytes = 0U;
  size_t value_bytes = 0U;
  size_t plan_bytes = 0U;
  size_t min_row_nnz = SIZE_MAX;
  size_t max_row_nnz = 0U;
  size_t *row_counts = NULL;
  struct SymmetryMatvecBlock *block = NULL;
  struct SymmetryMatvecPlan *plan = NULL;
  int count_error = 0;
  int fill_error = 0;
  int local_error = 0;
  int mpi_active;
  int mode;
  int vector_exchange_mode;
  int halo_reference_mode;

  if (X == NULL || X->Sym == NULL || X->Sym->enabled != TRUE) return -1;
  if (X->Sym->basis_layout != SYMMETRY_BASIS_REPLICATED) {
    fprintf(stdoutMPI,
            "Error: distributed symmetry basis plan construction is "
            "staged for B4 and cannot allocate legacy/allgather state.\n");
    return -1;
  }
  FreeSymmetryMatvecPlan(X->Sym->matvec_plan);
  X->Sym->matvec_plan = NULL;

  mode = select_matvec_mode();
  if (mode < 0) return -1;
  vector_exchange_mode = select_vector_exchange_mode(mode);
  if (vector_exchange_mode < 0) return -1;
  halo_reference_mode = select_halo_reference_mode();
  if (halo_reference_mode < 0) return -1;
  if (mode == SYMMETRY_MATVEC_MODE_LEGACY &&
      vector_exchange_mode == SYMMETRY_VECTOR_EXCHANGE_HALO) {
    fprintf(stdoutMPI,
            "Error: symmetry legacy matvec requires "
            "HPHI_SYMMETRY_VECTOR_EXCHANGE=allgather.\n");
    return -1;
  }
  if (halo_reference_mode == TRUE &&
      (mode != SYMMETRY_MATVEC_MODE_PLAN ||
       vector_exchange_mode != SYMMETRY_VECTOR_EXCHANGE_ALLGATHER)) {
    fprintf(stdoutMPI,
            "Error: HPHI_SYMMETRY_HALO_REFERENCE requires "
            "plan/allgather mode.\n");
    return -1;
  }
  X->Sym->matvec_mode = mode;
  X->Sym->vector_exchange_mode = vector_exchange_mode;
  mpi_active = SymmetryMpiCollectivesActive();
  if (configure_full_input_vector(
          X, vector_exchange_mode == SYMMETRY_VECTOR_EXCHANGE_ALLGATHER,
          mpi_active) != 0) {
    return -1;
  }
  if (mode == SYMMETRY_MATVEC_MODE_LEGACY) {
    fprintf(stdoutMPI,
            "Symmetry matvec: mode=legacy vector_exchange=allgather "
            "(replicated beta scan).\n");
    return 0;
  }

  plan = (struct SymmetryMatvecPlan *)calloc(1, sizeof(*plan));
  local_error = plan == NULL ? 1 : 0;
  if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;
  plan->dim = X->Sym->dim;
  plan->local_offset = X->Sym->local_offset;
  plan->local_dim = X->Sym->local_dim;
  plan->block_count = 1U;
  plan->blocks =
      (struct SymmetryMatvecBlock *)calloc(1U, sizeof(*plan->blocks));
  if (plan->blocks == NULL) local_error = 1;
  if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;
  block = &plan->blocks[0];
  block->local_row_begin = 0UL;
  block->local_row_count = plan->local_dim;

  if (plan->local_offset > plan->dim ||
      plan->local_dim > plan->dim - plan->local_offset ||
      plan->local_dim > (unsigned long int)(SIZE_MAX - 1U) ||
      (size_t)plan->local_dim + 1U > SIZE_MAX / sizeof(*block->row_ptr)) {
    local_error = 1;
  } else {
    row_ptr_bytes =
        ((size_t)plan->local_dim + 1U) * sizeof(*block->row_ptr);
    block->row_ptr =
        (size_t *)calloc((size_t)plan->local_dim + 1U,
                         sizeof(*block->row_ptr));
    if (block->row_ptr == NULL) local_error = 1;
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
  if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;

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
      if (block->row_ptr[local_row] > SIZE_MAX - row_count) {
        count_error = 1;
        break;
      }
      block->row_ptr[local_row + 1UL] =
          block->row_ptr[local_row] + row_count;
      if (row_count < min_row_nnz) min_row_nnz = row_count;
      if (row_count > max_row_nnz) max_row_nnz = row_count;
    }
  }
  StopTimer(1120);
  if (SymmetryMpiAgreeError(mpi_active, count_error) != 0) goto fail;
  free(row_counts);
  row_counts = NULL;
  block->nnz = block->row_ptr[block->local_row_count];
  plan->nnz = block->nnz;
  plan->row_nnz_max = max_row_nnz;

  StartTimer(1121);
  if (block->nnz > SIZE_MAX / sizeof(*block->global_columns) ||
      block->nnz > SIZE_MAX / sizeof(*block->values)) {
    local_error = 1;
  } else {
    col_bytes = block->nnz * sizeof(*block->global_columns);
    value_bytes = block->nnz * sizeof(*block->values);
  }
  if (local_error == 0 && block->nnz > 0U) {
    block->global_columns = (unsigned long int *)malloc(col_bytes);
    block->values = (double complex *)malloc(value_bytes);
    if (block->global_columns == NULL || block->values == NULL) {
      local_error = 1;
    }
  }
  StopTimer(1121);
  if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;

  StartTimer(1122);
#pragma omp parallel for default(none) schedule(static) reduction(|:fill_error) \
  shared(X, plan, block)
  for (local_row = 0UL; local_row < plan->local_dim; local_row++) {
    unsigned long int alpha = plan->local_offset + local_row + 1UL;
    struct FillEntriesContext fill;
    fill.block = block;
    fill.dim = plan->dim;
    fill.next = block->row_ptr[local_row];
    fill.end = block->row_ptr[local_row + 1UL];
    if (SymmetryEnumerateColumn(X, alpha, fill_transposed_entry, &fill) != 0 ||
        fill.next != fill.end) {
      fill_error = 1;
    }
  }
  StopTimer(1122);
  if (SymmetryMpiAgreeError(mpi_active, fill_error) != 0) goto fail;

  if (plan->dim - plan->local_dim > (unsigned long int)SIZE_MAX ||
      (size_t)(plan->dim - plan->local_dim) >
          SIZE_MAX / sizeof(double complex)) {
    local_error = 1;
  } else {
    plan->allgather_nonlocal_values_per_call =
        (size_t)(plan->dim - plan->local_dim);
    plan->allgather_payload_bytes_per_call =
        plan->allgather_nonlocal_values_per_call * sizeof(double complex);
  }
  if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;
  if (build_symmetry_matvec_plan_halo(plan, nproc, myrank) != 0) {
    goto fail;
  }
  plan->halo.reference_enabled = halo_reference_mode;
  if ((size_t)plan->local_dim >
      SIZE_MAX - plan->halo.ghost_count) {
    local_error = 1;
  } else {
    size_t slot_count = (size_t)plan->local_dim + plan->halo.ghost_count;
    plan->column_slot_width =
        slot_count <= (size_t)UINT32_MAX ? SYMMETRY_COLUMN_U32
                                         : SYMMETRY_COLUMN_U64;
  }
  if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;
  if (vector_exchange_mode == SYMMETRY_VECTOR_EXCHANGE_HALO) {
    if (plan->halo.ready != TRUE ||
        RemapSymmetryMatvecPlanColumns(plan) != 0) {
      local_error = 1;
    }
    if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;
    if (plan->column_slot_width == SYMMETRY_COLUMN_U32) {
      col_bytes = plan->nnz * sizeof(*block->column_slot32);
    } else {
      col_bytes = plan->nnz * sizeof(*block->column_slot64);
    }
    if (plan->halo.ghost_count >
        SIZE_MAX / sizeof(*plan->halo.ghost_global_index)) {
      local_error = 1;
    } else {
      ghost_index_bytes =
          plan->halo.ghost_count * sizeof(*plan->halo.ghost_global_index);
      if (ghost_index_bytes > plan->halo.schedule_bytes) local_error = 1;
    }
    if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;
    free(plan->halo.ghost_global_index);
    plan->halo.ghost_global_index = NULL;
    plan->halo.schedule_bytes -= ghost_index_bytes;
  }

  if (row_ptr_bytes > SIZE_MAX - col_bytes ||
      row_ptr_bytes + col_bytes > SIZE_MAX - value_bytes) {
    local_error = 1;
  } else {
    plan_bytes = row_ptr_bytes + col_bytes + value_bytes;
  }
  if (SymmetryMpiAgreeError(mpi_active, local_error) != 0) goto fail;
  plan->column_storage_bytes = col_bytes;
  plan->matrix_storage_bytes = plan_bytes;
  if (plan->local_dim == 0UL) min_row_nnz = 0U;
  sync_single_block_aliases(plan);
  plan->ready = TRUE;
  X->Sym->matvec_plan = plan;
  fprintf(stdoutMPI,
          "Symmetry matvec: mode=plan vector_exchange=%s "
          "local_rows=%lu local_nnz=%zu "
          "row_nnz_min=%zu row_nnz_max=%zu row_nnz_mean=%.3f bytes=%zu "
          "local_column_nnz=%zu remote_column_nnz=%zu ghost_count=%zu "
          "incoming_peers=%zu outgoing_peers=%zu halo_schedule=%s "
          "columns=%s column_slot_width=%u halo_reference=%s.\n",
          vector_exchange_mode == SYMMETRY_VECTOR_EXCHANGE_HALO
              ? "halo"
              : "allgather",
          plan->local_dim, plan->nnz, min_row_nnz, max_row_nnz,
          plan->local_dim > 0UL ? (double)plan->nnz / (double)plan->local_dim : 0.0,
          plan_bytes, plan->local_column_nnz, plan->remote_column_nnz,
          plan->halo.ghost_count, plan->halo.incoming_peer_count,
          plan->halo.outgoing_peer_count,
          plan->halo.ready == TRUE ? "ready" : "request-layout-only",
          plan->columns_remapped == TRUE ? "local/ghost-slots" : "global",
          (unsigned int)plan->column_slot_width,
          plan->halo.reference_enabled == TRUE ? "on" : "off");
  return 0;

fail:
  free(row_counts);
  FreeSymmetryMatvecPlan(plan);
  fprintf(stdoutMPI, "Error: failed to build symmetry matvec plan.\n");
  return -1;
}

static int apply_symmetry_matvec_plan_block(
    const struct SymmetryMatvecPlan *plan,
    const struct SymmetryMatvecBlockView *view,
    double complex *tmp_v0,
    const double complex *full_v1,
    double complex *local_prdct)
{
  unsigned long int block_row;
  double complex prdct = 0.0;
  /* unsigned long canonical loops require OpenMP 3.0 or later. */
#pragma omp parallel for default(none) schedule(static) reduction(+:prdct) \
  shared(plan, view, tmp_v0, full_v1)
  for (block_row = 0UL; block_row < view->local_row_count; block_row++) {
    size_t p;
    double complex sum = 0.0;
    unsigned long int local_row = view->local_row_begin + block_row;
    unsigned long int global_alpha =
        plan->local_offset + local_row + 1UL;
    for (p = view->row_ptr[block_row];
         p < view->row_ptr[block_row + 1UL]; p++) {
      sum += view->values[p] * full_v1[view->global_columns[p]];
    }
    tmp_v0[local_row + 1UL] += sum;
    prdct += conj(full_v1[global_alpha]) * sum;
  }
  *local_prdct = prdct;
  return 0;
}

int ApplySymmetryMatvecPlan(const struct BindStruct *X,
                            double complex *tmp_v0,
                            const double complex *full_v1,
                            double complex *local_prdct)
{
  const struct SymmetryMatvecPlan *plan;
  unsigned long int covered_rows = 0UL;
  double complex prdct = 0.0;
  size_t block_index;
  if (X == NULL || X->Sym == NULL || tmp_v0 == NULL || full_v1 == NULL ||
      local_prdct == NULL) {
    return -1;
  }
  plan = X->Sym->matvec_plan;
  if (plan == NULL || plan->ready != TRUE || plan->dim != X->Sym->dim ||
      plan->local_offset != X->Sym->local_offset ||
      plan->local_dim != X->Sym->local_dim ||
      plan->columns_remapped == TRUE || plan->blocks == NULL ||
      plan->block_count == 0U) {
    fprintf(stdoutMPI,
            "Error: symmetry matvec plan does not match the active sector.\n");
    return -1;
  }
  for (block_index = 0U; block_index < plan->block_count; block_index++) {
    struct SymmetryMatvecBlockView view;
    double complex block_prdct = 0.0;
    if (SymmetryMatvecPlanGetBlockView(plan, block_index, &view) != 0 ||
        view.local_row_begin != covered_rows ||
        (view.nnz > 0U && view.global_columns == NULL) ||
        apply_symmetry_matvec_plan_block(
            plan, &view, tmp_v0, full_v1, &block_prdct) != 0) {
      return -1;
    }
    covered_rows += view.local_row_count;
    prdct += block_prdct;
  }
  if (covered_rows != plan->local_dim) return -1;
  *local_prdct = prdct;
  return 0;
}

static int apply_symmetry_matvec_plan_halo_u32(
    const struct SymmetryMatvecPlan *plan,
    const struct SymmetryMatvecBlockView *view,
    double complex *tmp_v0,
    const double complex *local_v1,
    size_t slot_count,
    double complex *local_prdct)
{
  unsigned long int block_row;
  double complex prdct = 0.0;
  int apply_error = 0;
#pragma omp parallel for default(none) schedule(static) \
  reduction(+:prdct) reduction(|:apply_error) \
  shared(plan, view, tmp_v0, local_v1, slot_count)
  for (block_row = 0UL; block_row < view->local_row_count; block_row++) {
    size_t p;
    double complex sum = 0.0;
    unsigned long int local_row = view->local_row_begin + block_row;
    for (p = view->row_ptr[block_row];
         p < view->row_ptr[block_row + 1UL]; p++) {
      size_t slot = (size_t)view->column_slot32[p];
      double complex input_value;
      if (slot >= slot_count) {
        apply_error = 1;
        continue;
      }
      input_value =
          slot < (size_t)plan->local_dim
              ? local_v1[slot + 1U]
              : plan->halo.ghost_values[slot - (size_t)plan->local_dim];
      sum += view->values[p] * input_value;
    }
    tmp_v0[local_row + 1UL] += sum;
    prdct += conj(local_v1[local_row + 1UL]) * sum;
  }
  if (apply_error != 0) return -1;
  *local_prdct = prdct;
  return 0;
}

static int apply_symmetry_matvec_plan_halo_u64(
    const struct SymmetryMatvecPlan *plan,
    const struct SymmetryMatvecBlockView *view,
    double complex *tmp_v0,
    const double complex *local_v1,
    size_t slot_count,
    double complex *local_prdct)
{
  unsigned long int block_row;
  double complex prdct = 0.0;
  int apply_error = 0;
#pragma omp parallel for default(none) schedule(static) \
  reduction(+:prdct) reduction(|:apply_error) \
  shared(plan, view, tmp_v0, local_v1, slot_count)
  for (block_row = 0UL; block_row < view->local_row_count; block_row++) {
    size_t p;
    double complex sum = 0.0;
    unsigned long int local_row = view->local_row_begin + block_row;
    for (p = view->row_ptr[block_row];
         p < view->row_ptr[block_row + 1UL]; p++) {
      uint64_t raw_slot = view->column_slot64[p];
      size_t slot;
      double complex input_value;
      if (raw_slot >= (uint64_t)slot_count) {
        apply_error = 1;
        continue;
      }
      slot = (size_t)raw_slot;
      input_value =
          slot < (size_t)plan->local_dim
              ? local_v1[slot + 1U]
              : plan->halo.ghost_values[slot - (size_t)plan->local_dim];
      sum += view->values[p] * input_value;
    }
    tmp_v0[local_row + 1UL] += sum;
    prdct += conj(local_v1[local_row + 1UL]) * sum;
  }
  if (apply_error != 0) return -1;
  *local_prdct = prdct;
  return 0;
}

int ApplySymmetryMatvecPlanHalo(const struct BindStruct *X,
                                double complex *tmp_v0,
                                const double complex *local_v1,
                                double complex *local_prdct)
{
  const struct SymmetryMatvecPlan *plan;
  size_t slot_count;
  size_t block_index;
  unsigned long int covered_rows = 0UL;
  double complex prdct = 0.0;
  if (X == NULL || X->Sym == NULL || tmp_v0 == NULL || local_v1 == NULL ||
      local_prdct == NULL) {
    return -1;
  }
  plan = X->Sym->matvec_plan;
  if (plan == NULL || plan->ready != TRUE ||
      plan->columns_remapped != TRUE || plan->halo.ready != TRUE ||
      plan->dim != X->Sym->dim ||
      plan->local_offset != X->Sym->local_offset ||
      plan->local_dim != X->Sym->local_dim ||
      plan->blocks == NULL || plan->block_count == 0U ||
      (plan->column_slot_width != SYMMETRY_COLUMN_U32 &&
       plan->column_slot_width != SYMMETRY_COLUMN_U64) ||
      (size_t)plan->local_dim > SIZE_MAX - plan->halo.ghost_count) {
    fprintf(stdoutMPI,
            "Error: symmetry halo matvec plan does not match "
            "the active sector.\n");
    return -1;
  }
  slot_count = (size_t)plan->local_dim + plan->halo.ghost_count;
  for (block_index = 0U; block_index < plan->block_count; block_index++) {
    struct SymmetryMatvecBlockView view;
    double complex block_prdct = 0.0;
    int status;
    if (SymmetryMatvecPlanGetBlockView(plan, block_index, &view) != 0 ||
        view.local_row_begin != covered_rows) {
      return -1;
    }
    if (plan->column_slot_width == SYMMETRY_COLUMN_U32) {
      status = apply_symmetry_matvec_plan_halo_u32(
          plan, &view, tmp_v0, local_v1, slot_count, &block_prdct);
    } else {
      status = apply_symmetry_matvec_plan_halo_u64(
          plan, &view, tmp_v0, local_v1, slot_count, &block_prdct);
    }
    if (status != 0) return -1;
    covered_rows += view.local_row_count;
    prdct += block_prdct;
  }
  if (covered_rows != plan->local_dim) return -1;
  *local_prdct = prdct;
  return 0;
}
