#include <limits.h>
#include <stdint.h>
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
  size_t row_ptr_bytes;
  size_t col_bytes;
  size_t value_bytes;
  size_t plan_bytes;
  size_t min_row_nnz = SIZE_MAX;
  size_t max_row_nnz = 0U;
  struct SymmetryMatvecPlan *plan;
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

  plan = (struct SymmetryMatvecPlan *)calloc(1, sizeof(*plan));
  if (plan == NULL) return -1;
  plan->dim = X->Sym->dim;
  plan->local_offset = X->Sym->local_offset;
  plan->local_dim = X->Sym->local_dim;

  if (plan->local_dim > (unsigned long int)(SIZE_MAX - 1U) ||
      (size_t)plan->local_dim + 1U > SIZE_MAX / sizeof(*plan->row_ptr)) {
    goto fail;
  }
  row_ptr_bytes = ((size_t)plan->local_dim + 1U) * sizeof(*plan->row_ptr);
  plan->row_ptr = (size_t *)calloc((size_t)plan->local_dim + 1U,
                                  sizeof(*plan->row_ptr));
  if (plan->row_ptr == NULL) goto fail;

  StartTimer(1120);
  for (local_row = 0UL; local_row < plan->local_dim; local_row++) {
    unsigned long int alpha = plan->local_offset + local_row + 1UL;
    struct CountEntriesContext count = {0U};
    if (SymmetryEnumerateColumn(X, alpha, count_entry, &count) != 0) {
      StopTimer(1120);
      goto fail;
    }
    if (plan->row_ptr[local_row] > SIZE_MAX - count.count) {
      StopTimer(1120);
      goto fail;
    }
    plan->row_ptr[local_row + 1UL] = plan->row_ptr[local_row] + count.count;
    if (count.count < min_row_nnz) min_row_nnz = count.count;
    if (count.count > max_row_nnz) max_row_nnz = count.count;
  }
  StopTimer(1120);
  plan->nnz = plan->row_ptr[plan->local_dim];

  StartTimer(1121);
  if (plan->nnz > SIZE_MAX / sizeof(*plan->col_index) ||
      plan->nnz > SIZE_MAX / sizeof(*plan->values)) {
    StopTimer(1121);
    goto fail;
  }
  col_bytes = plan->nnz * sizeof(*plan->col_index);
  value_bytes = plan->nnz * sizeof(*plan->values);
  if (plan->nnz > 0U) {
    plan->col_index = (unsigned long int *)malloc(col_bytes);
    plan->values = (double complex *)malloc(value_bytes);
    if (plan->col_index == NULL || plan->values == NULL) {
      StopTimer(1121);
      goto fail;
    }
  }
  StopTimer(1121);

  StartTimer(1122);
  for (local_row = 0UL; local_row < plan->local_dim; local_row++) {
    unsigned long int alpha = plan->local_offset + local_row + 1UL;
    struct FillEntriesContext fill;
    fill.plan = plan;
    fill.next = plan->row_ptr[local_row];
    fill.end = plan->row_ptr[local_row + 1UL];
    if (SymmetryEnumerateColumn(X, alpha, fill_transposed_entry, &fill) != 0 ||
        fill.next != fill.end) {
      StopTimer(1122);
      goto fail;
    }
  }
  StopTimer(1122);

  if (row_ptr_bytes > SIZE_MAX - col_bytes ||
      row_ptr_bytes + col_bytes > SIZE_MAX - value_bytes) {
    goto fail;
  }
  plan_bytes = row_ptr_bytes + col_bytes + value_bytes;
  if (plan->local_dim == 0UL) min_row_nnz = 0U;
  plan->ready = TRUE;
  X->Sym->matvec_plan = plan;
  fprintf(stdoutMPI,
          "Symmetry matvec: mode=plan local_rows=%lu local_nnz=%zu "
          "row_nnz_min=%zu row_nnz_max=%zu row_nnz_mean=%.3f bytes=%zu.\n",
          plan->local_dim, plan->nnz, min_row_nnz, max_row_nnz,
          plan->local_dim > 0UL ? (double)plan->nnz / (double)plan->local_dim : 0.0,
          plan_bytes);
  return 0;

fail:
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
