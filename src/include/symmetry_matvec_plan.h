#ifndef HPHI_SYMMETRY_MATVEC_PLAN_H
#define HPHI_SYMMETRY_MATVEC_PLAN_H

#include <stdint.h>
#include "Common.h"
#include "symmetry_vector_halo.h"

struct BindStruct;

#define SYMMETRY_MATVEC_MODE_PLAN 0
#define SYMMETRY_MATVEC_MODE_LEGACY 1

#define SYMMETRY_VECTOR_EXCHANGE_ALLGATHER 0
#define SYMMETRY_VECTOR_EXCHANGE_HALO 1

enum SymmetryColumnWidth {
  SYMMETRY_COLUMN_U32 = 32,
  SYMMETRY_COLUMN_U64 = 64
};

struct SymmetryMatvecPlan {
  int ready;
  int columns_remapped;
  unsigned long int dim;
  unsigned long int local_offset;
  unsigned long int local_dim;
  size_t block_count;
  size_t nnz;
  size_t row_nnz_max;
  size_t *row_ptr;
  /* Global 1-origin columns before halo remap and in allgather mode. */
  unsigned long int *col_index;
  /* Exactly one adaptive 0-origin slot array owns halo-mode columns. */
  uint32_t *column_slot32;
  uint64_t *column_slot64;
  double complex *values;
  size_t column_storage_bytes;
  size_t matrix_storage_bytes;
  size_t local_column_nnz;
  size_t remote_column_nnz;
  size_t allgather_nonlocal_values_per_call;
  size_t allgather_payload_bytes_per_call;
  enum SymmetryColumnWidth column_slot_width;
  unsigned long long matvec_calls;
  unsigned long long input_allgather_calls;
  unsigned long long prdct_allreduce_calls;
  struct SymmetryVectorHaloPlan halo;
};

struct SymmetryMatvecBlockView {
  unsigned long int local_row_begin;
  unsigned long int local_row_count;
  size_t nnz;
  const size_t *row_ptr;
  const unsigned long int *global_columns;
  const uint32_t *column_slot32;
  const uint64_t *column_slot64;
  const double complex *values;
};

/**
 * Receive one matrix entry H(out_index, beta) from SymmetryEnumerateColumn().
 * coefficient includes the canonical phase and norm[out_index]/norm[beta],
 * but does not include an input-vector amplitude.
 */
typedef int (*SymmetryEntryCallback)(unsigned long int out_index,
                                     double complex coefficient,
                                     void *context);

int SymmetryEnumerateColumn(const struct BindStruct *X,
                            unsigned long int beta,
                            SymmetryEntryCallback callback,
                            void *context);
int BuildSymmetryMatvecPlan(struct BindStruct *X);
size_t SymmetryMatvecPlanBlockCount(
    const struct SymmetryMatvecPlan *plan);
int SymmetryMatvecPlanGetBlockView(
    const struct SymmetryMatvecPlan *plan,
    size_t block_index,
    struct SymmetryMatvecBlockView *view);
int ApplySymmetryMatvecPlan(const struct BindStruct *X,
                            double complex *tmp_v0,
                            const double complex *full_v1,
                            double complex *local_prdct);
int ApplySymmetryMatvecPlanHalo(const struct BindStruct *X,
                                double complex *tmp_v0,
                                const double complex *local_v1,
                                double complex *local_prdct);
int RemapSymmetryMatvecPlanColumns(struct SymmetryMatvecPlan *plan);
void FreeSymmetryMatvecPlan(struct SymmetryMatvecPlan *plan);

#endif /* HPHI_SYMMETRY_MATVEC_PLAN_H */
