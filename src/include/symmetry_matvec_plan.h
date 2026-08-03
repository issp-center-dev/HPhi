#ifndef HPHI_SYMMETRY_MATVEC_PLAN_H
#define HPHI_SYMMETRY_MATVEC_PLAN_H

#include <stdint.h>
#include "Common.h"
#include "symmetry_memory_policy.h"
#include "symmetry_vector_halo.h"

struct BindStruct;
struct SymmetryRepresentativeBatchOptions;

#define SYMMETRY_MATVEC_MODE_PLAN 0
#define SYMMETRY_MATVEC_MODE_LEGACY 1

#define SYMMETRY_VECTOR_EXCHANGE_ALLGATHER 0
#define SYMMETRY_VECTOR_EXCHANGE_HALO 1

#ifndef HPHI_SYMMETRY_PLAN_BLOCK_MEMORY_BYTES
#define HPHI_SYMMETRY_PLAN_BLOCK_MEMORY_BYTES \
  HPHI_SYMMETRY_MEMORY_LIMIT_BYTES
#endif

#ifndef HPHI_SYMMETRY_PLAN_LOCAL_ROWS_PER_BLOCK
#define HPHI_SYMMETRY_PLAN_LOCAL_ROWS_PER_BLOCK UINT64_C(65536)
#endif

enum SymmetryColumnWidth {
  SYMMETRY_COLUMN_U32 = 32,
  SYMMETRY_COLUMN_U64 = 64
};

struct SymmetryUnresolvedTransition {
  unsigned long int rep_state;
  /*
   * Off diagonal: hval * canonical phase, in production operation order.
   * Diagonal: the real owned diagonal promoted to double complex.
   */
  double complex phased_hval;
  int is_diagonal;
};

struct SymmetryUnresolvedMatvecBlock {
  unsigned long int local_row_begin;
  unsigned long int local_row_count;
  size_t transition_count;
  size_t offdiagonal_count;
  size_t *row_ptr;
  struct SymmetryUnresolvedTransition *transitions;
  size_t request_count;
  unsigned long int *request_keys;
  size_t storage_bytes;
  size_t temporary_peak_bytes;
  uint64_t memory_byte_limit;
};

struct SymmetryDistributedMatvecPlanOptions {
  unsigned long int local_rows_per_block;
  uint64_t block_memory_byte_limit;
  /* NULL selects the production directory batch path. */
  const struct SymmetryRepresentativeBatchOptions *directory_options;
};

struct SymmetryMatvecBlock {
  unsigned long int local_row_begin;
  unsigned long int local_row_count;
  size_t nnz;
  size_t *row_ptr;
  unsigned long int *global_columns;
  uint32_t *column_slot32;
  uint64_t *column_slot64;
  double complex *values;
};

struct SymmetryMatvecPlan {
  int ready;
  int columns_remapped;
  unsigned long int dim;
  unsigned long int local_offset;
  unsigned long int local_dim;
  size_t block_count;
  struct SymmetryMatvecBlock *blocks;
  size_t nnz;
  size_t row_nnz_max;
  /*
   * Non-owning compatibility aliases for a single-block plan.
   * Multi-block consumers must use SymmetryMatvecPlanGetBlockView().
   */
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
  size_t build_local_wave_count;
  size_t build_max_wave_count;
  size_t build_temporary_peak_bytes;
  uint64_t build_memory_warning_byte_threshold;
  uint64_t build_memory_byte_limit;
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
/**
 * Build a communication-free candidate block from distributed owned rows.
 *
 * Row order and Hamiltonian term order are preserved in transitions.
 * request_keys contains only off-diagonal representative keys, sorted and
 * deduplicated for one later collective directory batch.  memory_byte_limit
 * covers all block-owned storage plus the row-count construction workspace.
 * block_out must be empty and remains empty on failure.
 */
int BuildSymmetryUnresolvedMatvecBlock(
    const struct BindStruct *X,
    unsigned long int local_row_begin,
    unsigned long int local_row_count,
    uint64_t memory_byte_limit,
    struct SymmetryUnresolvedMatvecBlock *block_out);
void FreeSymmetryUnresolvedMatvecBlock(
    struct SymmetryUnresolvedMatvecBlock *block);
/**
 * Collectively resolve distributed transition blocks into global-column CSR.
 *
 * Every rank participates once per rank-local block wave up to the maximum
 * local block count, including ranks that have no block in a wave.  The
 * resulting plan has global columns and no halo/remap state; B4-C4 owns that
 * staged integration.
 */
int BuildSymmetryDistributedMatvecPlanWithOptions(
    struct BindStruct *X,
    const struct SymmetryDistributedMatvecPlanOptions *options);
int BuildSymmetryDistributedMatvecPlanForSolverWithOptions(
    struct BindStruct *X,
    const struct SymmetryDistributedMatvecPlanOptions *options);
int BuildSymmetryDistributedMatvecPlan(struct BindStruct *X);
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
