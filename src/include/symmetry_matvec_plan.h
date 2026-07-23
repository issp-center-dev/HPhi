#ifndef HPHI_SYMMETRY_MATVEC_PLAN_H
#define HPHI_SYMMETRY_MATVEC_PLAN_H

#include "Common.h"

struct BindStruct;

#define SYMMETRY_MATVEC_MODE_PLAN 0
#define SYMMETRY_MATVEC_MODE_LEGACY 1

struct SymmetryMatvecPlan {
  int ready;
  unsigned long int dim;
  unsigned long int local_offset;
  unsigned long int local_dim;
  size_t nnz;
  size_t row_nnz_max;
  size_t *row_ptr;
  unsigned long int *col_index;
  double complex *values;
  size_t local_column_nnz;
  size_t remote_column_nnz;
  size_t ghost_count;
  size_t send_value_count;
  size_t incoming_peer_count;
  size_t outgoing_peer_count;
  size_t max_recv_from_peer;
  size_t max_send_to_peer;
  size_t topology_scratch_bytes;
  size_t halo_schedule_bytes_estimate;
  size_t halo_runtime_buffer_bytes_estimate;
  size_t allgather_nonlocal_values_per_call;
  size_t allgather_payload_bytes_per_call;
  unsigned int column_slot_width;
  unsigned long long matvec_calls;
  unsigned long long input_allgather_calls;
  unsigned long long prdct_allreduce_calls;
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
int ApplySymmetryMatvecPlan(const struct BindStruct *X,
                            double complex *tmp_v0,
                            const double complex *full_v1,
                            double complex *local_prdct);
void FreeSymmetryMatvecPlan(struct SymmetryMatvecPlan *plan);

#endif /* HPHI_SYMMETRY_MATVEC_PLAN_H */
