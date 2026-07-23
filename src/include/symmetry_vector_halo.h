#ifndef HPHI_SYMMETRY_VECTOR_HALO_H
#define HPHI_SYMMETRY_VECTOR_HALO_H

#include "Common.h"

struct SymmetryVectorHaloPlan {
  int request_layout_ready;
  int ready;
  int reference_enabled;
  int nrank;
  int rank;
  unsigned long int dim;
  unsigned long int local_offset;
  unsigned long int local_dim;
  size_t ghost_count;
  size_t send_value_count;
  size_t incoming_peer_count;
  size_t outgoing_peer_count;
  size_t max_recv_from_peer;
  size_t max_send_to_peer;
  int *send_counts;
  int *send_displs;
  int *recv_counts;
  int *recv_displs;
  unsigned long int *ghost_global_index;
  unsigned long int *send_local_index;
  double complex *send_values;
  double complex *ghost_values;
  size_t topology_scratch_bytes;
  size_t schedule_bytes;
  size_t runtime_buffer_bytes;
  unsigned long long schedule_checksum;
  unsigned long long reference_exchange_calls;
};

int SymmetryVectorOwnerOfGlobalIndex(unsigned long int dim,
                                     int nrank,
                                     unsigned long int global_index);
int BuildSymmetryVectorHaloPlan(struct SymmetryVectorHaloPlan *halo,
                                unsigned long int dim,
                                unsigned long int local_offset,
                                unsigned long int local_dim,
                                const unsigned long int *global_columns,
                                size_t column_count,
                                int nrank,
                                int rank,
                                size_t *local_column_count,
                                size_t *remote_column_count);
int ExchangeSymmetryVectorHaloReference(
    struct SymmetryVectorHaloPlan *halo,
    const double complex *local_vector,
    const double complex *full_vector);
void FreeSymmetryVectorHaloPlan(struct SymmetryVectorHaloPlan *halo);

#endif /* HPHI_SYMMETRY_VECTOR_HALO_H */
