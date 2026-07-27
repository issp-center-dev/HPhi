/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*-------------------------------------------------------------*/
/*-------------------------------------------------------------
 * HPhi
 * timer program
 * "-lrt" option is needed for clock_gettime().
 *-------------------------------------------------------------
 * original code written by Satoshi Morita
 *-------------------------------------------------------------*/

#include "Common.h"
#include "FileIO.h"
#include "CalcTime.h"
#include "symmetry_basis.h"
#include "symmetry_directory.h"
#include "symmetry_matvec_plan.h"

#ifdef MPI
#include <mpi.h>
#endif

#ifdef MPI
static void OutputSymmetryRankStats(const struct BindStruct *X)
{
  static const int timer_ids[] = {
    1100, 1110, 1115, 1111, 1112, 1123, 1113, 1114,
    1101, 1120, 1121, 1122, 1130, 1131, 1132,
    1133, 1134, 1135, 1124, 1125, 4113,
    1, 1501, 1502, 1503, 1510, 1511, 1512, 1513
  };
  static const char *work_keys[] = {
    "basis_raw_states",
    "basis_representative_candidates",
    "basis_compatible_survivors",
    "basis_transform_calls",
    "basis_orbit_metadata_calls",
    "basis_thread_count",
    "basis_thread_raw_states_max",
    "basis_thread_representative_candidates_max",
    "basis_thread_compatible_survivors_max",
    "basis_thread_transform_calls_max",
    "basis_gather_entries",
    "basis_gather_bytes",
    "plan_local_rows",
    "plan_local_nnz",
    "plan_row_nnz_max",
    "plan_local_column_nnz",
    "plan_remote_column_nnz",
    "halo_ghost_count",
    "halo_send_value_count",
    "halo_incoming_peer_count",
    "halo_outgoing_peer_count",
    "halo_max_recv_from_peer",
    "halo_max_send_to_peer",
    "halo_schedule_bytes",
    "halo_runtime_buffer_bytes",
    "topology_scratch_bytes",
    "column_slot_width",
    "input_allgather_nonlocal_values_per_call",
    "input_allgather_payload_bytes_per_call",
    "symmetry_matvec_calls",
    "input_allgather_calls",
    "prdct_allreduce_calls",
    "halo_schedule_ready",
    "halo_reference_enabled",
    "halo_reference_exchange_calls",
    "columns_remapped",
    "full_input_vector_allocated",
    "halo_exchange_calls",
    "allocation_raw_basis_list_elements",
    "allocation_raw_diagonal_elements",
    "allocation_initial_vector_elements",
    "allocation_mpi_vector_buffer_elements",
    "allocation_auxiliary_vector_elements",
    "allocation_lobpcg_workspace_elements",
    "basis_state_enumerator_calls",
    "basis_diagonal_evaluator_calls",
    "local_basis_capacity_entries",
    "local_basis_bytes",
    "rank_offset_entries",
    "distribution_local_survivor_entries",
    "distribution_local_sample_entries",
    "distribution_global_sample_entries",
    "distribution_range_send_entries",
    "distribution_range_recv_entries",
    "distribution_rebalance_send_entries",
    "distribution_rebalance_recv_entries",
    "distribution_global_entries",
    "distribution_range_entries",
    "distribution_sample_gather_used_chunked",
    "distribution_sample_gather_message_byte_limit",
    "distribution_sample_gather_max_message_bytes",
    "distribution_range_exchange_used_chunked",
    "distribution_range_exchange_message_byte_limit",
    "distribution_range_exchange_max_message_bytes",
    "distribution_rebalance_exchange_used_chunked",
    "distribution_rebalance_exchange_message_byte_limit",
    "distribution_rebalance_exchange_max_message_bytes",
    "distribution_memory_byte_limit",
    "distribution_sort_temporary_peak_bytes",
    "distribution_rebalance_temporary_peak_bytes",
    "directory_nonempty_rank_count",
    "directory_splitter_entries",
    "directory_splitter_bytes",
    "directory_local_hash_entries",
    "directory_local_hash_table_entries",
    "directory_local_hash_bytes",
    "directory_hash_build_collisions",
    "directory_hash_build_max_probe",
    "directory_batch_calls",
    "directory_request_entries_sent",
    "directory_request_entries_received",
    "directory_found_entries",
    "directory_not_found_entries",
    "directory_lookup_probe_count",
    "directory_lookup_max_probe",
    "directory_owner_peer_count_max",
    "directory_requester_peer_count_max",
    "directory_exchange_used_chunked",
    "directory_exchange_message_byte_limit",
    "directory_exchange_max_message_bytes",
    "directory_exchange_send_messages",
    "directory_exchange_recv_messages",
    "directory_batch_temporary_peak_bytes",
    "directory_batch_memory_byte_limit",
    "plan_matrix_storage_bytes",
    "plan_column_storage_bytes",
    "directory_build_heavy_bytes",
    "directory_steady_heavy_bytes",
    "directory_heavy_storage_released",
    "plan_local_wave_count",
    "plan_max_wave_count"
  };
  static const char *metric_keys[] = {
    "plan_remote_column_nnz_ratio",
    "halo_ghost_global_ratio",
    "halo_ghost_nonlocal_ratio",
    "symmetry_matvec_seconds_per_call",
    "input_allgather_seconds_per_call",
    "input_allgather_effective_bandwidth_Bps",
    "plan_apply_seconds_per_call",
    "prdct_allreduce_seconds_per_call",
    "matvec_other_seconds",
    "matvec_other_seconds_per_call",
    "halo_reference_pack_seconds_per_call",
    "halo_reference_exchange_seconds_per_call",
    "halo_reference_validation_seconds_per_call",
    "halo_pack_seconds_per_call",
    "halo_exchange_seconds_per_call"
  };
  const size_t timer_count = sizeof(timer_ids) / sizeof(timer_ids[0]);
  const size_t work_count = sizeof(work_keys) / sizeof(work_keys[0]);
  const size_t metric_count = sizeof(metric_keys) / sizeof(metric_keys[0]);
  const struct SymmetryMatvecPlan *plan;
  double timer_local[sizeof(timer_ids) / sizeof(timer_ids[0])];
  double timer_min[sizeof(timer_ids) / sizeof(timer_ids[0])];
  double timer_max[sizeof(timer_ids) / sizeof(timer_ids[0])];
  double timer_sum[sizeof(timer_ids) / sizeof(timer_ids[0])];
  unsigned long long work_local[sizeof(work_keys) / sizeof(work_keys[0])];
  unsigned long long work_min[sizeof(work_keys) / sizeof(work_keys[0])];
  unsigned long long work_max[sizeof(work_keys) / sizeof(work_keys[0])];
  unsigned long long work_sum[sizeof(work_keys) / sizeof(work_keys[0])];
  double metric_local[sizeof(metric_keys) / sizeof(metric_keys[0])];
  double metric_min[sizeof(metric_keys) / sizeof(metric_keys[0])];
  double metric_max[sizeof(metric_keys) / sizeof(metric_keys[0])];
  double metric_sum[sizeof(metric_keys) / sizeof(metric_keys[0])];
  double row_mean_local;
  double row_mean_min;
  double row_mean_max;
  double row_mean_sum;
  struct SymmetryBasisDigest basis_digest;
  struct SymmetryRepresentativeDirectoryInfo directory_info;
  struct SymmetryLocalRepresentativeIndexStats directory_index_stats;
  struct SymmetryRepresentativeBatchStats directory_batch_stats;
  int basis_digest_status_local;
  int basis_digest_status_global;
  int basis_digest_algorithm_local;
  int basis_digest_algorithm_min;
  int basis_digest_algorithm_max;
  unsigned long long basis_digest_count_local;
  unsigned long long basis_digest_count_sum;
  unsigned long long basis_digest_local;
  unsigned long long basis_digest_min;
  unsigned long long basis_digest_max;
  unsigned long long basis_digest_xor_local;
  unsigned long long basis_digest_xor;
  unsigned long long basis_digest_sum_local;
  unsigned long long basis_digest_sum;
  unsigned long long halo_checksum_local;
  unsigned long long halo_checksum_xor;
  unsigned long long halo_checksum_sum;
  char fileName[D_FileNameMax];
  FILE *fp;
  size_t i;

  if (X == NULL || X->Def.iFlgSymmetryBasis == FALSE || X->Sym == NULL ||
      X->Sym->enabled != TRUE) {
    return;
  }
  plan = X->Sym->matvec_plan;
  memset(&directory_info, 0, sizeof(directory_info));
  memset(&directory_index_stats, 0, sizeof(directory_index_stats));
  memset(&directory_batch_stats, 0, sizeof(directory_batch_stats));
  if (X->Sym->representative_directory != NULL) {
    if (GetSymmetryRepresentativeDirectoryInfo(
            X->Sym->representative_directory, &directory_info) != 0 ||
        GetSymmetryRepresentativeDirectoryLocalIndexStats(
            X->Sym->representative_directory,
            &directory_index_stats) != 0 ||
        GetSymmetryRepresentativeDirectoryBatchStats(
            X->Sym->representative_directory,
            &directory_batch_stats) != 0) {
      memset(&directory_info, 0, sizeof(directory_info));
      memset(&directory_index_stats, 0, sizeof(directory_index_stats));
      memset(&directory_batch_stats, 0, sizeof(directory_batch_stats));
    }
  } else if (X->Sym->representative_directory_stats_ready == TRUE) {
    directory_info = X->Sym->representative_directory_info;
    directory_index_stats =
        X->Sym->representative_directory_index_stats;
    directory_batch_stats =
        X->Sym->representative_directory_batch_stats;
  }
  for (i = 0; i < timer_count; i++) timer_local[i] = Timer[timer_ids[i]];
  work_local[0] = X->Sym->basis_raw_states;
  work_local[1] = X->Sym->basis_representative_candidates;
  work_local[2] = X->Sym->basis_compatible_survivors;
  work_local[3] = X->Sym->basis_transform_calls;
  work_local[4] = X->Sym->basis_orbit_metadata_calls;
  work_local[5] = (unsigned long long)X->Sym->basis_thread_count;
  work_local[6] = X->Sym->basis_thread_raw_states_max;
  work_local[7] = X->Sym->basis_thread_representative_candidates_max;
  work_local[8] = X->Sym->basis_thread_compatible_survivors_max;
  work_local[9] = X->Sym->basis_thread_transform_calls_max;
  work_local[10] = X->Sym->basis_gather_entries;
  work_local[11] = X->Sym->basis_gather_bytes;
  work_local[12] = plan != NULL ? (unsigned long long)plan->local_dim : 0ULL;
  work_local[13] = plan != NULL ? (unsigned long long)plan->nnz : 0ULL;
  work_local[14] = plan != NULL ? (unsigned long long)plan->row_nnz_max : 0ULL;
  work_local[15] =
      plan != NULL ? (unsigned long long)plan->local_column_nnz : 0ULL;
  work_local[16] =
      plan != NULL ? (unsigned long long)plan->remote_column_nnz : 0ULL;
  work_local[17] =
      plan != NULL ? (unsigned long long)plan->halo.ghost_count : 0ULL;
  work_local[18] =
      plan != NULL ? (unsigned long long)plan->halo.send_value_count : 0ULL;
  work_local[19] =
      plan != NULL ? (unsigned long long)plan->halo.incoming_peer_count : 0ULL;
  work_local[20] =
      plan != NULL ? (unsigned long long)plan->halo.outgoing_peer_count : 0ULL;
  work_local[21] =
      plan != NULL ? (unsigned long long)plan->halo.max_recv_from_peer : 0ULL;
  work_local[22] =
      plan != NULL ? (unsigned long long)plan->halo.max_send_to_peer : 0ULL;
  work_local[23] = plan != NULL
                       ? (unsigned long long)plan->halo.schedule_bytes
                       : 0ULL;
  work_local[24] =
      plan != NULL
          ? (unsigned long long)plan->halo.runtime_buffer_bytes
          : 0ULL;
  work_local[25] =
      plan != NULL
          ? (unsigned long long)plan->halo.topology_scratch_bytes
          : 0ULL;
  work_local[26] =
      plan != NULL ? (unsigned long long)plan->column_slot_width : 0ULL;
  work_local[27] =
      plan != NULL
          ? (unsigned long long)plan->allgather_nonlocal_values_per_call
          : 0ULL;
  work_local[28] =
      plan != NULL
          ? (unsigned long long)plan->allgather_payload_bytes_per_call
          : 0ULL;
  work_local[29] =
      plan != NULL ? plan->matvec_calls : 0ULL;
  work_local[30] =
      plan != NULL ? plan->input_allgather_calls : 0ULL;
  work_local[31] =
      plan != NULL ? plan->prdct_allreduce_calls : 0ULL;
  work_local[32] =
      plan != NULL && plan->halo.ready == TRUE ? 1ULL : 0ULL;
  work_local[33] =
      plan != NULL && plan->halo.reference_enabled == TRUE ? 1ULL : 0ULL;
  work_local[34] =
      plan != NULL ? plan->halo.reference_exchange_calls : 0ULL;
  work_local[35] =
      plan != NULL && plan->columns_remapped == TRUE ? 1ULL : 0ULL;
  work_local[36] = X->Sym->mpi_full_v1 != NULL ? 1ULL : 0ULL;
  work_local[37] =
      plan != NULL ? plan->halo.exchange_calls : 0ULL;
  work_local[38] = X->Sym->allocation_raw_basis_list_elements;
  work_local[39] = X->Sym->allocation_raw_diagonal_elements;
  work_local[40] = X->Sym->allocation_initial_vector_elements;
  work_local[41] = X->Sym->allocation_mpi_vector_buffer_elements;
  work_local[42] = X->Sym->allocation_auxiliary_vector_elements;
  work_local[43] = X->Sym->allocation_lobpcg_workspace_elements;
  work_local[44] = X->Sym->basis_state_enumerator_calls;
  work_local[45] = X->Sym->basis_diagonal_evaluator_calls;
  work_local[46] = (unsigned long long)X->Sym->local_capacity;
  work_local[47] =
      (unsigned long long)X->Sym->local_capacity *
      (unsigned long long)sizeof(*X->Sym->local_basis);
  work_local[48] =
      X->Sym->rank_offsets != NULL ? (unsigned long long)nproc + 1ULL : 0ULL;
  work_local[49] = X->Sym->distribution_stats.local_survivor_entries;
  work_local[50] = X->Sym->distribution_stats.local_sample_entries;
  work_local[51] = X->Sym->distribution_stats.global_sample_entries;
  work_local[52] = X->Sym->distribution_stats.range_send_entries;
  work_local[53] = X->Sym->distribution_stats.range_recv_entries;
  work_local[54] = X->Sym->distribution_stats.rebalance_send_entries;
  work_local[55] = X->Sym->distribution_stats.rebalance_recv_entries;
  work_local[56] = X->Sym->distribution_stats.global_entries;
  work_local[57] = X->Sym->distribution_stats.range_entries;
  work_local[58] =
      X->Sym->distribution_stats.sample_gather_used_chunked != FALSE
          ? 1ULL
          : 0ULL;
  work_local[59] =
      X->Sym->distribution_stats.sample_gather_message_byte_limit;
  work_local[60] =
      X->Sym->distribution_stats.sample_gather_max_message_bytes;
  work_local[61] =
      X->Sym->distribution_stats.range_exchange_used_chunked != FALSE
          ? 1ULL
          : 0ULL;
  work_local[62] =
      X->Sym->distribution_stats.range_exchange_message_byte_limit;
  work_local[63] =
      X->Sym->distribution_stats.range_exchange_max_message_bytes;
  work_local[64] =
      X->Sym->distribution_stats.rebalance_exchange_used_chunked != FALSE
          ? 1ULL
          : 0ULL;
  work_local[65] =
      X->Sym->distribution_stats.rebalance_exchange_message_byte_limit;
  work_local[66] =
      X->Sym->distribution_stats.rebalance_exchange_max_message_bytes;
  work_local[67] =
      X->Sym->distribution_stats.distribution_memory_byte_limit;
  work_local[68] =
      (unsigned long long)
          X->Sym->distribution_stats.sort_temporary_peak_bytes;
  work_local[69] =
      (unsigned long long)
          X->Sym->distribution_stats.rebalance_temporary_peak_bytes;
  work_local[70] =
      (unsigned long long)directory_info.nonempty_rank_count;
  work_local[71] =
      (unsigned long long)directory_info.nonempty_rank_count;
  work_local[72] =
      (unsigned long long)directory_info.splitter_bytes;
  work_local[73] =
      (unsigned long long)directory_info.local_dim;
  work_local[74] =
      (unsigned long long)directory_index_stats.table_size;
  work_local[75] =
      (unsigned long long)directory_index_stats.table_bytes;
  work_local[76] =
      (unsigned long long)directory_index_stats.build_collisions;
  work_local[77] =
      (unsigned long long)directory_index_stats.build_max_probe;
  work_local[78] =
      (unsigned long long)directory_batch_stats.directory_batch_calls;
  work_local[79] =
      (unsigned long long)
          directory_batch_stats.directory_request_entries_sent;
  work_local[80] =
      (unsigned long long)
          directory_batch_stats.directory_request_entries_received;
  work_local[81] =
      (unsigned long long)directory_batch_stats.directory_found_entries;
  work_local[82] =
      (unsigned long long)
          directory_batch_stats.directory_not_found_entries;
  work_local[83] =
      (unsigned long long)
          directory_batch_stats.directory_lookup_probe_count;
  work_local[84] =
      (unsigned long long)
          directory_batch_stats.directory_lookup_max_probe;
  work_local[85] =
      (unsigned long long)
          directory_batch_stats.directory_owner_peer_count_max;
  work_local[86] =
      (unsigned long long)
          directory_batch_stats.directory_requester_peer_count_max;
  work_local[87] =
      directory_batch_stats.directory_exchange_used_chunked != FALSE
          ? 1ULL
          : 0ULL;
  work_local[88] =
      directory_batch_stats.directory_exchange_message_byte_limit;
  work_local[89] =
      directory_batch_stats.directory_exchange_max_message_bytes;
  work_local[90] =
      directory_batch_stats.directory_exchange_send_messages;
  work_local[91] =
      directory_batch_stats.directory_exchange_recv_messages;
  work_local[92] =
      (unsigned long long)
          directory_batch_stats.directory_batch_temporary_peak_bytes;
  work_local[93] =
      (unsigned long long)
          directory_batch_stats.directory_batch_memory_byte_limit;
  work_local[94] =
      plan != NULL
          ? (unsigned long long)plan->matrix_storage_bytes
          : 0ULL;
  work_local[95] =
      plan != NULL
          ? (unsigned long long)plan->column_storage_bytes
          : 0ULL;
  work_local[96] =
      (unsigned long long)directory_info.splitter_bytes +
      (unsigned long long)directory_index_stats.table_bytes;
  work_local[97] =
      X->Sym->representative_directory != NULL
          ? work_local[96] : 0ULL;
  work_local[98] =
      X->Sym->representative_directory_heavy_storage_released == TRUE
          ? 1ULL : 0ULL;
  work_local[99] =
      plan != NULL
          ? (unsigned long long)plan->build_local_wave_count
          : 0ULL;
  work_local[100] =
      plan != NULL
          ? (unsigned long long)plan->build_max_wave_count
          : 0ULL;
  row_mean_local = plan != NULL && plan->local_dim > 0UL
                       ? (double)plan->nnz / (double)plan->local_dim
                       : 0.0;
  metric_local[0] = plan != NULL && plan->nnz > 0U
                        ? (double)plan->remote_column_nnz / (double)plan->nnz
                        : 0.0;
  metric_local[1] = plan != NULL && plan->dim > 0UL
                        ? (double)plan->halo.ghost_count / (double)plan->dim
                        : 0.0;
  metric_local[2] =
      plan != NULL && plan->allgather_nonlocal_values_per_call > 0U
          ? (double)plan->halo.ghost_count /
                (double)plan->allgather_nonlocal_values_per_call
          : 0.0;
  metric_local[3] =
      plan != NULL && plan->matvec_calls > 0ULL
          ? Timer[1] / (double)plan->matvec_calls
          : 0.0;
  metric_local[4] =
      plan != NULL && plan->input_allgather_calls > 0ULL
          ? Timer[1501] / (double)plan->input_allgather_calls
          : 0.0;
  metric_local[5] =
      plan != NULL && plan->input_allgather_calls > 0ULL &&
              Timer[1501] > 0.0
          ? ((double)plan->allgather_payload_bytes_per_call *
             (double)plan->input_allgather_calls) /
                Timer[1501]
          : 0.0;
  metric_local[6] =
      plan != NULL && plan->matvec_calls > 0ULL
          ? Timer[1503] / (double)plan->matvec_calls
          : 0.0;
  metric_local[7] =
      plan != NULL && plan->prdct_allreduce_calls > 0ULL
          ? Timer[1513] / (double)plan->prdct_allreduce_calls
          : 0.0;
  metric_local[8] =
      Timer[1] - Timer[1501] - Timer[1502] - Timer[1503] -
      Timer[1510] - Timer[1511] - Timer[1512] - Timer[1513];
  metric_local[9] =
      plan != NULL && plan->matvec_calls > 0ULL
          ? metric_local[8] / (double)plan->matvec_calls
          : 0.0;
  metric_local[10] =
      plan != NULL && plan->halo.reference_exchange_calls > 0ULL
          ? Timer[1510] /
                (double)plan->halo.reference_exchange_calls
          : 0.0;
  metric_local[11] =
      plan != NULL && plan->halo.reference_exchange_calls > 0ULL
          ? Timer[1511] /
                (double)plan->halo.reference_exchange_calls
          : 0.0;
  metric_local[12] =
      plan != NULL && plan->halo.reference_exchange_calls > 0ULL
          ? Timer[1512] /
                (double)plan->halo.reference_exchange_calls
          : 0.0;
  metric_local[13] =
      plan != NULL && plan->halo.exchange_calls > 0ULL
          ? Timer[1510] / (double)plan->halo.exchange_calls
          : 0.0;
  metric_local[14] =
      plan != NULL && plan->halo.exchange_calls > 0ULL
          ? Timer[1511] / (double)plan->halo.exchange_calls
          : 0.0;
  basis_digest_status_local =
      ComputeSymmetryBasisDigest(X->Sym, &basis_digest) == 0 ? 0 : 1;
  basis_digest_algorithm_local = (int)basis_digest.algorithm;
  basis_digest_count_local = basis_digest.count;
  basis_digest_local = basis_digest.fnv1a64;
  basis_digest_xor_local = basis_digest.xor_hash;
  basis_digest_sum_local = basis_digest.sum_hash;
  halo_checksum_local =
      plan != NULL ? plan->halo.schedule_checksum : 0ULL;

  MPI_Allreduce(timer_local, timer_min, (int)timer_count, MPI_DOUBLE,
                MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(timer_local, timer_max, (int)timer_count, MPI_DOUBLE,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(timer_local, timer_sum, (int)timer_count, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(work_local, work_min, (int)work_count, MPI_UNSIGNED_LONG_LONG,
                MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(work_local, work_max, (int)work_count, MPI_UNSIGNED_LONG_LONG,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(work_local, work_sum, (int)work_count, MPI_UNSIGNED_LONG_LONG,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(metric_local, metric_min, (int)metric_count, MPI_DOUBLE,
                MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(metric_local, metric_max, (int)metric_count, MPI_DOUBLE,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(metric_local, metric_sum, (int)metric_count, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&row_mean_local, &row_mean_min, 1, MPI_DOUBLE,
                MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&row_mean_local, &row_mean_max, 1, MPI_DOUBLE,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&row_mean_local, &row_mean_sum, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&basis_digest_status_local, &basis_digest_status_global,
                1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&basis_digest_algorithm_local, &basis_digest_algorithm_min,
                1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&basis_digest_algorithm_local, &basis_digest_algorithm_max,
                1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&basis_digest_count_local, &basis_digest_count_sum, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&basis_digest_local, &basis_digest_min, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&basis_digest_local, &basis_digest_max, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&basis_digest_xor_local, &basis_digest_xor, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_BXOR, MPI_COMM_WORLD);
  MPI_Allreduce(&basis_digest_sum_local, &basis_digest_sum, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&halo_checksum_local, &halo_checksum_xor, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_BXOR, MPI_COMM_WORLD);
  MPI_Allreduce(&halo_checksum_local, &halo_checksum_sum, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

  sprintf(fileName, "CalcTimerRankStats.dat");
  if (childfopenMPI(fileName, "w", &fp) != 0) return;
  fprintf(fp,
          "format=HPhiCalcTimerRankStats version=9 ranks=%d "
          "basis_layout=%s matvec_mode=%s vector_exchange=%s\n",
          nproc,
          X->Sym->basis_layout == SYMMETRY_BASIS_DISTRIBUTED
              ? "distributed"
              : "replicated",
          X->Sym->matvec_mode == SYMMETRY_MATVEC_MODE_PLAN
              ? "plan"
              : "legacy",
          X->Sym->vector_exchange_mode == SYMMETRY_VECTOR_EXCHANGE_HALO
              ? "halo"
              : "allgather");
  for (i = 0; i < timer_count; i++) {
    fprintf(fp, "timer id=%d ranks=%d min=%.17g max=%.17g mean=%.17g\n",
            timer_ids[i], nproc, timer_min[i], timer_max[i],
            timer_sum[i] / (double)nproc);
  }
  for (i = 0; i < work_count; i++) {
    fprintf(fp, "work key=%s ranks=%d min=%llu max=%llu mean=%.17g\n",
            work_keys[i], nproc, work_min[i], work_max[i],
            (double)work_sum[i] / (double)nproc);
  }
  for (i = 0; i < metric_count; i++) {
    fprintf(fp, "metric key=%s ranks=%d min=%.17g max=%.17g mean=%.17g\n",
            metric_keys[i], nproc, metric_min[i], metric_max[i],
            metric_sum[i] / (double)nproc);
  }
  fprintf(fp,
          "work key=plan_row_nnz_mean ranks=%d min=%.17g max=%.17g mean=%.17g\n",
          nproc, row_mean_min, row_mean_max, row_mean_sum / (double)nproc);
  if (basis_digest_status_global != 0 ||
      basis_digest_algorithm_min != basis_digest_algorithm_max) {
    fprintf(fp,
            "basis_digest algorithm=invalid ranks=%d status=error\n",
            nproc);
  } else if (basis_digest_algorithm_min ==
             SYMMETRY_BASIS_DIGEST_REPLICATED_FNV1A64) {
    fprintf(fp,
            "basis_digest algorithm=fnv1a64-fields ranks=%d "
            "min=%016llx max=%016llx count=%llu status=ok\n",
            nproc, basis_digest_min, basis_digest_max,
            basis_digest_count_local);
  } else if (basis_digest_algorithm_min ==
             SYMMETRY_BASIS_DIGEST_DISTRIBUTED_GLOBAL_BETA) {
    fprintf(fp,
            "basis_digest algorithm=fnv1a64-global-beta-fields-xor-sum "
            "ranks=%d xor=%016llx sum=%016llx count=%llu status=ok\n",
            nproc, basis_digest_xor, basis_digest_sum,
            basis_digest_count_sum);
  } else {
    fprintf(fp,
            "basis_digest algorithm=invalid ranks=%d status=error\n",
            nproc);
  }
  fprintf(fp,
          "halo_schedule_digest algorithm=fnv1a64-layout ranks=%d "
          "xor=%016llx sum=%016llx\n",
          nproc, halo_checksum_xor, halo_checksum_sum);
  fclose(fp);
}
#endif
/** 
 * 
 * @brief function for displaying elapse time
 * 
 * @version 2.0
 */

void StampTime(FILE *fp, char *str, int num){
#ifdef MPI
  char str1[256];
  sprintf(str1, "%-50s [%04d] %12.5lf\n", str, num, Timer[num]);
  fprintf(fp, "%s", str1);
#endif
}

/** 
 * 
 * @brief function for initializing Timer[]
 * 
 * @version 2.0
 */
void InitTimer() {
#ifdef MPI
  int i;
  int NTimer=10000;
  Timer       = (double*)malloc((NTimer)*sizeof(double));
  TimerStart  = (double*)malloc((NTimer)*sizeof(double));
  for(i=0;i<NTimer;i++) Timer[i]=0.0;
  for(i=0;i<NTimer;i++) TimerStart[i]=0.0;
#endif
  return;
}
/** 
 * 
 * @brief function for initializing elapse time [start]
 * 
 * @version 2.0
 */

void StartTimer(int n) {
#ifdef MPI
  TimerStart[n]=MPI_Wtime();
#endif
  return;
}
/** 
 * 
 * @brief function for calculating elapse time [elapse time=StartTimer-StopTimer]
 * 
 * @version 2.0
 */
void StopTimer(int n) {
#ifdef MPI
  Timer[n] += MPI_Wtime() - TimerStart[n];
#endif
  return;
}
/** 
 * 
 * @brief function for outputting elapse time for each function
 * 
 * @version 2.0
 */
void OutputTimer(struct BindStruct *X) {

#ifdef MPI
  char fileName[D_FileNameMax];
  FILE *fp;
  sprintf(fileName, "CalcTimer.dat"); //TBC
  childfopenMPI(fileName,"w", &fp);
  //fp = fopen(fileName, "w");
  //fp = childfopenMPI(fileName, "w");
  StampTime(fp, "All", 0);
  StampTime(fp, "  sz", 1000);
  StampTime(fp, "  symmetry basis build/activate", 1100);
  StampTime(fp, "    symmetry basis raw enumeration", 1110);
  StampTime(fp, "      symmetry basis MPI gather/reduction", 1115);
  StampTime(fp, "    symmetry basis sort/merge", 1111);
  StampTime(fp, "    symmetry representative hash build", 1112);
  StampTime(fp, "    symmetry distributed sample sort", 1133);
  StampTime(fp, "    symmetry distributed exact rebalance", 1134);
  StampTime(fp, "    symmetry distributed storage validation", 1135);
  StampTime(fp, "    symmetry representative batch directory", 1123);
  StampTime(fp, "    symmetry plan unresolved block build", 1124);
  StampTime(fp, "    symmetry plan block finalize", 1125);
  StampTime(fp, "  symmetry solver storage allocation", 1113);
  StampTime(fp, "    symmetry dimension activation/validation", 1114);
  StampTime(fp, "  symmetry matvec plan build", 1101);
  StampTime(fp, "    symmetry plan count/prefix", 1120);
  StampTime(fp, "    symmetry plan storage allocation", 1121);
  StampTime(fp, "    symmetry plan fill", 1122);
  StampTime(fp, "    symmetry plan topology extraction", 1130);
  StampTime(fp, "    symmetry halo request schedule", 1131);
  StampTime(fp, "    symmetry CSR column remap", 1132);
  StampTime(fp, "  diagonalcalc", 2000);
  if(X->Def.iFlgCalcSpec == CALCSPEC_NOT){
    if(X->Def.iCalcType==TPQCalc || X->Def.iCalcType==cTPQ) {
      StampTime(fp, "  CalcByTPQ", 3000);
      StampTime(fp, "    FirstMultiply", 3100);
      StampTime(fp, "      rand   in FirstMultiply", 3101);
      StampTime(fp, "      mltply in FirstMultiply", 3102);
      StampTime(fp, "    expec_energy_flct        ", 3200);
      StampTime(fp, "      calc flctuation in expec_energy_flct ", 3201);
      StampTime(fp, "      mltply in expec_energy_flct ", 3202);
      StampTime(fp, "    expec_onebody            ", 3300);
      StampTime(fp, "    expec_twobody            ", 3400);
      StampTime(fp, "    Multiply                 ", 3500);
      StampTime(fp, "    FileIO                   ", 3600);
    }
    else if(X->Def.iCalcType==Lanczos){
      StampTime(fp, "  CalcByLanczos", 4000);
      StampTime(fp, "    LanczosEigenValue", 4100);
      StampTime(fp, "      mltply      in LanczosEigenValue", 4101);
      StampTime(fp, "      vec12       in LanczosEigenValue", 4102);
      StampTime(fp, "      DSEVvalue   in LanczosEigenValue", 4103);
      StampTime(fp, "      initial vector zero fill", 4110);
      StampTime(fp, "      initial vector random fill", 4111);
      StampTime(fp, "      initial vector local norm", 4112);
      StampTime(fp, "      initial vector MPI reduction", 4113);
      StampTime(fp, "      initial vector normalization", 4114);
      StampTime(fp, "    LanczosEigenVector", 4200);
      StampTime(fp, "      mltply      in LanczosEigenVector", 4201);
      StampTime(fp, "    expec_energy_flct", 4300);
      StampTime(fp, "      calc flctuation in expec_energy_flct ", 4301);
      StampTime(fp, "      mltply in expec_energy_flct ", 4302);
      StampTime(fp, "    CGEigenVector", 4400);
      StampTime(fp, "      mltply in CGEigenVector ", 4401);
      StampTime(fp, "    expec_onebody            ", 4500);
      StampTime(fp, "    expec_twobody            ", 4600);
      StampTime(fp, "    expec_TotalSz            ", 4700);
      StampTime(fp, "    FileIO                   ", 4800);
      StampTime(fp, "      Read Input Eigenvec ", 4801);    
    }
    else if(X->Def.iCalcType==FullDiag){
      StampTime(fp, "  CalcByFullDiag", 5000);
      StampTime(fp, "    MakeHam", 5100);
      StampTime(fp, "    LapackDiag", 5200);
      StampTime(fp, "    CalcPhys", 5300);
    StampTime(fp, "      calc flctuation in expec_energy_flct ", 5301);
    StampTime(fp, "      mltply in expec_energy_flct ", 5302);
        StampTime(fp, "    Output", 5400);
      StampTime(fp, "    OutputHam", 5500);
    }
  }
  else{ 
    StampTime(fp, "  CalcSpectrum by Lanczos method", 6000);
    StampTime(fp, "    Make excited state", 6100);
    StampTime(fp, "      Read origin state", 6101);
    StampTime(fp, "      Multiply excited operator", 6102);
    StampTime(fp, "    Calculate spectrum", 6200);
    if(X->Def.iCalcType==Lanczos){
      StampTime(fp, "      Read vector for recalculation", 6201);
      StampTime(fp, "      Read tridiagonal components for recalculation", 6202);
      StampTime(fp, "      Calculate tridiagonal components", 6203);
      StampTime(fp, "      Output tridiagonal components", 6204);
      StampTime(fp, "      Calculate spectrum by Lanczos method", 6205);
      StampTime(fp, "      Output vectors for recalculation", 6206);
    }
    else if(X->Def.iCalcType==FullDiag){
      StampTime(fp, "      MakeHam", 6301);
      StampTime(fp, "      lapackdiag", 6302);
      StampTime(fp, "      Calculate v1", 6303);
      StampTime(fp, "      Calculate spectrum", 6304);
    }
  }
  
  fprintf(fp,"================================================\n");
  
  StampTime(fp,"All mltply",1);
  StampTime(fp,"  symmetry input Allgatherv",1501);
  StampTime(fp,"  symmetry legacy beta scan",1502);
  StampTime(fp,"  symmetry local-row plan apply",1503);
  StampTime(fp,"  symmetry halo pack",1510);
  StampTime(fp,"  symmetry halo exchange",1511);
  StampTime(fp,"  symmetry halo reference validation",1512);
  StampTime(fp,"  symmetry prdct scalar Allreduce",1513);
  StampTime(fp,"  diagonal", 100);

  switch(X->Def.iCalcModel){
  case HubbardGC:
    StampTime(fp,"  HubbardGC", 200);
    StampTime(fp,"    trans    in HubbardGC", 210);
    StampTime(fp,"      double", 211);
    StampTime(fp,"      single", 212);
    StampTime(fp,"      inner", 213);
    StampTime(fp,"    interall in HubbardGC", 220);
    StampTime(fp,"      interPE", 221);
    StampTime(fp,"      inner", 222);
    StampTime(fp,"    pairhopp in HubbardGC", 230);
    StampTime(fp,"      interPE", 231);
    StampTime(fp,"      inner", 232);
    StampTime(fp,"    exchange in HubbardGC", 240);
    StampTime(fp,"      interPE", 241);
    StampTime(fp,"      inner", 242);
    break;
    
  case Hubbard:
  case tJ:
  case tJGC:
    StampTime(fp,"  Hubbard", 300);
    StampTime(fp,"    trans    in Hubbard", 310);
    StampTime(fp,"      double", 311);
    StampTime(fp,"      single", 312);
    StampTime(fp,"      inner", 313);
    StampTime(fp,"    interall in Hubbard", 320);
    StampTime(fp,"      interPE", 321);
    StampTime(fp,"      inner", 322);
    StampTime(fp,"    pairhopp in Hubbard", 330);
    StampTime(fp,"      interPE", 331);
    StampTime(fp,"      inner", 332);
    StampTime(fp,"    exchange in Hubbard", 340);
    StampTime(fp,"      interPE", 341);
    StampTime(fp,"      inner", 342);
    break;
    
  case Spin:
    fprintf(fp,"\n");
    StampTime(fp,"  Spin", 400);
    StampTime(fp,"    interall in Spin", 410);
    StampTime(fp,"      double", 411);
    StampTime(fp,"      single1", 412);
    StampTime(fp,"      single2", 413);
    StampTime(fp,"      inner", 414);
    StampTime(fp,"    exchange in Spin", 420);
    StampTime(fp,"      double", 421);
    StampTime(fp,"      single1", 422);
    StampTime(fp,"      single2", 423);
    StampTime(fp,"      inner", 424);
    break;
    
  case SpinGC:
    StampTime(fp,"  SpinGC", 500);
    StampTime(fp,"    trans    in SpinGC", 510);
    StampTime(fp,"      double", 511);
    StampTime(fp,"      inner", 512);
    StampTime(fp,"    interall in SpinGC", 520);
    StampTime(fp,"      double", 521);
    StampTime(fp,"      single", 522);
    StampTime(fp,"      inner", 523);
    StampTime(fp,"    exchange in SpinGC", 530);
    StampTime(fp,"      double", 531);
    StampTime(fp,"      single", 532);
    StampTime(fp,"      inner", 533);
    StampTime(fp,"    pairlift in SpinGC", 540);
    StampTime(fp,"      double", 541);
    StampTime(fp,"      single", 542);
    StampTime(fp,"      inner", 543);
    break;

  default:
    break;
  }
  fprintf(fp,"================================================\n");

  fclose(fp);
  OutputSymmetryRankStats(X);
  free(Timer);
  free(TimerStart);
#endif
  return;
}

/**
@page page_time Compute elapsed time for new functions

 Using StartTimer and StopTimer functions defined in time.c, we can measure the elapsed time for computation.

 1. Define an index and an output message by using StampTime function in time.c.
 For example, the index and the output message for the elapsed time of TPQ calculation is defined as follows.
 ```
       StampTime(fp, "  CalcByTPQ", 3000);
 ```

 2. Include CalcTime.h in the target source file.

 3. Set StartTimer and StopTimer functions in the region where you want to measure the time.
    It is noted that both functions must have the same index defined in time.c.
 For example, the elapsed time of TPQ calculation can be measured as follows.
 ```
       case TPQCalc:
        StartTimer(3000);
        if (CalcByTPQ(NumAve, X.Bind.Def.Param.ExpecInterval, &X) != TRUE) {
          FinalizeMPI();
          StopTimer(3000);
          return 0;
        }
        StopTimer(3000);
      break;

 ```

When above procedures were done, after calculation, you can see the elapsed time in CalcTimer.dat file as follows. The time unit is second.

```
 All                                                [0000]     37.94046
  sz                                               [1000]      0.00058
  diagonalcalc                                     [2000]      0.00046
  CalcByTPQ                                        [3000]     37.93129
```


*/
