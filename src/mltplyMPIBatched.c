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

/**
 * @file mltplyMPIBatched.c
 * @brief Implementation of batched MPI communication optimization
 *
 * This module groups multiple transfer terms that share the same MPI communication
 * partner (origin) to reduce the number of MPI_Sendrecv calls.
 *
 * Indexing convention note:
 * - Wavefunction buffers (`tmp_v*`, `v1buf`, `list_1*`) are 1-based in many legacy paths.
 * - Bit-operation state ids (`jreal`, `ioff` before +1 adjustment) are 0-based.
 * Preserve each function's existing `j` loop range and `+1` offsets when modifying loops.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */

#ifdef MPI
#include "mpi.h"
#endif

#include <stdlib.h>
#include <string.h>
#include "Common.h"
#include "mltply.h"
#include "mltplyCommon.h"
#include "bitcalc.h"
#include "wrapperMPI.h"
#include "mltplyMPIHubbardCore.h"
#include "mltplyMPIBatched.h"

/* Returns 1 if MPI communication batching is enabled, 0 if disabled
   (HPHI_MPI_NOBATCH=1). Available in both MPI and non-MPI builds. */
int MPIBatchingEnabled(void) {
    return iFlgMPIBatch;
}

#ifdef MPI

/**
 * @brief Initialize batched transfers for SpinlessFermionGC MPIsingle mode
 *
 * Groups MPIsingle transfers (one site local, one inter-process) by their
 * MPI communication partner (origin). Each group shares a single MPI_Sendrecv.
 *
 * MPI origin computation:
 * - origin = myrank XOR mask2, where mask2 = Tpow[inter-process site]
 * - XOR flips the bit corresponding to the inter-process site, giving the
 *   rank that "owns" the complementary bit pattern.
 *
 * @param X       Bind structure with transfer definitions
 * @param batched Output: populated MPIBatchedTransfers structure
 * @return 0 on success, -1 on allocation failure
 */
int InitializeMPIBatchedTransfers_SpinlessFermionGC(
    struct BindStruct *X,
    MPIBatchedTransfers *batched
) {
    int i, g, t;
    int mask2, origin;
    int *origin_count;
    int *unique_origins;
    int num_unique = 0;
    int max_transfers;

    if (batched == NULL) return -1;

    batched->is_initialized = 0;
    batched->num_groups = 0;
    batched->groups = NULL;

    // Count maximum possible transfers for MPIsingle
    max_transfers = 0;
    for (i = 0; i < (int)X->Def.EDNTransfer; i += 2) {
        // Check if this is an MPIsingle transfer (one site local, one site inter-process)
        int site1_local = (X->Def.EDGeneralTransfer[i][0] + 1 <= (int)X->Def.Nsite);
        int site2_local = (X->Def.EDGeneralTransfer[i][2] + 1 <= (int)X->Def.Nsite);

        if (site1_local != site2_local) {
            max_transfers += 2;  // Both directions
        }
    }

    if (max_transfers == 0) {
        batched->is_initialized = 1;
        return 0;
    }

    // Allocate temporary arrays for counting
    origin_count = (int *)calloc(nproc, sizeof(int));
    unique_origins = (int *)malloc(nproc * sizeof(int));
    if (origin_count == NULL || unique_origins == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    // First pass: count transfers per origin and identify unique origins
    for (i = 0; i < (int)X->Def.EDNTransfer; i += 2) {
        int site1 = X->Def.EDGeneralTransfer[i][0];
        int site2 = X->Def.EDGeneralTransfer[i][2];
        int site1_local = (site1 + 1 <= (int)X->Def.Nsite);
        int site2_local = (site2 + 1 <= (int)X->Def.Nsite);

        // Only process MPIsingle transfers (exactly one site in inter-process region)
        if (site1_local == site2_local) continue;

        // Determine the inter-process site
        int inter_site = site1_local ? site2 : site1;

        mask2 = (int)X->Def.Tpow[inter_site];
        // Partner rank differs by the PE-bit of the inter-process site.
        // XOR flips exactly that bit to select the communication peer.
        origin = myrank ^ mask2;

        if (origin_count[origin] == 0) {
            unique_origins[num_unique++] = origin;
        }
        origin_count[origin]++;
    }

    if (num_unique == 0) {
        free(origin_count);
        free(unique_origins);
        batched->is_initialized = 1;
        return 0;
    }

    // Allocate groups
    batched->num_groups = num_unique;
    batched->groups = (MPITransferGroup *)calloc(num_unique, sizeof(MPITransferGroup));
    if (batched->groups == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    // Initialize each group
    for (g = 0; g < num_unique; g++) {
        origin = unique_origins[g];
        int count = origin_count[origin];

        batched->groups[g].origin = origin;
        batched->groups[g].num_transfers = count;
        batched->groups[g].mask = 0;  // Will be set below
        batched->groups[g].transfer_indices = (int *)malloc(count * sizeof(int));
        batched->groups[g].coefficients = (double complex *)malloc(count * sizeof(double complex));
        batched->groups[g].Fsgn = (int *)malloc(count * sizeof(int));
        batched->groups[g].is_conj = (int *)malloc(count * sizeof(int));
        batched->groups[g].local_mask = (unsigned long int *)malloc(count * sizeof(unsigned long int));
        batched->groups[g].state1check = (unsigned long int *)malloc(count * sizeof(unsigned long int));
        batched->groups[g].bit1diff = (unsigned long int *)malloc(count * sizeof(unsigned long int));

        if (batched->groups[g].transfer_indices == NULL ||
            batched->groups[g].coefficients == NULL ||
            batched->groups[g].Fsgn == NULL ||
            batched->groups[g].is_conj == NULL ||
            batched->groups[g].local_mask == NULL ||
            batched->groups[g].state1check == NULL ||
            batched->groups[g].bit1diff == NULL) {
            // Cleanup on error
            for (int j = 0; j <= g; j++) {
                free(batched->groups[j].transfer_indices);
                free(batched->groups[j].coefficients);
                free(batched->groups[j].Fsgn);
                free(batched->groups[j].is_conj);
                free(batched->groups[j].local_mask);
                free(batched->groups[j].state1check);
                free(batched->groups[j].bit1diff);
            }
            free(batched->groups);
            free(origin_count);
            free(unique_origins);
            return -1;
        }
    }

    // Reset counts for second pass
    memset(origin_count, 0, nproc * sizeof(int));

    // Second pass: populate transfer information
    for (i = 0; i < (int)X->Def.EDNTransfer; i += 2) {
        int site1 = X->Def.EDGeneralTransfer[i][0];
        int site2 = X->Def.EDGeneralTransfer[i][2];
        int site1_local = (site1 + 1 <= (int)X->Def.Nsite);
        int site2_local = (site2 + 1 <= (int)X->Def.Nsite);

        // Only process MPIsingle transfers
        if (site1_local == site2_local) continue;

        // Determine local and inter-process sites
        int local_site, inter_site, trans_idx;
        double complex trans_coeff;

        if (site2_local) {
            // site1 is inter-process (use entry i+1 which has sites swapped)
            local_site = site2;
            inter_site = site1;
            trans_idx = i + 1;
            trans_coeff = X->Def.EDParaGeneralTransfer[i + 1];
        } else {
            // site2 is inter-process
            local_site = site1;
            inter_site = site2;
            trans_idx = i;
            trans_coeff = X->Def.EDParaGeneralTransfer[i];
        }

        mask2 = (int)X->Def.Tpow[inter_site];
        origin = myrank ^ mask2;

        // Find the group for this origin
        for (g = 0; g < num_unique; g++) {
            if (batched->groups[g].origin == origin) {
                break;
            }
        }

        t = origin_count[origin]++;

        // Store the mask (same for all transfers in this group)
        batched->groups[g].mask = mask2;

        // Compute Fsgn for the inter-process part
        int bit2diff = mask2 - 1;
        int Fsgn;
        SgnBit((unsigned long int)(origin & bit2diff), &Fsgn);

        // Store transfer information
        batched->groups[g].transfer_indices[t] = trans_idx;
        batched->groups[g].local_mask[t] = X->Def.Tpow[local_site];
        batched->groups[g].bit1diff[t] = X->Def.Tpow[X->Def.Nsite - 1] * 2 -
                                          batched->groups[g].local_mask[t] * 2;

        // Determine state1check and coefficient based on origin state
        // Store Fsgn/is_conj for dynamic coefficient computation during time evolution.
        // NOTE: `coefficients[t]` is kept for compatibility/debug, but runtime paths
        // intentionally re-read EDParaGeneralTransfer so time-dependent coefficients
        // (e.g., in TE mode) are reflected without re-initializing batches.
        int state2 = origin & mask2;
        if (state2 == mask2) {
            // Inter-process site has electron, so local site should be empty
            batched->groups[g].state1check[t] = 0;
            batched->groups[g].Fsgn[t] = Fsgn;
            batched->groups[g].is_conj[t] = 0;
            batched->groups[g].coefficients[t] = -(double)Fsgn * trans_coeff;
        } else {
            // Inter-process site is empty, so local site should have electron
            batched->groups[g].state1check[t] = batched->groups[g].local_mask[t];
            batched->groups[g].Fsgn[t] = Fsgn;
            batched->groups[g].is_conj[t] = 1;
            batched->groups[g].coefficients[t] = -(double)Fsgn * conj(trans_coeff);
        }
    }

    free(origin_count);
    free(unique_origins);

    // Debug output: show batching statistics (parity with the Hubbard/Spin
    // models, which already print this; lets tests confirm the batched
    // Spinless MPIsingle path actually fired).
    {
        int total_transfers = 0, gg;
        for (gg = 0; gg < batched->num_groups; gg++) {
            total_transfers += batched->groups[gg].num_transfers;
        }
        if (myrank == 0) {
            fprintf(stdout, "  [MPI Batching] Spinless: %d MPIsingle transfers -> %d groups (%.1fx reduction)\n",
                    total_transfers, batched->num_groups,
                    batched->num_groups > 0 ? (double)total_transfers / batched->num_groups : 0.0);
        }
    }

    batched->is_initialized = 1;
    return 0;
}

/**
 * @brief Free memory allocated for batched transfers
 */
void FinalizeMPIBatchedTransfers(MPIBatchedTransfers *batched) {
    int g;

    if (batched == NULL || !batched->is_initialized) return;

    for (g = 0; g < batched->num_groups; g++) {
        free(batched->groups[g].transfer_indices);
        free(batched->groups[g].coefficients);
        free(batched->groups[g].Fsgn);
        free(batched->groups[g].is_conj);
        free(batched->groups[g].local_mask);
        free(batched->groups[g].state1check);
        free(batched->groups[g].bit1diff);
    }
    free(batched->groups);

    batched->groups = NULL;
    batched->num_groups = 0;
    batched->is_initialized = 0;
}

/**
 * @brief Perform batched MPI hopping for SpinlessFermionGC MPIsingle mode
 *
 * Computes contribution from inter-process hopping terms to H|psi> (M_MLTPLY)
 * or correlation functions (M_CORR). All transfers in `group` share the same
 * MPI partner, so a single MPI_Sendrecv suffices.
 *
 * Parallel safety:
 * - M_MLTPLY/M_CALCSPEC: t-outer loop ensures each OMP region processes one
 *   transfer at a time; within that region, distinct j values map to distinct
 *   ioff (since ioff = jreal ^ mask1 is bijective for fixed mask1), so no race.
 * - M_CORR: No tmp_v0 write; dam_pr is OMP reduction. Loop fusion (j-outer)
 *   is safe and reduces OMP overhead.
 *
 * Index convention:
 * - j: 1-based buffer index (matches tmp_v1/v1buf layout)
 * - jreal = j-1: 0-based bit-state id
 * - ioff = (jreal ^ mask1) + 1: 1-based destination index
 *
 * Arithmetic equivalence: loop reorder changes only accumulation order
 * (floating-point roundoff at ~1e-15 level), not algebraic contributions.
 *
 * @param group  Transfer group (all share same MPI origin)
 * @param X      Bind structure
 * @param tmp_v0 Output vector (updated in M_MLTPLY/M_CALCSPEC only)
 * @param tmp_v1 Input vector
 * @return dam_pr = sum of <tmp_v1|H|tmp_v1> contributions (always computed)
 */
double complex X_child_GC_general_hopp_SpinlessFermion_MPIsingle_batched(
    MPITransferGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
) {
    unsigned long int idim_max_buf, j;
    int ierr, t, Fsgn;
    int num_transfers;
    MPI_Status statusMPI;
    double complex dam_pr = 0.0;

    if (group == NULL || group->num_transfers == 0) return 0.0;

    num_transfers = group->num_transfers;

    // Single MPI exchange for all transfers in this group
    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    // Process all transfers using the received data
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
        // t-outer loop: each transfer processed sequentially to avoid tmp_v0 race
        for (t = 0; t < num_transfers; t++) {
            unsigned long int mask1 = group->local_mask[t];
            unsigned long int state1check = group->state1check[t];
            unsigned long int bit1diff = group->bit1diff[t];
            int trans_idx = group->transfer_indices[t];
            int trans_Fsgn = group->Fsgn[t];
            double complex trans_coeff = X->Def.EDParaGeneralTransfer[trans_idx];
            double complex trans = group->is_conj[t] ?
                -(double)trans_Fsgn * conj(trans_coeff) :
                -(double)trans_Fsgn * trans_coeff;

            // Drop numerically negligible transfer amplitude (legacy threshold).
            if (cabs(trans) < 1e-15) continue;

#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, Fsgn) \
    firstprivate(idim_max_buf, trans, mask1, state1check, bit1diff) \
    shared(v1buf, tmp_v1, tmp_v0)
            for (j = 1; j <= idim_max_buf; j++) {
                // j: 1-based buffer index, jreal: 0-based bit-state id.
                unsigned long int jreal = j - 1;
                unsigned long int state1 = jreal & mask1;

                if (state1 == state1check) {
                    SgnBit((unsigned long int)(jreal & bit1diff), &Fsgn);
                    // `jreal ^ mask1` is 0-based destination state; +1 for buffer index.
                    unsigned long int ioff = (jreal ^ mask1) + 1;
                    double complex dmv = (double)Fsgn * trans * v1buf[j];
                    tmp_v0[ioff] += dmv;
                    dam_pr += conj(tmp_v1[ioff]) * dmv;
                }
            }
        }
    } else {
        // M_CORR: loop fusion safe (only dam_pr reduction, no tmp_v0 write)
#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, t, Fsgn) \
    firstprivate(idim_max_buf, num_transfers) \
    shared(v1buf, tmp_v1, group, X)
        for (j = 1; j <= idim_max_buf; j++) {
            unsigned long int jreal = j - 1;

            for (t = 0; t < num_transfers; t++) {
                unsigned long int mask1 = group->local_mask[t];
                unsigned long int state1check = group->state1check[t];
                unsigned long int state1 = jreal & mask1;

                if (state1 == state1check) {
                    unsigned long int bit1diff = group->bit1diff[t];
                    int trans_idx = group->transfer_indices[t];
                    int trans_Fsgn = group->Fsgn[t];
                    double complex trans_coeff = X->Def.EDParaGeneralTransfer[trans_idx];
                    double complex trans = group->is_conj[t] ?
                        -(double)trans_Fsgn * conj(trans_coeff) :
                        -(double)trans_Fsgn * trans_coeff;

                    SgnBit(jreal & bit1diff, &Fsgn);
                    unsigned long int ioff = (jreal ^ mask1) + 1;
                    double complex dmv = (double)Fsgn * trans * v1buf[j];
                    dam_pr += conj(tmp_v1[ioff]) * dmv;
                }
            }
        }
    }

    return dam_pr;
}

/**
 * @brief Initialize batched transfers for SpinlessFermion (canonical) MPIsingle mode
 *
 * Similar to GC version but requires list_1/list_1buf exchange.
 */
int InitializeMPIBatchedTransfers_SpinlessFermion(
    struct BindStruct *X,
    MPIBatchedTransfers *batched
) {
    // The initialization is identical to the GC version
    // The difference is in the processing function which uses list_1buf
    return InitializeMPIBatchedTransfers_SpinlessFermionGC(X, batched);
}

/**
 * @brief Perform batched MPI hopping for SpinlessFermion (canonical) MPIsingle mode
 *
 * Canonical version: Hilbert space is restricted by particle number, so basis
 * states are enumerated in list_1. GetOffComp maps bit-state to list index.
 *
 * Parallel safety:
 * - M_MLTPLY: t-outer ensures single transfer per OMP region; GetOffComp returns
 *   unique ioff per valid jreal, so no race on tmp_v0.
 * - M_CORR: No tmp_v0 write; loop fusion (j-outer) is safe.
 *
 * Index convention:
 * - j: 1-based index into list_1buf (received basis states)
 * - jreal = list_1buf[j]: actual bit-state id
 * - ioff: 1-based index into tmp_v0/tmp_v1 (from GetOffComp)
 *
 * @param group  Transfer group (all share same MPI origin)
 * @param X      Bind structure
 * @param tmp_v0 Output vector (updated in M_MLTPLY only)
 * @param tmp_v1 Input vector
 * @return dam_pr = sum of <tmp_v1|H|tmp_v1> contributions
 */
double complex X_child_general_hopp_Spinless_MPIsingle_batched(
    MPITransferGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
) {
    unsigned long int idim_max_buf, j, ioff, jreal, state1;
    int ierr, t, Fsgn;
    int num_transfers;
    MPI_Status statusMPI;
    double complex dam_pr = 0.0;

    if (group == NULL || group->num_transfers == 0) return 0.0;

    num_transfers = group->num_transfers;

    // Single MPI exchange for all transfers in this group
    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    // Process all transfers using the received data
    if (X->Large.mode == M_MLTPLY) {
        // t-outer loop: each transfer processed sequentially to avoid tmp_v0 race
        for (t = 0; t < num_transfers; t++) {
            unsigned long int mask1 = group->local_mask[t];
            unsigned long int state1check = group->state1check[t];
            unsigned long int bit1diff = group->bit1diff[t];
            int trans_idx = group->transfer_indices[t];
            int trans_Fsgn = group->Fsgn[t];
            double complex trans_coeff = X->Def.EDParaGeneralTransfer[trans_idx];
            double complex trans = group->is_conj[t] ?
                -(double)trans_Fsgn * conj(trans_coeff) :
                -(double)trans_Fsgn * trans_coeff;

            // Drop numerically negligible transfer amplitude (legacy threshold).
            if (cabs(trans) < 1e-15) continue;

#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, Fsgn, ioff, jreal, state1) \
    firstprivate(idim_max_buf, trans, mask1, state1check, bit1diff) \
    shared(list_1buf, list_2_1, list_2_2, v1buf, tmp_v1, tmp_v0, X)
            for (j = 1; j <= idim_max_buf; j++) {
                // j is 1-based basis index in exchanged list; jreal is bit-state id.
                jreal = list_1buf[j];
                state1 = jreal & mask1;

                if (state1 == state1check) {
                    SgnBit((unsigned long int)(jreal & bit1diff), &Fsgn);
                    if (GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
                                   X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == FALSE) {
                        continue;
                    }
                    double complex dmv = (double)Fsgn * trans * v1buf[j];
                    tmp_v0[ioff] += dmv;
                    dam_pr += conj(tmp_v1[ioff]) * dmv;
                }
            }
        }
    } else {
        // M_CORR: loop fusion safe (only dam_pr reduction, no tmp_v0 write)
#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, t, Fsgn, ioff, jreal, state1) \
    firstprivate(idim_max_buf, num_transfers) \
    shared(list_1buf, list_2_1, list_2_2, v1buf, tmp_v1, group, X)
        for (j = 1; j <= idim_max_buf; j++) {
            jreal = list_1buf[j];

            for (t = 0; t < num_transfers; t++) {
                unsigned long int mask1 = group->local_mask[t];
                unsigned long int state1check = group->state1check[t];
                state1 = jreal & mask1;

                if (state1 == state1check) {
                    unsigned long int bit1diff = group->bit1diff[t];
                    int trans_idx = group->transfer_indices[t];
                    int trans_Fsgn = group->Fsgn[t];
                    double complex trans_coeff = X->Def.EDParaGeneralTransfer[trans_idx];
                    double complex trans = group->is_conj[t] ?
                        -(double)trans_Fsgn * conj(trans_coeff) :
                        -(double)trans_Fsgn * trans_coeff;

                    SgnBit(jreal & bit1diff, &Fsgn);
                    if (GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
                                   X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == FALSE) {
                        continue;
                    }
                    double complex dmv = (double)Fsgn * trans * v1buf[j];
                    dam_pr += conj(tmp_v1[ioff]) * dmv;
                }
            }
        }
    }

    return dam_pr;
}

/*******************************************************************************
 * HubbardGC batched MPI functions
 ******************************************************************************/

/**
 * @brief Initialize batched transfers for HubbardGC MPIsingle mode
 *
 * Similar to SpinlessFermionGC but with spin handling (mask = Tpow[2*site + spin]).
 */
int InitializeMPIBatchedTransfers_HubbardGC(
    struct BindStruct *X,
    MPIBatchedTransfers *batched
) {
    int i, g, t;
    int mask2, origin;
    int *origin_count;
    int *unique_origins;
    int num_unique = 0;
    int max_transfers;

    if (batched == NULL) return -1;

    batched->is_initialized = 0;
    batched->num_groups = 0;
    batched->groups = NULL;

    // Count maximum possible transfers for MPIsingle
    max_transfers = 0;
    for (i = 0; i < (int)X->Def.EDNTransfer; i += 2) {
        int site1 = X->Def.EDGeneralTransfer[i][0];
        int site2 = X->Def.EDGeneralTransfer[i][2];
        int site1_local = (site1 + 1 <= (int)X->Def.Nsite);
        int site2_local = (site2 + 1 <= (int)X->Def.Nsite);

        if (site1_local != site2_local) {
            max_transfers += 2;
        }
    }

    if (max_transfers == 0) {
        batched->is_initialized = 1;
        return 0;
    }

    // Allocate temporary arrays for counting
    origin_count = (int *)calloc(nproc, sizeof(int));
    unique_origins = (int *)malloc(nproc * sizeof(int));
    if (origin_count == NULL || unique_origins == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    // First pass: count transfers per origin and identify unique origins
    for (i = 0; i < (int)X->Def.EDNTransfer; i += 2) {
        int site1 = X->Def.EDGeneralTransfer[i][0];
        int spin1 = X->Def.EDGeneralTransfer[i][1];
        int site2 = X->Def.EDGeneralTransfer[i][2];
        int spin2 = X->Def.EDGeneralTransfer[i][3];
        int site1_local = (site1 + 1 <= (int)X->Def.Nsite);
        int site2_local = (site2 + 1 <= (int)X->Def.Nsite);

        if (site1_local == site2_local) continue;

        // Determine the inter-process site
        int inter_site, inter_spin;
        if (site1_local) {
            inter_site = site2;
            inter_spin = spin2;
        } else {
            inter_site = site1;
            inter_spin = spin1;
        }

        mask2 = (int)X->Def.Tpow[2 * inter_site + inter_spin];
        // Partner rank differs by the PE-bit of (site,spin); XOR selects the peer rank.
        origin = myrank ^ mask2;

        if (origin_count[origin] == 0) {
            unique_origins[num_unique++] = origin;
        }
        origin_count[origin]++;
    }

    if (num_unique == 0) {
        free(origin_count);
        free(unique_origins);
        batched->is_initialized = 1;
        return 0;
    }

    // Allocate groups
    batched->num_groups = num_unique;
    batched->groups = (MPITransferGroup *)calloc(num_unique, sizeof(MPITransferGroup));
    if (batched->groups == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    // Initialize each group
    for (g = 0; g < num_unique; g++) {
        origin = unique_origins[g];
        int count = origin_count[origin];

        batched->groups[g].origin = origin;
        batched->groups[g].num_transfers = count;
        batched->groups[g].mask = 0;
        batched->groups[g].transfer_indices = (int *)malloc(count * sizeof(int));
        batched->groups[g].coefficients = (double complex *)malloc(count * sizeof(double complex));
        batched->groups[g].Fsgn = (int *)malloc(count * sizeof(int));
        batched->groups[g].is_conj = (int *)malloc(count * sizeof(int));
        batched->groups[g].local_mask = (unsigned long int *)malloc(count * sizeof(unsigned long int));
        batched->groups[g].state1check = (unsigned long int *)malloc(count * sizeof(unsigned long int));
        batched->groups[g].bit1diff = (unsigned long int *)malloc(count * sizeof(unsigned long int));

        if (batched->groups[g].transfer_indices == NULL ||
            batched->groups[g].coefficients == NULL ||
            batched->groups[g].Fsgn == NULL ||
            batched->groups[g].is_conj == NULL ||
            batched->groups[g].local_mask == NULL ||
            batched->groups[g].state1check == NULL ||
            batched->groups[g].bit1diff == NULL) {
            for (int j = 0; j <= g; j++) {
                free(batched->groups[j].transfer_indices);
                free(batched->groups[j].coefficients);
                free(batched->groups[j].Fsgn);
                free(batched->groups[j].is_conj);
                free(batched->groups[j].local_mask);
                free(batched->groups[j].state1check);
                free(batched->groups[j].bit1diff);
            }
            free(batched->groups);
            free(origin_count);
            free(unique_origins);
            return -1;
        }
    }

    // Reset counts for second pass
    memset(origin_count, 0, nproc * sizeof(int));

    // Second pass: populate transfer information
    for (i = 0; i < (int)X->Def.EDNTransfer; i += 2) {
        int site1 = X->Def.EDGeneralTransfer[i][0];
        int spin1 = X->Def.EDGeneralTransfer[i][1];
        int site2 = X->Def.EDGeneralTransfer[i][2];
        int spin2 = X->Def.EDGeneralTransfer[i][3];
        int site1_local = (site1 + 1 <= (int)X->Def.Nsite);
        int site2_local = (site2 + 1 <= (int)X->Def.Nsite);

        if (site1_local == site2_local) continue;

        int local_site, local_spin, inter_site, inter_spin, trans_idx;
        double complex trans_coeff;

        if (site2_local) {
            local_site = site2;
            local_spin = spin2;
            inter_site = site1;
            inter_spin = spin1;
            trans_idx = i + 1;
            trans_coeff = X->Def.EDParaGeneralTransfer[i + 1];
        } else {
            local_site = site1;
            local_spin = spin1;
            inter_site = site2;
            inter_spin = spin2;
            trans_idx = i;
            trans_coeff = X->Def.EDParaGeneralTransfer[i];
        }

        mask2 = (int)X->Def.Tpow[2 * inter_site + inter_spin];
        origin = myrank ^ mask2;

        for (g = 0; g < num_unique; g++) {
            if (batched->groups[g].origin == origin) {
                break;
            }
        }

        t = origin_count[origin]++;

        batched->groups[g].mask = mask2;

        int bit2diff = mask2 - 1;
        int Fsgn;
        SgnBit((unsigned long int)(origin & bit2diff), &Fsgn);

        batched->groups[g].transfer_indices[t] = trans_idx;
        batched->groups[g].local_mask[t] = X->Def.Tpow[2 * local_site + local_spin];
        batched->groups[g].bit1diff[t] = X->Def.Tpow[2 * X->Def.Nsite - 1] * 2 -
                                          batched->groups[g].local_mask[t] * 2;

        int state2 = origin & mask2;
        if (state2 == mask2) {
            batched->groups[g].state1check[t] = 0;
            batched->groups[g].Fsgn[t] = Fsgn;
            batched->groups[g].is_conj[t] = 0;
            batched->groups[g].coefficients[t] = -(double)Fsgn * trans_coeff;
        } else {
            batched->groups[g].state1check[t] = batched->groups[g].local_mask[t];
            batched->groups[g].Fsgn[t] = Fsgn;
            batched->groups[g].is_conj[t] = 1;
            batched->groups[g].coefficients[t] = -(double)Fsgn * conj(trans_coeff);
        }
    }

    free(origin_count);
    free(unique_origins);

    // Debug output: show batching statistics
    {
        int total_transfers = 0;
        for (g = 0; g < batched->num_groups; g++) {
            total_transfers += batched->groups[g].num_transfers;
        }
        if (myrank == 0) {
            fprintf(stdout, "  [MPI Batching] HubbardGC: %d MPIsingle transfers -> %d groups (%.1fx reduction)\n",
                    total_transfers, batched->num_groups,
                    batched->num_groups > 0 ? (double)total_transfers / batched->num_groups : 0.0);
        }
    }

    batched->is_initialized = 1;
    return 0;
}

/**
 * @brief Perform batched MPI hopping for HubbardGC MPIsingle mode
 *
 * HubbardGC variant: Hilbert space dimension is 4^Nsite (up/down per site).
 * Uses same parallel-safety strategy as SpinlessFermionGC version.
 *
 * Parallel safety:
 * - M_MLTPLY/M_CALCSPEC: t-outer ensures no concurrent tmp_v0[ioff] writes.
 * - M_CORR: dam_pr reduction only; loop fusion (j-outer) safe.
 *
 * Index convention:
 * - j: 0-based bit-state index (GC direct indexing)
 * - ioff = j ^ mask1: 0-based destination; +1 for 1-based buffer access
 *
 * @return dam_pr = sum of <tmp_v1|H|tmp_v1> contributions
 */
double complex X_child_GC_general_hopp_MPIsingle_batched(
    MPITransferGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
) {
    unsigned long int idim_max_buf, j;
    int ierr, t, Fsgn;
    int num_transfers;
    MPI_Status statusMPI;
    double complex dam_pr = 0.0;

    if (group == NULL || group->num_transfers == 0) return 0.0;

    num_transfers = group->num_transfers;

    // Single MPI exchange for all transfers in this group
    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    // Process all transfers using the received data
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
        // t-outer loop: each transfer processed sequentially to avoid tmp_v0 race
        for (t = 0; t < num_transfers; t++) {
            unsigned long int mask1 = group->local_mask[t];
            unsigned long int state1check = group->state1check[t];
            unsigned long int bit1diff = group->bit1diff[t];
            int trans_idx = group->transfer_indices[t];
            int trans_Fsgn = group->Fsgn[t];
            double complex trans_coeff = X->Def.EDParaGeneralTransfer[trans_idx];
            double complex trans = group->is_conj[t] ?
                -(double)trans_Fsgn * conj(trans_coeff) :
                -(double)trans_Fsgn * trans_coeff;

            // Drop numerically negligible transfer amplitude (legacy threshold).
            if (cabs(trans) < 1e-15) continue;

#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, Fsgn) \
    firstprivate(idim_max_buf, trans, mask1, state1check, bit1diff) \
    shared(v1buf, tmp_v1, tmp_v0)
            // Here j is 0-based state id; v1buf/tmp_v* access remains 1-based (+1).
            for (j = 0; j < idim_max_buf; j++) {
                unsigned long int state1 = j & mask1;

                if (state1 == state1check) {
                    SgnBit(j & bit1diff, &Fsgn);
                    unsigned long int ioff = j ^ mask1;
                    double complex dmv = (double)Fsgn * trans * v1buf[j + 1];
                    tmp_v0[ioff + 1] += dmv;
                    dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
                }
            }
        }
    } else {
        // M_CORR: loop fusion safe (only dam_pr reduction, no tmp_v0 write)
#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, t, Fsgn) \
    firstprivate(idim_max_buf, num_transfers) \
    shared(v1buf, tmp_v1, group, X)
        for (j = 0; j < idim_max_buf; j++) {
            for (t = 0; t < num_transfers; t++) {
                unsigned long int mask1 = group->local_mask[t];
                unsigned long int state1check = group->state1check[t];
                unsigned long int state1 = j & mask1;

                if (state1 == state1check) {
                    unsigned long int bit1diff = group->bit1diff[t];
                    int trans_idx = group->transfer_indices[t];
                    int trans_Fsgn = group->Fsgn[t];
                    double complex trans_coeff = X->Def.EDParaGeneralTransfer[trans_idx];
                    double complex trans = group->is_conj[t] ?
                        -(double)trans_Fsgn * conj(trans_coeff) :
                        -(double)trans_Fsgn * trans_coeff;

                    SgnBit(j & bit1diff, &Fsgn);
                    unsigned long int ioff = j ^ mask1;
                    double complex dmv = (double)Fsgn * trans * v1buf[j + 1];
                    dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
                }
            }
        }
    }

    return dam_pr;
}

/*******************************************************************************
 * Hubbard (canonical) batched MPI functions
 ******************************************************************************/

/**
 * @brief Initialize batched transfers for Hubbard (canonical) MPIsingle mode
 */
int InitializeMPIBatchedTransfers_Hubbard(
    struct BindStruct *X,
    MPIBatchedTransfers *batched
) {
    // The initialization is identical to the HubbardGC version
    return InitializeMPIBatchedTransfers_HubbardGC(X, batched);
}

/**
 * @brief Perform batched MPI hopping for Hubbard (canonical) MPIsingle mode
 *
 * Canonical Hubbard: Hilbert space restricted by (Nup, Ndown). Uses list_1buf
 * for basis state mapping and GetOffComp for index lookup.
 *
 * Parallel safety:
 * - M_MLTPLY/M_CALCSPEC: t-outer ensures single transfer per OMP region;
 *   GetOffComp returns unique ioff per valid state.
 * - M_CORR: dam_pr reduction only; loop fusion safe.
 *
 * Index convention:
 * - j: 1-based index into list_1buf
 * - jreal = list_1buf[j]: actual bit-state
 * - ioff: 1-based tmp_v0/tmp_v1 index from GetOffComp
 *
 * @return dam_pr = sum of <tmp_v1|H|tmp_v1> contributions
 */
double complex X_child_general_hopp_MPIsingle_batched(
    MPITransferGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
) {
    unsigned long int idim_max_buf, j, ioff, jreal, state1;
    int ierr, t, Fsgn;
    int num_transfers;
    MPI_Status statusMPI;
    double complex dam_pr = 0.0;

    if (group == NULL || group->num_transfers == 0) return 0.0;

    num_transfers = group->num_transfers;

    // Single MPI exchange for all transfers in this group
    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    // Process all transfers using the received data
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
        // t-outer loop: each transfer processed sequentially to avoid tmp_v0 race
        for (t = 0; t < num_transfers; t++) {
            unsigned long int mask1 = group->local_mask[t];
            unsigned long int state1check = group->state1check[t];
            unsigned long int bit1diff = group->bit1diff[t];
            int trans_idx = group->transfer_indices[t];
            int trans_Fsgn = group->Fsgn[t];
            double complex trans_coeff = X->Def.EDParaGeneralTransfer[trans_idx];
            double complex trans = group->is_conj[t] ?
                -(double)trans_Fsgn * conj(trans_coeff) :
                -(double)trans_Fsgn * trans_coeff;

            // Drop numerically negligible transfer amplitude (legacy threshold).
            if (cabs(trans) < 1e-15) continue;

#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, Fsgn, ioff, jreal, state1) \
    firstprivate(idim_max_buf, trans, mask1, state1check, bit1diff) \
    shared(list_1buf, list_2_1, list_2_2, v1buf, tmp_v1, tmp_v0, X)
            for (j = 1; j <= idim_max_buf; j++) {
                // j is 1-based basis index in exchanged list; jreal is bit-state id.
                jreal = list_1buf[j];
                state1 = jreal & mask1;

                if (state1 == state1check) {
                    SgnBit(jreal & bit1diff, &Fsgn);
                    if (GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
                                   X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == FALSE) {
                        continue;
                    }
                    double complex dmv = (double)Fsgn * trans * v1buf[j];
                    tmp_v0[ioff] += dmv;
                    dam_pr += conj(tmp_v1[ioff]) * dmv;
                }
            }
        }
    } else {
        // M_CORR: loop fusion safe (only dam_pr reduction, no tmp_v0 write)
#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, t, Fsgn, ioff, jreal, state1) \
    firstprivate(idim_max_buf, num_transfers) \
    shared(list_1buf, list_2_1, list_2_2, v1buf, tmp_v1, group, X)
        for (j = 1; j <= idim_max_buf; j++) {
            jreal = list_1buf[j];

            for (t = 0; t < num_transfers; t++) {
                unsigned long int mask1 = group->local_mask[t];
                unsigned long int state1check = group->state1check[t];
                state1 = jreal & mask1;

                if (state1 == state1check) {
                    unsigned long int bit1diff = group->bit1diff[t];
                    int trans_idx = group->transfer_indices[t];
                    int trans_Fsgn = group->Fsgn[t];
                    double complex trans_coeff = X->Def.EDParaGeneralTransfer[trans_idx];
                    double complex trans = group->is_conj[t] ?
                        -(double)trans_Fsgn * conj(trans_coeff) :
                        -(double)trans_Fsgn * trans_coeff;

                    SgnBit(jreal & bit1diff, &Fsgn);
                    if (GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
                                   X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == FALSE) {
                        continue;
                    }
                    double complex dmv = (double)Fsgn * trans * v1buf[j];
                    dam_pr += conj(tmp_v1[ioff]) * dmv;
                }
            }
        }
    }

    return dam_pr;
}

/*******************************************************************************
 * MPIdouble batched MPI functions (HubbardGC)
 ******************************************************************************/

/**
 * @brief Initialize batched MPIdouble transfers for HubbardGC
 */
int InitializeMPIBatchedDoubleTransfers_HubbardGC(
    struct BindStruct *X,
    MPIBatchedDoubleTransfers *batched
) {
    int i, g, t;
    int mask1, mask2, origin;
    int nproc, num_unique;
    int *origin_count;
    int *unique_origins;

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    batched->num_groups = 0;
    batched->groups = NULL;
    batched->is_initialized = 0;

    // First pass: count MPIdouble transfers per origin
    origin_count = (int *)calloc(nproc, sizeof(int));
    unique_origins = (int *)malloc(nproc * sizeof(int));
    if (origin_count == NULL || unique_origins == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    num_unique = 0;

    for (i = 0; i < (int)X->Def.EDNTransfer; i += 2) {
        int site1 = X->Def.EDGeneralTransfer[i][0];
        int spin1 = X->Def.EDGeneralTransfer[i][1];
        int site2 = X->Def.EDGeneralTransfer[i][2];
        int spin2 = X->Def.EDGeneralTransfer[i][3];

        // Check if both sites are inter-process (MPIdouble)
        if (site1 + 1 <= (int)X->Def.Nsite || site2 + 1 <= (int)X->Def.Nsite) {
            continue;  // Not MPIdouble
        }

        mask1 = (int)X->Def.Tpow[2 * site1 + spin1];
        mask2 = (int)X->Def.Tpow[2 * site2 + spin2];
        origin = myrank ^ (mask1 + mask2);

        // Don't filter by state here - all ranks must participate in communication
        // State validity is checked during processing

        if (origin_count[origin] == 0) {
            unique_origins[num_unique++] = origin;
        }
        origin_count[origin]++;
    }

    if (num_unique == 0) {
        if (myrank == 0) {
            fprintf(stdout, "  [MPI Batching] HubbardGC MPIdouble: No inter-process transfers\n");
        }
        free(origin_count);
        free(unique_origins);
        batched->is_initialized = 1;
        return 0;
    }

    // Allocate groups
    batched->num_groups = num_unique;
    batched->groups = (MPIDoubleTransferGroup *)calloc(num_unique, sizeof(MPIDoubleTransferGroup));
    if (batched->groups == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    // Initialize each group
    for (g = 0; g < num_unique; g++) {
        origin = unique_origins[g];
        int count = origin_count[origin];

        batched->groups[g].origin = origin;
        batched->groups[g].num_transfers = count;
        batched->groups[g].transfer_indices = (int *)malloc(count * sizeof(int));
        batched->groups[g].coefficients = (double complex *)malloc(count * sizeof(double complex));
        batched->groups[g].Fsgn = (int *)malloc(count * sizeof(int));
        batched->groups[g].is_conj = (int *)malloc(count * sizeof(int));

        if (batched->groups[g].transfer_indices == NULL ||
            batched->groups[g].coefficients == NULL ||
            batched->groups[g].Fsgn == NULL ||
            batched->groups[g].is_conj == NULL) {
            for (int j = 0; j <= g; j++) {
                free(batched->groups[j].transfer_indices);
                free(batched->groups[j].coefficients);
                free(batched->groups[j].Fsgn);
                free(batched->groups[j].is_conj);
            }
            free(batched->groups);
            free(origin_count);
            free(unique_origins);
            return -1;
        }
    }

    // Reset counts for second pass
    memset(origin_count, 0, nproc * sizeof(int));

    // Second pass: populate transfer information
    for (i = 0; i < (int)X->Def.EDNTransfer; i += 2) {
        int site1 = X->Def.EDGeneralTransfer[i][0];
        int spin1 = X->Def.EDGeneralTransfer[i][1];
        int site2 = X->Def.EDGeneralTransfer[i][2];
        int spin2 = X->Def.EDGeneralTransfer[i][3];
        double complex trans_coeff = X->Def.EDParaGeneralTransfer[i];

        if (site1 + 1 <= (int)X->Def.Nsite || site2 + 1 <= (int)X->Def.Nsite) {
            continue;
        }

        mask1 = (int)X->Def.Tpow[2 * site1 + spin1];
        mask2 = (int)X->Def.Tpow[2 * site2 + spin2];
        origin = myrank ^ (mask1 + mask2);

        int state1 = origin & mask1;
        int state2 = origin & mask2;

        int Fsgn;
        int bitdiff;
        if (mask2 > mask1) bitdiff = mask2 - mask1 * 2;
        else bitdiff = mask1 - mask2 * 2;
        SgnBit((unsigned long int)(origin & bitdiff), &Fsgn);

        int is_conj_flag;
        double complex trans;
        if (state1 == 0 && state2 == mask2) {
            trans = -(double)Fsgn * trans_coeff;
            is_conj_flag = 0;
        } else if (state1 == mask1 && state2 == 0) {
            trans = -(double)Fsgn * conj(trans_coeff);
            is_conj_flag = 1;
        } else {
            // Invalid state for this rank - set coefficient to 0
            // Communication still happens but no contribution
            trans = 0.0;
            is_conj_flag = -1;  // Special flag for invalid state
        }

        // Find the group for this origin
        for (g = 0; g < num_unique; g++) {
            if (batched->groups[g].origin == origin) break;
        }

        t = origin_count[origin]++;
        batched->groups[g].transfer_indices[t] = i;
        batched->groups[g].coefficients[t] = trans;
        batched->groups[g].Fsgn[t] = Fsgn;
        batched->groups[g].is_conj[t] = is_conj_flag;
    }

    free(origin_count);
    free(unique_origins);

    // Debug output
    {
        int total_transfers = 0;
        for (g = 0; g < batched->num_groups; g++) {
            total_transfers += batched->groups[g].num_transfers;
        }
        if (myrank == 0) {
            fprintf(stdout, "  [MPI Batching] HubbardGC MPIdouble: %d transfers -> %d groups (%.1fx reduction)\n",
                    total_transfers, batched->num_groups,
                    batched->num_groups > 0 ? (double)total_transfers / batched->num_groups : 0.0);
        }
    }

    batched->is_initialized = 1;
    return 0;
}

/**
 * @brief Free memory allocated for batched MPIdouble transfers
 */
void FinalizeMPIBatchedDoubleTransfers(MPIBatchedDoubleTransfers *batched) {
    int g;

    if (batched == NULL || !batched->is_initialized) return;

    for (g = 0; g < batched->num_groups; g++) {
        free(batched->groups[g].transfer_indices);
        free(batched->groups[g].coefficients);
        free(batched->groups[g].Fsgn);
        free(batched->groups[g].is_conj);
    }
    free(batched->groups);

    batched->groups = NULL;
    batched->num_groups = 0;
    batched->is_initialized = 0;
}

/**
 * @brief Perform batched MPI hopping for HubbardGC MPIdouble mode
 *
 * MPIdouble: Both hopping sites are inter-process, so destination index
 * equals source index (j -> j). This means each thread owns unique tmp_v0[j].
 *
 * Parallel safety:
 * - tmp_v0[j] += dmv: Each OMP thread owns distinct j, so no race.
 * - t-outer/j-inner chosen for efficiency: trans computed once per transfer
 *   (not O(idim_max_buf) times as in loop fusion).
 *
 * Threshold note:
 * - cabs(trans) < 1e-15: Skip numerically negligible amplitudes to avoid
 *   accumulating floating-point noise.
 *
 * @return dam_pr = sum of <tmp_v1|H|tmp_v1> contributions
 */
double complex X_child_GC_general_hopp_MPIdouble_batched(
    MPIDoubleTransferGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
) {
    unsigned long int idim_max_buf, j;
    int ierr, t;
    int num_transfers;
    MPI_Status statusMPI;
    double complex dam_pr = 0.0;

    if (group == NULL || group->num_transfers == 0) return 0.0;

    num_transfers = group->num_transfers;

    // Single MPI exchange for all transfers in this group
    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    // Process all transfers: t-outer to compute trans once per transfer
    // For MPIdouble, all j indices map 1:1 (no state filtering), so j-inner is efficient
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
        for (t = 0; t < num_transfers; t++) {
            int is_conj_flag = group->is_conj[t];
            if (is_conj_flag == -1) continue;  // Invalid state for this rank

            int trans_idx = group->transfer_indices[t];
            int trans_Fsgn = group->Fsgn[t];
            double complex trans_coeff = X->Def.EDParaGeneralTransfer[trans_idx];
            double complex trans = is_conj_flag ?
                -(double)trans_Fsgn * conj(trans_coeff) :
                -(double)trans_Fsgn * trans_coeff;

            // Drop numerically negligible transfer amplitude (legacy threshold).
            if (cabs(trans) < 1e-15) continue;

#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j) firstprivate(idim_max_buf, trans) shared(v1buf, tmp_v1, tmp_v0)
            for (j = 1; j <= idim_max_buf; j++) {
                double complex dmv = trans * v1buf[j];
                tmp_v0[j] += dmv;
                dam_pr += conj(tmp_v1[j]) * dmv;
            }
        }
    } else {
        for (t = 0; t < num_transfers; t++) {
            int is_conj_flag = group->is_conj[t];
            if (is_conj_flag == -1) continue;  // Invalid state for this rank

            int trans_idx = group->transfer_indices[t];
            int trans_Fsgn = group->Fsgn[t];
            double complex trans_coeff = X->Def.EDParaGeneralTransfer[trans_idx];
            double complex trans = is_conj_flag ?
                -(double)trans_Fsgn * conj(trans_coeff) :
                -(double)trans_Fsgn * trans_coeff;

            if (cabs(trans) < 1e-15) continue;

#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j) firstprivate(idim_max_buf, trans) shared(v1buf, tmp_v1)
            for (j = 1; j <= idim_max_buf; j++) {
                double complex dmv = trans * v1buf[j];
                dam_pr += conj(tmp_v1[j]) * dmv;
            }
        }
    }

    return dam_pr;
}

/*******************************************************************************
 * InterAll batched MPI functions (HubbardGC)
 ******************************************************************************/

/**
 * @brief Check if InterAll term requires inter-process communication
 *        and compute the origin (communication partner)
 *
 * This is similar to CheckBit_InterAllPE but simplified for origin computation.
 */
static int ComputeInterAllOrigin(
    int org_isite1, int org_ispin1,
    int org_isite2, int org_ispin2,
    int org_isite3, int org_ispin3,
    int org_isite4, int org_ispin4,
    struct BindStruct *X,
    unsigned long int *origin,
    int *is_hermite
) {
    unsigned long int tmp_org, tmp_off;
    unsigned long int tmp_ispin;
    int iflgBitExist = TRUE;
    int any_interPE = FALSE;

    // Check forward direction: c†_1 c_2 c†_3 c_4
    // Return contract:
    //   0  -> inter-process contribution exists; outputs {origin,is_hermite}
    //  -1  -> all sites are local (no MPI communication needed)
    //  -2  -> no valid occupancy path in either direction (zero contribution)
    tmp_org = (unsigned long int)myrank;
    tmp_off = 0;

    if (CheckPE(org_isite1, X) == TRUE) {
        any_interPE = TRUE;
        tmp_ispin = X->Def.Tpow[2 * org_isite1 + org_ispin1];
        if (CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
            iflgBitExist = FALSE;
        }
        tmp_org = tmp_off;
    }

    if (CheckPE(org_isite2, X) == TRUE) {
        any_interPE = TRUE;
        tmp_ispin = X->Def.Tpow[2 * org_isite2 + org_ispin2];
        if (CheckBit_Cis(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
            iflgBitExist = FALSE;
        }
        tmp_org = tmp_off;
    }

    if (CheckPE(org_isite3, X) == TRUE) {
        any_interPE = TRUE;
        tmp_ispin = X->Def.Tpow[2 * org_isite3 + org_ispin3];
        if (CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
            iflgBitExist = FALSE;
        }
        tmp_org = tmp_off;
    }

    if (CheckPE(org_isite4, X) == TRUE) {
        any_interPE = TRUE;
        tmp_ispin = X->Def.Tpow[2 * org_isite4 + org_ispin4];
        if (CheckBit_Cis(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
            iflgBitExist = FALSE;
        }
        tmp_org = tmp_off;
    }

    if (!any_interPE) {
        // All sites are local - no MPI needed
        return -1;
    }

    if (iflgBitExist == TRUE) {
        *origin = tmp_org;
        *is_hermite = 0;
        return 0;
    }

    // Try reverse direction (hermitian conjugate): c†_4 c_3 c†_2 c_1
    iflgBitExist = TRUE;
    tmp_org = (unsigned long int)myrank;

    if (CheckPE(org_isite4, X) == TRUE) {
        tmp_ispin = X->Def.Tpow[2 * org_isite4 + org_ispin4];
        if (CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
            iflgBitExist = FALSE;
        }
        tmp_org = tmp_off;
    }

    if (CheckPE(org_isite3, X) == TRUE) {
        tmp_ispin = X->Def.Tpow[2 * org_isite3 + org_ispin3];
        if (CheckBit_Cis(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
            iflgBitExist = FALSE;
        }
        tmp_org = tmp_off;
    }

    if (CheckPE(org_isite2, X) == TRUE) {
        tmp_ispin = X->Def.Tpow[2 * org_isite2 + org_ispin2];
        if (CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
            iflgBitExist = FALSE;
        }
        tmp_org = tmp_off;
    }

    if (CheckPE(org_isite1, X) == TRUE) {
        tmp_ispin = X->Def.Tpow[2 * org_isite1 + org_ispin1];
        if (CheckBit_Cis(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
            iflgBitExist = FALSE;
        }
        tmp_org = tmp_off;
    }

    if (iflgBitExist == TRUE) {
        *origin = tmp_org;
        *is_hermite = 1;
        return 0;
    }

    // Neither direction works - term doesn't contribute
    return -2;
}

/**
 * @brief Initialize batched InterAll for HubbardGC
 */
int InitializeMPIBatchedInterAll_HubbardGC(
    struct BindStruct *X,
    MPIBatchedInterAll *batched
) {
    int i, g;
    int nproc, num_unique;
    int *origin_count;
    unsigned long int *unique_origins;

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    batched->num_groups = 0;
    batched->groups = NULL;
    batched->is_initialized = 0;

    if (X->Def.NInterAll_OffDiagonal == 0) {
        if (myrank == 0) {
            fprintf(stdout, "  [MPI Batching] HubbardGC InterAll: No off-diagonal terms\n");
        }
        batched->is_initialized = 1;
        return 0;
    }

    // First pass: count terms per origin
    origin_count = (int *)calloc(nproc, sizeof(int));
    unique_origins = (unsigned long int *)malloc(nproc * sizeof(unsigned long int));
    if (origin_count == NULL || unique_origins == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    num_unique = 0;

    // Process pairs (i, i+1) for hermitian terms
    for (i = 0; i < (int)X->Def.NInterAll_OffDiagonal; i += 2) {
        int isite1 = X->Def.InterAll_OffDiagonal[i][0];
        int ispin1 = X->Def.InterAll_OffDiagonal[i][1];
        int isite2 = X->Def.InterAll_OffDiagonal[i][2];
        int ispin2 = X->Def.InterAll_OffDiagonal[i][3];
        int isite3 = X->Def.InterAll_OffDiagonal[i][4];
        int ispin3 = X->Def.InterAll_OffDiagonal[i][5];
        int isite4 = X->Def.InterAll_OffDiagonal[i][6];
        int ispin4 = X->Def.InterAll_OffDiagonal[i][7];

        unsigned long int origin;
        int is_hermite;
        int ret = ComputeInterAllOrigin(isite1, ispin1, isite2, ispin2,
                                        isite3, ispin3, isite4, ispin4,
                                        X, &origin, &is_hermite);

        if (ret != 0) continue;  // Skip local or non-contributing terms

        if ((int)origin == myrank) continue;  // Skip self-communication

        if (origin_count[origin] == 0) {
            unique_origins[num_unique++] = origin;
        }
        origin_count[origin]++;
    }

    if (num_unique == 0) {
        if (myrank == 0) {
            fprintf(stdout, "  [MPI Batching] HubbardGC InterAll: %d off-diagonal terms (all local)\n",
                    (int)X->Def.NInterAll_OffDiagonal / 2);
        }
        free(origin_count);
        free(unique_origins);
        batched->is_initialized = 1;
        return 0;
    }

    // Allocate groups
    batched->num_groups = num_unique;
    batched->groups = (MPIInterAllGroup *)calloc(num_unique, sizeof(MPIInterAllGroup));
    if (batched->groups == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    // Initialize each group
    for (g = 0; g < num_unique; g++) {
        unsigned long int origin = unique_origins[g];
        int count = origin_count[origin];

        batched->groups[g].origin = (int)origin;
        batched->groups[g].num_interall = count;
        batched->groups[g].interall_indices = (int *)malloc(count * sizeof(int));
        batched->groups[g].coefficients = (double complex *)malloc(count * sizeof(double complex));
        batched->groups[g].is_hermite = (int *)malloc(count * sizeof(int));
        batched->groups[g].tmp_isite1 = (unsigned long int *)malloc(count * sizeof(unsigned long int));
        batched->groups[g].tmp_isite2 = (unsigned long int *)malloc(count * sizeof(unsigned long int));
        batched->groups[g].tmp_isite3 = (unsigned long int *)malloc(count * sizeof(unsigned long int));
        batched->groups[g].tmp_isite4 = (unsigned long int *)malloc(count * sizeof(unsigned long int));
        batched->groups[g].isite1 = (unsigned long int *)malloc(count * sizeof(unsigned long int));
        batched->groups[g].isite2 = (unsigned long int *)malloc(count * sizeof(unsigned long int));
        batched->groups[g].isite3 = (unsigned long int *)malloc(count * sizeof(unsigned long int));
        batched->groups[g].isite4 = (unsigned long int *)malloc(count * sizeof(unsigned long int));

        if (batched->groups[g].interall_indices == NULL ||
            batched->groups[g].coefficients == NULL ||
            batched->groups[g].is_hermite == NULL ||
            batched->groups[g].tmp_isite1 == NULL ||
            batched->groups[g].tmp_isite2 == NULL ||
            batched->groups[g].tmp_isite3 == NULL ||
            batched->groups[g].tmp_isite4 == NULL ||
            batched->groups[g].isite1 == NULL ||
            batched->groups[g].isite2 == NULL ||
            batched->groups[g].isite3 == NULL ||
            batched->groups[g].isite4 == NULL) {
            // Cleanup on error
            for (int j = 0; j <= g; j++) {
                free(batched->groups[j].interall_indices);
                free(batched->groups[j].coefficients);
                free(batched->groups[j].is_hermite);
                free(batched->groups[j].tmp_isite1);
                free(batched->groups[j].tmp_isite2);
                free(batched->groups[j].tmp_isite3);
                free(batched->groups[j].tmp_isite4);
                free(batched->groups[j].isite1);
                free(batched->groups[j].isite2);
                free(batched->groups[j].isite3);
                free(batched->groups[j].isite4);
            }
            free(batched->groups);
            free(origin_count);
            free(unique_origins);
            return -1;
        }
    }

    // Reset counts for second pass
    memset(origin_count, 0, nproc * sizeof(int));

    // Second pass: populate InterAll information
    for (i = 0; i < (int)X->Def.NInterAll_OffDiagonal; i += 2) {
        int isite1 = X->Def.InterAll_OffDiagonal[i][0];
        int ispin1 = X->Def.InterAll_OffDiagonal[i][1];
        int isite2 = X->Def.InterAll_OffDiagonal[i][2];
        int ispin2 = X->Def.InterAll_OffDiagonal[i][3];
        int isite3 = X->Def.InterAll_OffDiagonal[i][4];
        int ispin3 = X->Def.InterAll_OffDiagonal[i][5];
        int isite4 = X->Def.InterAll_OffDiagonal[i][6];
        int ispin4 = X->Def.InterAll_OffDiagonal[i][7];
        double complex tmp_V = X->Def.ParaInterAll_OffDiagonal[i];

        unsigned long int origin;
        int is_hermite;
        int ret = ComputeInterAllOrigin(isite1, ispin1, isite2, ispin2,
                                        isite3, ispin3, isite4, ispin4,
                                        X, &origin, &is_hermite);

        if (ret != 0) continue;
        if ((int)origin == myrank) continue;

        // Find the group for this origin
        for (g = 0; g < num_unique; g++) {
            if (batched->groups[g].origin == (int)origin) break;
        }

        int t = origin_count[origin]++;

        batched->groups[g].interall_indices[t] = i;
        batched->groups[g].is_hermite[t] = is_hermite;

        if (is_hermite) {
            batched->groups[g].coefficients[t] = conj(tmp_V);
            batched->groups[g].tmp_isite1[t] = X->Def.OrgTpow[2 * isite4 + ispin4];
            batched->groups[g].tmp_isite2[t] = X->Def.OrgTpow[2 * isite3 + ispin3];
            batched->groups[g].tmp_isite3[t] = X->Def.OrgTpow[2 * isite2 + ispin2];
            batched->groups[g].tmp_isite4[t] = X->Def.OrgTpow[2 * isite1 + ispin1];
        } else {
            batched->groups[g].coefficients[t] = tmp_V;
            batched->groups[g].tmp_isite1[t] = X->Def.OrgTpow[2 * isite1 + ispin1];
            batched->groups[g].tmp_isite2[t] = X->Def.OrgTpow[2 * isite2 + ispin2];
            batched->groups[g].tmp_isite3[t] = X->Def.OrgTpow[2 * isite3 + ispin3];
            batched->groups[g].tmp_isite4[t] = X->Def.OrgTpow[2 * isite4 + ispin4];
        }

        batched->groups[g].isite1[t] = X->Def.Tpow[2 * isite1 + ispin1];
        batched->groups[g].isite2[t] = X->Def.Tpow[2 * isite2 + ispin2];
        batched->groups[g].isite3[t] = X->Def.Tpow[2 * isite3 + ispin3];
        batched->groups[g].isite4[t] = X->Def.Tpow[2 * isite4 + ispin4];
    }

    free(origin_count);
    free(unique_origins);

    // Debug output: show batching statistics
    {
        int total_interall = 0;
        for (g = 0; g < batched->num_groups; g++) {
            total_interall += batched->groups[g].num_interall;
        }
        if (myrank == 0) {
            fprintf(stdout, "  [MPI Batching] HubbardGC InterAll: %d MPI terms -> %d groups (%.1fx reduction)\n",
                    total_interall, batched->num_groups,
                    batched->num_groups > 0 ? (double)total_interall / batched->num_groups : 0.0);
        }
    }

    batched->is_initialized = 1;
    return 0;
}

/**
 * @brief Free memory allocated for batched InterAll
 */
void FinalizeMPIBatchedInterAll(MPIBatchedInterAll *batched) {
    int g;

    if (batched == NULL || !batched->is_initialized) return;

    for (g = 0; g < batched->num_groups; g++) {
        free(batched->groups[g].interall_indices);
        free(batched->groups[g].coefficients);
        free(batched->groups[g].is_hermite);
        free(batched->groups[g].tmp_isite1);
        free(batched->groups[g].tmp_isite2);
        free(batched->groups[g].tmp_isite3);
        free(batched->groups[g].tmp_isite4);
        free(batched->groups[g].isite1);
        free(batched->groups[g].isite2);
        free(batched->groups[g].isite3);
        free(batched->groups[g].isite4);
    }
    free(batched->groups);

    batched->groups = NULL;
    batched->num_groups = 0;
    batched->is_initialized = 0;
}

/**
 * @brief Perform batched MPI InterAll for HubbardGC
 */
double complex X_child_GC_InterAll_Hubbard_MPI_batched(
    MPIInterAllGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
) {
    unsigned long int idim_max_buf, j;
    int ierr, t, Fsgn;
    MPI_Status statusMPI;
    double complex dam_pr = 0.0;

    if (group == NULL || group->num_interall == 0) return 0.0;

    // Single MPI exchange for all InterAll terms in this group
    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    unsigned long int org_rankbit = X->Def.OrgTpow[2 * X->Def.Nsite] * group->origin;

    // Process all InterAll terms using the received data
    for (t = 0; t < group->num_interall; t++) {
        unsigned long int tmp_isite1 = group->tmp_isite1[t];
        unsigned long int tmp_isite2 = group->tmp_isite2[t];
        unsigned long int tmp_isite3 = group->tmp_isite3[t];
        unsigned long int tmp_isite4 = group->tmp_isite4[t];
        double complex tmp_V = group->coefficients[t];

        if (cabs(tmp_V) < 1e-15) continue;

        if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, Fsgn) \
    firstprivate(idim_max_buf, tmp_V, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit, X) \
    shared(v1buf, tmp_v1, tmp_v0)
            for (j = 1; j <= idim_max_buf; j++) {
                unsigned long int tmp_off;
                if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1,
                                   &Fsgn, X, (j - 1) + org_rankbit, &tmp_off) == TRUE) {
                    double complex dmv = tmp_V * v1buf[j] * Fsgn;
                    tmp_v0[tmp_off + 1] += dmv;
                    dam_pr += conj(tmp_v1[tmp_off + 1]) * dmv;
                }
            }
        } else {
#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, Fsgn) \
    firstprivate(idim_max_buf, tmp_V, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit, X) \
    shared(v1buf, tmp_v1, tmp_v0)
            for (j = 1; j <= idim_max_buf; j++) {
                unsigned long int tmp_off;
                if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1,
                                   &Fsgn, X, (j - 1) + org_rankbit, &tmp_off) == TRUE) {
                    double complex dmv = tmp_V * v1buf[j] * Fsgn;
                    dam_pr += conj(tmp_v1[tmp_off + 1]) * dmv;
                }
            }
        }
    }

    return dam_pr;
}

/*******************************************************************************
 * MPIdouble batched MPI functions (Hubbard canonical)
 ******************************************************************************/

/**
 * @brief Initialize batched MPIdouble transfers for Hubbard (canonical)
 *
 * Similar to HubbardGC but requires list_1/list_1buf exchange.
 */
int InitializeMPIBatchedDoubleTransfers_Hubbard(
    struct BindStruct *X,
    MPIBatchedDoubleTransfers *batched
) {
    int i, g, t;
    int mask1, mask2, origin;
    int nproc, num_unique;
    int *origin_count;
    int *unique_origins;

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    batched->num_groups = 0;
    batched->groups = NULL;
    batched->is_initialized = 0;

    // First pass: count MPIdouble transfers per origin
    origin_count = (int *)calloc(nproc, sizeof(int));
    unique_origins = (int *)malloc(nproc * sizeof(int));
    if (origin_count == NULL || unique_origins == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    num_unique = 0;

    for (i = 0; i < (int)X->Def.EDNTransfer; i += 2) {
        int site1 = X->Def.EDGeneralTransfer[i][0];
        int spin1 = X->Def.EDGeneralTransfer[i][1];
        int site2 = X->Def.EDGeneralTransfer[i][2];
        int spin2 = X->Def.EDGeneralTransfer[i][3];

        // Check if both sites are inter-process (MPIdouble)
        if (site1 + 1 <= (int)X->Def.Nsite || site2 + 1 <= (int)X->Def.Nsite) {
            continue;  // Not MPIdouble
        }

        mask1 = (int)X->Def.Tpow[2 * site1 + spin1];
        mask2 = (int)X->Def.Tpow[2 * site2 + spin2];
        origin = myrank ^ (mask1 + mask2);

        // Don't filter by state here - all ranks must participate in communication
        // State validity is checked during processing

        if (origin_count[origin] == 0) {
            unique_origins[num_unique++] = origin;
        }
        origin_count[origin]++;
    }

    if (num_unique == 0) {
        if (myrank == 0) {
            fprintf(stdout, "  [MPI Batching] Hubbard MPIdouble: No inter-process transfers\n");
        }
        free(origin_count);
        free(unique_origins);
        batched->is_initialized = 1;
        return 0;
    }

    // Allocate groups
    batched->num_groups = num_unique;
    batched->groups = (MPIDoubleTransferGroup *)calloc(num_unique, sizeof(MPIDoubleTransferGroup));
    if (batched->groups == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    // Initialize each group
    for (g = 0; g < num_unique; g++) {
        origin = unique_origins[g];
        int count = origin_count[origin];

        batched->groups[g].origin = origin;
        batched->groups[g].num_transfers = count;
        batched->groups[g].transfer_indices = (int *)malloc(count * sizeof(int));
        batched->groups[g].coefficients = (double complex *)malloc(count * sizeof(double complex));
        batched->groups[g].Fsgn = (int *)malloc(count * sizeof(int));
        batched->groups[g].is_conj = (int *)malloc(count * sizeof(int));

        if (batched->groups[g].transfer_indices == NULL ||
            batched->groups[g].coefficients == NULL ||
            batched->groups[g].Fsgn == NULL ||
            batched->groups[g].is_conj == NULL) {
            for (int j = 0; j <= g; j++) {
                free(batched->groups[j].transfer_indices);
                free(batched->groups[j].coefficients);
                free(batched->groups[j].Fsgn);
                free(batched->groups[j].is_conj);
            }
            free(batched->groups);
            free(origin_count);
            free(unique_origins);
            return -1;
        }
    }

    // Reset counts for second pass
    memset(origin_count, 0, nproc * sizeof(int));

    // Second pass: populate transfer information
    for (i = 0; i < (int)X->Def.EDNTransfer; i += 2) {
        int site1 = X->Def.EDGeneralTransfer[i][0];
        int spin1 = X->Def.EDGeneralTransfer[i][1];
        int site2 = X->Def.EDGeneralTransfer[i][2];
        int spin2 = X->Def.EDGeneralTransfer[i][3];
        double complex trans_coeff = X->Def.EDParaGeneralTransfer[i];

        if (site1 + 1 <= (int)X->Def.Nsite || site2 + 1 <= (int)X->Def.Nsite) {
            continue;
        }

        mask1 = (int)X->Def.Tpow[2 * site1 + spin1];
        mask2 = (int)X->Def.Tpow[2 * site2 + spin2];
        origin = myrank ^ (mask1 + mask2);

        int state1 = origin & mask1;
        int state2 = origin & mask2;

        int Fsgn;
        int bitdiff;
        if (mask2 > mask1) bitdiff = mask2 - mask1 * 2;
        else bitdiff = mask1 - mask2 * 2;
        SgnBit((unsigned long int)(origin & bitdiff), &Fsgn);

        int is_conj_flag;
        double complex trans;
        if (state1 == 0 && state2 == mask2) {
            trans = -(double)Fsgn * trans_coeff;
            is_conj_flag = 0;
        } else if (state1 == mask1 && state2 == 0) {
            trans = -(double)Fsgn * conj(trans_coeff);
            is_conj_flag = 1;
        } else {
            // Invalid state for this rank - set coefficient to 0
            // Communication still happens but no contribution
            trans = 0.0;
            is_conj_flag = -1;  // Special flag for invalid state
        }

        // Find the group for this origin
        for (g = 0; g < num_unique; g++) {
            if (batched->groups[g].origin == origin) break;
        }

        t = origin_count[origin]++;
        batched->groups[g].transfer_indices[t] = i;
        batched->groups[g].coefficients[t] = trans;
        batched->groups[g].Fsgn[t] = Fsgn;
        batched->groups[g].is_conj[t] = is_conj_flag;
    }

    free(origin_count);
    free(unique_origins);

    // Debug output
    {
        int total_transfers = 0;
        for (g = 0; g < batched->num_groups; g++) {
            total_transfers += batched->groups[g].num_transfers;
        }
        if (myrank == 0) {
            fprintf(stdout, "  [MPI Batching] Hubbard MPIdouble: %d transfers -> %d groups (%.1fx reduction)\n",
                    total_transfers, batched->num_groups,
                    batched->num_groups > 0 ? (double)total_transfers / batched->num_groups : 0.0);
        }
    }

    batched->is_initialized = 1;
    return 0;
}

/**
 * @brief Perform batched MPI hopping for Hubbard (canonical) MPIdouble mode
 *
 * Performs a single MPI_Sendrecv for all transfers in the group (including list_1),
 * then applies all transfer operations locally using the received data.
 */
double complex X_child_general_hopp_MPIdouble_batched(
    MPIDoubleTransferGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
) {
    unsigned long int idim_max_buf, j, ioff;
    int ierr, t;
    MPI_Status statusMPI;
    double complex dam_pr = 0.0;

    if (group == NULL || group->num_transfers == 0) return 0.0;

    // Single MPI exchange for all transfers in this group
    // For canonical Hubbard, we need to exchange list_1 as well
    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    // Process all transfers using the received data
    // Read coefficients from current EDParaGeneralTransfer to support time evolution
    // For canonical Hubbard, use list_1buf to find target index via GetOffComp
    for (t = 0; t < group->num_transfers; t++) {
        int is_conj_flag = group->is_conj[t];
        if (is_conj_flag == -1) continue;  // Invalid state for this rank

        int trans_idx = group->transfer_indices[t];
        int trans_Fsgn = group->Fsgn[t];
        double complex trans_coeff = X->Def.EDParaGeneralTransfer[trans_idx];
        double complex trans = is_conj_flag ?
            -(double)trans_Fsgn * conj(trans_coeff) :
            -(double)trans_Fsgn * trans_coeff;

        if (cabs(trans) < 1e-15) continue;

        if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, ioff) firstprivate(idim_max_buf, trans, X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
            for (j = 1; j <= idim_max_buf; j++) {
                GetOffComp(list_2_1, list_2_2, list_1buf[j],
                           X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
                double complex dmv = trans * v1buf[j];
                tmp_v0[ioff] += dmv;
                dam_pr += conj(tmp_v1[ioff]) * dmv;
            }
        } else {
#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, ioff) firstprivate(idim_max_buf, trans, X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
            for (j = 1; j <= idim_max_buf; j++) {
                GetOffComp(list_2_1, list_2_2, list_1buf[j],
                           X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
                double complex dmv = trans * v1buf[j];
                dam_pr += conj(tmp_v1[ioff]) * dmv;
            }
        }
    }

    return dam_pr;
}

/*******************************************************************************
 * Spin model batched MPI functions
 ******************************************************************************/

/**
 * @brief Initialize batched Exchange for SpinGC
 *
 * Groups Exchange and PairLift terms by their MPI communication partner.
 * For Spin 1/2, mask = Tpow[site] (not Tpow[2*site + spin] as in Hubbard).
 */
int InitializeMPIBatchedExchange_SpinGC(
    struct BindStruct *X,
    MPIBatchedSpinExchange *batched
) {
    int i, g, t;
    int mask, origin;
    int nproc, num_unique;
    int *origin_count;
    int *unique_origins;

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    batched->num_groups = 0;
    batched->groups = NULL;
    batched->is_initialized = 0;

    // Count total MPI terms (Exchange + PairLift with one site inter-process)
    int total_terms = 0;

    // Count Exchange terms requiring MPIsingle
    for (i = 0; i < (int)X->Def.NExchangeCoupling; i++) {
        int site0 = X->Def.ExchangeCoupling[i][0];
        int site1 = X->Def.ExchangeCoupling[i][1];
        // MPIsingle: exactly one site is inter-process
        if ((site0 + 1 > (int)X->Def.Nsite) != (site1 + 1 > (int)X->Def.Nsite)) {
            total_terms++;
        }
    }

    // Count PairLift terms requiring MPIsingle
    for (i = 0; i < (int)X->Def.NPairLiftCoupling; i++) {
        int site0 = X->Def.PairLiftCoupling[i][0];
        int site1 = X->Def.PairLiftCoupling[i][1];
        if ((site0 + 1 > (int)X->Def.Nsite) != (site1 + 1 > (int)X->Def.Nsite)) {
            total_terms++;
        }
    }

    if (total_terms == 0) {
        if (myrank == 0) {
            fprintf(stdout, "  [MPI Batching] SpinGC Exchange: No MPIsingle terms\n");
        }
        batched->is_initialized = 1;
        return 0;
    }

    // First pass: count terms per origin
    origin_count = (int *)calloc(nproc, sizeof(int));
    unique_origins = (int *)malloc(nproc * sizeof(int));
    if (origin_count == NULL || unique_origins == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    num_unique = 0;

    // Exchange terms
    for (i = 0; i < (int)X->Def.NExchangeCoupling; i++) {
        int site0 = X->Def.ExchangeCoupling[i][0];
        int site1 = X->Def.ExchangeCoupling[i][1];

        if (site1 + 1 > (int)X->Def.Nsite && site0 + 1 <= (int)X->Def.Nsite) {
            mask = (int)X->Def.Tpow[site1];
            origin = myrank ^ mask;
            if (origin_count[origin] == 0) {
                unique_origins[num_unique++] = origin;
            }
            origin_count[origin]++;
        }
        else if (site0 + 1 > (int)X->Def.Nsite && site1 + 1 <= (int)X->Def.Nsite) {
            mask = (int)X->Def.Tpow[site0];
            origin = myrank ^ mask;
            if (origin_count[origin] == 0) {
                unique_origins[num_unique++] = origin;
            }
            origin_count[origin]++;
        }
    }

    // PairLift terms
    for (i = 0; i < (int)X->Def.NPairLiftCoupling; i++) {
        int site0 = X->Def.PairLiftCoupling[i][0];
        int site1 = X->Def.PairLiftCoupling[i][1];

        if (site1 + 1 > (int)X->Def.Nsite && site0 + 1 <= (int)X->Def.Nsite) {
            mask = (int)X->Def.Tpow[site1];
            origin = myrank ^ mask;
            if (origin_count[origin] == 0) {
                unique_origins[num_unique++] = origin;
            }
            origin_count[origin]++;
        }
        else if (site0 + 1 > (int)X->Def.Nsite && site1 + 1 <= (int)X->Def.Nsite) {
            mask = (int)X->Def.Tpow[site0];
            origin = myrank ^ mask;
            if (origin_count[origin] == 0) {
                unique_origins[num_unique++] = origin;
            }
            origin_count[origin]++;
        }
    }

    // Allocate groups
    batched->num_groups = num_unique;
    batched->groups = (MPISpinExchangeGroup *)calloc(num_unique, sizeof(MPISpinExchangeGroup));
    if (batched->groups == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    // Initialize each group
    for (g = 0; g < num_unique; g++) {
        origin = unique_origins[g];
        int count = origin_count[origin];

        batched->groups[g].origin = origin;
        batched->groups[g].num_terms = count;
        batched->groups[g].term_indices = (int *)malloc(count * sizeof(int));
        batched->groups[g].coefficients = (double complex *)malloc(count * sizeof(double complex));
        batched->groups[g].org_isite1 = (int *)malloc(count * sizeof(int));
        batched->groups[g].org_ispin1 = (int *)malloc(count * sizeof(int));
        batched->groups[g].org_ispin2 = (int *)malloc(count * sizeof(int));
        batched->groups[g].state1check = (int *)malloc(count * sizeof(int));

        if (batched->groups[g].term_indices == NULL ||
            batched->groups[g].coefficients == NULL ||
            batched->groups[g].org_isite1 == NULL ||
            batched->groups[g].org_ispin1 == NULL ||
            batched->groups[g].org_ispin2 == NULL ||
            batched->groups[g].state1check == NULL) {
            // Cleanup on error
            for (int j = 0; j <= g; j++) {
                free(batched->groups[j].term_indices);
                free(batched->groups[j].coefficients);
                free(batched->groups[j].org_isite1);
                free(batched->groups[j].org_ispin1);
                free(batched->groups[j].org_ispin2);
                free(batched->groups[j].state1check);
            }
            free(batched->groups);
            free(origin_count);
            free(unique_origins);
            return -1;
        }
    }

    // Reset counts for second pass
    memset(origin_count, 0, nproc * sizeof(int));

    // Second pass: populate transfer information
    // Exchange terms (sigma1=0->1, sigma2=1->0 for exchange)
    for (i = 0; i < (int)X->Def.NExchangeCoupling; i++) {
        int site0 = X->Def.ExchangeCoupling[i][0];
        int site1 = X->Def.ExchangeCoupling[i][1];
        double complex J = X->Def.ParaExchangeCoupling[i];

        int local_site, inter_site;
        int is_reverse = 0;

        if (site1 + 1 > (int)X->Def.Nsite && site0 + 1 <= (int)X->Def.Nsite) {
            local_site = site0;
            inter_site = site1;
        }
        else if (site0 + 1 > (int)X->Def.Nsite && site1 + 1 <= (int)X->Def.Nsite) {
            local_site = site1;
            inter_site = site0;
            is_reverse = 1;
            J = conj(J);
        }
        else {
            continue;  // Not MPIsingle
        }

        mask = (int)X->Def.Tpow[inter_site];
        origin = myrank ^ mask;
        int state2 = (origin & mask) / mask;

        // Exchange: (0,1) <-> (1,0)
        // sigma1=0, sigma2=1 for local site
        // sigma3=1, sigma4=0 for inter site
        int org_ispin3 = 1;  // inter site before (up)
        int org_ispin4 = 0;  // inter site after (down)
        if (is_reverse) {
            org_ispin3 = 0;
            org_ispin4 = 1;
        }

        double complex Jint;
        int state1check;
        if (state2 == org_ispin4) {
            state1check = is_reverse ? 0 : 1;  // org_ispin2
            Jint = J;
        }
        else if (state2 == org_ispin3) {
            state1check = is_reverse ? 1 : 0;  // org_ispin1
            Jint = conj(J);
        }
        else {
            Jint = 0.0;
            state1check = 0;
        }

        // Find group
        for (g = 0; g < num_unique; g++) {
            if (batched->groups[g].origin == origin) break;
        }

        t = origin_count[origin]++;
        batched->groups[g].term_indices[t] = i;
        batched->groups[g].coefficients[t] = Jint;
        batched->groups[g].org_isite1[t] = local_site;
        batched->groups[g].org_ispin1[t] = is_reverse ? 1 : 0;
        batched->groups[g].org_ispin2[t] = is_reverse ? 0 : 1;
        batched->groups[g].state1check[t] = state1check;
    }

    // PairLift terms (similar pattern but sigma transitions are different)
    for (i = 0; i < (int)X->Def.NPairLiftCoupling; i++) {
        int site0 = X->Def.PairLiftCoupling[i][0];
        int site1 = X->Def.PairLiftCoupling[i][1];
        double complex J = X->Def.ParaPairLiftCoupling[i];

        int local_site, inter_site;
        int is_reverse = 0;

        if (site1 + 1 > (int)X->Def.Nsite && site0 + 1 <= (int)X->Def.Nsite) {
            local_site = site0;
            inter_site = site1;
        }
        else if (site0 + 1 > (int)X->Def.Nsite && site1 + 1 <= (int)X->Def.Nsite) {
            local_site = site1;
            inter_site = site0;
            is_reverse = 1;
            J = conj(J);
        }
        else {
            continue;
        }

        mask = (int)X->Def.Tpow[inter_site];
        origin = myrank ^ mask;
        int state2 = (origin & mask) / mask;

        // PairLift: (0,0) -> (1,1) or (1,1) -> (0,0)
        // For both sites: sigma1=0, sigma2=1 (up to down transition)
        int org_ispin3 = 0;  // inter site before
        int org_ispin4 = 1;  // inter site after

        double complex Jint;
        int state1check;
        if (state2 == org_ispin4) {
            state1check = 1;  // org_ispin2
            Jint = J;
        }
        else if (state2 == org_ispin3) {
            state1check = 0;  // org_ispin1
            Jint = conj(J);
        }
        else {
            Jint = 0.0;
            state1check = 0;
        }

        for (g = 0; g < num_unique; g++) {
            if (batched->groups[g].origin == origin) break;
        }

        t = origin_count[origin]++;
        batched->groups[g].term_indices[t] = i + (int)X->Def.NExchangeCoupling;  // Mark as PairLift
        batched->groups[g].coefficients[t] = Jint;
        batched->groups[g].org_isite1[t] = local_site;
        batched->groups[g].org_ispin1[t] = 0;
        batched->groups[g].org_ispin2[t] = 1;
        batched->groups[g].state1check[t] = state1check;
    }

    free(origin_count);
    free(unique_origins);

    // Debug output
    if (myrank == 0) {
        fprintf(stdout, "  [MPI Batching] SpinGC Exchange: %d MPIsingle terms -> %d groups (%.1fx reduction)\n",
                total_terms, batched->num_groups,
                batched->num_groups > 0 ? (double)total_terms / batched->num_groups : 0.0);
    }

    batched->is_initialized = 1;
    return 0;
}

/**
 * @brief Free memory allocated for batched Spin Exchange
 */
void FinalizeMPIBatchedSpinExchange(MPIBatchedSpinExchange *batched) {
    int g;

    if (batched == NULL || !batched->is_initialized) return;

    for (g = 0; g < batched->num_groups; g++) {
        free(batched->groups[g].term_indices);
        free(batched->groups[g].coefficients);
        free(batched->groups[g].org_isite1);
        free(batched->groups[g].org_ispin1);
        free(batched->groups[g].org_ispin2);
        free(batched->groups[g].state1check);
    }
    free(batched->groups);

    batched->groups = NULL;
    batched->num_groups = 0;
    batched->is_initialized = 0;
}

/**
 * @brief Perform batched MPI Exchange for SpinGC
 */
double complex X_child_GC_CisAitCiuAiv_spin_MPIsingle_batched(
    MPISpinExchangeGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
) {
    unsigned long int idim_max_buf, j, ioff;
    int ierr, t;
    long int state1;
    MPI_Status statusMPI;
    double complex dam_pr = 0.0;

    if (group == NULL || group->num_terms == 0) return 0.0;

    // Single MPI exchange for all terms in this group
    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    // Process all terms using the received data
    for (t = 0; t < group->num_terms; t++) {
        double complex Jint = group->coefficients[t];
        int org_isite1 = group->org_isite1[t];
        int state1check = group->state1check[t];
        unsigned long int mask1 = X->Def.Tpow[org_isite1];

        if (cabs(Jint) < 1e-15) continue;

        if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, state1, ioff) \
    firstprivate(idim_max_buf, Jint, mask1, state1check, X) \
    shared(v1buf, tmp_v1, tmp_v0)
            for (j = 0; j < idim_max_buf; j++) {
                // child_SpinGC_CisAit logic inline
                long unsigned int ibit_tmp_1 = j & mask1;
                if ((ibit_tmp_1 == 0 && state1check == 0) ||
                    (ibit_tmp_1 != 0 && state1check == 1)) {
                    if (state1check == 0) {
                        ioff = j + mask1;
                    } else {
                        ioff = j - mask1;
                    }
                    double complex dmv = Jint * v1buf[j + 1];
                    tmp_v0[ioff + 1] += dmv;
                    dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
                }
            }
        } else {
#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, state1, ioff) \
    firstprivate(idim_max_buf, Jint, mask1, state1check, X) \
    shared(v1buf, tmp_v1, tmp_v0)
            for (j = 0; j < idim_max_buf; j++) {
                long unsigned int ibit_tmp_1 = j & mask1;
                if ((ibit_tmp_1 == 0 && state1check == 0) ||
                    (ibit_tmp_1 != 0 && state1check == 1)) {
                    if (state1check == 0) {
                        ioff = j + mask1;
                    } else {
                        ioff = j - mask1;
                    }
                    double complex dmv = Jint * v1buf[j + 1];
                    dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
                }
            }
        }
    }

    return dam_pr;
}

/**
 * @brief Initialize batched Exchange for Spin (canonical)
 */
int InitializeMPIBatchedExchange_Spin(
    struct BindStruct *X,
    MPIBatchedSpinExchange *batched
) {
    int i, g, t;
    int mask, origin;
    int nproc, num_unique;
    int *origin_count;
    int *unique_origins;

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    batched->num_groups = 0;
    batched->groups = NULL;
    batched->is_initialized = 0;

    // Count Exchange terms requiring MPIsingle
    int total_terms = 0;
    for (i = 0; i < (int)X->Def.NExchangeCoupling; i++) {
        int site0 = X->Def.ExchangeCoupling[i][0];
        int site1 = X->Def.ExchangeCoupling[i][1];
        if ((site0 + 1 > (int)X->Def.Nsite) != (site1 + 1 > (int)X->Def.Nsite)) {
            total_terms++;
        }
    }

    if (total_terms == 0) {
        if (myrank == 0) {
            fprintf(stdout, "  [MPI Batching] Spin Exchange: No MPIsingle terms\n");
        }
        batched->is_initialized = 1;
        return 0;
    }

    origin_count = (int *)calloc(nproc, sizeof(int));
    unique_origins = (int *)malloc(nproc * sizeof(int));
    if (origin_count == NULL || unique_origins == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    num_unique = 0;

    for (i = 0; i < (int)X->Def.NExchangeCoupling; i++) {
        int site0 = X->Def.ExchangeCoupling[i][0];
        int site1 = X->Def.ExchangeCoupling[i][1];

        if (site1 + 1 > (int)X->Def.Nsite && site0 + 1 <= (int)X->Def.Nsite) {
            mask = (int)X->Def.Tpow[site1];
            origin = myrank ^ mask;
            if (origin_count[origin] == 0) {
                unique_origins[num_unique++] = origin;
            }
            origin_count[origin]++;
        }
        else if (site0 + 1 > (int)X->Def.Nsite && site1 + 1 <= (int)X->Def.Nsite) {
            mask = (int)X->Def.Tpow[site0];
            origin = myrank ^ mask;
            if (origin_count[origin] == 0) {
                unique_origins[num_unique++] = origin;
            }
            origin_count[origin]++;
        }
    }

    // Allocate and populate groups (similar to SpinGC)
    batched->num_groups = num_unique;
    batched->groups = (MPISpinExchangeGroup *)calloc(num_unique, sizeof(MPISpinExchangeGroup));
    if (batched->groups == NULL) {
        free(origin_count);
        free(unique_origins);
        return -1;
    }

    for (g = 0; g < num_unique; g++) {
        origin = unique_origins[g];
        int count = origin_count[origin];

        batched->groups[g].origin = origin;
        batched->groups[g].num_terms = count;
        batched->groups[g].term_indices = (int *)malloc(count * sizeof(int));
        batched->groups[g].coefficients = (double complex *)malloc(count * sizeof(double complex));
        batched->groups[g].org_isite1 = (int *)malloc(count * sizeof(int));
        batched->groups[g].org_ispin1 = (int *)malloc(count * sizeof(int));
        batched->groups[g].org_ispin2 = (int *)malloc(count * sizeof(int));
        batched->groups[g].state1check = (int *)malloc(count * sizeof(int));

        if (batched->groups[g].term_indices == NULL ||
            batched->groups[g].coefficients == NULL ||
            batched->groups[g].org_isite1 == NULL ||
            batched->groups[g].org_ispin1 == NULL ||
            batched->groups[g].org_ispin2 == NULL ||
            batched->groups[g].state1check == NULL) {
            for (int j = 0; j <= g; j++) {
                free(batched->groups[j].term_indices);
                free(batched->groups[j].coefficients);
                free(batched->groups[j].org_isite1);
                free(batched->groups[j].org_ispin1);
                free(batched->groups[j].org_ispin2);
                free(batched->groups[j].state1check);
            }
            free(batched->groups);
            free(origin_count);
            free(unique_origins);
            return -1;
        }
    }

    memset(origin_count, 0, nproc * sizeof(int));

    // Populate Exchange terms
    for (i = 0; i < (int)X->Def.NExchangeCoupling; i++) {
        int site0 = X->Def.ExchangeCoupling[i][0];
        int site1 = X->Def.ExchangeCoupling[i][1];
        double complex J = X->Def.ParaExchangeCoupling[i];

        int local_site, inter_site;
        int is_reverse = 0;

        if (site1 + 1 > (int)X->Def.Nsite && site0 + 1 <= (int)X->Def.Nsite) {
            local_site = site0;
            inter_site = site1;
        }
        else if (site0 + 1 > (int)X->Def.Nsite && site1 + 1 <= (int)X->Def.Nsite) {
            local_site = site1;
            inter_site = site0;
            is_reverse = 1;
            J = conj(J);
        }
        else {
            continue;
        }

        mask = (int)X->Def.Tpow[inter_site];
        origin = myrank ^ mask;
        int state2 = (origin & mask) / mask;

        int org_ispin3 = is_reverse ? 0 : 1;
        int org_ispin4 = is_reverse ? 1 : 0;

        double complex Jint;
        int state1check;
        if (state2 == org_ispin4) {
            state1check = is_reverse ? 0 : 1;
            Jint = J;
        }
        else if (state2 == org_ispin3) {
            state1check = is_reverse ? 1 : 0;
            Jint = conj(J);
        }
        else {
            Jint = 0.0;
            state1check = 0;
        }

        for (g = 0; g < num_unique; g++) {
            if (batched->groups[g].origin == origin) break;
        }

        t = origin_count[origin]++;
        batched->groups[g].term_indices[t] = i;
        batched->groups[g].coefficients[t] = Jint;
        batched->groups[g].org_isite1[t] = local_site;
        batched->groups[g].org_ispin1[t] = is_reverse ? 1 : 0;
        batched->groups[g].org_ispin2[t] = is_reverse ? 0 : 1;
        batched->groups[g].state1check[t] = state1check;
    }

    free(origin_count);
    free(unique_origins);

    if (myrank == 0) {
        fprintf(stdout, "  [MPI Batching] Spin Exchange: %d MPIsingle terms -> %d groups (%.1fx reduction)\n",
                total_terms, batched->num_groups,
                batched->num_groups > 0 ? (double)total_terms / batched->num_groups : 0.0);
    }

    batched->is_initialized = 1;
    return 0;
}

/**
 * @brief Perform batched MPI Exchange for Spin (canonical)
 *
 * For canonical Spin, we need to exchange list_1/list_1buf.
 */
double complex X_child_general_int_spin_MPIsingle_batched(
    MPISpinExchangeGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
) {
    unsigned long int idim_max_buf, j, ioff;
    int ierr, t;
    MPI_Status statusMPI;
    double complex dam_pr = 0.0;

    if (group == NULL || group->num_terms == 0) return 0.0;

    // Exchange dimensions, list_1, and wavefunction
    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, group->origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    // Process all terms using the received data
    for (t = 0; t < group->num_terms; t++) {
        double complex Jint = group->coefficients[t];
        int org_isite1 = group->org_isite1[t];
        int state1check = group->state1check[t];
        unsigned long int mask1 = X->Def.Tpow[org_isite1];

        if (cabs(Jint) < 1e-15) continue;

        if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, ioff) \
    firstprivate(idim_max_buf, Jint, mask1, state1check, X) \
    shared(list_1buf, list_2_1, list_2_2, v1buf, tmp_v1, tmp_v0)
            for (j = 1; j <= idim_max_buf; j++) {
                unsigned long int list_1_j = list_1buf[j];
                unsigned long int ibit_tmp_1 = list_1_j & mask1;
                if ((ibit_tmp_1 == 0 && state1check == 0) ||
                    (ibit_tmp_1 != 0 && state1check == 1)) {
                    unsigned long int off_state;
                    if (state1check == 0) {
                        off_state = list_1_j + mask1;
                    } else {
                        off_state = list_1_j - mask1;
                    }
                    GetOffComp(list_2_1, list_2_2, off_state,
                               X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
                    double complex dmv = Jint * v1buf[j];
                    tmp_v0[ioff] += dmv;
                    dam_pr += conj(tmp_v1[ioff]) * dmv;
                }
            }
        } else {
#pragma omp parallel for default(none) reduction(+:dam_pr) \
    private(j, ioff) \
    firstprivate(idim_max_buf, Jint, mask1, state1check, X) \
    shared(list_1buf, list_2_1, list_2_2, v1buf, tmp_v1, tmp_v0)
            for (j = 1; j <= idim_max_buf; j++) {
                unsigned long int list_1_j = list_1buf[j];
                unsigned long int ibit_tmp_1 = list_1_j & mask1;
                if ((ibit_tmp_1 == 0 && state1check == 0) ||
                    (ibit_tmp_1 != 0 && state1check == 1)) {
                    unsigned long int off_state;
                    if (state1check == 0) {
                        off_state = list_1_j + mask1;
                    } else {
                        off_state = list_1_j - mask1;
                    }
                    GetOffComp(list_2_1, list_2_2, off_state,
                               X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
                    double complex dmv = Jint * v1buf[j];
                    dam_pr += conj(tmp_v1[ioff]) * dmv;
                }
            }
        }
    }

    return dam_pr;
}

#endif /* MPI */
