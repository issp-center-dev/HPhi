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
 * @file mltplyMPIBatched.h
 * @brief Data structures and functions for batched MPI communication optimization
 *
 * This module groups multiple transfer terms that share the same MPI communication
 * partner (origin) to reduce the number of MPI_Sendrecv calls.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */

#pragma once

#include <complex.h>
#include "struct.h"

/**
 * @brief Structure representing a group of transfers to the same MPI origin
 *
 * When multiple transfer terms require communication with the same MPI rank,
 * they are grouped together to share a single MPI_Sendrecv call.
 */
typedef struct {
    int origin;                      /**< MPI rank of communication partner */
    int mask;                        /**< Bit mask for this origin (Tpow[org_isite2]) */
    int num_transfers;               /**< Number of transfers in this group */
    int *transfer_indices;           /**< Indices into EDGeneralTransfer array */
    double complex *coefficients;    /**< Pre-computed coefficients (Fsgn * trans) - may be stale for time evolution */
    int *Fsgn;                       /**< Fermion sign for each transfer */
    int *is_conj;                    /**< Whether to conjugate the coefficient */
    unsigned long int *local_mask;   /**< Local site masks for each transfer */
    unsigned long int *state1check;  /**< Expected state of local site for each transfer */
    unsigned long int *bit1diff;     /**< Bit difference mask for fermion sign */
} MPITransferGroup;

/**
 * @brief Container for all batched transfer groups
 */
typedef struct {
    int num_groups;                  /**< Number of unique origins (groups) */
    MPITransferGroup *groups;        /**< Array of transfer groups */
    int is_initialized;              /**< Flag indicating initialization status */
} MPIBatchedTransfers;

/** @brief Returns 1 if MPI communication batching is enabled, 0 if disabled
 *  (HPHI_MPI_NOBATCH=1). Parsed on rank 0 and broadcast in InitializeMPI;
 *  declared outside #ifdef MPI so non-MPI builds can call it too. */
int MPIBatchingEnabled(void);

#ifdef MPI

/**
 * @brief Initialize batched transfers for SpinlessFermionGC MPIsingle mode
 *
 * Scans all transfer terms and groups them by their MPI communication partner (origin).
 * Pre-computes coefficients including fermion signs.
 *
 * @param X Pointer to BindStruct containing transfer definitions
 * @param batched Pointer to MPIBatchedTransfers to initialize
 * @return 0 on success, -1 on error
 */
int InitializeMPIBatchedTransfers_SpinlessFermionGC(
    struct BindStruct *X,
    MPIBatchedTransfers *batched
);

/**
 * @brief Free memory allocated for batched transfers
 *
 * @param batched Pointer to MPIBatchedTransfers to finalize
 */
void FinalizeMPIBatchedTransfers(MPIBatchedTransfers *batched);

/**
 * @brief Perform batched MPI hopping for SpinlessFermionGC MPIsingle mode
 *
 * Performs a single MPI_Sendrecv for all transfers in the group, then
 * applies all transfer operations locally using the received data.
 *
 * @param group Pointer to the transfer group to process
 * @param X Pointer to BindStruct
 * @param tmp_v0 Output vector (updated in M_MLTPLY mode)
 * @param tmp_v1 Input vector
 * @return Sum of dam_pr contributions from all transfers
 */
double complex X_child_GC_general_hopp_SpinlessFermion_MPIsingle_batched(
    MPITransferGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
);

/**
 * @brief Initialize batched transfers for SpinlessFermion (canonical) MPIsingle mode
 *
 * Similar to GC version but requires list_1/list_1buf exchange.
 *
 * @param X Pointer to BindStruct containing transfer definitions
 * @param batched Pointer to MPIBatchedTransfers to initialize
 * @return 0 on success, -1 on error
 */
int InitializeMPIBatchedTransfers_SpinlessFermion(
    struct BindStruct *X,
    MPIBatchedTransfers *batched
);

/**
 * @brief Perform batched MPI hopping for SpinlessFermion (canonical) MPIsingle mode
 *
 * Performs a single MPI_Sendrecv for all transfers in the group (including list_1),
 * then applies all transfer operations locally using the received data.
 *
 * @param group Pointer to the transfer group to process
 * @param X Pointer to BindStruct
 * @param tmp_v0 Output vector (updated in M_MLTPLY mode)
 * @param tmp_v1 Input vector
 * @return Sum of dam_pr contributions from all transfers
 */
double complex X_child_general_hopp_Spinless_MPIsingle_batched(
    MPITransferGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
);

/**
 * @brief Initialize batched transfers for HubbardGC MPIsingle mode
 *
 * Similar to SpinlessFermionGC but with spin handling (mask = Tpow[2*site + spin]).
 *
 * @param X Pointer to BindStruct containing transfer definitions
 * @param batched Pointer to MPIBatchedTransfers to initialize
 * @return 0 on success, -1 on error
 */
int InitializeMPIBatchedTransfers_HubbardGC(
    struct BindStruct *X,
    MPIBatchedTransfers *batched
);

/**
 * @brief Perform batched MPI hopping for HubbardGC MPIsingle mode
 *
 * Performs a single MPI_Sendrecv for all transfers in the group, then
 * applies all transfer operations locally using the received data.
 *
 * @param group Pointer to the transfer group to process
 * @param X Pointer to BindStruct
 * @param tmp_v0 Output vector (updated in M_MLTPLY mode)
 * @param tmp_v1 Input vector
 * @return Sum of dam_pr contributions from all transfers
 */
double complex X_child_GC_general_hopp_MPIsingle_batched(
    MPITransferGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
);

/**
 * @brief Initialize batched transfers for Hubbard (canonical) MPIsingle mode
 *
 * Similar to HubbardGC but requires list_1/list_1buf exchange.
 *
 * @param X Pointer to BindStruct containing transfer definitions
 * @param batched Pointer to MPIBatchedTransfers to initialize
 * @return 0 on success, -1 on error
 */
int InitializeMPIBatchedTransfers_Hubbard(
    struct BindStruct *X,
    MPIBatchedTransfers *batched
);

/**
 * @brief Perform batched MPI hopping for Hubbard (canonical) MPIsingle mode
 *
 * @param group Pointer to the transfer group to process
 * @param X Pointer to BindStruct
 * @param tmp_v0 Output vector (updated in M_MLTPLY mode)
 * @param tmp_v1 Input vector
 * @return Sum of dam_pr contributions from all transfers
 */
double complex X_child_general_hopp_MPIsingle_batched(
    MPITransferGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
);

/**
 * @brief Structure representing a group of MPIdouble transfers to the same MPI origin
 *
 * For MPIdouble, both sites are inter-process. After MPI_Sendrecv,
 * processing is simple vector scaling (no per-element bit manipulation).
 */
typedef struct {
    int origin;                      /**< MPI rank of communication partner */
    int num_transfers;               /**< Number of transfers in this group */
    int *transfer_indices;           /**< Indices into EDGeneralTransfer array */
    double complex *coefficients;    /**< Pre-computed coefficients (Fsgn * trans) - may be stale for time evolution */
    int *Fsgn;                       /**< Fermion sign for each transfer */
    int *is_conj;                    /**< Whether to conjugate the coefficient */
} MPIDoubleTransferGroup;

/**
 * @brief Container for all batched MPIdouble transfer groups
 */
typedef struct {
    int num_groups;                  /**< Number of unique origins (groups) */
    MPIDoubleTransferGroup *groups;  /**< Array of transfer groups */
    int is_initialized;              /**< Flag indicating initialization status */
} MPIBatchedDoubleTransfers;

/**
 * @brief Initialize batched MPIdouble transfers for HubbardGC
 *
 * @param X Pointer to BindStruct containing transfer definitions
 * @param batched Pointer to MPIBatchedDoubleTransfers to initialize
 * @return 0 on success, -1 on error
 */
int InitializeMPIBatchedDoubleTransfers_HubbardGC(
    struct BindStruct *X,
    MPIBatchedDoubleTransfers *batched
);

/**
 * @brief Free memory allocated for batched MPIdouble transfers
 *
 * @param batched Pointer to MPIBatchedDoubleTransfers to finalize
 */
void FinalizeMPIBatchedDoubleTransfers(MPIBatchedDoubleTransfers *batched);

/**
 * @brief Perform batched MPI hopping for HubbardGC MPIdouble mode
 *
 * @param group Pointer to the transfer group to process
 * @param X Pointer to BindStruct
 * @param tmp_v0 Output vector (updated in M_MLTPLY mode)
 * @param tmp_v1 Input vector
 * @return Sum of dam_pr contributions from all transfers
 */
double complex X_child_GC_general_hopp_MPIdouble_batched(
    MPIDoubleTransferGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
);

/**
 * @brief Initialize batched MPIdouble transfers for Hubbard (canonical)
 *
 * Similar to HubbardGC but requires list_1/list_1buf exchange.
 *
 * @param X Pointer to BindStruct containing transfer definitions
 * @param batched Pointer to MPIBatchedDoubleTransfers to initialize
 * @return 0 on success, -1 on error
 */
int InitializeMPIBatchedDoubleTransfers_Hubbard(
    struct BindStruct *X,
    MPIBatchedDoubleTransfers *batched
);

/**
 * @brief Perform batched MPI hopping for Hubbard (canonical) MPIdouble mode
 *
 * Performs a single MPI_Sendrecv for all transfers in the group (including list_1),
 * then applies all transfer operations locally using the received data.
 *
 * @param group Pointer to the transfer group to process
 * @param X Pointer to BindStruct
 * @param tmp_v0 Output vector (updated in M_MLTPLY mode)
 * @param tmp_v1 Input vector
 * @return Sum of dam_pr contributions from all transfers
 */
double complex X_child_general_hopp_MPIdouble_batched(
    MPIDoubleTransferGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
);

/**
 * @brief Structure representing a group of InterAll terms to the same MPI origin
 *
 * When multiple InterAll terms require communication with the same MPI rank,
 * they are grouped together to share a single MPI_Sendrecv call.
 */
typedef struct {
    int origin;                      /**< MPI rank of communication partner */
    int num_interall;                /**< Number of InterAll terms in this group */
    int *interall_indices;           /**< Indices into InterAll_OffDiagonal array */
    double complex *coefficients;    /**< Pre-computed coefficients (tmp_V) */
    int *is_hermite;                 /**< Flag for hermitian conjugate term */
    unsigned long int *tmp_isite1;   /**< Pre-computed OrgTpow values */
    unsigned long int *tmp_isite2;
    unsigned long int *tmp_isite3;
    unsigned long int *tmp_isite4;
    unsigned long int *isite1;       /**< Pre-computed Tpow values */
    unsigned long int *isite2;
    unsigned long int *isite3;
    unsigned long int *isite4;
} MPIInterAllGroup;

/**
 * @brief Container for all batched InterAll groups
 */
typedef struct {
    int num_groups;                  /**< Number of unique origins (groups) */
    MPIInterAllGroup *groups;        /**< Array of InterAll groups */
    int is_initialized;              /**< Flag indicating initialization status */
} MPIBatchedInterAll;

/**
 * @brief Initialize batched InterAll for HubbardGC
 *
 * Scans all InterAll terms and groups them by their MPI communication partner (origin).
 *
 * @param X Pointer to BindStruct containing InterAll definitions
 * @param batched Pointer to MPIBatchedInterAll to initialize
 * @return 0 on success, -1 on error
 */
int InitializeMPIBatchedInterAll_HubbardGC(
    struct BindStruct *X,
    MPIBatchedInterAll *batched
);

/**
 * @brief Free memory allocated for batched InterAll
 *
 * @param batched Pointer to MPIBatchedInterAll to finalize
 */
void FinalizeMPIBatchedInterAll(MPIBatchedInterAll *batched);

/**
 * @brief Perform batched MPI InterAll for HubbardGC
 *
 * Performs a single MPI_Sendrecv for all InterAll terms in the group, then
 * applies all operations locally using the received data.
 *
 * @param group Pointer to the InterAll group to process
 * @param X Pointer to BindStruct
 * @param tmp_v0 Output vector (updated in M_MLTPLY mode)
 * @param tmp_v1 Input vector
 * @return Sum of dam_pr contributions from all InterAll terms
 */
double complex X_child_GC_InterAll_Hubbard_MPI_batched(
    MPIInterAllGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
);

/*******************************************************************************
 * Spin model batched MPI structures and functions
 ******************************************************************************/

/**
 * @brief Structure representing a group of Spin Exchange/PairLift terms to the same MPI origin
 *
 * For Spin 1/2 models, mask = Tpow[site] (no spin multiplication)
 */
typedef struct {
    int origin;                      /**< MPI rank of communication partner */
    int num_terms;                   /**< Number of Exchange/PairLift terms in this group */
    int *term_indices;               /**< Indices into Exchange/PairLift array */
    double complex *coefficients;    /**< Pre-computed coefficients */
    int *org_isite1;                 /**< Local site indices */
    int *org_ispin1;                 /**< Local spin indices (before flip) */
    int *org_ispin2;                 /**< Local spin indices (after flip) */
    int *state1check;                /**< State check for local site */
} MPISpinExchangeGroup;

/**
 * @brief Container for all batched Spin Exchange groups
 */
typedef struct {
    int num_groups;                  /**< Number of unique origins (groups) */
    MPISpinExchangeGroup *groups;    /**< Array of Exchange groups */
    int is_initialized;              /**< Flag indicating initialization status */
} MPIBatchedSpinExchange;

/**
 * @brief Initialize batched Exchange for SpinGC
 *
 * @param X Pointer to BindStruct containing Exchange definitions
 * @param batched Pointer to MPIBatchedSpinExchange to initialize
 * @return 0 on success, -1 on error
 */
int InitializeMPIBatchedExchange_SpinGC(
    struct BindStruct *X,
    MPIBatchedSpinExchange *batched
);

/**
 * @brief Initialize batched Exchange for Spin (canonical)
 *
 * @param X Pointer to BindStruct containing Exchange definitions
 * @param batched Pointer to MPIBatchedSpinExchange to initialize
 * @return 0 on success, -1 on error
 */
int InitializeMPIBatchedExchange_Spin(
    struct BindStruct *X,
    MPIBatchedSpinExchange *batched
);

/**
 * @brief Free memory allocated for batched Spin Exchange
 *
 * @param batched Pointer to MPIBatchedSpinExchange to finalize
 */
void FinalizeMPIBatchedSpinExchange(MPIBatchedSpinExchange *batched);

/**
 * @brief Perform batched MPI Exchange for SpinGC
 *
 * @param group Pointer to the Exchange group to process
 * @param X Pointer to BindStruct
 * @param tmp_v0 Output vector
 * @param tmp_v1 Input vector
 * @return Sum of dam_pr contributions
 */
double complex X_child_GC_CisAitCiuAiv_spin_MPIsingle_batched(
    MPISpinExchangeGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
);

/**
 * @brief Perform batched MPI Exchange for Spin (canonical)
 *
 * @param group Pointer to the Exchange group to process
 * @param X Pointer to BindStruct
 * @param tmp_v0 Output vector
 * @param tmp_v1 Input vector
 * @return Sum of dam_pr contributions
 */
double complex X_child_general_int_spin_MPIsingle_batched(
    MPISpinExchangeGroup *group,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
);

#endif /* MPI */
