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
 * @file phys_distributed.h
 *
 * @brief State-task-parallel (ExpecMode 1) FullDiag observables driver.
 *
 * Two layers, split across two translation units so that the guard
 * test/check_expec_local_calls.sh can scan the MPI-free local loop:
 *
 * - phys_stateparallel() (src/phys_distributed.c, MPI orchestration):
 *   state-panel redistribution, the single MIN-Allreduce rendezvous, the
 *   all_* Gatherv to rank 0, and the collective GreenOutputMergePartials
 *   call. Compiled only under _SCALAPACK (empty TU otherwise, like
 *   matrixscalapack.c). Permanently NOT scanned by the guard.
 *
 * - phys_stateparallel_local_loop() (src/phys_distributed_local.c, MPI-free
 *   local loop): the ExpecLocal-wrapped per-rank state loop that calls the
 *   existing expec_* observable evaluators. Contains ZERO raw MPI_* /
 *   exitMPI / non-allowlisted wrapperMPI calls and does not include mpi.h;
 *   it is the file added to check_expec_local_calls.sh's FILES. It references
 *   only build-independent globals (v0/v1, X->Phys.*) and always-available
 *   APIs (expec_*, ExpecLocal*, GreenOutput* partial helpers), so it compiles
 *   in every build (including the non-MPI / non-_SCALAPACK default build);
 *   it is simply never *called* there because the only caller,
 *   phys_stateparallel(), is compiled out.
 */
#pragma once
#include "struct.h"

/**
 * @brief MPI-free per-rank observable state loop (ExpecMode 1 local layer).
 *
 * Enters ExpecLocal mode, opens a green_output partial session for @p X's
 * rank, evaluates all FullDiag observables for the states this rank owns
 * (columns of @p panel), and records them into X->Phys.all_*[n-1]. Always
 * pairs GreenOutputClearPartialSuffix()/ExpecLocalLeave() before returning,
 * on both the success and the early-break (rc=-1) paths.
 *
 * @param[in,out] X    calculation parameters / result struct
 * @param[in] panel    this rank's column-major state panel (ld = NN); column
 *                     (n-jb) holds the full eigenvector of 1-based state n
 * @param[in] jb       first owned state (1-based)
 * @param[in] je       last owned state (1-based, inclusive); je<jb ==> no
 *                     states owned (zero-owner rank), loop body is skipped
 * @param[in] NN       full eigenvector length (= X->Check.idim_max)
 * @return 0 on success, -1 if any observable evaluation (or a deferred
 *         ExpecLocal error) failed for an owned state (rank-local verdict)
 */
int phys_stateparallel_local_loop(struct BindStruct *X,
                                  double complex *panel,
                                  long int jb, long int je, long int NN);

#ifdef _SCALAPACK
/**
 * @brief State-task-parallel FullDiag observables driver (ExpecMode 1).
 *
 * Redistributes the distributed eigenvector matrix (Z_vec/descZ_vec) into a
 * per-rank state-column panel, runs phys_stateparallel_local_loop() on the
 * owned states, then (collectively) shares the failure verdict, gathers the
 * all_* summary arrays to rank 0, re-renders the progress lines in state
 * order on rank 0, and merges the green_output partial files.
 *
 * @param[in,out] X   calculation parameters / result struct
 * @param[in] neig    number of eigenstates (= X->Check.idim_max = N)
 * @return 0 on success, -1 on any (globally-shared) failure
 */
int phys_stateparallel(struct BindStruct *X, unsigned long int neig);
#endif /* _SCALAPACK */
