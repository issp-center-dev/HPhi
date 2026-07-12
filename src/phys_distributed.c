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
 * @file phys_distributed.c
 *
 * @brief ExpecMode 1 (state-task-parallel) FullDiag observables driver:
 * MPI orchestration layer.
 *
 * This translation unit owns every genuine MPI collective of the
 * state-parallel path: the state-panel redistribution, the single
 * MIN-Allreduce failure rendezvous, the all_* Gatherv to rank 0, and the
 * collective GreenOutputMergePartials(). The MPI-free per-rank observable
 * loop lives in src/phys_distributed_local.c (which the ExpecLocal guard
 * scans); this file is permanently NOT scanned. The whole TU is compiled
 * only under _SCALAPACK (an empty TU otherwise, mirroring matrixscalapack.c),
 * because it references the distributed eigenvector globals Z_vec/descZ_vec
 * and the ScaLAPACK redistribution primitive.
 */
#include "phys_distributed.h"
#ifdef _SCALAPACK
#include "matrixscalapack.h"   /* mpi.h, global.h, Z_vec/descZ_vec, use_scalapack,
                                  RedistBlockCyclicToStatePanel */
#include "green_output.h"
#include "wrapperMPI.h"
#include "DefCommon.h"
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

int phys_stateparallel(struct BindStruct *X, unsigned long int neig) {
  long int NN = (long int)neig;
  int P = nproc;
  long int NC = (NN + P - 1) / P;
  long int jb = (long int)myrank * NC + 1;
  long int je = ((long int)myrank + 1) * NC;
  long int ncols;
  int rc_local = 0;
  double complex *panel = NULL;
  int ok, gok, g;
  long int r;

  if (je > NN) je = NN;
  ncols = (je >= jb) ? (je - jb + 1) : 0;

  /* --- allocate the state panel; sync the malloc verdict with a MIN-Allreduce
     BEFORE any collective touches it (same pattern as lapack_diag_elpa). --- */
  panel = malloc((((NN * ncols) > 0) ? (size_t)(NN * ncols) : 1)
                 * sizeof(double complex));
  ok = (panel != NULL) ? 0 : -1;
  if (ok != 0)
    fprintf(stdout, "  Error: malloc failed for state panel (rank %d).\n", myrank);
  MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (gok != 0) {
    free(panel);
    return -1;
  }

  /* Redistribute the 2D block-cyclic eigenvectors into this rank's state
     column panel, then release Z_vec (peak memory drops back to O(N^2/P)). */
  if (RedistBlockCyclicToStatePanel(NN, Z_vec, descZ_vec, jb, ncols, NN, panel)
      != 0) {
    free(panel);
    return -1;
  }
  free(Z_vec);
  Z_vec = NULL;

  /* --- MPI-free per-rank observable loop (partial green_output session). --- */
  rc_local = phys_stateparallel_local_loop(X, panel, jb, je, NN);
  free(panel);

  /* --- single rendezvous: share the failure verdict across all ranks. --- */
  MPI_Allreduce(&rc_local, &g, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (g != 0)
    return -1;

  /* --- Gatherv the all_* summary arrays to rank 0. Each rank owns the
     block-contiguous segment [jb-1, je-1]; recvcounts/displs are derived from
     the deterministic block ownership. All all_* arrays are double, so the
     same recvcounts/displs are reused for every array. --- */
  {
    int *recvcounts = NULL, *displs = NULL;
    double *tmp = NULL;
    int sendcount;

    /* long->int width check: the largest int-valued quantity fed to
       MPI_Gatherv is a displ/count bounded by NN. */
    ok = 0;
    if (NN > (long int)INT_MAX) {
      ok = -1;
      if (myrank == 0)
        fprintf(stdout,
                "  Error: N=%ld exceeds INT_MAX; cannot Gatherv all_* arrays.\n",
                NN);
    }
    /* rank-0 receive buffers; sync malloc success before the collective. */
    if (ok == 0 && myrank == 0) {
      recvcounts = malloc((size_t)P * sizeof(int));
      displs = malloc((size_t)P * sizeof(int));
      tmp = malloc((size_t)NN * sizeof(double));
      if (recvcounts == NULL || displs == NULL || tmp == NULL) {
        ok = -1;
        fprintf(stdout,
                "  Error: malloc failed for all_* Gatherv buffers (rank 0).\n");
      }
    }
    MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (gok != 0) {
      free(recvcounts);
      free(displs);
      free(tmp);
      return -1;
    }

    if (myrank == 0) {
      for (r = 0; r < P; r++) {
        long int jb_r = r * NC + 1;
        long int je_r = (r + 1) * NC;
        long int nc_r;
        if (je_r > NN) je_r = NN;
        nc_r = (je_r >= jb_r) ? (je_r - jb_r + 1) : 0;
        recvcounts[r] = (int)nc_r;
        displs[r] = (int)(r * NC);  /* 0-based start; unused when nc_r==0 */
      }
    }
    sendcount = (int)ncols;

    /* rank 0 receives into a separate temp buffer, then copies into
       X->Phys.all_* (avoids MPI_IN_PLACE aliasing with its own segment).
       Senders pass &all_xxx[jb-1] with sendcount=ncols (0 for zero-owner). */
#define PHYS_GATHER_ONE(arr)                                                   \
    MPI_Gatherv((ncols > 0) ? &X->Phys.arr[jb - 1] : X->Phys.arr, sendcount,   \
                MPI_DOUBLE, tmp, recvcounts, displs, MPI_DOUBLE, 0,            \
                MPI_COMM_WORLD);                                               \
    if (myrank == 0) {                                                         \
      long int k;                                                             \
      for (k = 0; k < NN; k++) X->Phys.arr[k] = tmp[k];                        \
    }
    PHYS_GATHER_ONE(all_energy)
    PHYS_GATHER_ONE(all_doublon)
    PHYS_GATHER_ONE(all_sz)
    PHYS_GATHER_ONE(all_s2)
    PHYS_GATHER_ONE(all_num_up)
    PHYS_GATHER_ONE(all_num_down)
#undef PHYS_GATHER_ONE

    free(recvcounts);
    free(displs);
    free(tmp);
  }

  /* --- rank 0 re-renders the progress lines in state order, in the serial
     format (S2 column included). Mode 1 is only ever entered for FullDiag. --- */
  if (myrank == 0 && X->Def.iCalcType == FullDiag) {
    long int i;
    double tmp_N;
    for (i = 0; i < NN; i++) {
      if (X->Def.iCalcModel == Spin || X->Def.iCalcModel == SpinGC) {
        tmp_N = X->Def.NsiteMPI;
      } else {
        tmp_N = X->Phys.all_num_up[i] + X->Phys.all_num_down[i];
      }
      fprintf(stdoutMPI,
              "i=%5ld Energy=%10lf N=%10lf Sz=%10lf S2=%10lf Doublon=%10lf \n",
              i, X->Phys.all_energy[i], tmp_N, X->Phys.all_sz[i],
              X->Phys.all_s2[i], X->Phys.all_doublon[i]);
    }
  }

  /* --- publish the aggregate green_output files: rank 0 concatenates the
     rank-ordered (= state-ordered) partials. Collective; identical rc on all
     ranks; publishes/deletes only on global success. --- */
  if (GreenOutputMergePartials(X) != 0)
    return -1;

  return 0;
}

#endif /* _SCALAPACK */
