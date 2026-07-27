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

/* Verify RedistPanelToBlockCyclic: build the same deterministic matrix
   (i) replicated + DivMat/pzelset and (ii) as a 1D column panel +
   pzgemr2d, then compare the resulting 2D block-cyclic local arrays
   element-wise. N=97 exercises non-divisible N vs both the process
   count and the (capped) block size. Run with any np. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>
#include "matrixscalapack.h"
#include "matrixlapack_elpa.h"

#define NDIM 97

/* RedistPanelToBlockCyclic (src/matrixscalapack.c) checks its 1D BLACS
   grid's column against the process-global `myrank` (extern in global.h,
   defined in src/global.c) to guard against a scrambled ownership
   mapping. This test links only matrixscalapack.c + matrixlapack_elpa.c
   (mirroring elpa_eigen_check's registration), so `myrank` is not
   otherwise defined; provide it here and set it from MPI_Comm_rank
   below, the same way wrapperMPI.c's InitializeMPI() does. */
int myrank = 0;

static double complex MatElem(long int i, long int j) { /* 0-based */
  return (1.0 / (1.0 + labs(i - j))) + ((double)(i - j) / (NDIM * NDIM)) * I;
}

int main(int argc, char **argv) {
  int i_negone = -1, i_zero_i = 0;
  int rank, size, ictxt, iam, nprocs, info, ok = 1;
  int nprow, npcol, myrow, mycol;
  long int n = NDIM, mb, mp, nq, i, j, k;
  int lld, dims[2] = {0, 0};
  int descA[9], descB[9];
  double complex *A_ref, *B_panel2d, *panel;
  long int NC, jb, je, ncols, lde;
  int redist_rc;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  myrank = rank; /* see comment on the `myrank` definition above */
  MPI_Dims_create(size, 2, dims);
  nprow = dims[0]; npcol = dims[1];
  mb = ElpaBlockSize(n, nprow, npcol);
  blacs_pinfo_(&iam, &nprocs);
  blacs_get_(&i_negone, &i_zero_i, &ictxt);
  blacs_gridinit_(&ictxt, "R", &nprow, &npcol);
  blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);

  mp = numroc_(&n, &mb, &myrow, &i_zero_i, &nprow);
  nq = numroc_(&n, &mb, &mycol, &i_zero_i, &npcol);
  lld = (mp > 0) ? (int)mp : 1;
  descinit_(descA, &n, &n, &mb, &mb, &i_zero_i, &i_zero_i, &ictxt, &lld, &info);
  descinit_(descB, &n, &n, &mb, &mb, &i_zero_i, &i_zero_i, &ictxt, &lld, &info);

  A_ref = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));
  B_panel2d = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));

  /* (i) replicated fill */
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      DivMat(i, j, MatElem(i, j), A_ref, descA);

  /* (ii) 1D column panel fill: rank owns 0-based columns [jb, je) */
  NC = (n + size - 1) / size;
  jb = (long int)rank * NC;
  je = jb + NC; if (je > n) je = n;
  ncols = (je > jb) ? (je - jb) : 0;
  lde = n;
  panel = malloc((((n * ncols) > 0) ? n * ncols : 1) * sizeof(double complex));
  for (j = jb; j < je; j++)
    for (i = 0; i < n; i++)
      panel[(j - jb) * lde + i] = MatElem(i, j);

  redist_rc = RedistPanelToBlockCyclic(n, jb + 1, ncols, lde, panel, B_panel2d, descB);
  if (redist_rc != 0) {
    if (rank == 0)
      fprintf(stderr, "RedistPanelToBlockCyclic failed: rc=%d\n", redist_rc);
    ok = 0;
  }

  if (ok) {
    for (k = 0; k < mp * nq; k++) {
      if (cabs(A_ref[k] - B_panel2d[k]) > 1e-14) {
        fprintf(stderr, "rank %d: mismatch at local index %ld: %g\n",
                rank, k, cabs(A_ref[k] - B_panel2d[k]));
        ok = 0; break;
      }
    }
  }
  /* Failure-propagation phase: rank 0 passes an invalid panel leading
     dimension (1 < LOCr = n), so descinit_ fails with info != 0 on rank 0
     ONLY. The synchronized verdict inside RedistPanelToBlockCyclic must
     make EVERY rank return -1 without entering the collective pzgemr2d_
     (a hang here means the verdict Allreduce is broken). Requires n > 1
     so that lde=1 is genuinely invalid on rank 0. Mirrors the analogous
     phase in elpa_statepanel_check.c for RedistBlockCyclicToStatePanel. */
  if (n > 1 && size >= 2) {
    long int lde_bad = (rank == 0) ? 1 : lde;
    int rc2 = RedistPanelToBlockCyclic(n, jb + 1, ncols, lde_bad, panel,
                                       B_panel2d, descB);
    if (rc2 != -1) {
      fprintf(stderr,
              "rank %d: expected rc=-1 from the descinit_ failure-propagation "
              "phase, got %d\n", rank, rc2);
      ok = 0;
    }
  }

  { int gok; MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); ok = gok; }
  if (rank == 0) printf("elpa_redist_check: %s\n", ok ? "OK" : "FAILED");
  free(A_ref); free(B_panel2d); free(panel);
  blacs_gridexit_(&ictxt);
  MPI_Finalize();
  return ok ? 0 : 1;
}
