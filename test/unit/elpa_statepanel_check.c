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

/* Verify RedistBlockCyclicToStatePanel: the mirror image of
   elpa_redist_check.c's RedistPanelToBlockCyclic check. Build a
   deterministic matrix (i) in the 2D block-cyclic layout via
   DivMat/pzelset_ (playing the role of an eigenvector matrix Z), then
   (ii) redistribute it with RedistBlockCyclicToStatePanel into a 1D
   state-column panel. For each 1-based state this rank owns
   (jbegin..jbegin+ncols-1), compare the panel column directly against
   MatElem(i, state-1) -- i.e. the corresponding column of the original
   (undistributed) matrix. No diagonalization is involved: this checks
   pure data movement. NDIM defaults to 97 (non-divisible vs both the
   process count and the capped block size) but can be overridden at
   compile time (-DNDIM=...) for small deliberate edge cases such as the
   zero-state-owner scenario. Run with any np. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>
#include "matrixscalapack.h"
#include "matrixlapack_elpa.h"

#ifndef NDIM
#define NDIM 97
#endif

/* RedistBlockCyclicToStatePanel (src/matrixscalapack.c), like
   RedistPanelToBlockCyclic, checks its 1D BLACS grid's column against
   the process-global `myrank` (extern in global.h, defined in
   src/global.c) to guard against a scrambled ownership mapping. This
   test links only matrixscalapack.c + matrixlapack_elpa.c (mirroring
   elpa_redist_check's registration), so `myrank` is not otherwise
   defined; provide it here and set it from MPI_Comm_rank below, the
   same way wrapperMPI.c's InitializeMPI() does. */
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
  int descZ[9];
  double complex *Z_ref, *panel;
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
  descinit_(descZ, &n, &n, &mb, &mb, &i_zero_i, &i_zero_i, &ictxt, &lld, &info);

  Z_ref = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));

  /* fill the 2D block-cyclic "eigenvector matrix" Z: column j plays the
     role of eigenvector for state j */
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      DivMat(i, j, MatElem(i, j), Z_ref, descZ);

  /* this rank's owned 0-based state range [jb, je) per the frozen
     ownership formula NC=ceil(N/P), first_state(r)=r*NC+1 (1-based) */
  NC = (n + size - 1) / size;
  jb = (long int)rank * NC;
  je = jb + NC; if (je > n) je = n;
  ncols = (je > jb) ? (je - jb) : 0;
  lde = n;
  panel = malloc((((n * ncols) > 0) ? n * ncols : 1) * sizeof(double complex));

  redist_rc = RedistBlockCyclicToStatePanel(n, Z_ref, descZ, jb + 1, ncols, lde, panel);
  if (redist_rc != 0) {
    if (rank == 0)
      fprintf(stderr, "RedistBlockCyclicToStatePanel failed: rc=%d\n", redist_rc);
    ok = 0;
  }

  if (ok) {
    for (j = jb; j < je && ok; j++) {
      for (i = 0; i < n; i++) {
        double complex expect = MatElem(i, j);
        double complex got = panel[(j - jb) * lde + i];
        if (cabs(expect - got) > 1e-14) {
          fprintf(stderr, "rank %d: mismatch at state %ld row %ld: %g\n",
                  rank, j, i, cabs(expect - got));
          ok = 0; break;
        }
      }
    }
  }
  { int gok; MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); ok = gok; }
  if (rank == 0) printf("elpa_statepanel_check: %s\n", ok ? "OK" : "FAILED");
  free(Z_ref); free(panel);
  blacs_gridexit_(&ictxt);
  MPI_Finalize();
  return ok ? 0 : 1;
}
