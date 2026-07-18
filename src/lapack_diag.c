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
#include "lapack_diag.h"
#include "matrixlapack.h"
#include "FileIO.h"
#include "DefCommon.h"
#ifdef _MAGMA
#include "matrixlapack_magma.h"
#endif
#ifdef _SCALAPACK
#include "matrixscalapack.h"
#endif
#ifdef _ELPA
#include "matrixlapack_elpa.h"

/**
 * @brief FullDiag via ELPA (phase 1: fill the 2D block-cyclic matrix from
 * the replicated Ham with pzelset, then call diag_elpa_cmp).
 * Eigenvalues land in v0 on all ranks; eigenvectors stay in Z_vec.
 * NOTE (phase 1): the replicated Ham plus A_distr plus Z_vec coexist in
 * memory, so verification runs must stay at small N (design doc sec. 3).
 */
static int lapack_diag_elpa(struct BindStruct *X, long int xMsize) {
  int i_negone = -1, i_zero = 0;
  int size;
  int nprow, npcol, myrow, mycol;
  int ictxt;
  long int mb, mp, nq, i, j;
  int lld, dims[2] = {0, 0};
  int iam, nprocs, info;
  double complex *A_distr;
  double *w;
  int descA[9];
  int ierr;

  fprintf(stdoutMPI, "Using ELPA (%s)\n\n",
          X->Def.iNGPU >= 1 ? "GPU" : "CPU");

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size, 2, dims);
  nprow = dims[0]; npcol = dims[1];

  if (xMsize < ((nprow > npcol) ? nprow : npcol)) {
    fprintf(stdoutMPI,
            "Error: matrix dimension (%ld) is smaller than the process grid (%d x %d).\n"
            "       Reduce the number of MPI ranks for this problem size.\n",
            xMsize, nprow, npcol);
    /* Pre-existing early return, predating the descinit-verdict block
       below (PR #276 review round 2): this runs before any BLACS grid
       exists, but iHamPanelActive/Ham_local are set by the caller
       (xsetmem.c) well before lapack_diag_elpa() is entered, so the
       distributed-panel mode leak fixed in the descinit-failure branch
       applies here too. Same one-line cleanup idiom, minus
       blacs_gridexit_ (no grid to tear down yet). */
    if (iHamPanelActive) {
      free(Ham_local);
      Ham_local = NULL;
      iHamPanelActive = 0;
    }
    return -1;
  }
  mb = ElpaBlockSize(xMsize, nprow, npcol);

#ifdef _ELPA_GPU
  /* Startup consistency warning (design doc sec. 2): ranks per node
     should be a multiple of NGPU (ideally equal: 1 rank per GPU). */
  if (X->Def.iNGPU >= 1) {
    MPI_Comm comm_node;
    int nrank_node;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                        MPI_INFO_NULL, &comm_node);
    MPI_Comm_size(comm_node, &nrank_node);
    MPI_Comm_free(&comm_node);
    if (nrank_node % X->Def.iNGPU != 0) {
      fprintf(stdoutMPI,
              "Warning: ranks per node (%d) is not a multiple of NGPU (%d):\n"
              "         GPUs may idle or be shared unevenly. Recommended: 1 rank per GPU.\n",
              nrank_node, X->Def.iNGPU);
    }
  }
#endif

  blacs_pinfo_(&iam, &nprocs);
  blacs_get_(&i_negone, &i_zero, &ictxt);
  blacs_gridinit_(&ictxt, "R", &nprow, &npcol);
  blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);

  mp = numroc_(&xMsize, &mb, &myrow, &i_zero, &nprow);
  nq = numroc_(&xMsize, &mb, &mycol, &i_zero, &npcol);
  lld = (mp > 0) ? mp : 1;

  descinit_(descA, &xMsize, &xMsize, &mb, &mb, &i_zero, &i_zero, &ictxt, &lld, &info);
  {
    /* descinit_ overwrites info per call: fold both verdicts, then
       synchronize so no rank proceeds into the collectives below with an
       invalid descriptor while others abort (same pattern as
       RedistBlockCyclicToStatePanel). Capture descA's info before the
       second descinit_ call overwrites it, so a descA-only failure is
       not misreported as info=0 (PR #276 review round 2). */
    int info_a = info;
    int ok = (info_a == 0) ? 0 : -1, gok;
    descinit_(descZ_vec, &xMsize, &xMsize, &mb, &mb, &i_zero, &i_zero, &ictxt, &lld, &info);
    if (info != 0) ok = -1;
    if (ok != 0) {
      fprintf(stdout,
              "  Error: descinit_ failed (info_A=%d, info_Z=%d) for the ELPA\n"
              "         descriptors on rank %d; aborting the ELPA diagonalization.\n",
              info_a, info, myrank);
    }
    MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (gok != 0) {
      /* Distributed-panel mode (phase 2): this return happens before the
         panel-mode branch below consumes Ham_local, so free it here too
         (same idiom as the panel-mode allocation-failure paths further
         down in this function; PR #276 review round 2). */
      if (iHamPanelActive) {
        free(Ham_local);
        Ham_local = NULL;
        iHamPanelActive = 0;
      }
      blacs_gridexit_(&ictxt);
      return -1;
    }
  }

  if (iHamPanelActive) {
    /* Distributed generation (phase 2): allocate only A_distr + w first,
       redistribute the 1D panel into A_distr, then free the panel and
       only THEN allocate Z_vec. This keeps the redistribution-time peak
       at 2 coexisting matrices (design doc sec. 3), not 3. Every
       rank-local malloc (including failure) is synced with
       MPI_Allreduce(MIN) before any collective touches the buffer
       (design doc sec. 4), same pattern as SyncError in
       matrixlapack_elpa.c / the ownership guard in
       RedistPanelToBlockCyclic. */
    int rerr, ok, gok;

    A_distr = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));
    w = malloc(xMsize * sizeof(double));
    ok = (A_distr != NULL && w != NULL) ? 0 : -1;
    if (ok != 0) {
      fprintf(stdout, "  Error: malloc failed for A_distr/w (rank %d).\n", myrank);
    }
    MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (gok != 0) {
      free(A_distr);
      free(w);
      free(Ham_local);
      Ham_local = NULL;
      iHamPanelActive = 0;
      return -1;
    }

    rerr = RedistPanelToBlockCyclic(xMsize, HamColBegin,
                                    (HamColEnd >= HamColBegin)
                                      ? (HamColEnd - HamColBegin + 1) : 0,
                                    HamPanelLd, Ham_local, A_distr, descA);
    free(Ham_local);
    Ham_local = NULL;
    iHamPanelActive = 0; /* panel consumed; phys.c uses Z_vec only */
    if (rerr != 0) {
      free(A_distr);
      free(w);
      return -1;
    }

    Z_vec = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));
    ok = (Z_vec != NULL) ? 0 : -1;
    if (ok != 0) {
      fprintf(stdout, "  Error: malloc failed for Z_vec (rank %d).\n", myrank);
    }
    MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (gok != 0) {
      free(A_distr);
      free(Z_vec);
      Z_vec = NULL;
      free(w);
      return -1;
    }
  } else {
    int ok, gok;

    A_distr = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));
    Z_vec = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));
    w = malloc(xMsize * sizeof(double));
    ok = (A_distr != NULL && Z_vec != NULL && w != NULL) ? 0 : -1;
    if (ok != 0) {
      fprintf(stdout, "  Error: malloc failed for A_distr/Z_vec/w (rank %d).\n", myrank);
    }
    MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (gok != 0) {
      free(A_distr);
      free(Z_vec);
      Z_vec = NULL;
      free(w);
      return -1;
    }

    for (i = 0; i < xMsize; i++) {
      for (j = 0; j < xMsize; j++) {
        DivMat(i, j, Ham[i][j], A_distr, descA);
      }
    }
  }

  ierr = diag_elpa_cmp((int)xMsize, A_distr, Z_vec, w,
                       (int)mp, (int)nq, (int)myrow, (int)mycol,
                       (int)mb, X->Def.iNGPU);
  free(A_distr);
  if (ierr != 0) {
    /* Failed diagonalization must not leave a dangling non-NULL Z_vec
       (PR #276 review): downstream code treats Z_vec != NULL as "valid
       distributed eigenvectors exist". */
    free(Z_vec);
    Z_vec = NULL;
    free(w);
    return -1;
  }

  for (i = 0; i < xMsize; i++) {
    v0[i] = w[i];
  }
  free(w);
  use_scalapack = 1;
  return 0;
}
#endif /* _ELPA */

/**
 *
 * @brief performing full diagonalization using lapack
 * @param[in,out] X
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return
 */
int lapack_diag(
struct BindStruct *X//!<[inout]
) {

  FILE *fp;
  char sdt[D_FileNameMax] = "";
  long int i, j, i_max, xMsize;
#ifdef _SCALAPACK
  int rank, size, nprocs, nprow, npcol, myrow, mycol, ictxt;
  int i_negone=-1, i_zero=0, iam;
  long int mb, nb, mp, nq;
  int dims[2]={0,0};
#endif

  i_max = X->Check.idim_max;
  if (!iHamPanelActive) {
    /* Distributed-panel mode (phase 2): Ham is NULL and the panel is
       already 0-based-packed by rows / start-packed by columns at
       generation time (Task 4), so this shift must be skipped. */
    for (i = 0; i < i_max; i++) {
      for (j = 0; j < i_max; j++) {
        Ham[i][j] = Ham[i + 1][j + 1];
      }
    }
  }
  xMsize = i_max;
  switch (X->Def.iSolver) {
  case SOLVER_SCALAPACK:
#ifdef _SCALAPACK
    if (nproc > 1) {
      fprintf(stdoutMPI, "Using SCALAPACK\n\n");
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      MPI_Dims_create(size, 2, dims);
      nprow = dims[0]; npcol = dims[1];

      blacs_pinfo_(&iam, &nprocs);
      blacs_get_(&i_negone, &i_zero, &ictxt);
      blacs_gridinit_(&ictxt, "R", &nprow, &npcol);
      blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);

      mb = GetBlockSize(xMsize, size);
      mp = numroc_(&xMsize, &mb, &myrow, &i_zero, &nprow);
      nq = numroc_(&xMsize, &mb, &mycol, &i_zero, &npcol);
      Z_vec = malloc(mp * nq * sizeof(complex double));
      diag_scalapack_cmp(xMsize, Ham, v0, Z_vec, descZ_vec);
    } else {
      ZHEEVall(xMsize, Ham, v0, L_vec);
    }
#endif
    break;

  case SOLVER_MAGMA:
#ifdef _MAGMA
    if (myrank == 0) {
      if (diag_magma_cmp(xMsize, Ham, v0, L_vec, X->Def.iNGPU) != 0) {
        return -1;
      }
    }
#endif
    break;

  case SOLVER_ELPA:
#ifdef _ELPA
    if (lapack_diag_elpa(X, xMsize) != 0) {
      return -1;
    }
#endif
    break;

  default: /* SOLVER_LAPACK */
    ZHEEVall(xMsize, Ham, v0, L_vec);
    break;
  }
  strcpy(sdt, cFileNameEigenvalue_Lanczos);
  if (childfopenMPI(sdt, "w", &fp) != 0) {
    return -1;
  }
  for (i = 0; i < i_max; i++) {
    fprintf(fp, " %ld %.10lf \n", i, creal(v0[i]));
  }
  fclose(fp);
  return 0;
}
