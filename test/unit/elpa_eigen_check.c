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

/* Mathematical verification of the ELPA FullDiag path (design doc section 6):
   residual ||A z - w z||, orthogonality ||Z^H Z - I||, and eigenvalue
   agreement with LAPACK zheev, all checked against the scale-rule thresholds
   c*N*eps*||A|| (residual, eigenvalues) and c*N*eps (orthogonality), c = 50.
   Run with any rank count (a 2D BLACS grid is built from MPI_Dims_create). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <mpi.h>
#include "matrixscalapack.h"
#include "matrixlapack_elpa.h"

#define NDIM 97
#define TOL_C 50.0

extern void zheev_(char *jobz, char *uplo, int *n, double complex *a,
                   int *lda, double *w, double complex *work, int *lwork,
                   double *rwork, int *info);

static double complex MatElem(int i, int j) {
  double re = 1.0 / (1.0 + fabs((double)(i - j)));
  double im = (double)(i - j) / (double)(NDIM * NDIM);
  return re + im * I;   /* MatElem(j,i) == conj(MatElem(i,j)) by construction */
}

int main(int argc, char **argv) {
  /* Single int zero/negative-one, matching the verified idiom in
     diag_scalapack_cmp (src/matrixscalapack.c) and lapack_diag_elpa
     (src/lapack_diag.c, Task 5): all of blacs_get_'s "request",
     numroc_'s "isrcproc", and descinit_'s "irsrc"/"icsrc" take
     `const int *`, not `const long int *`. */
  int i_negone = -1, i_zero = 0;
  int rank, size, ictxt, iam, nprocs, info, ok = 1;
  int nprow, npcol, myrow, mycol;
  long int n = NDIM, mb, mp, nq, i, j, k;
  int lld, dims[2] = {0, 0};
  int descA[9], descZ[9];
  double complex *A_distr, *Z_distr, *vecs, *vec_tmp;
  double *w, anorm = 0.0, eps = DBL_EPSILON;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size, 2, dims);
  nprow = dims[0]; npcol = dims[1];
  /* Cap the block size so every process row/column owns a block
     (e.g. np=3 -> 3x1 grid needs nblk <= 32 for N=97). */
  mb = ElpaBlockSize(n, nprow, npcol);
  blacs_pinfo_(&iam, &nprocs);
  blacs_get_(&i_negone, &i_zero, &ictxt);
  blacs_gridinit_(&ictxt, "R", &nprow, &npcol);
  blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);

  mp = numroc_(&n, &mb, &myrow, &i_zero, &nprow);
  nq = numroc_(&n, &mb, &mycol, &i_zero, &npcol);
  lld = (mp > 0) ? (int)mp : 1;
  descinit_(descA, &n, &n, &mb, &mb, &i_zero, &i_zero, &ictxt, &lld, &info);
  descinit_(descZ, &n, &n, &mb, &mb, &i_zero, &i_zero, &ictxt, &lld, &info);

  A_distr = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));
  Z_distr = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));
  w = malloc(n * sizeof(double));
  vec_tmp = malloc(n * sizeof(double complex));
  vecs = malloc(n * n * sizeof(double complex)); /* rank 0: gathered Z */

  for (i = 0; i < n; i++) {
    double colsum = 0.0;
    for (j = 0; j < n; j++) {
      DivMat(i, j, MatElem((int)i, (int)j), A_distr, descA);
      colsum += cabs(MatElem((int)i, (int)j));
    }
    if (colsum > anorm) anorm = colsum;  /* max abs row sum; == ||A||_1 since A is Hermitian */
  }

  if (diag_elpa_cmp((int)n, A_distr, Z_distr, w,
                    (int)mp, (int)nq, myrow, mycol, (int)mb, 0) != 0) {
    if (rank == 0) fprintf(stderr, "diag_elpa_cmp failed\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  for (k = 0; k < n; k++) {
    GetEigenVectorBlock(k, n, Z_distr, descZ, vec_tmp);
    if (rank == 0) for (i = 0; i < n; i++) vecs[k * n + i] = vec_tmp[i];
  }
  FreeEigenVectorGatherContext();

  if (rank == 0) {
    /* (a) eigenvalues vs LAPACK zheev */
    double complex *a_full = malloc(n * n * sizeof(double complex));
    double *w_ref = malloc(n * sizeof(double));
    double *rwork = malloc((3 * n - 2) * sizeof(double));
    int lwork = 4 * NDIM, n_int = NDIM;
    double complex *work = malloc(lwork * sizeof(double complex));
    for (j = 0; j < n; j++) for (i = 0; i < n; i++)
      a_full[j * n + i] = MatElem((int)i, (int)j);
    zheev_("N", "U", &n_int, a_full, &n_int, w_ref, work, &lwork, rwork, &info);
    for (k = 0; k < n; k++) {
      if (fabs(w[k] - w_ref[k]) > TOL_C * n * eps * anorm) {
        fprintf(stderr, "eigenvalue %ld mismatch: %e vs %e\n", k, w[k], w_ref[k]);
        ok = 0;
      }
    }
    /* (b) residual ||A z_k - w_k z_k||_inf */
    for (k = 0; k < n && ok; k++) {
      for (i = 0; i < n; i++) {
        double complex r = -w[k] * vecs[k * n + i];
        for (j = 0; j < n; j++) r += MatElem((int)i, (int)j) * vecs[k * n + j];
        if (cabs(r) > TOL_C * n * eps * anorm) {
          fprintf(stderr, "residual too large: state %ld row %ld: %e\n", k, i, cabs(r));
          ok = 0; break;
        }
      }
    }
    /* (c) orthogonality |z_k^H z_l - delta_kl| */
    for (k = 0; k < n && ok; k++) {
      for (j = k; j < n; j++) {
        double complex dot = 0.0;
        for (i = 0; i < n; i++) dot += conj(vecs[k * n + i]) * vecs[j * n + i];
        if (cabs(dot - (k == j ? 1.0 : 0.0)) > TOL_C * n * eps) {
          fprintf(stderr, "orthogonality violated: (%ld,%ld) %e\n", k, j, cabs(dot));
          ok = 0; break;
        }
      }
    }
    printf("elpa_eigen_check: %s\n", ok ? "OK" : "FAILED");
    free(a_full); free(w_ref); free(rwork); free(work);
  }
  MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  free(A_distr); free(Z_distr); free(w); free(vec_tmp); free(vecs);
  MPI_Finalize();
  return ok ? 0 : 1;
}
