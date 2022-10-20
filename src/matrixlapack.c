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
/*-------------------------------------------------------------
 *[ver.2009.05.25]
 *-------------------------------------------------------------
 * Copyright (C) 2007-2009 Daisuke Tahara. All rights reserved.
 * Copyright (C) 2009-     Takahiro Misawa. All rights reserved.
 * Some functions are added by TM.
 *-------------------------------------------------------------*/
/*=================================================================================================*/

#include "matrixlapack.h"
#include <stdlib.h>
/**
 * @file   matrixlapack.c
 * @version 0.1, 0.2
 * @author Takahiro Misawa (The University of Tokyo)
 * 
 * @brief  wrapper for linear algebra operations using lapack
 * 
 * 
 */

#ifdef SR
int dsyevd_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *iwork, int *liwork, int *info);
int zheevd_(char *jobz, char *uplo, int *n, double complex *a, int *lda, double *w, double complex *work, int *lwork, double *rwork, int *iwork, int *liwork, int *info);
#else
int dsyev_(char *jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info);
int zheev_(char *jobz, char *uplo, int *n, double complex *a, int *lda, double *w, double complex *work, int *lwork, double *rwork, int *info);
#endif

int dsyevx_(char *jobz, char *range, char *uplo, int *n, double *a, int *lda, double *vl, double *vu, 
        int *il, int *iu, double *abstol, int *m, double *w, double *z__, int *ldz, 
        double *work, int *lwork, int *iwork, int *ifail, int *info);

/**
 *
 * @brief function for transforming Row-major matrix (C) to Column-major matrix (Fortran)
 * @param[in] N
 * @param[in] M
 * @param[in] A  Row-major matrix
 * @param[out] a  Column-major matrix
 * @author Takahiro Misawa (The University of Tokyo)
 */
void to_f(
  int N,  //!<[in]
  int M,   //!<[in]
  double** A,   //!<[in]
  double* a   //!<[out]
) {
  int i, j, k;
  k = 0;
  for (j = 0; j < M; j++) {
    for (i = 0; i < N; i++) {
      a[k] = A[i][j];
      k++;
    }
  }
}

/**
 *
 * @brief obtain eigenvalues of real symmetric A
 * @param[in] xNsize
 * @param[in] A  matrix
 * @param[out] r  eigenvalues
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int DSEVvalue(int xNsize, double** A, double* r) {

  int k;
  char jobz, uplo;
  int n, lda, lwork, info;
  double* a, * w, * work;
#ifdef SR
  int* iwork, liwork;
  liwork = 5 * xNsize + 3;
  iwork = (int*)malloc(liwork * sizeof(double));
#endif

  n = lda = xNsize;
  lwork = 4 * xNsize; /* 3*xNsize OK?*/

  a = (double*)malloc(xNsize * xNsize * sizeof(double));
  w = (double*)malloc(xNsize * sizeof(double));
  work = (double*)malloc(lwork * sizeof(double));

  to_f(xNsize, xNsize, A, a);

  jobz = 'N';
  uplo = 'U';

#ifdef SR
  M_DSYEV(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &iwork, &liwork, &info);
  free(iwork);
#else
  M_DSYEV(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
#endif

  if (info != 0) {
    free(a);
    free(w);
    free(work);
    return 0;
  }

  for (k = 0; k < xNsize; k++) {
    r[k] = w[k];
  }

  free(a);
  free(w);
  free(work);

  return 1;
}

// added by Misawa 20090309
// obtain eigen vectors
/**
 *
 *
 * @brief obtain eigenvalues and eigenvectors of real symmetric A
 * @param xNsize  size of matrix
 * @param A  matrix
 * @param r  eigenvalues
 * @param vec  eignevectos
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int DSEVvector(int xNsize, double** A, double* r, double** vec) {

  int i, j, k;
  char jobz, uplo;
  int n, lda, lwork, info;
  double* a, * w, * work;
#ifdef SR
  int* iwork, liwork;
  liwork = 5 * xNsize + 3;
  iwork = (int*)malloc(liwork * sizeof(double));
  lwork = 2 * xNsize * xNsize + 6 * xNsize + 1; /* 3*xNsize OK?*/
#else
  lwork = 4 * xNsize; /* 3*xNsize OK?*/
#endif

  n = lda = xNsize;

  a = (double*)malloc(xNsize * xNsize * sizeof(double));
  w = (double*)malloc(xNsize * sizeof(double));
  work = (double*)malloc(lwork * sizeof(double));

  k = 0;
  for (j = 0; j < xNsize; j++) {
    for (i = 0; i < xNsize; i++) {
      a[k] = A[i][j];
      k++;
    }
  }

  jobz = 'V';
  uplo = 'U';

#ifdef SR
  M_DSYEV(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info);
  free(iwork);
#else
  dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
#endif

  k = 0;
  for (i = 0; i < xNsize; i++) {
    for (j = 0; j < xNsize; j++) {
      vec[i][j] = a[k];
      k++;
    }
  }

  if (info != 0) {
    free(a);
    free(w);
    free(work);
    return 0;
  }

  for (k = 0; k < xNsize; k++) {
    r[k] = w[k];
  }

  free(a);
  free(w);
  free(work);

  return 1;
}

//added by Misawa 130121
//For complex Hermite matrix
/** 
 * 
 * @brief obtain eigenvalues and eigenvectors of Hermite matrix A
 * @param xNsize size of matrix
 * @param A matrix
 * @param r eigenvalues
 * @param vec eigenvectors
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int ZHEEVall(int xNsize, double complex **A, double *r, double complex **vec) {

  int i, j, k;
  char jobz, uplo;
  int n, lda, lwork, info, lrwork;
  double *rwork;
  double *w;
  double complex *a, *work;
#ifdef SR
  int *iwork, liwork;
  liwork = 5 * xNsize + 3;
  iwork = (int*)malloc(liwork * sizeof(double));
#endif

  n = lda = xNsize;
#ifdef SR
  lwork = xNsize * xNsize + 2 * xNsize; /* 3*xNsize OK?*/
  lrwork = 2 * xNsize*xNsize + 5 * xNsize + 1;
#else
  lwork = 4 * xNsize; /* 3*xNsize OK?*/
  lrwork = lwork;
#endif

  a = (double complex*)malloc(xNsize*xNsize * sizeof(double complex));
  w = (double*)malloc(xNsize * sizeof(double));
  work = (double complex*)malloc(lwork * sizeof(double complex));
  rwork = (double*)malloc(lrwork * sizeof(double));

  k = 0;
  for (j = 0; j < xNsize; j++) {
    for (i = 0; i < xNsize; i++) {
      a[k] = A[i][j];
      k++;
    }
  }

  jobz = 'V';
  uplo = 'U';

#ifdef SR
  int zheevd_(char *jobz, char *uplo, int *n, double complex *a, int *lda, double *w, double complex *work, int *lwork, double *rwork, int *iwork, int *liwork, int *info);
  free(iwork);
#else
  zheev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &info);
#endif

  if (info != 0) {
    free(a);
    free(w);
    free(work);
    free(rwork);
    return 0;
  }

  k = 0;
  for (i = 0; i < xNsize; i++) {
    for (j = 0; j < xNsize; j++) {
      vec[j][i] = a[k];
      k++;
    }
  }


  for (k = 0; k < xNsize; k++) {
    r[k] = w[k];
  }

  free(a);
  free(w);
  free(work);
  free(rwork);

  return 1;
}
