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

#ifdef _MAGMA
#include "matrixlapack_magma.h"


int diag_magma_cmp(int xNsize, double complex **A,
                   double complex *r,double complex **vec, int ngpu){

    if(ngpu <= 0) {
      fprintf(stdout, "  Error: ngpu is less than or equal to 0. ngpu = %d\n", ngpu);
      fprintf(stdout, "         exit in this function.\n\n");
      return -1;
    }

    magma_init();

    // local variable
    int i,j,k;
    magma_vec_t jobz;
    magma_uplo_t uplo;
    magma_int_t n, lda, lwork, lrwork;
    magma_int_t *iwork, liwork;
    double *rwork;
    double *w;
    magmaDoubleComplex *a, *work;
    magma_int_t info;

    {
      magma_int_t maxngpu = 0;
      magma_device_t *dummy;
      magma_getdevices(dummy, 9999, &maxngpu);
      if(maxngpu <= 0) {
          fprintf(stdout, "  Error: Your HPhi was built with a GPGPU library but NO GPU device is found.\n");
          fprintf(stdout, "         To perform fulldiag calculation without GPU, set \"ngpu 0\" in \"calcmod.def\".\n");
          exit(1);
      }
      if(ngpu > maxngpu) {
        fprintf(stdout, "  Warning: ngpu is beyond maxgpu. ngpu = %d, maxngpu =  %d\n", ngpu, maxngpu);
        fprintf(stdout, "           MAGMA library executes maxngpu.\n\n");
        ngpu = maxngpu;
      }
    }

    // auxiliary variable
    magmaDoubleComplex aux_work[1];
    double aux_rwork[1];
    magma_int_t aux_iwork[1];

    n = lda = xNsize;
    jobz = MagmaVec;
    uplo = MagmaUpper;

    // get optimal work size
    magma_zheevd(jobz, uplo,
                 n, NULL, lda, NULL,  // A, w
                 aux_work,  -1,
                 aux_rwork, -1,
                 aux_iwork, -1,
                 &info);

    lwork = (magma_int_t) MAGMA_Z_REAL(aux_work[0]);
    lrwork = (magma_int_t) aux_rwork[0];
    liwork = aux_iwork[0];

    // allocate
    magma_zmalloc_cpu(&a, xNsize * xNsize);
    magma_dmalloc_cpu(&w, xNsize);
    magma_zmalloc_pinned(&work, lwork);
    magma_dmalloc_cpu(&rwork, lrwork);
    magma_imalloc_cpu(&iwork, liwork);

    k=0;
    for(j=0;j<xNsize;j++){
        for(i=0;i<xNsize;i++){
            MAGMA_Z_REAL(a[k]) = creal(A[i][j]);
            MAGMA_Z_IMAG(a[k]) = cimag(A[i][j]);
            k++;
        }
    }

    // calculation eigen value and eigen vector
    magma_zheevd_m(ngpu, jobz, uplo,
                   n, a, lda, w,
                   work, lwork,
                   rwork, lrwork,
                   iwork, liwork,
                   &info);

    magma_free_cpu(iwork);
    if(info != 0){
        magma_free_cpu(a);
        magma_free_cpu(w);
        magma_free_pinned(work);
        magma_free_cpu(rwork);
        return -1;
    }

    k=0;
    for(i=0;i<xNsize;i++){
        for(j=0;j<xNsize;j++){
            vec[i][j] = MAGMA_Z_REAL(a[k]) + MAGMA_Z_IMAG(a[k]) * I;
            k++;
        }
    }

    for(k=0;k<xNsize;k++){
        r[k] = w[k];
    }

    magma_free_cpu(a);
    magma_free_cpu(w);
    magma_free_pinned(work);
    magma_free_cpu(rwork);

    magma_finalize();

    return 0;
}

int diag_magma_real(int xNsize, double **A,
                    double *r, double **vec, int ngpu){

    if(ngpu <= 0) {
      fprintf(stdout, "  Error: ngpu is less than or equal to 0. ngpu = %d\n", ngpu);
      fprintf(stdout, "         exit in this function.\n\n");
      return -1;
    }

    magma_init();

    // local variable
    int i,j,k;
    magma_vec_t jobz;
    magma_uplo_t uplo;
    magma_int_t n, lda, lwork;
    magma_int_t *iwork, liwork;
    double *w;
    double *a, *work;
    magma_int_t info;

    {
      magma_int_t maxngpu = 0;
      magma_device_t *dummy;
      magma_getdevices(dummy, 9999, &maxngpu);
      if(ngpu > maxngpu) {
        fprintf(stdout, "  Warning: ngpu is beyond maxgpu. ngpu = %d, maxngpu =  %d\n", ngpu, maxngpu);
        fprintf(stdout, "           MAGMA library executes maxngpu.\n\n");
        ngpu = maxngpu;
      }
    }

    // auxiliary variable
    double aux_work[1];
    magma_int_t aux_iwork[1];

    n = lda = xNsize;
    jobz = MagmaVec;
    uplo = MagmaUpper;

    // get optimal work size
    magma_dsyevd(jobz, uplo,
                 n, NULL, lda, NULL,  // A, w
                 aux_work,  -1,
                 aux_iwork, -1,
                 &info);

    lwork = (magma_int_t) MAGMA_D_REAL(aux_work[0]);
    liwork = aux_iwork[0];

    // allocate
    magma_dmalloc_cpu(&a, xNsize * xNsize);
    magma_dmalloc_cpu(&w, xNsize);
    magma_dmalloc_pinned(&work, lwork);
    magma_imalloc_cpu(&iwork, liwork);

    k=0;
    for(j=0;j<xNsize;j++){
        for(i=0;i<xNsize;i++){
            a[k] = A[i][j];
            k++;
        }
    }

    // calculation eigen value and eigen vector
    magma_dsyevd_m(ngpu, jobz, uplo,
                   n, a, lda, w,
                   work, lwork,
                   iwork, liwork,
                   &info);

    magma_free_cpu(iwork);
    if(info != 0){
        magma_free_cpu(a);
        magma_free_cpu(w);
        magma_free_pinned(work);
        return -1;
    }

    k=0;
    for(i=0;i<xNsize;i++){
        for(j=0;j<xNsize;j++){
            vec[i][j]=a[k];
            k++;
        }
    }

    for(k=0;k<xNsize;k++){
        r[k] = w[k];
    }

    magma_free_cpu(a);
    magma_free_cpu(w);
    magma_free_pinned(work);

    magma_finalize();

    return 0;
}
#endif
