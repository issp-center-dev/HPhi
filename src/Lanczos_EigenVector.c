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
 * @file Lanczos_EigenVector.c
 *
 * @brief Reconstruct eigenvector from Lanczos tridiagonal matrix
 *
 * After Lanczos_EigenValue() finds eigenvalues and stores the tridiagonal
 * matrix elements (alpha[], beta[]), this function reconstructs the
 * corresponding eigenvector.
 *
 * Algorithm:
 * The eigenvector |psi> in the original basis is:
 *   |psi> = sum_j c_j |q_j>
 * where |q_j> are Lanczos vectors and c_j are components of the
 * eigenvector of the tridiagonal matrix T.
 *
 * Since storing all Lanczos vectors |q_j> would require O(N*M) memory
 * (N = Hilbert space dim, M = Lanczos steps), we regenerate them:
 *
 * 1. Diagonalize tridiagonal matrix T to get eigenvector components c_j
 * 2. Re-run Lanczos iteration with same initial vector
 * 3. At each step j, accumulate: |psi> += c_j * |q_j>
 *
 * This requires O(N) memory but O(M) H*v multiplications.
 *
 * Initial vector:
 * - initial_mode = 0: Deterministic (one element = 1)
 * - initial_mode = 1: Random (for grand canonical ensembles)
 *
 * @version 0.2 Added real/complex initial vector option
 * @version 0.1
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
#include "Common.h"
#include "mltply.h"
#include "CalcTime.h"
#include "Lanczos_EigenVector.h"
#include "wrapperMPI.h"

/**
 * @brief Reconstruct eigenvector by regenerating Lanczos vectors
 *
 * Computes the eigenvector corresponding to the k_exct-th eigenvalue
 * found by Lanczos_EigenValue(). The eigenvector is stored in v0[].
 *
 * Process:
 * 1. Diagonalize stored tridiagonal matrix (alpha[], beta[])
 * 2. Extract eigenvector components vec_j for target eigenvalue
 * 3. Initialize Lanczos iteration with same initial vector
 * 4. For each step j = 0 to Lanczos_restart:
 *    - vg[] += vec_j * (current Lanczos vector)
 *    - Perform Lanczos step to get next vector
 * 5. Normalize result: v0[] = vg[] / |vg|
 *
 * Memory usage:
 * - v0, v1: Lanczos vectors (size idim_max each)
 * - vg: Accumulated eigenvector (size idim_max)
 * - vec: Tridiagonal eigenvector components (size Lanczos_restart)
 *
 * @param X Struct with Lanczos parameters and tridiagonal matrix [in,out]
 *          Uses: alpha[], beta[], Lanczos_restart, k_exct
 *          Sets: v0[] (final eigenvector)
 *
 * @version 0.2 Added real/complex initial vector option
 * @version 0.1
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void Lanczos_EigenVector(struct BindStruct *X){

  fprintf(stdoutMPI, "%s", cLogLanczos_EigenVectorStart);
  TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenVectorStart, "a");

  long int i,j,i_max,iv;
  int k_exct, iproc;
  double beta1,alpha1,dnorm, dnorm_inv;
  double complex temp1,temp2,cdnorm;
  int mythread;

// for GC
  long unsigned int u_long_i, sum_i_max, i_max_tmp;
  dsfmt_t dsfmt;

  k_exct = X->Def.k_exct;

  iv=X->Large.iv;
  i_max=X->Check.idim_max;

  if(initial_mode == 0){

    sum_i_max = SumMPI_li(X->Check.idim_max);
    X->Large.iv = (sum_i_max / 2 + X->Def.initial_iv) % sum_i_max + 1;
    iv=X->Large.iv;
#pragma omp parallel for default(none) private(i) shared(v0, v1,vg) firstprivate(i_max)
    for(i = 1; i <= i_max; i++){
      v0[i]=0.0;
      v1[i]=0.0;
      vg[i]=0.0;
    }

    sum_i_max = 0;
    for (iproc = 0; iproc < nproc; iproc++) {

      i_max_tmp = BcastMPI_li(iproc, i_max);
      if (sum_i_max <= iv && iv < sum_i_max + i_max_tmp) {

        if (myrank == iproc) {
          v1[iv - sum_i_max+1] = 1.0;
          if (X->Def.iInitialVecType == 0) {
            v1[iv - sum_i_max+1] += 1.0*I;
            v1[iv - sum_i_max+1] /= sqrt(2.0);
          }
          vg[iv - sum_i_max+1]=conj(vec[k_exct][1])*v1[iv - sum_i_max+1];

        }/*if (myrank == iproc)*/
      }/*if (sum_i_max <= iv && iv < sum_i_max + i_max_tmp)*/

      sum_i_max += i_max_tmp;

    }/*for (iproc = 0; iproc < nproc; iproc++)*/

  }/*if(initial_mode == 0)*/
  else if(initial_mode==1){
    iv = X->Def.initial_iv;
    //fprintf(stdoutMPI, "  initial_mode=%d (random): iv = %ld i_max=%ld k_exct =%d \n",initial_mode,iv,i_max,k_exct);
    #pragma omp parallel default(none) private(i, u_long_i, mythread, dsfmt) \
            shared(v0, v1, iv, X, nthreads, myrank) firstprivate(i_max)
    {

#pragma omp for
      for (i = 1; i <= i_max; i++) {
        v0[i] = 0.0;
      }
      /*
       Initialize MT
      */
#ifdef _OPENMP
      mythread = omp_get_thread_num();
#else
      mythread = 0;
#endif
      u_long_i = 123432 + labs(iv) + mythread + nthreads * myrank;
      dsfmt_init_gen_rand(&dsfmt, u_long_i);

      if (X->Def.iInitialVecType == 0) {
#pragma omp for
        for (i = 1; i <= i_max; i++)
          v1[i] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5) + 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5)*I;
      }
      else {
#pragma omp for
        for (i = 1; i <= i_max; i++)
          v1[i] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5);
      }
    }/*#pragma omp parallel*/
    /*
     Normalize
    */
    cdnorm=0.0;
#pragma omp parallel for default(none) private(i) shared(v1, i_max) reduction(+: cdnorm) 
    for(i=1;i<=i_max;i++){
     cdnorm += conj(v1[i])*v1[i];
    }
    cdnorm = SumMPI_dc(cdnorm);
    dnorm=creal(cdnorm);
    dnorm=sqrt(dnorm);
#pragma omp parallel for default(none) private(i) shared(v1, vec, vg) firstprivate(i_max, dnorm, k_exct)
    for(i=1;i<=i_max;i++){
      v1[i] = v1[i]/dnorm;
      vg[i] = v1[i]*conj(vec[k_exct][1]);
    }
  }/*else if(initial_mode==1)*/
  StartTimer(4201);
  mltply(X, v0, v1);
  StopTimer(4201);

  alpha1=alpha[1];
  beta1=beta[1];

#pragma omp parallel for default(none) private(j) shared(vec, v0, v1, vg) firstprivate(alpha1, beta1, i_max, k_exct)
  for(j=1;j<=i_max;j++){
    vg[j]+=conj(vec[k_exct][2])*(v0[j]-alpha1*v1[j])/beta1;
  }

  //iteration
  for(i=2;i<=X->Large.itr-1;i++) {
      /*
    if (abs(beta[i]) < pow(10.0, -15)) {
      break;
    }
*/
#pragma omp parallel for default(none) private(j, temp1, temp2) shared(v0, v1) firstprivate(i_max, alpha1, beta1)
    for (j = 1; j <= i_max; j++) {
      temp1 = v1[j];
      temp2 = (v0[j] - alpha1 * v1[j]) / beta1;
      v0[j] = -beta1 * temp1;
      v1[j] = temp2;
    }
    StartTimer(4201);
    mltply(X, v0, v1);
    StopTimer(4201);
    alpha1 = alpha[i];
    beta1 = beta[i];
#pragma omp parallel for default(none) private(j) shared(vec, v0, v1, vg) firstprivate(alpha1, beta1, i_max, k_exct, i)
    for (j = 1; j <= i_max; j++) {
      vg[j] += conj(vec[k_exct][i + 1]) * (v0[j] - alpha1 * v1[j]) / beta1;
    }
  }

#pragma omp parallel for default(none) private(j) shared(v0, vg) firstprivate(i_max)
    for(j=1;j<=i_max;j++){
      v0[j] = vg[j];
    } 
      
  //normalization
  dnorm=0.0;
#pragma omp parallel for default(none) reduction(+:dnorm) private(j) shared(v0) firstprivate(i_max)
  for(j=1;j<=i_max;j++){
    dnorm += conj(v0[j])*v0[j];
  }
  dnorm = SumMPI_d(dnorm);
  dnorm=sqrt(dnorm);
  dnorm_inv=1.0/dnorm;
#pragma omp parallel for default(none) private(j) shared(v0) firstprivate(i_max, dnorm_inv)
  for(j=1;j<=i_max;j++){
    v0[j] = v0[j]*dnorm_inv;
  }
  
  TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenVectorFinish, "a");
  fprintf(stdoutMPI, "%s", cLogLanczos_EigenVectorEnd);
}
