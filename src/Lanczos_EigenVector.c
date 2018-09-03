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

#include "Common.h"
#include "mltply.h"
#include "CalcTime.h"
#include "Lanczos_EigenVector.h"
#include "wrapperMPI.h"

/**
 *
 * @file   Lanczos_EigenVector.c
 * @version 0.1, 0.2
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 * 
 * @brief  Calculate eigenvectors by the Lanczos method.
 * 
 */

/** 
 * @brief Calculate eigenvectors by the Lanczos method. \n
 * The calculated tridiagonal matrix components @f$ \alpha_i, \beta_i@f$ are stored in each array @f$ \verb|alpha| @f$ and @f$\verb|beta|@f$\n
 * (@f$ i = 0\cdots N_c@f$, where @f$ N_c@f$ is the step where the calculated energy satisfies the convergence condition).
 * 
 * @param X [in,out] Struct for getting information to calculate eigenvectors.
 * @version 0.2
 * @details add an option to choose a type of initial vectors from complex or real types. 
 * @version 0.1
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
