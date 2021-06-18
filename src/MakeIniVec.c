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
#include "FirstMultiply.h"
#include "expec_energy_flct.h"
#include "common/setmemory.h"
#include "wrapperMPI.h"
#include "CalcTime.h"

/**
 * @file   MakeIniVec.c
 * @author Takahiro Misawa (BAQIS)
 * @version 3.4
 * @brief  Generating the random initial vectors for TPQ  mode (@f$ v_1=v_0 @f$ is the random or inputted vector).
 *
 */

///
/// \brief  Generating the random initial vectors for TPQ  mode (@f$ v_1=v_0 @f$ is the random or inputted vector).
/// \param rand_i [in] A random number seed for giving the initial vector @f$ v_1=v_0 @f$.
/// \param X [in] Struct to get necessary information.
/// \retval -1 fail generating the random initial vectors.
/// \retval 0 succeed in generating the random initial vectors.
/// \version 3.4
/// \author Takahiro Misawa (BAQIS)
/*Note: X->Def.iInitialVecType == 0: All components are given by random complex numbers x+i*y, x = [-1,1),y=[-1,1]*/
/*Note: X->Def.iInitialVecType ==-1: random vectors on the 2*i_max complex sphere*/
/*Note: X->Def.iInitialVecType == others: All components are given by random real numbers x, x=[-1,1)*/
int MakeIniVec(int rand_i, struct BindStruct *X) {

  long int i, i_max;
  double complex dnorm;
  double Ns;
  long unsigned int u_long_i;
  dsfmt_t dsfmt;
  int mythread;
  double rand_X,rand_Y;
  double complex rand_Z1,rand_Z2;

  Ns = 1.0*X->Def.NsiteMPI;
  i_max = X->Check.idim_max;

#pragma omp parallel default(none) private(i, mythread, u_long_i, dsfmt,rand_X,rand_Y,rand_Z1,rand_Z2) \
        shared(v0, v1, nthreads, myrank, rand_i, X, stdoutMPI, cLogCheckInitComplex, cLogCheckInitReal) \
        firstprivate(i_max)
  {
#pragma omp for
    for (i = 1; i <= i_max; i++) {
      v0[i] = 0.0;
      v1[i] = 0.0;
    }
    /*
    Initialise MT
    */
#ifdef _OPENMP
    mythread = omp_get_thread_num();
#else
    mythread = 0;
#endif
    u_long_i = 123432 + (rand_i+1)*labs(X->Def.initial_iv) + mythread + nthreads * myrank;
    dsfmt_init_gen_rand(&dsfmt, u_long_i);

    if (X->Def.iInitialVecType == 0) {
      StartTimer(3101);
      #pragma omp for
      for (i = 1; i <= i_max; i++)
        v1[i] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5) + 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5)*I;
      /*if (X->Def.iInitialVecType == 0)*/
    }else if (X->Def.iInitialVecType == -1) {
      StartTimer(3101);
      #pragma omp for 
      for (i = 1; i <= i_max; i++){
        rand_X   = dsfmt_genrand_close_open(&dsfmt);
        rand_Y   = dsfmt_genrand_close_open(&dsfmt);
        rand_Z1  = sqrt(-2.0*log(rand_X))*cos(2.0*M_PI*rand_Y);
        rand_Z2  = sqrt(-2.0*log(rand_X))*sin(2.0*M_PI*rand_Y);
        v1[i]    = rand_Z1+I*rand_Z2;
      } 
    /*if (X->Def.iInitialVecType == -1)*/
    }else {
      #pragma omp for
      for (i = 1; i <= i_max; i++)
         v1[i] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5);
    }
    StopTimer(3101);
  }/*#pragma omp parallel*/
  /*Normalize v*/
  dnorm=0.0;
#pragma omp parallel for default(none) private(i) shared(v1, i_max) reduction(+: dnorm) 
  for(i=1;i<=i_max;i++){
    dnorm += conj(v1[i])*v1[i];
  }
  dnorm = SumMPI_dc(dnorm);
  dnorm=sqrt(dnorm);
  global_1st_norm = dnorm;
#pragma omp parallel for default(none) private(i) shared(v0,v1) firstprivate(i_max, dnorm)
  for(i=1;i<=i_max;i++){
    v1[i] = v1[i]/dnorm;
    v0[i] = v1[i];
  }
  
  TimeKeeperWithRandAndStep(X, cFileNameTimeKeep, cTPQStep, "a", rand_i, step_i);
   
  return 0;
}
