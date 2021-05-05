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
 * @file   FirstMultiply.c
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @version 0.1
 * @brief  Multiplication @f$ v_0 = H v_1 @f$ at the first step for TPQ mode (@f$ v_1 @f$ is the random or inputted vector).
 *
 */

///
/// \brief Multiplication @f$ v_0 = H v_1 @f$ at the first step for TPQ mode (@f$ v_1 @f$ is the random or inputted vector).
/// \param rand_i [in] A rundom number seed for giving the initial vector @f$ v_1 @f$.
/// \param X [in] Struct to get information of the vector @f$ v_1 @f$ for the first step calculation.
/// \retval -1 fail the multiplication @f$ v_0 = H v_1 @f$.
/// \retval 0 succeed the multiplication @f$ v_0 = H v_1 @f$.
/// \version 0.1
/// \author Takahiro Misawa (The University of Tokyo)
/// \author Kazuyoshi Yoshimi (The University of Tokyo)
int MakeIniVec(int rand_i, struct BindStruct *X) {

  long int i, i_max;
  double complex dnorm;
  double Ns;
  long unsigned int u_long_i;
  dsfmt_t dsfmt;
  int mythread;

  Ns = 1.0*X->Def.NsiteMPI;
  i_max = X->Check.idim_max;

#pragma omp parallel default(none) private(i, mythread, u_long_i, dsfmt) \
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
    u_long_i = 123432 + (rand_i + 1)*labs(X->Def.initial_iv) + mythread + nthreads * myrank;
    dsfmt_init_gen_rand(&dsfmt, u_long_i);

    if (X->Def.iInitialVecType == 0) {
  
    StartTimer(3101);
#pragma omp for
      for (i = 1; i <= i_max; i++)
        v1[i] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5) + 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5)*I;
    }/*if (X->Def.iInitialVecType == 0)*/
    else {
#pragma omp for
      for (i = 1; i <= i_max; i++)
          v1[i] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5);
    }
    StopTimer(3101);

  }/*#pragma omp parallel*/
  /*
    Normalize v
  */
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
