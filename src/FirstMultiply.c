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
#include "MakeIniVec.h"
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
int FirstMultiply(int rand_i, struct BindStruct *X) {

  long int i, i_max;
  double complex dnorm;
  double Ns;
  long unsigned int u_long_i;
  dsfmt_t dsfmt;
  int mythread;

  Ns = 1.0*X->Def.NsiteMPI;
  i_max = X->Check.idim_max;

  /**@brief
  Initialize v1 and v0 = v1
  */
  MakeIniVec(rand_i,X);
  
  TimeKeeperWithRandAndStep(X, cFileNameTimeKeep, cTPQStep, "a", rand_i, step_i);
   
  StartTimer(3102);
  if(expec_energy_flct(X) !=0){ //v1 <- v0 and v0 = H*v1
    StopTimer(3102);
    return -1;
  }
  StopTimer(3102);
#pragma omp parallel for default(none) private(i) shared(v0, v1, list_1) firstprivate(i_max, Ns, LargeValue, myrank)
  for(i = 1; i <= i_max; i++){
    v0[i]=LargeValue*v1[i]-v0[i]/Ns;
  }

  dnorm=0.0;
#pragma omp parallel for default(none) private(i) shared(v0) firstprivate(i_max) reduction(+: dnorm)
  for(i=1;i<=i_max;i++){
    dnorm += conj(v0[i])*v0[i];
  }
  dnorm = SumMPI_dc(dnorm);
  dnorm=sqrt(dnorm);
  global_norm = dnorm;
#pragma omp parallel for default(none) private(i) shared(v0) firstprivate(i_max, dnorm)
  for(i=1;i<=i_max;i++){
    v0[i] = v0[i]/dnorm;
  }
  TimeKeeperWithRandAndStep(X, cFileNameTimeKeep, cTPQStepEnd, "a", rand_i, step_i);
  return 0;
}
