/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima */

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

/** 
 * 
 * 
 * @param dsfmt 
 * @param X 
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 *
 * @return 
 */
int FirstMultiply(dsfmt_t *dsfmt,struct BindStruct *X){
  
  long int i,i_max;      
  double dnorm;
  double Ns;
      
  Ns = 1.0*X->Def.Nsite;
  i_max=X->Check.idim_max;      

  if(X->Def.iInitialVecType==0){
    //For getting random numbers without any dependencies of threads,
    //we do not adopt omp in this part.
    for(i = 1; i <= i_max; i++){
      v0[i]=0.0;  
      v1[i]=2.0*(dsfmt_genrand_close_open(dsfmt)-0.5)+I*2.0*(dsfmt_genrand_close_open(dsfmt)-0.5);
    }
  }
  else{
    for(i = 1; i <= i_max; i++){
      v0[i]=0.0;
      v1[i]=2.0*(dsfmt_genrand_close_open(dsfmt)-0.5);
    }
  }
  
  dnorm=0.0;
#pragma omp parallel for default(none) private(i) shared(v1, i_max) reduction(+: dnorm) 
  for(i=1;i<=i_max;i++){
    dnorm += conj(v1[i])*v1[i];
  }
  dnorm=sqrt(dnorm);
  global_1st_norm = dnorm;
#pragma omp parallel for default(none) private(i) shared(v1) firstprivate(i_max, dnorm)
  for(i=1;i<=i_max;i++){
    v1[i] = v1[i]/dnorm;
  }

  TimeKeeperWithStep(X, cFileNameTimeKeep, cTPQStep, "a", step_i);
   
  mltply(X, v0, v1);
#pragma omp parallel for default(none) private(i) shared(v0, v1) firstprivate(i_max, Ns, LargeValue)
  for(i = 1; i <= i_max; i++){
    v0[i]=LargeValue*v1[i]-v0[i]/Ns;
  }
  dnorm=0.0;
#pragma omp parallel for default(none) private(i) shared(v0) firstprivate(i_max) reduction(+: dnorm)
  for(i=1;i<=i_max;i++){
    dnorm += conj(v0[i])*v0[i];
  }
  dnorm=sqrt(dnorm);
  global_norm = dnorm;
#pragma omp parallel for default(none) private(i) shared(v0) firstprivate(i_max, dnorm)
  for(i=1;i<=i_max;i++){
    v0[i] = v0[i]/dnorm;
  }

  TimeKeeperWithStep(X, cFileNameTimeKeep, cTPQStepEnd,"a", step_i);

  return 0;
}
