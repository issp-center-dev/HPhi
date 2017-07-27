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
#include "Multiply.h"
#include "wrapperMPI.h"
#include "mltply.h"

/** 
 * @brief  Function of multiplying Hamiltonian for TPQ calculation
 * 
 * @param X data list for calculation
 * @retval 0  normally finished
 * @retval -1 unnormally finished
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int Multiply
(
 struct BindStruct *X
 )
{
  
  long int i,i_max;
  double complex dnorm;
  double Ns;

  i_max=X->Check.idim_max;      
  Ns = 1.0*X->Def.NsiteMPI;
 // mltply is in expec_energy.c v0=H*v1
  dnorm=0.0;
#pragma omp parallel for default(none) reduction(+: dnorm) private(i) shared(v0, v1) firstprivate(i_max, Ns, LargeValue)
  for(i = 1; i <= i_max; i++){
    v0[i]=LargeValue*v1[i]-v0[i]/Ns;  //v0=(l-H/Ns)*v1
    dnorm += conj(v0[i])*v0[i];
  }
  dnorm=SumMPI_dc(dnorm);
  dnorm=sqrt(dnorm);
  global_norm = dnorm;
#pragma omp parallel for default(none) private(i) shared(v0) firstprivate(i_max, dnorm)
  for(i=1;i<=i_max;i++){
    v0[i] = v0[i]/dnorm;
  }
  return 0;
}

/**
 * @brief  Function of multiplying Hamiltonian for Time Evolution
 *
 * @param X data list for calculation
 *
 * @retval 0  normally finished
 * @retval -1 unnormally finished
 */
int MultiplyForTEM
        (
                struct BindStruct *X
        )
{

  long int i,i_max;
  int coef;
  double complex dnorm;
  double  complex tmp1 = 1.0,tmp2;
  double dt=X->Def.Param.TimeSlice;

  //Make |v0> = |psi(t+dt)> from |v1> = |psi(t)> and |v0> = H |psi(t)>
  i_max=X->Check.idim_max;
  // mltply is in expec_energy.c v0=H*v1
  tmp1 *= -I*dt;
#pragma omp parallel for default(none) reduction(+: dnorm) private(i) shared(v0, v1, v2) firstprivate(i_max, dt, tmp1,tmp2)
  for(i = 1; i <= i_max; i++){
    tmp2 = v0[i];
    v0[i]=v1[i]+tmp1*tmp2;  //v0=(1-i*dt*H)*v1
    v1[i]=tmp2;
    v2[i]=0.0 + I*0.0;
  }

  for(coef = 2; coef <= X->Def.Param.ExpandCoef; coef++){
    tmp1 *= -I*dt/(double complex)coef;
    //v2 = H*v1 = H^coef |psi(t)>
    mltply(X, v2, v1);
    //[TODO] mltply diagonal term

#pragma omp parallel for default(none) private(i) shared(v0, v1, v2) firstprivate(i_max, tmp1, myrank)
    for(i = 1; i <= i_max; i++){
      v0[i] += tmp1*v2[i];
      v1[i]=v2[i];
      v2[i]=0.0 + I*0.0;
    }
  }

  dnorm=0.0;
#pragma omp parallel for default(none) reduction(+: dnorm) private(i) shared(v0) firstprivate(i_max, dt)
  for(i = 1; i <= i_max; i++){
    dnorm += conj(v0[i])*v0[i];
  }
  dnorm=SumMPI_dc(dnorm);
  dnorm=sqrt(dnorm);
  global_norm = dnorm;
#pragma omp parallel for default(none) private(i) shared(v0) firstprivate(i_max, dnorm)
  for(i=1;i<=i_max;i++){
    v0[i] = v0[i]/dnorm;
  }
  return 0;
}
