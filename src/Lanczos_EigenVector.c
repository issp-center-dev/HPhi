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

#include "Common.h"
#include "Lanczos_EigenVector.h"
#include "wrapperMPI.h"

/** 
 * 
 * 
 * @param X 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
void Lanczos_EigenVector(struct BindStruct *X){

  fprintf(stdoutMPI, "%s", cLogLanczos_EigenVectorStart);
  
  int i,j,i_max,iv;  	 
  int k_exct;
  double beta1,alpha1,dnorm, dnorm_inv;
  double complex temp1,temp2;

// for GC
  long unsigned int u_long_i;
  dsfmt_t dsfmt;

  k_exct = X->Def.k_exct;
	
  iv=X->Large.iv;
  i_max=X->Check.idim_max;
 
  //Eigenvectors by Lanczos method
  //initialization: initialization should be identical to that of Lanczos_EigenValue.c
#pragma omp parallel for default(none) private(i) shared(v0, v1, vg) firstprivate(i_max)
  for(i=1;i<=i_max;i++){
    v0[i]=0.0+0.0*I;
    v1[i]=0.0+0.0*I;
    vg[i]=0.0+0.0*I;
  }
    
  if(initial_mode == 0){
    v1[iv]=1.0;
    vg[iv]=vec[k_exct][1];
  }else if(initial_mode==1){      
    iv = X->Def.initial_iv;
    u_long_i = 123432 + abs(iv);
    dsfmt_init_gen_rand(&dsfmt, u_long_i);    
    for(i = 1; i <= i_max; i++){
      v1[i]=2.0*(dsfmt_genrand_close_open(&dsfmt)-0.5)+2.0*(dsfmt_genrand_close_open(&dsfmt)-0.5)*I;
    }
    dnorm=0;
    #pragma omp parallel for default(none) private(i) shared(v1, i_max) reduction(+: dnorm) 
    for(i=1;i<=i_max;i++){
      dnorm += conj(v1[i])*v1[i];
    }    
    dnorm=sqrt(dnorm);
    dnorm_inv=1.0/dnorm;
#pragma omp parallel for default(none) private(i) shared(v1,vg,vec,k_exct) firstprivate(i_max, dnorm_inv)
    for(i=1;i<=i_max;i++){
      v1[i]        = v1[i]*dnorm_inv;
      vg[i]        = v1[i]*vec[k_exct][1];
    }
  }
  
  mltply(X, v0, v1);
  
  alpha1=alpha[1];
  beta1=beta[1];

#pragma omp parallel for default(none) private(j) shared(vec, v0, v1, vg) firstprivate(alpha1, beta1, i_max, k_exct)
  for(j=1;j<=i_max;j++){
    vg[j]+=vec[k_exct][2]*(v0[j]-alpha1*v1[j])/beta1;
  }
    
  //iteration
  for(i=2;i<=X->Large.itr-1;i++){
#pragma omp parallel for default(none) private(j, temp1, temp2) shared(v0, v1) firstprivate(i_max, alpha1, beta1)
    for(j=1;j<=i_max;j++){
      temp1=v1[j];
      temp2=(v0[j]-alpha1*v1[j])/beta1;
      v0[j]=-beta1*temp1;
      v1[j]=temp2;        
    }
    mltply(X, v0, v1);   
	
    alpha1 = alpha[i];
    beta1  = beta[i];

#pragma omp parallel for default(none) private(j) shared(vec, v0, v1, vg) firstprivate(alpha1, beta1, i_max, k_exct, i)
    for(j=1;j<=i_max;j++){
      vg[j] += vec[k_exct][i+1]*(v0[j]-alpha1*v1[j])/beta1;
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
  dnorm=sqrt(dnorm);
  dnorm_inv=dnorm;
#pragma omp parallel for default(none) private(j) shared(v0) firstprivate(i_max, dnorm_inv)
  for(j=1;j<=i_max;j++){
    v0[j] = v0[j]*dnorm_inv;
  }
  
  TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenVectorFinish, "a");
  fprintf(stdoutMPI, "%s", cLogLanczos_EigenVectorEnd);
}
