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

#include "matrixlapack.h"
#include "Common.h"
#include "wrapperMPI.h"
#include "mfmemory.h"
#include "xsetmem.h"

/** 
 * 
 * 
 * @param alpha 
 * @param beta 
 * @param ndim 
 * @param E 
 * @param X 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Youhei Yamaji (The University of Tokyo)
 */
void vec12(double alpha[],double beta[],int ndim,
	   double E[],struct BindStruct *X){
  
  int j,k,nvec;

  double **tmpA, **tmpvec;
  double *tmpr;
    
  nvec = X->Def.nvec;

  d_malloc2(tmpA,ndim,ndim); 
  d_malloc2(tmpvec,ndim,ndim); 
  d_malloc1(tmpr,ndim);

#pragma omp parallel for default(none) firstprivate(ndim) private(j,k) shared(tmpA)
  for(k=0;k<=ndim-1;k++){
    for(j=0;j<=ndim-1;j++){
      tmpA[k][j]=0.0;
    }
  }
#pragma omp parallel for default(none) firstprivate(ndim, nvec) private(j,k) shared(vec)
  for(k=1;k<=nvec;k++){
    for(j=1;j<=ndim;j++){
      vec[k][j]=0.0;
    }
  }

#pragma omp parallel for default(none) firstprivate(ndim, alpha, beta) private(j) shared(tmpA)
  for(j=0;j<=ndim-2;j++){
    tmpA[j][j]=alpha[j+1];
    tmpA[j][j+1]=beta[j+1];
    tmpA[j+1][j]=beta[j+1];
  }
  tmpA[ndim-1][ndim-1]=alpha[ndim];

  DSEVvector( ndim, tmpA, tmpr, tmpvec );
  fprintf(stdoutMPI, "  Lanczos EigenValue in vec12 = %.10lf \n ",tmpr[0]);
  
 
  if(nvec<=ndim){ 
#pragma omp parallel for default(none) firstprivate(ndim, nvec) private(j,k) shared(tmpvec, vec)
    for(k=1;k<=nvec;k++){
      for(j=1;j<=ndim;j++){
        vec[k][j]=tmpvec[k-1][j-1];
      } 
    }
  }
  else{
#pragma omp parallel for default(none) firstprivate(ndim, nvec) private(j,k) shared(tmpvec, vec)
    for(k=1;k<=ndim;k++){
      for(j=1;j<=ndim;j++){
        vec[k][j]=tmpvec[k-1][j-1];
      } 
    }
  }

  free(tmpA);
  free(tmpr);
  free(tmpvec);

}   
