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
/**@file
@brief Functions to Diagonalize a tri-diagonal matrix and store eigenvectors 
into ::vec
*/
#include "matrixlapack.h"
#include "Common.h"
#include "wrapperMPI.h"
#include "common/setmemory.h"
#include "xsetmem.h"
/**
@brief Diagonalize a tri-diagonal matrix and store eigenvectors 
into ::vec
@author Takahiro Misawa (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
void vec12(
  double alpha[],
  double beta[],
  unsigned int ndim,
  double tmp_E[], 
  struct BindStruct *X
) {
  unsigned int j,k,nvec;

  double **tmpA, **tmpvec;
  double *tmpr;
    
  nvec = X->Def.nvec;
  tmpA = d_2d_allocate(ndim, ndim);
  tmpvec = d_2d_allocate(ndim, ndim);
  tmpr = d_1d_allocate(ndim);

#pragma omp parallel for default(none) firstprivate(ndim) private(j,k) shared(tmpA)
  for(k=0;k<=ndim-1;k++)
    for(j=0;j<=ndim-1;j++) tmpA[k][j]=0.0;
#pragma omp parallel for default(none) firstprivate(ndim, nvec) private(j,k) shared(vec)
  for(k=1;k<=nvec;k++)
    for(j=1;j<=ndim;j++) vec[k][j]=0.0;

#pragma omp parallel for default(none) firstprivate(ndim, alpha, beta) private(j) shared(tmpA)
  for(j=0;j<=ndim-2;j++){
    tmpA[j][j]=alpha[j+1];
    tmpA[j][j+1]=beta[j+1];
    tmpA[j+1][j]=beta[j+1];
  }/*for(j=0;j<=ndim-2;j++)*/
  tmpA[ndim-1][ndim-1]=alpha[ndim];

  DSEVvector( ndim, tmpA, tmpr, tmpvec );
  if(X->Def.iCalcType==Lanczos && X->Def.iFlgCalcSpec == 0)
    fprintf(stdoutMPI, "  Lanczos EigenValue in vec12 = %.10lf \n ",tmpr[0]);
 
  if (nvec <= ndim) {
    if (nvec < X->Def.LanczosTarget) nvec = X->Def.LanczosTarget;
    
#pragma omp parallel for default(none) firstprivate(ndim, nvec) private(j,k) shared(tmpvec, vec, tmp_E, tmpr)
    for(k=1;k<=nvec;k++){
      tmp_E[k]=tmpr[k-1];
      for (j = 1; j <= ndim; j++) vec[k][j] = tmpvec[k - 1][j - 1];
    }/*for(k=1;k<=nvec;k++)*/
  }/*if(nvec<=ndim)*/
  else{
#pragma omp parallel for default(none) firstprivate(ndim, nvec) private(j,k) shared(tmpvec, vec, tmp_E, tmpr)
    for(k=1;k<=ndim;k++){
      tmp_E[k]=tmpr[k-1];
      for (j = 1; j <= ndim; j++) vec[k][j] = tmpvec[k - 1][j - 1];
    }/*for(k=1;k<=ndim;k++)*/
  }/*if(nvec>ndim)*/
  free_d_2d_allocate(tmpA);
  free_d_1d_allocate(tmpr);
  free_d_2d_allocate(tmpvec);
}/*void vec12*/
