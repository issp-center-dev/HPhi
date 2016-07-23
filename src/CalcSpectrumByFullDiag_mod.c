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

#include <complex.h>
#include <math.h>
#include "global.h"
#include "time.h"
#include "struct.h"
#include "lapack_diag.h"
#include "makeHam.h"
#include "wrapperMPI.h"
#include "mfmemory.h"
#include "CalcSpectrum.h"

void zcopy_(int *n, double complex *x, int *incx, double complex *y, int *incy);
void zdotc_(double complex *xy, int *n, double complex *x, int *incx, double complex *y, int *incy);
/**
*
* Compute the Green function with the Lehmann representation and FD
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
int CalcSpectrumByFullDiag(
  struct EDMainCalStruct *X,
  int Nomega,
  double complex *dcSpectrum,
  double complex *dcomega)
{
  int idim, iomega;
  int idim_max_int;
  int incr=1;

  double complex *tmpv0;
  double complex *tmpv1;
  double complex tmpnorm;
  double tmpdnorm; 
 
  idim_max_int = (int)X->Bind.Check.idim_max;
  zcopy_(&idim_max_int, &v0[1], &incr, &vg[0], &incr);
  /*
   In generating Hamiltonian, v0 & v1 are overwritten.
  */
  makeHam(&(X->Bind));
  /*
   v0 : becomes eigenvalues in lapack_diag()
   L_vec : Eigenvectors
  */
  lapack_diag(&(X->Bind));
  
  c_malloc1(tmpv0, (unsigned long int) (idim_max_int+1));
  c_malloc1(tmpv1, (unsigned long int) (idim_max_int+1));

//  for (idim = 0; idim < idim_max_int; idim++) {
//    tmpv1[idim+1] = L_vec[0][idim];
//  }
  zcopy_(&idim_max_int, &L_vec[0][0], &incr, &tmpv1[1], &incr);

  for (idim = 0; idim < idim_max_int; idim++) {
    if(isnan(creal(tmpv1[(unsigned long int) (idim+1)]))) {
      fprintf(stdoutMPI, " #CSFD nan %d \n",idim);
      break;
    }
    if(isnan(cimag(tmpv1[(unsigned long int) (idim+1)]))) {
      fprintf(stdoutMPI, " #CSFD nan %d \n",idim);
      break;
    }
  }

  fprintf(stdoutMPI, " #CSFD1  %.10lf %.10lf \n",creal(tmpv1[9]),cimag(tmpv1[9]));
  GetExcitedState(&(X->Bind), tmpv0, tmpv1);
  fprintf(stdoutMPI, " #CSFD2  %.10lf %.10lf \n",creal(tmpv0[9]),cimag(tmpv0[9]));
  tmpnorm=0;

  for (idim = 0; idim < idim_max_int; idim++) {
  //for (idim = 0; idim < 1; idim++) {
    tmpnorm += conj(tmpv0[(unsigned long int) (idim+1)])*tmpv0[(unsigned long int) (idim+1)];
    //tmpnorm = (double) idim;
    if(isnan(creal(tmpv0[(unsigned long int) (idim+1)]))) {
      fprintf(stdoutMPI, " #CSFD nan %d \n",idim);
      break;
    }
  }
  fprintf(stdoutMPI, " #CSFD2.0  %.10lf %.10lf \n",creal(tmpnorm),cimag(tmpnorm));
  tmpdnorm = creal(tmpnorm);
  tmpdnorm = sqrt(tmpdnorm);
  fprintf(stdoutMPI, " #CSFD2.1  %.10lf \n",tmpdnorm);
  for (idim = 0; idim < idim_max_int; idim++) {
    tmpv0[(unsigned long int) (idim+1)] = tmpv0[(unsigned long int) (idim+1)] / tmpdnorm; 
  }
//  #CSFD1  -0.0002369629 0.0000000000 
// #CSFD2  0.0000000000 0.0001077884 
// #CSFD2.1  -nan 
// #CSFD3.0  -nan -nan 
// #CSFD3  -nan -nan 
// #CSFD3.1  -nan -nan 
// #CSFD4  -nan nan 
//  End:  Calculating a spectrum.

//  for (idim = 0; idim < idim_max_int; idim++) {
//    vg[idim] = tmpv0[idim+1];
//  }
  zcopy_(&idim_max_int, &tmpv0[1], &incr, &vg[0], &incr);
  fprintf(stdoutMPI, " #CSFD3.0  %.10lf %.10lf \n",creal(tmpv0[1]),cimag(tmpv0[1]));

  c_free1(tmpv0, (unsigned long int) (idim_max_int+1)); 
  c_free1(tmpv1, (unsigned long int) (idim_max_int+1)); 

  fprintf(stdoutMPI, " #CSFD3  %.10lf %.10lf \n",creal(vg[9]),cimag(vg[9]));
  fprintf(stdoutMPI, " #CSFD3.1  %.10lf %.10lf \n",creal(vg[0]),cimag(vg[0]));

  /*
   v1 = |<n|c|0>|^2 
  */
  for (idim = 0; idim < idim_max_int; idim++) {
    zdotc_(&v1[idim], &idim_max_int, &vg[0], &incr, &L_vec[idim][0], &incr);
    /*
    v1[idim] = 0;
    for (iomega=0; iomega < idim_max_int; iomega++) {
      v1[idim] += conj(vg[iomega])*L_vec[idim][iomega];
    }
    */
    v1[idim] = conj(v1[idim]) * v1[idim];
  }/*for (idim = 0; idim < idim_max_int; idim++)*/
  fprintf(stdoutMPI, " #CSFD4  %.10lf %.10lf \n",creal(v1[9]),cimag(v1[9]));
  /*
   Sum_n |<n|c|0>|^2 / (z - E_n)
  */
  for (iomega = 0; iomega < Nomega; iomega++) {
    dcSpectrum[iomega] = 0.0;
    for (idim = 0; idim < idim_max_int; idim++) {
      dcSpectrum[iomega] += v1[idim] / (dcomega[iomega] - v0[idim]);
    }/*for (idim = 0; idim < idim_max_int; idim++)*/
  }/*for (iomega = 0; iomega < Nomega; iomega++)*/

  return TRUE;

}/*CalcSpectrumByFullDiag*/
