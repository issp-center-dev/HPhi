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
#include "global.h"
#include <time.h>
#include "struct.h"
#include "lapack_diag.h"
#include "makeHam.h"
#include "CalcTime.h"

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
  int idim, jdim, iomega;
  int idim_max_int;
  int incr=1;
  
  idim_max_int = (int)X->Bind.Check.idim_max;
  zcopy_(&idim_max_int, &v0[1], &incr, &vg[0], &incr);
  /*
   In generating Hamiltonian, v0 & v1 are overwritten.
  */
  StartTimer(6301);
  makeHam(&(X->Bind));
  StopTimer(6301);
  /*
   v0 : becomes eigenvalues in lapack_diag()
   L_vec : Eigenvectors
  */
  StartTimer(6302);
  lapack_diag(&(X->Bind));
  StopTimer(6302);
  /*
   v1 = |<n|c|0>|^2 
  */
  StartTimer(6303);

  for (idim = 0; idim < idim_max_int; idim++) {
    v1[idim] = 0.0;
    for (jdim = 0; jdim < idim_max_int; jdim++) v1[idim] += conj(vg[jdim]) * L_vec[idim][jdim];
    //zdotc_(&v1[idim], &idim_max_int, &vg[0], &incr, &L_vec[idim][0], &incr);
    v1[idim] = conj(v1[idim]) * v1[idim];
  }/*for (idim = 0; idim < idim_max_int; idim++)*/
  StopTimer(6303);
  /*
   Sum_n |<n|c|0>|^2 / (z - E_n)
  */
  StartTimer(6304);
  for (iomega = 0; iomega < Nomega; iomega++) {
    dcSpectrum[iomega] = 0.0;
    for (idim = 0; idim < idim_max_int; idim++) {
      dcSpectrum[iomega] += v1[idim] / (dcomega[iomega] - v0[idim]);
    }/*for (idim = 0; idim < idim_max_int; idim++)*/
  }/*for (iomega = 0; iomega < Nomega; iomega++)*/
  StopTimer(6304);
  return TRUE;
}/*CalcSpectrumByFullDiag*/

