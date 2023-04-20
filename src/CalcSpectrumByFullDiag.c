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
/**@file
@brief Functions to perform spectrum calculations with the
full-diagonalization method.
*/
#include <complex.h>
#include "global.h"
#include <time.h>
#include "struct.h"
#include "lapack_diag.h"
#include "mltply.h"
#include "mltplyCommon.h"
#include "CalcTime.h"
#include "common/setmemory.h"
#include "CalcSpectrum.h"
/**
@brief Compute the Green function with the Lehmann representation and FD
@f[
G(z) = \sum_n \frac{|\langle n|c|0\rangle|^2}{z - E_n}
@f]
@author Mitsuaki Kawamura (The University of Tokyo)
*/
int CalcSpectrumByFullDiag(
  struct EDMainCalStruct *X,//!<[inout]
  int Nomega,//!<[in] Number of frequencies
  int NdcSpectrum,
  double complex ***dcSpectrum,//!<[out] [Nomega] Spectrum
  double complex **dcomega,//!<[in] [Nomega] Frequency
  double complex **v1Org
)
{
  int idim0, idim1, idim2, iomega;
  int idim_max_int, idim_maxorg_int;
  int idcSpectrum;
  double complex **vR, **vL, **vRv, **vLv;
  /**
  <ul>
  <li>Generate fully stored Hamiltonian. Because ::v0 & ::v1 are overwritten,
  copy ::v0 into ::vg.</li>
  */
  idim_max_int = (int)X->Bind.Check.idim_max;
  idim_maxorg_int = (int)X->Bind.Check.idim_maxOrg;
  vR = cd_2d_allocate(idim_max_int+1, idim_maxorg_int);
  vL = cd_2d_allocate(idim_max_int+1, idim_maxorg_int);
  vLv = cd_2d_allocate(idim_max_int, idim_maxorg_int);
  vRv = cd_2d_allocate(idim_max_int, idim_maxorg_int);

  StartTimer(6301);
  v0[0][0] = 0.0;
  zclear((X->Bind.Check.idim_max + 1)*X->Bind.Check.idim_max, &v0[0][0]);
  zclear((X->Bind.Check.idim_max + 1)*X->Bind.Check.idim_max, &v1[0][0]);
  for (idim0 = 1; idim0 <= X->Bind.Check.idim_max; idim0++) v1[idim0][idim0] = 1.0;
  mltply(&(X->Bind), X->Bind.Check.idim_max, v0, v1);
  StopTimer(6301);
  /**
  <li>::v0 becomes eigenvalues in lapack_diag(), and
   ::v1 becomes eigenvectors</li>
  */
  StartTimer(6302);
  lapack_diag(&(X->Bind));
  StopTimer(6302);
  /**
  <li>Compute @f$|\langle n|c|0\rangle|^2@f$ for all @f$n@f$ and store them into ::v1,
  where @f$c|0\rangle@f$ is ::vg.</li>
  */
  zclear(X->Bind.Check.idim_max, &vR[1][0]);
  GetExcitedState(&(X->Bind), X->Bind.Check.idim_maxOrg, vR, v1Org, 0);
  for (idim0 = 1; idim0 < idim_max_int+1; idim0++)
    for (idim1 = 0; idim1 < idim_max_int; idim1++)
      for (idim2 = 0; idim2 < idim_maxorg_int; idim2++)
        vRv[idim1][idim2] += conj(v0[idim0][idim1]) * vR[idim0][idim2];
  for (idcSpectrum = 0; idcSpectrum < NdcSpectrum; idcSpectrum++) {
    StartTimer(6303);
    zclear(X->Bind.Check.idim_max, &vL[1][0]);
    GetExcitedState(&(X->Bind), X->Bind.Check.idim_maxOrg, vL, v1Org, idcSpectrum + 1);
    zclear(X->Bind.Check.idim_max* X->Bind.Check.idim_max, &vLv[0][0]);
    for (idim0 = 1; idim0 < idim_max_int + 1; idim0++)
      for (idim1 = 0; idim1 < idim_max_int; idim1++)
        for (idim2 = 0; idim2 < idim_maxorg_int; idim2++)
          vLv[idim1][idim2] += conj(v0[idim0][idim1]) * vL[idim0][idim2];
    StopTimer(6303);
    /**
    <li>Compute spectrum
    @f[
    \sum_n \frac{|\langle n|c|0\rangle|^2}{z - E_n}
    @f]
    </li>
    </ul>
    */
    StartTimer(6304);
    for (idim0 = 0; idim0 < idim_maxorg_int; idim0++) {
      for (iomega = 0; iomega < Nomega; iomega++) {
        dcSpectrum[idim0][iomega][idcSpectrum] = 0.0;
        for (idim1 = 0; idim1 < idim_max_int; idim1++) {
          dcSpectrum[idim0][iomega][idcSpectrum] += conj(vLv[idim1][idim0]) * vRv[idim1][idim0]
            / (dcomega[idim0][iomega] - X->Bind.Phys.energy[idim1]);
        }/*for (idim = 0; idim < idim_max_int; idim++)*/
      }/*for (iomega = 0; iomega < Nomega; iomega++)*/
    }
    StopTimer(6304);
  }/*for (idcSpectrum = 1; idcSpectrum < NdcSpectrum; idcSpectrum++)*/
  free_cd_2d_allocate(vL);
  free_cd_2d_allocate(vR);
  free_cd_2d_allocate(vLv);
  free_cd_2d_allocate(vRv);
  return TRUE;
}/*CalcSpectrumByFullDiag*/

