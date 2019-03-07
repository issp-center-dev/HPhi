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

// Define Mode for mltply
// complex version
#include <bitcalc.h>
#include "mltply.h"
#include "mltplySpin.h"
#include "mltplyHubbard.h"
#include "wrapperMPI.h"
#include "CalcTime.h"
#include "mltplyCommon.h"
#include "diagonalcalc.h"

/**
 * @file   mltply.c
 *
 * @brief  Multiplying the wavefunction by the Hamiltonian. @f$ H v_1@f$.
 *
 * @version 0.2
 * @details add function to treat the case of generalspin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */


/**
 * @brief Parent function of multiplying the wavefunction by the Hamiltonian. @f$ H v_1@f$.\n
 * First, the calculation of diagonal term is done by using the list @f$ \verb|list_diaognal| @f$. \n
 * Next, the calculation of off-diagonal term is done.\n
 * @note If @f$ \verb|mode| @f$ in BindStruct X is @f$ \verb|M_CORR| @f$, the wave function is not updated. The expected values are only calculated.\n
 * Otherwise, the wavefunction @f$ v_0 @f$ is updated as @f$ v_0 += H v_1@f$.
 *
 * @param X [in] Struct for getting the information of the operators.
 * @param tmp_v0 [in, out]
 * @param tmp_v1 [in]
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int mltply(struct BindStruct *X, int nstate, double complex **tmp_v0,double complex **tmp_v1) {
  int one = 1;
  long unsigned int j=0;
  long unsigned int irght=0;
  long unsigned int ilft=0;
  long unsigned int ihfbit=0;


  long unsigned int i_max;

  StartTimer(1);
  i_max = X->Check.idim_max;
  X->Large.prdct = 0.0;

  if(i_max!=0){
    if (X->Def.iFlgGeneralSpin == FALSE) {
      if (GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit) != 0) {
        return -1;
      }
    }
    else{
      if(X->Def.iCalcModel==Spin){
        if (GetSplitBitForGeneralSpin(X->Def.Nsite, &ihfbit, X->Def.SiteToBit) != 0) {
          return -1;
        }
      }
    }
  }
  else{
    irght=0;
    ilft=0;
    ihfbit=0;
  }
  X->Large.i_max = i_max;
  X->Large.irght = irght;
  X->Large.ilft = ilft;
  X->Large.ihfbit = ihfbit;
  X->Large.mode = M_MLTPLY;

  StartTimer(100);
#pragma omp parallel for default(none)  firstprivate(i_max) shared(tmp_v0, tmp_v1, list_Diagonal)
  for (j = 1; j <= i_max; j++) {
    zaxpy_(&nstate, &list_Diagonal[j], &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
  }
  StopTimer(100);
  if (X->Def.iCalcType == TimeEvolution) diagonalcalcForTE(step_i, X, nstate, &tmp_v0[0][0], &tmp_v1[0][0]);
  
  switch (X->Def.iCalcModel) {
  case HubbardGC:
    mltplyHubbardGC(X, nstate, tmp_v0, tmp_v1);
    break;
      
  case KondoGC:
  case Hubbard:
  case Kondo:
    mltplyHubbard(X, nstate, tmp_v0, tmp_v1);
    break;
      
  case Spin:
    mltplySpin(X, nstate, tmp_v0, tmp_v1);
    break;
      
  case SpinGC:
    mltplySpinGC(X, nstate, tmp_v0, tmp_v1);
    break;
      
  default:
    return -1;
  }
  
  X->Large.prdct = SumMPI_dc(X->Large.prdct);  
  StopTimer(1);
  return 0;
}
/**
@brief Wrapper of zaxpy.
*/
void zaxpy_long(
  unsigned long int n, 
  double complex a, 
  double complex *x, 
  double complex *y
) {
  unsigned long int i;

#pragma omp parallel for default(none) private(i) shared(n, a, x, y)
  for (i = 0; i < n; i++) 
    y[i] += a * x[i] + y[i];
}
/**
@brief Wrapper of zswap.
*/
void zswap_long(
  unsigned long int n,
  double complex *x,
  double complex *y
) {
  unsigned long int i;
  double complex x0;

#pragma omp parallel for default(none) private(i,x0) shared(n, x, y)
  for (i = 0; i < n; i++) {
    x0 = x[i];
    x[i] = y[i];
    y[i] = x0;
  }
}
/**
@brief Wrapper of zswap.
*/
void zclear(
  unsigned long int n,
  double complex *x
) {
  unsigned long int i;

#pragma omp parallel for default(none) private(i) shared(n, x)
  for (i = 0; i < n; i++) 
    x[i] = 0.0;
}
