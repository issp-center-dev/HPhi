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
@brief Functions for spin Hamiltonian (Core)
*/

#include <bitcalc.h>
#include "mfmemory.h"
#include "xsetmem.h"
#include "wrapperMPI.h"
#include "mltplyCommon.h"
#include "mltplySpinCore.h"

/******************************************************************************/
//[s] GetInfo functions
/******************************************************************************/

/**
 *
 *
 * @param iExchange
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int child_exchange_spin_GetInfo(
  int iExchange,
  struct BindStruct *X
) {
  int isite1 = X->Def.ExchangeCoupling[iExchange][0] + 1;
  int isite2 = X->Def.ExchangeCoupling[iExchange][1] + 1;
  X->Large.tmp_J = X->Def.ParaExchangeCoupling[iExchange];
  X->Large.is1_up = X->Def.Tpow[isite1 - 1];
  X->Large.is2_up = X->Def.Tpow[isite2 - 1];
  X->Large.isA_spin = X->Large.is1_up + X->Large.is2_up;
  return 0;
}/*int child_exchange_spin_GetInfo*/
/**
 *
 *
 * @param iPairLift
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int child_pairlift_spin_GetInfo(
  int iPairLift,
  struct BindStruct *X
) {
  int isite1 = X->Def.PairLiftCoupling[iPairLift][0] + 1;
  int isite2 = X->Def.PairLiftCoupling[iPairLift][1] + 1;
  X->Large.tmp_J = X->Def.ParaPairLiftCoupling[iPairLift];
  X->Large.is1_up = X->Def.Tpow[isite1 - 1];
  X->Large.is2_up = X->Def.Tpow[isite2 - 1];
  X->Large.isA_spin = X->Large.is1_up + X->Large.is2_up;
  return 0;
}/*int child_pairlift_spin_GetInfo*/
/**
 *
 *
 * @param X
 * @param isite1
 * @param isite2
 * @param sigma1
 * @param sigma2
 * @param sigma3
 * @param sigma4
 * @param tmp_V
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int child_general_int_spin_GetInfo(
  struct BindStruct *X,
  long unsigned int isite1,
  long unsigned int isite2,
  long unsigned int sigma1,
  long unsigned int sigma2,
  long unsigned int sigma3,
  long unsigned int sigma4,
  double complex tmp_V
) {
  X->Large.tmp_V = tmp_V;
  X->Large.isite1 = isite1;
  X->Large.isite2 = isite2;
  X->Large.is1_up = X->Def.Tpow[isite1 - 1];
  X->Large.is2_up = X->Def.Tpow[isite2 - 1];
  X->Large.is1_spin = sigma1;
  X->Large.is2_spin = sigma2;
  X->Large.is3_spin = sigma3;
  X->Large.is4_spin = sigma4;
  return 0;
}/*int child_general_int_spin_GetInfo*/

/******************************************************************************/
//[e] GetInfo functions
/******************************************************************************/

/******************************************************************************/
//[s] core routines
/******************************************************************************/

/**
 * @return XXX
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int X_Spin_CisAit(
  long unsigned int j,//!<[in]
  struct BindStruct *X,//!<[in]
  long unsigned int is1_spin,//!<[in]
  long unsigned int sigma2,//!<[in]
  long unsigned int *list_1_Org_,//!<[in]
  long unsigned int *list_2_1_,//!<[in]
  long unsigned int *list_2_2_,//!<[in]
  long unsigned int *tmp_off//!<[in]
) {
  long unsigned int list_1_j;
  long unsigned int off;
  list_1_j = list_1_Org_[j];
  if (X_SpinGC_CisAit(list_1_j + 1, X, is1_spin, sigma2, &off) != 0) {
    GetOffComp(list_2_1_, list_2_2_, off, X->Large.irght, X->Large.ilft, X->Large.ihfbit, tmp_off);
    return 1;
  }
  else {
    *tmp_off = 1;
    return 0;
  }
}/*int X_Spin_CisAit*/
/**
 *
 *
 * @param j
 * @param X
 * @param is1_spin
 * @param sigma1
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int X_Spin_CisAis(
  long unsigned int j,
  struct BindStruct *X,
  long unsigned int is1_spin,
  long unsigned int sigma1
) {
  int A_ibit_tmp;
  // off = j
  A_ibit_tmp = ((list_1[j] & is1_spin) / is1_spin) ^ (1 - sigma1);
  return A_ibit_tmp;
}/*int X_Spin_CisAis*/
/**
 *
 *
 * @param j
 * @param X
 * @param is1_spin
 * @param is2_spin
 * @param sigma1
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int X_Spin_CisAjs(
  long unsigned int j,
  struct BindStruct *X,
  long unsigned int is1_spin,
  long unsigned int is2_spin,
  long unsigned int sigma1,
  long unsigned int *tmp_off
) {
  int ibit_tmp;
  int jbit_tmp;
  long unsigned int iexchg;
  long unsigned int irght = X->Large.irght;
  long unsigned int ilft = X->Large.ilft;
  long unsigned int ihfbit = X->Large.ihfbit;
  long unsigned int is_up;
  is_up = is1_spin + is2_spin;

  //bit_tmp =0 -> down spin, 1 -> up spin
  ibit_tmp = ((list_1[j] & is1_spin) / is1_spin) ^ (1 - sigma1);//create: 1=(0^1, 1^0) is OK
  jbit_tmp = ((list_1[j] & is2_spin) / is2_spin) ^ (1 - sigma1);//anihilate: 0=(0^0, 1^1) is OK
    //1-sigma1 = 0 -> down spin, 1->up spin
  if (ibit_tmp != 0 && jbit_tmp == 0) {
    iexchg = list_1[j] ^ is_up;
    GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, tmp_off);
    return 1;
  }
  *tmp_off = 1;
  return 0;
}/*int X_Spin_CisAjs*/
/**
 *
 *
 * @param j
 * @param X
 * @param is1_spin
 * @param sigma1
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int X_SpinGC_CisAis(
  long unsigned int j,
  struct BindStruct *X,
  long unsigned int is1_spin,
  long unsigned int sigma1
) {
  int A_ibit_tmp;
  long unsigned int list_1_j;
  // off = j
  list_1_j = j - 1;
  A_ibit_tmp = ((list_1_j & is1_spin) / is1_spin) ^ (1 - sigma1);
  return A_ibit_tmp;
}/*int X_SpinGC_CisAis*/
/**
 *
 *
 * @param j
 * @param X
 * @param is1_spin
 * @param sigma2
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int X_SpinGC_CisAit(
  long unsigned int j,
  struct BindStruct *X,
  long unsigned int is1_spin,
  long unsigned int sigma2,
  long unsigned int *tmp_off
) {
  long unsigned int list_1_j, ibit_tmp_1;

  list_1_j = j - 1;

  ibit_tmp_1 = list_1_j & is1_spin;
  if (ibit_tmp_1 == 0 && sigma2 == 0) {    // down -> up
    *tmp_off = list_1_j + is1_spin;
    return 1;
  }
  else if (ibit_tmp_1 != 0 && sigma2 == 1) { // up -> down
    *tmp_off = list_1_j - is1_spin;
    return 1;
  }
  else {
    *tmp_off = 1;
    return 0;
  }
}/*int X_SpinGC_CisAit*/

/******************************************************************************/
//[e] core routines
/******************************************************************************/

/**
 *
 *
 * @param j
 * @param X
 * @param isA_up
 * @param isB_up
 * @param sigmaA
 * @param sigmaB
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int X_child_exchange_spin_element(
  long unsigned int j,
  struct BindStruct *X,
  long unsigned int isA_up,
  long unsigned int isB_up,
  long unsigned int sigmaA,
  long unsigned int sigmaB,
  long unsigned int *tmp_off
) {
  long unsigned int iexchg, off;
  long unsigned int irght = X->Large.irght;
  long unsigned int ilft = X->Large.ilft;
  long unsigned int ihfbit = X->Large.ihfbit;
  long unsigned int ibit_tmp_A, ibit_tmp_B;

  ibit_tmp_A = ((list_1[j] & isA_up) / isA_up);
  ibit_tmp_B = ((list_1[j] & isB_up) / isB_up);
  if (ibit_tmp_A == sigmaA && ibit_tmp_B == sigmaB) {
    iexchg = list_1[j] ^ (isA_up + isB_up);
    GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
    *tmp_off = off;
    return 1;
  }
  else {
    *tmp_off = 1; // just tentative
    return 0;
  }
}/*int X_child_exchange_spin_element*/
/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex child_exchange_spin_element(
  long unsigned int j,
  double complex *tmp_v0,
  double complex *tmp_v1,
  struct BindStruct *X,
  long unsigned int *tmp_off
) {
  long unsigned int off;
  double complex dmv;
  long unsigned int iexchg;
  long unsigned int is_up = X->Large.isA_spin;
  long unsigned int irght = X->Large.irght;
  long unsigned int ilft = X->Large.ilft;
  long unsigned int ihfbit = X->Large.ihfbit;
  double complex tmp_J = X->Large.tmp_J;
  int mode = X->Large.mode;
  double complex dam_pr = 0;
  long unsigned int ibit_tmp;

  ibit_tmp = (list_1[j] & is_up);
  if (ibit_tmp == 0 || ibit_tmp == is_up) {
    return dam_pr;
  }
  else {
    iexchg = list_1[j] ^ is_up;
    GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
    *tmp_off = off;
    dmv = tmp_J * tmp_v1[j];
    if (mode == M_MLTPLY) {
      tmp_v0[off] += dmv;
    }
    dam_pr += dmv * conj(tmp_v1[off]);
    return dam_pr;
  }
}/*double complex child_exchange_spin_element*/
/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex GC_child_exchange_spin_element(
  long unsigned int j,
  double complex *tmp_v0,
  double complex *tmp_v1,
  struct BindStruct *X,
  long unsigned int *tmp_off
) {
  double complex dmv;
  long unsigned int is_up = X->Large.isA_spin;
  double complex tmp_J = X->Large.tmp_J;
  int mode = X->Large.mode;
  long unsigned int list_1_j, list_1_off;

  double complex dam_pr = 0;
  list_1_j = j - 1;

  long unsigned int ibit_tmp;
  ibit_tmp = (list_1_j & is_up);
  if (ibit_tmp == 0 || ibit_tmp == is_up) {
    return dam_pr;
  }
  else {
    list_1_off = list_1_j ^ is_up;
    *tmp_off = list_1_off;
    dmv = tmp_J * tmp_v1[j];
    if (mode == M_MLTPLY) {
      tmp_v0[list_1_off + 1] += dmv;
    }
    dam_pr += dmv * conj(tmp_v1[list_1_off + 1]);
    return dam_pr;
  }
}/*double complex GC_child_exchange_spin_element*/
/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex child_pairlift_spin_element(
  long unsigned int j,
  double complex *tmp_v0,
  double complex *tmp_v1,
  struct BindStruct *X,
  long unsigned int *tmp_off
) {
  long unsigned int off;
  double complex dmv;
  long unsigned int iexchg;
  long unsigned int is1_up = X->Large.is1_up;
  long unsigned int is2_up = X->Large.is2_up;
  long unsigned int is_up = X->Large.isA_spin;
  long unsigned int irght = X->Large.irght;
  long unsigned int ilft = X->Large.ilft;
  long unsigned int ihfbit = X->Large.ihfbit;
  double complex tmp_J = X->Large.tmp_J;
  int mode = X->Large.mode;
  double complex dam_pr = 0;

  long unsigned int ibit_tmp;
  ibit_tmp = ((list_1[j] & is1_up) / is1_up) ^ ((list_1[j] & is2_up) / is2_up);
  if (ibit_tmp == 0) {
    iexchg = list_1[j] ^ is_up; //Change: ++ -> -- or -- -> ++
    GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
    dmv = tmp_J * tmp_v1[j] * ibit_tmp;
    *tmp_off = off;
    if (mode == M_MLTPLY) {
      tmp_v0[off] += dmv;
    }
    dam_pr += dmv * conj(tmp_v1[off]);
  }
  return dam_pr;
}/*double complex child_pairlift_spin_element*/
/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex GC_child_pairlift_spin_element(
  long unsigned int j,
  double complex *tmp_v0,
  double complex *tmp_v1,
  struct BindStruct *X,
  long unsigned int *tmp_off
) {
  double complex dmv;
  long unsigned int is_up = X->Large.isA_spin;
  double complex tmp_J = X->Large.tmp_J;
  int mode = X->Large.mode;
  double complex dam_pr = 0;
  long unsigned int list_1_off;
  long unsigned int list_1_j = j - 1;
  long unsigned int ibit_tmp;
  //ibit_tmp = ((list_1_j & is1_up) / is1_up) ^ ((list_1_j & is2_up) / is2_up);
  ibit_tmp = (list_1_j & is_up);
  if (ibit_tmp == 0 || ibit_tmp == is_up) {
    list_1_off = list_1_j ^ is_up; //Change: ++ -> -- or -- -> ++
    *tmp_off = list_1_off;
    dmv = tmp_J * tmp_v1[j];//* ibit_tmp;
    if (mode == M_MLTPLY) {
      tmp_v0[list_1_off + 1] += dmv;
    }
    dam_pr += dmv * conj(tmp_v1[list_1_off + 1]);
    return dam_pr;
  }
  else {
    return dam_pr;
  }
}/*double complex GC_child_pairlift_spin_element*/
//[s]Spin
/**
 *
 *
 * @param j
 * @param isA_up
 * @param isB_up
 * @param org_sigma2
 * @param org_sigma4
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex child_CisAisCisAis_spin_element(
  long unsigned int j,
  long unsigned int isA_up,
  long unsigned int isB_up,
  long unsigned int org_sigma2,
  long unsigned int org_sigma4,
  double complex tmp_V,
  double complex *tmp_v0,
  double complex *tmp_v1,
  struct BindStruct *X
) {
  int tmp_sgn;
  double complex dmv;
  double complex dam_pr = 0;

  tmp_sgn = X_Spin_CisAis(j, X, isB_up, org_sigma4);
  tmp_sgn *= X_Spin_CisAis(j, X, isA_up, org_sigma2);
  dmv = tmp_v1[j] * tmp_sgn * tmp_V;
  if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
    tmp_v0[j] += dmv;
  }
  dam_pr = conj(tmp_v1[j]) * dmv;
  return dam_pr;
}/*double complex child_CisAisCisAis_spin_element*/
/**
 *
 *
 * @param j
 * @param isA_up
 * @param isB_up
 * @param org_sigma2
 * @param org_sigma4
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex child_CisAjsCjtAit_spin_element(
  long unsigned int j,
  long unsigned int isA_up,
  long unsigned int isB_up,
  long unsigned int org_sigma2,
  long unsigned int org_sigma4,
  double complex tmp_V,
  double complex *tmp_v0,
  double complex *tmp_v1,
  struct BindStruct *X,
  long unsigned int *tmp_off
) {
  int tmp_sgn;
  long unsigned int tmp_off_1;
  double complex dmv = 0;
  double complex dam_pr = 0;
  tmp_sgn = X_Spin_CisAjs(j, X, isB_up, isA_up, org_sigma4, &tmp_off_1);
  if (tmp_sgn != 0) {
    tmp_sgn *= X_Spin_CisAjs(tmp_off_1, X, isA_up, isB_up, org_sigma2, tmp_off);
    if (tmp_sgn != 0) {
      dmv = tmp_v1[j] * tmp_V;
      if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
        tmp_v0[*tmp_off] += dmv;
      }
      dam_pr = conj(tmp_v1[*tmp_off]) * dmv;
    }
  }
  return dam_pr;
}/*double complex child_CisAjsCjtAit_spin_element*/
//[e]Spin

//[s]GC Spin
/**
 *
 *
 * @param j
 * @param isA_up
 * @param isB_up
 * @param org_sigma2
 * @param org_sigma4
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex GC_child_CisAisCisAis_spin_element(
  long unsigned int j,
  long unsigned int isA_up,
  long unsigned int isB_up,
  long unsigned int org_sigma2,
  long unsigned int org_sigma4,
  double complex tmp_V,
  double complex *tmp_v0,
  double complex *tmp_v1,
  struct BindStruct *X
) {
  int tmp_sgn;
  double complex dmv = 0;
  double complex dam_pr = 0;

  tmp_sgn = X_SpinGC_CisAis(j, X, isB_up, org_sigma4);
  tmp_sgn *= X_SpinGC_CisAis(j, X, isA_up, org_sigma2);
  if (tmp_sgn != 0) {
    dmv = tmp_v1[j] * tmp_sgn * tmp_V;
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[j] += dmv;
    }
    dam_pr = conj(tmp_v1[j]) * dmv;
  }
  return dam_pr;
}/*double complex GC_child_CisAisCisAis_spin_element*/
/**
 *
 *
 * @param j
 * @param org_sigma2
 * @param org_sigma4
 * @param isA_up
 * @param isB_up
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex GC_child_CisAisCitAiu_spin_element(
  long unsigned int j,
  long unsigned int org_sigma2,
  long unsigned int org_sigma4,
  long unsigned int isA_up,
  long unsigned int isB_up,
  double complex tmp_V,
  double complex *tmp_v0,
  double complex *tmp_v1,
  struct BindStruct *X,
  long unsigned int *tmp_off
) {
  int tmp_sgn;
  double complex dmv;
  double complex dam_pr = 0 + 0 * I;
  tmp_sgn = X_SpinGC_CisAit(j, X, isB_up, org_sigma4, tmp_off);
  if (tmp_sgn != 0) {
    tmp_sgn *= X_SpinGC_CisAis((*tmp_off + 1), X, isA_up, org_sigma2);
    if (tmp_sgn != 0) {
      dmv = tmp_v1[j] * tmp_sgn * tmp_V;
      if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
        tmp_v0[*tmp_off + 1] += dmv;
      }
      dam_pr = conj(tmp_v1[*tmp_off + 1]) * dmv;
    }/*if (tmp_sgn != 0)*/
  }/*if (tmp_sgn != 0)*/
  return dam_pr;
}/*double complex GC_child_CisAisCitAiu_spin_element*/
/**
 *
 *
 * @param j
 * @param org_sigma2
 * @param org_sigma4
 * @param isA_up
 * @param isB_up
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex GC_child_CisAitCiuAiu_spin_element(
  long unsigned int j,
  long unsigned int org_sigma2,
  long unsigned int org_sigma4,
  long unsigned int isA_up,
  long unsigned int isB_up,
  double complex tmp_V,
  double complex *tmp_v0,
  double complex *tmp_v1,
  struct BindStruct *X,
  long unsigned int *tmp_off
) {
  int tmp_sgn;
  double complex dmv;
  double complex dam_pr = 0 + 0 * I;
  tmp_sgn = X_SpinGC_CisAis(j, X, isB_up, org_sigma4);
  if (tmp_sgn != 0) {
    tmp_sgn *= X_SpinGC_CisAit(j, X, isA_up, org_sigma2, tmp_off);
    if (tmp_sgn != 0) {
      dmv = tmp_v1[j] * tmp_sgn * tmp_V;
      if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
        tmp_v0[*tmp_off + 1] += dmv;
      }
      dam_pr = conj(tmp_v1[*tmp_off + 1]) * dmv;
    }/*if (tmp_sgn != 0)*/
  }/*if (tmp_sgn != 0)*/
  return dam_pr;
}/*double complex GC_child_CisAitCiuAiu_spin_element*/
/**
 *
 *
 * @param j
 * @param org_sigma2
 * @param org_sigma4
 * @param isA_up
 * @param isB_up
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off_2
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex GC_child_CisAitCiuAiv_spin_element(
  long unsigned int j,
  long unsigned int org_sigma2,
  long unsigned int org_sigma4,
  long unsigned int isA_up,
  long unsigned int isB_up,
  double complex tmp_V,
  double complex *tmp_v0,
  double complex *tmp_v1,
  struct BindStruct *X,
  long unsigned int *tmp_off_2
) {
  int tmp_sgn;
  long unsigned int tmp_off_1;
  double complex dmv;
  double complex dam_pr = 0 + 0 * I;
  tmp_sgn = X_SpinGC_CisAit(j, X, isB_up, org_sigma4, &tmp_off_1);
  if (tmp_sgn != 0) {
    tmp_sgn *= X_SpinGC_CisAit((tmp_off_1 + 1), X, isA_up, org_sigma2, tmp_off_2);
    if (tmp_sgn != 0) {
      dmv = tmp_v1[j] * tmp_sgn * tmp_V;
      if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
        tmp_v0[*tmp_off_2 + 1] += dmv;
      }
      dam_pr = conj(tmp_v1[*tmp_off_2 + 1]) * dmv;
    }/*if (tmp_sgn != 0)*/
  }/*if (tmp_sgn != 0)*/
  return dam_pr;
}/*double complex GC_child_CisAitCiuAiv_spin_element*/
//[e]GC Spin
