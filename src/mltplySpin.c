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
@brief Functions for spin Hamiltonian

- mltplySpin() : Main routine of spin Hamiltonian (canonical)
  - mltplyHalfSpin() : 1/2 spin
  - mltplyGeneralSpin() : general spin
- mltplySpinGC() : Main routine of spin Hamiltonian (grandcanonical)
  - mltplyHalfSpinGC() : 1/2 spin
  - mltplyGeneralSpinGC() : general spin
  - mltplySpinGCBoost() : 

Hub routines
<table>
  <tr>
    <td></td>
    <td>Get info</td>
    <td>Canonical</td>
    <td>Grandcanonical</td>
  </tr>
  <tr>
    <td>Exchange</td>
    <td>::exchange_spin_GetInfo</td>
    <td>::exchange_spin, ::exchange_spin_element</td>
    <td>::GC_exchange_spin, ::GC_exchange_spin_element</td>
  </tr>
  <tr>
    <td>Pair lift</td>
    <td>::pairlift_spin_GetInfo</td>
    <td></td>
    <td>::GC_pairlift_spin, ::GC_pairlift_spin_element</td>
  </tr>
  <tr>
    <td>General int.</td>
    <td>::general_int_spin_GetInfo</td>
    <td>::general_int_spin, ::general_int_spin_MPIsingle
    ::child_general_int_spin_MPIsingle, ::general_int_spin_MPIdouble,
    ::child_general_int_spin_MPIdouble</td>
    <td>::GC_general_int_spin, ::GC_general_int_spin_MPIsingle,
    ::GC_general_int_spin_MPIdouble</td>
  </tr>
  <tr>
    <td>General int for 1/2 spin</td>
    <td>::general_int_spin_GetInfo</td>
    <td>::general_int_spin, ::general_int_spin_MPIsingle
    ::child_general_int_spin_MPIsingle, ::general_int_spin_MPIdouble,
    ::child_general_int_spin_MPIdouble</td>
    <td>::GC_general_int_spin, ::GC_general_int_spin_MPIsingle,
    ::GC_general_int_spin_MPIdouble</td>
  </tr>
  <tr>
    <td>General int for general spin</td>
    <td></td>
    <td>::general_int_GeneralSpin_MPIsingle,
    ::general_int_GeneralSpin_MPIdouble</td>
    <td>::GC_general_int_GeneralSpin_MPIsingle,
    ::GC_general_int_GeneralSpin_MPIdouble</td>
  </tr>
</table>

General on-site term
<table>
  <tr>
    <td></td>
    <td>1/2 spin</td>
    <td>1/2 spin</td>
    <td>1/2 spin</td>
    <td>1/2 spin</td>
    <td>General spin</td>
    <td>General spin</td>
  </tr>
  <tr>
    <td></td>
    <td>Canonical</td>
    <td>Canonical</td>
    <td>Grand canonical</td>
    <td>Grand canonical</td>
    <td>Canonical</td>
    <td>Grand canonical</td>
  </tr>
  <tr>
    <td></td>
    <td>In process</td>
    <td>Across process</td>
    <td>In process</td>
    <td>Across process</td>
    <td>Across process</td>
    <td>Across process</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i s}@f$</td>
    <td>::child_Spin_CisAis</td>
    <td>::child_CisAis_spin_MPIdouble</td>
    <td>::child_SpinGC_CisAis</td>
    <td>::child_GC_CisAis_spin_MPIdouble</td>
    <td>::child_CisAis_GeneralSpin_MPIdouble</td>
    <td>::child_GC_CisAis_GeneralSpin_MPIdouble</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i t}@f$</td>
    <td>::child_Spin_CisAit</td>
    <td>::child_CisAit_spin_MPIdouble</td>
    <td>::child_SpinGC_CisAit</td>
    <td>::child_GC_CisAit_spin_MPIdouble</td>
    <td>::child_CisAit_GeneralSpin_MPIdouble</td>
    <td>::child_GC_CisAit_GeneralSpin_MPIdouble</td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i s} c_{i s}^\dagger c_{i s}@f$</td>
    <td>::CisAisCisAis_spin_element</td>
    <td></td>
    <td>::GC_CisAisCisAis_spin_element</td>
    <td></td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i s} c_{i t}^\dagger c_{i u}@f$</td>
    <td></td>
    <td></td>
    <td>::GC_CisAisCitAiu_spin_element</td>
    <td></td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i t} c_{i u}^\dagger c_{i u}@f$</td>
    <td></td>
    <td></td>
    <td>::GC_CisAitCiuAiu_spin_element</td>
    <td></td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>@f$c_{i s}^\dagger c_{i t} c_{i u}^\dagger c_{i v}@f$</td>
    <td></td>
    <td></td>
    <td>::GC_CisAitCiuAiv_spin_element</td>
    <td>::GC_CisAitCiuAiv_spin_MPIsingle, ::child_GC_CisAitCiuAiv_spin_MPIsingle,
    ::GC_CisAitCiuAiv_spin_MPIdouble, ::child_GC_CisAitCiuAiv_spin_MPIdouble</td>
    <td></td>
    <td></td>
  </tr>
</table>
*/
#include <bitcalc.h>
#include "common/setmemory.h"
#include "mltplyCommon.h"
#include "mltplySpin.h"
#include "CalcTime.h"
#include "mltplySpinCore.h"
#include "mltplyHubbardCore.h"
#include "mltplyMPISpin.h"
#include "mltplyMPISpinCore.h"
#include "mltplyMPIBoost.h"
/**
@brief Driver function for Spin hamiltonian
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltplySpin(
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1//!<[in] Input producted vector
) {
  int iret=0;
  if (X->Def.iFlgGeneralSpin == FALSE)
    iret = mltplyHalfSpin(X, tmp_v0, tmp_v1);
  else
    iret = mltplyGeneralSpin(X, tmp_v0, tmp_v1);
  return iret;
}/*int mltplySpin*/
/**
@brief Driver function for Spin 1/2 hamiltonian
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltplyHalfSpin(
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1//!<[in] Input producted vector
) {
  long unsigned int i;
  long unsigned int isite1, isite2, sigma1, sigma2;
  long unsigned int sigma3, sigma4;

  double complex dam_pr;
  /*[s] For InterAll */
  double complex tmp_V;
  /*[e] For InterAll */
  int ihermite=0;
  int idx=0;

  StartTimer(400);
  /**
  Transfer absorbed in Diagonal term.
  InterAll
  */
  StartTimer(410);
  for (i = 0; i < X->Def.NInterAll_OffDiagonal; i+=2) {
    if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite &&
        X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
      StartTimer(411);
      general_int_spin_MPIdouble(i, X, tmp_v0, tmp_v1);
      StopTimer(411);
    }
    else if (X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
      StartTimer(412);
      general_int_spin_MPIsingle(i, X, tmp_v0, tmp_v1);
      StopTimer(412);
    }
    else if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite) {
      StartTimer(413);
      general_int_spin_MPIsingle(i + 1, X, tmp_v0, tmp_v1);
      StopTimer(413);
    }
    else {
      StartTimer(414);
      for (ihermite = 0; ihermite<2; ihermite++) {
        idx = i + ihermite;
        isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
        isite2 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
        sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
        sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
        sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
        sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
        tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];
        general_int_spin_GetInfo(X, isite1, isite2, sigma1, sigma2, sigma3, sigma4, tmp_V);
        dam_pr = general_int_spin(tmp_v0, tmp_v1, X);
        X->Large.prdct += dam_pr;
      }/*for (ihermite = 0; ihermite<2; ihermite++)*/
      StopTimer(414);
    }
  }/*for (i = 0; i < X->Def.NInterAll_OffDiagonal; i+=2)*/
  StopTimer(410);
  /**
  Exchange 
  */
  StartTimer(420);   
  for (i = 0; i < X->Def.NExchangeCoupling; i++) {
    sigma1=0; sigma2=1;
    if (X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite &&
        X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite) {
      StartTimer(421);
      dam_pr = child_general_int_spin_MPIdouble(
        X->Def.ExchangeCoupling[i][0], sigma1, sigma2, 
        X->Def.ExchangeCoupling[i][1], sigma2, sigma1, 
        X->Def.ParaExchangeCoupling[i], X, tmp_v0, tmp_v1);
      StopTimer(421);
    }
    else if (X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite) {
      StartTimer(422);
      dam_pr = child_general_int_spin_MPIsingle(
        X->Def.ExchangeCoupling[i][0], sigma1, sigma2, 
        X->Def.ExchangeCoupling[i][1], sigma2, sigma1,
        X->Def.ParaExchangeCoupling[i], X, tmp_v0, tmp_v1);
      StopTimer(422);
    }
    else if (X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite) {
      StartTimer(423);
      dam_pr = child_general_int_spin_MPIsingle(
        X->Def.ExchangeCoupling[i][1], sigma2, sigma1, 
        X->Def.ExchangeCoupling[i][0], sigma1, sigma2, 
        conj(X->Def.ParaExchangeCoupling[i]), X, tmp_v0, tmp_v1);
      StopTimer(423);
    }
    else {
      StartTimer(424);
      exchange_spin_GetInfo(i, X);
      dam_pr = exchange_spin(tmp_v0, tmp_v1, X);
      StopTimer(424);
    }
    X->Large.prdct += dam_pr;
  }/*for (i = 0; i < X->Def.NExchangeCoupling; i += 2)*/
  StopTimer(420);

  StopTimer(400);
  return 0;
}/*int mltplyHalfSpin*/
/**
@brief Driver function for General Spin hamiltonian
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltplyGeneralSpin(
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1//!<[in] Input producted vector
){
  long unsigned int j;
  long unsigned int i;
  long unsigned int off = 0;
  long unsigned int tmp_off = 0;
  long unsigned int tmp_off2 = 0;
  long unsigned int ihfbit=0;
  long unsigned int isite1, isite2, sigma1, sigma2;
  long unsigned int sigma3, sigma4;

  double complex dam_pr;
  long int tmp_sgn;
  /*[s] For InterAll */
  double complex tmp_V;
  double complex dmv=0;
  /*[e] For InterAll */

  long unsigned int i_max;
  i_max = X->Check.idim_max;
  int ihermite=0;
  int idx=0;

  StartTimer(400);
  /**
  Transfer absorbed in Diagonal term.
  InterAll
  */
  StartTimer(410);
  ihfbit =X->Check.sdim;
  for (i = 0; i < X->Def.NInterAll_OffDiagonal; i += 2) {
    if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite &&
        X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
      StartTimer(411);
      general_int_GeneralSpin_MPIdouble(i, X, tmp_v0, tmp_v1);
      StopTimer(411);
    }
    else if (X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
      StartTimer(412);
      general_int_GeneralSpin_MPIsingle(i, X, tmp_v0, tmp_v1);
      StopTimer(412);
    }
    else if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite) {
      StartTimer(412);
      general_int_GeneralSpin_MPIsingle(i + 1, X, tmp_v0, tmp_v1);
      StopTimer(412);
    }
    else {
      StartTimer(413);
      for (ihermite = 0; ihermite < 2; ihermite++) {
        idx = i + ihermite;
        isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
        isite2 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
        sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
        sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
        sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
        sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
        tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];
        dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) \
  private(j, tmp_sgn, dmv, off, tmp_off, tmp_off2) \
  firstprivate(i_max, isite1, isite2, sigma1, sigma2, sigma3, sigma4, X, tmp_V, ihfbit) \
  shared(tmp_v0, tmp_v1, list_1, list_2_1, list_2_2)
        for (j = 1; j <= i_max; j++) {
          tmp_sgn = GetOffCompGeneralSpin(list_1[j], isite2, sigma4, sigma3, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
          if (tmp_sgn == TRUE) {
            tmp_sgn = GetOffCompGeneralSpin(tmp_off, isite1, sigma2, sigma1, &tmp_off2, X->Def.SiteToBit, X->Def.Tpow);
            if (tmp_sgn == TRUE) {
              ConvertToList1GeneralSpin(tmp_off2, ihfbit, &off);
              dmv = tmp_v1[j] * tmp_V;
              if (X->Large.mode == M_MLTPLY) { // for multply
                tmp_v0[off] += dmv;
              }
              dam_pr += conj(tmp_v1[off]) * dmv;
            }
          }/*if (tmp_sgn == TRUE)*/
        }/*for (j = 1; j <= i_max; j++)*/
        X->Large.prdct += dam_pr;
      }/*for (ihermite = 0; ihermite < 2; ihermite++)*/
      StopTimer(413);
    }
  }/*for (i = 0; i < X->Def.NInterAll_OffDiagonal; i += 2)*/
  StopTimer(410);
  StopTimer(400);
  return 0;  
}/*int mltplyGeneralSpin*/
/**
@brief Driver function for Spin hamiltonian
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltplySpinGC(
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1//!<[in] Input producted vector
) {
  int iret=0;
  if (X->Def.iFlgGeneralSpin == FALSE) 
    iret = mltplyHalfSpinGC(X, tmp_v0, tmp_v1);
  else 
    iret = mltplyGeneralSpinGC(X, tmp_v0, tmp_v1);

  if(iret != 0) return iret;
  
  if(X->Boost.flgBoost == 1)
    iret = mltplySpinGCBoost(X, tmp_v0, tmp_v1);

  return iret;
}/*int mltplySpinGC*/
/**
@brief Driver function for Spin 1/2 Hamiltonian (grandcanonical)
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void mltplyHalfSpinGC_mini(
  struct BindStruct *X,//!<[inout]
  int site_i,
  int spin_i,
  int site_j,
  int spin_j,
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1//!<[in] Input producted vector
) {
  long unsigned int j;
  long unsigned int i;
  long unsigned int off = 0;
  long unsigned int is1_spin = 0;
  /**/
  long unsigned int isite1;
  long unsigned int org_isite1, org_isite2;
  long unsigned int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long unsigned int isA_up, isB_up;
  long unsigned int tmp_off = 0;
  double complex dam_pr;
  double complex tmp_trans;
  long int tmp_sgn;
  /*[s] For InterAll */
  double complex tmp_V;
  /*[e] For InterAll */

  long unsigned int i_max;
  i_max = X->Check.idim_max;

  int ihermite=0;
  int idx=0;

  //EDGeneralTransfer[i][0] -> site_i  
  //EDGeneralTransfer[i][1] -> spin_i  
  //EDGeneralTransfer[i][2] -> site_j  
  //EDGeneralTransfer[i][3] -> spin_j  
  //X->Large.mode == M_MLTPLY 
  org_isite1 = site_i+1;
  org_isite2 = site_j+1;
  org_sigma1 = spin_i;
  org_sigma2 = spin_j;
  dam_pr=0.0;
  if(org_isite1 == org_isite2){
    if(org_isite1 > X->Def.Nsite){
      if(org_sigma1==org_sigma2){  // longitudinal magnetic field
        dam_pr += child_GC_CisAis_spin_MPIdouble(org_isite1-1, org_sigma1, 1.0, X, tmp_v0, tmp_v1);
      }else{  // transverse magnetic field
        X->Large.mode = M_MLTPLY2;
        dam_pr += child_GC_CisAit_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, 1.0, X, tmp_v0, tmp_v1);
        X->Large.mode = M_MLTPLY;
      }
    }else{
      isite1 = X->Def.Tpow[org_isite1-1];
      if(org_sigma1==org_sigma2){
       // longitudinal magnetic field
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn) firstprivate(i_max, isite1, org_sigma1, X) shared(tmp_v0,tmp_v1)
        for(j=1;j<=i_max;j++){
          tmp_v0[j] += child_SpinGC_CisAis(j, X, isite1, org_sigma1)*tmp_v1[j]; 
        }
      }else{
        // transverse magnetic field
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, tmp_off) firstprivate(i_max, isite1, org_sigma2, X) shared(tmp_v0,tmp_v1)
        for(j=1;j <= i_max;j++){
          tmp_sgn  =  child_SpinGC_CisAit(j,X,isite1,org_sigma2,&tmp_off);
          if(tmp_sgn !=0){
            tmp_v0[tmp_off+1]  +=  tmp_sgn*tmp_v1[j];
          }
        }
      }
    }
  }//else{
     // hopping is not allowed in localized spin system
     //dam_pr=0.0;
  //}
  //return 0;
}/*int mltplyHalfSpinGC*/



/**
@brief Driver function for Spin 1/2 Hamiltonian (grandcanonical)
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltplyHalfSpinGC(
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1//!<[in] Input producted vector
) {
  long unsigned int j;
  long unsigned int i;
  long unsigned int off = 0;
  long unsigned int is1_spin = 0;
  long unsigned int isite1, isite2, sigma1, sigma2;
  long unsigned int sigma3, sigma4;
  double complex dam_pr;
  double complex tmp_trans;
  long int tmp_sgn;
  /*[s] For InterAll */
  double complex tmp_V;
  /*[e] For InterAll */

  long unsigned int i_max;
  i_max = X->Check.idim_max;

  int ihermite=0;
  int idx=0;

  StartTimer(500);

  StartTimer(510);
  for (i = 0; i < X->Def.EDNTransfer; i+=2 ) {
    if(X->Def.EDGeneralTransfer[i][0]+1 > X->Def.Nsite){
      dam_pr=0;
      if(X->Def.EDGeneralTransfer[i][1]==X->Def.EDGeneralTransfer[i][3]){
        fprintf(stderr, "Transverse_OffDiagonal component is illegal.\n");
      }
      else{
        StartTimer(511);
        dam_pr += child_GC_CisAit_spin_MPIdouble(
          X->Def.EDGeneralTransfer[i][0], X->Def.EDGeneralTransfer[i][1], 
          X->Def.EDGeneralTransfer[i][3], -X->Def.EDParaGeneralTransfer[i], 
          X, tmp_v0, tmp_v1);
        StopTimer(511);
      }
    }/*if(X->Def.EDGeneralTransfer[i][0]+1 > X->Def.Nsite)*/
    else{
      StartTimer(512);
      dam_pr = 0;
      for(ihermite=0; ihermite<2; ihermite++){
        idx=i+ihermite;
        isite1 = X->Def.EDGeneralTransfer[idx][0] + 1;
        isite2 = X->Def.EDGeneralTransfer[idx][2] + 1;
        sigma1 = X->Def.EDGeneralTransfer[idx][1];
        sigma2 = X->Def.EDGeneralTransfer[idx][3];
        tmp_trans = -X->Def.EDParaGeneralTransfer[idx];
        if (general_hopp_GetInfo(X, isite1, isite2, sigma1, sigma2) != 0) {
          return -1;
        }
       
        if(sigma1==sigma2){
          fprintf(stderr, "Transverse_OffDiagonal component is illegal.\n");
        }
        else{
          // longitudinal magnetic field (considerd in diagonalcalc.c)
          // transverse magnetic field
          is1_spin = X->Def.Tpow[isite1 - 1];
#pragma omp parallel for default(none) reduction(+:dam_pr) \
private(j, tmp_sgn) firstprivate(i_max, is1_spin, sigma2, X,off, tmp_trans) \
shared(tmp_v0, tmp_v1)
          for (j = 1; j <= i_max; j++) {
            tmp_sgn = child_SpinGC_CisAit(j, X, is1_spin, sigma2, &off);
            if(tmp_sgn !=0){
              tmp_v0[off+1] += tmp_v1[j]*tmp_trans;
              dam_pr += tmp_trans * conj(tmp_v1[off + 1]) * tmp_v1[j];
            }/*if(tmp_sgn !=0)*/
          }/*for (j = 1; j <= i_max; j++)*/
        }//sigma1 != sigma2
      }/*for(ihermite=0; ihermite<2; ihermite++)*/
      StopTimer(512);
    }
    X->Large.prdct += dam_pr;
  }/*for (i = 0; i < X->Def.EDNTransfer; i+=2 )*/
  StopTimer(510);
  /**
  InterAll 
  */
  StartTimer(520);
  for (i = 0; i < X->Def.NInterAll_OffDiagonal; i += 2) {
    if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite &&
        X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
      StartTimer(521);
      GC_general_int_spin_MPIdouble(i, X, tmp_v0, tmp_v1);
      StopTimer(521);
    }
    else if (X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
      StartTimer(522);
      GC_general_int_spin_MPIsingle(i, X, tmp_v0, tmp_v1);
      StopTimer(522);
    }
    else if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite) {
      StartTimer(522);
      GC_general_int_spin_MPIsingle(i + 1, X, tmp_v0, tmp_v1);
      StopTimer(522);
    }
    else {
      StartTimer(523);
      for (ihermite = 0; ihermite < 2; ihermite++) {
        idx = i + ihermite;
        isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
        isite2 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
        sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
        sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
        sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
        sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
        tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];
        general_int_spin_GetInfo(X, isite1, isite2, sigma1, sigma2, sigma3, sigma4, tmp_V);
        dam_pr = GC_general_int_spin(tmp_v0, tmp_v1, X);
        X->Large.prdct += dam_pr;
      }
      StopTimer(523);
    }
  }/*for (i = 0; i < X->Def.NInterAll_OffDiagonal; i += 2)*/
  StopTimer(520);
  /**
  Exchange
  */
  StartTimer(530);
  for (i = 0; i < X->Def.NExchangeCoupling; i++) {
    sigma1=0; sigma2=1;
    if (X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite &&
        X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite){
      StartTimer(531);
      dam_pr = child_GC_CisAitCiuAiv_spin_MPIdouble(
        X->Def.ExchangeCoupling[i][0], sigma1, sigma2, 
        X->Def.ExchangeCoupling[i][1], sigma2, sigma1, 
        X->Def.ParaExchangeCoupling[i], X, tmp_v0, tmp_v1);
      StopTimer(531);
    }
    else if (X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite) {
      StartTimer(532);
      dam_pr=child_GC_CisAitCiuAiv_spin_MPIsingle(
        X->Def.ExchangeCoupling[i][0], sigma1, sigma2,
        X->Def.ExchangeCoupling[i][1], sigma2, sigma1,
        X->Def.ParaExchangeCoupling[i], X, tmp_v0, tmp_v1);
      StopTimer(532);
    }
    else if (X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite) {
      StartTimer(532);
      dam_pr=child_GC_CisAitCiuAiv_spin_MPIsingle(
        X->Def.ExchangeCoupling[i][1], sigma2, sigma1,
        X->Def.ExchangeCoupling[i][0], sigma1, sigma2,
        conj(X->Def.ParaExchangeCoupling[i]), X, tmp_v0, tmp_v1);
      StopTimer(532);
    }
    else {
      StartTimer(533);
      exchange_spin_GetInfo(i, X);
      dam_pr = GC_exchange_spin(tmp_v0, tmp_v1, X);
      StopTimer(533);
    }
    X->Large.prdct += dam_pr;
  }/* for (i = 0; i < X->Def.NExchangeCoupling; i ++) */
  StopTimer(530);
  /**
  PairLift
  */
  StartTimer(540);
  for (i = 0; i < X->Def.NPairLiftCoupling; i++) {
    sigma1 =0; sigma2=1;
    if (X->Def.PairLiftCoupling[i][0] + 1 > X->Def.Nsite &&
        X->Def.PairLiftCoupling[i][1] + 1 > X->Def.Nsite) {
      StartTimer(541);
      dam_pr = child_GC_CisAitCiuAiv_spin_MPIdouble(
        X->Def.PairLiftCoupling[i][0], sigma1, sigma2, 
        X->Def.PairLiftCoupling[i][1], sigma1, sigma2,
        X->Def.ParaPairLiftCoupling[i], X, tmp_v0, tmp_v1);
      StopTimer(541);
    }
    else if (X->Def.PairLiftCoupling[i][1] + 1 > X->Def.Nsite) {
      StartTimer(542);
      dam_pr = child_GC_CisAitCiuAiv_spin_MPIsingle(
        X->Def.PairLiftCoupling[i][0], sigma1, sigma2, 
        X->Def.PairLiftCoupling[i][1], sigma1, sigma2, 
        X->Def.ParaPairLiftCoupling[i], X, tmp_v0, tmp_v1);
      StopTimer(542);
    }
    else if (X->Def.PairLiftCoupling[i][0] + 1 > X->Def.Nsite) {
      StartTimer(542);
      dam_pr = child_GC_CisAitCiuAiv_spin_MPIsingle(
        X->Def.PairLiftCoupling[i][1], sigma1, sigma2,
        X->Def.PairLiftCoupling[i][0], sigma1, sigma2,
        conj(X->Def.ParaPairLiftCoupling[i]), X, tmp_v0, tmp_v1);
      StopTimer(542);
    }
    else {
      StartTimer(543);
      pairlift_spin_GetInfo(i, X);
      dam_pr = GC_pairlift_spin(tmp_v0, tmp_v1, X);
      StopTimer(543);
    }
    X->Large.prdct += dam_pr;
  }/*for (i = 0; i < X->Def.NPairLiftCoupling; i += 2)*/
  StopTimer(540);

  StopTimer(500);
  return 0;
}/*int mltplyHalfSpinGC*/
/**
@brief Driver function for General Spin hamiltonian (grandcanonical)
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltplyGeneralSpinGC(
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1//!<[in] Input producted vector
) {
  long unsigned int j;
  long unsigned int i;
  long unsigned int off = 0;
  long unsigned int tmp_off = 0;
  long unsigned int isite1, isite2, sigma1, sigma2;
  long unsigned int sigma3, sigma4;
  double complex dam_pr;
  double complex tmp_trans;
  long int tmp_sgn;
  double num1 = 0;
  /*[s] For InterAll */
  double complex tmp_V;
  double complex dmv=0;
  /*[e] For InterAll */

  long unsigned int i_max;
  i_max = X->Check.idim_max;

  int ihermite=0;
  int idx=0;

  StartTimer(500);

  StartTimer(510);
  for (i = 0; i < X->Def.EDNTransfer; i += 2) {
    isite1 = X->Def.EDGeneralTransfer[i][0] + 1;
    isite2 = X->Def.EDGeneralTransfer[i][2] + 1;
    sigma1 = X->Def.EDGeneralTransfer[i][1];
    sigma2 = X->Def.EDGeneralTransfer[i][3];
    tmp_trans = -X->Def.EDParaGeneralTransfer[idx];
    dam_pr = 0.0;
    if (isite1 == isite2) {
      if (sigma1 != sigma2) {
        if (isite1 > X->Def.Nsite) {
          dam_pr = child_GC_CisAit_GeneralSpin_MPIdouble(
            isite1 - 1, sigma1, sigma2, tmp_trans, X, tmp_v0, tmp_v1);
          X->Large.prdct += dam_pr;
        }/*if (isite1 > X->Def.Nsite)*/
        else {
          for (ihermite = 0; ihermite<2; ihermite++) {
            idx = i + ihermite;
            isite1 = X->Def.EDGeneralTransfer[idx][0] + 1;
            isite2 = X->Def.EDGeneralTransfer[idx][2] + 1;
            sigma1 = X->Def.EDGeneralTransfer[idx][1];
            sigma2 = X->Def.EDGeneralTransfer[idx][3];
            tmp_trans = -X->Def.EDParaGeneralTransfer[idx];
    
            // transverse magnetic field
            dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) \
private(j, tmp_sgn, num1) firstprivate(i_max, isite1, sigma1, sigma2, X, off, tmp_trans) \
shared(tmp_v0, tmp_v1)
            for (j = 1; j <= i_max; j++) {
              num1 = GetOffCompGeneralSpin(
                j - 1, isite1, sigma2, sigma1, &off, X->Def.SiteToBit, X->Def.Tpow);
              if (num1 != 0) { // for multply
                tmp_v0[off + 1] += tmp_v1[j] * tmp_trans;
                dam_pr += conj(tmp_v1[off + 1]) * tmp_v1[j] * tmp_trans;
              }/*if (num1 != 0)*/
            }/*for (j = 1; j <= i_max; j++)*/
            X->Large.prdct += dam_pr;
          }/*for (ihermite = 0; ihermite<2; ihermite++)*/
        }
      }// sigma1 != sigma2          
      else{ // sigma1 = sigma2
        fprintf(stderr, "Error: Transverse_Diagonal component must be absorbed !");
      }
    }//isite1 = isite2
    else { // isite1 != isite2
      // hopping is not allowed in localized spin system
      return -1;
    }
  }/*for (i = 0; i < X->Def.EDNTransfer; i += 2)*/
  StopTimer(510);
  /**
  InterAll
  */
  StartTimer(520);
  for (i = 0; i< X->Def.NInterAll_OffDiagonal; i += 2) {
    if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite &&
        X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
      StartTimer(521);
      GC_general_int_GeneralSpin_MPIdouble(i, X, tmp_v0, tmp_v1);
      StopTimer(521);
    }
    else if (X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
      StartTimer(522);
      GC_general_int_GeneralSpin_MPIsingle(i, X, tmp_v0, tmp_v1);
      StopTimer(522);
    }
    else if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite) {
      StartTimer(522);
      GC_general_int_GeneralSpin_MPIsingle(i + 1, X, tmp_v0, tmp_v1);
      StopTimer(522);
    }
    else {
      StartTimer(523);
      for (ihermite = 0; ihermite < 2; ihermite++) {
        idx = i + ihermite;
        isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
        isite2 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
        sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
        sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
        sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
        sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
        tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];

        dam_pr = 0.0;
        if (sigma1 == sigma2) {
          if (sigma3 == sigma4) {
            fprintf(stderr, "InterAll_OffDiagonal component is illegal.\n");
            return -1;
          }/*if (sigma3 == sigma4)*/
          else {
            //sigma3=sigma4 term is considerd as a diagonal term.
#pragma omp parallel for default(none) reduction(+:dam_pr)              \
  private(j, tmp_sgn, dmv, off)                                         \
  firstprivate(i_max, isite1, isite2, sigma1, sigma3, sigma4, X, tmp_V) \
  shared(tmp_v0, tmp_v1)
            for (j = 1; j <= i_max; j++) {
              tmp_sgn = GetOffCompGeneralSpin(
                j - 1, isite2, sigma4, sigma3, &off, X->Def.SiteToBit, X->Def.Tpow);
              if (tmp_sgn == TRUE) {
                tmp_sgn = BitCheckGeneral(off, isite1, sigma1, X->Def.SiteToBit, X->Def.Tpow);
                if (tmp_sgn == TRUE) {
                  dmv = tmp_v1[j] * tmp_V;
                  if (X->Large.mode == M_MLTPLY) { // for multply
                    tmp_v0[off + 1] += dmv;
                  }
                  dam_pr += conj(tmp_v1[off + 1]) * dmv;
                }/*if (tmp_sgn == TRUE)*/
              }/*if (tmp_sgn == TRUE)*/
            }/*for (j = 1; j <= i_max; j++)*/
          }
        }/*if (sigma1 == sigma2)*/
        else if (sigma3 == sigma4) {
          //sigma1=sigma2 term is considerd as a diagonal term.
#pragma omp parallel for default(none) reduction(+:dam_pr)              \
  private(j, tmp_sgn, dmv, off, tmp_off)                                \
  firstprivate(i_max, isite1, isite2, sigma1, sigma2, sigma3, sigma4, X, tmp_V) \
  shared(tmp_v0, tmp_v1)
          for (j = 1; j <= i_max; j++) {
            tmp_sgn = BitCheckGeneral(j - 1, isite2, sigma3, X->Def.SiteToBit, X->Def.Tpow);
            if (tmp_sgn == TRUE) {
              tmp_sgn = GetOffCompGeneralSpin(
                j - 1, isite1, sigma2, sigma1, &off, X->Def.SiteToBit, X->Def.Tpow);
              if (tmp_sgn == TRUE) {
                dmv = tmp_v1[j] * tmp_V;
                if (X->Large.mode == M_MLTPLY) { // for multply
                  tmp_v0[off + 1] += dmv;
                }
                dam_pr += conj(tmp_v1[off + 1]) * dmv;
              }/*if (tmp_sgn == TRUE)*/
            }/*if (tmp_sgn == TRUE)*/
          }/*for (j = 1; j <= i_max; j++)*/
        }/*else if (sigma3 == sigma4)*/
        else {
#pragma omp parallel for default(none) reduction(+:dam_pr)              \
  private(j, tmp_sgn, dmv, off, tmp_off)                                \
  firstprivate(i_max, isite1, isite2, sigma1, sigma2, sigma3, sigma4, X, tmp_V) \
  shared(tmp_v0, tmp_v1)
          for (j = 1; j <= i_max; j++) {
            tmp_sgn = GetOffCompGeneralSpin(
              j - 1, isite2, sigma4, sigma3, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
            if (tmp_sgn == TRUE) {
              tmp_sgn = GetOffCompGeneralSpin(
                tmp_off, isite1, sigma2, sigma1, &off, X->Def.SiteToBit, X->Def.Tpow);
              if (tmp_sgn == TRUE) {
                dmv = tmp_v1[j] * tmp_V;
                if (X->Large.mode == M_MLTPLY) { // for multply
                  tmp_v0[off + 1] += dmv;
                }
                dam_pr += conj(tmp_v1[off + 1]) * dmv;
              }/*if (tmp_sgn == TRUE)*/
            }/*if (tmp_sgn == TRUE)*/
          }/*for (j = 1; j <= i_max; j++)*/
        }
        X->Large.prdct += dam_pr;
      }
      StopTimer(523);
    }
  }/*for (i = 0; i< X->Def.NInterAll_OffDiagonal; i += 2)*/
  StopTimer(520);

  StopTimer(500);
  return 0;
}/*int mltplyGeneralSpinGC*/
/**
@brief Driver function for Spin hamiltonian (Boost)
@return error code
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int mltplySpinGCBoost(
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1//!<[in] Input producted vector
)
{
  long unsigned int j;

  double complex dam_pr;

  /* SpinGCBoost */
  double complex* tmp_v2;
  double complex* tmp_v3;
  /* SpinGCBoost */
  
  long unsigned int i_max;
  i_max = X->Check.idim_max;

  StartTimer(500);
  tmp_v2 = cd_1d_allocate(i_max+1);
  tmp_v3 = cd_1d_allocate(i_max+1);

  general_int_spin_MPIBoost(X, tmp_v0, tmp_v1, tmp_v2, tmp_v3);
  dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) \
private(j) shared(tmp_v1,tmp_v0) firstprivate(i_max) 
  for(j=1;j<=i_max;j++)
    dam_pr += conj(tmp_v1[j])*tmp_v0[j]; // <H>=<v1|H|v1>
  X->Large.prdct += dam_pr;  
    
  /* SpinGCBoost */
  free_cd_1d_allocate(tmp_v2);
  free_cd_1d_allocate(tmp_v3);
  /* SpinGCBoost */
  StopTimer(500);
  return 0;
}/*int mltplySpinGCBoost*/

/******************************************************************************/
//[s] child functions
/******************************************************************************/

/**
@brief Compute exchange term of spin Hamiltonian (canonical)
@return Fragment of @f$\langle v_1|{\hat H}|v_1\rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex exchange_spin(
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1,//!<[in] Input producted vector
  struct BindStruct *X//!<[inout]
) {
  long unsigned int j;
  long unsigned int i_max = X->Large.i_max;
  long unsigned int off = 0;
  double complex dam_pr = 0;

#pragma omp parallel for default(none) reduction(+:dam_pr) \
firstprivate(i_max, X,off) private(j) shared(tmp_v0, tmp_v1)
  for (j = 1; j <= i_max; j++) 
    dam_pr += exchange_spin_element(j, tmp_v0, tmp_v1, X, &off);
  return dam_pr;
}/*double complex exchange_spin*/
/**
@brief Compute exchange term of spin Hamiltonian (grandcanonical)
@return Fragment of @f$\langle v_1|{\hat H}|v_1\rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex GC_exchange_spin(
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1,//!<[in] Input producted vector
  struct BindStruct *X//!<[inout]
) {
  long unsigned int j;
  long unsigned int i_max = X->Large.i_max;
  long unsigned int off = 0;
  double complex dam_pr = 0;

#pragma omp parallel for default(none) reduction(+:dam_pr) \
firstprivate(i_max, X,off) private(j) shared(tmp_v0, tmp_v1)
  for (j = 1; j <= i_max; j++)
    dam_pr += GC_exchange_spin_element(j, tmp_v0, tmp_v1, X, &off);
  return dam_pr;
}/*double complex GC_exchange_spin*/
/**
@brief Compute pair-lift term of spin Hamiltonian (grandcanonical)
@return Fragment of @f$\langle v_1|{\hat H}|v_1\rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex GC_pairlift_spin(
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1,//!<[in] Input producted vector
  struct BindStruct *X//!<[inout]
) {
  long unsigned int j;
  long unsigned int i_max = X->Large.i_max;
  long unsigned int off = 0;
  double complex dam_pr = 0;

#pragma omp parallel for default(none) reduction(+:dam_pr) \
firstprivate(i_max, X,off) private(j) shared(tmp_v0, tmp_v1)
  for (j = 1; j <= i_max; j++) 
    dam_pr += GC_pairlift_spin_element(j, tmp_v0, tmp_v1, X, &off);
  return dam_pr;
}/*double complex GC_pairlift_spin*/
/**
@brief Compute Inter-All term of spin Hamiltonian (canonical)
@return Fragment of @f$\langle v_1|{\hat H}|v_1\rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex general_int_spin(
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1,//!<[in] Input producted vector
  struct BindStruct *X//!<[inout]
) {
  double complex dam_pr, tmp_V, dmv;
  long unsigned int j, i_max;
  long unsigned int org_sigma2, org_sigma4;
  long unsigned int isA_up, isB_up;
  long unsigned int tmp_off = 0;
  int tmp_sgn;

  i_max = X->Large.i_max;
  org_sigma2 = X->Large.is2_spin;
  org_sigma4 = X->Large.is4_spin;
  tmp_V = X->Large.tmp_V;
  isA_up = X->Large.is1_up;
  isB_up = X->Large.is2_up;
  dam_pr = 0.0;

#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(tmp_v1, tmp_v0)
  for (j = 1; j <= i_max; j++) {
    tmp_sgn = child_exchange_spin_element(j, X, isA_up, isB_up, org_sigma2, org_sigma4, &tmp_off);
    if (tmp_sgn != 0) {
      dmv = tmp_v1[j] * tmp_sgn * tmp_V;
      tmp_v0[tmp_off] += dmv;
      dam_pr += conj(tmp_v1[tmp_off]) * dmv;
    }/*if (tmp_sgn != 0)*/
  }/*for (j = 1; j <= i_max; j++)*/
  return dam_pr;
}/*double complex general_int_spin*/
/**
@brief Compute Inter-All term of spin Hamiltonian (grandcanonical)
@return Fragment of @f$\langle v_1|{\hat H}|v_1\rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex GC_general_int_spin(
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1,//!<[in] Input producted vector
  struct BindStruct *X//!<[inout]
) {
  double complex dam_pr, tmp_V;
  long unsigned int j, i_max;
  long unsigned int org_isite1, org_isite2;
  long unsigned int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long unsigned int isA_up, isB_up;
  long unsigned int tmp_off = 0;

  i_max = X->Large.i_max;
  org_isite1 = X->Large.isite1;
  org_isite2 = X->Large.isite2;
  org_sigma1 = X->Large.is1_spin;
  org_sigma2 = X->Large.is2_spin;
  org_sigma3 = X->Large.is3_spin;
  org_sigma4 = X->Large.is4_spin;
  tmp_V = X->Large.tmp_V;
  dam_pr = 0.0;
  isA_up = X->Def.Tpow[org_isite1 - 1];
  isB_up = X->Def.Tpow[org_isite2 - 1];

#pragma omp parallel default(none) reduction(+:dam_pr) \
private(j) shared(tmp_v0, tmp_v1) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma1,org_sigma2,org_sigma3,org_sigma4,tmp_off, tmp_V) 
  {
    if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp for
      for (j = 1; j <= i_max; j++)
        dam_pr += GC_CisAisCisAis_spin_element(
          j, isA_up, isB_up, org_sigma2, org_sigma4, tmp_V, tmp_v0, tmp_v1, X);
    }
    else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp for
      for (j = 1; j <= i_max; j++)
        dam_pr += GC_CisAisCitAiu_spin_element(
          j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, tmp_v0, tmp_v1, X, &tmp_off);
    }
    else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
#pragma omp for
      for (j = 1; j <= i_max; j++)
        dam_pr += GC_CisAitCiuAiu_spin_element(
          j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, tmp_v0, tmp_v1, X, &tmp_off);
    }
    else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp for
      for (j = 1; j <= i_max; j++)
        dam_pr += GC_CisAitCiuAiv_spin_element(
          j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, tmp_v0, tmp_v1, X, &tmp_off);
    }
  }/*End of parallel region*/
  return dam_pr;
}/*double complex GC_general_int_spin*/
/******************************************************************************/
//[e] child functions
/******************************************************************************/
