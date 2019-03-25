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
@brief Functions for spin Hamiltonian + MPI
*/

#include "Common.h"
#include "mltplyCommon.h"
#include "bitcalc.h"
#include "wrapperMPI.h"
#include "mltplyMPISpin.h"
#include "mltplyMPISpinCore.h"

/**
@brief Exchange term in Spin model
  When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void child_general_int_spin_MPIdouble(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[out] Result v0 = H v1
  double complex **tmp_v1//!<[in] v0 = H v1
){
  X_child_general_int_spin_MPIdouble(
    (int)X->Def.InterAll_OffDiagonal[i_int][0], (int)X->Def.InterAll_OffDiagonal[i_int][1],
    (int)X->Def.InterAll_OffDiagonal[i_int][3], (int)X->Def.InterAll_OffDiagonal[i_int][4],
    (int)X->Def.InterAll_OffDiagonal[i_int][5], (int)X->Def.InterAll_OffDiagonal[i_int][7],
    X->Def.ParaInterAll_OffDiagonal[i_int], X, nstate, tmp_v0, tmp_v1);
}/*void child_general_int_spin_MPIdouble*/
/**
@brief Exchange term in Spin model
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void X_child_general_int_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_J,//!<[in] Copupling constatnt
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex **tmp_v1//!<[in] Vector to be producted
) {
  int mask1, mask2, state1, state2, origin;
  unsigned long int idim_max_buf, j, ioff;
  double complex Jint;
  int one = 1;

  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  origin = myrank ^ (mask1 + mask2);

  state1 = (origin & mask1) / mask1;
  state2 = (origin & mask2) / mask2;

  if (state1 == org_ispin2 && state2 == org_ispin4) {
    Jint = tmp_J;
  }
  else if (state1 == org_ispin1 && state2 == org_ispin3) {
    Jint = conj(tmp_J);
    if (X->Large.mode == M_CORR || X->Large.mode == M_CALCSPEC) {
      Jint = 0;
    }
  }
  else return;

  idim_max_buf = SendRecv_i(origin, X->Check.idim_max);
  SendRecv_iv(origin, X->Check.idim_max + 1, idim_max_buf + 1, list_1, list_1buf);
  SendRecv_cv(origin, X->Check.idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);

#pragma omp parallel for default(none)  private(j, ioff) \
  firstprivate(idim_max_buf,Jint,X) shared(list_2_1,list_2_2,list_1buf,v1buf,tmp_v1,tmp_v0,nstate,one)
  for (j = 1; j <= idim_max_buf; j++) {
    GetOffComp(list_2_1, list_2_2, list_1buf[j],
        X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
    zaxpy_(&nstate, &Jint, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
  }/*for (j = 1; j <= idim_max_buf; j++)*/
}/*double complex X_child_general_int_spin_MPIdouble*/
/**
@brief Exchange term in Spin model
  When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void X_child_general_int_spin_TotalS_MPIdouble(
  int org_isite1,//!<[in] site 1
  int org_isite3,//!<[in] site 3
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex **tmp_v1//!<[in] Vector to be producted
){
  int mask1, mask2, num1_up, num2_up, origin, one = 1;
  unsigned long int idim_max_buf, j, ioff, ibit_tmp;
  double complex dmv;

  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  if (mask1 == mask2) origin = myrank ^ mask1;
  else origin = myrank ^ (mask1 + mask2);
  num1_up = (origin & mask1) / mask1;
  num2_up = (origin & mask2) / mask2;

  ibit_tmp = (num1_up) ^ (num2_up);
  if (ibit_tmp == 0) return;

  idim_max_buf = SendRecv_i(origin, X->Check.idim_max);
  SendRecv_iv(origin, X->Check.idim_max + 1, idim_max_buf + 1, list_1, list_1buf);
  SendRecv_cv(origin, X->Check.idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);

#pragma omp parallel for default(none)  private(j, dmv, ioff) \
    firstprivate(idim_max_buf, X) \
  shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0,nstate,one)
  for (j = 1; j <= idim_max_buf; j++) {
    GetOffComp(list_2_1, list_2_2, list_1buf[j],
      X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
    zaxpy_(&nstate, &dmv, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
  }/*for (j = 1; j <= idim_max_buf; j++)*/
  return;
}/*double complex X_child_general_int_spin_MPIdouble*/
/**
@brief Exchange term in Spin model
  When only site2 is in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void child_general_int_spin_MPIsingle(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[out] Result v0 = H v1
  double complex **tmp_v1//!<[in] v0 = H v1
){

  X_child_general_int_spin_MPIsingle(
    (int)X->Def.InterAll_OffDiagonal[i_int][0], (int)X->Def.InterAll_OffDiagonal[i_int][1], 
    (int)X->Def.InterAll_OffDiagonal[i_int][3], (int)X->Def.InterAll_OffDiagonal[i_int][4],
    (int)X->Def.InterAll_OffDiagonal[i_int][5], (int)X->Def.InterAll_OffDiagonal[i_int][7],
    X->Def.ParaInterAll_OffDiagonal[i_int], X, nstate, tmp_v0, tmp_v1);
}/*void child_general_int_spin_MPIsingle*/
/*
@brief General interaction term of canonical spin system.
site 3 is in the inter process region
*/
void X_child_general_int_spin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_J,//!<[in] Copupling constatnt
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex **tmp_v1//!<[in] Vector to be producted
) {
  int mask2, state2, origin;
  unsigned long int mask1, idim_max_buf, j, ioff, state1, jreal, state1check;
  double complex Jint;
  int one = 1;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[org_isite3];
  origin = myrank ^ mask2;
  state2 = (origin & mask2) / mask2;

  if (state2 == org_ispin4) {
    state1check = (unsigned long int) org_ispin2;
    Jint = tmp_J;
  }
  else if (state2 == org_ispin3) {
    state1check = (unsigned long int) org_ispin1;
    Jint = conj(tmp_J);
    if (X->Large.mode == M_CORR || X->Large.mode == M_CALCSPEC) {
      Jint = 0;
    }
  }
  else return;

  idim_max_buf = SendRecv_i(origin, X->Check.idim_max);
  SendRecv_iv(origin, X->Check.idim_max + 1, idim_max_buf + 1, list_1, list_1buf);
  SendRecv_cv(origin, X->Check.idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);
  /*
  Index in the intra PE
  */
  mask1 = X->Def.Tpow[org_isite1];
  //printf("debug1 %ld\n", idim_max_buf);

#pragma omp parallel for default(none)  private(j, ioff, jreal, state1) \
firstprivate(idim_max_buf, Jint, X, mask1, state1check, org_isite1) \
shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0,nstate,one)
  for (j = 1; j <= idim_max_buf; j++) {
    jreal = list_1buf[j];

    state1 = (jreal & mask1) / mask1;
    if (state1 == state1check) {
      GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
        X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
      zaxpy_(&nstate, &Jint, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
    }
  }
}/*double complex X_child_general_int_spin_MPIsingle*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_child_general_int_spin_MPIdouble(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[out] Result v0 = H v1
  double complex **tmp_v1//!<[in] v0 = H v1
){
  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
      X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
    GC_child_CisAisCjuAjv_spin_MPIdouble(i_int, X, nstate, tmp_v0, tmp_v1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
    GC_child_CisAitCjuAju_spin_MPIdouble(i_int, X, nstate, tmp_v0, tmp_v1);
  }
  else {
    GC_child_CisAitCiuAiv_spin_MPIdouble(i_int, X, nstate, tmp_v0, tmp_v1);
  }
}/*void GC_child_general_int_spin_MPIdouble*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_child_general_int_spin_MPIsingle(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[out] Result v0 = H v1
  double complex **tmp_v1//!<[in] v0 = H v1
){
  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
      X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
    GC_child_CisAisCjuAjv_spin_MPIsingle(i_int, X, nstate, tmp_v0, tmp_v1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
    GC_child_CisAitCjuAju_spin_MPIsingle(i_int, X, nstate, tmp_v0, tmp_v1);
  }
  else {
    GC_child_CisAitCiuAiv_spin_MPIsingle(i_int, X, nstate, tmp_v0, tmp_v1);
  }
}/*void GC_child_general_int_spin_MPIsingle*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_child_general_int_GeneralSpin_MPIdouble(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[out] Result v0 = H v1
  double complex **tmp_v1//!<[in] v0 = H v1
){

  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
      X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
    X_GC_child_CisAisCjuAjv_GeneralSpin_MPIdouble(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
      X->Def.InterAll_OffDiagonal[i_int][4], X->Def.InterAll_OffDiagonal[i_int][5],
      X->Def.InterAll_OffDiagonal[i_int][7], X->Def.ParaInterAll_OffDiagonal[i_int], X, nstate, tmp_v0, tmp_v1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
    X_GC_child_CisAitCjuAju_GeneralSpin_MPIdouble(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
      X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
      X->Def.InterAll_OffDiagonal[i_int][5], X->Def.ParaInterAll_OffDiagonal[i_int], X, nstate, tmp_v0, tmp_v1);
  }
  else {
    X_GC_child_CisAitCjuAjv_GeneralSpin_MPIdouble(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1], 
      X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
      X->Def.InterAll_OffDiagonal[i_int][5], X->Def.InterAll_OffDiagonal[i_int][7],
      X->Def.ParaInterAll_OffDiagonal[i_int], X, nstate, tmp_v0, tmp_v1);
  }
}/*void GC_child_general_int_spin_MPIdouble*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_child_general_int_GeneralSpin_MPIsingle(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[out] Result v0 = H v1
  double complex **tmp_v1//!<[in] v0 = H v1
){

  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
      X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
    X_GC_child_CisAisCjuAjv_GeneralSpin_MPIsingle(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
      X->Def.InterAll_OffDiagonal[i_int][4], X->Def.InterAll_OffDiagonal[i_int][5],
      X->Def.InterAll_OffDiagonal[i_int][7], X->Def.ParaInterAll_OffDiagonal[i_int], X, nstate, tmp_v0, tmp_v1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
    X_GC_child_CisAitCjuAju_GeneralSpin_MPIsingle(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1], 
      X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
      X->Def.InterAll_OffDiagonal[i_int][5], X->Def.ParaInterAll_OffDiagonal[i_int], X, nstate, tmp_v0, tmp_v1);
    }
  else {
    X_GC_child_CisAitCjuAjv_GeneralSpin_MPIsingle(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1], 
      X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
      X->Def.InterAll_OffDiagonal[i_int][5], X->Def.InterAll_OffDiagonal[i_int][7],
      X->Def.ParaInterAll_OffDiagonal[i_int], X, nstate, tmp_v0, tmp_v1);
  }

}/*void GC_child_general_int_spin_MPIsingle*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void child_general_int_GeneralSpin_MPIdouble(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[out] Result v0 = H v1
  double complex **tmp_v1//!<[in] v0 = H v1
){
    X_child_CisAitCjuAjv_GeneralSpin_MPIdouble(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1], 
      X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
      X->Def.InterAll_OffDiagonal[i_int][5], X->Def.InterAll_OffDiagonal[i_int][7],
      X->Def.ParaInterAll_OffDiagonal[i_int], X, nstate, tmp_v0, tmp_v1);

}/*void GC_child_general_int_spin_MPIdouble*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void child_general_int_GeneralSpin_MPIsingle(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[out] Result v0 = H v1
  double complex **tmp_v1//!<[in] v0 = H v1
){

  X_child_CisAitCjuAjv_GeneralSpin_MPIsingle(
    X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
    X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
    X->Def.InterAll_OffDiagonal[i_int][5], X->Def.InterAll_OffDiagonal[i_int][7],
    X->Def.ParaInterAll_OffDiagonal[i_int], X, nstate, tmp_v0, tmp_v1);

}/*void GC_child_general_int_spin_MPIsingle*/
