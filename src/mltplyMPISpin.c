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

#ifdef MPI
#include "mpi.h"
#endif
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
void general_int_spin_MPIdouble(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
#ifdef MPI
  double complex dam_pr = 0;
  dam_pr = child_general_int_spin_MPIdouble(
    (int)X->Def.InterAll_OffDiagonal[i_int][0], (int)X->Def.InterAll_OffDiagonal[i_int][1],
    (int)X->Def.InterAll_OffDiagonal[i_int][3], (int)X->Def.InterAll_OffDiagonal[i_int][4],
    (int)X->Def.InterAll_OffDiagonal[i_int][5], (int)X->Def.InterAll_OffDiagonal[i_int][7],
    X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  /**
  Add @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
  to LargeList::prdct
  */
  X->Large.prdct += dam_pr;
#endif
}/*void general_int_spin_MPIdouble*/
/**
@brief Exchange term in Spin model
 When both site1 and site2 are in the inter process region.
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double complex child_general_int_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_J,//!<[in] Copupling constatnt
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex *tmp_v1//!<[in] Vector to be producted
) {
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin;
  unsigned long int idim_max_buf, j, ioff;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;

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
  else return 0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                      list_1buf,   idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                       v1buf,      idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
  if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff) \
  firstprivate(idim_max_buf, Jint, X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {
      GetOffComp(list_2_1, list_2_2, list_1buf[j],
        X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
      dmv = Jint * v1buf[j];
      tmp_v0[ioff] += dmv;
      dam_pr += conj(tmp_v1[ioff]) * dmv;
    }/*for (j = 1; j <= idim_max_buf; j++)*/
  }/*if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC)*/
  else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff) \
  firstprivate(idim_max_buf, Jint, X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {
      GetOffComp(list_2_1, list_2_2, list_1buf[j],
        X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
      dmv = Jint * v1buf[j];
      dam_pr += conj(tmp_v1[ioff]) * dmv;
    }/*for (j = 1; j <= idim_max_buf; j++)*/
  }/*if (! (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC))*/
  return dam_pr;
#else
  return 0.0;
#endif
}/*double complex child_general_int_spin_MPIdouble*/
/**
@brief Exchange term in Spin model
  When both site1 and site2 are in the inter process region.
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double complex child_general_int_spin_TotalS_MPIdouble(
  int org_isite1,//!<[in] site 1
  int org_isite3,//!<[in] site 3
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex *tmp_v1//!<[in] Vector to be producted
){
#ifdef MPI
  int mask1, mask2, num1_up, num2_up, ierr, origin;
  unsigned long int idim_max_buf, j, ioff, ibit_tmp;
  MPI_Status statusMPI;
  double complex dmv, dam_pr;

  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  if (mask1 == mask2) origin = myrank ^ mask1;
  else origin = myrank ^ (mask1 + mask2);
  num1_up = (origin & mask1) / mask1;
  num2_up = (origin & mask2) / mask2;

  ibit_tmp = (num1_up) ^ (num2_up);
  if (ibit_tmp == 0) return 0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                      list_1buf,   idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,       idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff) \
    firstprivate(idim_max_buf,  X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
  for (j = 1; j <= idim_max_buf; j++) {
    GetOffComp(list_2_1, list_2_2, list_1buf[j],
      X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
    dmv = 0.5 * v1buf[j];
    dam_pr += conj(tmp_v1[ioff]) * dmv;
  }/*for (j = 1; j <= idim_max_buf; j++)*/
  return dam_pr;
#else
  return 0.0;
#endif
}/*double complex child_general_int_spin_MPIdouble*/
/**
@brief Exchange term in Spin model
  When only site2 is in the inter process region.
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void general_int_spin_MPIsingle(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
#ifdef MPI
  double complex dam_pr = 0;

  dam_pr = child_general_int_spin_MPIsingle(
    (int)X->Def.InterAll_OffDiagonal[i_int][0], (int)X->Def.InterAll_OffDiagonal[i_int][1], 
    (int)X->Def.InterAll_OffDiagonal[i_int][3], (int)X->Def.InterAll_OffDiagonal[i_int][4],
    (int)X->Def.InterAll_OffDiagonal[i_int][5], (int)X->Def.InterAll_OffDiagonal[i_int][7],
    X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  /**
  Add @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
  to LargeList::prdct
  */
  X->Large.prdct += dam_pr;
#endif
}/*void general_int_spin_MPIsingle*/
/*
@brief General interaction term of canonical spin system.
site 3 is in the inter process region
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
*/
double complex child_general_int_spin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_J,//!<[in] Copupling constatnt
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex *tmp_v1//!<[in] Vector to be producted
) {
#ifdef MPI
  int mask2, state2, ierr, origin;
  unsigned long int mask1, idim_max_buf, j, ioff, state1, jreal, state1check;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;
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
  else return 0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                      list_1buf,   idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,       idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  /*
  Index in the intra PE
  */
  mask1 = X->Def.Tpow[org_isite1];

  dam_pr = 0.0;
  if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff, jreal, state1) \
firstprivate(idim_max_buf, Jint, X, mask1, state1check, org_isite1) \
shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {

      jreal = list_1buf[j];

      state1 = (jreal & mask1) / mask1;
      if (state1 == state1check) {
        GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
          X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);

        dmv = Jint * v1buf[j];
        tmp_v0[ioff] += dmv;
        dam_pr += conj(tmp_v1[ioff]) * dmv;
      }
    }
  }
  else if (X->Large.mode == M_TOTALS) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff, jreal, state1) \
firstprivate(idim_max_buf, Jint, X, mask1, state1check, org_isite1) \
shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {

      jreal = list_1buf[j];

      state1 = (jreal & mask1) / mask1;
      if (state1 == state1check) {
        GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
          X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);

        dmv = Jint * v1buf[j];
        dmv = 0.5 * v1buf[j];
        dam_pr += conj(tmp_v1[ioff]) * dmv;
      }/*if (state1 == state1check)*/
    }/*for (j = 1; j <= idim_max_buf; j++)*/
  }/*if (X->Large.mode == M_TOTALS)*/
  else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff, jreal, state1) \
firstprivate(idim_max_buf, Jint, X, mask1, state1check, org_isite1) \
shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {

      jreal = list_1buf[j];

      state1 = (jreal & mask1) / mask1;
      if (state1 == state1check) {
        GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
          X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
        dmv = Jint * v1buf[j];
        dam_pr += conj(tmp_v1[ioff]) * dmv;
      }/*if (state1 == state1check)*/
    }/*for (j = 1; j <= idim_max_buf; j++)*/
  }
  return dam_pr;
#else
  return 0.0;
#endif
}/*double complex child_general_int_spin_MPIsingle*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_general_int_spin_MPIdouble(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
      X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
    GC_CisAisCjuAjv_spin_MPIdouble(i_int, X, tmp_v0, tmp_v1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
    GC_CisAitCjuAju_spin_MPIdouble(i_int, X, tmp_v0, tmp_v1);
  }
  else {
    GC_CisAitCiuAiv_spin_MPIdouble(i_int, X, tmp_v0, tmp_v1);
  }
}/*void GC_general_int_spin_MPIdouble*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_general_int_spin_MPIsingle(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
      X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
    GC_CisAisCjuAjv_spin_MPIsingle(i_int, X, tmp_v0, tmp_v1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
    GC_CisAitCjuAju_spin_MPIsingle(i_int, X, tmp_v0, tmp_v1);
  }
  else {
    GC_CisAitCiuAiv_spin_MPIsingle(i_int, X, tmp_v0, tmp_v1);
  }
}/*void GC_general_int_spin_MPIsingle*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_general_int_GeneralSpin_MPIdouble(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
#ifdef MPI
  double complex dam_pr;
  // MPI_Status statusMPI;

  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
      X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
    dam_pr = child_GC_CisAisCjuAjv_GeneralSpin_MPIdouble(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
      X->Def.InterAll_OffDiagonal[i_int][4], X->Def.InterAll_OffDiagonal[i_int][5],
      X->Def.InterAll_OffDiagonal[i_int][7], X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
    dam_pr = child_GC_CisAitCjuAju_GeneralSpin_MPIdouble(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
      X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
      X->Def.InterAll_OffDiagonal[i_int][5], X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  }
  else {
    dam_pr = child_GC_CisAitCjuAjv_GeneralSpin_MPIdouble(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1], 
      X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
      X->Def.InterAll_OffDiagonal[i_int][5], X->Def.InterAll_OffDiagonal[i_int][7],
      X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  }
  X->Large.prdct += dam_pr;
#endif
}/*void GC_general_int_spin_MPIdouble*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_general_int_GeneralSpin_MPIsingle(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
#ifdef MPI
  double complex dam_pr;

  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
      X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
    dam_pr = child_GC_CisAisCjuAjv_GeneralSpin_MPIsingle(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
      X->Def.InterAll_OffDiagonal[i_int][4], X->Def.InterAll_OffDiagonal[i_int][5],
      X->Def.InterAll_OffDiagonal[i_int][7], X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
    dam_pr = child_GC_CisAitCjuAju_GeneralSpin_MPIsingle(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1], 
      X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
      X->Def.InterAll_OffDiagonal[i_int][5], X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
    }
  else {
    dam_pr = child_GC_CisAitCjuAjv_GeneralSpin_MPIsingle(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1], 
      X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
      X->Def.InterAll_OffDiagonal[i_int][5], X->Def.InterAll_OffDiagonal[i_int][7],
      X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  }

  X->Large.prdct += dam_pr;
#endif
}/*void GC_general_int_spin_MPIsingle*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void general_int_GeneralSpin_MPIdouble(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
    double complex dam_pr;
    dam_pr = child_CisAitCjuAjv_GeneralSpin_MPIdouble(
      X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1], 
      X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
      X->Def.InterAll_OffDiagonal[i_int][5], X->Def.InterAll_OffDiagonal[i_int][7],
      X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
    X->Large.prdct += dam_pr;

}/*void GC_general_int_spin_MPIdouble*/
/**
@brief General interaction term in the Spin model + GC
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void general_int_GeneralSpin_MPIsingle(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
  double complex dam_pr;

  dam_pr = child_CisAitCjuAjv_GeneralSpin_MPIsingle(
    X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
    X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
    X->Def.InterAll_OffDiagonal[i_int][5], X->Def.InterAll_OffDiagonal[i_int][7],
    X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);

  X->Large.prdct += dam_pr;
}/*void GC_general_int_spin_MPIsingle*/
