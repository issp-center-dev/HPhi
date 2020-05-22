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
@brief Functions for spin Hamiltonian + MPI (Core)

General two body term:
<table>
  <tr>
    <td></td>
    <td>1/2 spin</td>
    <td>1/2 spin</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td></td>
    <td>MPI single</td>
    <td>MPI double</td>
    <td>MPI single</td>
    <td>MPI double</td>
  </tr>
  <tr>
    <td>@f$c_{is}^\dagger c_{is} c_{ju}^\dagger c_{ju}@f$</td>
    <td>::GC_CisAisCjuAju_spin_MPIsingle, ::child_GC_CisAisCjuAjv_spin_MPIsingle</td>
    <td>::GC_CisAisCjuAju_spin_MPIdouble, ::child_GC_CisAisCjuAjv_spin_MPIdouble</td>
    <td>::child_CisAisCjuAju_GeneralSpin_MPIsingle, ::child_GC_CisAisCjuAjv_GeneralSpin_MPIsingle</td>
    <td>::child_CisAisCjuAju_GeneralSpin_MPIsingle, ::child_GC_CisAisCjuAjv_GeneralSpin_MPIsingle</td>
  </tr>
  <tr>
    <td>@f$c_{is}^\dagger c_{is} c_{ju}^\dagger c_{jv}@f$</td>
    <td>::GC_CisAisCjuAjv_spin_MPIsingle, ::child_GC_CisAisCjuAjv_spin_MPIsingle</td>
    <td>::GC_CisAisCjuAjv_spin_MPIdouble, ::child_GC_CisAisCjuAjv_spin_MPIdouble</td>
    <td>::child_CisAisCjuAjv_GeneralSpin_MPIsingle, ::child_GC_CisAisCjuAjv_GeneralSpin_MPIsingle</td>
    <td>::child_CisAisCjuAjv_GeneralSpin_MPIsingle, ::child_GC_CisAisCjuAjv_GeneralSpin_MPIsingle</td>
  </tr>
  <tr>
    <td>@f$c_{is}^\dagger c_{it} c_{ju}^\dagger c_{ju}@f$</td>
    <td>::GC_CisAitCjuAju_spin_MPIsingle, ::child_GC_CisAisCjuAjv_spin_MPIsingle</td>
    <td>::GC_CisAitCjuAju_spin_MPIdouble, ::child_GC_CisAisCjuAjv_spin_MPIdouble</td>
    <td>::child_CisAitCjuAju_GeneralSpin_MPIsingle, ::child_GC_CisAitCjuAju_GeneralSpin_MPIsingle</td>
    <td>::child_CisAitCjuAju_GeneralSpin_MPIsingle, ::child_GC_CisAitCjuAju_GeneralSpin_MPIsingle</td>
  </tr>
  <tr>
    <td>@f$c_{is}^\dagger c_{it} c_{ju}^\dagger c_{jv}@f$</td>
    <td>::GC_CisAitCjuAjv_spin_MPIsingle, ::child_GC_CisAisCjuAjv_spin_MPIsingle</td>
    <td>::GC_CisAitCjuAjv_spin_MPIdouble, ::child_GC_CisAisCjuAjv_spin_MPIdouble</td>
    <td>::child_CisAitCjuAjv_GeneralSpin_MPIsingle, ::child_GC_CisAitCjuAjv_GeneralSpin_MPIsingle</td>
    <td>::child_CisAitCjuAjv_GeneralSpin_MPIsingle, ::child_GC_CisAitCjuAjv_GeneralSpin_MPIsingle</td>
  </tr>
</table>
*/
#ifdef MPI
#include "mpi.h"
#endif
#include "Common.h"
#include "mltplyCommon.h"
#include "mltplySpinCore.h"
#include "mltplyMPISpinCore.h"
#include "bitcalc.h"
#include "wrapperMPI.h"
/**
@brief Exchange and Pairlifting term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_CisAitCiuAiv_spin_MPIdouble(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  double complex dam_pr;  
  dam_pr =  child_GC_CisAitCiuAiv_spin_MPIdouble(
    X->Def.InterAll_OffDiagonal[i_int][0],  X->Def.InterAll_OffDiagonal[i_int][1], 
    X->Def.InterAll_OffDiagonal[i_int][3],  X->Def.InterAll_OffDiagonal[i_int][4], 
    X->Def.InterAll_OffDiagonal[i_int][5],  X->Def.InterAll_OffDiagonal[i_int][7],
    X->Def.ParaInterAll_OffDiagonal[i_int],X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}/*void GC_CisAitCiuAiv_spin_MPIdouble*/
/**
@brief @f$c_{is}^\dagger c_{it} c_{iu}^\dagger c_{iv}@f$ term
       in Spin model + GC.
       When both site1 and site2 are in the inter process region.
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double complex child_GC_CisAitCiuAiv_spin_MPIdouble(
  int org_isite1,//!<[in] site i
  int org_ispin1,//!<[in] spin s
  int org_ispin2,//!<[in] spin t
  int org_isite3,//!<[in] site i?
  int org_ispin3,//!<[in] spin u
  int org_ispin4,//!<[in] spin v
  double complex tmp_J,//!<[in] Copupling constatnt
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex *tmp_v1//!<[in] Vector to be producted
) {
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin;
  unsigned long int idim_max_buf, j;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;

  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  if (org_isite1 != org_isite3) {
    origin = myrank ^ (mask1 + mask2);
  }
  else {
    if (org_ispin1 == org_ispin4 && org_ispin2 == org_ispin3) { //CisAitCitAis=CisAis
      dam_pr = child_GC_CisAis_spin_MPIdouble(org_isite1, org_ispin1, tmp_J, X, tmp_v0, tmp_v1);
      return (dam_pr);
    }
    else { //CisAitCisAit=0
      return 0.0;
    }
  }

  state1 = (origin & mask1) / mask1;
  state2 = (origin & mask2) / mask2;

  if (state1 == org_ispin2 && state2 == org_ispin4) {
    Jint = tmp_J;
  }
  else if (state1 == org_ispin1 && state2 == org_ispin3) {
    Jint = conj(tmp_J);
    if (X->Large.mode == M_CORR ||X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC) {
      Jint = 0;
    }
  }
  else {
    return 0;
  }

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {
        dmv = Jint * v1buf[j];
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }else if(X->Large.mode == H_CORR){
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {
        dmv = Jint * v1buf[j];
        dam_pr += conj(tmp_v0[j]) * dmv;
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }else {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {
        dmv = Jint * v1buf[j];
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }
  }/*End of parallel region*/
  return dam_pr;
#else
 return 0.0;
#endif
}/*void GC_CisAitCiuAiv_spin_MPIdouble*/
/**
@brief Wrapper for calculating CisAisCjuAjv term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_CisAisCjuAjv_spin_MPIdouble(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/
){
#ifdef MPI
  double complex dam_pr;
  dam_pr = child_GC_CisAisCjuAjv_spin_MPIdouble(
    X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
    X->Def.InterAll_OffDiagonal[i_int][4], X->Def.InterAll_OffDiagonal[i_int][5],
    X->Def.InterAll_OffDiagonal[i_int][7], X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}/*void GC_CisAitCiuAiv_spin_MPIdouble*/
/**
@brief CisAisCjuAjv term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@return @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex child_GC_CisAisCjuAjv_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_J,//!<[in] Copupling constatnt
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex *tmp_v1//!<[in] Vector to be producted
) {
#ifdef MPI
  int mask1, mask2, state2, ierr;
  long int origin, num1;
  unsigned long int idim_max_buf, j;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;

  if (org_isite1 == org_isite3 && org_ispin1 == org_ispin4) {//CisAisCitAis
      return 0.0;
  }

  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  origin = myrank ^ mask2;
  state2 = (origin & mask2) / mask2;
  num1 = child_SpinGC_CisAis((unsigned long int) myrank + 1, X, mask1, org_ispin1);
  if (num1 != 0 && state2 == org_ispin4) {
    Jint = tmp_J;
  }
  else if (child_SpinGC_CisAis(origin + 1, X, mask1, org_ispin1) == TRUE && state2 == org_ispin3) {
    Jint = conj(tmp_J);
    if (X->Large.mode == M_CORR ||X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC) Jint = 0;
  }
  else {
    return 0.0;
  }

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,       idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
  if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {
      dmv = Jint * v1buf[j];
      tmp_v0[j] += dmv;
      dam_pr += conj(tmp_v1[j]) * dmv;
    }
  }else if(X->Large.mode == H_CORR){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {
      dmv = Jint * v1buf[j];
      dam_pr += conj(tmp_v0[j]) * dmv;
    }
  }else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {
      dmv = Jint * v1buf[j];
      dam_pr += conj(tmp_v1[j]) * dmv;
    }
  }
  return (dam_pr);
#else
 return 0.0;
#endif
}/*double complex child_GC_CisAisCjuAjv_spin_MPIdouble*/
/**
@brief Wrapper for calculating CisAitCjuAju term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_CisAitCjuAju_spin_MPIdouble(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
)
{
#ifdef MPI
  double complex dam_pr;
  dam_pr = child_GC_CisAitCjuAju_spin_MPIdouble(
    X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
    X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4], 
    X->Def.InterAll_OffDiagonal[i_int][5], X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}/*void GC_CisAitCiuAiv_spin_MPIdouble*/
/**
@brief CisAisCjuAjv term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@return @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex child_GC_CisAitCjuAju_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_J,//!<[in] Copupling constatnt
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex *tmp_v1//!<[in] Vector to be producted
) {
#ifdef MPI
  int mask1, mask2, state1, ierr, num1;
  long int origin;
  unsigned long int idim_max_buf, j;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;

  if (org_isite1 == org_isite3 && org_ispin1 == org_ispin3) {//cisaitcisais
    return 0.0;
  }

  mask1 = (int)X->Def.Tpow[org_isite1];
  origin = myrank ^ mask1;
  state1 = (origin & mask1) / mask1;
  mask2 = (int)X->Def.Tpow[org_isite3];
  num1 = child_SpinGC_CisAis(origin + 1, X, mask2, org_ispin3);
  if (state1 == org_ispin2) {
    if (num1 != 0) {
      Jint = tmp_J;
    }
    else {
      return 0.0;
    }
  }/*if (state1 == org_ispin2)*/
  else {//state1 = org_ispin1
    num1 = child_SpinGC_CisAis((unsigned long int) myrank + 1, X, mask2, org_ispin3);
    if (num1 != 0) {
      Jint = conj(tmp_J);
      if (X->Large.mode == M_CORR ||X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC) {
        Jint = 0;
      }
    }
    else {
      return 0.0;
    }
  }

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,       idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {
        dmv = Jint * v1buf[j];
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }else if(X->Large.mode == H_CORR){
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {
        dmv = Jint * v1buf[j];
        dam_pr += conj(tmp_v0[j]) * dmv;
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }else {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {
        dmv = Jint * v1buf[j];
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }
  }/*End of parallel region*/
  return (dam_pr);
#else
 return 0.0;
#endif
}/*double complex child_GC_CisAisCjuAjv_spin_MPIdouble*/
/**
@brief CisAisCjuAjv term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@return @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex child_GC_CisAisCjuAju_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_J,//!<[in] Copupling constatnt
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex *tmp_v1//!<[in] Vector to be producted
){
#ifdef MPI
  long unsigned int mask1, mask2, num1,num2;
  unsigned long int  j;
//  MPI_Status statusMPI;
  double complex dmv, dam_pr;
  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  num1 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, mask1, org_ispin1);
  num2 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, mask2, org_ispin3);
  
  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(tmp_J, X, num1, num2) shared(tmp_v1, tmp_v0)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = num1*num2*tmp_v1[j] * tmp_J;
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++) */
    }else if(X->Large.mode == H_CORR){
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = num1 * num2 * tmp_v1[j] * tmp_J;
        dam_pr += conj(tmp_v0[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = num1 * num2 * tmp_v1[j] * tmp_J;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
  }/*End of parallel region*/
  return(dam_pr);
#else
 return 0.0;
#endif
}/*double complex child_GC_CisAisCjuAju_spin_MPIdouble*/
/**
@brief CisAisCjuAjv term in Spin model + GC
       When both site1 and site2 are in the inter process region.
@return @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex child_GC_CisAisCjuAju_spin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_J,//!<[in] Copupling constatnt
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex *tmp_v1//!<[in] Vector to be producted
) {
#ifdef MPI
  long unsigned int mask1, mask2, num1, num2;
  unsigned long int j;
//  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;
  Jint = tmp_J;
  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  num2 = child_SpinGC_CisAis((unsigned long int) myrank + 1, X, mask2, org_ispin3);

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) private(j, dmv, num1) \
  firstprivate(Jint, X, num2, mask1, org_ispin1) shared(tmp_v1, tmp_v0)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        num1 = child_SpinGC_CisAis(j, X, mask1, org_ispin1);
        dmv = Jint * num1 * num2 * tmp_v1[j];
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }else if(X->Large.mode == H_CORR){
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        num1 = child_SpinGC_CisAis(j, X, mask1, org_ispin1);
        dmv = Jint * num1 * num2 * tmp_v1[j];
        dam_pr += conj(tmp_v0[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        num1 = child_SpinGC_CisAis(j, X, mask1, org_ispin1);
        dmv = Jint * num1 * num2 * tmp_v1[j];
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
  }/*End of parallel region*/
  return (dam_pr);
#else
 return 0.0;
#endif
}/*double complex child_GC_CisAisCjuAju_spin_MPIdouble*/
/**
@brief Exchange and Pairlifting term in Spin model + GC
       When only site2 is in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_CisAitCiuAiv_spin_MPIsingle(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
#ifdef MPI
  double complex dam_pr;  
  dam_pr =child_GC_CisAitCiuAiv_spin_MPIsingle(
    X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
    X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
    X->Def.InterAll_OffDiagonal[i_int][5], X->Def.InterAll_OffDiagonal[i_int][7],
    X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}/*void GC_CisAitCiuAiv_spin_MPIsingle*/

/**
@brief Exchange and Pairlifting term in Spin model + GC
       When only site2 is in the inter process region.
@return @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double complex child_GC_CisAitCiuAiv_spin_MPIsingle(
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
  unsigned long int mask1, idim_max_buf, j, ioff, state1, state1check;
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
    if (X->Large.mode == M_CORR ||X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC) {
      Jint = 0;
    }
  }
  else return 0.0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
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
#pragma omp parallel default(none) reduction(+:dam_pr) private(j, dmv, state1, ioff) \
  firstprivate(idim_max_buf, Jint, X, state1check, mask1) shared(v1buf, tmp_v1, tmp_v0)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 0; j < idim_max_buf; j++) {
        state1 = child_SpinGC_CisAit(j + 1, X, mask1, state1check, &ioff);
        if (state1 != 0) {
          dmv = Jint * v1buf[j + 1];
          tmp_v0[ioff + 1] += dmv;
          dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
        }/*if (state1 != 0)*/
      }/*for (j = 0; j < idim_max_buf; j++)*/
    }else if (X->Large.mode == H_CORR) {
#pragma omp for
      for (j = 0; j < idim_max_buf; j++) {
        state1 = child_SpinGC_CisAit(j + 1, X, mask1, state1check, &ioff);
        if (state1 != 0) {
          dmv = Jint * v1buf[j + 1];
          dam_pr += conj(tmp_v0[ioff + 1]) * dmv;
        }/*if (state1 != 0)*/
      }/*for (j = 0; j < idim_max_buf; j++)*/
    }else {
#pragma omp for
      for (j = 0; j < idim_max_buf; j++) {
        state1 = child_SpinGC_CisAit(j + 1, X, mask1, state1check, &ioff);
        if (state1 != 0) {
          dmv = Jint * v1buf[j + 1];
          dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
        }/*if (state1 != 0)*/
      }/*for (j = 0; j < idim_max_buf; j++)*/
    }
  }/*End of parallel region*/
  return (dam_pr);
#else
 return 0.0;
#endif
}/*void GC_CisAitCiuAiv_spin_MPIsingle*/
/**
@brief Wrapper for CisAisCjuAjv term in Spin model + GC
       When only site2 is in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_CisAisCjuAjv_spin_MPIsingle(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
#ifdef MPI
  double complex dam_pr;  
  dam_pr =child_GC_CisAisCjuAjv_spin_MPIsingle(
    X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
    X->Def.InterAll_OffDiagonal[i_int][4], X->Def.InterAll_OffDiagonal[i_int][5],
    X->Def.InterAll_OffDiagonal[i_int][7], X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}/*void GC_CisAisCjuAjv_spin_MPIsingle*/
/**
@brief CisAisCjuAjv term in Spin model + GC
       When only site2 is in the inter process region.
@return @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex child_GC_CisAisCjuAjv_spin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 1
  int org_ispin3,//!<[in] Spin 2
  int org_ispin4,//!<[in] Spin 2
  double complex tmp_J,//!<[in] Copupling constatnt
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex *tmp_v1//!<[in] Vector to be producted
) {
#ifdef MPI
  int mask2, state2, ierr, origin;
  unsigned long int mask1, idim_max_buf, j, state1, state1check;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[org_isite3];
  origin = myrank ^ mask2;
  state2 = (origin & mask2) / mask2;
  if (state2 == org_ispin4) {
    state1check = (unsigned long int) org_ispin1;
    Jint = tmp_J;
  }
  else if (state2 == org_ispin3) {
    state1check = (unsigned long int) org_ispin1;
    Jint = conj(tmp_J);
    if (X->Large.mode == M_CORR ||X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC) {
      Jint = 0;
    }
  }
  else return 0.0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
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
#pragma omp parallel default(none) reduction(+:dam_pr) private(j, dmv, state1) \
  firstprivate(idim_max_buf, Jint, X, state1check, mask1) shared(v1buf, tmp_v1, tmp_v0)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 0; j < idim_max_buf; j++) {
        state1 = (j & mask1) / mask1;
        if (state1 == state1check) {
          dmv = Jint * v1buf[j + 1];
          tmp_v0[j + 1] += dmv;
          dam_pr += conj(tmp_v1[j + 1]) * dmv;
        }/*if (state1 == state1check)*/
      }/*for (j = 0; j < idim_max_buf; j++)*/
    }else if(X->Large.mode == H_CORR){
#pragma omp for
      for (j = 0; j < idim_max_buf; j++) {
        state1 = (j & mask1) / mask1;
        if (state1 == state1check) {
          dmv = Jint * v1buf[j + 1];
          dam_pr += conj(tmp_v0[j + 1]) * dmv;
        }/*if (state1 == state1check)*/
      }/*for (j = 0; j < idim_max_buf; j++)*/
    }else {
#pragma omp for
      for (j = 0; j < idim_max_buf; j++) {
        state1 = (j & mask1) / mask1;
        if (state1 == state1check) {
          dmv = Jint * v1buf[j + 1];
          dam_pr += conj(tmp_v1[j + 1]) * dmv;
        }/*if (state1 == state1check)*/
      }/*for (j = 0; j < idim_max_buf; j++)*/
    }
  }/*End of parallel region*/
  return (dam_pr);
#else
 return 0.0;
#endif
}/*void GC_CisAitCiuAiv_spin_MPIsingle*/
/**
@brief Wrapper for CisAisCjuAjv term in Spin model + GC
       When only site2 is in the inter process region.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void GC_CisAitCjuAju_spin_MPIsingle(
  unsigned long int i_int,//!<[in] Interaction ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
#ifdef MPI
  double complex dam_pr;  
  dam_pr =child_GC_CisAitCjuAju_spin_MPIsingle(
    X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
    X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4],
    X->Def.InterAll_OffDiagonal[i_int][5], X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}/*void GC_CisAisCjuAjv_spin_MPIsingle*/
/**
@brief CisAisCjuAjv term in Spin model + GC
       When only site2 is in the inter process region.
@return @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex child_GC_CisAitCjuAju_spin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_J,//!<[in] Copupling constatnt
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex *tmp_v1//!<[in] Vector to be producted
) {
#ifdef MPI
  int mask2, state2;
  unsigned long int mask1, j, ioff, state1, state1check;
  //MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[org_isite3];
  state2 = (myrank & mask2) / mask2;

  if (state2 == org_ispin3) {
    state1check = org_ispin2;
    Jint = tmp_J;
  }
  else {
    return 0.0;
  }

  mask1 = (int)X->Def.Tpow[org_isite1];

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) private(j, dmv, state1, ioff) \
  firstprivate(Jint, X, state1check, mask1) shared(tmp_v1, tmp_v0)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 0; j < X->Check.idim_max; j++) {

        state1 = (j & mask1) / mask1;
        ioff = j ^ mask1;
        if (state1 == state1check) {
          dmv = Jint * tmp_v1[j + 1];
        }
        else {
          dmv = conj(Jint) * tmp_v1[j + 1];
        }
        tmp_v0[ioff + 1] += dmv;
        dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
      }/*for (j = 0; j < X->Check.idim_max; j++)*/
    }else if (X->Large.mode == H_CORR) {
#pragma omp for
      for (j = 0; j < X->Check.idim_max; j++) {

        state1 = (j & mask1) / mask1;
        ioff = j ^ mask1;
        if (state1 == state1check) {
          dmv = Jint * tmp_v1[j + 1];
        }
        else {
          dmv = 0.0;
        }
        dam_pr += conj(tmp_v0[ioff + 1]) * dmv;
      }/*for (j = 0; j < X->Check.idim_max; j++)*/
    }


    else if (X->Large.mode == M_CORR) {
#pragma omp for
      for (j = 0; j < X->Check.idim_max; j++) {

        state1 = (j & mask1) / mask1;
        ioff = j ^ mask1;
        if (state1 == state1check) {
          dmv = Jint * tmp_v1[j + 1];
        }
        else {
          dmv = 0.0;
        }
        dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
      }/*for (j = 0; j < X->Check.idim_max; j++)*/
    }
    else {
#pragma omp for
      for (j = 0; j < X->Check.idim_max; j++) {
        state1 = (j & mask1) / mask1;
        ioff = j ^ mask1;
        if (state1 == state1check) {
          dmv = Jint * tmp_v1[j + 1];
        }
        else {
          dmv = conj(Jint) * tmp_v1[j + 1];
        }
        dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
      }/*for (j = 0; j < X->Check.idim_max; j++)*/
    }
  }/*End of parallel region*/
  return (dam_pr);
#else
 return 0.0;
#endif
}/*void GC_CisAitCiuAiv_spin_MPIsingle*/
/**
@brief @f$c_{is}^\dagger c_{is} c_{ju}^\dagger c_{jv}@f$ term in Spin model.
 When both site1 and site3 are in the inter process region.
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
*/
double complex child_GC_CisAisCjuAjv_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_J,//!<[in] Copupling constatnt
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex *tmp_v1//!<[in] Vector to be producted
) {
#ifdef MPI
  unsigned long int off, j;
  int origin, ierr;
  double complex tmp_V, dmv, dam_pr;
  MPI_Status statusMPI;
  if (org_isite1 == org_isite3 && org_ispin1 == org_ispin4) {//cisaisciuais=0 && cisaiucisais=0
    return 0.0;
  }

  if (BitCheckGeneral(myrank, org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow) == TRUE
    && GetOffCompGeneralSpin((unsigned long int) myrank, org_isite3 + 1, org_ispin3, org_ispin4,
      &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE)
    tmp_V = tmp_J;
  else {
    if (GetOffCompGeneralSpin((unsigned long int) myrank, org_isite3 + 1, org_ispin4, org_ispin3,
      &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE)
    {
      if (BitCheckGeneral((unsigned long int)off, org_isite1 + 1, org_ispin1, X->Def.SiteToBit,
        X->Def.Tpow) == TRUE)
      {
        tmp_V = conj(tmp_J);
        if(X->Large.mode == M_CORR || X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC) tmp_V = 0.0;
      }/*BitCheckGeneral(off, org_isite1 + 1, org_ispin1)*/
      else return 0.0;
    }/*GetOffCompGeneralSpin(myrank, org_isite3 + 1, org_ispin4, org_ispin3, &off)*/
    else return 0.0;
  }
  origin = (int)off;
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,  X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) \
private(j, dmv) shared (tmp_v0, tmp_v1, v1buf)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = v1buf[j] * tmp_V;
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }
    }
    else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = v1buf[j] * tmp_V;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }
    }
  }/*End of parallel region*/
  return dam_pr;
#else
 return 0.0;
#endif
}/*double complex child_GC_CisAisCjuAjv_GeneralSpin_MPIdouble*/
/**
@brief @f$c_{is}^\dagger c_{it} c_{ju}^\dagger c_{ju}@f$ term in Spin model.
When both site1 and site3 are in the inter process region.
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
*/
double complex child_GC_CisAitCjuAju_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_J,//!<[in] Copupling constatnt
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] @f${\bf v}_0=H {\bf v}_1@f$
  double complex *tmp_v1//!<[in] Vector to be producted
) {
#ifdef MPI
  unsigned long int j, off;
  int origin, ierr;
  double complex tmp_V, dmv, dam_pr;
  MPI_Status statusMPI;

  if (org_isite1 == org_isite3 && org_ispin1 == org_ispin3) {//cisaitcisais=0 && cisaiscitais=0
    return 0.0;
  }

  if (BitCheckGeneral(myrank, org_isite3 + 1, org_ispin3, X->Def.SiteToBit, X->Def.Tpow) == TRUE
    && GetOffCompGeneralSpin((unsigned long int) myrank, org_isite1 + 1, org_ispin2, org_ispin1, &off,
      X->Def.SiteToBit, X->Def.Tpow) == TRUE)
  {
    tmp_V = conj(tmp_J);
    if (X->Large.mode == M_CORR || X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC) tmp_V = 0.0;
  }
  else if (GetOffCompGeneralSpin((unsigned long int) myrank, org_isite1 + 1, org_ispin1, org_ispin2,
    &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE)
  {
    if (BitCheckGeneral((unsigned long int)off, org_isite3 + 1, org_ispin3,
      X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
      tmp_V = tmp_J;
    }
    else return 0.0;
  }
  else return 0.0;

  origin = (int)off;

  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,  X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) \
shared (tmp_v0, tmp_v1, v1buf)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = v1buf[j] * tmp_V;
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }
    }
    else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = v1buf[j] * tmp_V;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }
    }
  }/*End of parallel region*/
  return dam_pr;
#else
 return 0.0;
#endif
}/*double complex child_GC_CisAitCjuAju_GeneralSpin_MPIdouble*/
/**
@brief Compute @f$c_{is}^\dagger c_{it} c_{ju}^\dagger c_{jv}@f$ term in the
grandcanonical general spin system when both site is in the inter process region
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
*/
double complex child_GC_CisAitCjuAjv_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_J,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1//!<[in] Input wavefunction
) {
#ifdef MPI
  unsigned long int tmp_off, off, j;
  int origin, ierr, ihermite;
  double complex tmp_V, dmv, dam_pr;
  MPI_Status statusMPI;

  ihermite = TRUE;

  if (org_isite1 == org_isite3 && org_ispin1 == org_ispin4 &&
    org_ispin2 == org_ispin3) { //cisaitcitais=cisais && cisaitcitais =cisais
    dam_pr = child_GC_CisAis_GeneralSpin_MPIdouble(org_isite1, org_ispin1, tmp_J, X, tmp_v0, tmp_v1);
    return (dam_pr);
  }
  //cisaitcisait
  if (GetOffCompGeneralSpin((unsigned long int) myrank, org_isite1 + 1, org_ispin1, org_ispin2,
    &tmp_off, X->Def.SiteToBit, X->Def.Tpow) == TRUE) {

    if (GetOffCompGeneralSpin(tmp_off, org_isite3 + 1, org_ispin3, org_ispin4,
      &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE) {

      tmp_V = tmp_J;
    }
    else ihermite = FALSE;
  }
  else {
    ihermite = FALSE;
  }

  if (ihermite == FALSE) {
    if (GetOffCompGeneralSpin((unsigned long int) myrank, org_isite3 + 1, org_ispin4, org_ispin3, &tmp_off,
      X->Def.SiteToBit, X->Def.Tpow) == TRUE) {

      if (GetOffCompGeneralSpin(tmp_off, org_isite1 + 1, org_ispin2, org_ispin1, &off, X->Def.SiteToBit,
                                      X->Def.Tpow) == TRUE) {
        tmp_V = conj(tmp_J);
        if (X->Large.mode == M_CORR || X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC) tmp_V = 0.0;
      }
      else return 0.0;
    }
    else return 0.0;
  }

  origin = (int)off;

  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,  X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) \
  shared (tmp_v0, tmp_v1, v1buf)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = v1buf[j] * tmp_V;
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }
    }
    else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = v1buf[j] * tmp_V;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }
    }
  }/*End of parallel region*/
  return dam_pr;
#else
 return 0.0;
#endif
}/*double complex child_GC_CisAitCjuAjv_GeneralSpin_MPIdouble*/
 /**
 @brief Compute @f$c_{is}^\dagger c_{is} c_{ju}^\dagger c_{ju}@f$ term in the
 grandcanonical general spin system when both site is in the inter process region
 @return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
 */
double complex child_GC_CisAisCjuAju_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_J,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1//!<[in] Input wavefunction
) {
#ifdef MPI
  unsigned long int j, num1;
  double complex tmp_V, dmv, dam_pr;
  //MPI_Status statusMPI;

  num1 = BitCheckGeneral((unsigned long int) myrank, org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);

  if (num1 == TRUE) {
    num1 = BitCheckGeneral((unsigned long int) myrank, org_isite3 + 1, org_ispin3, X->Def.SiteToBit, X->Def.Tpow);
    if (num1 == TRUE) {
      tmp_V = tmp_J;
    }
    else return 0.0;
  }
  else return 0.0;

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) \
shared (tmp_v0, tmp_v1)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = tmp_v1[j] * tmp_V;
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
    else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = tmp_v1[j] * tmp_V;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
  }/*End of parallel region*/
  return dam_pr;
#else
  return 0.0;
#endif
}/*double complex child_GC_CisAisCjuAju_GeneralSpin_MPIdouble*/
 /**
 @brief Compute @f$c_{is}^\dagger c_{it}@f$ term in the
 grandcanonical general spin system when both site is in the inter process region
 @return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
 */
double complex child_GC_CisAit_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  double complex tmp_trans,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1//!<[in] Input wavefunction
) {
#ifdef MPI
  unsigned long int off, j;
  int origin, ierr;
  double complex tmp_V, dmv, dam_pr;
  MPI_Status statusMPI;

  if (GetOffCompGeneralSpin((unsigned long int) myrank, org_isite1 + 1, org_ispin1, org_ispin2,
    &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
    tmp_V = tmp_trans;
  }
  else if (GetOffCompGeneralSpin((unsigned long int) myrank,
    org_isite1 + 1, org_ispin2, org_ispin1, &off,
    X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
    tmp_V = conj(tmp_trans);
    if (X->Large.mode == M_CORR || X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC) tmp_V = 0.0;
  }
  else return 0.0;

  origin = (int)off;

  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,  X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) \
shared (tmp_v0, tmp_v1, v1buf)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = v1buf[j] * tmp_V;
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
    else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = v1buf[j] * tmp_V;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
  }/*End of parallel region*/
  return dam_pr;
#else
  return 0.0;
#endif
}/*double complex child_GC_CisAit_GeneralSpin_MPIdouble*/
 /**
 @brief Compute @f$c_{is}^\dagger c_{is}@f$ term in the
 grandcanonical general spin system when both site is in the inter process region
 @return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
 */
double complex child_GC_CisAis_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  double complex tmp_trans,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1//!<[in] Input wavefunction
) {
#ifdef MPI
  unsigned long int j, num1;
  double complex tmp_V, dmv, dam_pr;
  //MPI_Status statusMPI;

  num1 = BitCheckGeneral((unsigned long int) myrank,
    org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);
  if (num1 != 0) {
    tmp_V = tmp_trans;
  }
  else return 0.0;

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) \
shared (tmp_v0, tmp_v1)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = tmp_v1[j] * tmp_V;
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
    else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = tmp_v1[j] * tmp_V;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
  }/*End of parallel region*/
  return dam_pr;
#else
  return 0.0;
#endif
}/*double complex child_GC_CisAis_GeneralSpin_MPIdouble*/
 /**
 @brief Compute @f$c_{is} c_{is}^\dagger@f$ term in the
 grandcanonical general spin system when both site is in the inter process region
 @return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
 */
double complex child_GC_AisCis_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  double complex tmp_trans,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1//!<[in] Input wavefunction
) {
#ifdef MPI
  unsigned long int j, num1;
  double complex tmp_V, dmv, dam_pr;
  //MPI_Status statusMPI;

  num1 = BitCheckGeneral((unsigned long int) myrank,
    org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);
  if (num1 == 0) {
    tmp_V = tmp_trans;
  }
  else return 0.0;

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) \
shared (tmp_v0, tmp_v1)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = tmp_v1[j] * tmp_V;
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
    else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = tmp_v1[j] * tmp_V;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
  }/*End of Parallel region*/
  return dam_pr;
#else
  return 0.0;
#endif
}/*double complex child_GC_AisCis_GeneralSpin_MPIdouble*/
/**
@brief Compute @f$c_{is}^\dagger c_{it}@f$ term in the
canonical general spin system when both site is in the inter process region
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
*/
double complex child_CisAit_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  double complex tmp_trans,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Input wavefunction
  double complex *tmp_v1buf,//!<[inout] buffer for wavefunction
  unsigned long int idim_max,//!<[in] Similar to CheckList::idim_max
  long unsigned int *list_1_org,//!<[in] Similar to ::list_1
  long unsigned int *list_1buf_org,//!<[in] Similar to ::list_1buf
  long unsigned int _ihfbit//!<[in] Similer to LargeList::ihfbit
)
{
#ifdef MPI
  unsigned long int off, j, tmp_off,idim_max_buf;
  int origin, ierr;
  double complex tmp_V, dmv;
  MPI_Status statusMPI;
  
  if (GetOffCompGeneralSpin((unsigned long int) myrank, org_isite1 + 1, org_ispin1, org_ispin2,
                            &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
    tmp_V = tmp_trans;
  }
  else if (GetOffCompGeneralSpin((unsigned long int) myrank,
                                 org_isite1 + 1, org_ispin2, org_ispin1, &off,
                                 X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
    tmp_V = conj(tmp_trans);
    if (X->Large.mode == M_CORR || X->Large.mode == H_CORR || X->Large.mode ==M_CALCSPEC) tmp_V = 0.0;
  }
  else return 0.0;
  
  origin = (int) off;

  ierr = MPI_Sendrecv(&idim_max,     1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, 
                      MPI_COMM_WORLD, &statusMPI);
  if(ierr != 0) exitMPI(-1);
  
  ierr = MPI_Sendrecv(list_1_org,        idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                      list_1buf_org, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  
  ierr = MPI_Sendrecv(tmp_v1,    idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none)\
firstprivate(X, tmp_V, idim_max_buf, list_1buf_org) private(j, dmv, tmp_off) \
shared (tmp_v0, tmp_v1, v1buf)
    for (j = 1; j <= idim_max_buf; j++) {
      ConvertToList1GeneralSpin(list_1buf_org[j], X->Large.ihfbit, &tmp_off);
      dmv = v1buf[j] * tmp_V;
      tmp_v0[tmp_off] += dmv;
    }/*for (j = 1; j <= idim_max_buf; j++)*/
  }
  else {
    tmp_off = 0;
    return 0;
  }
  return 1;
#else
 return 0.0;
#endif
}/*double complex child_CisAit_GeneralSpin_MPIdouble*/

/**
@brief Compute @f$c_{is}^\dagger c_{is}c_{ju}^\dagger c_{jv}@f$ term in the
grandcanonical general spin system when one of these site is in the inter process region
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
*/
double complex child_GC_CisAisCjuAjv_GeneralSpin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_J,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1//!<[in] Input wavefunction
){
#ifdef MPI
  unsigned long int off, j, num1;
  int origin, ierr, isite, IniSpin;
  double complex tmp_V, dmv, dam_pr;
  MPI_Status statusMPI;

  if (GetOffCompGeneralSpin((unsigned long int)myrank,
    org_isite3 + 1, org_ispin3, org_ispin4, &off,
    X->Def.SiteToBit, X->Def.Tpow) == TRUE)
  {
    tmp_V = tmp_J;
    isite = org_isite1 + 1;
    IniSpin = org_ispin1;
  }
  else if (GetOffCompGeneralSpin((unsigned long int)myrank,
    org_isite3 + 1, org_ispin4, org_ispin3, &off,
    X->Def.SiteToBit, X->Def.Tpow) == TRUE)
  {
    tmp_V = conj(tmp_J);
    if (X->Large.mode == M_CORR || X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC) tmp_V = 0.0;
    isite = org_isite1 + 1;
    IniSpin = org_ispin1;
  }
  else return 0.0;
  
  origin = (int)off;
  
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,  X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, isite, IniSpin) \
private(j, dmv, num1) shared (tmp_v0, tmp_v1, v1buf)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        num1 = BitCheckGeneral(j - 1, isite, IniSpin, X->Def.SiteToBit, X->Def.Tpow);
        if (num1 != 0) {
          dmv = v1buf[j] * tmp_V;
          tmp_v0[j] += dmv;
          dam_pr += conj(tmp_v1[j]) * dmv;
        }/*if (num1 != 0)*/
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
    else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        num1 = BitCheckGeneral(j - 1, isite, IniSpin, X->Def.SiteToBit, X->Def.Tpow);
        if (num1 != 0) {
          dmv = v1buf[j] * tmp_V;
          dam_pr += conj(tmp_v1[j]) * dmv;
        }/*if (num1 != 0)*/
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
  }/*End of parallel region*/
  return dam_pr;
#else
 return 0.0;
#endif
}/*double complex child_GC_CisAisCjuAjv_GeneralSpin_MPIsingle*/
/**
@brief Compute @f$c_{is}^\dagger c_{it}c_{ju}^\dagger c_{ju}@f$ term in the
grandcanonical general spin system when one of these site is in the inter process region
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
*/
double complex child_GC_CisAitCjuAju_GeneralSpin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_J,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1//!<[in] Input wavefunction
){
#ifdef MPI
  unsigned long int num1, j, off;
  int isite, IniSpin, FinSpin;
  double complex tmp_V, dmv, dam_pr;
  //MPI_Status statusMPI;

  num1 = BitCheckGeneral((unsigned long int)myrank, 
    org_isite3+1, org_ispin3, X->Def.SiteToBit, X->Def.Tpow);
  if(num1 != 0){
    tmp_V = tmp_J;
    isite = org_isite1 + 1;
    IniSpin = org_ispin2;
    FinSpin = org_ispin1;
  }
  else return 0.0;

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) \
firstprivate(X, tmp_V, isite, IniSpin, FinSpin) private(j, dmv, num1, off) \
shared (tmp_v0, tmp_v1, v1buf)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        if (GetOffCompGeneralSpin(j - 1, isite, IniSpin, FinSpin, &off,
          X->Def.SiteToBit, X->Def.Tpow) == TRUE)
        {
          dmv = tmp_v1[j] * tmp_V;
          tmp_v0[off + 1] += dmv;
          dam_pr += conj(tmp_v1[off + 1]) * dmv;
        }
        else if (GetOffCompGeneralSpin(j - 1, isite, FinSpin, IniSpin, &off,
          X->Def.SiteToBit, X->Def.Tpow) == TRUE)
        {
          dmv = tmp_v1[j] * conj(tmp_V);
          tmp_v0[off + 1] += dmv;
          dam_pr += conj(tmp_v1[off + 1]) * dmv;
        }
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
    else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        if (GetOffCompGeneralSpin(j - 1, isite, IniSpin, FinSpin, &off,
          X->Def.SiteToBit, X->Def.Tpow) == TRUE) 
        {
          dmv = tmp_v1[j] * tmp_V;
          dam_pr += conj(tmp_v1[off + 1]) * dmv;
        }
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
  }/*End of parallel region*/
  return dam_pr;
#else
 return 0.0;
#endif
}/*double complex child_GC_CisAitCjuAju_GeneralSpin_MPIsingle*/
/**
@brief Compute @f$c_{is}^\dagger c_{is}c_{ju}^\dagger c_{jv}@f$ term in the
grandcanonical general spin system when one of these site is in the inter process region
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
*/
double complex child_GC_CisAitCjuAjv_GeneralSpin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_J,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1//!<[in] Input wavefunction
){
#ifdef MPI
  unsigned long int off, j;
  int origin, ierr, isite, IniSpin, FinSpin;
  double complex tmp_V, dmv, dam_pr;
  MPI_Status statusMPI;

  if (GetOffCompGeneralSpin((unsigned long int)myrank,
    org_isite3 + 1, org_ispin3, org_ispin4, &off,
    X->Def.SiteToBit, X->Def.Tpow) == TRUE)
  {
    tmp_V = tmp_J;
    isite = org_isite1 + 1;
    IniSpin = org_ispin2;
    FinSpin = org_ispin1;
  }
  else if (GetOffCompGeneralSpin((unsigned long int)myrank,
    org_isite3 + 1, org_ispin4, org_ispin3, &off,
    X->Def.SiteToBit, X->Def.Tpow) == TRUE)
  {
    tmp_V = conj(tmp_J);
    if (X->Large.mode == M_CORR || X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC) tmp_V = 0.0;
    isite = org_isite1 + 1;
    IniSpin = org_ispin1;
    FinSpin = org_ispin2;
  }
  else return 0.0;

  origin = (int)off;

  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,  X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) \
firstprivate(X, tmp_V, isite, IniSpin, FinSpin) private(j, dmv, off) shared (tmp_v0, tmp_v1, v1buf)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        if (GetOffCompGeneralSpin(j - 1, isite, IniSpin, FinSpin, &off,
          X->Def.SiteToBit, X->Def.Tpow) == TRUE)
        {
          dmv = v1buf[j] * tmp_V;
          tmp_v0[off + 1] += dmv;
          dam_pr += conj(tmp_v1[off + 1]) * dmv;
        }
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
    else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        if (GetOffCompGeneralSpin(j - 1, isite, IniSpin, FinSpin, &off,
          X->Def.SiteToBit, X->Def.Tpow) == TRUE) 
        {
          dmv = v1buf[j] * tmp_V;
          dam_pr += conj(tmp_v1[off + 1]) * dmv;
        }
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
  }/*End of parallel region*/
  return dam_pr;
#else
 return 0.0;
#endif
}/*double complex child_GC_CisAitCjuAjv_GeneralSpin_MPIsingle*/
/**
@brief Compute @f$c_{is}^\dagger c_{is}c_{ju}^\dagger c_{ju}@f$ term in the
grandcanonical general spin system when one of these site is in the inter process region
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
*/
double complex child_GC_CisAisCjuAju_GeneralSpin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_J,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1//!<[in] Input wavefunction
){
#ifdef MPI
  unsigned long int j, num1;
  double complex tmp_V, dmv, dam_pr;
  //MPI_Status statusMPI;

  num1 = BitCheckGeneral((unsigned long int)myrank, org_isite3+1, org_ispin3, X->Def.SiteToBit, X->Def.Tpow);
  if (num1 != FALSE) {
    tmp_V = tmp_J;
  }
  else return 0.0;
  
  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) \
firstprivate(X, tmp_V, org_isite1, org_ispin1) private(j, dmv, num1) shared (tmp_v0, tmp_v1)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        num1 = BitCheckGeneral(j - 1, org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);

        dmv = tmp_v1[j] * tmp_V * num1;
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
    else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        num1 = BitCheckGeneral(j - 1, org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);
        dmv = tmp_v1[j] * tmp_V * num1;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
  }/*End of parallel region*/
  return dam_pr;
#else
 return 0.0;
#endif
}/*double complex child_GC_CisAisCjuAju_GeneralSpin_MPIsingle*/
/**
@brief Compute @f$c_{is}^\dagger c_{it}c_{ju}^\dagger c_{jv}@f$ term in the
canonical general spin system when both sites are in the inter process region
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
*/
double complex child_CisAitCjuAjv_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_J,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1//!<[in] Input wavefunction
){
#ifdef MPI
  unsigned long int tmp_off, off, j, idim_max_buf;
  int origin, ierr;
  double complex tmp_V, dmv, dam_pr;
  MPI_Status statusMPI;
  int ihermite=TRUE;

  if (GetOffCompGeneralSpin((unsigned long int)myrank, org_isite1 + 1, org_ispin1, org_ispin2, &tmp_off, X->Def.SiteToBit, X->Def.Tpow) == TRUE)
  {
    if (GetOffCompGeneralSpin(tmp_off, org_isite3 + 1, org_ispin3, org_ispin4, &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE)
    {
      tmp_V = tmp_J;
    }
    else{
      ihermite =FALSE;
    }
  }
  else{
    ihermite=FALSE;
  }
  
  if(ihermite==FALSE){
    if(GetOffCompGeneralSpin((unsigned long int)myrank, org_isite3 + 1, org_ispin4, org_ispin3, &tmp_off, X->Def.SiteToBit, X->Def.Tpow) == TRUE)
      {
 if (GetOffCompGeneralSpin(tmp_off, org_isite1 + 1, org_ispin2, org_ispin1, &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE)
   {
     tmp_V = conj(tmp_J);
     if(X->Large.mode == M_CORR || X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC){
       tmp_V=0.0;
     }
   }
 else return 0.0;
      }
    else return 0.0;
  }
  
  
  origin = (int)off;

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
#pragma omp parallel default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, idim_max_buf) \
private(j, dmv, off) shared (tmp_v0, tmp_v1, list_1buf, v1buf)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {
        ConvertToList1GeneralSpin(list_1buf[j], X->Check.sdim, &off);
        dmv = v1buf[j] * tmp_V;
        tmp_v0[off] += dmv;
        dam_pr += conj(tmp_v1[off]) * dmv;
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }
    else {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {
        ConvertToList1GeneralSpin(list_1buf[j], X->Check.sdim, &off);
        dmv = v1buf[j] * tmp_V;
        dam_pr += conj(tmp_v1[off]) * dmv;
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }
  }/*End of parallel region*/
  return dam_pr;
#else
 return 0.0;
#endif
}/*double complex child_CisAitCjuAjv_GeneralSpin_MPIdouble*/
/**
@brief Compute @f$c_{is}^\dagger c_{is}c_{ju}^\dagger c_{ju}@f$ term in the
canonical general spin system when both sites are in the inter process region
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
*/
double complex child_CisAisCjuAju_GeneralSpin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_J,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1//!<[in] Input wavefunction
) {
#ifdef MPI
  unsigned long int j, num1;
  double complex tmp_V, dmv, dam_pr;

  if (org_isite1 == org_isite3 && org_ispin1 == org_ispin3) {
    num1 = BitCheckGeneral((unsigned long int) myrank, org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);
    if (num1 != FALSE) {
      tmp_V = tmp_J;
    }
    else {
      return 0.0;
    }
  }
  else {
    num1 = BitCheckGeneral((unsigned long int) myrank, org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);
    if (num1 != FALSE) {
      num1 = BitCheckGeneral((unsigned long int) myrank, org_isite3 + 1, org_ispin3, X->Def.SiteToBit,
        X->Def.Tpow);
      if (num1 != FALSE) {
        tmp_V = tmp_J;
      }
      else {
        return 0.0;
      }
    }
    else {
      return 0.0;
    }
  }
  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) \
shared (tmp_v0, tmp_v1)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = tmp_v1[j] * tmp_V;
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
    else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = tmp_v1[j] * tmp_V;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
  }/*End of parallel region*/
  return dam_pr;
#else
 return 0.0;
#endif
}/*double complex child_CisAisCjuAju_GeneralSpin_MPIdouble*/
/**
@brief Compute @f$c_{is}^\dagger c_{is}c_{ju}^\dagger c_{ju}@f$ term in the
canonical general spin system when one of these sites is in the inter process region
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
*/
double complex child_CisAisCjuAju_GeneralSpin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_J,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1//!<[in] Input wavefunction
)
{
#ifdef MPI
  unsigned long int j, num1;
  double complex tmp_V, dmv, dam_pr;
  //MPI_Status statusMPI;

  num1 = BitCheckGeneral((unsigned long int) myrank, org_isite3 + 1, org_ispin3, X->Def.SiteToBit, X->Def.Tpow);
  if (num1 != FALSE) {
    tmp_V = tmp_J;
  }
  else return 0.0;

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) \
firstprivate(X, tmp_V, org_isite1, org_ispin1) private(j, dmv, num1) shared (tmp_v0, tmp_v1, list_1)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        num1 = BitCheckGeneral(list_1[j], org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);

        dmv = tmp_v1[j] * tmp_V * num1;
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
    else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        num1 = BitCheckGeneral(list_1[j], org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);

        dmv = tmp_v1[j] * tmp_V * num1;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
  }/*End of parallel region*/
  return dam_pr;
#else
 return 0.0;
#endif
}/*double complex child_CisAisCjuAju_GeneralSpin_MPIsingle*/
 /**
 @brief Compute @f$c_{is}^\dagger c_{it}c_{ju}^\dagger c_{jv}@f$ term in the
 canonical general spin system when one of these sites is in the inter process region
 @return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
 */
double complex child_CisAitCjuAjv_GeneralSpin_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_J,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1//!<[in] Input wavefunction
){
#ifdef MPI
  unsigned long int tmp_off, off, j, idim_max_buf;
  int origin, ierr, isite, IniSpin, FinSpin;
  double complex tmp_V, dmv, dam_pr;
  MPI_Status statusMPI;
  
  if (GetOffCompGeneralSpin((unsigned long int)myrank,
    org_isite3 + 1, org_ispin3, org_ispin4, &off,
    X->Def.SiteToBit, X->Def.Tpow) == TRUE)
  {
    tmp_V = tmp_J;
    isite = org_isite1 + 1;
    IniSpin = org_ispin2;
    FinSpin = org_ispin1;
  }
  else if (GetOffCompGeneralSpin((unsigned long int)myrank,
    org_isite3 + 1, org_ispin4, org_ispin3, &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE)
  {
    tmp_V = conj(tmp_J);
    if (X->Large.mode == M_CORR || X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC) tmp_V = 0.0;
    isite = org_isite1 + 1;
    IniSpin = org_ispin1;
    FinSpin = org_ispin2;
  }
  else return 0.0;

  origin = (int)off;
  
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
#pragma omp parallel default(none) reduction(+:dam_pr) \
firstprivate(X, tmp_V, idim_max_buf, IniSpin, FinSpin, isite) \
private(j, dmv, off, tmp_off) shared (tmp_v0, tmp_v1, list_1buf, v1buf)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {

        if (GetOffCompGeneralSpin(list_1buf[j], isite, IniSpin, FinSpin, &tmp_off,
          X->Def.SiteToBit, X->Def.Tpow) == TRUE) 
        {
          ConvertToList1GeneralSpin(tmp_off, X->Check.sdim, &off);
          dmv = v1buf[j] * tmp_V;
          tmp_v0[off] += dmv;
          dam_pr += conj(tmp_v1[off]) * dmv;
        }
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }
    else {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {

        if (GetOffCompGeneralSpin(list_1buf[j], isite, IniSpin, FinSpin, &tmp_off,
          X->Def.SiteToBit, X->Def.Tpow) == TRUE) 
        {
          ConvertToList1GeneralSpin(tmp_off, X->Check.sdim, &off);
          dmv = v1buf[j] * tmp_V;
          dam_pr += conj(tmp_v1[off]) * dmv;
        }
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }
  }/*End of parallel region*/
  return dam_pr;
#else
  return 0.0;
#endif
}/*double complex child_CisAitCjuAjv_GeneralSpin_MPIsingle*/
/**
@brief Hopping term in Spin + GC
       When both site1 and site2 are in the inter process region.
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex child_GC_CisAit_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_ispin2,//!<[in] Spin 2
  double complex tmp_trans,//!<[in] Coupling constant
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  int mask1, state1, ierr, origin;
  unsigned long int idim_max_buf, j;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;
  
  mask1 = (int)X->Def.Tpow[org_isite1];
  origin = myrank ^ mask1;
  state1 = (origin & mask1)/mask1;

  //fprintf(stdout, "Debug: myrank=%d, origin=%d, state1=%d\n", myrank, origin, state1);

  if(state1 ==  org_ispin2){
    trans = tmp_trans;
  }
  else if(state1 == org_ispin1) {
    trans = conj(tmp_trans);
    if(X->Large.mode == M_CORR || X->Large.mode == H_CORR|| X->Large.mode ==M_CALCSPEC || X->Large.mode == M_MLTPLY2){
      trans = 0.0;
    }
  }
  else{
    return 0.0;
  }

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,       idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) private(j, dmv) \
firstprivate(idim_max_buf, trans, X) shared(v1buf, tmp_v1, tmp_v0)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC || X->Large.mode == M_MLTPLY2) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = trans * v1buf[j];
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }else if (X->Large.mode == H_CORR) {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = trans * v1buf[j];
        dam_pr += conj(tmp_v0[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }else {
#pragma omp for
      for (j = 1; j <= X->Check.idim_max; j++) {
        dmv = trans * v1buf[j];
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= X->Check.idim_max; j++)*/
    }
  }/*End of parallel region*/
  return (dam_pr);
#else
 return 0.0;
#endif
}/*double complex  child_GC_CisAit_spin_MPIdouble*/
/**
@brief Hopping term in Spin + Canonical for CalcSpectrum
       When both site1 and site2 are in the inter process region.
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex child_CisAit_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin2,//!<[in] Spin 2
  double complex tmp_trans,//!<[in] Coupling constant
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1, /**< [in] v0 = H v1*/
  double complex *tmp_v1buf,//!<[in] buffer for wavefunction
  unsigned long int idim_max,//!<[in] Similar to CheckList::idim_max
  long unsigned int *Tpow,//!<[in] Similar to DefineList::Tpow
  long unsigned int *list_1_org,//!<[in] Similar to ::list_1
  long unsigned int *list_1buf_org,//!<[in] Similar to ::list_1buf
  long unsigned int *list_2_1_target,//!<[in] Similar to ::list_2_1
  long unsigned int *list_2_2_target,//!<[in] Similar to ::list_2_2
  long unsigned int _irght,//!<[in] Similer to LargeList::irght
  long unsigned int _ilft,//!<[in] Similer to LargeList::ilft
  long unsigned int _ihfbit//!<[in] Similer to LargeList::ihfbit
){
#ifdef MPI
  int mask1, state1, ierr, origin;
  unsigned long int idim_max_buf, j;
  unsigned long int tmp_off;
  MPI_Status statusMPI;
  double complex trans, dmv;
  
  mask1 = (int)X->Def.Tpow[org_isite1];
  origin = myrank ^ mask1;
  state1 = (origin & mask1)/mask1;

  if(state1 ==  org_ispin2){
    trans = tmp_trans;
  }
  else{
    trans =0.0;
  }

  //  fprintf(stdout, "Debug: myrank=%d, origin=%d, trans=%lf\n", myrank, origin, trans);
  
  ierr = MPI_Sendrecv(&idim_max,     1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  ierr = MPI_Sendrecv(list_1_org,        idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                      list_1buf_org, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  ierr = MPI_Sendrecv(tmp_v1,    idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
    
  if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) private(j, dmv, tmp_off) \
firstprivate(idim_max_buf, trans, X, list_1buf_org, list_2_1_target, list_2_2_target) \
shared(v1buf, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {
      GetOffComp(list_2_1_target, list_2_2_target, list_1buf_org[j], X->Large.irght, X->Large.ilft, X->Large.ihfbit, &tmp_off);
      dmv = trans * v1buf[j];
      tmp_v0[tmp_off] += dmv;
    }
  }
  else {
    tmp_off = 0;
    return 0;
  }
  return 1;
#else
 return 0.0;
#endif
}/*double complex  child_CisAit_spin_MPIdouble*/
/**
@brief Hopping term in Spin + GC
       When both site1 and site2 are in the inter process region.
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex child_GC_CisAis_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  double complex tmp_trans,//!<[in] Coupling constant
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
 double complex *tmp_v1 /**< [in] v0 = H v1*/
){
#ifdef MPI
  long unsigned int j;
  int mask1;
  int ibit1;
  double complex dam_pr;
  mask1 = (int)X->Def.Tpow[org_isite1];
  ibit1 = (((unsigned long int)myrank& mask1)/mask1)^(1-org_ispin1);

  dam_pr = 0.0;
#pragma omp parallel reduction(+:dam_pr)default(none) shared(tmp_v1, tmp_v0, ibit1) \
  firstprivate(X, tmp_trans) private(j)
  {
    if (ibit1 != 0) {
      if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
#pragma omp for
        for (j = 1; j <= X->Check.idim_max; j++) {
          tmp_v0[j] += tmp_v1[j] * tmp_trans;
          dam_pr += tmp_trans * conj(tmp_v1[j]) * tmp_v1[j];
        }/*for (j = 1; j <= X->Check.idim_max; j++)*/
      }else if (X->Large.mode == H_CORR){
#pragma omp for
        for (j = 1; j <= X->Check.idim_max; j++) {
          dam_pr += tmp_trans * conj(tmp_v0[j]) * tmp_v1[j];
        }/*for (j = 1; j <= X->Check.idim_max; j++)*/
      }else {
#pragma omp for
        for (j = 1; j <= X->Check.idim_max; j++) {
          dam_pr += tmp_trans * conj(tmp_v1[j]) * tmp_v1[j];
        }/*for (j = 1; j <= X->Check.idim_max; j++)*/
      }
    }/*if (ibit1 != 0)*/
  }/*End of parallel region*/
  return dam_pr;
#else
 return 0.0;
#endif
}/*double complex child_GC_CisAis_spin_MPIdouble*/
/**
@brief Hopping term in Spin + GC
       When both site1 and site2 are in the inter process region.
@return @f$\langle v_1| H_{\rm this} | v_1 \rangle@f$
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex child_GC_AisCis_spin_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  double complex tmp_trans,//!<[in] Coupling constant
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/
){
#ifdef MPI
  long unsigned int j;
  int mask1;
  int ibit1;
  double complex dam_pr;
  mask1 = (int)X->Def.Tpow[org_isite1];
  ibit1 = (((unsigned long int)myrank& mask1) / mask1) ^ (1 - org_ispin1);

  dam_pr = 0.0;
#pragma omp parallel reduction(+:dam_pr)default(none) shared(tmp_v1, tmp_v0, ibit1) \
  firstprivate(X, tmp_trans) private(j)
  {
    if (ibit1 == 0) {
      if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
#pragma omp for
        for (j = 1; j <= X->Check.idim_max; j++) {
          tmp_v0[j] += tmp_v1[j] * tmp_trans;
          dam_pr += tmp_trans * conj(tmp_v1[j]) * tmp_v1[j];
        }/*for (j = 1; j <= X->Check.idim_max; j++)*/
      } else if (X->Large.mode == H_CORR) { //
#pragma omp for
        for (j = 1; j <= X->Check.idim_max; j++) {
          dam_pr    += tmp_trans * conj(tmp_v0[j]) * tmp_v1[j];
        }/*for (j = 1; j <= X->Check.idim_max; j++)*/
      }else {
#pragma omp for
        for (j = 1; j <= X->Check.idim_max; j++) {
          dam_pr += tmp_trans * conj(tmp_v1[j]) * tmp_v1[j];
        }/*for (j = 1; j <= X->Check.idim_max; j++)*/
      }
    }/*if (ibit1 == 0)*/
  }/*End of parallel region*/
  return dam_pr;
#else
 return 0.0;
#endif
}/*double complex child_GC_AisCis_spin_MPIdouble*/
