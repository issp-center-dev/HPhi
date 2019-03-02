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
@brief Functions for Hubbar + MPI (Core)
*/
#include "Common.h"
#include "mltplyCommon.h"
#include "mltplyHubbardCore.h"
#include "mltplyMPIHubbard.h"
#include "mltplyMPIHubbardCore.h"
#include "bitcalc.h"
#include "wrapperMPI.h"
/**
@brief Check whether this site is in the inter process region or not
@return 1 if it is inter-process region, 0 if not.
*/
int CheckPE(
  int org_isite,//!<[in] Site index
  struct BindStruct *X//!<[inout]
){
  if (org_isite + 1 > X->Def.Nsite) {
    return TRUE;
  }
  else {
    return FALSE;
  }
}/*int CheckPE*/
/**
@brief Check the occupation of @f$(i,s)@f$ state,
and compute the index of final wavefunction associated to 
@f$c^\dagger_{is}@f$
@return 1 if unoccupied, 0 if occupied
*/
int CheckBit_Cis(
  long unsigned int is1_spin,//!<[in] Index of site+spin
  long unsigned int orgbit,//!<[in] Index of initial wavefunction
  long unsigned int *offbit//!<[out] Index of final wavefunction
) {
  long unsigned int ibit_tmp;
  ibit_tmp = orgbit & is1_spin;
  if (ibit_tmp == 0) {
    *offbit = orgbit + is1_spin;
    return TRUE;
  }
  *offbit = 0;
  return FALSE;
}/*int CheckBit_Cis*/
/**
@brief Check the occupation of @f$(i,s)@f$ state,
and compute the index of final wavefunction associated to
@f$c_{jt}@f$
@return 1 if occupied, 0 if unoccupied
*/
int CheckBit_Ajt(
  long unsigned int is1_spin,//!<[in] Index of site+spin
  long unsigned int orgbit,//!<[in] Index of initial wavefunction
  long unsigned int *offbit//!<[out] Index of final wavefunction
) {
  long unsigned int ibit_tmp;
  ibit_tmp = orgbit & is1_spin;
  if (ibit_tmp != 0) {
    *offbit = orgbit - is1_spin;
    return TRUE;
  }
  *offbit = 0;
  return FALSE;
}/*int CheckBit_Ajt*/
/**
@brief Compute the index of final wavefunction associated to
@f$c_{4}^\dagger c_{3}c_{2}^\dagger c_{1}@f$, and
check whether this operator is relevant or not
@return 1 if relevant, 0 if irrelevant
*/
int CheckBit_InterAllPE(
  int org_isite1,//!<[in] Site 1
  int org_isigma1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_isigma2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_isigma3,//!<[in] Spin 3
  int org_isite4,//!<[in] Site 4
  int org_isigma4,//!<[in] Spin 4
  struct BindStruct *X,//!<[inout]
  long unsigned int orgbit,//!<[in] Index of initial wavefunction
  long unsigned int *offbit//!<[out] Index of final wavefunction
){
  long unsigned int tmp_ispin;
  long unsigned int tmp_org, tmp_off;
  int iflgBitExist = TRUE;
  tmp_org=orgbit;
  tmp_off=0;

  if (CheckPE(org_isite1, X) == TRUE) {
    tmp_ispin = X->Def.Tpow[2 * org_isite1 + org_isigma1];
    if (CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
    tmp_org = tmp_off;
  }

  if (CheckPE(org_isite2, X) == TRUE ) {
    tmp_ispin = X->Def.Tpow[2 * org_isite2 + org_isigma2];
    if (CheckBit_Cis(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
    tmp_org = tmp_off;
  }

  if (CheckPE(org_isite3, X) == TRUE) {
    tmp_ispin = X->Def.Tpow[2 * org_isite3 + org_isigma3];
    if (CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
    tmp_org = tmp_off;
  }

  if (CheckPE(org_isite4, X) == TRUE) {
    tmp_ispin = X->Def.Tpow[2 * org_isite4 + org_isigma4];
    if (CheckBit_Cis(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
    tmp_org = tmp_off;
  }

  if(iflgBitExist != TRUE){
    *offbit=0;
    return FALSE;
  }
  
  *offbit=tmp_org;
  return TRUE;
}/*int CheckBit_InterAllPE*/
/**
@brief Check the occupation of both site 1 and site 3
@return 1 if both sites are occupied, 0 if not
*/
int CheckBit_PairPE(
  int org_isite1,//!<[in] Site 1
  int org_isigma1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_isigma3,//!<[in] Spin 4
  struct BindStruct *X,//!<[inout]
  long unsigned int orgbit//!<[in] Index pf intial wavefunction
){
  long unsigned int tmp_ispin;
  long unsigned int tmp_org, tmp_off;
  int iflgBitExist = TRUE;
  tmp_org=orgbit;
  
  if(CheckPE(org_isite1, X)==TRUE){
    tmp_ispin = X->Def.Tpow[2 * org_isite1 + org_isigma1];
    if (CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist=FALSE;
    }
  }
  
  if (CheckPE(org_isite3, X) == TRUE) {
    tmp_ispin = X->Def.Tpow[2 * org_isite3 + org_isigma3];
    if (CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE) {
      iflgBitExist = FALSE;
    }
  }
  
  if(iflgBitExist != TRUE){
    return FALSE;
  }

  return TRUE;
}/*int CheckBit_PairPE*/
/**
@brief Compute the index of final wavefunction associated to
@f$c_{4}^\dagger c_{3}c_{2}^\dagger c_{1}@f$, and
Fermion sign
@return 1 if relevant, 0 if irrelevant
*/
int GetSgnInterAll(
  unsigned long int isite1,//!<[in] Site 1
  unsigned long int isite2,//!<[in] Site 2
  unsigned long int isite3,//!<[in] Site 3
  unsigned long int isite4,//!<[in] Site 4
  int *Fsgn,//!<[out] Fermion sign
  struct BindStruct *X,//!<[inout]
  unsigned long int orgbit,//!<[in] Index of the initial state
  unsigned long int *offbit//!<[out] Index of the final state
){
  long unsigned int diffA;
  long unsigned int tmp_off;
  long unsigned int tmp_ispin1, tmp_ispin2;
  int tmp_sgn=0;

  tmp_ispin1=isite2;
  tmp_ispin2=isite1;

  if (tmp_ispin1 == tmp_ispin2) {
    if ((orgbit & tmp_ispin1) == 0) {
      *offbit = 0;
      *Fsgn = tmp_sgn;
      return FALSE;
    }
    tmp_sgn=1;
    tmp_off = orgbit;
  }
  else {
    if (tmp_ispin2 > tmp_ispin1) diffA = tmp_ispin2 - tmp_ispin1 * 2;
    else diffA = tmp_ispin1-tmp_ispin2*2;

    tmp_sgn=X_GC_CisAjt(orgbit, X, tmp_ispin1, tmp_ispin2, tmp_ispin1+tmp_ispin2, diffA, &tmp_off);
    if(tmp_sgn ==0){
      *offbit =0;
      *Fsgn = 0;
      return FALSE;
    }
  }

  tmp_ispin1 = isite4;
  tmp_ispin2 = isite3;

  if(tmp_ispin1 == tmp_ispin2){
    if( (tmp_off & tmp_ispin1) == 0){
      *offbit =0;
      *Fsgn = 0;
      return FALSE;
    }
    *offbit=tmp_off;
  }
  else{
    if(tmp_ispin2 > tmp_ispin1) diffA = tmp_ispin2 - tmp_ispin1*2;
    else diffA = tmp_ispin1-tmp_ispin2*2;
    tmp_sgn *=X_GC_CisAjt(tmp_off, X, tmp_ispin1, tmp_ispin2, tmp_ispin1+tmp_ispin2, diffA, offbit);
    
    if(tmp_sgn ==0){
      *offbit =0;
      *Fsgn = 0;
      return FALSE;
    }
  }
  
  *Fsgn =tmp_sgn;
  *offbit = *offbit%X->Def.OrgTpow[2*X->Def.Nsite];
  return TRUE;
}/*int GetSgnInterAll*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{jt}@f$
term of grandcanonical Hubbard system
*/
void X_GC_child_CisAisCjtAjt_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_V,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] Resulting wavefunction
  double complex **tmp_v1//!<[inout] Initial wavefunction
) {
  int iCheck;
  unsigned long int tmp_ispin1;
  unsigned long int i_max = X->Check.idim_max;
  unsigned long int tmp_off, j;
  int one = 1;

  iCheck=CheckBit_PairPE(org_isite1, org_ispin1, org_isite3, org_ispin3, X, (long unsigned int) myrank);
  if(iCheck != TRUE){
    return;
  }

  if (org_isite1 + 1 > X->Def.Nsite && org_isite3 + 1 > X->Def.Nsite) {
    zaxpy_long(i_max*nstate, tmp_V, &tmp_v1[1][0], &tmp_v0[1][0]);
  }/*if (org_isite1 + 1 > X->Def.Nsite && org_isite3 + 1 > X->Def.Nsite)*/
  else if (org_isite1 + 1 > X->Def.Nsite || org_isite3 + 1 > X->Def.Nsite) {
    if (org_isite1 > org_isite3) tmp_ispin1 = X->Def.Tpow[2 * org_isite3 + org_ispin3];
    else                         tmp_ispin1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];

#pragma omp parallel  default(none) shared(org_isite1, org_ispin1, org_isite3, org_ispin3, nstate, tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X) private(j, tmp_off, tmp_ispin1)
#pragma omp for
    for (j = 1; j <= i_max; j++) {
      if (CheckBit_Ajt(tmp_ispin1, j - 1, &tmp_off) == TRUE) {
        zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
      }
    }/*for (j = 1; j <= i_max; j++)*/
  }
}/*double complex X_GC_child_CisAisCjtAjt_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{ku}@f$
term of grandcanonical Hubbard system
*/
void X_GC_child_CisAjtCkuAku_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_V,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] Resulting wavefunction
  double complex **tmp_v1//!<[inout] Initial wavefunction
) {
  unsigned long int i_max = X->Check.idim_max;
  unsigned long int idim_max_buf;
  int iCheck, ierr, Fsgn;
  unsigned long int isite1, isite2, isite3;
  unsigned long int tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4;
  unsigned long int j, Asum, Adiff;
  double complex dmv;
  unsigned long int origin, tmp_off;
  unsigned long int org_rankbit;
  int one = 1;

  iCheck = CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite3, org_ispin3, X, (long unsigned int) myrank, &origin);
  isite1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
  isite2 = X->Def.Tpow[2 * org_isite2 + org_ispin2];
  isite3 = X->Def.Tpow[2 * org_isite3 + org_ispin3];

  if (iCheck == TRUE) {
    tmp_isite1 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
    tmp_isite2 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
    tmp_isite3 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
    tmp_isite4 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
    Asum = tmp_isite1 + tmp_isite2;
    if (tmp_isite2 > tmp_isite1) Adiff = tmp_isite2 - tmp_isite1 * 2;
    else Adiff = tmp_isite1 - tmp_isite2 * 2;
  }
  else {
    iCheck = CheckBit_InterAllPE(org_isite3, org_ispin3, org_isite3, org_ispin3, org_isite2, org_ispin2, org_isite1, org_ispin1, X, (long unsigned int) myrank, &origin);
    if (iCheck == TRUE) {
      tmp_V = conj(tmp_V);
      tmp_isite4 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
      tmp_isite3 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
      tmp_isite2 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
      tmp_isite1 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
      Asum = tmp_isite3 + tmp_isite4;
      if (tmp_isite4 > tmp_isite3) Adiff = tmp_isite4 - tmp_isite3 * 2;
      else Adiff = tmp_isite3 - tmp_isite4 * 2;
      if (X->Large.mode == M_CORR || X->Large.mode == M_CALCSPEC) {
        tmp_V = 0;
      }
    }
    else {
      return;
    }
  }

  if (myrank == origin) {// only k is in PE

    if (CheckBit_Ajt(isite3, myrank, &tmp_off) == FALSE) return;

#pragma omp parallel default(none)  \
firstprivate(i_max,X,Asum,Adiff,isite1,isite2, tmp_V) private(j,tmp_off) shared(tmp_v0, tmp_v1)
    {
#pragma omp for
      for (j = 1; j <= i_max; j++) 
        GC_CisAjt(j, nstate, tmp_v0, tmp_v1, X, isite2, isite1, Asum, Adiff, tmp_V, &tmp_off);

      if (X->Large.mode != M_CORR) {
#pragma omp for
        for (j = 1; j <= i_max; j++) 
          GC_CisAjt(j, nstate, tmp_v0, tmp_v1, X, isite1, isite2, Asum, Adiff, tmp_V, &tmp_off);
      }/*if (X->Large.mode != M_CORR)*/
    }/*End of paralle region*/
    return;
  }//myrank =origin
  else {
    idim_max_buf = SendRecv_i(origin, X->Check.idim_max);
    SendRecv_cv(origin, X->Check.idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);

#pragma omp parallel default(none)  private(j, dmv, tmp_off, Fsgn, org_rankbit, Adiff) \
shared(v1buf, tmp_v1, nstate, tmp_v0, myrank, origin, isite3, org_isite3, isite1, isite2, org_isite2, org_isite1) \
firstprivate(idim_max_buf, tmp_V, X, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4)
    {
      if (org_isite1 + 1 > X->Def.Nsite && org_isite2 + 1 > X->Def.Nsite) {
        if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
        else Adiff = isite1 - isite2 * 2;
        SgnBit(((long unsigned int) myrank & Adiff), &Fsgn);
        tmp_V *= Fsgn;

        if (org_isite3 + 1 > X->Def.Nsite) {
#pragma omp for
          for (j = 1; j <= idim_max_buf; j++) {
            zaxpy_(&nstate, &tmp_V, &v1buf[j][0], &one, &tmp_v0[j][0], &one);
          }/*for (j = 1; j <= idim_max_buf; j++)*/
        }
        else { //org_isite3 <= X->Def.Nsite
#pragma omp for
          for (j = 1; j <= idim_max_buf; j++) {
            if (CheckBit_Ajt(isite3, j - 1, &tmp_off) == TRUE) {
              zaxpy_(&nstate, &tmp_V, &v1buf[j][0], &one, &tmp_v0[j][0], &one);
            }
          }/*for (j = 1; j <= idim_max_buf; j++)*/
        }
      }/*if (org_isite1 + 1 > X->Def.Nsite && org_isite2 + 1 > X->Def.Nsite)*/
      else {
        org_rankbit = X->Def.OrgTpow[2 * X->Def.Nsite] * origin;
#pragma omp for
        for (j = 1; j <= idim_max_buf; j++) {
          if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, X, (j - 1) + org_rankbit, &tmp_off) == TRUE) {
            dmv = tmp_V * Fsgn;
            zaxpy_(&nstate, &dmv, &v1buf[j][0], &one, &tmp_v0[tmp_off + 1][0], &one);
          }
        }/*for (j = 1; j <= idim_max_buf; j++)*/
      }
    }/*End of parallel region*/
  }/*myrank != origin*/
}/*double complex X_GC_child_CisAjtCkuAku_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{ku}@f$
term of grandcanonical Hubbard system
*/
void X_GC_child_CisAisCjtAku_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_isite4,//!<[in] Site 4
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_V,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] Resulting wavefunction
  double complex **tmp_v1//!<[inout] Initial wavefunction
) {
  X_GC_child_CisAjtCkuAku_Hubbard_MPI(
    org_isite4, org_ispin4, org_isite3, org_ispin3,
    org_isite1, org_ispin1, conj(tmp_V), X, nstate, tmp_v0, tmp_v1);
}/*double complex X_GC_child_CisAisCjtAku_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{lv}@f$
term of grandcanonical Hubbard system
*/
void X_GC_child_CisAjtCkuAlv_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_isite4,//!<[in] Site 4
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_V,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] Resulting wavefunction
  double complex **tmp_v1//!<[inout] Initial wavefunction
) {
  unsigned long int i_max = X->Check.idim_max;
  unsigned long int idim_max_buf;
  int iCheck, ierr, Fsgn;
  unsigned long int isite1, isite2, isite3, isite4;
  unsigned long int tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4;
  unsigned long int j, Adiff, Bdiff;
  double complex dmv;
  unsigned long int origin, tmp_off, tmp_off2;
  unsigned long int org_rankbit;
  int iFlgHermite = FALSE;
  int one = 1;

  iCheck = CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2,
                               org_isite3, org_ispin3, org_isite4, org_ispin4,
                               X, (long unsigned int) myrank, &origin);
  isite1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
  isite2 = X->Def.Tpow[2 * org_isite2 + org_ispin2];
  isite3 = X->Def.Tpow[2 * org_isite3 + org_ispin3];
  isite4 = X->Def.Tpow[2 * org_isite4 + org_ispin4];

  if (iCheck == TRUE) {
    tmp_isite1 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
    tmp_isite2 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
    tmp_isite3 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
    tmp_isite4 = X->Def.OrgTpow[2 * org_isite4 + org_ispin4];
  }
  else {
    iCheck = CheckBit_InterAllPE(org_isite4, org_ispin4, org_isite3, org_ispin3,
                                 org_isite2, org_ispin2, org_isite1, org_ispin1,
                                 X, (long unsigned int) myrank, &origin);
    if (iCheck == TRUE) {
      tmp_V = conj(tmp_V);
      tmp_isite4 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
      tmp_isite3 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
      tmp_isite2 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
      tmp_isite1 = X->Def.OrgTpow[2 * org_isite4 + org_ispin4];
      iFlgHermite = TRUE;
      if (X->Large.mode == M_CORR || X->Large.mode == M_CALCSPEC) {
        tmp_V = 0;
      }
    }
    else {
      return;
    }
  }

  if (myrank == origin) {
    if (isite1 == isite4 && isite2 == isite3) { // CisAjvCjvAis =Cis(1-njv)Ais=nis-nisnjv
            //calc nis
      X_GC_child_CisAis_Hubbard_MPI(org_isite1, org_ispin1, tmp_V, X, nstate, tmp_v0, tmp_v1);
      //calc -nisniv
      X_GC_child_CisAisCjtAjt_Hubbard_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3, -tmp_V, X, nstate, tmp_v0, tmp_v1);
    }/*if (isite1 == isite4 && isite2 == isite3)*/
    else if (isite2 == isite3) { // CisAjvCjvAku= Cis(1-njv)Aku=-CisAkunjv+CisAku: j is in PE
            //calc CisAku
      if (isite4 > isite1) Adiff = isite4 - isite1 * 2;
      else Adiff = isite1 - isite4 * 2;

#pragma omp parallel for default(none)  private(j, tmp_off) \
firstprivate(i_max, tmp_V, X, isite1, isite4, Adiff) shared(tmp_v1, tmp_v0)
      for (j = 1; j <= i_max; j++) 
        GC_CisAjt(j - 1, nstate, tmp_v0, tmp_v1, X, isite1, isite4, (isite1 + isite4), Adiff, tmp_V, &tmp_off);
      
      //calc -CisAku njv
      X_GC_child_CisAjtCkuAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite4, org_ispin4, 
                                          org_isite2, org_ispin2, -tmp_V, X, nstate, tmp_v0, tmp_v1);
      if (X->Large.mode != M_CORR) { //for hermite
#pragma omp parallel for default(none)  private(j, tmp_off) \
firstprivate(i_max, tmp_V, X, isite1, isite4, Adiff) shared(tmp_v1, tmp_v0)
        for (j = 1; j <= i_max; j++) 
          GC_CisAjt(j - 1, nstate, tmp_v0, tmp_v1, X, isite4, isite1, (isite1 + isite4), Adiff, tmp_V, &tmp_off);
        
        //calc -njvCkuAis
        X_GC_child_CisAisCjtAku_Hubbard_MPI(org_isite2, org_ispin2, org_isite4, org_ispin4,
                                            org_isite1, org_ispin1, -tmp_V, X, nstate, tmp_v0, tmp_v1);
      }/*if (X->Large.mode != M_CORR)*/
    }/*if (isite2 == isite3)*/
    else {// CisAjtCkuAis = -CisAisCkuAjt: i is in PE
      X_GC_child_CisAisCjtAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3,
                                          org_isite2, org_ispin2, -tmp_V, X, nstate, tmp_v0, tmp_v1);
      if (X->Large.mode != M_CORR) { //for hermite
        X_GC_child_CisAisCjtAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite2, org_ispin2,
                                            org_isite3, org_ispin3, -tmp_V, X, nstate, tmp_v0, tmp_v1);
      }/*if (X->Large.mode != M_CORR)*/
    }/*if (isite2 != isite3)*/
    return;
  }//myrank =origin
  else {
    idim_max_buf = SendRecv_i(origin, X->Check.idim_max);
    SendRecv_cv(origin, X->Check.idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);

    if (org_isite1 + 1 > X->Def.Nsite && org_isite2 + 1 > X->Def.Nsite
     && org_isite3 + 1 > X->Def.Nsite && org_isite4 + 1 > X->Def.Nsite) {

      if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
      else Adiff = isite1 - isite2 * 2;
      if (isite4 > isite3) Bdiff = isite4 - isite3 * 2;
      else Bdiff = isite3 - isite4 * 2;

      if (iFlgHermite == FALSE) {
        Fsgn = X_GC_CisAjt((long unsigned int) myrank, X, isite2, isite1, (isite1 + isite2), Adiff, &tmp_off2);
        Fsgn *= X_GC_CisAjt(tmp_off2, X, isite4, isite3, (isite3 + isite4), Bdiff, &tmp_off);
        tmp_V *= Fsgn;
      }/*if (iFlgHermite == FALSE)*/
      else {
        Fsgn = X_GC_CisAjt((long unsigned int) myrank, X, isite3, isite4, (isite3 + isite4), Bdiff, &tmp_off2);
        Fsgn *= X_GC_CisAjt(tmp_off2, X, isite1, isite2, (isite1 + isite2), Adiff, &tmp_off);
        tmp_V *= Fsgn;
      }/*if (iFlgHermite == TRUE)*/

      zaxpy_long(i_max*nstate, tmp_V, &v1buf[1][0], &tmp_v0[1][0]);
    }
    else {
      org_rankbit = X->Def.OrgTpow[2 * X->Def.Nsite] * origin;
#pragma omp parallel for default(none)  private(j, dmv, tmp_off, Fsgn) firstprivate(idim_max_buf, tmp_V, X, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit) shared(v1buf, tmp_v1, tmp_v0)
      for (j = 1; j <= idim_max_buf; j++) {
        if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, X, (j - 1) + org_rankbit, &tmp_off) == TRUE) {
          dmv = tmp_V * Fsgn;
          zaxpy_(&nstate, &dmv, &v1buf[j][0], &one, &tmp_v0[tmp_off + 1][0], &one);
        }
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }
  }/*myrank != origin*/
}/*double complex X_GC_child_CisAjtCkuAlv_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{is}@f$
term of grandcanonical Hubbard system
*/
void X_GC_child_CisAis_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  double complex tmp_V,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] Resulting wavefunction
  double complex **tmp_v1//!<[inout] Initial wavefunction
) {
  unsigned long int i_max = X->Check.idim_max;
  unsigned long int j, isite1, tmp_off;
  int one = 1;

  isite1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
  if (org_isite1 + 1 > X->Def.Nsite) {
    if (CheckBit_Ajt(isite1, (unsigned long int) myrank, &tmp_off) == FALSE) return;

    zaxpy_long(i_max*nstate, tmp_V, &tmp_v1[1][0], &tmp_v0[1][0]);
  }/*if (org_isite1 + 1 > X->Def.Nsite)*/
  else {
#pragma omp parallel  default(none) shared(tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X, isite1) private(j, tmp_off)
    {
#pragma omp for
      for (j = 1; j <= i_max; j++) {
        if (CheckBit_Ajt(isite1, j - 1, &tmp_off) == TRUE) {
          zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
        }/*if (CheckBit_Ajt(isite1, j - 1, &tmp_off) == TRUE)*/
      }/*for (j = 1; j <= i_max; j++)*/
    }/*End of parallel region*/
  }/*if (org_isite1 + 1 <= X->Def.Nsite)*/
}/*double complex X_GC_child_CisAis_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt}@f$
term of grandcanonical Hubbard system
*/
void X_GC_child_CisAjt_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  double complex tmp_trans,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] Resulting wavefunction
  double complex **tmp_v1//!<[inout] Initial wavefunction
) {

  if (org_isite1 + 1 > X->Def.Nsite && org_isite2 + 1 > X->Def.Nsite) {
    X_GC_child_general_hopp_MPIdouble(org_isite1, org_ispin1, org_isite2, org_ispin2, tmp_trans, X, nstate, tmp_v0, tmp_v1);
  }
  else if (org_isite1 + 1 > X->Def.Nsite || org_isite2 + 1 > X->Def.Nsite) {
    X_GC_child_general_hopp_MPIsingle(org_isite1, org_ispin1, org_isite2, org_ispin2, tmp_trans, X, nstate, tmp_v0, tmp_v1);
  }
  else {
    //error message will be added.
    exitMPI(-1);
  }
}/*double complex X_GC_child_CisAjt_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{jt}@f$
term of canonical Hubbard system
*/
void X_child_CisAisCjtAjt_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_V,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] Resulting wavefunction
  double complex **tmp_v1//!<[inout] Initial wavefunction
) {
  int iCheck;
  unsigned long int tmp_ispin1;
  unsigned long int i_max = X->Check.idim_max;
  unsigned long int tmp_off, j;
  int one = 1;

  iCheck = CheckBit_PairPE(org_isite1, org_ispin1, org_isite3, org_ispin3, X, (long unsigned int) myrank);
  if (iCheck != TRUE) return;
  
  if (org_isite1 + 1 > X->Def.Nsite && org_isite3 + 1 > X->Def.Nsite) {
#pragma omp for
    zaxpy_long(i_max*nstate, tmp_V, &tmp_v1[1][0], &tmp_v0[1][0]);
  }/*if (org_isite1 + 1 > X->Def.Nsite && org_isite3 + 1 > X->Def.Nsite)*/
  else if (org_isite1 + 1 > X->Def.Nsite || org_isite3 + 1 > X->Def.Nsite) {
    if (org_isite1 > org_isite3) tmp_ispin1 = X->Def.Tpow[2 * org_isite3 + org_ispin3];
    else                         tmp_ispin1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];

#pragma omp parallel for default(none) \
shared(tmp_v0, tmp_v1, list_1, org_isite1, org_ispin1, org_isite3, org_ispin3) \
  firstprivate(i_max, tmp_V, X, tmp_ispin1) private(j, tmp_off)
    for (j = 1; j <= i_max; j++) {
      if (CheckBit_Ajt(tmp_ispin1, list_1[j], &tmp_off) == TRUE) {
        zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
      }
    }/*for (j = 1; j <= i_max; j++)*/
  }/*if (org_isite1 + 1 > X->Def.Nsite || org_isite3 + 1 > X->Def.Nsite)*/
}/*double complex X_child_CisAisCjtAjt_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{lv}@f$
term of canonical Hubbard system
*/
void X_child_CisAjtCkuAlv_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_isite4,//!<[in] Site 4
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_V,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] Resulting wavefunction
  double complex **tmp_v1//!<[inout] Initial wavefunction
) {
  unsigned long int i_max = X->Check.idim_max;
  unsigned long int idim_max_buf;
  int iCheck, ierr, Fsgn;
  unsigned long int isite1, isite2, isite3, isite4;
  unsigned long int tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4;
  unsigned long int j, Adiff, Bdiff;
  double complex dmv;
  unsigned long int origin, tmp_off, tmp_off2;
  unsigned long int org_rankbit, ioff;
  int iFlgHermite = FALSE;
  int one = 1;

  iCheck = CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2,
                               org_isite3, org_ispin3, org_isite4, org_ispin4,
                               X, (long unsigned int) myrank, &origin);
  //printf("iCheck=%d, myrank=%d, origin=%d\n", iCheck, myrank, origin);
  isite1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
  isite2 = X->Def.Tpow[2 * org_isite2 + org_ispin2];
  isite3 = X->Def.Tpow[2 * org_isite3 + org_ispin3];
  isite4 = X->Def.Tpow[2 * org_isite4 + org_ispin4];

  if (iCheck == TRUE) {
    tmp_isite1 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
    tmp_isite2 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
    tmp_isite3 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
    tmp_isite4 = X->Def.OrgTpow[2 * org_isite4 + org_ispin4];
  }/*if (iCheck == TRUE)*/
  else {
    iCheck = CheckBit_InterAllPE(org_isite4, org_ispin4, org_isite3, org_ispin3,
                                 org_isite2, org_ispin2, org_isite1, org_ispin1,
                                 X, (long unsigned int) myrank, &origin);
    if (iCheck == TRUE) {
      tmp_V = conj(tmp_V);
      tmp_isite4 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
      tmp_isite3 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
      tmp_isite2 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
      tmp_isite1 = X->Def.OrgTpow[2 * org_isite4 + org_ispin4];
      iFlgHermite = TRUE;
      if (X->Large.mode == M_CORR || X->Large.mode == M_CALCSPEC) tmp_V = 0;     
    }/*if (iCheck == TRUE)*/
    else return;
  }/*if (iCheck == FALSE)*/

  if (myrank == origin) {
    if (isite1 == isite4 && isite2 == isite3) { // CisAjvCjvAis =Cis(1-njv)Ais=nis-nisnjv
            //calc nis
      X_child_CisAis_Hubbard_MPI(org_isite1, org_ispin1, tmp_V, X, nstate, tmp_v0, tmp_v1);
      //calc -nisniv
      X_child_CisAisCjtAjt_Hubbard_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3, -tmp_V, X, nstate, tmp_v0, tmp_v1);
    }/* if (isite1 == isite4 && isite2 == isite3)*/
    else if (isite2 == isite3) { // CisAjvCjvAku= Cis(1-njv)Aku=-CisAkunjv+CisAku: j is in PE
      if (isite4 > isite1) Adiff = isite4 - isite1 * 2;
      else Adiff = isite1 - isite4 * 2;

      //calc CisAku
#pragma omp parallel for default(none)  private(j, tmp_off) \
firstprivate(i_max, tmp_V, X, isite1, isite4, Adiff) shared(tmp_v1, nstate, tmp_v0, list_1)
      for (j = 1; j <= i_max; j++)
        CisAjt(j, nstate, tmp_v0, tmp_v1, X, isite1, isite4, (isite1 + isite4), Adiff, tmp_V);
      
      //calc -CisAku njv
      X_child_CisAjtCkuAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite4, org_ispin4,
                                       org_isite2, org_ispin2, -tmp_V, X, nstate, tmp_v0, tmp_v1);

      if (X->Large.mode != M_CORR) {  //for hermite
#pragma omp parallel for default(none)  private(j, tmp_off) \
firstprivate(i_max, tmp_V, X, isite1, isite4, Adiff) shared(tmp_v1, tmp_v0)
        for (j = 1; j <= i_max; j++) 
          CisAjt(j, nstate, tmp_v0, tmp_v1, X, isite4, isite1, (isite1 + isite4), Adiff, tmp_V);
                
        //calc -njvCkuAis
        X_child_CisAisCjtAku_Hubbard_MPI(org_isite2, org_ispin2, org_isite4, org_ispin4, 
                                         org_isite1, org_ispin1, -tmp_V, X, nstate, tmp_v0, tmp_v1);
      }/*if (X->Large.mode != M_CORR)*/
    }/*if (isite2 == isite3)*/
    else {// CisAjtCkuAis = -CisAisCkuAjt: i is in PE
      X_child_CisAisCjtAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3, 
                                       org_isite2, org_ispin2, -tmp_V, X, nstate, tmp_v0, tmp_v1);

      if (X->Large.mode != M_CORR) //for hermite: CisAkuCjtAis=-CisAisCjtAku
        X_child_CisAisCjtAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite2, org_ispin2,
                                         org_isite3, org_ispin3, -tmp_V, X, nstate, tmp_v0, tmp_v1);    
    }/*if (isite2 != isite3)*/
    return;
  }//myrank =origin
  else {
    idim_max_buf = SendRecv_i(origin, X->Check.idim_max);
    SendRecv_iv(origin, X->Check.idim_max + 1, idim_max_buf + 1, list_1, list_1buf);

    SendRecv_cv(origin, X->Check.idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);
    if (org_isite1 + 1 > X->Def.Nsite && org_isite2 + 1 > X->Def.Nsite
     && org_isite3 + 1 > X->Def.Nsite && org_isite4 + 1 > X->Def.Nsite)
    {
      if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
      else Adiff = isite1 - isite2 * 2;
      if (isite4 > isite3) Bdiff = isite4 - isite3 * 2;
      else Bdiff = isite3 - isite4 * 2;

      if (iFlgHermite == FALSE) {
        Fsgn = X_GC_CisAjt((long unsigned int) myrank, X, isite2, isite1, (isite1 + isite2), Adiff, &tmp_off2);
        Fsgn *= X_GC_CisAjt(tmp_off2, X, isite4, isite3, (isite3 + isite4), Bdiff, &tmp_off);
        tmp_V *= Fsgn;
      }/*if (iFlgHermite == FALSE)*/
      else {
        Fsgn = X_GC_CisAjt((long unsigned int) myrank, X, isite3, isite4, (isite3 + isite4), Bdiff, &tmp_off2);
        Fsgn *= X_GC_CisAjt(tmp_off2, X, isite1, isite2, (isite1 + isite2), Adiff, &tmp_off);
        tmp_V *= Fsgn;
      }/*if (iFlgHermite == TRUE)*/
#pragma omp parallel default(none)  private(j, ioff) \
firstprivate(idim_max_buf, tmp_V, X) shared(v1buf, tmp_v1, nstate, tmp_v0, list_2_1, list_2_2, list_1buf)
      {
#pragma omp for
        for (j = 1; j <= idim_max_buf; j++) {
          if (GetOffComp(list_2_1, list_2_2, list_1buf[j],
            X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == TRUE)
          {
            zaxpy_(&nstate, &tmp_V, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
          }
        }/*for (j = 1; j <= idim_max_buf; j++)*/
      }/*End of parallel region*/
    }//org_isite1+1 > X->Def.Nsite && org_isite2+1 > X->Def.Nsite
            // && org_isite3+1 > X->Def.Nsite && org_isite4+1 > X->Def.Nsite
    else {
      org_rankbit = X->Def.OrgTpow[2 * X->Def.Nsite] * origin;

#pragma omp parallel default(none)  private(j, dmv, tmp_off, Fsgn, ioff) \
firstprivate(myrank, idim_max_buf, tmp_V, X, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit, \
org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite4, org_ispin4) \
shared(v1buf, tmp_v1, nstate, tmp_v0, list_1buf, list_2_1, list_2_2)
      {
#pragma omp for
        for (j = 1; j <= idim_max_buf; j++) {
          if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, X,
            list_1buf[j] + org_rankbit, &tmp_off) == TRUE)
          {
            if (GetOffComp(list_2_1, list_2_2, tmp_off, X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == TRUE)
            {
              dmv = tmp_V * Fsgn;
              zaxpy_(&nstate, &dmv, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
            }
          }
        }/*for (j = 1; j <= idim_max_buf; j++)*/
      }/*End of parallel region*/
    }
  }/*if (myrank != origin)*/
}/*double complex X_child_CisAjtCkuAlv_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{ku}@f$
term of canonical Hubbard system
*/
void X_child_CisAjtCkuAku_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  double complex tmp_V,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] Resulting wavefunction
  double complex **tmp_v1//!<[inout] Initial wavefunction
) {
  unsigned long int i_max = X->Check.idim_max;
  unsigned long int idim_max_buf, ioff;
  int iCheck, ierr, Fsgn;
  unsigned long int isite1, isite2, isite3;
  unsigned long int tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4;
  unsigned long int j, Asum, Adiff;
  double complex dmv;
  unsigned long int origin, tmp_off;
  unsigned long int org_rankbit;
  int one = 1;
  //printf("Deubg0-0: org_isite1=%d, org_ispin1=%d, org_isite2=%d, org_ispin2=%d, org_isite3=%d, org_ispin3=%d\n", org_isite1, org_ispin1,org_isite2, org_ispin2,org_isite3, org_ispin3);
  iCheck = CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite3, org_ispin3, X, (long unsigned int) myrank, &origin);
  //printf("iCheck=%d, myrank=%d, origin=%d\n", iCheck, myrank, origin);

  isite1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
  isite2 = X->Def.Tpow[2 * org_isite2 + org_ispin2];
  isite3 = X->Def.Tpow[2 * org_isite3 + org_ispin3];

  if (iCheck == TRUE) {
    tmp_isite1 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
    tmp_isite2 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
    tmp_isite3 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
    tmp_isite4 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
    Asum = tmp_isite1 + tmp_isite2;
    if (tmp_isite2 > tmp_isite1) Adiff = tmp_isite2 - tmp_isite1 * 2;
    else Adiff = tmp_isite1 - tmp_isite2 * 2;
  }/*if (iCheck == TRUE)*/
  else {
    iCheck = CheckBit_InterAllPE(org_isite3, org_ispin3, org_isite3, org_ispin3, org_isite2, org_ispin2, org_isite1, org_ispin1, X, (long unsigned int) myrank, &origin);
    if (iCheck == TRUE) {
      tmp_V = conj(tmp_V);
      tmp_isite4 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
      tmp_isite3 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
      tmp_isite2 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
      tmp_isite1 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
      Asum = tmp_isite3 + tmp_isite4;
      if (tmp_isite4 > tmp_isite3) Adiff = tmp_isite4 - tmp_isite3 * 2;
      else Adiff = tmp_isite3 - tmp_isite4 * 2;
      if (X->Large.mode == M_CORR || X->Large.mode == M_CALCSPEC) tmp_V = 0;
      //printf("tmp_isite1=%ld, tmp_isite2=%ld, Adiff=%ld\n", tmp_isite1, tmp_isite2, Adiff);
    }/*if (iCheck == TRUE)*/
    else return;   
  }/*if (iCheck == FALSE)*/

  if (myrank == origin) {// only k is in PE
    //for hermite
#pragma omp parallel default(none)  \
firstprivate(i_max, Asum, Adiff, isite1, isite2, tmp_V, X) private(j) shared(tmp_v0, tmp_v1)
    {
#pragma omp for
      for (j = 1; j <= i_max; j++)
        CisAjt(j, nstate, tmp_v0, tmp_v1, X, isite1, isite2, Asum, Adiff, tmp_V);

      if (X->Large.mode != M_CORR) {
#pragma omp for
        for (j = 1; j <= i_max; j++)
          CisAjt(j, nstate, tmp_v0, tmp_v1, X, isite2, isite1, Asum, Adiff, tmp_V);
      }/*if (X->Large.mode != M_CORR)*/
    }/*End of parallel region*/
    return;
  }//myrank =origin
  else {
    idim_max_buf = SendRecv_i(origin, X->Check.idim_max);
    SendRecv_iv(origin, X->Check.idim_max + 1, idim_max_buf + 1, list_1, list_1buf);
    SendRecv_cv(origin, X->Check.idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &v1buf[1][0]);

#pragma omp parallel default(none)  private(j, dmv, ioff, tmp_off, Fsgn, Adiff) \
firstprivate(idim_max_buf, tmp_V, X, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit, isite3) \
shared(v1buf, tmp_v1, nstate, tmp_v0, list_1buf, list_2_1, list_2_2, origin, org_isite3, myrank, isite1, isite2, org_isite1, org_isite2)
    {

      if (org_isite1 + 1 > X->Def.Nsite && org_isite2 + 1 > X->Def.Nsite) {
        if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
        else Adiff = isite1 - isite2 * 2;
        SgnBit(((long unsigned int) myrank & Adiff), &Fsgn);
        tmp_V *= Fsgn;

        if (org_isite3 + 1 > X->Def.Nsite) {
#pragma omp for
          for (j = 1; j <= idim_max_buf; j++) {
            GetOffComp(list_2_1, list_2_2, list_1buf[j],
              X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
            zaxpy_(&nstate, &tmp_V, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
          }/*for (j = 1; j <= idim_max_buf; j++)*/
        }/*if (org_isite3 + 1 > X->Def.Nsite)*/
        else { //org_isite3 <= X->Def.Nsite
#pragma omp for
          for (j = 1; j <= idim_max_buf; j++) {
            if (CheckBit_Ajt(isite3, list_1buf[j], &tmp_off) == TRUE) {
              GetOffComp(list_2_1, list_2_2, list_1buf[j],
                X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
              zaxpy_(&nstate, &tmp_V, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
            }
          }/*for (j = 1; j <= idim_max_buf; j++)*/
        }/*if (org_isite3 + 1 <= X->Def.Nsite)*/
      }/*if (org_isite1 + 1 > X->Def.Nsite && org_isite2 + 1 > X->Def.Nsite)*/
      else {
        org_rankbit = X->Def.OrgTpow[2 * X->Def.Nsite] * origin;
#pragma omp for
        for (j = 1; j <= idim_max_buf; j++) {
          if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, X,
            list_1buf[j] + org_rankbit, &tmp_off) == TRUE) {
            dmv = tmp_V * Fsgn;
            GetOffComp(list_2_1, list_2_2, tmp_off,
              X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
            zaxpy_(&nstate, &dmv, &v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
          }
        }/*for (j = 1; j <= idim_max_buf; j++)*/
      }
    }/*End of parallel region*/
  }/*if (myrank != origin)*/
}/*double complex X_child_CisAjtCkuAku_Hubbard_MPI*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{ku}@f$
term of canonical Hubbard system
*/
void X_child_CisAisCjtAku_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite3,//!<[in] Site 3
  int org_ispin3,//!<[in] Spin 3
  int org_isite4,//!<[in] Site 4
  int org_ispin4,//!<[in] Spin 4
  double complex tmp_V,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] Resulting wavefunction
  double complex **tmp_v1//!<[inout] Initial wavefunction
) {
  X_child_CisAjtCkuAku_Hubbard_MPI(
    org_isite4, org_ispin4, org_isite3, org_ispin3,
    org_isite1, org_ispin1, conj(tmp_V), X, nstate, tmp_v0, tmp_v1);
}/*double complex X_child_CisAisCjtAku_Hubbard_MPI*/

void X_child_CisAis_Hubbard_MPI(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  double complex tmp_V,//!<[in] Coupling constant
  struct BindStruct *X,//!<[inout]
  int nstate, double complex **tmp_v0,//!<[inout] Resulting wavefunction
  double complex **tmp_v1//!<[inout] Initial wavefunction
) {
  unsigned long int i_max = X->Check.idim_max;
  unsigned long int j, isite1, tmp_off;
  int one = 1;

  isite1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
  if (org_isite1 + 1 > X->Def.Nsite) {
    if (CheckBit_Ajt(isite1, (unsigned long int) myrank, &tmp_off) == FALSE)
      return;

    zaxpy_long(i_max*nstate, tmp_V, &tmp_v1[1][0], &tmp_v0[1][0]);
  }/*if (org_isite1 + 1 > X->Def.Nsite)*/
  else {
#pragma omp parallel  default(none) shared(tmp_v0, tmp_v1, list_1) \
  firstprivate(i_max, tmp_V, X, isite1) private(j, tmp_off)
    {
#pragma omp for
      for (j = 1; j <= i_max; j++) {
        if (X_CisAis(list_1[j], X, isite1) != 0) {
          zaxpy_(&nstate, &tmp_V, &tmp_v1[j][0], &one, &tmp_v0[j][0], &one);
        }/*if (X_CisAis(list_1[j], X, isite1) != 0)*/
      }/*for (j = 1; j <= i_max; j++)*/
    }/*End of parallel region*/
  }/*if (org_isite1 + 1 <= X->Def.Nsite)*/
}/*double complex X_child_CisAis_Hubbard_MPI*/
/**
@brief Single creation/annihilation operator
 in the inter process region for HubbardGC.
@author Mitsuaki Kawamura (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
void X_GC_Cis_MPI(
  int org_isite,//!<[in] Site i
  int org_ispin,//!<[in] Spin s
  double complex tmp_trans,//!<[in] Coupling constant//!<[in]
  int nstate, double complex **tmp_v0,//!<[out] Result v0 += H v1*/,
  double complex **tmp_v1,//!<[in] v0 += H v1*/,
  unsigned long int idim_max,//!<[in] Similar to CheckList::idim_max
  double complex **tmp_v1buf,//!<[in] buffer for wavefunction
  unsigned long int *Tpow//!<[in] Similar to DefineList::Tpow
) {
  int mask2, state2, ierr, origin, bit2diff, Fsgn;
  unsigned long int idim_max_buf, j;
  double complex trans;

  // org_isite >= Nsite
  mask2 = (int)Tpow[2 * org_isite + org_ispin];

  origin = myrank ^ mask2; // XOR
  state2 = origin & mask2;

  //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in myrank
  //origin: if the state (org_isite, org_ispin) is occupied in myrank, the state is not occupied in origin.

  bit2diff = myrank - ((2 * mask2 - 1) & myrank);

  //SgnBit((unsigned long int) (origin & bit2diff), &Fsgn); // Fermion sign
  SgnBit((unsigned long int) (bit2diff), &Fsgn); // Fermion sign

  idim_max_buf = SendRecv_i(origin, idim_max);
  SendRecv_cv(origin, idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &tmp_v1buf[1][0]);

  if (state2 == mask2) {
    trans = 0;
  }
  else if (state2 == 0) {
    trans = (double)Fsgn * tmp_trans;
  }
  else return;

  zaxpy_long(idim_max_buf*nstate, trans, &tmp_v1buf[1][0], &tmp_v0[1][0]);
}/*double complex X_GC_Cis_MPI*/
/**
@brief Single creation/annihilation operator
  in the inter process region for HubbardGC.
@author Mitsuaki Kawamura (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
void X_GC_Ajt_MPI(
  int org_isite,//!<[in] Site j
  int org_ispin,//!<[in] Spin t
  double complex tmp_trans,//!<[in] Coupling constant//!<[in]
  int nstate, double complex **tmp_v0,//!<[out] Result v0 += H v1*/,
  double complex **tmp_v1,//!<[in] v0 += H v1*/,
  unsigned long int idim_max,//!<[in] Similar to CheckList::idim_max
  double complex **tmp_v1buf,//!<[in] buffer for wavefunction
  unsigned long int *Tpow//!<[in] Similar to DefineList::Tpow
) {
  int mask2, state2, ierr, origin, bit2diff, Fsgn;
  unsigned long int idim_max_buf, j;
  double complex trans;

  // org_isite >= Nsite
  mask2 = (int)Tpow[2 * org_isite + org_ispin];

  origin = myrank ^ mask2; // XOR
  state2 = origin & mask2;

  //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in myrank
  //origin: if the state (org_isite, org_ispin) is occupied in myrank, the state is not occupied in origin.

  bit2diff = myrank - ((2 * mask2 - 1) & myrank);

  //SgnBit((unsigned long int) (origin & bit2diff), &Fsgn); // Fermion sign
  SgnBit((unsigned long int) (bit2diff), &Fsgn); // Fermion sign

  idim_max_buf = SendRecv_i(origin, idim_max);

  SendRecv_cv(origin, idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &tmp_v1buf[1][0]);

  if (     state2 == 0    ) trans = 0;
  else if (state2 == mask2) trans = (double)Fsgn * tmp_trans;
  else return;

  zaxpy_long(idim_max_buf*nstate, trans, &tmp_v1buf[1][0], &tmp_v0[1][0]);
}/*double complex X_GC_Ajt_MPI*/
/**
@brief Compute @f$c_{is}^\dagger@f$
term of canonical Hubbard system
*/
void X_Cis_MPI(
  int org_isite,//!<[in] Site i
  unsigned int org_ispin,//!<[in] Spin s
  double complex tmp_trans,//!<[in] Coupling constant
  int nstate, double complex **tmp_v0,//!<[inout] Resulting wavefunction
  double complex **tmp_v1,//!<[inout] Initial wavefunction
  double complex **tmp_v1buf,//!<[in] buffer for wavefunction
  unsigned long int idim_max,//!<[in] Similar to CheckList::idim_max
  long unsigned int *Tpow,//!<[in] Similar to DefineList::Tpow
  long unsigned int *list_1_org,//!<[in] Similar to ::list_1
  long unsigned int *list_1buf_org,//!<[in] Similar to ::list_1buf
  long unsigned int *list_2_1_target,//!<[in] Similar to ::list_2_1
  long unsigned int *list_2_2_target,//!<[in] Similar to ::list_2_2
  long unsigned int _irght,//!<[in] Similer to LargeList::irght
  long unsigned int _ilft,//!<[in] Similer to LargeList::ilft
  long unsigned int _ihfbit//!<[in] Similer to LargeList::ihfbit
) {
  int mask2, state2, ierr, origin, bit2diff, Fsgn;
  unsigned long int idim_max_buf, j, ioff;
  double complex trans;
  int one = 1;

  // org_isite >= Nsite
  mask2 = (int)Tpow[2 * org_isite + org_ispin];

  origin = myrank ^ mask2; // XOR
  state2 = origin & mask2;

  //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in myrank
  //origin: if the state (org_isite, org_ispin) is occupied in myrank, the state is not occupied in origin.

  bit2diff = myrank - ((2 * mask2 - 1) & myrank);

  SgnBit((unsigned long int) (bit2diff), &Fsgn); // Fermion sign

  idim_max_buf = SendRecv_i(origin, idim_max);

  SendRecv_iv(origin, idim_max + 1, idim_max_buf + 1, list_1_org, list_1buf_org);

  SendRecv_cv(origin, idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &tmp_v1buf[1][0]);

  if (state2 == mask2) {
    trans = 0;
  }
  else if (state2 == 0) {
    trans = (double)Fsgn * tmp_trans;
  }
  else retur;

#pragma omp parallel for default(none) private(j) \
firstprivate(idim_max_buf, trans, ioff, _irght, _ilft, _ihfbit, list_2_1_target, list_2_2_target) \
shared(tmp_v1buf, tmp_v1, nstate, tmp_v0, list_1buf_org)
  for (j = 1; j <= idim_max_buf; j++) {//idim_max_buf -> original
    GetOffComp(list_2_1_target, list_2_2_target, list_1buf_org[j],
      _irght, _ilft, _ihfbit, &ioff);
    zaxpy_(&nstate, &trans, &tmp_v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
  }/*for (j = 1; j <= idim_max_buf; j++)*/
}/*double complex X_GC_Cis_MPI*/
/**
@brief Compute @f$c_{jt}@f$
term of canonical Hubbard system
*/
void X_Ajt_MPI(
  int org_isite,//!<[in] Site j
  unsigned int org_ispin,//!<[in] Spin t
  double complex tmp_trans,//!<[in] Coupling constant
  int nstate, double complex **tmp_v0,//!<[inout] Resulting wavefunction
  double complex **tmp_v1,//!<[inout] Initial wavefunction
  double complex **tmp_v1buf,//!<[in] buffer for wavefunction
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
  int mask2, state2, ierr, origin, bit2diff, Fsgn;
  unsigned long int idim_max_buf, j, ioff;
  double complex trans;
  int one = 1;

  // org_isite >= Nsite
  mask2 = (int)Tpow[2 * org_isite + org_ispin];

  origin = myrank ^ mask2; // XOR
  state2 = origin & mask2;

  //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in myrank
  //origin: if the state (org_isite, org_ispin) is occupied in myrank, the state is not occupied in origin.

  bit2diff = myrank - ((2 * mask2 - 1) & myrank);

  SgnBit((unsigned long int) (bit2diff), &Fsgn); // Fermion sign
  idim_max_buf = SendRecv_i(origin, idim_max);
  SendRecv_iv(origin, idim_max + 1, idim_max_buf + 1, list_1_org, list_1buf_org);
  SendRecv_cv(origin, idim_max*nstate, idim_max_buf*nstate, &tmp_v1[1][0], &tmp_v1buf[1][0]);

  if (state2 == 0) {
    trans = 0;
  }
  else if (state2 == mask2) {
    trans = (double)Fsgn * tmp_trans;
  }
  else return;

#pragma omp parallel for default(none) private(j) \
firstprivate(idim_max_buf, trans, ioff, _irght, _ilft, _ihfbit, list_2_1_target, list_2_2_target) \
shared(tmp_v1buf, tmp_v1, nstate, tmp_v0, list_1buf_org)
  for (j = 1; j <= idim_max_buf; j++) {
    GetOffComp(list_2_1_target, list_2_2_target, list_1buf_org[j],
      _irght, _ilft, _ihfbit, &ioff);
    zaxpy_(&nstate, &trans, &tmp_v1buf[j][0], &one, &tmp_v0[ioff][0], &one);
  }
}/*double complex X_Ajt_MPI*/
