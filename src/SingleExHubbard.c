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
#include "bitcalc.h"
#include "wrapperMPI.h"
#include "common/setmemory.h"
#include "mltplyHubbardCore.h"
#include "mltplyMPIHubbardCore.h"
#include "mltplyCommon.h"
#ifdef MPI
#include "mpi.h"
#endif
/**@file
@brief Functions to compute singly excited state in Hubbard model
*/
/**
@brief Calculation of Single excited state for Hubbard canonical system
@returns TRUE: Normally finished
@returns FALSE: Abnormally finished
@author Kazuyoshi Yoshimi
@version 1.2
*/
int GetSingleExcitedStateHubbard(
  struct BindStruct *X,//!<define list to get and put information of calculation
  int nstate, 
  double complex **tmp_v0,//!<[out] Result v0 = H v1
  double complex **tmp_v1,//!<[in] v0 = H v1
  int iEx
) {
  long unsigned int idim_max;
  long unsigned int i, j;
  long unsigned int org_isite, ispin, itype;
  long unsigned int is1_spin;
  int isgn = 1, one = 1;
  double complex tmpphi, dmv;
  long unsigned int tmp_off = 0;
  //tmp_v0
  if (X->Def.NSingleExcitationOperator[iEx] == 0) {
    return TRUE;
  }

  idim_max = X->Check.idim_maxOrg;
  for (i = 0; i < X->Def.NSingleExcitationOperator[iEx]; i++) {
    org_isite = X->Def.SingleExcitationOperator[iEx][i][0];
    ispin = X->Def.SingleExcitationOperator[iEx][i][1];
    itype = X->Def.SingleExcitationOperator[iEx][i][2];
    tmpphi = X->Def.ParaSingleExcitationOperator[iEx][i];
    is1_spin = X->Def.Tpow[2 * org_isite + ispin];
    if (itype == 1) {
      if (org_isite >= X->Def.Nsite) {
        X_Cis_MPI(org_isite, ispin, tmpphi, nstate, tmp_v0, tmp_v1, idim_max, 
          X->Def.Tpow, 
          X->Large.irght, X->Large.ilft, X->Large.ihfbit);
      }
      else {
#pragma omp parallel for default(none) shared(nstate,tmp_v0, tmp_v1, X, list_1_org,one) \
firstprivate(idim_max, tmpphi, org_isite, ispin, list_2_1, list_2_2, is1_spin) \
private(j,  isgn,tmp_off,dmv)
        for (j = 1; j <= idim_max; j++) {//idim_max -> original dimension
          isgn = X_Cis(j, is1_spin, &tmp_off, list_1_org, list_2_1, list_2_2,
            X->Large.irght, X->Large.ilft, X->Large.ihfbit);
          dmv = isgn * tmpphi;
          zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[tmp_off], &one);
        }
      }
    }
    else if (itype == 0) {
      if (org_isite >= X->Def.Nsite) {
        X_Ajt_MPI(org_isite, ispin, tmpphi, nstate, tmp_v0, tmp_v1, 
          idim_max, X->Def.Tpow, X->Large.irght, X->Large.ilft, X->Large.ihfbit);
      }
      else {
#pragma omp parallel for default(none) shared(tmp_v0,tmp_v1,X,list_1_org,list_1,one,nstate) \
firstprivate(idim_max,tmpphi,org_isite,ispin,list_2_1,list_2_2,is1_spin,myrank) \
private(j, isgn, tmp_off,dmv)
        for (j = 1; j <= idim_max; j++) {//idim_max -> original dimension
          isgn = X_Ajt(j, is1_spin, &tmp_off, list_1_org, list_2_1, list_2_2, 
            X->Large.irght, X->Large.ilft, X->Large.ihfbit);
          dmv = isgn * tmpphi;
          zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[tmp_off], &one);
        }
      }
    }
  }
  return TRUE;
}/*int GetSingleExcitedStateHubbard*/
/**
@brief Calculation of Single excited state for Hubbard Grand canonical system
@returns TRUE: Normally finished
@returns FALSE: Abnormally finished
@author Kazuyoshi Yoshimi
@version 1.2
*/
int GetSingleExcitedStateHubbardGC(
  struct BindStruct *X,//!<define list to get and put information of calculation
  int nstate, 
  double complex **tmp_v0,//!<[out] Result v0 = H v1
  double complex **tmp_v1,//!<[in] v0 = H v1
  int iEx
) {
  long unsigned int idim_max;
  long unsigned int i, j;
  long unsigned int org_isite, ispin, itype;
  long unsigned int is1_spin;
  double complex tmpphi;
  long unsigned int tmp_off = 0;
  //idim_max = X->Check.idim_max;
  idim_max = X->Check.idim_maxOrg;
  //tmp_v0
  if (X->Def.NSingleExcitationOperator[iEx] == 0) {
    return TRUE;
  }

  // SingleEx
  for (i = 0; i < X->Def.NSingleExcitationOperator[iEx]; i++) {
    org_isite = X->Def.SingleExcitationOperator[iEx][i][0];
    ispin = X->Def.SingleExcitationOperator[iEx][i][1];
    itype = X->Def.SingleExcitationOperator[iEx][i][2];
    tmpphi = X->Def.ParaSingleExcitationOperator[iEx][i];
    if (itype == 1) {
      if (org_isite >= X->Def.Nsite) {
        X_GC_Cis_MPI(org_isite, ispin, tmpphi, nstate, tmp_v0, tmp_v1, 
          idim_max, X->Def.Tpow);
      }
      else {
#pragma omp parallel for default(none) shared(tmp_v0,tmp_v1,X,nstate)         \
firstprivate(idim_max, tmpphi, org_isite, ispin) private(j, is1_spin, tmp_off)
        for (j = 1; j <= idim_max; j++) {
          is1_spin = X->Def.Tpow[2 * org_isite + ispin];
          GC_Cis(j, nstate, tmp_v0, tmp_v1, is1_spin, tmpphi, &tmp_off);
        }/*for (j = 1; j <= idim_max; j++)*/
      }
    }
    else if (itype == 0) {
      if (org_isite >= X->Def.Nsite) {
        X_GC_Ajt_MPI(org_isite, ispin, tmpphi, nstate, tmp_v0, tmp_v1, 
          idim_max, X->Def.Tpow);
      }
      else {
#pragma omp parallel for default(none) shared(tmp_v0,tmp_v1,X,nstate)         \
firstprivate(idim_max, tmpphi, org_isite, ispin) private(j, is1_spin, tmp_off)
        for (j = 1; j <= idim_max; j++) {
          is1_spin = X->Def.Tpow[2 * org_isite + ispin];
          GC_Ajt(j, nstate, tmp_v0, tmp_v1, is1_spin, tmpphi, &tmp_off);
        }/*for (j = 1; j <= idim_max; j++)*/
      }
    }
  }
  return TRUE;
}/*int GetSingleExcitedStateHubbardGC*/
