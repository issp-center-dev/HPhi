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
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
) {
  long int idim_max, idim_maxMPI;
  long unsigned int i, j;
  long unsigned int org_isite, ispin, itype;
  long unsigned int is1_spin;
  int isgn = 1;
  double complex tmpphi;
  long unsigned int tmp_off = 0;
  //tmp_v0
  if (X->Def.NSingleExcitationOperator == 0) {
    return TRUE;
  }
  double complex *tmp_v1bufOrg;
  //set size
#ifdef MPI
  idim_maxMPI = MaxMPI_li(X->Check.idim_maxOrg);
  tmp_v1bufOrg=cd_1d_allocate(idim_maxMPI + 1);
#endif // MPI

  idim_max = X->Check.idim_maxOrg;
  for (i = 0; i < X->Def.NSingleExcitationOperator; i++) {
    org_isite = X->Def.SingleExcitationOperator[i][0];
    ispin = X->Def.SingleExcitationOperator[i][1];
    itype = X->Def.SingleExcitationOperator[i][2];
    tmpphi = X->Def.ParaSingleExcitationOperator[i];
    is1_spin = X->Def.Tpow[2 * org_isite + ispin];
    if (itype == 1) {
      if (org_isite >= X->Def.Nsite) {
        child_Cis_MPI(org_isite, ispin, tmpphi, tmp_v0, tmp_v1, tmp_v1bufOrg, idim_max, \
          X->Def.Tpow, list_1_org, list_1buf_org, list_2_1, list_2_2, \
          X->Large.irght, X->Large.ilft, X->Large.ihfbit);
      }
      else {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1, X, list_1_org) \
  firstprivate(idim_max, tmpphi, org_isite, ispin, list_2_1, list_2_2, is1_spin) private(j,  isgn,tmp_off)
        for (j = 1; j <= idim_max; j++) {//idim_max -> original dimension
          isgn = child_Cis(j, is1_spin, &tmp_off, list_1_org, list_2_1, list_2_2, X->Large.irght, X->Large.ilft, X->Large.ihfbit);
          tmp_v0[tmp_off] += tmp_v1[j] * isgn*tmpphi;
        }
      }
    }
    else if (itype == 0) {
      if (org_isite >= X->Def.Nsite) {
        child_Ajt_MPI(org_isite, ispin, tmpphi, tmp_v0, tmp_v1, tmp_v1bufOrg, \
          idim_max, X->Def.Tpow, list_1_org, list_1buf_org, \
          list_2_1, list_2_2, X->Large.irght, X->Large.ilft, X->Large.ihfbit);
      }
      else {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1, X, list_1_org, list_1) \
  firstprivate(idim_max, tmpphi, org_isite, ispin, list_2_1, list_2_2, is1_spin, myrank) private(j, isgn, tmp_off)
        for (j = 1; j <= idim_max; j++) {//idim_max -> original dimension
          isgn = child_Ajt(j, is1_spin, &tmp_off, list_1_org, list_2_1, list_2_2, X->Large.irght, X->Large.ilft, X->Large.ihfbit);
          tmp_v0[tmp_off] += tmp_v1[j] * isgn*tmpphi;
        }
      }
    }
  }
#ifdef MPI
  free_cd_1d_allocate(tmp_v1bufOrg);
#endif
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
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
  long int idim_max, idim_maxMPI;
  long unsigned int i, j;
  long unsigned int org_isite, ispin, itype;
  long unsigned int is1_spin;
  double complex tmpphi;
  long unsigned int tmp_off = 0;
  //idim_max = X->Check.idim_max;
  idim_max = X->Check.idim_maxOrg;
  //tmp_v0
  if (X->Def.NSingleExcitationOperator == 0) {
    return TRUE;
  }
  double complex *tmp_v1bufOrg;
  //set size
#ifdef MPI
  idim_maxMPI = MaxMPI_li(X->Check.idim_maxOrg);
  tmp_v1bufOrg=cd_1d_allocate(idim_maxMPI + 1);
#endif // MPI

  // SingleEx
  for (i = 0; i < X->Def.NSingleExcitationOperator; i++) {
    org_isite = X->Def.SingleExcitationOperator[i][0];
    ispin = X->Def.SingleExcitationOperator[i][1];
    itype = X->Def.SingleExcitationOperator[i][2];
    tmpphi = X->Def.ParaSingleExcitationOperator[i];
    if (itype == 1) {
      if (org_isite >= X->Def.Nsite) {
        child_GC_Cis_MPI(org_isite, ispin, tmpphi, tmp_v0, tmp_v1, idim_max, tmp_v1bufOrg, X->Def.Tpow);
      }
      else {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1, X) \
  firstprivate(idim_max, tmpphi, org_isite, ispin) private(j, is1_spin, tmp_off)
        for (j = 1; j <= idim_max; j++) {
          is1_spin = X->Def.Tpow[2 * org_isite + ispin];
          GC_Cis(j, tmp_v0, tmp_v1, is1_spin, tmpphi, &tmp_off);
        }/*for (j = 1; j <= idim_max; j++)*/
      }
    }
    else if (itype == 0) {
      if (org_isite >= X->Def.Nsite) {
        child_GC_Ajt_MPI(org_isite, ispin, tmpphi, tmp_v0, tmp_v1, idim_max, tmp_v1bufOrg, X->Def.Tpow);
      }
      else {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1, X) \
  firstprivate(idim_max, tmpphi, org_isite, ispin) private(j, is1_spin, tmp_off)
        for (j = 1; j <= idim_max; j++) {
          is1_spin = X->Def.Tpow[2 * org_isite + ispin];
          GC_Ajt(j, tmp_v0, tmp_v1, is1_spin, tmpphi, &tmp_off);
        }/*for (j = 1; j <= idim_max; j++)*/
      }
    }
  }
#ifdef MPI
  free_cd_1d_allocate(tmp_v1bufOrg);
#endif
  return TRUE;
}/*int GetSingleExcitedStateHubbardGC*/
