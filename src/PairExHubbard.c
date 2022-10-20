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
/*-------------------------------------------------------------*/
#include "bitcalc.h"
#include "wrapperMPI.h"
#include "mltplyCommon.h"
#include "mltplyHubbard.h"
#include "mltplyHubbardCore.h"
#include "mltplyMPIHubbard.h"
#include "mltplyMPIHubbardCore.h"
#ifdef MPI
#include "common/setmemory.h"
#endif

///
/// Calculation of pair excited state for Hubbard Grand canonical system
/// \param X [in,out] define list to get and put information of calculation
/// \param tmp_v0 [out] Result v0 = H v1
/// \param tmp_v1 [in] v0 = H v1
/// \returns TRUE: Normally finished
/// \returns FALSE: Abnormally finished
/// \author Kazuyoshi Yoshimi
/// \version 1.2
int GetPairExcitedStateHubbardGC(
  struct BindStruct *X,/**< [inout] define list to get and put information of calculation*/
  int nstate, 
  double complex **tmp_v0, /**< [out] Result v0 = H v1*/
  double complex **tmp_v1, /**< [in] v0 = H v1*/
  int iEx
) {
  long unsigned int i, j;
  long unsigned int isite1;
  long unsigned int org_isite1, org_isite2, org_sigma1, org_sigma2;

  double complex tmp_trans = 0;
  long unsigned int i_max;
  long unsigned int ibit;
  long unsigned int is;
  i_max = X->Check.idim_maxOrg;
  for (i = 0; i < X->Def.NPairExcitationOperator[iEx]; i++) {
    org_isite1 = X->Def.PairExcitationOperator[iEx][i][0] + 1;
    org_isite2 = X->Def.PairExcitationOperator[iEx][i][2] + 1;
    org_sigma1 = X->Def.PairExcitationOperator[iEx][i][1];
    org_sigma2 = X->Def.PairExcitationOperator[iEx][i][3];
    tmp_trans = X->Def.ParaPairExcitationOperator[iEx][i];

    if (org_isite1 > X->Def.Nsite &&
      org_isite2 > X->Def.Nsite) {
      if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {
        if (X->Def.PairExcitationOperator[iEx][i][4] == 0) {
          if (org_sigma1 == 0) {
            is = X->Def.Tpow[2 * org_isite1 - 2];
          }
          else {
            is = X->Def.Tpow[2 * org_isite1 - 1];
          }
          ibit = (unsigned long int) myrank & is;
          if (ibit != is) {
            //minus sign comes from negative tmp_trans due to readdef
            zaxpy_long(i_max*nstate, -tmp_trans, &tmp_v1[1][0], &tmp_v0[1][0]);
          }
        }
        else {//X->Def.PairExcitationOperator[iEx][i][4]==1
          if (org_sigma1 == 0) {
            is = X->Def.Tpow[2 * org_isite1 - 2];
          }
          else {
            is = X->Def.Tpow[2 * org_isite1 - 1];
          }
          ibit = (unsigned long int) myrank & is;
          if (ibit == is) {
            zaxpy_long(i_max*nstate, tmp_trans, &tmp_v1[1][0], &tmp_v0[1][0]);
          }
        }
      }
      else {
            child_GC_general_hopp_MPIdouble(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          -tmp_trans, X, nstate, tmp_v0, tmp_v1);
      }
    }
    else if (org_isite2 > X->Def.Nsite || org_isite1 > X->Def.Nsite) {
      if (org_isite1 < org_isite2) {
        child_GC_general_hopp_MPIsingle(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          -tmp_trans, X, nstate, tmp_v0, tmp_v1);
      }
      else {
        child_GC_general_hopp_MPIsingle(org_isite2 - 1, org_sigma2, org_isite1 - 1, org_sigma1,
          -conj(tmp_trans), X, nstate, tmp_v0, tmp_v1);
      }
    }
    else {

      if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2 && X->Def.PairExcitationOperator[iEx][i][4] == 0) {
        isite1 = X->Def.Tpow[2 * org_isite1 - 2 + org_sigma1];
#pragma omp parallel for default(none) private(j) \
firstprivate(i_max,X,isite1, tmp_trans) shared(tmp_v0,tmp_v1,nstate)
        for (j = 1; j <= i_max; j++) {
          GC_AisCis(j, nstate, tmp_v0, tmp_v1, isite1, -tmp_trans);
        }
      }
      else {
        if (general_hopp_GetInfo(X, org_isite1, org_isite2, org_sigma1, org_sigma2) != 0) {
          return -1;
        }
        GC_general_hopp(nstate, tmp_v0, tmp_v1, X, tmp_trans);
      }
    }
  }
  return TRUE;
}
///
/// Calculation of pair excited state for Hubbard canonical system
/// \param X [in,out] define list to get and put information of calculation
/// \param tmp_v0 [out] Result v0 = H v1
/// \param tmp_v1 [in] v0 = H v1
/// \returns TRUE: Normally finished
/// \returns FALSE: Abnormally finished
/// \author Kazuyoshi Yoshimi
/// \version 1.2
int GetPairExcitedStateHubbard(
  struct BindStruct *X,/**< [inout] define list to get and put information of calculation*/
  int nstate,
  double complex **tmp_v0, /**< [out] Result v0 = H v1*/
  double complex **tmp_v1, /**< [in] v0 = H v1*/
  int iEx
) {
  long unsigned int i, j;
  long unsigned int irght, ilft, ihfbit;
  long unsigned int org_isite1, org_isite2, org_sigma1, org_sigma2;
  long unsigned int tmp_off = 0;

  double complex tmp_trans = 0, dmv;
  long unsigned int i_max;
  int tmp_sgn, num1, one = 1;
  long unsigned int ibit;
  long unsigned int is, Asum, Adiff;
  long unsigned int ibitsite1, ibitsite2;

  //  i_max = X->Check.idim_max;
  i_max = X->Check.idim_maxOrg;
  if (GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit) != 0) {
    return -1;
  }
  X->Large.i_max = i_max;
  X->Large.irght = irght;
  X->Large.ilft = ilft;
  X->Large.ihfbit = ihfbit;
  X->Large.mode = M_CALCSPEC;

  for (i = 0; i < X->Def.NPairExcitationOperator[iEx]; i++) {
    org_isite1 = X->Def.PairExcitationOperator[iEx][i][0] + 1;
    org_isite2 = X->Def.PairExcitationOperator[iEx][i][2] + 1;
    org_sigma1 = X->Def.PairExcitationOperator[iEx][i][1];
    org_sigma2 = X->Def.PairExcitationOperator[iEx][i][3];
    tmp_trans = X->Def.ParaPairExcitationOperator[iEx][i];
    ibitsite1 = X->Def.OrgTpow[2 * org_isite1 - 2 + org_sigma1];
    ibitsite2 = X->Def.OrgTpow[2 * org_isite2 - 2 + org_sigma2];
    general_hopp_GetInfo(X, org_isite1, org_isite2, org_sigma1, org_sigma2);
    Asum = X->Large.isA_spin;
    Adiff = X->Large.A_spin;

    if (X->Def.iFlagListModified == TRUE // Not to adopt HubbrdNConserved
      && org_sigma1 != org_sigma2) {
      if (org_isite1 > X->Def.Nsite &&
        org_isite2 > X->Def.Nsite)
      {
        child_CisAjt_MPIdouble(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          -tmp_trans, X, nstate, tmp_v0, tmp_v1);
      }
      else if (org_isite2 > X->Def.Nsite
        || org_isite1 > X->Def.Nsite)
      {
        if (org_isite1 < org_isite2) {
          child_CisAjt_MPIsingle(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
            -tmp_trans, X, nstate, tmp_v0, tmp_v1);
        }
        else {
          child_CisAjt_MPIsingle(org_isite2 - 1, org_sigma2, org_isite1 - 1, org_sigma1,
            -conj(tmp_trans), X, nstate, tmp_v0, tmp_v1);
        }
      }
      else {
#pragma omp parallel for default(none) shared(tmp_v0,tmp_v1,one,nstate)       \
firstprivate(i_max,tmp_trans,Asum,Adiff,ibitsite1,ibitsite2,X,list_1_org,list_1,myrank) \
  private(j,tmp_sgn,tmp_off,dmv)
        for (j = 1; j <= i_max; j++) {
          tmp_sgn = child_CisAjt(list_1_org[j], X, ibitsite1, ibitsite2, Asum, Adiff, &tmp_off);
          dmv = tmp_trans * tmp_sgn;
          zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[tmp_off], &one);
        }
      }
    }
    else {
      if (org_isite1 > X->Def.Nsite &&
        org_isite2 > X->Def.Nsite) {
        if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {//diagonal
          is = X->Def.Tpow[2 * org_isite1 - 2 + org_sigma1];
          ibit = (unsigned long int) myrank & is;
          if (X->Def.PairExcitationOperator[iEx][i][4] == 0) {
            if (ibit != is) {
              dmv = -tmp_trans;
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1,one,dmv,nstate) \
firstprivate(i_max, tmp_trans) private(j)
              for (j = 1; j <= i_max; j++) {
                zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
              }
            }
          }
          else {
            if (ibit == is) {
              zaxpy_long(i_max*nstate, tmp_trans, &tmp_v1[1][0], &tmp_v0[1][0]);
            }
          }
        }
        else {
          child_general_hopp_MPIdouble(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
            -tmp_trans, X, nstate, tmp_v0, tmp_v1);
        }
      }
      else if (org_isite2 > X->Def.Nsite || org_isite1 > X->Def.Nsite) {
        if (org_isite1 < org_isite2) {
          child_general_hopp_MPIsingle(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
            -tmp_trans, X, nstate, tmp_v0, tmp_v1);
        }
        else {
          child_general_hopp_MPIsingle(org_isite2 - 1, org_sigma2, org_isite1 - 1, org_sigma1,
            -conj(tmp_trans), X, nstate, tmp_v0, tmp_v1);
        }
      }
      else {
        if (general_hopp_GetInfo(X, org_isite1, org_isite2, org_sigma1, org_sigma2) != 0) {
          return -1;
        }
        if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {
          is = X->Def.Tpow[2 * org_isite1 - 2 + org_sigma1];
          if (X->Def.PairExcitationOperator[iEx][i][4] == 0) {
#pragma omp parallel for default(none) shared(list_1,nstate,tmp_v0,tmp_v1,one) \
firstprivate(i_max,is,tmp_trans) private(num1,ibit,dmv)
            for (j = 1; j <= i_max; j++) {
              ibit = list_1[j] & is;
              num1 = (1 - ibit / is);
              dmv = -tmp_trans * num1;
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
            }
          }
          else {
#pragma omp parallel for default(none) shared(list_1,nstate,tmp_v0,tmp_v1,one) \
  firstprivate(i_max,is,tmp_trans) private(num1,ibit,dmv)
            for (j = 1; j <= i_max; j++) {
              ibit = list_1[j] & is;
              num1 = ibit / is;
              dmv = tmp_trans * num1;
              zaxpy_(&nstate, &dmv, tmp_v1[j], &one, tmp_v0[j], &one);
            }
          }
        }
        else {
          general_hopp(nstate, tmp_v0, tmp_v1, X, tmp_trans);
        }
      }
    }
  }
  return TRUE;
}
