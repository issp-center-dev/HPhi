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

#include "mltply.h"
#include "mltplyCommon.h"
#include "FileIO.h"
#include "bitcalc.h"
#include "expec_cisajscktaltdc.h"
#include "mltplySpin.h"
#include "mltplySpinCore.h"
#include "mltplyHubbardCore.h"
#include "wrapperMPI.h"
#include "mltplyMPISpin.h"
#include "mltplyMPISpinCore.h"
#include "mltplyMPIHubbardCore.h"
#include "common/setmemory.h"
/**
 * @file   expec_cisajscktaltdc.c
 * 
 * @brief  File for calculating two-body green's functions
 * 
 * @version 0.2
 * @details add function to treat the case of general spin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 */
///
/// \brief Rearray interactions
/// \param i
/// \param org_isite1 a site number on the site 1.
/// \param org_isite2 a site number on the site 2.
/// \param org_isite3 a site number on the site 3.
/// \param org_isite4 a site number on the site 4.
/// \param org_sigma1 a spin index on the site 1.
/// \param org_sigma2 a spin index on the site 2.
/// \param org_sigma3 a spin index on the site 3.
/// \param org_sigma4 a spin index on the site 4.
/// \param tmp_V a value of interaction
/// \param X  data list for calculation
/// \return 0 normally finished
/// \return -1 unnormally finished
int Rearray_Interactions(
                         int i,
                         long unsigned int *org_isite1,
                         long unsigned int *org_isite2,
                         long unsigned int *org_isite3,
                         long unsigned int *org_isite4,
                         long unsigned int *org_sigma1,
                         long unsigned int *org_sigma2,
                         long unsigned int *org_sigma3,
                         long unsigned int *org_sigma4,
                         double complex *tmp_V,
                         struct BindStruct *X,
                         int type
                         )
{
  long unsigned int tmp_org_isite1,tmp_org_isite2,tmp_org_isite3,tmp_org_isite4;
  long unsigned int tmp_org_sigma1,tmp_org_sigma2,tmp_org_sigma3,tmp_org_sigma4;

  if(type==6){
    tmp_org_isite1   = X->Def.SBody[i][0]+1;
    tmp_org_sigma1   = X->Def.SBody[i][1];
    tmp_org_isite2   = X->Def.SBody[i][2]+1;
    tmp_org_sigma2   = X->Def.SBody[i][3];
    tmp_org_isite3   = X->Def.SBody[i][4]+1;
    tmp_org_sigma3   = X->Def.SBody[i][5];
    tmp_org_isite4   = X->Def.SBody[i][6]+1;
    tmp_org_sigma4   = X->Def.SBody[i][7];
  }else if(type==4){
    tmp_org_isite1   = X->Def.FBody[i][0]+1;
    tmp_org_sigma1   = X->Def.FBody[i][1];
    tmp_org_isite2   = X->Def.FBody[i][2]+1;
    tmp_org_sigma2   = X->Def.FBody[i][3];
    tmp_org_isite3   = X->Def.FBody[i][4]+1;
    tmp_org_sigma3   = X->Def.FBody[i][5];
    tmp_org_isite4   = X->Def.FBody[i][6]+1;
    tmp_org_sigma4   = X->Def.FBody[i][7];
  }else if(type==3){
    tmp_org_isite1   = X->Def.TBody[i][0]+1;
    tmp_org_sigma1   = X->Def.TBody[i][1];
    tmp_org_isite2   = X->Def.TBody[i][2]+1;
    tmp_org_sigma2   = X->Def.TBody[i][3];
    tmp_org_isite3   = X->Def.TBody[i][4]+1;
    tmp_org_sigma3   = X->Def.TBody[i][5];
    tmp_org_isite4   = X->Def.TBody[i][6]+1;
    tmp_org_sigma4   = X->Def.TBody[i][7];
  }else{
    tmp_org_isite1   = X->Def.CisAjtCkuAlvDC[i][0]+1;
    tmp_org_sigma1   = X->Def.CisAjtCkuAlvDC[i][1];
    tmp_org_isite2   = X->Def.CisAjtCkuAlvDC[i][2]+1;
    tmp_org_sigma2   = X->Def.CisAjtCkuAlvDC[i][3];
    tmp_org_isite3   = X->Def.CisAjtCkuAlvDC[i][4]+1;
    tmp_org_sigma3   = X->Def.CisAjtCkuAlvDC[i][5];
    tmp_org_isite4   = X->Def.CisAjtCkuAlvDC[i][6]+1;
    tmp_org_sigma4   = X->Def.CisAjtCkuAlvDC[i][7];
  }

  if(tmp_org_isite1==tmp_org_isite2 && tmp_org_isite3==tmp_org_isite4){
    if(tmp_org_isite1 > tmp_org_isite3){
      *org_isite1   = tmp_org_isite3;
      *org_sigma1   = tmp_org_sigma3;
      *org_isite2   = tmp_org_isite4;
      *org_sigma2   = tmp_org_sigma4;
      *org_isite3   = tmp_org_isite1;
      *org_sigma3   = tmp_org_sigma1;
      *org_isite4   = tmp_org_isite2;
      *org_sigma4   = tmp_org_sigma2;
    }
    else{
      *org_isite1   = tmp_org_isite1;
      *org_sigma1   = tmp_org_sigma1;
      *org_isite2   = tmp_org_isite2;
      *org_sigma2   = tmp_org_sigma2;
      *org_isite3   = tmp_org_isite3;
      *org_sigma3   = tmp_org_sigma3;
      *org_isite4   = tmp_org_isite4;
      *org_sigma4   = tmp_org_sigma4;
    }
    *tmp_V = 1.0;

  }
  else if(tmp_org_isite1==tmp_org_isite4 && tmp_org_isite3==tmp_org_isite2){
    if(tmp_org_isite1 > tmp_org_isite3){
      *org_isite1   = tmp_org_isite3;
      *org_sigma1   = tmp_org_sigma3;
      *org_isite2   = tmp_org_isite2;
      *org_sigma2   = tmp_org_sigma2;
      *org_isite3   = tmp_org_isite1;
      *org_sigma3   = tmp_org_sigma1;
      *org_isite4   = tmp_org_isite4;
      *org_sigma4   = tmp_org_sigma4;
    }
    else{
      *org_isite1   = tmp_org_isite1;
      *org_sigma1   = tmp_org_sigma1;
      *org_isite2   = tmp_org_isite4;
      *org_sigma2   = tmp_org_sigma4;
      *org_isite3   = tmp_org_isite3;
      *org_sigma3   = tmp_org_sigma3;
      *org_isite4   = tmp_org_isite2;
      *org_sigma4   = tmp_org_sigma2;
    }
    *tmp_V =-1.0;
  }
  else{
    return -1;
  }
  return 0;
}
/**
 * @brief Child function to calculate two-body green's functions for Hubbard GC model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
int expec_cisajscktalt_HubbardGC(
  struct BindStruct *X, 
  int nstate, 
  double complex **Xvec,
  double complex **vec, 
  double complex **prod
) {
  long unsigned int i, j;
  long unsigned int isite1, isite2, isite3, isite4;
  long unsigned int org_isite1, org_isite2, org_isite3, org_isite4;
  long unsigned int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long unsigned int Asum, Bsum, Adiff, Bdiff;
  long unsigned int tmp_off = 0;
  long unsigned int tmp_off_2 = 0;
  double complex tmp_V = 1.0 + 0.0*I;
  long unsigned int i_max;

  for (i = 0; i < X->Def.NCisAjtCkuAlvDC; i++) {
    zclear(X->Large.i_max*nstate, &Xvec[1][0]);
    org_isite1 = X->Def.CisAjtCkuAlvDC[i][0] + 1;
    org_sigma1 = X->Def.CisAjtCkuAlvDC[i][1];
    org_isite2 = X->Def.CisAjtCkuAlvDC[i][2] + 1;
    org_sigma2 = X->Def.CisAjtCkuAlvDC[i][3];
    org_isite3 = X->Def.CisAjtCkuAlvDC[i][4] + 1;
    org_sigma3 = X->Def.CisAjtCkuAlvDC[i][5];
    org_isite4 = X->Def.CisAjtCkuAlvDC[i][6] + 1;
    org_sigma4 = X->Def.CisAjtCkuAlvDC[i][7];

    if (CheckPE(org_isite1 - 1, X) == TRUE || CheckPE(org_isite2 - 1, X) == TRUE ||
      CheckPE(org_isite3 - 1, X) == TRUE || CheckPE(org_isite4 - 1, X) == TRUE) {
      isite1 = X->Def.OrgTpow[2 * org_isite1 - 2 + org_sigma1];
      isite2 = X->Def.OrgTpow[2 * org_isite2 - 2 + org_sigma2];
      isite3 = X->Def.OrgTpow[2 * org_isite3 - 2 + org_sigma3];
      isite4 = X->Def.OrgTpow[2 * org_isite4 - 2 + org_sigma4];
      if (isite1 == isite2 && isite3 == isite4) {
        child_GC_CisAisCjtAjt_Hubbard_MPI(org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3,
          1.0, X, nstate, Xvec, vec);
      }
      else if (isite1 == isite2 && isite3 != isite4) {
        child_GC_CisAisCjtAku_Hubbard_MPI(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, org_isite4 - 1, org_sigma4,
          1.0, X, nstate, Xvec, vec);
      }
      else if (isite1 != isite2 && isite3 == isite4) {
        child_GC_CisAjtCkuAku_Hubbard_MPI(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          org_isite3 - 1, org_sigma3, 1.0, X, nstate, Xvec, vec);
      }
      else if (isite1 != isite2 && isite3 != isite4) {
        child_GC_CisAjtCkuAlv_Hubbard_MPI(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          org_isite3 - 1, org_sigma3, org_isite4 - 1, org_sigma4, 1.0, X, nstate, Xvec, vec);
      }
    }//InterPE
    else {
      general_int_GetInfo(i, X, org_isite1, org_isite2, org_isite3, org_isite4,
        org_sigma1, org_sigma2, org_sigma3, org_sigma4, tmp_V);

      i_max = X->Large.i_max;
      isite1 = X->Large.is1_spin;
      isite2 = X->Large.is2_spin;
      Asum = X->Large.isA_spin;
      Adiff = X->Large.A_spin;

      isite3 = X->Large.is3_spin;
      isite4 = X->Large.is4_spin;
      Bsum = X->Large.isB_spin;
      Bdiff = X->Large.B_spin;

      if (isite1 == isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) private(j) shared(vec,Xvec,nstate)     \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V)
        for (j = 1; j <= i_max; j++) {
          GC_CisAisCisAis_element(j, isite1, isite3, tmp_V, nstate, Xvec, vec, X, &tmp_off);
        }
      }
      else if (isite1 == isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) private(j) shared(vec,Xvec,nstate) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V)
        for (j = 1; j <= i_max; j++) {
          GC_CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff,
            tmp_V, nstate, Xvec, vec, X, &tmp_off);
        }
      }
      else if (isite1 != isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) private(j) shared(vec,Xvec,nstate) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) 
        for (j = 1; j <= i_max; j++) {
          GC_CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff,
            tmp_V, nstate, Xvec, vec, X, &tmp_off);
        }
      }
      else if (isite1 != isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) private(j) shared(vec,Xvec,nstate) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) 
        for (j = 1; j <= i_max; j++) {
          GC_CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff,
            tmp_V, nstate, Xvec, vec, X, &tmp_off_2);
        }
      }
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }//Intra PE
  return 0;
}
/**
 * @brief Child function to calculate two-body green's functions for Hubbard model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
int expec_cisajscktalt_Hubbard(
  struct BindStruct *X,
  int nstate, 
  double complex **Xvec, 
  double complex **vec,
  double complex **prod
){
  long unsigned int i, j;
  long unsigned int isite1, isite2, isite3, isite4;
  long unsigned int org_isite1, org_isite2, org_isite3, org_isite4;
  long unsigned int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long unsigned int Asum, Bsum, Adiff, Bdiff;
  long unsigned int tmp_off = 0;
  long unsigned int tmp_off_2 = 0;
  double complex tmp_V;
  long unsigned int i_max;

  for (i = 0; i < X->Def.NCisAjtCkuAlvDC; i++) {
    zclear(X->Large.i_max*nstate, &Xvec[1][0]);
    org_isite1 = X->Def.CisAjtCkuAlvDC[i][0] + 1;
    org_sigma1 = X->Def.CisAjtCkuAlvDC[i][1];
    org_isite2 = X->Def.CisAjtCkuAlvDC[i][2] + 1;
    org_sigma2 = X->Def.CisAjtCkuAlvDC[i][3];
    org_isite3 = X->Def.CisAjtCkuAlvDC[i][4] + 1;
    org_sigma3 = X->Def.CisAjtCkuAlvDC[i][5];
    org_isite4 = X->Def.CisAjtCkuAlvDC[i][6] + 1;
    org_sigma4 = X->Def.CisAjtCkuAlvDC[i][7];
    tmp_V = 1.0;

    if (X->Def.iFlgSzConserved == TRUE) {
      if (org_sigma1 + org_sigma3 != org_sigma2 + org_sigma4) {
        zclear(nstate, prod[i]);
        continue;
      }
    }

    if (CheckPE(org_isite1 - 1, X) == TRUE || CheckPE(org_isite2 - 1, X) == TRUE ||
      CheckPE(org_isite3 - 1, X) == TRUE || CheckPE(org_isite4 - 1, X) == TRUE) {
      isite1 = X->Def.OrgTpow[2 * org_isite1 - 2 + org_sigma1];
      isite2 = X->Def.OrgTpow[2 * org_isite2 - 2 + org_sigma2];
      isite3 = X->Def.OrgTpow[2 * org_isite3 - 2 + org_sigma3];
      isite4 = X->Def.OrgTpow[2 * org_isite4 - 2 + org_sigma4];
      if (isite1 == isite2 && isite3 == isite4) {
        child_CisAisCjtAjt_Hubbard_MPI(org_isite1 - 1, org_sigma1,
          org_isite3 - 1, org_sigma3, 1.0, X, nstate, Xvec, vec);
      }
      else if (isite1 == isite2 && isite3 != isite4) {
        //printf("org_isite1=%d, org_isite2=%d, org_isite3=%d, org_isite4=%d\n", org_isite1, org_isite2, org_isite3, org_isite4);
        child_CisAisCjtAku_Hubbard_MPI(org_isite1 - 1, org_sigma1,
          org_isite3 - 1, org_sigma3, org_isite4 - 1, org_sigma4, 1.0, X, nstate, Xvec, vec);
      }
      else if (isite1 != isite2 && isite3 == isite4) {
        child_CisAjtCkuAku_Hubbard_MPI(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          org_isite3 - 1, org_sigma3, 1.0, X, nstate, Xvec, vec);

      }
      else if (isite1 != isite2 && isite3 != isite4) {
        child_CisAjtCkuAlv_Hubbard_MPI(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          org_isite3 - 1, org_sigma3, org_isite4 - 1, org_sigma4, 1.0, X, nstate, Xvec, vec);
      }
    }//InterPE
    else {
      general_int_GetInfo(
        i, X, org_isite1, org_isite2, org_isite3, org_isite4,
        org_sigma1, org_sigma2, org_sigma3, org_sigma4, tmp_V
      );

      i_max = X->Large.i_max;
      isite1 = X->Large.is1_spin;
      isite2 = X->Large.is2_spin;
      Asum = X->Large.isA_spin;
      Adiff = X->Large.A_spin;

      isite3 = X->Large.is3_spin;
      isite4 = X->Large.is4_spin;
      Bsum = X->Large.isB_spin;
      Bdiff = X->Large.B_spin;

      tmp_V = 1.0;
      if (isite1 == isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) private(j) shared(vec,tmp_V,Xvec,nstate) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2)
        for (j = 1; j <= i_max; j++) {
          CisAisCisAis_element(j, isite1, isite3, tmp_V, nstate, Xvec, vec, X, &tmp_off);
        }
      }
      else if (isite1 == isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) private(j) shared(vec,tmp_V,Xvec,nstate) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2)
        for (j = 1; j <= i_max; j++) {
          CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, 
            tmp_V, nstate, Xvec, vec, X, &tmp_off);
        }
      }
      else if (isite1 != isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) private(j) shared(vec,tmp_V,Xvec,nstate) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2)
        for (j = 1; j <= i_max; j++) {
          CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, 
            tmp_V, nstate, Xvec, vec, X, &tmp_off);
        }
      }
      else if (isite1 != isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) private(j) shared(vec,tmp_V,Xvec,nstate) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2)
        for (j = 1; j <= i_max; j++) {
          CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, 
            tmp_V, nstate, Xvec, vec, X, &tmp_off_2);
        }
      }
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief Child function to calculate two-body green's functions for 1/2 Spin model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
int expec_cisajscktalt_SpinHalf(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec,
  double complex **prod
){
  long unsigned int i, j;
  long unsigned int org_isite1, org_isite2, org_isite3, org_isite4;
  long unsigned int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long unsigned int isA_up, isB_up;
  long unsigned int is1_up, is2_up;
  long unsigned int tmp_off = 0;
  int tmp_sgn, num1, num2, one = 1;
  double complex tmp_V;
  long unsigned int i_max;
  double complex dmv;

  i_max = X->Check.idim_max;
  X->Large.mode = M_CORR;

  for (i = 0; i < X->Def.NCisAjtCkuAlvDC; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    if (Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X, 2) != 0) {
      //error message will be added
      zclear(nstate, prod[i]);
      continue;
    }

    if (org_isite1 > X->Def.Nsite && org_isite3 > X->Def.Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        is1_up = X->Def.Tpow[org_isite1 - 1];
        is2_up = X->Def.Tpow[org_isite3 - 1];
        num1 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is1_up, org_sigma1);
        num2 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is2_up, org_sigma3);
        zaxpy_long(i_max*nstate, tmp_V * num1*num2, &vec[1][0], &Xvec[1][0]);
      }
      else if (org_isite1 == org_isite3 && org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) {
        is1_up = X->Def.Tpow[org_isite1 - 1];
        num1 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is1_up, org_sigma1);
        zaxpy_long(i_max*nstate, tmp_V * num1, &vec[1][0], &Xvec[1][0]);
      }
      else if (org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) {//exchange
        child_general_int_spin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else {  // other process is not allowed
                // error message will be added
      }
    }
    else if (org_isite1 > X->Def.Nsite || org_isite3 > X->Def.Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        is1_up = X->Def.Tpow[org_isite1 - 1];
        is2_up = X->Def.Tpow[org_isite3 - 1];
        num2 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is2_up, org_sigma3);
#pragma omp parallel for default(none)shared(vec,Xvec,nstate,one)      \
  firstprivate(i_max, tmp_V, is1_up, org_sigma1, X, num2) private(j, num1,dmv)
        for (j = 1; j <= i_max; j++) {
          num1 = child_Spin_CisAis(j, X, is1_up, org_sigma1);
          dmv = tmp_V * num1*num2;
          zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
        }
      }
      else if (org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) {//exchange
        child_general_int_spin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else {  // other process is not allowed
                // error message will be added
      }
    }
    else {
      isA_up = X->Def.Tpow[org_isite1 - 1];
      isB_up = X->Def.Tpow[org_isite3 - 1];
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp parallel for default(none) private(j) shared(vec,Xvec,nstate) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off, tmp_V)
        for (j = 1; j <= i_max; j++) {
          CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4,
            tmp_V, nstate, Xvec, vec, X);
        }
      }
      else if (org_isite1 == org_isite3 && org_sigma1 == org_sigma4 && org_sigma3 == org_sigma2) {
#pragma omp parallel for default(none) private(j, dmv) \
  firstprivate(i_max,X,isA_up,org_sigma1, tmp_V) shared(vec, list_1,Xvec,nstate,one)
        for (j = 1; j <= i_max; j++) {
          dmv = tmp_V * child_Spin_CisAis(j, X, isA_up, org_sigma1);
          zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
        }
      }
      else if (org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) { // exchange
#pragma omp parallel for default(none) private(j, tmp_sgn, dmv) shared(vec,Xvec,nstate,one) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V)
        for (j = 1; j <= i_max; j++) {
          tmp_sgn = child_exchange_spin_element(j, X, isA_up, isB_up, org_sigma2, org_sigma4, &tmp_off);
          dmv = tmp_sgn;
          zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[tmp_off][0], &one);
        }
      }
      else {  // other process is not allowed
                // error message will be added
      }
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief Child function to calculate two-body green's functions for General Spin model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
int expec_cisajscktalt_SpinGeneral(
  struct BindStruct *X,
  int nstate, 
  double complex **Xvec, 
  double complex **vec, 
  double complex **prod
){
  long unsigned int i, j;
  long unsigned int org_isite1, org_isite2, org_isite3, org_isite4;
  long unsigned int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long unsigned int tmp_off = 0;
  long unsigned int tmp_off_2 = 0;
  long unsigned int list1_off = 0;
  int num1, one = 1;
  double complex tmp_V;
  long unsigned int i_max;
  int tmp_Sz;
  long unsigned int tmp_org = 0;
  i_max = X->Check.idim_max;
  X->Large.mode = M_CORR;

  for (i = 0; i < X->Def.NCisAjtCkuAlvDC; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);

    if (Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X, 2) != 0) {
      zclear(nstate, prod[i]);
      continue;
    }
    tmp_Sz = 0;

    for (j = 0; j < 2; j++) {
      tmp_org = X->Def.CisAjtCkuAlvDC[i][4 * j + 1] * X->Def.Tpow[X->Def.CisAjtCkuAlvDC[i][4 * j]];
      tmp_Sz += GetLocal2Sz(X->Def.CisAjtCkuAlvDC[i][4 * j] + 1, tmp_org, X->Def.SiteToBit, X->Def.Tpow);
      tmp_org = X->Def.CisAjtCkuAlvDC[i][4 * j + 3] * X->Def.Tpow[X->Def.CisAjtCkuAlvDC[i][4 * j + 2]];
      tmp_Sz -= GetLocal2Sz(X->Def.CisAjtCkuAlvDC[i][4 * j + 2] + 1, tmp_org, X->Def.SiteToBit, X->Def.Tpow);
    }
    if (tmp_Sz != 0) { // not Sz conserved
      zclear(nstate, prod[i]);
      continue;
    }

    if (org_isite1 > X->Def.Nsite && org_isite3 > X->Def.Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        child_CisAisCjuAju_GeneralSpin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        child_CisAitCjuAjv_GeneralSpin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4,
          tmp_V, X, nstate, Xvec, vec);
      }
      else {
      }
    }
    else if (org_isite3 > X->Def.Nsite || org_isite1 > X->Def.Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        child_CisAisCjuAju_GeneralSpin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        child_CisAitCjuAjv_GeneralSpin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4,
          tmp_V, X, nstate, Xvec, vec);
      }
      else {
      }
    }
    else {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp parallel for default(none) private(j, num1) shared(vec,list_1,Xvec,nstate,one) \
firstprivate(i_max,X,org_isite1, org_sigma1,org_isite3, org_sigma3, tmp_V)
        for (j = 1; j <= i_max; j++) {
          num1 = BitCheckGeneral(list_1[j], org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
          if (num1 != FALSE) {
            num1 = BitCheckGeneral(list_1[j], org_isite3, org_sigma3, X->Def.SiteToBit, X->Def.Tpow);
            if (num1 != FALSE) {
              zaxpy_(&nstate, &tmp_V, &vec[j][0], &one, &Xvec[j][0], &one);
            }
          }
        }
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) private(j,num1) \
firstprivate(i_max,X,org_isite1,org_isite3,org_sigma1,org_sigma2,org_sigma3,org_sigma4,tmp_off,tmp_off_2,list1_off,myrank,tmp_V) \
  shared(vec,list_1,Xvec,nstate,one)
        for (j = 1; j <= i_max; j++) {
          num1 = GetOffCompGeneralSpin(list_1[j], org_isite3, org_sigma4, org_sigma3, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
          if (num1 != FALSE) {
            num1 = GetOffCompGeneralSpin(tmp_off, org_isite1, org_sigma2, org_sigma1, &tmp_off_2,
              X->Def.SiteToBit, X->Def.Tpow);
            if (num1 != FALSE) {
              ConvertToList1GeneralSpin(tmp_off_2, X->Check.sdim, &list1_off);
              zaxpy_(&nstate, &tmp_V, &vec[j][0], &one, &Xvec[list1_off][0], &one);
            }
          }
        }
      }
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief Child function to calculate two-body green's functions for 1/2 Spin GC model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
int expec_cisajscktalt_SpinGCHalf(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec, 
  double complex **prod
){
  long unsigned int i, j;
  long unsigned int org_isite1, org_isite2, org_isite3, org_isite4;
  long unsigned int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long unsigned int isA_up, isB_up;
  long unsigned int tmp_off = 0;
  double complex tmp_V;
  long unsigned int i_max;
  i_max = X->Check.idim_max;

  for (i = 0; i < X->Def.NCisAjtCkuAlvDC; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);

    if (Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X, 2) != 0) {
      //error message will be added
      zclear(nstate, prod[i]);
      continue;
    }

    if (org_isite1 > X->Def.Nsite && org_isite3 > X->Def.Nsite) { //org_isite3 >= org_isite1 > Nsite

      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        child_GC_CisAisCjuAju_spin_MPIdouble(
          org_isite1 - 1, org_sigma1, (org_isite3 - 1), org_sigma3, tmp_V, X, nstate, Xvec, vec);

      }
      else if (org_isite1 == org_isite3 && org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) { //diagonal (for spin: cuadcdau=cuau)
        child_GC_CisAis_spin_MPIdouble(
          org_isite1 - 1, org_sigma1, tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
        child_GC_CisAisCjuAjv_spin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, org_sigma4, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
        child_GC_CisAitCjuAju_spin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3,
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        child_GC_CisAitCiuAiv_spin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4, 
          tmp_V, X, nstate, Xvec, vec);
      }
    }
    else if (org_isite3 > X->Def.Nsite || org_isite1 > X->Def.Nsite) { //org_isite3 > Nsite >= org_isite1
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        child_GC_CisAisCjuAju_spin_MPIsingle(
          org_isite1 - 1, org_sigma1, (org_isite3 - 1), org_sigma3, tmp_V, X, nstate, Xvec, vec);

      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
        child_GC_CisAisCjuAjv_spin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, org_sigma4,
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
        child_GC_CisAitCjuAju_spin_MPIsingle(
          org_isite1 - 1, org_sigma2, org_isite3 - 1, org_sigma3, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        child_GC_CisAitCiuAiv_spin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4,
          tmp_V, X, nstate, Xvec, vec);
      }
    }
    else {
      if (org_isite1 == org_isite2 && org_isite3 == org_isite4) {
        isA_up = X->Def.Tpow[org_isite2 - 1];
        isB_up = X->Def.Tpow[org_isite4 - 1];
        if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp parallel for default(none) private(j) shared(vec,Xvec,nstate) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V)
          for (j = 1; j <= i_max; j++) {
            GC_CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4,
              tmp_V, nstate, Xvec, vec, X);
          }
        }
        else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) private(j) shared(vec,Xvec,nstate) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V)
          for (j = 1; j <= i_max; j++) {
            GC_CisAisCitAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up,
              tmp_V, nstate, Xvec, vec, X, &tmp_off);
          }
        }
        else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
#pragma omp parallel for default(none) private(j) shared(vec,Xvec,nstate) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V)
          for (j = 1; j <= i_max; j++) {
            GC_CisAitCiuAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up,
              tmp_V, nstate, Xvec, vec, X, &tmp_off);
          }
        }
        else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) private(j) shared(vec,Xvec,nstate) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V)
          for (j = 1; j <= i_max; j++) {
            GC_CisAitCiuAiv_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up,
              tmp_V, nstate, Xvec, vec, X, &tmp_off);
          }
        }
      }
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief Child function to calculate two-body green's functions for General Spin GC model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
int expec_cisajscktalt_SpinGCGeneral(
  struct BindStruct *X,
  int nstate, 
  double complex **Xvec, 
  double complex **vec, 
  double complex **prod
){
  long unsigned int i, j;
  long unsigned int org_isite1, org_isite2, org_isite3, org_isite4;
  long unsigned int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
  long unsigned int tmp_off = 0;
  long unsigned int tmp_off_2 = 0;
  int  num1, one = 1;
  double complex tmp_V;
  long unsigned int i_max;
  i_max = X->Check.idim_max;
  X->Large.mode = M_CORR;

  for(i=0;i<X->Def.NCisAjtCkuAlvDC;i++){
    zclear(i_max*nstate, &Xvec[1][0]);

    if (Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X, 2) != 0) {
      //error message will be added
      zclear(nstate, prod[i]);
      continue;
    }

    if (org_isite1 > X->Def.Nsite && org_isite3 > X->Def.Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        child_GC_CisAisCjuAju_GeneralSpin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
        child_GC_CisAisCjuAjv_GeneralSpin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, org_sigma4, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
        child_GC_CisAitCjuAju_GeneralSpin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        child_GC_CisAitCjuAjv_GeneralSpin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4,
          tmp_V, X, nstate, Xvec, vec);
      }
    }
    else if (org_isite3 > X->Def.Nsite || org_isite1 > X->Def.Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        child_GC_CisAisCjuAju_GeneralSpin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
        child_GC_CisAisCjuAjv_GeneralSpin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, org_sigma4, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
        child_GC_CisAitCjuAju_GeneralSpin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        child_GC_CisAitCjuAjv_GeneralSpin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4, 
          tmp_V, X, nstate, Xvec, vec);
      }
    }
    else {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp parallel for default(none) private(j, num1) shared(vec,Xvec,nstate,one) \
firstprivate(i_max,X,org_isite1, org_sigma1,org_isite3, org_sigma3, tmp_V)
        for (j = 1; j <= i_max; j++) {
          num1 = BitCheckGeneral(j - 1, org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
          if (num1 != FALSE) {
            num1 = BitCheckGeneral(j - 1, org_isite3, org_sigma3, X->Def.SiteToBit, X->Def.Tpow);
            if (num1 != FALSE) {
              zaxpy_(&nstate, &tmp_V, &vec[j][0], &one, &Xvec[j][0], &one);
            }
          }
        }
      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) private(j, num1) shared(vec,Xvec,nstate,one) \
firstprivate(i_max,X, org_isite1, org_isite3, org_sigma1,org_sigma3,org_sigma4, tmp_off, tmp_V)
        for (j = 1; j <= i_max; j++) {
          num1 = GetOffCompGeneralSpin(j - 1, org_isite3, org_sigma4, org_sigma3,
            &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
          if (num1 != FALSE) {
            num1 = BitCheckGeneral(tmp_off, org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
            if (num1 != FALSE) {
              zaxpy_(&nstate, &tmp_V, &vec[j][0], &one, &Xvec[tmp_off + 1][0], &one);
            }
          }
        }
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
#pragma omp parallel for default(none) private(j, num1) shared(vec,Xvec,nstate,one) \
firstprivate(i_max,X, org_isite1, org_isite3, org_sigma1,org_sigma2, org_sigma3, tmp_off, tmp_V)
        for (j = 1; j <= i_max; j++) {
          num1 = BitCheckGeneral(j - 1, org_isite3, org_sigma3, X->Def.SiteToBit, X->Def.Tpow);
          if (num1 != FALSE) {
            num1 = GetOffCompGeneralSpin(j - 1, org_isite1, org_sigma2, org_sigma1,
              &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
            if (num1 != FALSE) {
              zaxpy_(&nstate, &tmp_V, &vec[j][0], &one, &Xvec[tmp_off + 1][0], &one);
            }
          }
        }
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) private(j, num1) \
firstprivate(i_max,X,org_isite1,org_isite3,org_sigma1,org_sigma2,org_sigma3,org_sigma4,tmp_off,tmp_off_2,tmp_V) \
  shared(vec,Xvec,nstate,one)
        for (j = 1; j <= i_max; j++) {
          num1 = GetOffCompGeneralSpin(j - 1, org_isite3, org_sigma4, org_sigma3,
            &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
          if (num1 != FALSE) {
            num1 = GetOffCompGeneralSpin(tmp_off, org_isite1, org_sigma2, org_sigma1, 
              &tmp_off_2, X->Def.SiteToBit, X->Def.Tpow);
            if (num1 != FALSE) {
              zaxpy_(&nstate, &tmp_V, &vec[j][0], &one, &Xvec[tmp_off_2 + 1][0], &one);
            }
          }
        }
      }
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief Child function to calculate six-body green's functions for 1/2 Spin GC model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
int expec_Sixbody_SpinGCHalf(
  struct BindStruct *X, 
  int nstate,
  double complex** Xvec,
  double complex**vec,
  double complex**prod
){
    long unsigned int i,j;
    long unsigned int tmp_org_isite1,tmp_org_isite2,tmp_org_isite3,tmp_org_isite4,tmp_org_isite5,tmp_org_isite6,tmp_org_isite7,tmp_org_isite8,tmp_org_isite9,tmp_org_isite10,tmp_org_isite11,tmp_org_isite12;
    long unsigned int tmp_org_sigma1,tmp_org_sigma2,tmp_org_sigma3,tmp_org_sigma4,tmp_org_sigma5,tmp_org_sigma6,tmp_org_sigma7,tmp_org_sigma8,tmp_org_sigma9,tmp_org_sigma10,tmp_org_sigma11,tmp_org_sigma12;
    long unsigned int org_isite1,org_isite2,org_isite3,org_isite4,org_isite5,org_isite6,org_isite7,org_isite8,org_isite9,org_isite10,org_isite11,org_isite12;
    long unsigned int org_sigma1,org_sigma2,org_sigma3,org_sigma4,org_sigma5,org_sigma6,org_sigma7,org_sigma8,org_sigma9,org_sigma10,org_sigma11,org_sigma12;
    long unsigned int isA_up, isB_up;
    long unsigned int tmp_off=0;
    double complex tmp_V;
    double complex dam_pr;
    double complex **vec_pr;

    long int i_max;
    i_max=X->Check.idim_max;
    vec_pr = cd_2d_allocate(i_max + 1, nstate);

    for(i=0;i<X->Def.NSBody;i++){
        //printf("%d %d \n",i,X->Def.NSBody);

        tmp_org_isite1   = X->Def.SBody[i][0]+1;
        tmp_org_sigma1   = X->Def.SBody[i][1];
        tmp_org_isite2   = X->Def.SBody[i][2]+1;
        tmp_org_sigma2   = X->Def.SBody[i][3];
        /**/
        tmp_org_isite3   = X->Def.SBody[i][4]+1;
        tmp_org_sigma3   = X->Def.SBody[i][5];
        tmp_org_isite4   = X->Def.SBody[i][6]+1;
        tmp_org_sigma4   = X->Def.SBody[i][7];
        /**/
        tmp_org_isite5   = X->Def.SBody[i][8]+1;
        tmp_org_sigma5   = X->Def.SBody[i][9];
        tmp_org_isite6   = X->Def.SBody[i][10]+1;
        tmp_org_sigma6   = X->Def.SBody[i][11];
        /**/
        tmp_org_isite7   = X->Def.SBody[i][12]+1;
        tmp_org_sigma7   = X->Def.SBody[i][13];
        tmp_org_isite8   = X->Def.SBody[i][14]+1;
        tmp_org_sigma8   = X->Def.SBody[i][15];
        /**/
        tmp_org_isite9   = X->Def.SBody[i][16]+1;
        tmp_org_sigma9   = X->Def.SBody[i][17];
        tmp_org_isite10  = X->Def.SBody[i][18]+1;
        tmp_org_sigma10  = X->Def.SBody[i][19];
        /**/
        tmp_org_isite11  = X->Def.SBody[i][20]+1;
        tmp_org_sigma11  = X->Def.SBody[i][21];
        tmp_org_isite12  = X->Def.SBody[i][22]+1;
        tmp_org_sigma12  = X->Def.SBody[i][23];
 
        /**/
        org_isite5   = X->Def.SBody[i][8]+1;
        org_sigma5   = X->Def.SBody[i][9];
        org_isite6   = X->Def.SBody[i][10]+1;
        org_sigma6   = X->Def.SBody[i][11];
        /**/
        org_isite7   = X->Def.SBody[i][12]+1;
        org_sigma7   = X->Def.SBody[i][13];
        org_isite8   = X->Def.SBody[i][14]+1;
        org_sigma8   = X->Def.SBody[i][15];
        /**/
        org_isite9   = X->Def.SBody[i][16]+1;
        org_sigma9   = X->Def.SBody[i][17];
        org_isite10  = X->Def.SBody[i][18]+1;
        org_sigma10  = X->Def.SBody[i][19];
        /**/
        org_isite11  = X->Def.SBody[i][20]+1;
        org_sigma11  = X->Def.SBody[i][21];
        org_isite12  = X->Def.SBody[i][22]+1;
        org_sigma12  = X->Def.SBody[i][23];

        X->Large.mode = M_MLTPLY;
        /* |vec_pr_0>= c11a12|vec>*/
        zclear((i_max + 1) * nstate, &Xvec[0][0]);
        zclear((i_max + 1) * nstate, &vec_pr[0][0]);
        mltplyHalfSpinGC_mini(X,tmp_org_isite11-1,tmp_org_sigma11,tmp_org_isite12-1,tmp_org_sigma12,nstate, Xvec,vec);
        /* |vec_pr_1>= c9a10|vec_pr_0>*/
        mltplyHalfSpinGC_mini(X,tmp_org_isite9-1,tmp_org_sigma9,tmp_org_isite10-1,tmp_org_sigma10,nstate, vec_pr, Xvec);
        zclear((i_max + 1) * nstate, &Xvec[0][0]);
        /* |vec_pr_2>= c7a8|vec_pr_1>*/
        mltplyHalfSpinGC_mini(X,tmp_org_isite7-1,tmp_org_sigma7,tmp_org_isite8-1,tmp_org_sigma8, nstate, Xvec, vec_pr);
        zclear((i_max + 1) * nstate, &vec_pr[0][0]);
        /* |vec_pr>= c5a6|vec_pr_2>*/
        mltplyHalfSpinGC_mini(X,tmp_org_isite5-1,tmp_org_sigma5,tmp_org_isite6-1,tmp_org_sigma6, nstate,vec_pr,Xvec);
        zclear((i_max + 1) * nstate, &Xvec[0][0]);
        X->Large.mode = H_CORR;

        if(Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X,6)!=0){
            //error message will be added
            zclear(nstate, prod[i]);
            continue;
        }
        /*
        printf("check: %d %d %d %d %d %d %d %d %d %d %d %d \n",
        org_isite1-1,org_sigma1 ,org_isite2-1,org_sigma2,
        org_isite3-1,org_sigma3 ,org_isite4-1,org_sigma4,
        org_isite5-1,org_sigma5 ,org_isite6-1,org_sigma6
        );
        */
       
        dam_pr=0.0;
        if(org_isite1>X->Def.Nsite && org_isite3>X->Def.Nsite){ //org_isite3 >= org_isite1 > Nsite
           //printf("D-MPI \n");
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr += child_GC_CisAisCjuAju_spin_MPIdouble( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_isite1 ==org_isite3 && org_sigma1 ==org_sigma4 && org_sigma2 ==org_sigma3){ //diagonal (for spin: cuadcdau=cuau)
                dam_pr += child_GC_CisAis_spin_MPIdouble((org_isite1-1), org_sigma1, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIdouble(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, nstate, Xvec, vec_pr);
            }
        }
        else if(org_isite3>X->Def.Nsite || org_isite1>X->Def.Nsite){ //org_isite3 > Nsite >= org_isite1
           //printf("S-MPI \n");
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr += child_GC_CisAisCjuAju_spin_MPIsingle( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIsingle(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIsingle(org_isite1-1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, nstate, Xvec, vec_pr);
            }
        }
        else{
            if(org_isite1==org_isite2 && org_isite3==org_isite4){
                isA_up = X->Def.Tpow[org_isite2-1];
                isB_up = X->Def.Tpow[org_isite4-1];
                if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j,i) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V,nstate) shared(Xvec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr +=GC_CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4, tmp_V, nstate, Xvec, vec_pr, X);
                    }
                }else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V,nstate) shared(Xvec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAisCitAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, nstate, Xvec, vec_pr, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V,nstate) shared(Xvec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, nstate, Xvec, vec_pr, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V,nstate) shared(Xvec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiv_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, nstate, Xvec, vec_pr, X, &tmp_off);
                    }
                }
            }
        }
        MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
    }
    return 0;
}
/**
 * @brief Parent function to calculate two-body green's functions for Spin model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
int expec_cisajscktalt_Spin(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec,
  double complex **prod
) {
  int info = 0;
  if (X->Def.iFlgGeneralSpin == FALSE) {
    info = expec_cisajscktalt_SpinHalf(X, nstate, Xvec, vec, prod);
  }
  else {
    info = expec_cisajscktalt_SpinGeneral(X, nstate, Xvec, vec, prod);
  }
  return info;
}
/**
 * @brief Child function to calculate four-body green's functions for 1/2 Spin GC model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
int expec_Fourbody_SpinGCHalf(
  struct BindStruct *X, 
  int nstate,
  double complex** Xvec,
  double complex**vec,
  double complex**prod
){
    long unsigned int i,j;
    long unsigned int tmp_org_isite1,tmp_org_isite2,tmp_org_isite3,tmp_org_isite4,tmp_org_isite5,tmp_org_isite6,tmp_org_isite7,tmp_org_isite8;
    long unsigned int tmp_org_sigma1,tmp_org_sigma2,tmp_org_sigma3,tmp_org_sigma4,tmp_org_sigma5,tmp_org_sigma6,tmp_org_sigma7,tmp_org_sigma8;
    long unsigned int org_isite1,org_isite2,org_isite3,org_isite4,org_isite5,org_isite6,org_isite7,org_isite8;
    long unsigned int org_sigma1,org_sigma2,org_sigma3,org_sigma4,org_sigma5,org_sigma6,org_sigma7,org_sigma8;
    long unsigned int isA_up, isB_up;
    long unsigned int tmp_off=0;
    double complex tmp_V;
    double complex dam_pr;
    double complex **vec_pr;

    long int i_max;
    i_max=X->Check.idim_max;
    vec_pr = cd_2d_allocate(i_max + 1, nstate);

    for(i=0;i<X->Def.NFBody;i++){

        tmp_org_isite1   = X->Def.FBody[i][0]+1;
        tmp_org_sigma1   = X->Def.FBody[i][1];
        tmp_org_isite2   = X->Def.FBody[i][2]+1;
        tmp_org_sigma2   = X->Def.FBody[i][3];
        /**/
        tmp_org_isite3   = X->Def.FBody[i][4]+1;
        tmp_org_sigma3   = X->Def.FBody[i][5];
        tmp_org_isite4   = X->Def.FBody[i][6]+1;
        tmp_org_sigma4   = X->Def.FBody[i][7];
        /**/
        tmp_org_isite5   = X->Def.FBody[i][8]+1;
        tmp_org_sigma5   = X->Def.FBody[i][9];
        tmp_org_isite6   = X->Def.FBody[i][10]+1;
        tmp_org_sigma6   = X->Def.FBody[i][11];
        /**/
        tmp_org_isite7   = X->Def.FBody[i][12]+1;
        tmp_org_sigma7   = X->Def.FBody[i][13];
        tmp_org_isite8   = X->Def.FBody[i][14]+1;
        tmp_org_sigma8   = X->Def.FBody[i][15];
        /**/
        org_isite5   = X->Def.FBody[i][8]+1;
        org_sigma5   = X->Def.FBody[i][9];
        org_isite6   = X->Def.FBody[i][10]+1;
        org_sigma6   = X->Def.FBody[i][11];
        /**/
        org_isite7   = X->Def.FBody[i][12]+1;
        org_sigma7   = X->Def.FBody[i][13];
        org_isite8   = X->Def.FBody[i][14]+1;
        org_sigma8   = X->Def.FBody[i][15];

        X->Large.mode = M_MLTPLY;
        /* |vec_pr_tmp>= c7a8|vec>*/
        zclear((i_max + 1) * nstate, &Xvec[0][0]);
        zclear((i_max + 1) * nstate, &vec_pr[0][0]);
        mltplyHalfSpinGC_mini(X,tmp_org_isite7-1,tmp_org_sigma7,tmp_org_isite8-1,tmp_org_sigma8,nstate,Xvec,vec);
        /* |vec_pr>= c5a6|vec_pr_tmp>*/
        mltplyHalfSpinGC_mini(X,tmp_org_isite5-1,tmp_org_sigma5,tmp_org_isite6-1,tmp_org_sigma6,nstate,vec_pr, Xvec);
        zclear((i_max + 1) * nstate, &Xvec[0][0]);
        X->Large.mode = H_CORR;

        if(Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X,4)!=0){
          zclear(nstate, prod[i]);
          continue;
        }
        /*
        printf("check: %d %d %d %d %d %d %d %d %d %d %d %d \n",
        org_isite1-1,org_sigma1 ,org_isite2-1,org_sigma2,
        org_isite3-1,org_sigma3 ,org_isite4-1,org_sigma4,
        org_isite5-1,org_sigma5 ,org_isite6-1,org_sigma6
        );
        */
       
        dam_pr=0.0;
        if(org_isite1>X->Def.Nsite && org_isite3>X->Def.Nsite){ //org_isite3 >= org_isite1 > Nsite
           //printf("D-MPI \n");
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr += child_GC_CisAisCjuAju_spin_MPIdouble( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_isite1 ==org_isite3 && org_sigma1 ==org_sigma4 && org_sigma2 ==org_sigma3){ //diagonal (for spin: cuadcdau=cuau)
                dam_pr += child_GC_CisAis_spin_MPIdouble((org_isite1-1), org_sigma1, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIdouble(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, nstate, Xvec, vec_pr);
            }
        }
        else if(org_isite3>X->Def.Nsite || org_isite1>X->Def.Nsite){ //org_isite3 > Nsite >= org_isite1
           //printf("S-MPI \n");
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr += child_GC_CisAisCjuAju_spin_MPIsingle( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIsingle(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIsingle(org_isite1-1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, nstate, Xvec, vec_pr);
            }
        }
        else{
            if(org_isite1==org_isite2 && org_isite3==org_isite4){
                isA_up = X->Def.Tpow[org_isite2-1];
                isB_up = X->Def.Tpow[org_isite4-1];
                if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j,i) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V,nstate) shared(Xvec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr +=GC_CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4, tmp_V, nstate, Xvec, vec_pr, X);
                    }
                }else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V,nstate) shared(Xvec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAisCitAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, nstate, Xvec, vec_pr, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V,nstate) shared(Xvec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, nstate, Xvec, vec_pr, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V,nstate) shared(Xvec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiv_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, nstate, Xvec, vec_pr, X, &tmp_off);
                    }
                }
            }
        }
        MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
    }
    free_cd_2d_allocate(vec_pr);
    return 0;
}




/**
 * @brief Child function to calculate three-body green's functions for 1/2 Spin GC model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
int expec_Threebody_SpinGCHalf(
  struct BindStruct *X, 
  int nstate, 
  double complex **Xvec,
  double complex **vec, 
  double complex **prod
){
    long unsigned int i,j;
    long unsigned int tmp_org_isite1,tmp_org_isite2,tmp_org_isite3,tmp_org_isite4,tmp_org_isite5,tmp_org_isite6;
    long unsigned int tmp_org_sigma1,tmp_org_sigma2,tmp_org_sigma3,tmp_org_sigma4,tmp_org_sigma5,tmp_org_sigma6;
    long unsigned int org_isite1,org_isite2,org_isite3,org_isite4,org_isite5,org_isite6;
    long unsigned int org_sigma1,org_sigma2,org_sigma3,org_sigma4,org_sigma5,org_sigma6;
    long unsigned int isA_up, isB_up;
    long unsigned int tmp_off=0;
    double complex tmp_V;
    double complex dam_pr;
    double complex **vec_pr;

    long int i_max;
    i_max=X->Check.idim_max;

    vec_pr = cd_2d_allocate(i_max + 1, nstate);
    for(i=0;i<X->Def.NTBody;i++){
        tmp_org_isite1   = X->Def.TBody[i][0]+1;
        tmp_org_sigma1   = X->Def.TBody[i][1];
        tmp_org_isite2   = X->Def.TBody[i][2]+1;
        tmp_org_sigma2   = X->Def.TBody[i][3];
        /**/
        tmp_org_isite3   = X->Def.TBody[i][4]+1;
        tmp_org_sigma3   = X->Def.TBody[i][5];
        tmp_org_isite4   = X->Def.TBody[i][6]+1;
        tmp_org_sigma4   = X->Def.TBody[i][7];
        /**/
        tmp_org_isite5   = X->Def.TBody[i][8]+1;
        tmp_org_sigma5   = X->Def.TBody[i][9];
        tmp_org_isite6   = X->Def.TBody[i][10]+1;
        tmp_org_sigma6   = X->Def.TBody[i][11];
        /**/
        org_isite5   = X->Def.TBody[i][8]+1;
        org_sigma5   = X->Def.TBody[i][9];
        org_isite6   = X->Def.TBody[i][10]+1;
        org_sigma6   = X->Def.TBody[i][11];

        X->Large.mode = M_MLTPLY;
        /* |vec_pr>= c5a6|phi>*/
        zclear((i_max + 1) * nstate, &Xvec[0][0]);
        zclear((i_max + 1) * nstate, &vec_pr[0][0]);
        mltplyHalfSpinGC_mini(X,tmp_org_isite5-1,tmp_org_sigma5,tmp_org_isite6-1,tmp_org_sigma6,nstate,vec_pr,vec);
        X->Large.mode = H_CORR;

        if(Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X,3)!=0){
            //error message will be added
          zclear(nstate, prod[i]);
          continue;
        }
        /*
        printf("check: %d %d %d %d %d %d %d %d %d %d %d %d \n",
        org_isite1-1,org_sigma1 ,org_isite2-1,org_sigma2,
        org_isite3-1,org_sigma3 ,org_isite4-1,org_sigma4,
        org_isite5-1,org_sigma5 ,org_isite6-1,org_sigma6
        );
        */
       
        dam_pr=0.0;
        if(org_isite1>X->Def.Nsite && org_isite3>X->Def.Nsite){ //org_isite3 >= org_isite1 > Nsite
           //printf("D-MPI \n");
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr += child_GC_CisAisCjuAju_spin_MPIdouble( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_isite1 ==org_isite3 && org_sigma1 ==org_sigma4 && org_sigma2 ==org_sigma3){ //diagonal (for spin: cuadcdau=cuau)
                dam_pr += child_GC_CisAis_spin_MPIdouble((org_isite1-1), org_sigma1, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIdouble(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, nstate, Xvec, vec_pr);
            }
        }
        else if(org_isite3>X->Def.Nsite || org_isite1>X->Def.Nsite){ //org_isite3 > Nsite >= org_isite1
           //printf("S-MPI \n");
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr += child_GC_CisAisCjuAju_spin_MPIsingle( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIsingle(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIsingle(org_isite1-1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, nstate, Xvec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, nstate, Xvec, vec_pr);
            }
        }
        else{
            if(org_isite1==org_isite2 && org_isite3==org_isite4){
                isA_up = X->Def.Tpow[org_isite2-1];
                isB_up = X->Def.Tpow[org_isite4-1];
                if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j,i) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V,nstate) shared(Xvec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr +=GC_CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4, tmp_V, nstate, Xvec, vec_pr, X);
                    }
                }else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V,nstate) shared(Xvec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAisCitAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, nstate, Xvec, vec_pr, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V,nstate) shared(Xvec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, nstate, Xvec, vec_pr, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V,nstate) shared(Xvec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiv_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, nstate, Xvec, vec_pr, X, &tmp_off);
                    }
                }
            }
        }
        MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
    }
    free_cd_2d_allocate(vec_pr);
    return 0;
}


/**
 * @brief Parent function to calculate two-body green's functions for Spin GC model
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 * @param _fp [in] output file name
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 *
 */
int expec_cisajscktalt_SpinGC(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec,
  double complex **prod,
  double complex** prod_2,
  double complex** prod_3,
  double complex** prod_4
) {
  int info = 0;
  if (X->Def.iFlgGeneralSpin == FALSE) {
    info = expec_cisajscktalt_SpinGCHalf(X, nstate, Xvec, vec, prod);
    if (X->Def.NTBody > 0) {
      info = expec_Threebody_SpinGCHalf(X, nstate,Xvec, vec, prod_2);
    }
    if (X->Def.NFBody > 0) {
      info = expec_Fourbody_SpinGCHalf(X, nstate, Xvec, vec, prod_3);
    }
    if (X->Def.NSBody > 0) {
      info = expec_Sixbody_SpinGCHalf(X, nstate, Xvec, vec, prod_4);
    }
  }
  else {
    info = expec_cisajscktalt_SpinGCGeneral(X, nstate, Xvec, vec, prod);
  }
  return info;
}
/**
 * @brief Parent function to calculate two-body green's functions
 *
 * @param X [in] data list for calculation
 * @param vec [in] eigenvectors
 *
 * @retval 0 normally finished
 * @retval -1 abnormally finished
 * @note The origin of function's name cisajscktalt comes from c=creation, i=ith site, s=spin, a=annihiration, j=jth site and so on.
 *
 * @version 0.2
 * @details add function to treat the case of general spin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int expec_cisajscktaltdc
(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec
) {
  FILE *fp;
  char sdt[D_FileNameMax], sdt_2[D_FileNameMax], sdt_3[D_FileNameMax], sdt_4[D_FileNameMax], * tmp_char;
  long unsigned int irght, ilft, ihfbit, icaca;
  double complex **prod, ** prod_2, ** prod_3, ** prod_4;
  //For TPQ
  int step = 0, rand_i = 0, istate;

  if (X->Def.NCisAjtCkuAlvDC < 1) return 0;
  X->Large.mode = M_CORR;

  if (GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit) != 0) {
    return -1;
  }

  //Make File Name for output
  prod = cd_2d_allocate(X->Def.NCisAjtCkuAlvDC, nstate);
  prod_2 = cd_2d_allocate(X->Def.NTBody, nstate);
  prod_3 = cd_2d_allocate(X->Def.NFBody, nstate);
  prod_4 = cd_2d_allocate(X->Def.NSBody, nstate);
  switch (X->Def.iCalcType) {
  case TPQCalc:
  case cTPQ:
    step = X->Def.istep;
    TimeKeeperWithRandAndStep(X, cFileNameTimeKeep, cTPQExpecTwoBodyGStart, "a", 0, step);
    break;
  case TimeEvolution:
    step = X->Def.istep;
    TimeKeeperWithStep(X, cFileNameTimeKeep, cTEExpecTwoBodyGStart, "a", step);
    break;
  case FullDiag:
  case CG:
    break;
  }

  switch (X->Def.iCalcModel) {
  case HubbardGC:
    if (expec_cisajscktalt_HubbardGC(X, nstate, Xvec, vec, prod) != 0) {
      return -1;
    }
    break;

  case KondoGC:
  case Hubbard:
  case Kondo:
    if (expec_cisajscktalt_Hubbard(X, nstate, Xvec, vec, prod) != 0) {
      return -1;
    }
    break;

  case Spin:
    if (expec_cisajscktalt_Spin(X, nstate, Xvec, vec, prod) != 0) {
      return -1;
    }
    break;

  case SpinGC:
    if (expec_cisajscktalt_SpinGC(X, nstate, Xvec, vec, prod, prod_2, prod_3, prod_4) != 0) {
      return -1;
    }
    break;

  default:
    return -1;
  }

  for (istate = 0; istate < nstate; istate++) {
    switch (X->Def.iCalcType) {
    case Lanczos:
      if (X->Def.St == 0) {
        sprintf(sdt, cFileName2BGreen_Lanczos, X->Def.CDataFileHead);
        sprintf(sdt_2, cFileName3BGreen_Lanczos, X->Def.CDataFileHead);
        sprintf(sdt_3, cFileName4BGreen_Lanczos, X->Def.CDataFileHead);
        sprintf(sdt_4, cFileName6BGreen_Lanczos, X->Def.CDataFileHead);
        TimeKeeper(X, cFileNameTimeKeep, cLanczosExpecTwoBodyGStart, "a");
        fprintf(stdoutMPI, "%s", cLogLanczosExpecTwoBodyGStart);
      }
      else if (X->Def.St == 1) {
        sprintf(sdt, cFileName2BGreen_CG, X->Def.CDataFileHead);
        sprintf(sdt_2, cFileName3BGreen_Lanczos, X->Def.CDataFileHead);
        sprintf(sdt_3, cFileName4BGreen_Lanczos, X->Def.CDataFileHead);
        sprintf(sdt_4, cFileName6BGreen_Lanczos, X->Def.CDataFileHead);
        TimeKeeper(X, cFileNameTimeKeep, cCGExpecTwoBodyGStart, "a");
        fprintf(stdoutMPI, "%s", cLogLanczosExpecTwoBodyGStart);
      }
      break;
    case TPQCalc:
    case cTPQ:
      step = X->Def.istep;
      sprintf(sdt, cFileName2BGreen_TPQ, X->Def.CDataFileHead, istate, step);
      sprintf(sdt_2, cFileName3BGreen_TPQ, X->Def.CDataFileHead, istate, step);
      sprintf(sdt_3, cFileName4BGreen_TPQ, X->Def.CDataFileHead, istate, step);
      sprintf(sdt_4, cFileName6BGreen_TPQ, X->Def.CDataFileHead, istate, step);
      break; 
    case TimeEvolution:
      step = X->Def.istep;
      sprintf(sdt, cFileName2BGreen_TE, X->Def.CDataFileHead, step);
      sprintf(sdt_2, cFileName3BGreen_TE, X->Def.CDataFileHead, step);
      sprintf(sdt_3, cFileName4BGreen_TE, X->Def.CDataFileHead, step);
      sprintf(sdt_4, cFileName6BGreen_TE, X->Def.CDataFileHead, step);
      break;
    case FullDiag:
    case CG:
      sprintf(sdt, cFileName2BGreen_FullDiag, X->Def.CDataFileHead, istate);
      sprintf(sdt, cFileName2BGreen_FullDiag, X->Def.CDataFileHead, istate);
      sprintf(sdt_2, cFileName3BGreen_FullDiag, X->Def.CDataFileHead, istate);
      sprintf(sdt_3, cFileName4BGreen_FullDiag, X->Def.CDataFileHead, istate);
      sprintf(sdt_4, cFileName6BGreen_FullDiag, X->Def.CDataFileHead, istate);
      break;
    }
    if (childfopenMPI(sdt, "w", &fp) == 0) {
      for (icaca = 0; icaca < X->Def.NCisAjtCkuAlvDC; icaca++) {
        fprintf(fp, " %4d %4d %4d %4d %4d %4d %4d %4d %.10lf %.10lf\n",
          X->Def.CisAjtCkuAlvDC[icaca][0], X->Def.CisAjtCkuAlvDC[icaca][1],
          X->Def.CisAjtCkuAlvDC[icaca][2], X->Def.CisAjtCkuAlvDC[icaca][3],
          X->Def.CisAjtCkuAlvDC[icaca][4], X->Def.CisAjtCkuAlvDC[icaca][5],
          X->Def.CisAjtCkuAlvDC[icaca][6], X->Def.CisAjtCkuAlvDC[icaca][7],
          creal(prod[icaca][istate]), cimag(prod[icaca][istate]));
      }
      fclose(fp);
    }
    if (X->Def.NTBody > 0) {
      if (childfopenMPI(sdt_2, "w", &fp) == 0) {
        for (icaca = 0; icaca < X->Def.NTBody; icaca++) {
          fprintf(fp, " %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",
            X->Def.TBody[icaca][0], X->Def.TBody[icaca][1], X->Def.TBody[icaca][2], X->Def.TBody[icaca][3],
            X->Def.TBody[icaca][4], X->Def.TBody[icaca][5], X->Def.TBody[icaca][6], X->Def.TBody[icaca][7],
            X->Def.TBody[icaca][8], X->Def.TBody[icaca][9], X->Def.TBody[icaca][10], X->Def.TBody[icaca][11],
            creal(prod_2[icaca][istate]), cimag(prod_2[icaca][istate]));
        }
        fclose(fp);
      }
    }
    if (X->Def.NFBody > 0) {
      if (childfopenMPI(sdt_3, "w", &fp) == 0) {
        for (icaca = 0; icaca < X->Def.NFBody; icaca++) {
          fprintf(fp, " %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld  %4ld %4ld %4ld %4ld %.10lf %.10lf \n",
            X->Def.FBody[icaca][0], X->Def.FBody[icaca][1], X->Def.FBody[icaca][2], X->Def.FBody[icaca][3],
            X->Def.FBody[icaca][4], X->Def.FBody[icaca][5], X->Def.FBody[icaca][6], X->Def.FBody[icaca][7],
            X->Def.FBody[icaca][8], X->Def.FBody[icaca][9], X->Def.FBody[icaca][10], X->Def.FBody[icaca][11],
            X->Def.FBody[icaca][12], X->Def.FBody[icaca][13], X->Def.FBody[icaca][14], X->Def.FBody[icaca][15],
            creal(prod_3[icaca][istate]), cimag(prod_3[icaca][istate]));
        }
        fclose(fp);
      }
    }
    if (X->Def.NSBody > 0) {
      if (childfopenMPI(sdt_4, "w", &fp) == 0) {
        for (icaca = 0; icaca < X->Def.NSBody; icaca++) {
          fprintf(fp, " %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld  %4ld %4ld %4ld %4ld  %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",
            X->Def.SBody[icaca][0], X->Def.SBody[icaca][1], X->Def.SBody[icaca][2], X->Def.SBody[icaca][3],
            X->Def.SBody[icaca][4], X->Def.SBody[icaca][5], X->Def.SBody[icaca][6], X->Def.SBody[icaca][7],
            X->Def.SBody[icaca][8], X->Def.SBody[icaca][9], X->Def.SBody[icaca][10], X->Def.SBody[icaca][11],
            X->Def.SBody[icaca][12], X->Def.SBody[icaca][13], X->Def.SBody[icaca][14], X->Def.SBody[icaca][15],
            X->Def.SBody[icaca][16], X->Def.SBody[icaca][17], X->Def.SBody[icaca][18], X->Def.SBody[icaca][19],
            X->Def.SBody[icaca][20], X->Def.SBody[icaca][21], X->Def.SBody[icaca][22], X->Def.SBody[icaca][23],
            creal(prod_4[icaca][istate]), cimag(prod_4[icaca][istate]));
        }
        fclose(fp);
      }
    }
  }/*for (istate = 0; istate < nstate; istate++)*/

  if (X->Def.iCalcType == TPQCalc || X->Def.iCalcType == cTPQ) {
    TimeKeeperWithRandAndStep(X, cFileNameTimeKeep, cTPQExpecTwoBodyGFinish, "a", rand_i, step);
  }
  else if (X->Def.iCalcType == TimeEvolution) {
    TimeKeeperWithStep(X, cFileNameTimeKeep, cTEExpecTwoBodyGFinish, "a", step);
  }
  //[s] this part will be added
  /* For FullDiag, it is convinient to calculate the total spin for each vector.
     Such functions will be added
     if(X->Def.iCalcType==FullDiag){
     if(X->Def.iCalcModel==Spin){
     expec_cisajscktaltdc_alldiag_spin(X,vec);
     }else if(X->Def.iCalcModel==Hubbard || X->Def.iCalcModel==Kondo){
     expec_cisajscktaltdc_alldiag(X,vec);
     }else{//
     X->Phys.s2=0.0;
     }
     }
  */
  //[e]
  free_cd_2d_allocate(prod);
  free_cd_2d_allocate(prod_2);
  free_cd_2d_allocate(prod_3);
  free_cd_2d_allocate(prod_4);
  return 0;
}
