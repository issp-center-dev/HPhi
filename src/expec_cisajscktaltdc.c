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
                         struct BindStruct *X
                         )
{
  long unsigned int tmp_org_isite1,tmp_org_isite2,tmp_org_isite3,tmp_org_isite4;
  long unsigned int tmp_org_sigma1,tmp_org_sigma2,tmp_org_sigma3,tmp_org_sigma4;

  tmp_org_isite1   = X->Def.CisAjtCkuAlvDC[i][0]+1;
  tmp_org_sigma1   = X->Def.CisAjtCkuAlvDC[i][1];
  tmp_org_isite2   = X->Def.CisAjtCkuAlvDC[i][2]+1;
  tmp_org_sigma2   = X->Def.CisAjtCkuAlvDC[i][3];
  tmp_org_isite3   = X->Def.CisAjtCkuAlvDC[i][4]+1;
  tmp_org_sigma3   = X->Def.CisAjtCkuAlvDC[i][5];
  tmp_org_isite4   = X->Def.CisAjtCkuAlvDC[i][6]+1;
  tmp_org_sigma4   = X->Def.CisAjtCkuAlvDC[i][7];

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
  long int i_max;

  for (i = 0; i < X->Def.NCisAjtCkuAlvDC; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
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

        X_GC_child_CisAisCjtAjt_Hubbard_MPI(org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3,
          1.0, X, nstate, Xvec, vec);
      }
      else if (isite1 == isite2 && isite3 != isite4) {

        X_GC_child_CisAisCjtAku_Hubbard_MPI(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, org_isite4 - 1, org_sigma4,
          1.0, X, nstate, Xvec, vec);

      }
      else if (isite1 != isite2 && isite3 == isite4) {

        X_GC_child_CisAjtCkuAku_Hubbard_MPI(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          org_isite3 - 1, org_sigma3, 1.0, X, nstate, Xvec, vec);

      }
      else if (isite1 != isite2 && isite3 != isite4) {
        X_GC_child_CisAjtCkuAlv_Hubbard_MPI(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          org_isite3 - 1, org_sigma3, org_isite4 - 1, org_sigma4, 1.0, X, nstate, Xvec, vec);
      }

    }//InterPE
    else {
      child_general_int_GetInfo(i, X, org_isite1, org_isite2, org_isite3, org_isite4,
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
#pragma omp parallel for default(none) private(j) shared(vec) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V)
        for (j = 1; j <= i_max; j++) {
          GC_child_CisAisCisAis_element(j, isite1, isite3, tmp_V, nstate, Xvec, vec, X, &tmp_off);
        }
      }
      else if (isite1 == isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) private(j) shared(vec) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V)
        for (j = 1; j <= i_max; j++) {
          GC_child_CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, 
            tmp_V, nstate, Xvec, vec, X, &tmp_off);
        }
      }
      else if (isite1 != isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) private(j) shared(vec) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) 
        for (j = 1; j <= i_max; j++) {
          GC_child_CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, 
            tmp_V, nstate, Xvec, vec, X, &tmp_off);
        }
      }
      else if (isite1 != isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) private(j) shared(vec) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) 
        for (j = 1; j <= i_max; j++) {
          GC_child_CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, 
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
  long int i_max;

  for (i = 0; i < X->Def.NCisAjtCkuAlvDC; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
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
        X_child_CisAisCjtAjt_Hubbard_MPI(org_isite1 - 1, org_sigma1,
          org_isite3 - 1, org_sigma3, 1.0, X, nstate, Xvec, vec);
      }
      else if (isite1 == isite2 && isite3 != isite4) {
        //printf("org_isite1=%d, org_isite2=%d, org_isite3=%d, org_isite4=%d\n", org_isite1, org_isite2, org_isite3, org_isite4);
        X_child_CisAisCjtAku_Hubbard_MPI(org_isite1 - 1, org_sigma1,
          org_isite3 - 1, org_sigma3, org_isite4 - 1, org_sigma4, 1.0, X, nstate, Xvec, vec);
      }
      else if (isite1 != isite2 && isite3 == isite4) {
        X_child_CisAjtCkuAku_Hubbard_MPI(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          org_isite3 - 1, org_sigma3, 1.0, X, nstate, Xvec, vec);

      }
      else if (isite1 != isite2 && isite3 != isite4) {
        X_child_CisAjtCkuAlv_Hubbard_MPI(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          org_isite3 - 1, org_sigma3, org_isite4 - 1, org_sigma4, 1.0, X, nstate, Xvec, vec);
      }
    }//InterPE
    else {
      child_general_int_GetInfo(
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
#pragma omp parallel for default(none) private(j) shared(vec,tmp_V) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2)
        for (j = 1; j <= i_max; j++) {
          child_CisAisCisAis_element(j, isite1, isite3, tmp_V, nstate, Xvec, vec, X, &tmp_off);
        }
      }
      else if (isite1 == isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) private(j) shared(vec,tmp_V) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2)
        for (j = 1; j <= i_max; j++) {
          child_CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, 
            tmp_V, nstate, Xvec, vec, X, &tmp_off);
        }
      }
      else if (isite1 != isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) private(j) shared(vec,tmp_V) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2)
        for (j = 1; j <= i_max; j++) {
          child_CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, 
            tmp_V, nstate, Xvec, vec, X, &tmp_off);
        }
      }
      else if (isite1 != isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) private(j) shared(vec,tmp_V) \
firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2)
        for (j = 1; j <= i_max; j++) {
          child_CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, 
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
  long unsigned int tmp_org_isite1, tmp_org_isite2, tmp_org_isite3, tmp_org_isite4;
  long unsigned int tmp_org_sigma1, tmp_org_sigma2, tmp_org_sigma3, tmp_org_sigma4;
  long unsigned int isA_up, isB_up;
  long unsigned int is1_up, is2_up;
  long unsigned int tmp_off = 0;
  int tmp_sgn, num1, num2, one = 1;
  double complex tmp_V;
  long int i_max;
  double complex dmv;

  i_max = X->Check.idim_max;
  X->Large.mode = M_CORR;

  for (i = 0; i < X->Def.NCisAjtCkuAlvDC; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    tmp_org_isite1 = X->Def.CisAjtCkuAlvDC[i][0] + 1;
    tmp_org_sigma1 = X->Def.CisAjtCkuAlvDC[i][1];
    tmp_org_isite2 = X->Def.CisAjtCkuAlvDC[i][2] + 1;
    tmp_org_sigma2 = X->Def.CisAjtCkuAlvDC[i][3];
    tmp_org_isite3 = X->Def.CisAjtCkuAlvDC[i][4] + 1;
    tmp_org_sigma3 = X->Def.CisAjtCkuAlvDC[i][5];
    tmp_org_isite4 = X->Def.CisAjtCkuAlvDC[i][6] + 1;
    tmp_org_sigma4 = X->Def.CisAjtCkuAlvDC[i][7];
    if (Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X) != 0) {
      //error message will be added
      zclear(nstate, prod[i]);
      continue;
    }

    if (org_isite1 > X->Def.Nsite && org_isite3 > X->Def.Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        is1_up = X->Def.Tpow[org_isite1 - 1];
        is2_up = X->Def.Tpow[org_isite3 - 1];
        num1 = X_SpinGC_CisAis((unsigned long int)myrank + 1, X, is1_up, org_sigma1);
        num2 = X_SpinGC_CisAis((unsigned long int)myrank + 1, X, is2_up, org_sigma3);
        zaxpy_long(i_max*nstate, tmp_V * num1*num2, &vec[1][0], &Xvec[1][0]);
      }
      else if (org_isite1 == org_isite3 && org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) {
        is1_up = X->Def.Tpow[org_isite1 - 1];
        num1 = X_SpinGC_CisAis((unsigned long int)myrank + 1, X, is1_up, org_sigma1);
        zaxpy_long(i_max*nstate, tmp_V * num1, &vec[1][0], &Xvec[1][0]);
      }
      else if (org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) {//exchange
        X_child_general_int_spin_MPIdouble(
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
        num2 = X_SpinGC_CisAis((unsigned long int)myrank + 1, X, is2_up, org_sigma3);
#pragma omp parallel for default(none)shared(vec) \
  firstprivate(i_max, tmp_V, is1_up, org_sigma1, X, num2) private(j, num1)
        for (j = 1; j <= i_max; j++) {
          num1 = X_Spin_CisAis(j, X, is1_up, org_sigma1);
          dmv = tmp_V * num1*num2;
          zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
        }
      }
      else if (org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) {//exchange
        X_child_general_int_spin_MPIsingle(
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
#pragma omp parallel for default(none) private(j) shared(vec) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off, tmp_V)
        for (j = 1; j <= i_max; j++) {
          child_CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4,
            tmp_V, nstate, Xvec, vec, X);
        }
      }
      else if (org_isite1 == org_isite3 && org_sigma1 == org_sigma4 && org_sigma3 == org_sigma2) {
#pragma omp parallel for default(none) private(j, dmv) \
firstprivate(i_max,X,isA_up,org_sigma1, tmp_V) shared(vec, list_1)
        for (j = 1; j <= i_max; j++) {
          dmv = tmp_V * X_Spin_CisAis(j, X, isA_up, org_sigma1);
          zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
        }
      }
      else if (org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) { // exchange
#pragma omp parallel for default(none) private(j, tmp_sgn, dmv) shared(vec) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V)
        for (j = 1; j <= i_max; j++) {
          tmp_sgn = X_child_exchange_spin_element(j, X, isA_up, isB_up, org_sigma2, org_sigma4, &tmp_off);
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
  long unsigned int tmp_org_isite1, tmp_org_isite2, tmp_org_isite3, tmp_org_isite4;
  long unsigned int tmp_org_sigma1, tmp_org_sigma2, tmp_org_sigma3, tmp_org_sigma4;
  long unsigned int tmp_off = 0;
  long unsigned int tmp_off_2 = 0;
  long unsigned int list1_off = 0;
  int num1, one = 1;
  double complex tmp_V;
  long int i_max;
  int tmp_Sz;
  long unsigned int tmp_org = 0;
  i_max = X->Check.idim_max;
  X->Large.mode = M_CORR;

  for (i = 0; i < X->Def.NCisAjtCkuAlvDC; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    tmp_org_isite1 = X->Def.CisAjtCkuAlvDC[i][0] + 1;
    tmp_org_sigma1 = X->Def.CisAjtCkuAlvDC[i][1];
    tmp_org_isite2 = X->Def.CisAjtCkuAlvDC[i][2] + 1;
    tmp_org_sigma2 = X->Def.CisAjtCkuAlvDC[i][3];
    tmp_org_isite3 = X->Def.CisAjtCkuAlvDC[i][4] + 1;
    tmp_org_sigma3 = X->Def.CisAjtCkuAlvDC[i][5];
    tmp_org_isite4 = X->Def.CisAjtCkuAlvDC[i][6] + 1;
    tmp_org_sigma4 = X->Def.CisAjtCkuAlvDC[i][7];

    if (Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X) != 0) {
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
        X_child_CisAisCjuAju_GeneralSpin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        X_child_CisAitCjuAjv_GeneralSpin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4,
          tmp_V, X, nstate, Xvec, vec);
      }
      else {
      }
    }
    else if (org_isite3 > X->Def.Nsite || org_isite1 > X->Def.Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        X_child_CisAisCjuAju_GeneralSpin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        X_child_CisAitCjuAjv_GeneralSpin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4,
          tmp_V, X, nstate, Xvec, vec);
      }
      else {
      }
    }
    else {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp parallel for default(none) private(j, num1) shared(vec,list_1) \
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
#pragma omp parallel for default(none) private(j, num1) \
firstprivate(i_max,X, org_isite1, org_isite3, org_sigma1, org_sigma2, org_sigma3, org_sigma4, tmp_off, tmp_off_2, list1_off, myrank, tmp_V) shared(vec, list_1)
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
        //printf("DEBUG: rank=%d, dam_pr=%lf\n", myrank, creal(dam_pr));
      }
      else {
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
  long unsigned int tmp_org_isite1, tmp_org_isite2, tmp_org_isite3, tmp_org_isite4;
  long unsigned int tmp_org_sigma1, tmp_org_sigma2, tmp_org_sigma3, tmp_org_sigma4;
  long unsigned int isA_up, isB_up;
  long unsigned int tmp_off = 0;
  double complex tmp_V;
  long int i_max;
  i_max = X->Check.idim_max;

  for (i = 0; i < X->Def.NCisAjtCkuAlvDC; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    tmp_org_isite1 = X->Def.CisAjtCkuAlvDC[i][0] + 1;
    tmp_org_sigma1 = X->Def.CisAjtCkuAlvDC[i][1];
    tmp_org_isite2 = X->Def.CisAjtCkuAlvDC[i][2] + 1;
    tmp_org_sigma2 = X->Def.CisAjtCkuAlvDC[i][3];
    tmp_org_isite3 = X->Def.CisAjtCkuAlvDC[i][4] + 1;
    tmp_org_sigma3 = X->Def.CisAjtCkuAlvDC[i][5];
    tmp_org_isite4 = X->Def.CisAjtCkuAlvDC[i][6] + 1;
    tmp_org_sigma4 = X->Def.CisAjtCkuAlvDC[i][7];

    if (Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X) != 0) {
      //error message will be added
      zclear(nstate, prod[i]);
      continue;
    }

    if (org_isite1 > X->Def.Nsite && org_isite3 > X->Def.Nsite) { //org_isite3 >= org_isite1 > Nsite

      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        X_GC_child_CisAisCjuAju_spin_MPIdouble(
          org_isite1 - 1, org_sigma1, (org_isite3 - 1), org_sigma3, tmp_V, X, nstate, Xvec, vec);

      }
      else if (org_isite1 == org_isite3 && org_sigma1 == org_sigma4 && org_sigma2 == org_sigma3) { //diagonal (for spin: cuadcdau=cuau)
        X_GC_child_CisAis_spin_MPIdouble(
          org_isite1 - 1, org_sigma1, tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
        X_GC_child_CisAisCjuAjv_spin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, org_sigma4, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
        X_GC_child_CisAitCjuAju_spin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3,
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        X_GC_child_CisAitCiuAiv_spin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4, 
          tmp_V, X, nstate, Xvec, vec);
      }
    }
    else if (org_isite3 > X->Def.Nsite || org_isite1 > X->Def.Nsite) { //org_isite3 > Nsite >= org_isite1
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        X_GC_child_CisAisCjuAju_spin_MPIsingle(
          org_isite1 - 1, org_sigma1, (org_isite3 - 1), org_sigma3, tmp_V, X, nstate, Xvec, vec);

      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
        X_GC_child_CisAisCjuAjv_spin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, org_sigma4,
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
        X_GC_child_CisAitCjuAju_spin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        X_GC_child_CisAitCiuAiv_spin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4,
          tmp_V, X, nstate, Xvec, vec);
      }
    }
    else {
      if (org_isite1 == org_isite2 && org_isite3 == org_isite4) {
        isA_up = X->Def.Tpow[org_isite2 - 1];
        isB_up = X->Def.Tpow[org_isite4 - 1];
        if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp parallel for default(none) private(j) shared(vec) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V)
          for (j = 1; j <= i_max; j++) {
            GC_child_CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4,
              tmp_V, nstate, Xvec, vec, X);
          }
        }
        else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) private(j) shared(vec) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V)
          for (j = 1; j <= i_max; j++) {
            GC_child_CisAisCitAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up,
              tmp_V, nstate, Xvec, vec, X, &tmp_off);
          }
        }
        else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
#pragma omp parallel for default(none) private(j) shared(vec) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V)
          for (j = 1; j <= i_max; j++) {
            GC_child_CisAitCiuAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, 
              tmp_V, nstate, Xvec, vec, X, &tmp_off);
          }
        }
        else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) private(j) shared(vec) \
firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V)
          for (j = 1; j <= i_max; j++) {
            GC_child_CisAitCiuAiv_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up,
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
  long unsigned int tmp_org_isite1, tmp_org_isite2, tmp_org_isite3, tmp_org_isite4;
  long unsigned int tmp_org_sigma1, tmp_org_sigma2, tmp_org_sigma3, tmp_org_sigma4;
  long unsigned int tmp_off = 0;
  long unsigned int tmp_off_2 = 0;
  int  num1, one = 1;
  double complex tmp_V;
  long int i_max;
  i_max = X->Check.idim_max;
  X->Large.mode = M_CORR;

  for(i=0;i<X->Def.NCisAjtCkuAlvDC;i++){
    zclear(i_max*nstate, &Xvec[1][0]);
    tmp_org_isite1 = X->Def.CisAjtCkuAlvDC[i][0] + 1;
    tmp_org_sigma1 = X->Def.CisAjtCkuAlvDC[i][1];
    tmp_org_isite2 = X->Def.CisAjtCkuAlvDC[i][2] + 1;
    tmp_org_sigma2 = X->Def.CisAjtCkuAlvDC[i][3];
    tmp_org_isite3 = X->Def.CisAjtCkuAlvDC[i][4] + 1;
    tmp_org_sigma3 = X->Def.CisAjtCkuAlvDC[i][5];
    tmp_org_isite4 = X->Def.CisAjtCkuAlvDC[i][6] + 1;
    tmp_org_sigma4 = X->Def.CisAjtCkuAlvDC[i][7];

    if (Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X) != 0) {
      //error message will be added
      zclear(nstate, prod[i]);
      continue;
    }

    if (org_isite1 > X->Def.Nsite && org_isite3 > X->Def.Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        X_GC_child_CisAisCjuAju_GeneralSpin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
        X_GC_child_CisAisCjuAjv_GeneralSpin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, org_sigma4, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
        X_GC_child_CisAitCjuAju_GeneralSpin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        X_GC_child_CisAitCjuAjv_GeneralSpin_MPIdouble(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4,
          tmp_V, X, nstate, Xvec, vec);
      }
    }
    else if (org_isite3 > X->Def.Nsite || org_isite1 > X->Def.Nsite) {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
        X_GC_child_CisAisCjuAju_GeneralSpin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
        X_GC_child_CisAisCjuAjv_GeneralSpin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_isite3 - 1, org_sigma3, org_sigma4, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
        X_GC_child_CisAitCjuAju_GeneralSpin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, 
          tmp_V, X, nstate, Xvec, vec);
      }
      else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
        X_GC_child_CisAitCjuAjv_GeneralSpin_MPIsingle(
          org_isite1 - 1, org_sigma1, org_sigma2, org_isite3 - 1, org_sigma3, org_sigma4, 
          tmp_V, X, nstate, Xvec, vec);
      }
    }
    else {
      if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp parallel for default(none) private(j, num1) shared(vec) \
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
#pragma omp parallel for default(none) private(j, num1) shared(vec) \
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
#pragma omp parallel for default(none) private(j, num1) shared(vec) \
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
firstprivate(i_max,X, org_isite1, org_isite3, org_sigma1, org_sigma2, org_sigma3, org_sigma4, tmp_off, tmp_off_2, tmp_V) shared(vec)
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
  double complex **prod
) {
  int info = 0;
  if (X->Def.iFlgGeneralSpin == FALSE) {
    info = expec_cisajscktalt_SpinGCHalf(X, nstate, Xvec, vec, prod);
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
  char sdt[D_FileNameMax];
  long unsigned int irght, ilft, ihfbit, icaca;
  double complex **prod;
  //For TPQ
  int step = 0, rand_i = 0, istate;

  if (X->Def.NCisAjtCkuAlvDC < 1) return 0;
  X->Large.mode = M_CORR;

  if (GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit) != 0) {
    return -1;
  }

  //Make File Name for output
  prod = cd_2d_allocate(X->Def.NCisAjt, nstate);
  switch (X->Def.iCalcType) {
  case TPQCalc:
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
    if (expec_cisajscktalt_SpinGC(X, nstate, Xvec, vec, prod) != 0) {
      return -1;
    }
    break;

  default:
    return -1;
  }

  for (istate = 0; istate < nstate; istate++) {
    switch (X->Def.iCalcModel) {
    case TPQCalc:
      step = X->Def.istep;
      sprintf(sdt, cFileName2BGreen_TPQ, X->Def.CDataFileHead, istate, step);
      break;
    case TimeEvolution:
      step = X->Def.istep;
      sprintf(sdt, cFileName2BGreen_TE, X->Def.CDataFileHead, step);
      break;
    case FullDiag:
    case CG:
      sprintf(sdt, cFileName2BGreen_FullDiag, X->Def.CDataFileHead, istate);
      break;
    }
    if (childfopenMPI(sdt, "w", &fp) == 0) {
      for (icaca = 0; icaca < X->Def.NCisAjt; icaca++) {
        fprintf(fp, " %4d %4d %4d %4d %4d %4d %4d %4d %.10lf %.10lf\n",
          X->Def.CisAjtCkuAlvDC[icaca][0], X->Def.CisAjtCkuAlvDC[icaca][1],
          X->Def.CisAjtCkuAlvDC[icaca][2], X->Def.CisAjtCkuAlvDC[icaca][3],
          X->Def.CisAjtCkuAlvDC[icaca][4], X->Def.CisAjtCkuAlvDC[icaca][5],
          X->Def.CisAjtCkuAlvDC[icaca][6], X->Def.CisAjtCkuAlvDC[icaca][7],
          creal(prod[icaca][istate]), cimag(prod[icaca][istate]));
      }
      fclose(fp);
    }
    else return -1;
  }/*for (istate = 0; istate < nstate; istate++)*/

  if (X->Def.iCalcType == TPQCalc) {
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
  return 0;
}
