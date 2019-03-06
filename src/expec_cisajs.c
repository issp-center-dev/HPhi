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

#include "mltplyCommon.h"
#include "mltply.h"
#include "FileIO.h"
#include "bitcalc.h"
#include "wrapperMPI.h"
#include "mltplyHubbard.h"
#include "mltplyHubbardCore.h"
#include "mltplySpinCore.h"
#include "mltplyMPIHubbard.h"
#include "mltplyMPISpinCore.h"
#include "common/setmemory.h"

/**
 * @file   expec_cisajs.c
 * 
 * @brief  File for calculation of one body green's function
 *
 * @version 0.1, 0.2
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 */
/** 
 * @brief function of calculation for one body green's function
 * 
 * @param X [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvectors.
 * 
 * @version 0.2
 * @details add calculation one body green's functions for general spin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec
){
  FILE *fp;
  char sdt[D_FileNameMax];
  double complex **prod;
  long unsigned int irght, ilft, ihfbit, ica;
  long int i_max;
  //For TPQ
  int step = 0, rand_i = 0, istate;

  if(X->Def.NCisAjt <1) return 0;

  i_max = X->Check.idim_max;      
  if(GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit)!=0){
    return -1;
  }
  X->Large.i_max    = i_max;
  X->Large.irght    = irght;
  X->Large.ilft     = ilft;
  X->Large.ihfbit   = ihfbit;
  X->Large.mode     = M_CORR;
 
  switch(X->Def.iCalcType){
  case TPQCalc:
    step=X->Def.istep;
    TimeKeeperWithRandAndStep(X, cFileNameTimeKeep,  cTPQExpecOneBodyGStart, "a", 0, step);
    break;
  case TimeEvolution:
    step = X->Def.istep;
    TimeKeeperWithStep(X, cFileNameTimeKeep, cTEExpecOneBodyGStart, "a", step);
    break;
  case FullDiag:
  case CG:
    break;
  }
  
  prod = cd_2d_allocate(X->Def.NCisAjt, nstate);
  switch(X->Def.iCalcModel){
  case HubbardGC:
    if(expec_cisajs_HubbardGC(X, nstate, Xvec, vec, prod)!=0){
      return -1;
    }
    break;
    
  case KondoGC:
  case Hubbard:
  case Kondo:
    if (expec_cisajs_Hubbard(X, nstate, Xvec, vec, prod) != 0) {
      return -1;
    }
    break;

  case Spin: // for the Sz-conserved spin system
    if (expec_cisajs_Spin(X, nstate, Xvec, vec, prod) != 0) {
      return -1;
    }
    break;
    
  case SpinGC:
    if (expec_cisajs_SpinGC(X, nstate, Xvec, vec, prod) != 0) {
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
      sprintf(sdt, cFileName1BGreen_TPQ, X->Def.CDataFileHead, istate, step);
      break;
    case TimeEvolution:
      step = X->Def.istep;
      sprintf(sdt, cFileName1BGreen_TE, X->Def.CDataFileHead, step);
      break;
    case FullDiag:
    case CG:
      sprintf(sdt, cFileName1BGreen_FullDiag, X->Def.CDataFileHead, istate);
      break;
    }
    if (childfopenMPI(sdt, "w", &fp) == 0) {
      for (ica = 0; ica < X->Def.NCisAjt; ica++) {
        fprintf(fp, " %4ld %4ld %4ld %4ld %.10lf %.10lf\n",
          X->Def.CisAjt[ica][0], X->Def.CisAjt[ica][1], X->Def.CisAjt[ica][2], X->Def.CisAjt[ica][3],
          creal(prod[ica][istate]), cimag(prod[ica][istate]));
      }
      fclose(fp);
    }
    else return -1;
  }/*for (istate = 0; istate < nstate; istate++)*/

  if(X->Def.St==0){
    if(X->Def.iCalcType==TPQCalc){
      TimeKeeperWithRandAndStep(X, cFileNameTimeKeep, cTPQExpecOneBodyGFinish, "a", rand_i, step);     
    }
    else if(X->Def.iCalcType==TimeEvolution){
      TimeKeeperWithStep(X, cFileNameTimeKeep, cTEExpecOneBodyGFinish, "a", step);
    }
  }else if(X->Def.St==1){
    TimeKeeper(X, cFileNameTimeKeep, cCGExpecOneBodyGFinish, "a");
    fprintf(stdoutMPI, "%s", cLogCGExpecOneBodyGEnd);
  }
  free_cd_2d_allocate(prod);
  return 0;
}
/**
 * @brief function of calculation for one body green's function for Hubbard GC model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_HubbardGC(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec, 
  double complex **prod
){
  long unsigned int i, j;
  long unsigned int org_isite1, org_isite2, org_sigma1, org_sigma2;
  double complex dam_pr = 0;
  long int i_max;
  long int ibit;
  long unsigned int is;
  double complex tmp_OneGreen = 1.0;

  i_max = X->Check.idim_max;

  for (i = 0; i < X->Def.NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    org_isite1 = X->Def.CisAjt[i][0] + 1;
    org_isite2 = X->Def.CisAjt[i][2] + 1;
    org_sigma1 = X->Def.CisAjt[i][1];
    org_sigma2 = X->Def.CisAjt[i][3];
    dam_pr = 0;
    if (org_isite1 > X->Def.Nsite &&
      org_isite2 > X->Def.Nsite) {
      if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {
        if (org_sigma1 == 0) {
          is = X->Def.Tpow[2 * org_isite1 - 2];
        }
        else {
          is = X->Def.Tpow[2 * org_isite1 - 1];
        }
        ibit = (unsigned long int)myrank & is;
        if (ibit == is) {
          zaxpy_long(i_max*nstate, tmp_OneGreen, &vec[1][0], &Xvec[1][0]);
        }
      }
      else {
        X_GC_child_general_hopp_MPIdouble(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          -tmp_OneGreen, X, nstate, Xvec, vec);
      }
    }
    else if (org_isite2 > X->Def.Nsite || org_isite1 > X->Def.Nsite) {
      if (org_isite1 < org_isite2) {
        X_GC_child_general_hopp_MPIsingle(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          -tmp_OneGreen, X, nstate, Xvec, vec);
      }
      else {
        X_GC_child_general_hopp_MPIsingle(org_isite2 - 1, org_sigma2, org_isite1 - 1, org_sigma1, 
          -tmp_OneGreen, X, nstate, Xvec, vec);
        zswap_long(i_max*nstate, &vec[1][0], &Xvec[1][0]);
      }
    }
    else {
      if (child_general_hopp_GetInfo(X, org_isite1, org_isite2, org_sigma1, org_sigma2) != 0) {
        return -1;
      }
      GC_child_general_hopp(nstate, Xvec, vec, X, tmp_OneGreen);
    }

    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief function of calculation for one body green's function for Hubbard model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_Hubbard(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec, 
  double complex **prod
) {
  long unsigned int i, j;
  long unsigned int org_isite1, org_isite2, org_sigma1, org_sigma2;
  double complex dam_pr = 0;
  long int i_max;
  int num1;
  long int ibit;
  long unsigned int is;
  double complex tmp_OneGreen = 1.0;

  i_max = X->Check.idim_max;
  for (i = 0; i < X->Def.NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    org_isite1 = X->Def.CisAjt[i][0] + 1;
    org_isite2 = X->Def.CisAjt[i][2] + 1;
    org_sigma1 = X->Def.CisAjt[i][1];
    org_sigma2 = X->Def.CisAjt[i][3];
    dam_pr = 0.0;

    if (X->Def.iFlgSzConserved == TRUE) {
      if (org_sigma1 != org_sigma2) {
        zclear(nstate, prod[i]);
        continue;
      }
    }

    if (X->Def.iCalcModel == Kondo || X->Def.iCalcModel == KondoGC) {
      if ((X->Def.LocSpn[org_isite1 - 1] == 1 && X->Def.LocSpn[org_isite2 - 1] == 0) ||
          (X->Def.LocSpn[org_isite1 - 1] == 0 && X->Def.LocSpn[org_isite2 - 1] == 1)
        )
      {
        zclear(nstate, prod[i]);
        continue;
      }
    }

    if (org_isite1 > X->Def.Nsite &&
      org_isite2 > X->Def.Nsite) {
      if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {//diagonal

        is = X->Def.Tpow[2 * org_isite1 - 2 + org_sigma1];
        ibit = (unsigned long int)myrank & is;
        if (ibit == is) {
          zaxpy_long(i_max*nstate, tmp_OneGreen, &vec[1][0], &Xvec[1][0]);
        }
      }
      else {
        X_child_general_hopp_MPIdouble(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2, 
          -tmp_OneGreen, X, nstate, Xvec, vec);
      }
    }
    else if (org_isite2 > X->Def.Nsite || org_isite1 > X->Def.Nsite) {
      if (org_isite1 < org_isite2) {
        X_child_general_hopp_MPIsingle(org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
          -tmp_OneGreen, X, nstate, Xvec, vec);
      }
      else {
        X_child_general_hopp_MPIsingle(org_isite2 - 1, org_sigma2, org_isite1 - 1, org_sigma1, 
          -tmp_OneGreen, X, nstate, Xvec, vec);
        zswap_long(i_max*nstate, &vec[1][0], &Xvec[1][0]);
      }
    }
    else {
      if (child_general_hopp_GetInfo(X, org_isite1, org_isite2, org_sigma1, org_sigma2) != 0) {
        return -1;
      }
      if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {
        is = X->Def.Tpow[2 * org_isite1 - 2 + org_sigma1];

#pragma omp parallel for default(none) shared(list_1, vec) reduction(+:dam_pr) \
firstprivate(i_max, is) private(num1, ibit)
        for (j = 1; j <= i_max; j++) {
          ibit = list_1[j] & is;
          num1 = ibit / is;
          zaxpy_(&nstate, &tmp_OneGreen, vec[j], &one, Xvec[j], &one);
        }
      }
      else {
        child_general_hopp(nstate, Xvec, vec, X, tmp_OneGreen);
      }
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief function of calculation for one body green's function for Spin model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_Spin(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec, 
    double complex **prod
) {
  int info = 0;
  if (X->Def.iFlgGeneralSpin == FALSE) {
    info = expec_cisajs_SpinHalf(X, nstate, Xvec, vec, prod);
  }
  else {
    info = expec_cisajs_SpinGeneral(X, nstate, Xvec, vec, prod);
  }
  return info;
}
/**
 * @brief function of calculation for one body green's function for Half-Spin model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_SpinHalf(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec, 
  double complex **prod
) {
  long unsigned int i, j;
  long unsigned int isite1;
  long unsigned int org_isite1, org_isite2, org_sigma1, org_sigma2;
  double complex dam_pr = 0, dmv;
  long int i_max;
  long int ibit1;
  long unsigned int is1_up;
  int one = 1;

  i_max = X->Check.idim_max;

  for (i = 0; i < X->Def.NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    org_isite1 = X->Def.CisAjt[i][0] + 1;
    org_isite2 = X->Def.CisAjt[i][2] + 1;
    org_sigma1 = X->Def.CisAjt[i][1];
    org_sigma2 = X->Def.CisAjt[i][3];

    if (org_sigma1 == org_sigma2) {
      if (org_isite1 == org_isite2) {
        if (org_isite1 > X->Def.Nsite) {
          is1_up = X->Def.Tpow[org_isite1 - 1];
          ibit1 = X_SpinGC_CisAis((unsigned long int)myrank + 1, X, is1_up, org_sigma1);
          if (ibit1 != 0) {
            zaxpy_long(i_max*nstate, 1.0, &vec[1][0], &Xvec[1][0]);
          }
        }// org_isite1 > X->Def.Nsite
        else {
          isite1 = X->Def.Tpow[org_isite1 - 1];
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) \
firstprivate(i_max, isite1, org_sigma1, X) shared(vec)
          for (j = 1; j <= i_max; j++) {
            dmv = X_Spin_CisAis(j, X, isite1, org_sigma1);
            zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
          }
        }
      }
      else {
        dam_pr = 0.0;
      }
    }
    else {
      // for the canonical case
      dam_pr = 0.0;
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief function of calculation for one body green's function for General-Spin model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_SpinGeneral(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec, 
  double complex **prod
) {
  long unsigned int i, j;
  long unsigned int org_isite1, org_isite2, org_sigma1, org_sigma2;
  double complex dam_pr = 0, dmv;
  long int i_max;
  int num1, one = 1;
  i_max = X->Check.idim_max;

  for (i = 0; i < X->Def.NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    org_isite1 = X->Def.CisAjt[i][0] + 1;
    org_isite2 = X->Def.CisAjt[i][2] + 1;
    org_sigma1 = X->Def.CisAjt[i][1];
    org_sigma2 = X->Def.CisAjt[i][3];

    if (org_isite1 == org_isite2) {
      if (org_isite1 > X->Def.Nsite) {
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
          num1 = BitCheckGeneral((unsigned long int)myrank,
                                           org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
          dam_pr = 0.0;
          if (num1 != 0) {
            zaxpy_long(i_max*nstate, 1.0, &vec[1][0], &Xvec[1][0]);
          }
        }
        else {
          dam_pr = 0.0;
        }
      }
      else {//org_isite1 <= X->Def.Nsite
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
          dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) \
firstprivate(i_max, org_isite1, org_sigma1, X) shared(vec, list_1)
          for (j = 1; j <= i_max; j++) {
            dmv = BitCheckGeneral(list_1[j], org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
            zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
          }
        }
        else {
          dam_pr = 0.0;
        }
      }
    }
    else {
      // hopping is not allowed in localized spin system
      dam_pr = 0.0;
    }//org_isite1 != org_isite2

    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief function of calculation for one body green's function for SpinGC model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_SpinGC(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec, 
  double complex **prod
) {
  int info = 0;
  if (X->Def.iFlgGeneralSpin == FALSE) {
    info = expec_cisajs_SpinGCHalf(X, nstate, Xvec, vec, _fp);
  }
  else {
    info = expec_cisajs_SpinGCGeneral(X, nstate, Xvec, vec, _fp);
  }
  return info;
}
/**
 * @brief function of calculation for one body green's function for Half-SpinGC model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_SpinGCHalf(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec, 
    double complex **prod
) {
  long unsigned int i, j;
  long unsigned int isite1;
  long unsigned int org_isite1, org_isite2, org_sigma1, org_sigma2;
  double complex dam_pr = 0, dmv;
  long int i_max;
  int tmp_sgn, one = 1;
  long unsigned int tmp_off = 0;

  i_max = X->Check.idim_max;

  for (i = 0; i < X->Def.NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    org_isite1 = X->Def.CisAjt[i][0] + 1;
    org_isite2 = X->Def.CisAjt[i][2] + 1;
    org_sigma1 = X->Def.CisAjt[i][1];
    org_sigma2 = X->Def.CisAjt[i][3];
    dam_pr = 0.0;

    if (org_isite1 == org_isite2) {
      if (org_isite1 > X->Def.Nsite) {
        if (org_sigma1 == org_sigma2) {  // longitudinal magnetic field
          X_GC_child_CisAis_spin_MPIdouble(org_isite1 - 1, org_sigma1, 1.0, X, nstate, Xvec, vec);
        }
        else {  // transverse magnetic field
          X_GC_child_CisAit_spin_MPIdouble(org_isite1 - 1, org_sigma1, org_sigma2, 1.0, X, nstate, Xvec, vec);
        }
      }
      else {
        isite1 = X->Def.Tpow[org_isite1 - 1];

        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn) \
firstprivate(i_max, isite1, org_sigma1, X) shared(vec)
          for (j = 1; j <= i_max; j++) {
            dmv = X_SpinGC_CisAis(j, X, isite1, org_sigma1);
            zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
          }
        }
        else {
          // transverse magnetic field
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, tmp_off) \
firstprivate(i_max, isite1, org_sigma2, X) shared(vec)
          for (j = 1; j <= i_max; j++) {
            tmp_sgn = X_SpinGC_CisAit(j, X, isite1, org_sigma2, &tmp_off);
            if (tmp_sgn != 0) {
              dmv = (double complex)tmp_sgn;
              zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[tmp_off + 1][0], &one);
            }
          }
        }
      }
    }
    else {
      // hopping is not allowed in localized spin system
      dam_pr = 0.0;
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
/**
 * @brief function of calculation for one body green's function for General SpinGC model.
 *
 * @param X  [in] list for getting information to calculate one body green's function.
 * @param vec [in] eigenvector
 * @param _fp [in] pointer to output file
 * @retval 0 normally finished.
 * @retval -1 abnormally finished.
 */
int expec_cisajs_SpinGCGeneral(
  struct BindStruct *X,
  int nstate,
  double complex **Xvec,
  double complex **vec,
  double complex **prod
) {
  long unsigned int i, j;
  long unsigned int org_isite1, org_isite2, org_sigma1, org_sigma2;
  double complex dam_pr = 0, dmv;
  long int i_max;
  long unsigned int tmp_off = 0;
  int num1, one = 1;

  i_max = X->Check.idim_max;

  for (i = 0; i < X->Def.NCisAjt; i++) {
    zclear(i_max*nstate, &Xvec[1][0]);
    org_isite1 = X->Def.CisAjt[i][0] + 1;
    org_isite2 = X->Def.CisAjt[i][2] + 1;
    org_sigma1 = X->Def.CisAjt[i][1];
    org_sigma2 = X->Def.CisAjt[i][3];
    if (org_isite1 == org_isite2) {
      if (org_isite1 > X->Def.Nsite) {
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
          X_GC_child_CisAis_GeneralSpin_MPIdouble(org_isite1 - 1, org_sigma1,
            1.0, X, nstate, Xvec, vec);
        }
        else {
          // transverse magnetic field
          X_GC_child_CisAit_GeneralSpin_MPIdouble(
            org_isite1 - 1, org_sigma1, org_sigma2, 1.0, X, nstate, Xvec, vec);
        }
      }
      else {//org_isite1 <= X->Def.Nsite
        if (org_sigma1 == org_sigma2) {
          // longitudinal magnetic field
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) \
firstprivate(i_max, org_isite1, org_sigma1, X) shared(vec)
          for (j = 1; j <= i_max; j++) {
            num1 = BitCheckGeneral(j - 1, org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
            dmv = (double complex)num1;
            zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[j][0], &one);
          }
        }
        else {
          // transverse magnetic field
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) \
firstprivate(i_max, org_isite1, org_sigma1, org_sigma2, X,tmp_off) shared(vec)
          for (j = 1; j <= i_max; j++) {
            num1 = GetOffCompGeneralSpin(
              j - 1, org_isite1, org_sigma2, org_sigma1, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
            if (num1 != 0) {
              dmv = (double complex)num1;
              zaxpy_(&nstate, &dmv, &vec[j][0], &one, &Xvec[tmp_off + 1][0], &one);
            }
          }
        }
      }
    }
    MultiVecProdMPI(i_max, nstate, vec, Xvec, prod[i]);
  }
  return 0;
}
