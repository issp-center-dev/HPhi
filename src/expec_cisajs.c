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


int expec_cisajs_HubbardGC(struct BindStruct *X,double complex *vec, FILE **_fp);
int expec_cisajs_Hubbard(struct BindStruct *X,double complex *vec, FILE **_fp);

int expec_cisajs_Spin(struct BindStruct *X,double complex *vec, FILE **_fp);
int expec_cisajs_SpinHalf(struct BindStruct *X,double complex *vec, FILE **_fp);
int expec_cisajs_SpinGeneral(struct BindStruct *X,double complex *vec, FILE **_fp);

int expec_cisajs_SpinGC(struct BindStruct *X,double complex *vec, FILE **_fp);
int expec_cisajs_SpinGCHalf(struct BindStruct *X,double complex *vec, FILE **_fp);
int expec_cisajs_SpinGCGeneral(struct BindStruct *X,double complex *vec, FILE **_fp);



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
int expec_cisajs(struct BindStruct *X,double complex *vec){

  FILE *fp;
  char sdt[D_FileNameMax];

  long unsigned int irght,ilft,ihfbit;
  long int i_max;
  //For TPQ
  int step=0;
  int rand_i=0;

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
  case Lanczos:
    if(X->Def.St==0){
      sprintf(sdt, cFileName1BGreen_Lanczos, X->Def.CDataFileHead);
        fprintf(stdoutMPI, "%s", cLogLanczosExpecOneBodyGStart);
      TimeKeeper(X, cFileNameTimeKeep, cLanczosExpecOneBodyGStart, "a");
    }else if(X->Def.St==1){
      sprintf(sdt, cFileName1BGreen_CG, X->Def.CDataFileHead);
        TimeKeeper(X, cFileNameTimeKeep, cCGExpecOneBodyGStart, "a");
        fprintf(stdoutMPI, "%s", cLogCGExpecOneBodyGStart);
    }
    //vec=v0;
    break;
  case TPQCalc:
  case cTPQ:
    step=X->Def.istep;
    rand_i=X->Def.irand;
    TimeKeeperWithRandAndStep(X, cFileNameTimeKeep,  cTPQExpecOneBodyGStart, "a", rand_i, step);
    sprintf(sdt, cFileName1BGreen_TPQ, X->Def.CDataFileHead, rand_i, step);
    //vec=v0;
    break;
  case TimeEvolution:
      step=X->Def.istep;
      TimeKeeperWithStep(X, cFileNameTimeKeep,  cTEExpecOneBodyGStart, "a", step);
      sprintf(sdt, cFileName1BGreen_TE, X->Def.CDataFileHead, step);
      break;

  case FullDiag:
  case CG:
    sprintf(sdt, cFileName1BGreen_FullDiag, X->Def.CDataFileHead, X->Phys.eigen_num);
    //vec=v0;
    break;
  }
  
  if(childfopenMPI(sdt, "w", &fp)!=0){
    return -1;
  } 
  switch(X->Def.iCalcModel){
  case HubbardGC:
    if(expec_cisajs_HubbardGC(X, vec, &fp)!=0){
        return -1;
    }
    break;
    
  case KondoGC:
  case Hubbard:
  case Kondo:
      if(expec_cisajs_Hubbard(X, vec, &fp)!=0){
          return -1;
      }
    break;

  case Spin: // for the Sz-conserved spin system
      if(expec_cisajs_Spin(X, vec, &fp)!=0){
          return -1;
      }
    break;
    
  case SpinGC:
      if(expec_cisajs_SpinGC(X, vec, &fp)!=0){
          return -1;
      }
          break;
        
  default:
    return -1;
  }

  fclose(fp);
  if(X->Def.St==0){
    if(X->Def.iCalcType==Lanczos){
      TimeKeeper(X, cFileNameTimeKeep, cLanczosExpecOneBodyGFinish, "a");
      fprintf(stdoutMPI, "%s", cLogLanczosExpecOneBodyGEnd);
      TimeKeeper(X, cFileNameTimeKeep, cLanczosExpecOneBodyGFinish, "a");
    }
    else if(X->Def.iCalcType==TPQCalc){
      TimeKeeperWithRandAndStep(X, cFileNameTimeKeep, cTPQExpecOneBodyGFinish, "a", rand_i, step);     
    }
    else if(X->Def.iCalcType==TimeEvolution){
      TimeKeeperWithStep(X, cFileNameTimeKeep, cTEExpecOneBodyGFinish, "a", step);
    }
  }else if(X->Def.St==1){
    TimeKeeper(X, cFileNameTimeKeep, cCGExpecOneBodyGFinish, "a");
    fprintf(stdoutMPI, "%s", cLogCGExpecOneBodyGEnd);
  }
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
int expec_cisajs_HubbardGC(struct BindStruct *X, double complex *vec, FILE **_fp){
    long unsigned int i,j;
    long unsigned int org_isite1,org_isite2,org_sigma1,org_sigma2;
    double complex dam_pr=0;
    long int i_max;
    long int ibit;
    long unsigned int is;
    double complex tmp_OneGreen=1.0;

    i_max = X->Check.idim_max;

    for(i=0;i<X->Def.NCisAjt;i++){
        org_isite1 = X->Def.CisAjt[i][0]+1;
        org_isite2 = X->Def.CisAjt[i][2]+1;
        org_sigma1 = X->Def.CisAjt[i][1];
        org_sigma2 = X->Def.CisAjt[i][3];
        dam_pr=0;
        if (org_isite1  > X->Def.Nsite &&
            org_isite2  > X->Def.Nsite) {
            if(org_isite1==org_isite2 && org_sigma1==org_sigma2){
                if(org_sigma1==0){
                    is   = X->Def.Tpow[2 * org_isite1 - 2];
                }
                else{
                    is = X->Def.Tpow[2 * org_isite1 - 1];
                }
                ibit = (unsigned long int)myrank & is;
                if (ibit == is) {
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(vec)  \
  firstprivate(i_max) private(j)
                    for (j = 1; j <= i_max; j++) dam_pr += vec[j]*conj(vec[j]);
                }
            }
            else{
                dam_pr =child_GC_general_hopp_MPIdouble(org_isite1-1, org_sigma1, org_isite2-1, org_sigma2, -tmp_OneGreen, X, vec, vec);
            }
        }
        else if (org_isite2  > X->Def.Nsite || org_isite1  > X->Def.Nsite){
            if(org_isite1<org_isite2){
                dam_pr =child_GC_general_hopp_MPIsingle(org_isite1-1, org_sigma1, org_isite2-1, org_sigma2, -tmp_OneGreen, X, vec, vec);
            }
            else{
                dam_pr =child_GC_general_hopp_MPIsingle(org_isite2-1, org_sigma2, org_isite1-1, org_sigma1, -tmp_OneGreen, X, vec, vec);
                dam_pr = conj(dam_pr);
            }
        }
        else{
            if(general_hopp_GetInfo( X,org_isite1,org_isite2,org_sigma1,org_sigma2)!=0){
                return -1;
            }
            dam_pr = GC_general_hopp(vec,vec,X,tmp_OneGreen);
        }

        dam_pr= SumMPI_dc(dam_pr);
        fprintf(*_fp," %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,creal(dam_pr),cimag(dam_pr));
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
int expec_cisajs_Hubbard(struct BindStruct *X, double complex *vec, FILE **_fp) {
    long unsigned int i,j;
    long unsigned int org_isite1,org_isite2,org_sigma1,org_sigma2;
    double complex dam_pr=0;
    long int i_max;
    int num1;
    long int ibit;
    long unsigned int is;
    double complex tmp_OneGreen=1.0;

    i_max = X->Check.idim_max;
    for(i=0;i<X->Def.NCisAjt;i++){
        org_isite1 = X->Def.CisAjt[i][0]+1;
        org_isite2 = X->Def.CisAjt[i][2]+1;
        org_sigma1 = X->Def.CisAjt[i][1];
        org_sigma2 = X->Def.CisAjt[i][3];
        dam_pr=0.0;

        if(X->Def.iFlgSzConserved ==TRUE){
            if(org_sigma1 != org_sigma2){
                dam_pr =0.0;
                fprintf(*_fp," %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,creal(dam_pr),cimag(dam_pr));
                continue;
            }
        }

        if(X->Def.iCalcModel==Kondo || X->Def.iCalcModel==KondoGC) {
          if( (X->Def.LocSpn[org_isite1 - 1] == 1 && X->Def.LocSpn[org_isite2 - 1] == 0) ||
                  (X->Def.LocSpn[org_isite1 - 1] == 0 && X->Def.LocSpn[org_isite2 - 1] == 1)
                  )
          {
            dam_pr =0.0;
            fprintf(*_fp," %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,creal(dam_pr),cimag(dam_pr));
            continue;
          }
        }

        if (org_isite1  > X->Def.Nsite &&
            org_isite2  > X->Def.Nsite) {
          if(org_isite1==org_isite2 && org_sigma1==org_sigma2){//diagonal

                is   = X->Def.Tpow[2 * org_isite1 - 2+org_sigma1];
                ibit = (unsigned long int)myrank & is;
                if (ibit == is) {
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(vec) \
  firstprivate(i_max) private(j)
                    for (j = 1; j <= i_max; j++) dam_pr += vec[j]*conj(vec[j]);
                }

            }
            else{
                dam_pr =child_general_hopp_MPIdouble(org_isite1-1, org_sigma1, org_isite2-1, org_sigma2, -tmp_OneGreen, X, vec, vec);
            }
        }
        else if (org_isite2  > X->Def.Nsite || org_isite1  > X->Def.Nsite){
          if(org_isite1 < org_isite2){
             dam_pr =child_general_hopp_MPIsingle(org_isite1-1, org_sigma1,org_isite2-1, org_sigma2, -tmp_OneGreen, X, vec, vec);
            }
            else{
            dam_pr = child_general_hopp_MPIsingle(org_isite2-1, org_sigma2, org_isite1-1, org_sigma1, -tmp_OneGreen, X, vec, vec);
                dam_pr = conj(dam_pr);
            }
        }
        else{
          if(general_hopp_GetInfo( X,org_isite1,org_isite2,org_sigma1,org_sigma2)!=0){
                return -1;
            }
            if(org_isite1==org_isite2 && org_sigma1==org_sigma2){
              //fprintf(stdoutMPI,"DEBUG1-3-1\n");
              is   = X->Def.Tpow[2 * org_isite1 - 2 + org_sigma1];

#pragma omp parallel for default(none) shared(list_1, vec) reduction(+:dam_pr) firstprivate(i_max, is) private(num1, ibit)
                for(j = 1;j <= i_max;j++){
                    ibit = list_1[j]&is;
                    num1  = ibit/is;
                    dam_pr += num1*conj(vec[j])*vec[j];
                }
            }
            else{
                dam_pr = general_hopp(vec,vec,X,tmp_OneGreen);
            }
        }
        dam_pr= SumMPI_dc(dam_pr);
      //fprintf(stdoutMPI, "rank=%d, dam_pr=%lf\n", myrank, creal(dam_pr));
      fprintf(*_fp," %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,creal(dam_pr),cimag(dam_pr));
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
int expec_cisajs_Spin(struct BindStruct *X, double complex *vec, FILE **_fp) {
    int info=0;
    if (X->Def.iFlgGeneralSpin == FALSE) {
        info=expec_cisajs_SpinHalf(X,vec, _fp);
    } else {
        info=expec_cisajs_SpinGeneral(X,vec, _fp);
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
int expec_cisajs_SpinHalf(struct BindStruct *X, double complex *vec, FILE **_fp) {
    long unsigned int i,j;
    long unsigned int isite1;
    long unsigned int org_isite1,org_isite2,org_sigma1,org_sigma2;
    double complex dam_pr=0;
    long int i_max;
    long int ibit1;
    long unsigned int is1_up;

    i_max = X->Check.idim_max;

    for(i=0;i<X->Def.NCisAjt;i++){
        org_isite1 = X->Def.CisAjt[i][0]+1;
        org_isite2 = X->Def.CisAjt[i][2]+1;
        org_sigma1 = X->Def.CisAjt[i][1];
        org_sigma2 = X->Def.CisAjt[i][3];

        if(org_sigma1 == org_sigma2){
            if(org_isite1==org_isite2){
                if(org_isite1 > X->Def.Nsite){
                    is1_up = X->Def.Tpow[org_isite1 - 1];
                    ibit1 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is1_up, org_sigma1);
                    dam_pr=0;
                    if(ibit1 !=0){
#pragma omp parallel for reduction(+:dam_pr)default(none) shared(vec) \
  firstprivate(i_max) private(j)
                        for (j = 1; j <= i_max; j++) dam_pr += conj(vec[j])*vec[j];
                    }
                }// org_isite1 > X->Def.Nsite
                else{
                    isite1     = X->Def.Tpow[org_isite1-1];
                    dam_pr=0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max, isite1, org_sigma1, X) shared(vec)
                    for(j=1;j<=i_max;j++){
                        dam_pr+=child_Spin_CisAis(j,X, isite1,org_sigma1)*conj(vec[j])*vec[j];
                    }
                }
            }
            else{
                dam_pr=0.0;
            }
        }else{
            // for the canonical case
            dam_pr =0.0;
        }
        dam_pr = SumMPI_dc(dam_pr);
        fprintf(*_fp," %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1, org_sigma1, org_isite2-1, org_sigma2, creal(dam_pr), cimag(dam_pr));
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
int expec_cisajs_SpinGeneral(struct BindStruct *X, double complex *vec, FILE **_fp) {
    long unsigned int i,j;
    long unsigned int org_isite1,org_isite2,org_sigma1,org_sigma2;
    double complex dam_pr=0;
    long int i_max;
    int num1;
    i_max = X->Check.idim_max;

    for(i=0;i<X->Def.NCisAjt;i++){
        org_isite1 = X->Def.CisAjt[i][0]+1;
        org_isite2 = X->Def.CisAjt[i][2]+1;
        org_sigma1 = X->Def.CisAjt[i][1];
        org_sigma2 = X->Def.CisAjt[i][3];

        if(org_isite1 == org_isite2){
            if(org_isite1 >X->Def.Nsite){
                if(org_sigma1==org_sigma2){
                    // longitudinal magnetic field
                    num1 = BitCheckGeneral((unsigned long int)myrank,
                                           org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
                    dam_pr=0.0;
                    if (num1 != 0) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max, org_isite1, org_sigma1, X) shared(vec)
                        for(j=1;j<=i_max;j++){
                            dam_pr+=conj(vec[j])*vec[j];
                        }
                    }
                }else{
                    dam_pr=0.0;
                }
            }
            else {//org_isite1 <= X->Def.Nsite
                if(org_sigma1==org_sigma2){
                    // longitudinal magnetic field
                    dam_pr=0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, X) shared(vec, list_1)
                    for(j=1;j<=i_max;j++){
                        num1 = BitCheckGeneral(list_1[j], org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
                        dam_pr+=conj(vec[j])*vec[j]*num1;
                    }
                }else{
                    dam_pr=0.0;
                }
            }
        }else{
            // hopping is not allowed in localized spin system
            dam_pr=0.0;
        }//org_isite1 != org_isite2

        dam_pr = SumMPI_dc(dam_pr);
        fprintf(*_fp," %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1, org_sigma1, org_isite2-1, org_sigma2,creal(dam_pr),cimag(dam_pr));
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
int expec_cisajs_SpinGC(struct BindStruct *X, double complex *vec, FILE **_fp) {
    int info=0;
    if (X->Def.iFlgGeneralSpin == FALSE) {
        info=expec_cisajs_SpinGCHalf(X,vec, _fp);
    } else {
        info=expec_cisajs_SpinGCGeneral(X,vec, _fp);
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
int expec_cisajs_SpinGCHalf(struct BindStruct *X, double complex *vec, FILE **_fp) {
    long unsigned int i,j;
    long unsigned int isite1;
    long unsigned int org_isite1,org_isite2,org_sigma1,org_sigma2;
    double complex dam_pr=0;
    long int i_max;
    int tmp_sgn;
    long unsigned int tmp_off=0;

    i_max = X->Check.idim_max;

    for(i=0;i<X->Def.NCisAjt;i++){
        org_isite1 = X->Def.CisAjt[i][0]+1;
        org_isite2 = X->Def.CisAjt[i][2]+1;
        org_sigma1 = X->Def.CisAjt[i][1];
        org_sigma2 = X->Def.CisAjt[i][3];
        dam_pr=0.0;

        if(org_isite1 == org_isite2){
            if(org_isite1 > X->Def.Nsite){
                if(org_sigma1==org_sigma2){  // longitudinal magnetic field
                    dam_pr += child_GC_CisAis_spin_MPIdouble(org_isite1-1, org_sigma1, 1.0, X, vec, vec);
                }
                else{  // transverse magnetic field
                    dam_pr += child_GC_CisAit_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, 1.0, X, vec, vec);
                }
            }else{
                isite1 = X->Def.Tpow[org_isite1-1];

                if(org_sigma1==org_sigma2){
                    // longitudinal magnetic field
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn) firstprivate(i_max, isite1, org_sigma1, X) shared(vec)
                    for(j=1;j<=i_max;j++){
                        dam_pr += child_SpinGC_CisAis(j, X, isite1, org_sigma1)*conj(vec[j])*vec[j];
                    }
                }else{
                    // transverse magnetic field
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, tmp_off) firstprivate(i_max, isite1, org_sigma2, X) shared(vec)
                    for(j=1;j<=i_max;j++){
                        tmp_sgn  =  child_SpinGC_CisAit(j,X, isite1,org_sigma2,&tmp_off);
                        if(tmp_sgn !=0){
                            dam_pr  +=  tmp_sgn*conj(vec[tmp_off+1])*vec[j];
                        }
                    }
                }
            }
        }else{
            // hopping is not allowed in localized spin system
            dam_pr=0.0;
        }

        dam_pr = SumMPI_dc(dam_pr);
        fprintf(*_fp," %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1, org_sigma1, org_isite2-1, org_sigma2,creal(dam_pr),cimag(dam_pr));
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
int expec_cisajs_SpinGCGeneral(struct BindStruct *X, double complex *vec, FILE **_fp) {
    long unsigned int i, j;
    long unsigned int org_isite1, org_isite2, org_sigma1, org_sigma2;
    double complex dam_pr = 0;
    long int i_max;
    long unsigned int tmp_off = 0;
    int num1;

    i_max = X->Check.idim_max;

    for (i = 0; i < X->Def.NCisAjt; i++) {
        org_isite1 = X->Def.CisAjt[i][0] + 1;
        org_isite2 = X->Def.CisAjt[i][2] + 1;
        org_sigma1 = X->Def.CisAjt[i][1];
        org_sigma2 = X->Def.CisAjt[i][3];
        if (org_isite1 == org_isite2) {
            if (org_isite1 > X->Def.Nsite) {
                if (org_sigma1 == org_sigma2) {
// longitudinal magnetic field
                    dam_pr = child_GC_CisAis_GeneralSpin_MPIdouble(org_isite1 - 1, org_sigma1, 1.0, X, vec, vec);
                } else {
// transverse magnetic field
                    dam_pr = child_GC_CisAit_GeneralSpin_MPIdouble(org_isite1 - 1, org_sigma1, org_sigma2, 1.0, X,
                                                                     vec, vec);
                }
            } else {//org_isite1 <= X->Def.Nsite
                if (org_sigma1 == org_sigma2) {
// longitudinal magnetic field
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, X) shared(vec)
                    for (j = 1; j <= i_max; j++) {
                        num1 = BitCheckGeneral(j - 1, org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
                        dam_pr += conj(vec[j]) * vec[j] * num1;
                    }
                } else {
// transverse magnetic field
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, org_sigma2, X,tmp_off) shared(vec)
                    for (j = 1; j <= i_max; j++) {
                        num1 = GetOffCompGeneralSpin(j - 1, org_isite1, org_sigma2, org_sigma1, &tmp_off,
                                                     X->Def.SiteToBit, X->Def.Tpow);
                        if (num1 != 0) {
                            dam_pr += conj(vec[tmp_off + 1]) * vec[j] * num1;
                        }
                    }
                }
            }
        } else {
// hopping is not allowed in localized spin system
            dam_pr = 0.0;
        }
        dam_pr = SumMPI_dc(dam_pr);
        fprintf(*_fp, " %4ld %4ld %4ld %4ld %.10lf %.10lf\n", org_isite1 - 1, org_sigma1, org_isite2 - 1, org_sigma2,
                creal(dam_pr), cimag(dam_pr));
    }
    return 0;
}
