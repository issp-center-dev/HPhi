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

int expec_cisajscktalt_HubbardGC(struct BindStruct *X,double complex *vec, FILE **_fp);
int expec_cisajscktalt_Hubbard(struct BindStruct *X,double complex *vec, FILE **_fp);

int expec_cisajscktalt_Spin(struct BindStruct *X,double complex *vec, FILE **_fp);
int expec_cisajscktalt_SpinHalf(struct BindStruct *X,double complex *vec, FILE **_fp);
int expec_cisajscktalt_SpinGeneral(struct BindStruct *X,double complex *vec, FILE **_fp);

int expec_cisajscktalt_SpinGC(struct BindStruct *X,double complex *vec, FILE **_fp,FILE **_fp_2,FILE **_fp_3,FILE **_fp_4);
int expec_cisajscktalt_SpinGCHalf(struct BindStruct *X,double complex *vec, FILE **_fp);
int expec_cisajscktalt_SpinGCGeneral(struct BindStruct *X,double complex *vec, FILE **_fp);

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
);
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
 double complex *vec
 )
{

  FILE *fp,*fp_2,*fp_3,*fp_4;
  char sdt[D_FileNameMax],sdt_2[D_FileNameMax],sdt_3[D_FileNameMax],sdt_4[D_FileNameMax],*tmp_char;
  long unsigned int irght,ilft,ihfbit;

  //For TPQ
  int step=0;
  int rand_i=0;

  if(X->Def.NCisAjtCkuAlvDC < 1 && X->Def.NTBody < 1 && X->Def.NFBody < 1 && X->Def.NSBody < 1) return 0;
  X->Large.mode=M_CORR;

  if(GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit)!=0){
    return -1;
  }

  //Make File Name for output
  switch (X->Def.iCalcType){
  case Lanczos:
    if(X->Def.St==0){
      sprintf(sdt, cFileName2BGreen_Lanczos, X->Def.CDataFileHead);
      sprintf(sdt_2, cFileName3BGreen_Lanczos, X->Def.CDataFileHead);
      sprintf(sdt_3, cFileName4BGreen_Lanczos, X->Def.CDataFileHead);
      sprintf(sdt_4, cFileName6BGreen_Lanczos, X->Def.CDataFileHead);
      TimeKeeper(X, cFileNameTimeKeep, cLanczosExpecTwoBodyGStart,"a");
      fprintf(stdoutMPI, "%s", cLogLanczosExpecTwoBodyGStart);
    }else if(X->Def.St==1){
      sprintf(sdt, cFileName2BGreen_CG, X->Def.CDataFileHead);
      sprintf(sdt_2, cFileName3BGreen_Lanczos, X->Def.CDataFileHead);
      sprintf(sdt_3, cFileName4BGreen_Lanczos, X->Def.CDataFileHead);
      sprintf(sdt_4, cFileName6BGreen_Lanczos, X->Def.CDataFileHead);
        TimeKeeper(X, cFileNameTimeKeep, cCGExpecTwoBodyGStart,"a");
      fprintf(stdoutMPI, "%s", cLogLanczosExpecTwoBodyGStart);
    }
    break;

  case TPQCalc:
  case cTPQ:
    step=X->Def.istep;
    rand_i=X->Def.irand;
    TimeKeeperWithRandAndStep(X, cFileNameTimeKeep, cTPQExpecTwoBodyGStart, "a", rand_i, step);
    sprintf(sdt, cFileName2BGreen_TPQ, X->Def.CDataFileHead, rand_i, step);
    sprintf(sdt_2, cFileName3BGreen_TPQ, X->Def.CDataFileHead, rand_i, step);
    sprintf(sdt_3, cFileName4BGreen_TPQ, X->Def.CDataFileHead, rand_i, step);
    sprintf(sdt_4, cFileName6BGreen_TPQ, X->Def.CDataFileHead, rand_i, step);
    break;

  case TimeEvolution:
    step=X->Def.istep;
    TimeKeeperWithStep(X, cFileNameTimeKeep, cTEExpecTwoBodyGStart, "a", step);
    sprintf(sdt, cFileName2BGreen_TE, X->Def.CDataFileHead, step);
    sprintf(sdt_2, cFileName3BGreen_TE, X->Def.CDataFileHead, step);
    sprintf(sdt_3, cFileName4BGreen_TE, X->Def.CDataFileHead, step);
    sprintf(sdt_4, cFileName6BGreen_TE, X->Def.CDataFileHead, step);
    break;

  case FullDiag:
  case CG:
    sprintf(sdt,  cFileName2BGreen_FullDiag, X->Def.CDataFileHead, X->Phys.eigen_num);
    sprintf(sdt_2,cFileName3BGreen_FullDiag, X->Def.CDataFileHead, X->Phys.eigen_num);
    sprintf(sdt_3,cFileName4BGreen_FullDiag, X->Def.CDataFileHead, X->Phys.eigen_num);
    sprintf(sdt_4,cFileName6BGreen_FullDiag, X->Def.CDataFileHead, X->Phys.eigen_num);
    break;
  }

  if(childfopenMPI(sdt, "w", &fp)!=0){
    return -1;
  }
  if(X->Def.NTBody>0){
    if(childfopenMPI(sdt_2, "w", &fp_2)!=0){
      return -1;
    }
  }else{
    fp_2 = fp;
  }
  if(X->Def.NFBody>0){
    if(childfopenMPI(sdt_3, "w", &fp_3)!=0){
      return -1;
    }
  }else{
    fp_3 = fp;
  }
  if(X->Def.NSBody>0){
    if(childfopenMPI(sdt_4, "w", &fp_4)!=0){
      return -1;
    }
  }else{
    fp_4 = fp;
  }



  switch(X->Def.iCalcModel){
  case HubbardGC:
      if(expec_cisajscktalt_HubbardGC(X, vec, &fp)!=0){
          return -1;
      }
    break;
 
  case KondoGC:
  case Hubbard:
  case Kondo:
      if(expec_cisajscktalt_Hubbard(X, vec, &fp)!=0){
          return -1;
      }
    break;
  
  case Spin:
      if(expec_cisajscktalt_Spin(X, vec, &fp)!=0){
          return -1;
      }
    break;

  case SpinGC:
      if(expec_cisajscktalt_SpinGC(X, vec, &fp,&fp_2,&fp_3,&fp_4)!=0){
          return -1;
      }
    break;

  default:
    return -1;
  }
  
  fclose(fp);
  if(X->Def.NTBody>0){
    fclose(fp_2);
  }
  if(X->Def.NFBody>0){
    fclose(fp_3);
  }
  if(X->Def.NSBody>0){
    fclose(fp_4);
  }
  
  if(X->Def.iCalcType==Lanczos){
    if(X->Def.St==0){
      TimeKeeper(X, cFileNameTimeKeep, cLanczosExpecTwoBodyGFinish,"a");
      fprintf(stdoutMPI, "%s", cLogLanczosExpecTwoBodyGFinish);
    }else if(X->Def.St==1){
      TimeKeeper(X, cFileNameTimeKeep, cCGExpecTwoBodyGFinish,"a");
      fprintf(stdoutMPI, "%s", cLogCGExpecTwoBodyGFinish);
    }
  }
  else if(X->Def.iCalcType==TPQCalc){
    TimeKeeperWithRandAndStep(X, cFileNameTimeKeep, cTPQExpecTwoBodyGFinish, "a", rand_i, step);
  }
  else if(X->Def.iCalcType==TimeEvolution){
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
  return 0;
}

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
int expec_cisajscktalt_HubbardGC(struct BindStruct *X,double complex *vec, FILE **_fp){
    long unsigned int i,j;
    long unsigned int isite1,isite2,isite3,isite4;
    long unsigned int org_isite1,org_isite2,org_isite3,org_isite4;
    long unsigned int org_sigma1,org_sigma2,org_sigma3,org_sigma4;
    long unsigned int Asum,Bsum,Adiff,Bdiff;
    long unsigned int tmp_off=0;
    long unsigned int tmp_off_2=0;
    double complex tmp_V= 1.0+0.0*I;
    ;
    double complex dam_pr;
    long int i_max;

    for(i=0;i<X->Def.NCisAjtCkuAlvDC;i++){
        org_isite1   = X->Def.CisAjtCkuAlvDC[i][0]+1;
        org_sigma1   = X->Def.CisAjtCkuAlvDC[i][1];
        org_isite2   = X->Def.CisAjtCkuAlvDC[i][2]+1;
        org_sigma2   = X->Def.CisAjtCkuAlvDC[i][3];
        org_isite3   = X->Def.CisAjtCkuAlvDC[i][4]+1;
        org_sigma3   = X->Def.CisAjtCkuAlvDC[i][5];
        org_isite4   = X->Def.CisAjtCkuAlvDC[i][6]+1;
        org_sigma4   = X->Def.CisAjtCkuAlvDC[i][7];
        dam_pr=0.0;

        if(CheckPE(org_isite1-1, X)==TRUE || CheckPE(org_isite2-1, X)==TRUE ||
           CheckPE(org_isite3-1, X)==TRUE || CheckPE(org_isite4-1, X)==TRUE){
            isite1 = X->Def.OrgTpow[2*org_isite1-2+org_sigma1] ;
            isite2 = X->Def.OrgTpow[2*org_isite2-2+org_sigma2] ;
            isite3 = X->Def.OrgTpow[2*org_isite3-2+org_sigma3] ;
            isite4 = X->Def.OrgTpow[2*org_isite4-2+org_sigma4] ;
            if(isite1 == isite2 && isite3 == isite4){

                dam_pr = child_GC_CisAisCjtAjt_Hubbard_MPI(org_isite1-1, org_sigma1,
                                                             org_isite3-1, org_sigma3,
                                                             1.0, X, vec, vec);
            }
            else if(isite1 == isite2 && isite3 != isite4){

                dam_pr = child_GC_CisAisCjtAku_Hubbard_MPI(org_isite1-1, org_sigma1,
                                                             org_isite3-1, org_sigma3, org_isite4-1, org_sigma4,
                                                             1.0, X, vec, vec);

            }
            else if(isite1 != isite2 && isite3 == isite4){

                dam_pr = child_GC_CisAjtCkuAku_Hubbard_MPI(org_isite1-1, org_sigma1, org_isite2-1, org_sigma2,
                                                             org_isite3-1, org_sigma3,
                                                             1.0, X, vec, vec);

            }
            else if(isite1 != isite2 && isite3 != isite4){
                dam_pr = child_GC_CisAjtCkuAlv_Hubbard_MPI(org_isite1-1, org_sigma1, org_isite2-1, org_sigma2,
                                                             org_isite3-1, org_sigma3, org_isite4-1, org_sigma4,
                                                             1.0, X, vec, vec);
            }

        }//InterPE
        else{
            general_int_GetInfo
                    (i, X, org_isite1, org_isite2, org_isite3, org_isite4,
                     org_sigma1, org_sigma2, org_sigma3, org_sigma4, tmp_V
                    );

            i_max  = X->Large.i_max;
            isite1 = X->Large.is1_spin;
            isite2 = X->Large.is2_spin;
            Asum   = X->Large.isA_spin;
            Adiff  = X->Large.A_spin;

            isite3 = X->Large.is3_spin;
            isite4 = X->Large.is4_spin;
            Bsum   = X->Large.isB_spin;
            Bdiff  = X->Large.B_spin;

            if(isite1 == isite2 && isite3 == isite4){
                dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) shared(vec)
                for(j=1;j<=i_max;j++){
                    dam_pr += GC_CisAisCisAis_element(j, isite1, isite3, tmp_V, vec, vec, X, &tmp_off);
                }
            }else if(isite1 == isite2 && isite3 != isite4){
                dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) shared(vec)
                for(j=1;j<=i_max;j++){
                    dam_pr += GC_CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, tmp_V, vec, vec, X, &tmp_off);
                }
            }else if(isite1 != isite2 && isite3 == isite4){
                dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) shared(vec)
                for(j=1;j<=i_max;j++){
                    dam_pr +=GC_CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, tmp_V, vec, vec, X, &tmp_off);
                }

            }else if(isite1 != isite2 && isite3 != isite4){
                dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) shared(vec)
                for(j=1;j<=i_max;j++){
                    dam_pr +=GC_CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, tmp_V, vec, vec, X, &tmp_off_2);
                }
            }
        }
        dam_pr = SumMPI_dc(dam_pr);
        fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1, org_isite2-1,org_sigma2, org_isite3-1, org_sigma3, org_isite4-1,org_sigma4, creal(dam_pr), cimag(dam_pr));

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
int expec_cisajscktalt_Hubbard(struct BindStruct *X,double complex *vec, FILE **_fp){
    long unsigned int i,j;
    long unsigned int isite1,isite2,isite3,isite4;
    long unsigned int org_isite1,org_isite2,org_isite3,org_isite4;
    long unsigned int org_sigma1,org_sigma2,org_sigma3,org_sigma4;
    long unsigned int Asum,Bsum,Adiff,Bdiff;
    long unsigned int tmp_off=0;
    long unsigned int tmp_off_2=0;
    double complex tmp_V;
    double complex dam_pr;
    long int i_max;

    for(i=0;i<X->Def.NCisAjtCkuAlvDC;i++){
        org_isite1   = X->Def.CisAjtCkuAlvDC[i][0]+1;
        org_sigma1   = X->Def.CisAjtCkuAlvDC[i][1];
        org_isite2   = X->Def.CisAjtCkuAlvDC[i][2]+1;
        org_sigma2   = X->Def.CisAjtCkuAlvDC[i][3];
        org_isite3   = X->Def.CisAjtCkuAlvDC[i][4]+1;
        org_sigma3   = X->Def.CisAjtCkuAlvDC[i][5];
        org_isite4   = X->Def.CisAjtCkuAlvDC[i][6]+1;
        org_sigma4   = X->Def.CisAjtCkuAlvDC[i][7];
        tmp_V    = 1.0;

        dam_pr=0.0;
        if(X->Def.iFlgSzConserved ==TRUE){
            if(org_sigma1+org_sigma3 != org_sigma2+org_sigma4){
                dam_pr=SumMPI_dc(dam_pr);
                fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1, org_sigma1, org_isite2-1, org_sigma2, org_isite3-1, org_sigma3, org_isite4-1, org_sigma4, creal(dam_pr), cimag(dam_pr));
                continue;
            }
        }

        if(CheckPE(org_isite1-1, X)==TRUE || CheckPE(org_isite2-1, X)==TRUE ||
           CheckPE(org_isite3-1, X)==TRUE || CheckPE(org_isite4-1, X)==TRUE){
            isite1 = X->Def.OrgTpow[2*org_isite1-2+org_sigma1] ;
            isite2 = X->Def.OrgTpow[2*org_isite2-2+org_sigma2] ;
            isite3 = X->Def.OrgTpow[2*org_isite3-2+org_sigma3] ;
            isite4 = X->Def.OrgTpow[2*org_isite4-2+org_sigma4] ;
            if(isite1 == isite2 && isite3 == isite4){
                dam_pr = child_CisAisCjtAjt_Hubbard_MPI(org_isite1-1, org_sigma1,
                                                          org_isite3-1, org_sigma3,
                                                          1.0, X, vec, vec);
            }
            else if(isite1 == isite2 && isite3 != isite4){
              //printf("org_isite1=%d, org_isite2=%d, org_isite3=%d, org_isite4=%d\n", org_isite1, org_isite2, org_isite3, org_isite4);
              dam_pr = child_CisAisCjtAku_Hubbard_MPI(org_isite1-1, org_sigma1,
                                                          org_isite3-1, org_sigma3, org_isite4-1, org_sigma4,
                                                          1.0, X, vec, vec);
            }
            else if(isite1 != isite2 && isite3 == isite4){
                dam_pr = child_CisAjtCkuAku_Hubbard_MPI(org_isite1-1, org_sigma1, org_isite2-1, org_sigma2,
                                                          org_isite3-1, org_sigma3,
                                                          1.0, X, vec, vec);

            }
            else if(isite1 != isite2 && isite3 != isite4){
              dam_pr = child_CisAjtCkuAlv_Hubbard_MPI(org_isite1-1, org_sigma1, org_isite2-1, org_sigma2,
                                                          org_isite3-1, org_sigma3, org_isite4-1, org_sigma4,
                                                          1.0, X, vec, vec);
            }

        }//InterPE
        else{
            general_int_GetInfo(
                    i, X, org_isite1, org_isite2, org_isite3, org_isite4,
                    org_sigma1, org_sigma2, org_sigma3, org_sigma4, tmp_V
            );

            i_max  = X->Large.i_max;
            isite1 = X->Large.is1_spin;
            isite2 = X->Large.is2_spin;
            Asum   = X->Large.isA_spin;
            Adiff  = X->Large.A_spin;

            isite3 = X->Large.is3_spin;
            isite4 = X->Large.is4_spin;
            Bsum   = X->Large.isB_spin;
            Bdiff  = X->Large.B_spin;

            tmp_V  = 1.0;
            dam_pr = 0.0;
            if(isite1 == isite2 && isite3 == isite4){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2) shared(vec,tmp_V)
                for(j=1;j<=i_max;j++){
                    dam_pr += CisAisCisAis_element(j, isite1, isite3, tmp_V, vec, vec, X, &tmp_off);
                }
            }else if(isite1 == isite2 && isite3 != isite4){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2) shared(vec,tmp_V)
                for(j=1;j<=i_max;j++){
                    dam_pr += CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, tmp_V, vec, vec, X, &tmp_off);
                }
            }else if(isite1 != isite2 && isite3 == isite4){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2) shared(vec,tmp_V)
                for(j=1;j<=i_max;j++){
                    dam_pr +=CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, tmp_V, vec, vec, X, &tmp_off);
                }
            }else if(isite1 != isite2 && isite3 != isite4){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2) shared(vec,tmp_V)
                for(j=1;j<=i_max;j++){
                    dam_pr +=CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, tmp_V, vec, vec, X, &tmp_off_2);

                }
            }
        }
        dam_pr = SumMPI_dc(dam_pr);
        fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1, org_isite2-1,org_sigma2, org_isite3-1, org_sigma3, org_isite4-1,org_sigma4, creal(dam_pr), cimag(dam_pr));
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
int expec_cisajscktalt_Spin(struct BindStruct *X,double complex *vec, FILE **_fp) {
    int info=0;
    if (X->Def.iFlgGeneralSpin == FALSE) {
        info=expec_cisajscktalt_SpinHalf(X,vec, _fp);
    } else {
        info=expec_cisajscktalt_SpinGeneral(X,vec, _fp);
    }
    return info;
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
int expec_cisajscktalt_SpinHalf(struct BindStruct *X,double complex *vec, FILE **_fp){
    long unsigned int i,j;
    long unsigned int org_isite1,org_isite2,org_isite3,org_isite4;
    long unsigned int org_sigma1,org_sigma2,org_sigma3,org_sigma4;
    long unsigned int tmp_org_isite1,tmp_org_isite2,tmp_org_isite3,tmp_org_isite4;
    long unsigned int tmp_org_sigma1,tmp_org_sigma2,tmp_org_sigma3,tmp_org_sigma4;
    long unsigned int isA_up, isB_up;
    long unsigned int is1_up, is2_up;
    long unsigned int tmp_off=0;
    int tmp_sgn, num1, num2;
    double complex tmp_V;
    double complex dam_pr;
    long int i_max;
    double complex dmv;

    i_max=X->Check.idim_max;
    X->Large.mode=M_CORR;


    for(i=0;i<X->Def.NCisAjtCkuAlvDC;i++){
        tmp_org_isite1   = X->Def.CisAjtCkuAlvDC[i][0]+1;
        tmp_org_sigma1   = X->Def.CisAjtCkuAlvDC[i][1];
        tmp_org_isite2   = X->Def.CisAjtCkuAlvDC[i][2]+1;
        tmp_org_sigma2   = X->Def.CisAjtCkuAlvDC[i][3];
        tmp_org_isite3   = X->Def.CisAjtCkuAlvDC[i][4]+1;
        tmp_org_sigma3   = X->Def.CisAjtCkuAlvDC[i][5];
        tmp_org_isite4   = X->Def.CisAjtCkuAlvDC[i][6]+1;
        tmp_org_sigma4   = X->Def.CisAjtCkuAlvDC[i][7];
        if(Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X,2)!=0){
            //error message will be added
            fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, tmp_org_isite3-1,tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4,0.0,0.0);
            continue;
        }

        dam_pr = 0.0;
        if(org_isite1 >X->Def.Nsite && org_isite3>X->Def.Nsite){
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                is1_up = X->Def.Tpow[org_isite1 - 1];
                is2_up = X->Def.Tpow[org_isite3 - 1];
                num1 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is1_up, org_sigma1);
                num2 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is2_up, org_sigma3);
#pragma omp parallel for default(none) reduction (+:dam_pr) shared(vec) \
  firstprivate(i_max, num1, num2, tmp_V) private(j)
                for (j = 1; j <= i_max; j++) {
                    dam_pr += tmp_V*num1*num2*vec[j]*conj(vec[j]);
                }
            }
            else if(org_isite1==org_isite3 && org_sigma1==org_sigma4 && org_sigma2==org_sigma3){
                is1_up = X->Def.Tpow[org_isite1 - 1];
                num1 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is1_up, org_sigma1);
#pragma omp parallel for default(none) reduction (+:dam_pr) shared(vec) \
  firstprivate(i_max, num1, num2, tmp_V) private(j)
                for (j = 1; j <= i_max; j++) {
                    dam_pr += tmp_V*num1*vec[j]*conj(vec[j]);
                }
            }
            else if(org_sigma1==org_sigma4 && org_sigma2==org_sigma3){//exchange
                dam_pr += child_general_int_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec);
            }
            else{  // other process is not allowed
                // error message will be added
            }
        }
        else if(org_isite1 > X->Def.Nsite || org_isite3>X->Def.Nsite){
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                is1_up = X->Def.Tpow[org_isite1 - 1];
                is2_up = X->Def.Tpow[org_isite3 - 1];
                num2 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is2_up, org_sigma3);
                dam_pr=0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr)shared(vec) \
  firstprivate(i_max, tmp_V, is1_up, org_sigma1, X, num2) private(j, num1)
                for (j = 1; j <= i_max; j++) {
                    num1 = child_Spin_CisAis(j, X, is1_up, org_sigma1);
                    dam_pr += tmp_V*num1*num2*conj(vec[j])*vec[j];
                }
            }
            else if(org_sigma1==org_sigma4 && org_sigma2==org_sigma3){//exchange
                dam_pr += child_general_int_spin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec);
            }
            else{  // other process is not allowed
                // error message will be added
                dam_pr=0.0;
            }
        }
        else{
            isA_up = X->Def.Tpow[org_isite1-1];
            isB_up = X->Def.Tpow[org_isite3-1];
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off, tmp_V) shared(vec)
                for(j=1;j<=i_max;j++){
                    dam_pr +=CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4, tmp_V, vec, vec, X);
                }
            }else if(org_isite1==org_isite3 && org_sigma1==org_sigma4 && org_sigma3==org_sigma2){
                dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) firstprivate(i_max,X,isA_up,org_sigma1, tmp_V) shared(vec, list_1)
                for(j=1;j<=i_max;j++){
                    dmv = child_Spin_CisAis(j, X, isA_up, org_sigma1);
                    dam_pr += vec[j]*tmp_V*dmv*conj(vec[j]);
                }
            }
            else if(org_sigma1==org_sigma4 && org_sigma2==org_sigma3){ // exchange
                dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec)
                for(j=1;j<=i_max;j++){
                    tmp_sgn    =  child_exchange_spin_element(j,X,isA_up,isB_up,org_sigma2,org_sigma4,&tmp_off);
                    dmv        = vec[j]*tmp_sgn;
                    dam_pr    += conj(vec[tmp_off])*dmv;
                }
            }else{  // other process is not allowed
                // error message will be added
                dam_pr=0.0;
            }
        }
        dam_pr = SumMPI_dc(dam_pr);
        fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, tmp_org_isite3-1, tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4,creal(dam_pr),cimag(dam_pr));

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
int expec_cisajscktalt_SpinGeneral(struct BindStruct *X,double complex *vec, FILE **_fp){
    long unsigned int i,j;
    long unsigned int org_isite1,org_isite2,org_isite3,org_isite4;
    long unsigned int org_sigma1,org_sigma2,org_sigma3,org_sigma4;
    long unsigned int tmp_org_isite1,tmp_org_isite2,tmp_org_isite3,tmp_org_isite4;
    long unsigned int tmp_org_sigma1,tmp_org_sigma2,tmp_org_sigma3,tmp_org_sigma4;
    long unsigned int tmp_off=0;
    long unsigned int tmp_off_2=0;
    long unsigned int list1_off=0;
    int num1;
    double complex tmp_V;
    double complex dam_pr;
    long int i_max;
    int tmp_Sz;
    long unsigned int tmp_org=0;
    vec[0]=0;
    i_max=X->Check.idim_max;
    X->Large.mode=M_CORR;

    for(i=0;i<X->Def.NCisAjtCkuAlvDC;i++){
        tmp_org_isite1   = X->Def.CisAjtCkuAlvDC[i][0]+1;
        tmp_org_sigma1   = X->Def.CisAjtCkuAlvDC[i][1];
        tmp_org_isite2   = X->Def.CisAjtCkuAlvDC[i][2]+1;
        tmp_org_sigma2   = X->Def.CisAjtCkuAlvDC[i][3];
        tmp_org_isite3   = X->Def.CisAjtCkuAlvDC[i][4]+1;
        tmp_org_sigma3   = X->Def.CisAjtCkuAlvDC[i][5];
        tmp_org_isite4   = X->Def.CisAjtCkuAlvDC[i][6]+1;
        tmp_org_sigma4   = X->Def.CisAjtCkuAlvDC[i][7];

        if(Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X,2)!=0){
            fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, tmp_org_isite3-1,tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4,0.0,0.0);
            continue;
        }
        tmp_Sz=0;

      for(j=0;j<2; j++) {
        tmp_org = X->Def.CisAjtCkuAlvDC[i][4*j+1]*X->Def.Tpow[X->Def.CisAjtCkuAlvDC[i][4 * j]];
        tmp_Sz += GetLocal2Sz(X->Def.CisAjtCkuAlvDC[i][4 * j] + 1, tmp_org, X->Def.SiteToBit, X->Def.Tpow);
        tmp_org = X->Def.CisAjtCkuAlvDC[i][4*j+3]*X->Def.Tpow[X->Def.CisAjtCkuAlvDC[i][4 * j+2]];
        tmp_Sz -= GetLocal2Sz(X->Def.CisAjtCkuAlvDC[i][4 * j+2] + 1, tmp_org, X->Def.SiteToBit, X->Def.Tpow);
      }
      if(tmp_Sz !=0){ // not Sz conserved
        fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, tmp_org_isite3-1,tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4,0.0,0.0);
        continue;
      }

        dam_pr = 0.0;
        if(org_isite1 >X->Def.Nsite && org_isite3>X->Def.Nsite){
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr=child_CisAisCjuAju_GeneralSpin_MPIdouble(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, tmp_V, X, vec, vec);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr=child_CisAitCjuAjv_GeneralSpin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4,tmp_V, X, vec, vec);
            }
            else{
                dam_pr=0.0;
            }
        }
        else if(org_isite3 > X->Def.Nsite || org_isite1 > X->Def.Nsite){
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr=child_CisAisCjuAju_GeneralSpin_MPIsingle(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, tmp_V, X, vec, vec);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr=child_CisAitCjuAjv_GeneralSpin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4,tmp_V, X, vec, vec);
            }
            else{
                dam_pr=0.0;
            }
        }
        else{
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) firstprivate(i_max,X,org_isite1, org_sigma1,org_isite3, org_sigma3, tmp_V) shared(vec,list_1)
                for(j=1;j<=i_max;j++){
                    num1=BitCheckGeneral(list_1[j], org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
                    if(num1 != FALSE){
                        num1=BitCheckGeneral(list_1[j], org_isite3, org_sigma3, X->Def.SiteToBit, X->Def.Tpow);
                        if(num1 != FALSE){
                            dam_pr += tmp_V*conj(vec[j])*vec[j];
                        }
                    }
                }
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) firstprivate(i_max,X, org_isite1, org_isite3, org_sigma1, org_sigma2, org_sigma3, org_sigma4, tmp_off, tmp_off_2, list1_off, myrank, tmp_V) shared(vec, list_1)
                for(j=1;j<=i_max;j++){
                  num1 = GetOffCompGeneralSpin(list_1[j], org_isite3, org_sigma4, org_sigma3, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
                  if(num1 != FALSE) {
                    num1 = GetOffCompGeneralSpin(tmp_off, org_isite1, org_sigma2, org_sigma1, &tmp_off_2,
                                                 X->Def.SiteToBit, X->Def.Tpow);
                    if (num1 != FALSE) {
                      ConvertToList1GeneralSpin(tmp_off_2, X->Check.sdim, &list1_off);
                      dam_pr += tmp_V * conj(vec[list1_off]) * vec[j];
                    }
                  }
                }
              //printf("DEBUG: rank=%d, dam_pr=%lf\n", myrank, creal(dam_pr));
            }
            else{
                dam_pr=0.0;
            }
        }
        dam_pr = SumMPI_dc(dam_pr);
        fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, tmp_org_isite3-1, tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4, creal(dam_pr),cimag(dam_pr));
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
int expec_Sixbody_SpinGCHalf(struct BindStruct *X,double complex *vec, FILE **_fp){
    long unsigned int i,j;
    long unsigned int tmp_org_isite1,tmp_org_isite2,tmp_org_isite3,tmp_org_isite4,tmp_org_isite5,tmp_org_isite6,tmp_org_isite7,tmp_org_isite8,tmp_org_isite9,tmp_org_isite10,tmp_org_isite11,tmp_org_isite12;
    long unsigned int tmp_org_sigma1,tmp_org_sigma2,tmp_org_sigma3,tmp_org_sigma4,tmp_org_sigma5,tmp_org_sigma6,tmp_org_sigma7,tmp_org_sigma8,tmp_org_sigma9,tmp_org_sigma10,tmp_org_sigma11,tmp_org_sigma12;
    long unsigned int org_isite1,org_isite2,org_isite3,org_isite4,org_isite5,org_isite6,org_isite7,org_isite8,org_isite9,org_isite10,org_isite11,org_isite12;
    long unsigned int org_sigma1,org_sigma2,org_sigma3,org_sigma4,org_sigma5,org_sigma6,org_sigma7,org_sigma8,org_sigma9,org_sigma10,org_sigma11,org_sigma12;
    long unsigned int isA_up, isB_up;
    long unsigned int tmp_off=0;
    double complex tmp_V;
    double complex dam_pr;
    double complex *vec_pr,*vec_pr_0,*vec_pr_1,*vec_pr_2;

    long int i_max;
    i_max=X->Check.idim_max;

    for(i=0;i<X->Def.NSBody;i++){
        //printf("%d %d \n",i,X->Def.NSBody);
        vec_pr_0 = cd_1d_allocate(i_max + 1);
        vec_pr_1 = cd_1d_allocate(i_max + 1);
        vec_pr_2 = cd_1d_allocate(i_max + 1);
        vec_pr   = cd_1d_allocate(i_max + 1);

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
        mltplyHalfSpinGC_mini(X,tmp_org_isite11-1,tmp_org_sigma11,tmp_org_isite12-1,tmp_org_sigma12,vec_pr_0,vec);
        /* |vec_pr_1>= c9a10|vec_pr_0>*/
        mltplyHalfSpinGC_mini(X,tmp_org_isite9-1,tmp_org_sigma9,tmp_org_isite10-1,tmp_org_sigma10,vec_pr_1,vec_pr_0);
        /* |vec_pr_2>= c7a8|vec_pr_1>*/
        mltplyHalfSpinGC_mini(X,tmp_org_isite7-1,tmp_org_sigma7,tmp_org_isite8-1,tmp_org_sigma8,vec_pr_2,vec_pr_1);
        /* |vec_pr>= c5a6|vec_pr_2>*/
        mltplyHalfSpinGC_mini(X,tmp_org_isite5-1,tmp_org_sigma5,tmp_org_isite6-1,tmp_org_sigma6,vec_pr,vec_pr_2);
        X->Large.mode = H_CORR;

        if(Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X,6)!=0){
            //error message will be added
            fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, tmp_org_isite3-1,tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4,0.0,0.0);
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
                dam_pr += child_GC_CisAisCjuAju_spin_MPIdouble( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, vec, vec_pr);
            }
            else if(org_isite1 ==org_isite3 && org_sigma1 ==org_sigma4 && org_sigma2 ==org_sigma3){ //diagonal (for spin: cuadcdau=cuau)
                dam_pr += child_GC_CisAis_spin_MPIdouble((org_isite1-1), org_sigma1, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIdouble(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec_pr);
            }
        }
        else if(org_isite3>X->Def.Nsite || org_isite1>X->Def.Nsite){ //org_isite3 > Nsite >= org_isite1
           //printf("S-MPI \n");
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr += child_GC_CisAisCjuAju_spin_MPIsingle( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIsingle(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec_pr);
            }
        }
        else{
            if(org_isite1==org_isite2 && org_isite3==org_isite4){
                isA_up = X->Def.Tpow[org_isite2-1];
                isB_up = X->Def.Tpow[org_isite4-1];
                if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j,i) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr +=GC_CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4, tmp_V, vec, vec_pr, X);
                    }
                }else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAisCitAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec_pr, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec_pr, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiv_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec_pr, X, &tmp_off);
                    }
                }
            }
        }
        dam_pr = SumMPI_dc(dam_pr);
        fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld  %4ld %4ld %4ld %4ld  %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",
        tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, 
        tmp_org_isite3-1, tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4,
        tmp_org_isite5-1, tmp_org_sigma5, tmp_org_isite6-1, tmp_org_sigma6,
        tmp_org_isite7-1, tmp_org_sigma7, tmp_org_isite8-1, tmp_org_sigma8,
        tmp_org_isite9-1, tmp_org_sigma9, tmp_org_isite10-1, tmp_org_sigma10,
        tmp_org_isite11-1, tmp_org_sigma11, tmp_org_isite12-1, tmp_org_sigma12,
        creal(dam_pr),cimag(dam_pr));
        free_cd_1d_allocate(vec_pr);
        free_cd_1d_allocate(vec_pr_0);
        free_cd_1d_allocate(vec_pr_1);
        free_cd_1d_allocate(vec_pr_2);
    }
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
int expec_cisajscktalt_SpinGC(struct BindStruct *X,double complex *vec, FILE **_fp,FILE **_fp_2,FILE **_fp_3,FILE **_fp_4){
    int info=0;
    if (X->Def.iFlgGeneralSpin == FALSE) {
        info = expec_cisajscktalt_SpinGCHalf(X,vec, _fp);
        if(X->Def.NTBody>0){
          info = expec_Threebody_SpinGCHalf(X,vec, _fp_2);
        } 
        if(X->Def.NFBody>0){
          info = expec_Fourbody_SpinGCHalf(X,vec,_fp_3);
        }
        if(X->Def.NSBody>0){
          info = expec_Sixbody_SpinGCHalf(X,vec,_fp_4);
        }
    } else {
        info=expec_cisajscktalt_SpinGCGeneral(X,vec, _fp);
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
int expec_Fourbody_SpinGCHalf(struct BindStruct *X,double complex *vec, FILE **_fp){
    long unsigned int i,j;
    long unsigned int tmp_org_isite1,tmp_org_isite2,tmp_org_isite3,tmp_org_isite4,tmp_org_isite5,tmp_org_isite6,tmp_org_isite7,tmp_org_isite8;
    long unsigned int tmp_org_sigma1,tmp_org_sigma2,tmp_org_sigma3,tmp_org_sigma4,tmp_org_sigma5,tmp_org_sigma6,tmp_org_sigma7,tmp_org_sigma8;
    long unsigned int org_isite1,org_isite2,org_isite3,org_isite4,org_isite5,org_isite6,org_isite7,org_isite8;
    long unsigned int org_sigma1,org_sigma2,org_sigma3,org_sigma4,org_sigma5,org_sigma6,org_sigma7,org_sigma8;
    long unsigned int isA_up, isB_up;
    long unsigned int tmp_off=0;
    double complex tmp_V;
    double complex dam_pr;
    double complex *vec_pr,*vec_pr_tmp;

    long int i_max;
    i_max=X->Check.idim_max;

    for(i=0;i<X->Def.NFBody;i++){
        vec_pr     = cd_1d_allocate(i_max + 1);
        vec_pr_tmp = cd_1d_allocate(i_max + 1);

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
        mltplyHalfSpinGC_mini(X,tmp_org_isite7-1,tmp_org_sigma7,tmp_org_isite8-1,tmp_org_sigma8,vec_pr_tmp,vec);
        /* |vec_pr>= c5a6|vec_pr_tmp>*/
        mltplyHalfSpinGC_mini(X,tmp_org_isite5-1,tmp_org_sigma5,tmp_org_isite6-1,tmp_org_sigma6,vec_pr,vec_pr_tmp);
        X->Large.mode = H_CORR;

        if(Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X,4)!=0){
            //error message will be added
            fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, tmp_org_isite3-1,tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4,0.0,0.0);
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
                dam_pr += child_GC_CisAisCjuAju_spin_MPIdouble( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, vec, vec_pr);
            }
            else if(org_isite1 ==org_isite3 && org_sigma1 ==org_sigma4 && org_sigma2 ==org_sigma3){ //diagonal (for spin: cuadcdau=cuau)
                dam_pr += child_GC_CisAis_spin_MPIdouble((org_isite1-1), org_sigma1, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIdouble(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec_pr);
            }
        }
        else if(org_isite3>X->Def.Nsite || org_isite1>X->Def.Nsite){ //org_isite3 > Nsite >= org_isite1
           //printf("S-MPI \n");
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr += child_GC_CisAisCjuAju_spin_MPIsingle( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIsingle(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec_pr);
            }
        }
        else{
            if(org_isite1==org_isite2 && org_isite3==org_isite4){
                isA_up = X->Def.Tpow[org_isite2-1];
                isB_up = X->Def.Tpow[org_isite4-1];
                if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j,i) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr +=GC_CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4, tmp_V, vec, vec_pr, X);
                    }
                }else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAisCitAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec_pr, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec_pr, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiv_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec_pr, X, &tmp_off);
                    }
                }
            }
        }
        dam_pr = SumMPI_dc(dam_pr);
        fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld  %4ld %4ld %4ld %4ld %.10lf %.10lf \n",
        tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, 
        tmp_org_isite3-1, tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4,
        tmp_org_isite5-1, tmp_org_sigma5, tmp_org_isite6-1, tmp_org_sigma6,
        tmp_org_isite7-1, tmp_org_sigma7, tmp_org_isite8-1, tmp_org_sigma8,
        creal(dam_pr),cimag(dam_pr));
        free_cd_1d_allocate(vec_pr);
        free_cd_1d_allocate(vec_pr_tmp);
    }
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
int expec_Threebody_SpinGCHalf(struct BindStruct *X,double complex *vec, FILE **_fp){
    long unsigned int i,j;
    long unsigned int tmp_org_isite1,tmp_org_isite2,tmp_org_isite3,tmp_org_isite4,tmp_org_isite5,tmp_org_isite6;
    long unsigned int tmp_org_sigma1,tmp_org_sigma2,tmp_org_sigma3,tmp_org_sigma4,tmp_org_sigma5,tmp_org_sigma6;
    long unsigned int org_isite1,org_isite2,org_isite3,org_isite4,org_isite5,org_isite6;
    long unsigned int org_sigma1,org_sigma2,org_sigma3,org_sigma4,org_sigma5,org_sigma6;
    long unsigned int isA_up, isB_up;
    long unsigned int tmp_off=0;
    double complex tmp_V;
    double complex dam_pr;
    double complex *vec_pr;

    long int i_max;
    i_max=X->Check.idim_max;

    for(i=0;i<X->Def.NTBody;i++){
        vec_pr = cd_1d_allocate(i_max + 1);
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
        mltplyHalfSpinGC_mini(X,tmp_org_isite5-1,tmp_org_sigma5,tmp_org_isite6-1,tmp_org_sigma6,vec_pr,vec);
        X->Large.mode = H_CORR;

        if(Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X,3)!=0){
            //error message will be added
            fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, tmp_org_isite3-1,tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4,0.0,0.0);
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
                dam_pr += child_GC_CisAisCjuAju_spin_MPIdouble( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, vec, vec_pr);
            }
            else if(org_isite1 ==org_isite3 && org_sigma1 ==org_sigma4 && org_sigma2 ==org_sigma3){ //diagonal (for spin: cuadcdau=cuau)
                dam_pr += child_GC_CisAis_spin_MPIdouble((org_isite1-1), org_sigma1, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIdouble(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec_pr);
            }
        }
        else if(org_isite3>X->Def.Nsite || org_isite1>X->Def.Nsite){ //org_isite3 > Nsite >= org_isite1
           //printf("S-MPI \n");
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr += child_GC_CisAisCjuAju_spin_MPIsingle( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIsingle(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, vec, vec_pr);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec_pr);
            }
        }
        else{
            if(org_isite1==org_isite2 && org_isite3==org_isite4){
                isA_up = X->Def.Tpow[org_isite2-1];
                isB_up = X->Def.Tpow[org_isite4-1];
                if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j,i) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr +=GC_CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4, tmp_V, vec, vec_pr, X);
                    }
                }else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAisCitAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec_pr, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec_pr, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec,vec_pr)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiv_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec_pr, X, &tmp_off);
                    }
                }
            }
        }
        dam_pr = SumMPI_dc(dam_pr);
        fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",
        tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, 
        tmp_org_isite3-1, tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4,
        tmp_org_isite5-1, tmp_org_sigma5, tmp_org_isite6-1, tmp_org_sigma6,
        creal(dam_pr),cimag(dam_pr));
        free_cd_1d_allocate(vec_pr);
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
int expec_cisajscktalt_SpinGCHalf(struct BindStruct *X,double complex *vec, FILE **_fp){
    long unsigned int i,j;
    long unsigned int org_isite1,org_isite2,org_isite3,org_isite4;
    long unsigned int org_sigma1,org_sigma2,org_sigma3,org_sigma4;
    long unsigned int tmp_org_isite1,tmp_org_isite2,tmp_org_isite3,tmp_org_isite4;
    long unsigned int tmp_org_sigma1,tmp_org_sigma2,tmp_org_sigma3,tmp_org_sigma4;
    long unsigned int isA_up, isB_up;
    long unsigned int tmp_off=0;
    double complex tmp_V;
    double complex dam_pr;
    long int i_max;
    i_max=X->Check.idim_max;

    for(i=0;i<X->Def.NCisAjtCkuAlvDC;i++){
        tmp_org_isite1   = X->Def.CisAjtCkuAlvDC[i][0]+1;
        tmp_org_sigma1   = X->Def.CisAjtCkuAlvDC[i][1];
        tmp_org_isite2   = X->Def.CisAjtCkuAlvDC[i][2]+1;
        tmp_org_sigma2   = X->Def.CisAjtCkuAlvDC[i][3];
        tmp_org_isite3   = X->Def.CisAjtCkuAlvDC[i][4]+1;
        tmp_org_sigma3   = X->Def.CisAjtCkuAlvDC[i][5];
        tmp_org_isite4   = X->Def.CisAjtCkuAlvDC[i][6]+1;
        tmp_org_sigma4   = X->Def.CisAjtCkuAlvDC[i][7];

        if(Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X,2)!=0){
            //error message will be added
            fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, tmp_org_isite3-1,tmp_org_sigma3, tmp_org_isite4-1,tmp_org_sigma4,0.0,0.0);
            continue;
        }

        dam_pr=0.0;
        if(org_isite1>X->Def.Nsite && org_isite3>X->Def.Nsite){ //org_isite3 >= org_isite1 > Nsite

            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr += child_GC_CisAisCjuAju_spin_MPIdouble( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, vec, vec);

            }
            else if(org_isite1 ==org_isite3 && org_sigma1 ==org_sigma4 && org_sigma2 ==org_sigma3){ //diagonal (for spin: cuadcdau=cuau)
                dam_pr += child_GC_CisAis_spin_MPIdouble((org_isite1-1), org_sigma1, tmp_V, X, vec, vec);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIdouble(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, vec, vec);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec);
            }
        }
        else if(org_isite3>X->Def.Nsite || org_isite1>X->Def.Nsite){ //org_isite3 > Nsite >= org_isite1
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr += child_GC_CisAisCjuAju_spin_MPIsingle( (org_isite1-1), org_sigma1, (org_isite3-1), org_sigma3, tmp_V, X, vec, vec);

            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr += child_GC_CisAisCjuAjv_spin_MPIsingle(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr += child_GC_CisAitCjuAju_spin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, vec, vec);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr +=  child_GC_CisAitCiuAiv_spin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec);
            }
        }
        else{
            if(org_isite1==org_isite2 && org_isite3==org_isite4){
                isA_up = X->Def.Tpow[org_isite2-1];
                isB_up = X->Def.Tpow[org_isite4-1];
                if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec)
                    for(j=1;j<=i_max;j++){
                        dam_pr +=GC_CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4, tmp_V, vec, vec, X);
                    }
                }else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAisCitAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec, X, &tmp_off);
                    }
                }else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(vec)
                    for(j=1;j<=i_max;j++){
                        dam_pr += GC_CisAitCiuAiv_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec, X, &tmp_off);
                    }
                }
            }
        }
        dam_pr = SumMPI_dc(dam_pr);
        fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, tmp_org_isite3-1, tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4,creal(dam_pr),cimag(dam_pr));
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
int expec_cisajscktalt_SpinGCGeneral(struct BindStruct *X,double complex *vec, FILE **_fp){
    long unsigned int i,j;
    long unsigned int org_isite1,org_isite2,org_isite3,org_isite4;
    long unsigned int org_sigma1,org_sigma2,org_sigma3,org_sigma4;
    long unsigned int tmp_org_isite1,tmp_org_isite2,tmp_org_isite3,tmp_org_isite4;
    long unsigned int tmp_org_sigma1,tmp_org_sigma2,tmp_org_sigma3,tmp_org_sigma4;
    long unsigned int tmp_off=0;
    long unsigned int tmp_off_2=0;
    int  num1;
    double complex tmp_V;
    double complex dam_pr;
    long int i_max;
    i_max=X->Check.idim_max;
    X->Large.mode=M_CORR;


  for(i=0;i<X->Def.NCisAjtCkuAlvDC;i++){
        tmp_org_isite1   = X->Def.CisAjtCkuAlvDC[i][0]+1;
        tmp_org_sigma1   = X->Def.CisAjtCkuAlvDC[i][1];
        tmp_org_isite2   = X->Def.CisAjtCkuAlvDC[i][2]+1;
        tmp_org_sigma2   = X->Def.CisAjtCkuAlvDC[i][3];
        tmp_org_isite3   = X->Def.CisAjtCkuAlvDC[i][4]+1;
        tmp_org_sigma3   = X->Def.CisAjtCkuAlvDC[i][5];
        tmp_org_isite4   = X->Def.CisAjtCkuAlvDC[i][6]+1;
        tmp_org_sigma4   = X->Def.CisAjtCkuAlvDC[i][7];
        dam_pr = 0.0;

      if(Rearray_Interactions(i, &org_isite1, &org_isite2, &org_isite3, &org_isite4, &org_sigma1, &org_sigma2, &org_sigma3, &org_sigma4, &tmp_V, X,2)!=0){
          //error message will be added
          fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, tmp_org_isite3-1,tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4,0.0,0.0);
          continue;
      }

      if(org_isite1 > X->Def.Nsite && org_isite3 > X->Def.Nsite){
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr=child_GC_CisAisCjuAju_GeneralSpin_MPIdouble(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, tmp_V, X, vec, vec);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr=child_GC_CisAisCjuAjv_GeneralSpin_MPIdouble(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr=child_GC_CisAitCjuAju_GeneralSpin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, vec, vec);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr=child_GC_CisAitCjuAjv_GeneralSpin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4,tmp_V, X, vec, vec);
            }
        }
        else if(org_isite3 > X->Def.Nsite || org_isite1 > X->Def.Nsite){
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
                dam_pr=child_GC_CisAisCjuAju_GeneralSpin_MPIsingle(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, tmp_V, X, vec, vec);
            }
            else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr=child_GC_CisAisCjuAjv_GeneralSpin_MPIsingle(org_isite1-1, org_sigma1, org_isite3-1, org_sigma3, org_sigma4, tmp_V, X, vec, vec);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
                dam_pr=child_GC_CisAitCjuAju_GeneralSpin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, tmp_V, X, vec, vec);
            }
            else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
                dam_pr=child_GC_CisAitCjuAjv_GeneralSpin_MPIsingle(org_isite1-1, org_sigma1, org_sigma2, org_isite3-1, org_sigma3, org_sigma4,tmp_V, X, vec, vec);
            }
        }
        else{
            if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) firstprivate(i_max,X,org_isite1, org_sigma1,org_isite3, org_sigma3, tmp_V) shared(vec)
                for(j=1;j<=i_max;j++){
                    num1=BitCheckGeneral(j-1, org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
                    if(num1 != FALSE){
                        num1=BitCheckGeneral(j-1, org_isite3, org_sigma3, X->Def.SiteToBit, X->Def.Tpow);
                        if(num1 != FALSE){
                            dam_pr += tmp_V*conj(vec[j])*vec[j];
                        }
                    }
                }
            }else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) firstprivate(i_max,X, org_isite1, org_isite3, org_sigma1,org_sigma3,org_sigma4, tmp_off, tmp_V) shared(vec)
                for(j=1;j<=i_max;j++){
                    num1 = GetOffCompGeneralSpin(j-1, org_isite3, org_sigma4, org_sigma3, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
                    if(num1 != FALSE){
                        num1=BitCheckGeneral(tmp_off, org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
                        if(num1 != FALSE){
                            dam_pr += tmp_V*conj(vec[tmp_off+1])*vec[j];
                        }
                    }
                }
            }else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) firstprivate(i_max,X, org_isite1, org_isite3, org_sigma1,org_sigma2, org_sigma3, tmp_off, tmp_V) shared(vec)
                for(j=1;j<=i_max;j++){
                    num1 = BitCheckGeneral(j-1, org_isite3, org_sigma3, X->Def.SiteToBit, X->Def.Tpow);
                    if(num1 != FALSE){
                        num1 = GetOffCompGeneralSpin(j-1, org_isite1, org_sigma2, org_sigma1, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
                        if(num1 != FALSE){
                            dam_pr +=  tmp_V*conj(vec[tmp_off+1])*vec[j];
                        }
                    }
                }
            }else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) firstprivate(i_max,X, org_isite1, org_isite3, org_sigma1, org_sigma2, org_sigma3, org_sigma4, tmp_off, tmp_off_2, tmp_V) shared(vec)
                for(j=1;j<=i_max;j++){
                    num1 = GetOffCompGeneralSpin(j-1, org_isite3, org_sigma4, org_sigma3, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
                    if(num1 != FALSE){
                        num1 = GetOffCompGeneralSpin(tmp_off, org_isite1, org_sigma2, org_sigma1, &tmp_off_2, X->Def.SiteToBit, X->Def.Tpow);
                        if(num1 != FALSE){
                            dam_pr +=  tmp_V*conj(vec[tmp_off_2+1])*vec[j];
                        }
                    }

                }
            }
        }
        dam_pr = SumMPI_dc(dam_pr);
        fprintf(*_fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",tmp_org_isite1-1, tmp_org_sigma1, tmp_org_isite2-1, tmp_org_sigma2, tmp_org_isite3-1, tmp_org_sigma3, tmp_org_isite4-1, tmp_org_sigma4, creal(dam_pr),cimag(dam_pr));
    }
    return 0;
}
