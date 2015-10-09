/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima */

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
#include "expec_cisajscktaltdc.h"

/** 
 * 
 * 
 * @param X 
 * @param vec 
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return 
 */
int expec_cisajscktaltdc
(
 struct BindStruct *X,
 double complex *vec
 )
{

  FILE *fp;
  char sdt[D_FileNameMax];

  long unsigned int i,j;
  long unsigned int irght,ilft,ihfbit;
  long unsigned int isite1,isite2,isite3,isite4;
  long unsigned int org_isite1,org_isite2,org_isite3,org_isite4,org_sigma1,org_sigma2;
  long unsigned int org_sigma3,org_sigma4;
  long unsigned int isA_up, isB_up;
  long unsigned int Asum,Bsum,Adiff,Bdiff;
  long unsigned int tmp_off=0;
  long unsigned int tmp_off_2=0;
  int tmp_sgn;
  double complex tmp_V;
  double complex dam_pr;
  long int i_max;
  
  //For TPQ
  int step=0;
  int rand_i=0;
  //For Kond
  double complex dmv;
  
  i_max=X->Check.idim_max;
  X->Large.mode=M_CORR;
  tmp_V    = 1.0+0.0*I;
  
  if(GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit)!=0){
    return -1;
  }
 
  dam_pr=0.0;

  //Make File Name for output
  switch (X->Def.iCalcType){
  case Lanczos:
    if(X->Def.St==0){
      sprintf(sdt, cFileName2BGreen_Lanczos, X->Def.CDataFileHead);
      printf("Start: Calculate two bodies Green functions by Lanczos method.\n");
    }else if(X->Def.St==1){
      sprintf(sdt, cFileName2BGreen_CG, X->Def.CDataFileHead);
      printf("Start: Calculate two bodies Green functions by CG method.\n");
    }
    break;

  case TPQCalc:
    step=X->Def.istep;
    rand_i=X->Def.irand;
    TimeKeeperWithRandAndStep(X, cFileNameTimeKeep, cTPQExpecTwoBodyGStart, "a", rand_i, step);
    sprintf(sdt, cFileName2BGreen_TPQ, X->Def.CDataFileHead, rand_i, step);
    break;

  case FullDiag:
    sprintf(sdt, cFileName2BGreen_FullDiag, X->Def.CDataFileHead, X->Phys.eigen_num);
    break;
  }

  if(!childfopen(sdt, "w", &fp)==0){
    return -1;
  }


  switch(X->Def.iCalcModel){
  case HubbardGC:
    for(i=0;i<X->Def.NCisAjtCkuAlvDC;i++){
      org_isite1   = X->Def.CisAjtCkuAlvDC[i][0]+1;
      org_sigma1   = X->Def.CisAjtCkuAlvDC[i][1];
      org_isite2   = X->Def.CisAjtCkuAlvDC[i][2]+1;
      org_sigma2   = X->Def.CisAjtCkuAlvDC[i][3];
      org_isite3   = X->Def.CisAjtCkuAlvDC[i][4]+1;
      org_sigma3   = X->Def.CisAjtCkuAlvDC[i][5];
      org_isite4   = X->Def.CisAjtCkuAlvDC[i][6]+1;
      org_sigma4   = X->Def.CisAjtCkuAlvDC[i][7];
      child_general_int_GetInfo
	(
	 i, 
	 X,
	 org_isite1,
	 org_isite2,
	 org_isite3,
	 org_isite4,
	 org_sigma1,
	 org_sigma2,
	 org_sigma3,
	 org_sigma4,
	 tmp_V
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
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) shared(vec)
        for(j=1;j<=i_max;j++){
	  dam_pr += GC_child_CisAisCisAis_element(j, isite1, isite3, tmp_V, vec, vec, X, &tmp_off);
        }
        fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1, org_isite2-1,org_sigma2, org_isite3-1, org_sigma3, org_isite4-1,org_sigma4, creal(dam_pr), cimag(dam_pr));
      }else if(isite1 == isite2 && isite3 != isite4){
        dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) shared(vec)
        for(j=1;j<=i_max;j++){
	  dam_pr += GC_child_CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, tmp_V, vec, vec, X, &tmp_off);
        }
	fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1, org_isite2-1,org_sigma2, org_isite3-1, org_sigma3, org_isite4-1,org_sigma4, creal(dam_pr), cimag(dam_pr));
      }else if(isite1 != isite2 && isite3 == isite4){
	dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) shared(vec)
	for(j=1;j<=i_max;j++){
	  dam_pr +=GC_child_CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, tmp_V, vec, vec, X, &tmp_off);
	} 
        fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1, org_isite2-1,org_sigma2, org_isite3-1, org_sigma3, org_isite4-1,org_sigma4, creal(dam_pr), cimag(dam_pr));
      }else if(isite1 != isite2 && isite3 != isite4){
	dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) shared(vec)
	for(j=1;j<=i_max;j++){
	  dam_pr +=GC_child_CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, tmp_V, vec, vec, X, &tmp_off_2);   
	}
      } 
      fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1, org_isite2-1,org_sigma2, org_isite3-1, org_sigma3, org_isite4-1,org_sigma4, creal(dam_pr), cimag(dam_pr));
    }
    break;
 
  case KondoGC:
  case Hubbard:
  case Kondo:
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
           
      child_general_int_GetInfo(
				i, 
				X,
				org_isite1,
				org_isite2,
				org_isite3,
				org_isite4,
				org_sigma1,
				org_sigma2,
				org_sigma3,
				org_sigma4,
				tmp_V
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
      dam_pr = 0.0+I*0.0;
      if(isite1 == isite2 && isite3 == isite4){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2) shared(vec,tmp_V)
	for(j=1;j<=i_max;j++){
	  dam_pr += child_CisAisCisAis_element(j, isite1, isite3, tmp_V, vec, vec, X, &tmp_off);
	}
	fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1, org_isite2-1,org_sigma2, org_isite3-1, org_sigma3, org_isite4-1,org_sigma4, creal(dam_pr), cimag(dam_pr));
      }else if(isite1 == isite2 && isite3 != isite4){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2) shared(vec,tmp_V)
	for(j=1;j<=i_max;j++){
	  dam_pr += child_CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, tmp_V, vec, vec, X, &tmp_off);
	}
	fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1, org_isite2-1,org_sigma2, org_isite3-1, org_sigma3, org_isite4-1,org_sigma4, creal(dam_pr), cimag(dam_pr));
      }else if(isite1 != isite2 && isite3 == isite4){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2) shared(vec,tmp_V)
	for(j=1;j<=i_max;j++){
	  dam_pr +=child_CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, tmp_V, vec, vec, X, &tmp_off);
	} 
	fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1, org_isite2-1,org_sigma2, org_isite3-1, org_sigma3, org_isite4-1,org_sigma4, creal(dam_pr), cimag(dam_pr));
      }else if(isite1 != isite2 && isite3 != isite4){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2) shared(vec,tmp_V)
	for(j=1;j<=i_max;j++){
	  dam_pr +=child_CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, tmp_V, vec, vec, X, &tmp_off_2);
	  
	} 
	fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf\n",org_isite1-1,org_sigma1, org_isite2-1,org_sigma2, org_isite3-1, org_sigma3, org_isite4-1,org_sigma4, creal(dam_pr), cimag(dam_pr));
      }    
    }
    break;
  
  case Spin:
    for(i=0;i<X->Def.NCisAjtCkuAlvDC;i++){
      org_isite1   = X->Def.CisAjtCkuAlvDC[i][0]+1;
      org_sigma1   = X->Def.CisAjtCkuAlvDC[i][1];
      org_isite2   = X->Def.CisAjtCkuAlvDC[i][2]+1;
      org_sigma2   = X->Def.CisAjtCkuAlvDC[i][3];
      org_isite3   = X->Def.CisAjtCkuAlvDC[i][4]+1;
      org_sigma3   = X->Def.CisAjtCkuAlvDC[i][5];
      org_isite4   = X->Def.CisAjtCkuAlvDC[i][6]+1;
      org_sigma4   = X->Def.CisAjtCkuAlvDC[i][7];

      if(org_isite1==org_isite2 && org_isite3==org_isite4){
        isA_up = X->Def.Tpow[org_isite2-1];
        isB_up = X->Def.Tpow[org_isite4-1];
        if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
          dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_off_2, tmp_V) shared(vec)
          for(j=1;j<=i_max;j++){
	    dam_pr +=child_CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4, tmp_V, vec, vec, X);
          }
          fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,creal(dam_pr),cimag(dam_pr));
        }else if(org_sigma1==org_sigma4 && org_sigma2==org_sigma3){ // exchange
          dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_off_2,tmp_V) shared(vec)
          for(j=1;j<=i_max;j++){
            tmp_sgn    =  X_child_exchange_spin_element(j,X,isA_up,isB_up,org_sigma2,org_sigma4,&tmp_off);
            dmv        = vec[j]*tmp_sgn;
            dam_pr    += conj(vec[tmp_off])*dmv;
          }
          fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,creal(dam_pr),cimag(dam_pr));
        }else{  // other process is not allowed
          // error message will be added
          fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,0.0,0.0);
        }
      }else if(org_isite1==org_isite4 && org_isite2==org_isite3){
        isA_up = X->Def.Tpow[org_isite1-1];
        isB_up = X->Def.Tpow[org_isite2-1];
        if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4){ // exchange
	  dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_off_2,tmp_V) shared(vec)
	  for(j=1;j<=i_max;j++){
	    tmp_sgn    = X_child_exchange_spin_element(j,X,isA_up,isB_up,org_sigma4,org_sigma2,&tmp_off);
	    dmv        = vec[j]*tmp_sgn;
	    dam_pr    += conj(vec[tmp_off])*dmv;
	  }
	  dam_pr = -1.0*dam_pr;
	  fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,creal(dam_pr),cimag(dam_pr));
        }else{ // this process is not allowed
          //error message will be added 
          fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,0.0,0.0);
        }
      }else{
        //error message will be added 
        fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,0.0,0.0);
      }
    }
    break;

  case SpinGC:
    for(i=0;i<X->Def.NCisAjtCkuAlvDC;i++){
      org_isite1   = X->Def.CisAjtCkuAlvDC[i][0]+1;
      org_sigma1   = X->Def.CisAjtCkuAlvDC[i][1];
      org_isite2   = X->Def.CisAjtCkuAlvDC[i][2]+1;
      org_sigma2   = X->Def.CisAjtCkuAlvDC[i][3];
      org_isite3   = X->Def.CisAjtCkuAlvDC[i][4]+1;
      org_sigma3   = X->Def.CisAjtCkuAlvDC[i][5];
      org_isite4   = X->Def.CisAjtCkuAlvDC[i][6]+1;
      org_sigma4   = X->Def.CisAjtCkuAlvDC[i][7];

      if(org_isite1==org_isite2 && org_isite3==org_isite4){
        isA_up = X->Def.Tpow[org_isite2-1];
        isB_up = X->Def.Tpow[org_isite4-1];
        if(org_sigma1==org_sigma2 && org_sigma3==org_sigma4 ){ //diagonal
          dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_off_2,tmp_V) shared(vec)
          for(j=1;j<=i_max;j++){
	    dam_pr +=GC_child_CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4, tmp_V, vec, vec, X);
          }
          fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,creal(dam_pr),cimag(dam_pr));
        }else if(org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4){ 
          dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_off_2,tmp_V) shared(vec)
          for(j=1;j<=i_max;j++){
	    dam_pr += GC_child_CisAisCitAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec, X, &tmp_off);
          } 
          fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,creal(dam_pr),cimag(dam_pr));
        }else if(org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4){ 
          dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_off_2,tmp_V) shared(vec)
          for(j=1;j<=i_max;j++){
         dam_pr += GC_child_CisAitCiuAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec, X, &tmp_off);
          } 
          fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,creal(dam_pr),cimag(dam_pr));
        }else if(org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4){ 
          dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_off_2,tmp_V) shared(vec)
          for(j=1;j<=i_max;j++){
	    dam_pr += GC_child_CisAitCiuAiv_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, vec, vec, X, &tmp_off);
          }
          fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,creal(dam_pr),cimag(dam_pr));
        }
      }
      else{// hopping process is not allowed
        //error message will be added 
        fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,0.0,0.0);
      }
      if(org_isite1==org_isite4 && org_isite2==org_isite3){ // 
        // 2 <-> 4
        isA_up = X->Def.Tpow[org_isite4-1];
        isB_up = X->Def.Tpow[org_isite2-1];
        if(org_sigma1==org_sigma4 && org_sigma3==org_sigma2 ){ //diagonal
          dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_off_2) shared(vec)
          for(j=1;j<=i_max;j++){
            tmp_sgn      = X_SpinGC_CisAis(j,X,isB_up,org_sigma2);
            tmp_sgn     *= X_SpinGC_CisAis(j,X,isA_up,org_sigma4);
            dmv          = vec[j]*tmp_sgn;
            dam_pr      += conj(vec[j])*dmv;
          }
          dam_pr = -1.0*dam_pr;
          fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,creal(dam_pr),cimag(dam_pr));
        }else if(org_sigma1 == org_sigma4 && org_sigma3 != org_sigma2){ 
          dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_off_2) shared(vec)
          for(j=1;j<=i_max;j++){
            tmp_sgn   =  X_SpinGC_CisAit(j,X,isB_up,org_sigma2,&tmp_off);
            if(tmp_sgn != 0){
              tmp_sgn   *= X_SpinGC_CisAis(j,X,isA_up,org_sigma4);
              dmv        = vec[j]*tmp_sgn;
              dam_pr    += conj(vec[tmp_off]+1)*dmv;
            }
          } 
          dam_pr = -1.0*dam_pr;
          fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,creal(dam_pr),cimag(dam_pr));
        }else if(org_sigma1 != org_sigma4 && org_sigma3 == org_sigma2){ 
          dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_off_2) shared(vec)
          for(j=1;j<=i_max;j++){
            tmp_sgn      = X_SpinGC_CisAis(j,X,isB_up,org_sigma2);
            if(tmp_sgn != 0){
              tmp_sgn   *= X_SpinGC_CisAit(j,X,isA_up,org_sigma4,&tmp_off); 
              dmv        = vec[j]*tmp_sgn;
              dam_pr    += conj(vec[tmp_off]+1)*dmv;
            }
          } 
          dam_pr = -1.0*dam_pr;
          fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,creal(dam_pr),cimag(dam_pr));
        }else if(org_sigma1 != org_sigma4 && org_sigma3 != org_sigma2){ 
          dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_off_2) shared(vec)
          for(j=1;j<=i_max;j++){
            tmp_sgn      =  X_SpinGC_CisAit(j,X,isB_up,org_sigma2,&tmp_off);
            if(tmp_sgn != 0){
              tmp_sgn   *= X_SpinGC_CisAit(tmp_off+1,X,isA_up,org_sigma4,&tmp_off_2); 
              dmv        = vec[j]*tmp_sgn;
              dam_pr    += conj(vec[tmp_off_2]+1)*dmv;
            }
          }
          dam_pr = -1.0*dam_pr;
          fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,creal(dam_pr),cimag(dam_pr));
        }
      }else{// hopping process is not allowed
        //error message will be added 
        fprintf(fp," %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n",org_isite1-1,org_sigma1,org_isite2-1,org_sigma2,org_isite3-1,org_sigma3,org_isite4-1,org_sigma4,0.0,0.0);
      }
      
    }
	break;
  default:
    return -1;
  }
  
  fclose(fp);
  
  if(X->Def.iCalcType==Lanczos){
    if(X->Def.St==0){
      TimeKeeper(X, cFileNameTimeKeep, cLanczosExpecTwoBodyGFinish,"a");
      printf("%s", cLogLanczosExpecTwoBodyGFinish);
    }else if(X->Def.St==1){
      TimeKeeper(X, cFileNameTimeKeep, cCGExpecTwoBodyGFinish,"a");
      printf("%s", cLogCGExpecTwoBodyGFinish);
    }
  }
  else if(X->Def.iCalcType==TPQCalc){
    TimeKeeperWithRandAndStep(X, cFileNameTimeKeep, cTPQExpecTwoBodyGFinish, "a", rand_i, step);
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


/** 
 * 
 * 
 * @param X 
 * @param vec 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void expec_cisajscktaltdc_alldiag(struct BindStruct *X,double complex *vec){ // only for Spin and Hubbard

  long unsigned int j;
  long unsigned int irght,ilft,ihfbit;
  long unsigned int isite1,isite2;
  long unsigned int is1_up,is2_up,is1_down,is2_down;
  long unsigned int iexchg, off;
  int num1_up,num2_up;
  int num1_down,num2_down; 
  long unsigned int ibit1_up,ibit2_up,ibit1_down,ibit2_down; 
  double complex spn_z,chrg_z;
  double complex spn,chrg;
  double complex t_spn_z;  
  
  int N2,i_max;
    
  N2=2*X->Def.Nsite;
  i_max=X->Check.idim_max;
    
  irght=pow(2,((N2+1)/2))-1;
  ilft=pow(2,N2)-1;
  ilft=ilft ^ irght; 
  ihfbit=pow(2,((N2+1)/2));
  chrg=0.0;
  spn=0.0;
  t_spn_z=0.0;
  
  for(isite1=1;isite1<=X->Def.Nsite;isite1++){
    for(isite2=1;isite2<=X->Def.Nsite;isite2++){
            
      is1_up=X->Def.Tpow[2*isite1-2];
      is1_down=X->Def.Tpow[2*isite1-1];
      is2_up=X->Def.Tpow[2*isite2-2];
      is2_down=X->Def.Tpow[2*isite2-1];
      
#pragma omp parallel for reduction(+: t_spn_z,spn, chrg) firstprivate(i_max, is1_up, is2_up, is1_down, is2_down, irght, ilft, ihfbit) private(ibit1_up, num1_up, ibit2_up, num2_up, ibit1_down, num1_down, ibit2_down, num2_down, spn_z, chrg_z, iexchg, off)
      for(j=1;j<=i_max;j++){                    
	ibit1_up= list_1[j]&is1_up;
	num1_up=ibit1_up/is1_up;            
	ibit2_up= list_1[j]&is2_up;
	num2_up=ibit2_up/is2_up;
            
	ibit1_down= list_1[j]&is1_down;
	num1_down=ibit1_down/is1_down;
            
	ibit2_down= list_1[j]&is2_down;
	num2_down=ibit2_down/is2_down;
            
	spn_z=(num1_up-num1_down)*(num2_up-num2_down);
	chrg_z=(num1_up+num1_down)*(num2_up+num2_down);
            
	t_spn_z+= conj(vec[j])* vec[j]*(num1_up-num1_down);
	spn+= conj(vec[j])* vec[j]*spn_z;
	chrg+= conj(vec[j])* vec[j]*chrg_z;
            
	if(isite1==isite2){
	  spn+=2* conj(vec[j])* vec[j]*(num1_up+num1_down-2*num1_up*num1_down);
	}else{
	  if(ibit1_up!=0 && ibit1_down==0 && ibit2_up==0 &&ibit2_down!=0 ){
	    iexchg= list_1[j]-(is1_up+is2_down);
	    iexchg+=(is2_up+is1_down);
	    GetOffComp(list_2_1, list_2_2, iexchg,  irght, ilft, ihfbit, &off);                    
	    spn+=2* conj(vec[j])* vec[off];
	  }else if(ibit1_up==0 && ibit1_down!=0 && ibit2_up!=0 && ibit2_down==0){
	    iexchg= list_1[j]-(is1_down+is2_up);
	    iexchg+=(is2_down+is1_up);
	    GetOffComp(list_2_1, list_2_2, iexchg,  irght, ilft, ihfbit, &off);
	    spn+=2* conj(vec[j])* vec[off];
	  }
	}
      }
    }
  }
  spn = spn/X->Def.Nsite;
  X->Phys.s2=spn;
}

