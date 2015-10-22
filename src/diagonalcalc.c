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

/**
 * @file   diagonalcalc.c
 * @version 0.2
 * @details modify functions to calculate diagonal components for general spin
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  File to define functions for calculating diagonal components
 * 
 * 
 */

#include <bitcalc.h>
#include "diagonalcalc.h"
#include "mltply.h" 
#include "wrapperMPI.h" 


/** 
 * @fn function for calculating diagonal components
 * 
 * @param X 
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return 
 */
int diagonalcalc
(
 struct BindStruct *X
 ){
    
  FILE *fp;
  long unsigned int i,j;
  long unsigned int isite1,isite2;
  long unsigned int spin;
  double tmp_V;

  /*[s] For InterAll*/
  long unsigned int A_spin,B_spin;
  /*[e] For InterAll*/
  long unsigned int i_max=X->Check.idim_max;

  fprintf(stdoutMPI, "%s", cProStartCalcDiag);
  
#pragma omp parallel for default(none) private(j) shared(list_Diagonal) firstprivate(i_max)
  for(j = 1;j <= i_max; j++){
    list_Diagonal[j]=0.0;
  }
  
  if(X->Def.NCoulombIntra>0){
    if(childfopenMPI(cFileNameCheckCoulombIntra, "w", &fp)!=0){
      return -1;
    }
    for(i = 0; i < X->Def.NCoulombIntra; i++){
      isite1 = X->Def.CoulombIntra[i][0]+1;
      tmp_V  = X->Def.ParaCoulombIntra[i];     
      fprintf(fp,"i=%ld isite1=%ld tmp_V=%lf \n",i,isite1,tmp_V);    
      SetDiagonalCoulombIntra(isite1, tmp_V, X);
    }
    fclose(fp);
  }

  if(X->Def.EDNChemi>0){
    if(childfopenMPI(cFileNameCheckChemi,"w", &fp)!=0){
      return -1;
    }
    for(i = 0; i < X->Def.EDNChemi; i++){
      isite1 = X->Def.EDChemi[i]+1;
      spin   = X->Def.EDSpinChemi[i];
      tmp_V  = -X->Def.EDParaChemi[i];
      fprintf(fp,"i=%ld spin=%ld isite1=%ld tmp_V=%lf \n",i,spin,isite1,tmp_V);
      if(SetDiagonalChemi(isite1, tmp_V,spin,  X) !=0){
	return -1;
      }
    }
    fclose(fp);	
  }
   
  if(X->Def.NCoulombInter>0){
    if(childfopenMPI(cFileNameCheckInterU,"w", &fp)!=0){
      return -1;
    }
    for(i = 0; i < X->Def.NCoulombInter; i++){
      isite1 = X->Def.CoulombInter[i][0]+1;
      isite2 = X->Def.CoulombInter[i][1]+1;
      tmp_V  = X->Def.ParaCoulombInter[i];
      fprintf(fp,"i=%ld isite1=%ld isite2=%ld tmp_V=%lf \n",i,isite1,isite2,tmp_V);
      if(SetDiagonalCoulombInter(isite1, isite2, tmp_V,  X) !=0){
	return -1;
      }
    }
    fclose(fp);   
  }
  if(X->Def.NHundCoupling>0){
    if(childfopenMPI(cFileNameCheckHund,"w", &fp) !=0){
      return -1;
    }
    for(i = 0; i < X->Def.NHundCoupling; i++){
      isite1 = X->Def.HundCoupling[i][0]+1;
      isite2 = X->Def.HundCoupling[i][1]+1;
      tmp_V  = -X->Def.ParaHundCoupling[i];
      if(SetDiagonalHund(isite1, isite2, tmp_V,  X) !=0){
	return -1;
      }
      fprintf(fp,"i=%ld isite1=%ld isite2=%ld tmp_V=%lf \n",i,isite1,isite2,tmp_V);    
    }
    fclose(fp);   
  }

  if(X->Def.NInterAll_Diagonal>0){    
    if(childfopenMPI(cFileNameCheckInterAll,"w", &fp) !=0){
      return -1;
    }
    for(i = 0; i < X->Def.NInterAll_Diagonal; i++){
      isite1=X->Def.InterAll_Diagonal[i][0]+1;
      A_spin=X->Def.InterAll_Diagonal[i][1];
      isite2=X->Def.InterAll_Diagonal[i][2]+1;
      B_spin=X->Def.InterAll_Diagonal[i][3];
      tmp_V =  X->Def.ParaInterAll_Diagonal[i];
      fprintf(fp,"i=%ld isite1=%ld A_spin=%ld isite2=%ld B_spin=%ld tmp_V=%lf \n", i, isite1, A_spin, isite2, B_spin, tmp_V);
      SetDiagonalInterAll(isite1, isite2, A_spin, B_spin, tmp_V, X);
    }      
     fclose(fp);   
    }
  
  TimeKeeper(X, cFileNameTimeKeep, cDiagonalCalcFinish, "w");
  fprintf(stdoutMPI, "%s", cProEndCalcDiag);
  return 0;
}

/** 
 * 
 * 
 * @brief coulombintra interaction are given by \f$ U_i n_ {i \uparrow}n_{i \downarrow} \f$
 * @param isite1  a site number
 * @param dtmp_V A value of coulombintra interaction \f$ U_i \f$
 * @param X Define list to get dimesnion number 
 * @return 
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalCoulombIntra
(
 const long unsigned int isite1,
 double dtmp_V,
 struct BindStruct *X
 ){
  long unsigned int is;
  long unsigned int ibit;
  long unsigned int is1_up, is1_down;

  long unsigned int j;
  long unsigned int i_max=X->Check.idim_max;


  switch (X->Def.iCalcModel){
  case HubbardGC:
    is1_up   = X->Def.Tpow[2*isite1-2];
    is1_down = X->Def.Tpow[2*isite1-1];
    is=is1_up+is1_down;
#pragma omp parallel for default(none) shared(list_Diagonal, list_1) firstprivate(i_max, is, dtmp_V) private(ibit) 
    for(j = 1;j <= i_max;j++){
      ibit=(j-1)&is;
      if(ibit==is){
	list_Diagonal[j]+=dtmp_V;
      }
    }

   break;
  case KondoGC:
  case Hubbard:
  case Kondo:
    is1_up   = X->Def.Tpow[2*isite1-2];
    is1_down = X->Def.Tpow[2*isite1-1];
    is=is1_up+is1_down;
#pragma omp parallel for default(none) shared(list_Diagonal, list_1) firstprivate(i_max, is, dtmp_V) private(ibit) 
    for(j = 1;j <= i_max;j++){
      ibit=list_1[j]&is;
      if(ibit==is){
	list_Diagonal[j]+=dtmp_V;
      }
    }
    break;
    
  case Spin:
  case SpinGC:
    break;

  default:
    fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
    return -1;
    break;
  }
  return 0;
}


/** 
 * 
 * 
 * @param isite1 
 * @param dtmp_V 
 * @param spin 
 * @param X 
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return 
 */
int SetDiagonalChemi
(
 const long unsigned int isite1,
 double dtmp_V,
 long unsigned int spin,
 struct BindStruct *X
 ){
  long unsigned int is1_up;
  long unsigned int ibit1_up;
  long unsigned int num1;

  long unsigned int is1,ibit1;

  long unsigned int j;
  long unsigned int i_max=X->Check.idim_max;
    
  switch (X->Def.iCalcModel){
  case HubbardGC:
    if(spin==0){
      is1   = X->Def.Tpow[2*isite1-2];
    }else{
      is1 = X->Def.Tpow[2*isite1-1];
    }
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1) private(num1, ibit1)
    for(j = 1;j <= i_max;j++){

      ibit1 = (j-1)&is1;
      num1  = ibit1/is1;
      //fprintf(stdoutMPI, "DEBUG: spin=%ld  is1=%ld: isite1=%ld j=%ld num1=%ld \n",spin,is1,isite1,j,num1);
	      
      list_Diagonal[j]+=num1*dtmp_V;
    }
    break;
  case KondoGC:
  case Hubbard:
  case Kondo:
    if(spin==0){
      is1   = X->Def.Tpow[2*isite1-2];
    }else{
      is1 = X->Def.Tpow[2*isite1-1];
    }

#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1) private(num1, ibit1)
    for(j = 1;j <= i_max;j++){

      ibit1 = list_1[j]&is1;
      num1  = ibit1/is1;	      
      list_Diagonal[j]+=num1*dtmp_V;
    }
    break;
    
  case SpinGC:
    is1_up   = X->Def.Tpow[isite1-1];
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, spin) private(num1, ibit1_up)
    for(j = 1;j <= i_max;j++){
      ibit1_up=(((j-1)& is1_up)/is1_up)^(1-spin);
      list_Diagonal[j] += dtmp_V * ibit1_up;
    }
    break;

  case Spin:
    is1_up   = X->Def.Tpow[isite1-1];
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, spin) private(num1, ibit1_up)
    for(j = 1;j <= i_max;j++){
      ibit1_up=((list_1[j]& is1_up)/is1_up)^(1-spin);
      list_Diagonal[j] += dtmp_V * ibit1_up;
      
    }
    break;
  default:
    fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
    return -1;
  }

  return 0;
}

/** 
 * 
 * 
 * @param isite1 
 * @param isite2 
 * @param dtmp_V 
 * @param X 
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return 
 */
int SetDiagonalCoulombInter
(
 const long unsigned int isite1,
 const long unsigned int isite2,
 double dtmp_V,
 struct BindStruct *X
 )
{

  long unsigned int is1_up, is1_down;
  long unsigned int ibit1_up, ibit1_down;
  long unsigned int num1;
  long unsigned int is2_up, is2_down;
  long unsigned int ibit2_up, ibit2_down;
  long unsigned int num2;

  long unsigned int j;
  long unsigned int i_max=X->Check.idim_max;
  switch (X->Def.iCalcModel){

  case HubbardGC: //list_1[j] -> j-1
    is1_up   = X->Def.Tpow[2*isite1-2];
    is1_down = X->Def.Tpow[2*isite1-1];
    is2_up   = X->Def.Tpow[2*isite2-2];
    is2_down = X->Def.Tpow[2*isite2-1];
#pragma omp parallel for default(none) shared( list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is1_down, is2_up, is2_down) private(num1, ibit1_up, ibit1_down, num2, ibit2_up, ibit2_down)
    for(j = 1;j <= i_max;j++){
      num1=0;
      num2=0;              
      ibit1_up=(j-1)&is1_up;
      num1+=ibit1_up/is1_up;
      ibit1_down=(j-1)&is1_down;
      num1+=ibit1_down/is1_down;

      ibit2_up=(j-1)&is2_up;
      num2+=ibit2_up/is2_up;
      ibit2_down=(j-1)&is2_down;
      num2+=ibit2_down/is2_down;
             
      list_Diagonal[j]+=num1*num2*dtmp_V;
    } 
    break;
  case KondoGC:
  case Hubbard:
  case Kondo:
    is1_up   = X->Def.Tpow[2*isite1-2];
    is1_down = X->Def.Tpow[2*isite1-1];
    is2_up   = X->Def.Tpow[2*isite2-2];
    is2_down = X->Def.Tpow[2*isite2-1];

#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is1_down, is2_up, is2_down) private(num1, ibit1_up, ibit1_down, num2, ibit2_up, ibit2_down)
    for(j = 1;j <= i_max;j++){
      num1=0;
      num2=0;              
      ibit1_up=list_1[j]&is1_up;
      num1+=ibit1_up/is1_up;
      ibit1_down=list_1[j]&is1_down;
      num1+=ibit1_down/is1_down;

      ibit2_up=list_1[j]&is2_up;
      num2+=ibit2_up/is2_up;
      ibit2_down=list_1[j]&is2_down;
      num2+=ibit2_down/is2_down;
             
      list_Diagonal[j]+=num1*num2*dtmp_V;
    } 
    break;
    
  case Spin:
  case SpinGC:
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V)
    for(j = 1;j <= i_max; j++){
      list_Diagonal[j] += dtmp_V;
    } 
    break;
  default:
    fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
    return -1;
  }
  return 0;
}

/** 
 * 
 * 
 * @param isite1 
 * @param isite2 
 * @param dtmp_V 
 * @param X 
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return 
 */
int SetDiagonalHund
(
 const long unsigned int isite1,
 const long unsigned int isite2,
 double dtmp_V,
 struct BindStruct *X
 ){

  long unsigned int is1_up, is1_down;
  long unsigned int ibit1_up, ibit1_down;
  long unsigned int num1_up, num1_down;
  long unsigned int is2_up, is2_down;
  long unsigned int ibit2_up, ibit2_down;
  long unsigned int num2_up, num2_down;

  long unsigned int is_up;
  long unsigned int ibit;
  long unsigned int j;
  long unsigned int i_max=X->Check.idim_max;

  switch (X->Def.iCalcModel){
  case HubbardGC: // list_1[j] -> j-1
    is1_up   = X->Def.Tpow[2*isite1-2];
    is1_down = X->Def.Tpow[2*isite1-1];
    is2_up   = X->Def.Tpow[2*isite2-2];
    is2_down = X->Def.Tpow[2*isite2-1];
    
#pragma omp parallel for default(none) shared( list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is1_down, is2_up, is2_down) private(num1_up, num1_down, num2_up, num2_down, ibit1_up, ibit1_down, ibit2_up, ibit2_down)
    for(j = 1; j <= i_max;j++){
      num1_up=0;
      num1_down=0;
      num2_up=0;
      num2_down=0;
        
      ibit1_up=(j-1)&is1_up;
      num1_up=ibit1_up/is1_up;
      ibit1_down=(j-1)&is1_down;
      num1_down=ibit1_down/is1_down;
        
      ibit2_up=(j-1)&is2_up;
      num2_up=ibit2_up/is2_up;
      ibit2_down=(j-1)&is2_down;
      num2_down=ibit2_down/is2_down;
        
      list_Diagonal[j]+=dtmp_V*(num1_up*num2_up+num1_down*num2_down);
    }
    break;
  case KondoGC:
  case Hubbard:
  case Kondo:
    is1_up   = X->Def.Tpow[2*isite1-2];
    is1_down = X->Def.Tpow[2*isite1-1];
    is2_up   = X->Def.Tpow[2*isite2-2];
    is2_down = X->Def.Tpow[2*isite2-1];
    
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is1_down, is2_up, is2_down) private(num1_up, num1_down, num2_up, num2_down, ibit1_up, ibit1_down, ibit2_up, ibit2_down)
    for(j = 1; j <= i_max;j++){
      num1_up=0;
      num1_down=0;
      num2_up=0;
      num2_down=0;
        
      ibit1_up=list_1[j]&is1_up;
      num1_up=ibit1_up/is1_up;
      ibit1_down=list_1[j]&is1_down;
      num1_down=ibit1_down/is1_down;
        
      ibit2_up=list_1[j]&is2_up;
      num2_up=ibit2_up/is2_up;
      ibit2_down=list_1[j]&is2_down;
      num2_down=ibit2_down/is2_down;
        
      list_Diagonal[j]+=dtmp_V*(num1_up*num2_up+num1_down*num2_down);
    }
    break;

  case SpinGC:
    is1_up   = X->Def.Tpow[isite1-1];
    is2_up   = X->Def.Tpow[isite2-1];
    is_up    = is1_up+is2_up;
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is2_up, is_up) private(j, ibit) 
    for(j = 1;j <= i_max;j++){
      ibit = (j-1) & is_up;
      if(ibit == 0 || ibit == is_up){
	list_Diagonal[j]+= dtmp_V;
      }
    }
    break;

  case Spin:
    is1_up   = X->Def.Tpow[isite1-1];
    is2_up   = X->Def.Tpow[isite2-1];
    is_up    = is1_up+is2_up;
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is2_up, is_up) private(j, ibit) 
    for(j = 1;j <= i_max;j++){
      ibit = list_1[j] & is_up;
      if(ibit == 0 || ibit == is_up){
	list_Diagonal[j]+= dtmp_V;
      }
    }
    break;
  default:
    fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
    return -1;
  }
  return 0;
}

/** 
 * 
 * 
 * @param isite1 
 * @param isite2 
 * @param isigma1 
 * @param isigma2 
 * @param dtmp_V 
 * @param X 
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return 
 */
int SetDiagonalInterAll
(
 const long unsigned int isite1,
 const long unsigned int isite2,
 long unsigned int isigma1,
 long unsigned int isigma2,
 double dtmp_V,
 struct BindStruct *X
 )
{
  long unsigned int is1_spin;
  long unsigned int is2_spin;
  long unsigned int is1_up;
  long unsigned int is2_up;

  long unsigned int ibit1_spin;
  long unsigned int ibit2_spin;
  
  long unsigned int num1;
  long unsigned int num2;

  long unsigned int j;
  long unsigned int i_max=X->Check.idim_max;
  
  switch (X->Def.iCalcModel){
  case HubbardGC: //list_1[j] -> j-1
    is1_spin   = X->Def.Tpow[2*isite1-2+isigma1];
    is2_spin   = X->Def.Tpow[2*isite2-2+isigma2];
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V, is1_spin, is2_spin) private(num1, ibit1_spin, num2, ibit2_spin)
    for(j = 1;j <= i_max;j++){
      num1=0;
      num2=0;              
      ibit1_spin=(j-1)&is1_spin;
      num1+=ibit1_spin/is1_spin;
      ibit2_spin=(j-1)&is2_spin;
      num2+=ibit2_spin/is2_spin;
      list_Diagonal[j]+=num1*num2*dtmp_V;
    } 
    break;
  case KondoGC:
  case Hubbard:
  case Kondo:
    is1_spin  = X->Def.Tpow[2*isite1-2+isigma1];
    is2_spin = X->Def.Tpow[2*isite2-2+isigma2];

#pragma omp parallel for default(none) shared(list_Diagonal, list_1) firstprivate(i_max, dtmp_V, is1_spin, is2_spin) private(num1, ibit1_spin, num2, ibit2_spin)
    for(j = 1;j <= i_max;j++){
      num1=0;
      num2=0;              
      ibit1_spin=list_1[j]&is1_spin;
      num1+=ibit1_spin/is1_spin;

      ibit2_spin=list_1[j]&is2_spin;
      num2+=ibit2_spin/is2_spin;             
      list_Diagonal[j]+=num1*num2*dtmp_V;
    } 
    break;  
    
  case Spin:
    is1_up   = X->Def.Tpow[isite1-1];
    is2_up   = X->Def.Tpow[isite2-1];
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is2_up, isigma1, isigma2, X) private(j, num1, num2)
    for(j = 1;j <= i_max; j++){
      num1=X_Spin_CisAis(j, X, is1_up, isigma1);
      num2=X_Spin_CisAis(j, X, is2_up, isigma2);      
      list_Diagonal[j] += num1*num2*dtmp_V;
    }
    break;

 case SpinGC:
   if(X->Def.iFlgGeneralSpin==FALSE){
     is1_up   = X->Def.Tpow[isite1-1];
     is2_up   = X->Def.Tpow[isite2-1];
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is2_up, isigma1, isigma2, X) private(j, num1, num2)
     for(j = 1;j <= i_max; j++){
       num1=X_SpinGC_CisAis(j, X, is1_up, isigma1);
       num2=X_SpinGC_CisAis(j, X, is2_up, isigma2);      
       list_Diagonal[j] += num1*num2*dtmp_V;
     } 
   }
   else{
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is2_up, isigma1, isigma2, X) private(j, num1, num2)
     for(j = 1;j <= i_max; j++){
       num1=BitCheckGeneral (j, isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
       num2=BitCheckGeneral (j, isite2, isigma2, X->Def.SiteToBit, X->Def.Tpow);
       list_Diagonal[j] += num1*num2*dtmp_V;
     }
   }
   break;
    
  default:
    fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
    return -1;
  }
   
  return 0;

}
