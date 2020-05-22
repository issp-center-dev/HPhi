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

/**
 * @file   diagonalcalc.c
 * @version 2.1
 * @details add functions to calculate diagonal components for Time evolution.
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @version 0.2
 * @details modify functions to calculate diagonal components for general spin.
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  Calculate diagonal components, i.e. @f$ H_d |\phi_0> = E_d |\phi_0> @f$.
 * 
 * 
 */

#include <bitcalc.h>
#include "FileIO.h"
#include "diagonalcalc.h"
#include "mltplySpinCore.h"
#include "wrapperMPI.h"


int SetDiagonalTETransfer(
        long unsigned int isite1,
        double dtmp_V,
        long unsigned int spin,
        struct BindStruct *X,
        double complex *tmp_v0,
        double complex *tmp_v1
);

int SetDiagonalTEInterAll(
        long unsigned int isite1,
        long unsigned int isite2,
        long unsigned int isigma1,
        long unsigned int isigma2,
        double dtmp_V,
        struct BindStruct *X,
        double complex *tmp_v0,
        double complex *tmp_v1
);

int SetDiagonalTEChemi(
        long unsigned int isite1,
        long unsigned int spin,
        double dtmp_V,
        struct BindStruct *X,
        double complex *tmp_v0,
        double complex *tmp_v1
);

/**
 *
 *
 * @brief Calculate diagonal components and obtain the list, list_diagonal.
 *
 * @param X [in] Struct to get the information of the diagonal operators.
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @retval -1 fail to calculate diagonal components.
 * @retval 0 succeed to calculate diagonal components.
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
  TimeKeeper(X, cFileNameTimeKeep, cDiagonalCalcStart, "a");

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
  
  TimeKeeper(X, cFileNameTimeKeep, cDiagonalCalcFinish, "a");
  fprintf(stdoutMPI, "%s", cProEndCalcDiag);
  return 0;
}

///  @fn diagonalcalcForTE() Update the vector for diagonal operators ( using in Time Evolution mode).
/// \param X [in] Struct to get the information of the diagonal operators.
/// \param tmp_v0 [in,out] Result vector
/// \param tmp_v1 [in] Input produced vector
/// \retval -1 fail to update the vector.
/// \retval  0 succeed to update the vector.
/// \version 2.1
int diagonalcalcForTE
        (
                const int _istep,
                struct BindStruct *X,
                double complex *tmp_v0,
                double complex *tmp_v1
        ) {

  long unsigned int i;
  long unsigned int isite1, isite2;
  long unsigned int A_spin, B_spin;
  double tmp_V;

  if (X->Def.NTETransferDiagonal[_istep] > 0) {
    for (i = 0; i < X->Def.NTETransferDiagonal[_istep]; i++) {
      isite1 = X->Def.TETransferDiagonal[_istep][i][0] + 1;
      A_spin = X->Def.TETransferDiagonal[_istep][i][1];
      tmp_V = -X->Def.ParaTETransferDiagonal[_istep][i];
      SetDiagonalTETransfer(isite1, tmp_V, A_spin, X, tmp_v0, tmp_v1);
    }
  }
  else if (X->Def.NTEInterAllDiagonal[_istep] >0) {
    for (i = 0; i < X->Def.NTEInterAllDiagonal[_istep]; i++) {
      //Assume n_{1\sigma_1} n_{2\sigma_2}
      isite1 = X->Def.TEInterAllDiagonal[_istep][i][0] + 1;
      A_spin = X->Def.TEInterAllDiagonal[_istep][i][1];
      isite2 = X->Def.TEInterAllDiagonal[_istep][i][2] + 1;
      B_spin = X->Def.TEInterAllDiagonal[_istep][i][3];
      tmp_V = X->Def.ParaTEInterAllDiagonal[_istep][i];

      if (SetDiagonalTEInterAll(isite1, isite2, A_spin, B_spin, tmp_V, X, tmp_v0, tmp_v1) != 0) {
        return -1;
      }
    }

    if (X->Def.NTEChemi[_istep] > 0) {
      for(i=0; i< X->Def.NTEChemi[_istep]; i++) {
        isite1 = X->Def.TEChemi[_istep][i] + 1;
        A_spin = X->Def.SpinTEChemi[_istep][i];
        tmp_V = -X->Def.ParaTEChemi[_istep][i];
        if (SetDiagonalTEChemi(isite1, A_spin, tmp_V, X, tmp_v0, tmp_v1) != 0) {
          return -1;
        }
      }
    }
  }
  return 0;
}


/**
 * 
 * 
 * @brief Calculate the components for Coulombintra interaction, \f$ U_i n_ {i \uparrow}n_{i \downarrow} \f$
 * @param isite1 [in] a site number
 * @param dtmp_V [in] A value of coulombintra interaction \f$ U_i \f$
 * @param X [in] Define list to get dimension number
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalCoulombIntra
(
 long unsigned int isite1,
 double dtmp_V,
 struct BindStruct *X
 ){
  long unsigned int is;
  long unsigned int ibit;
  long unsigned int is1_up, is1_down;

  long unsigned int j;
  long unsigned int i_max=X->Check.idim_max;

  /*
   When isite1 is in the inter process region
  */
  if (isite1 > X->Def.Nsite){
    
    switch (X->Def.iCalcModel) {

    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:

      is1_up   = X->Def.Tpow[2 * isite1 - 2];
      is1_down = X->Def.Tpow[2 * isite1 - 1];
      is = is1_up + is1_down;
      ibit = (unsigned long int)myrank & is;
      if (ibit == is) {
#pragma omp parallel for default(none) shared(list_Diagonal) \
                                       firstprivate(i_max, dtmp_V) private(j) 
        for (j = 1; j <= i_max; j++) list_Diagonal[j] += dtmp_V;
      }

      break; /*case HubbardGC, KondoGC, Hubbard, Kondo:*/

    case Spin:
    case SpinGC:
      /*
       They do not have the Coulomb term
      */
      break;

    default:
      fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
      return -1;
      //break;

    }/*switch (X->Def.iCalcModel)*/

    return 0;

  }/*if (isite1 >= X->Def.Nsite*/
  else{
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
      //break;
    }
  }
  return 0;
}


/**
 *
 *
 * @brief Calculate the components for the chemical potential \f$ \mu_{i \sigma_1} n_ {i \sigma_1} \f$
 * @param isite1 [in] a site number
 * @param dtmp_V [in] A value of coulombintra interaction \f$ \mu_{i \sigma_1} \f$
 * @param spin [in] Spin index for the chemical potential
 * @param X [in] Define list to get dimension number
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalChemi
(
 long unsigned int isite1,
 double dtmp_V,
 long unsigned int spin,
 struct BindStruct *X
 ){
  long unsigned int is1_up;
  long unsigned int ibit1_up;
  long unsigned int num1;
  long unsigned int isigma1 =spin;
  long unsigned int is1,ibit1;

  long unsigned int j;
  long unsigned int i_max=X->Check.idim_max;

  /*
    When isite1 is in the inter process region
  */
  if (isite1 > X->Def.Nsite){

    switch (X->Def.iCalcModel) {

    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:

      if (spin == 0) {
        is1 = X->Def.Tpow[2 * isite1 - 2];
      }
      else {
        is1 = X->Def.Tpow[2 * isite1 - 1];
      }
      ibit1 = (unsigned long int)myrank & is1;
      num1 = ibit1 / is1;
#pragma omp parallel for default(none) shared(list_Diagonal) \
                     firstprivate(i_max, dtmp_V, num1) private(j)
      for (j = 1; j <= i_max; j++) list_Diagonal[j] += num1*dtmp_V;

      break;/*case HubbardGC, case KondoGC, Hubbard, Kondo:*/

    case SpinGC:
    case Spin:

      if (X->Def.iFlgGeneralSpin == FALSE) {
        is1_up = X->Def.Tpow[isite1 - 1];
        ibit1_up = (((unsigned long int)myrank& is1_up) / is1_up) ^ (1 - spin);
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, ibit1_up) private(j)
        for (j = 1; j <= i_max; j++) list_Diagonal[j] += dtmp_V * ibit1_up;
      } /*if (X->Def.iFlgGeneralSpin == FALSE)*/
      else /*if (X->Def.iFlgGeneralSpin == TRUE)*/ {
        num1 = BitCheckGeneral((unsigned long int)myrank, 
          isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
        if (num1 != 0) {
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V) private(j)
          for (j = 1; j <= i_max; j++) list_Diagonal[j] += dtmp_V;
        }/*if (num1 != 0)*/
      }/*if (X->Def.iFlgGeneralSpin == TRUE)*/
      break;/*case SpinGC, Spin:*/

    default:
      fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
      return -1;

    } /*switch (X->Def.iCalcModel)*/

    return 0;

  }/*if (isite1 >= X->Def.Nsite*/

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
    if(X->Def.iFlgGeneralSpin==FALSE){
      is1_up   = X->Def.Tpow[isite1-1];
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, spin) private(num1, ibit1_up)
      for(j = 1;j <= i_max;j++){
        ibit1_up=(((j-1)& is1_up)/is1_up)^(1-spin);
        list_Diagonal[j] += dtmp_V * ibit1_up;
      }
    }
    else{
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
     for(j = 1;j <= i_max; j++){
       num1=BitCheckGeneral (j-1, isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
       if(num1 != 0){
         list_Diagonal[j] += dtmp_V;
       }
     }
    }
    break;

  case Spin:
    if(X->Def.iFlgGeneralSpin==FALSE){
      is1_up   = X->Def.Tpow[isite1-1];
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, spin) private(num1, ibit1_up)
      for(j = 1;j <= i_max;j++){
      ibit1_up=((list_1[j]& is1_up)/is1_up)^(1-spin);
      list_Diagonal[j] += dtmp_V * ibit1_up;
      }
    }
    else{
#pragma omp parallel for default(none) shared(list_Diagonal, list_1) firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1) 
     for(j = 1;j <= i_max; j++){
       num1=BitCheckGeneral (list_1[j], isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
       if(num1 != 0){
         list_Diagonal[j] += dtmp_V;
       }
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
 * @brief Calculate the components for Coulombinter interaction, \f$ V_{ij} n_ {i}n_{j} \f$
 * @param isite1 [in] a site number \f$i \f$
 * @param isite2 [in] a site number \f$j \f$
 * @param dtmp_V [in] A value of coulombinter interaction \f$ V_{ij} \f$
 * @param X [in] Define list to get the operator information.
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalCoulombInter
(
  long unsigned int isite1,
  long unsigned int isite2,
 double dtmp_V,
 struct BindStruct *X
 )
{

  long unsigned int is1_up = 0;
  long unsigned int is1_down = 0;
  long unsigned int ibit1_up = 0;
  long unsigned int ibit1_down = 0;
  long unsigned int num1 = 0;
  long unsigned int is2_up = 0;
  long unsigned int is2_down = 0;
  long unsigned int ibit2_up = 0;
  long unsigned int ibit2_down = 0;
  long unsigned int num2 = 0;

  long unsigned int j;
  long unsigned int i_max=X->Check.idim_max;

  /*
   Force isite1 <= isite2
  */
  if (isite2 < isite1) {
    j = isite2;
    isite2 = isite1;
    isite1 = j;
  }/*if (isite2 < isite1)*/
  /*
    When isite1 & site2 are in the inter process region
  */
  if (/*isite2 => */ isite1 > X->Def.Nsite) {

    switch (X->Def.iCalcModel) {

    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:

      is1_up   = X->Def.Tpow[2 * isite1 - 2];
      is1_down = X->Def.Tpow[2 * isite1 - 1];
      is2_up   = X->Def.Tpow[2 * isite2 - 2];
      is2_down = X->Def.Tpow[2 * isite2 - 1];

      num1 = 0;
      num2 = 0;

      ibit1_up = (unsigned long int)myrank&is1_up;
      num1 += ibit1_up / is1_up;
      ibit1_down = (unsigned long int)myrank&is1_down;
      num1 += ibit1_down / is1_down;

      ibit2_up = (unsigned long int)myrank&is2_up;
      num2 += ibit2_up / is2_up;
      ibit2_down = (unsigned long int)myrank&is2_down;
      num2 += ibit2_down / is2_down;
      
#pragma omp parallel for default(none) shared(list_Diagonal) \
      firstprivate(i_max, dtmp_V, num1, num2) private(j)
      for (j = 1; j <= i_max; j++) list_Diagonal[j] += num1*num2*dtmp_V;

      break;/*case HubbardGC, KondoGC, Hubbard, Kondo:*/

    case Spin:
    case SpinGC:
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V)
      for (j = 1; j <= i_max; j++) {
        list_Diagonal[j] += dtmp_V;
      }
      break;/*case Spin, SpinGC*/

    default:
      fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
      return -1;

    }/*switch (X->Def.iCalcModel)*/

    return 0;

  }/*if (isite1 > X->Def.Nsite)*/
  else if (isite2 > X->Def.Nsite /* => isite1 */) {
    
    switch(X->Def.iCalcModel){
    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:
      is1_up   = X->Def.Tpow[2 * isite1 - 2];
      is1_down = X->Def.Tpow[2 * isite1 - 1];
      is2_up   = X->Def.Tpow[2 * isite2 - 2];
      is2_down = X->Def.Tpow[2 * isite2 - 1];      
      num2 = 0;
      ibit2_up = (unsigned long int)myrank&is2_up;
      num2 += ibit2_up / is2_up;
      ibit2_down = (unsigned long int)myrank&is2_down;
      num2 += ibit2_down / is2_down;
      break;
      
    case Spin:
    case SpinGC:
      break;

    default:
      fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
      return -1;
    }
    
    switch (X->Def.iCalcModel) {

    case HubbardGC:
      
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, num2, is1_up, is1_down) \
private(num1, ibit1_up, ibit1_down, j)
      for (j = 1; j <= i_max; j++) {
        num1 = 0;
        ibit1_up = (j - 1)&is1_up;
        num1 += ibit1_up / is1_up;
        ibit1_down = (j - 1)&is1_down;
        num1 += ibit1_down / is1_down;

        list_Diagonal[j] += num1*num2*dtmp_V;
      }

      break;/*case HubbardGC*/

    case KondoGC:
    case Hubbard:
    case Kondo:
      
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up, is1_down, num2) \
private(num1, ibit1_up, ibit1_down, j)
      for (j = 1; j <= i_max; j++) {
        num1 = 0;
        ibit1_up = list_1[j] & is1_up;
        num1 += ibit1_up / is1_up;
        ibit1_down = list_1[j] & is1_down;
        num1 += ibit1_down / is1_down;

        list_Diagonal[j] += num1*num2*dtmp_V;
      }
      break;/*case KondoGC, Hubbard, Kondo:*/

    case Spin:
    case SpinGC:
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V)
      for (j = 1; j <= i_max; j++) {
        list_Diagonal[j] += dtmp_V;
      }
      break;/* case Spin, SpinGC:*/

    default:
      fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
      return -1;

    }/*switch (X->Def.iCalcModel)*/

    return 0;

  }/*else if (isite2 > X->Def.Nsite)*/
  else{
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
  }
  
  return 0;
}

/**
 *
 * @brief Calculate the components for Hund interaction, \f$ H_{ij}(n_ {i\uparrow}n_{j\uparrow}+ n_ {i\downarrow}n_{j\downarrow})\f$
 * @param isite1 [in] a site number \f$i \f$
 * @param isite2 [in] a site number \f$j \f$
 * @param dtmp_V [in] A value of Hund interaction \f$ H_{ij} \f$
 * @param X [in] Define list to get the operator information.
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalHund
(
 long unsigned int isite1,
 long unsigned int isite2,
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
  /*
  Force isite1 <= isite2
  */
  if (isite2 < isite1) {
    j = isite2;
    isite2 = isite1;
    isite1 = j;
  }
  /*
  When isite1 & site2 are in the inter process region
  */
  if (/*isite2 >= */ isite1 > X->Def.Nsite){

    switch (X->Def.iCalcModel) {

    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:

      is1_up   = X->Def.Tpow[2 * isite1 - 2];
      is1_down = X->Def.Tpow[2 * isite1 - 1];
      is2_up   = X->Def.Tpow[2 * isite2 - 2];
      is2_down = X->Def.Tpow[2 * isite2 - 1];

      num1_up = 0;
      num1_down = 0;
      num2_up = 0;
      num2_down = 0;

      ibit1_up = (unsigned long int)myrank &is1_up;
      num1_up = ibit1_up / is1_up;
      ibit1_down = (unsigned long int)myrank &is1_down;
      num1_down = ibit1_down / is1_down;

      ibit2_up = (unsigned long int)myrank &is2_up;
      num2_up = ibit2_up / is2_up;
      ibit2_down = (unsigned long int)myrank &is2_down;
      num2_down = ibit2_down / is2_down;

#pragma omp parallel for default(none) shared(list_Diagonal) \
  firstprivate(i_max, dtmp_V, num1_up, num1_down, num2_up, num2_down) private(j)
      for (j = 1; j <= i_max; j++)
        list_Diagonal[j] += dtmp_V*(num1_up*num2_up + num1_down*num2_down);

      break;/*case HubbardGC, KondoGC, Hubbard, Kondo:*/

    case SpinGC:
    case Spin:

      is1_up = X->Def.Tpow[isite1 - 1];
      is2_up = X->Def.Tpow[isite2 - 1];
      is_up = is1_up + is2_up;
      ibit = (unsigned long int)myrank & is_up;
      if (ibit == 0 || ibit == is_up) {
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V) private(j) 
        for (j = 1; j <= i_max; j++) list_Diagonal[j] += dtmp_V;
      }
      break;/*case SpinGC, Spin:*/

    default:
      fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
      return -1;
    }

    return 0;

  }/*if (isite1 > X->Def.Nsite)*/
  else if (isite2 > X->Def.Nsite /* >= isite1 */) {

    switch (X->Def.iCalcModel) {

    case HubbardGC:

      is1_up   = X->Def.Tpow[2 * isite1 - 2];
      is1_down = X->Def.Tpow[2 * isite1 - 1];
      is2_up   = X->Def.Tpow[2 * isite2 - 2];
      is2_down = X->Def.Tpow[2 * isite2 - 1];

      num2_up = 0;
      num2_down = 0;

      ibit2_up =  (unsigned long int)myrank &is2_up;
      num2_up = ibit2_up / is2_up;
      ibit2_down =  (unsigned long int)myrank &is2_down;
      num2_down = ibit2_down / is2_down;

#pragma omp parallel for default(none) shared( list_Diagonal) \
firstprivate(i_max, dtmp_V, num2_up, num2_down, is1_up, is1_down) \
private(num1_up, num1_down, ibit1_up, ibit1_down, j)
      for (j = 1; j <= i_max; j++) {
        num1_up = 0;
        num1_down = 0;
 
        ibit1_up = (j - 1)&is1_up;
        num1_up = ibit1_up / is1_up;
        ibit1_down = (j - 1)&is1_down;
        num1_down = ibit1_down / is1_down;

        list_Diagonal[j] += dtmp_V*(num1_up*num2_up + num1_down*num2_down);
      }
      break;/*case HubbardGC:*/

    case KondoGC:
    case Hubbard:
    case Kondo:

      is1_up   = X->Def.Tpow[2 * isite1 - 2];
      is1_down = X->Def.Tpow[2 * isite1 - 1];
      is2_up   = X->Def.Tpow[2 * isite2 - 2];
      is2_down = X->Def.Tpow[2 * isite2 - 1];

      num2_up = 0;
      num2_down = 0;

      ibit2_up = (unsigned long int)myrank&is2_up;
      num2_up = ibit2_up / is2_up;
      ibit2_down = (unsigned long int)myrank&is2_down;
      num2_down = ibit2_down / is2_down;

#pragma omp parallel for default(none) shared(list_1, list_Diagonal) \
firstprivate(i_max, dtmp_V, num2_up, num2_down, is1_up, is1_down) \
private(num1_up, num1_down, ibit1_up, ibit1_down, j)
      for (j = 1; j <= i_max; j++) {
        num1_up = 0;
        num1_down = 0;

        ibit1_up = list_1[j] & is1_up;
        num1_up = ibit1_up / is1_up;
        ibit1_down = list_1[j] & is1_down;
        num1_down = ibit1_down / is1_down;

        list_Diagonal[j] += dtmp_V*(num1_up*num2_up + num1_down*num2_down);
      }
      break;/*case KondoGC, Hubbard, Kondo:*/

    case SpinGC:
      is1_up = X->Def.Tpow[isite1 - 1];
      is2_up = X->Def.Tpow[isite2 - 1];
      ibit2_up = (unsigned long int)myrank & is2_up;

      if (ibit2_up == is2_up) {
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up) private(j, ibit1_up) 
        for (j = 1; j <= i_max; j++) {
          ibit1_up = (j - 1) & is1_up;
          if (ibit1_up == is1_up) {
            list_Diagonal[j] += dtmp_V;
          }
        }
      }
      else if(ibit2_up == 0){
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up) private(j, ibit1_up) 
        for (j = 1; j <= i_max; j++) {
          ibit1_up = (j - 1) & is1_up;
          if (ibit1_up == 0) {
            list_Diagonal[j] += dtmp_V;
          }
        }
      }
      break;/*case SpinGC:*/

    case Spin:
      is1_up = X->Def.Tpow[isite1 - 1];
      is2_up = X->Def.Tpow[isite2 - 1];
      ibit2_up = (unsigned long int)myrank & is2_up;

      if (ibit2_up == is2_up) {
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up) private(j, ibit1_up) 
        for (j = 1; j <= i_max; j++) {
          ibit1_up = list_1[j] & is1_up;
          if (ibit1_up == is1_up) {
            list_Diagonal[j] += dtmp_V;
          }
        }
      }
      else if (ibit2_up == 0) {
#pragma omp parallel for default(none) shared(list_1, list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up) private(j, ibit1_up) 
        for (j = 1; j <= i_max; j++) {
          ibit1_up = list_1[j] & is1_up;
          if (ibit1_up == 0) {
            list_Diagonal[j] += dtmp_V;
          }
        }
      }
      break;/*case Spin:*/

    default:
      fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
      return -1;

    }/*switch (X->Def.iCalcModel)*/

    return 0;

  }/*else if (isite2 > X->Def.Nsite)*/
  else{
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
  }
  return 0;
}

/**
 *
 * @brief Calculate the components for general two-body diagonal interaction, \f$ H_{i\sigma_1 j\sigma_2} n_ {i\sigma_1}n_{j\sigma_2}\f$
 * @param isite1 [in] a site number \f$i \f$
 * @param isite2 [in] a site number \f$j \f$
 * @param isigma1 [in] a spin index at \f$i \f$ site.
 * @param isigma2 [in] a spin index at \f$j \f$ site.
 * @param dtmp_V [in] A value of general two-body diagonal interaction \f$ H_{i\sigma_1 j\sigma_2} \f$
 * @param X [in] Define list to get the operator information.
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalInterAll
(
 long unsigned int isite1,
 long unsigned int isite2,
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
  
  /*
  Forse isite1 <= isite2
  */
  if (isite2 < isite1) {
    j = isite2;
    isite2 = isite1;
    isite1 = j;
    j = isigma2;
    isigma2 = isigma1;
    isigma1 = j;
  }
  /*
  When isite1 & site2 are in the inter process regino
  */
  if (isite1 > X->Def.Nsite) {

    switch (X->Def.iCalcModel) {

    case HubbardGC:
    case KondoGC:
    case Hubbard:
    case Kondo:

      is1_spin = X->Def.Tpow[2 * isite1 - 2 + isigma1];
      is2_spin = X->Def.Tpow[2 * isite2 - 2 + isigma2];

      num1 = 0;
      ibit1_spin = (unsigned long int)myrank&is1_spin;
      num1 += ibit1_spin / is1_spin;

      num2 = 0;
      ibit2_spin = (unsigned long int)myrank&is2_spin;
      num2 += ibit2_spin / is2_spin;

#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, num2, num1) private(ibit1_spin, j)
      for (j = 1; j <= i_max; j++) list_Diagonal[j] += num1*num2*dtmp_V;

      break;/*case HubbardGC, KondoGC, Hubbard, Kondo:*/

    case SpinGC:
    case Spin:

      if (X->Def.iFlgGeneralSpin == FALSE) {
        is1_up = X->Def.Tpow[isite1 - 1];
        is2_up = X->Def.Tpow[isite2 - 1];
        num1 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is1_up, isigma1);
        num2 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is2_up, isigma2);

#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up, isigma1, X, num1, num2) private(j)
        for (j = 1; j <= i_max; j++) {
          list_Diagonal[j] += num1*num2*dtmp_V;
        }
      }/*if (X->Def.iFlgGeneralSpin == FALSE)*/
      else {//start:generalspin
        num1 = BitCheckGeneral((unsigned long int)myrank, isite1, isigma1, 
          X->Def.SiteToBit, X->Def.Tpow);
        num2 = BitCheckGeneral((unsigned long int)myrank, isite2, isigma2,
          X->Def.SiteToBit, X->Def.Tpow);
        if (num1 !=0 && num2 != 0) {
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, num1, X) private(j)
          for (j = 1; j <= i_max; j++) list_Diagonal[j] += dtmp_V*num1;
        }
      }/*if (X->Def.iFlgGeneralSpin == TRUE)*/

      break;/*case SpinGC, Spin:*/

    default:
      fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
      return -1;

    }/*if (isite1 > X->Def.Nsite)*/

    return 0;

  }/*if (isite1 > X->Def.Nsite)*/
  else if (isite2 > X->Def.Nsite) {

    switch (X->Def.iCalcModel) {

    case HubbardGC:

      is1_spin = X->Def.Tpow[2 * isite1 - 2 + isigma1];
      is2_spin = X->Def.Tpow[2 * isite2 - 2 + isigma2];

      num2 = 0;
      ibit2_spin = (unsigned long int)myrank&is2_spin;
      num2 += ibit2_spin / is2_spin;

#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_spin, num2) private(num1, ibit1_spin, j)
      for (j = 1; j <= i_max; j++) {
        num1 = 0;
        ibit1_spin = (j - 1)&is1_spin;
        num1 += ibit1_spin / is1_spin;
        list_Diagonal[j] += num1*num2*dtmp_V;
      }
      break;/*case HubbardGC:*/

    case KondoGC:
    case Hubbard:
    case Kondo:

      is1_spin = X->Def.Tpow[2 * isite1 - 2 + isigma1];
      is2_spin = X->Def.Tpow[2 * isite2 - 2 + isigma2];

      num2 = 0;
      ibit2_spin = (unsigned long int)myrank&is2_spin;
      num2 += ibit2_spin / is2_spin;

#pragma omp parallel for default(none) shared(list_Diagonal, list_1) \
firstprivate(i_max, dtmp_V, is1_spin, num2) private(num1, ibit1_spin, j)
      for (j = 1; j <= i_max; j++) {
        num1 = 0;
        ibit1_spin = list_1[j] & is1_spin;
        num1 += ibit1_spin / is1_spin;
        list_Diagonal[j] += num1*num2*dtmp_V;
      }
      break;/*case KondoGC, Hubbard, Kondo:*/

     case SpinGC:
   
      if (X->Def.iFlgGeneralSpin == FALSE) {
        is1_up = X->Def.Tpow[isite1 - 1];
        is2_up = X->Def.Tpow[isite2 - 1];
        num2 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is2_up, isigma2);

#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up, isigma1, X, num2) private(j, num1)
        for (j = 1; j <= i_max; j++) {
          num1 = child_SpinGC_CisAis(j, X, is1_up, isigma1);
          list_Diagonal[j] += num1*num2*dtmp_V;
        }
      }/* if (X->Def.iFlgGeneralSpin == FALSE)*/
      else {//start:generalspin
        num2 = BitCheckGeneral((unsigned long int)myrank, isite2, isigma2, 
          X->Def.SiteToBit, X->Def.Tpow);
        if (num2 != 0) {
#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
          for (j = 1; j <= i_max; j++) {
            num1 = BitCheckGeneral(j - 1, isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
            list_Diagonal[j] += dtmp_V*num1;
          }
        }
      }/* if (X->Def.iFlgGeneralSpin == TRUE)*/

      break;/*case SpinGC:*/

    case Spin:
   
      if (X->Def.iFlgGeneralSpin == FALSE) {
        is1_up = X->Def.Tpow[isite1 - 1];
        is2_up = X->Def.Tpow[isite2 - 1];
        num2 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is2_up, isigma2);

#pragma omp parallel for default(none) shared(list_Diagonal) \
firstprivate(i_max, dtmp_V, is1_up, isigma1, X, num2) private(j, num1)
        for (j = 1; j <= i_max; j++) {
          num1 = child_Spin_CisAis(j, X, is1_up, isigma1);
          list_Diagonal[j] += num1*num2*dtmp_V;
        }
      }/* if (X->Def.iFlgGeneralSpin == FALSE)*/
      else /* if (X->Def.iFlgGeneralSpin == TRUE)*/{
        num2 = BitCheckGeneral((unsigned long int)myrank, isite2, isigma2, \
          X->Def.SiteToBit, X->Def.Tpow);
        if (num2 != 0) {
#pragma omp parallel for default(none) shared(list_Diagonal, list_1) \
firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
          for (j = 1; j <= i_max; j++) {
            num1 = BitCheckGeneral(list_1[j], isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
            list_Diagonal[j] += dtmp_V*num1;
          }
        }
      } /* if (X->Def.iFlgGeneralSpin == TRUE)*/

      break;/*case Spin:*/

    default:
      fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
      return -1;

    }/*switch (X->Def.iCalcModel)*/

    return 0;

  }/*else if (isite2 > X->Def.Nsite)*/

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
   if(X->Def.iFlgGeneralSpin==FALSE){
     is1_up   = X->Def.Tpow[isite1-1];
     is2_up   = X->Def.Tpow[isite2-1];
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is2_up, isigma1, isigma2, X) private(j, num1, num2)
     for(j = 1;j <= i_max; j++){
       num1=child_Spin_CisAis(j, X, is1_up, isigma1);
       num2=child_Spin_CisAis(j, X, is2_up, isigma2);      
       list_Diagonal[j] += num1*num2*dtmp_V;
     }
   }
   else{
#pragma omp parallel for default(none) shared(list_Diagonal, list_1) firstprivate(i_max, dtmp_V, isite1, isite2, isigma1, isigma2, X) private(j, num1)
     for(j = 1;j <= i_max; j++){
       num1=BitCheckGeneral (list_1[j], isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
       if(num1 != 0){
         num1=BitCheckGeneral (list_1[j], isite2, isigma2, X->Def.SiteToBit, X->Def.Tpow);
         list_Diagonal[j] += dtmp_V*num1;
       }
     }

   }
    break;

 case SpinGC:
   if(X->Def.iFlgGeneralSpin==FALSE){
     is1_up   = X->Def.Tpow[isite1-1];
     is2_up   = X->Def.Tpow[isite2-1];
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V, is1_up, is2_up, isigma1, isigma2, X) private(j, num1, num2)
     for(j = 1;j <= i_max; j++){
       num1=child_SpinGC_CisAis(j, X, is1_up, isigma1);
       num2=child_SpinGC_CisAis(j, X, is2_up, isigma2);      
       list_Diagonal[j] += num1*num2*dtmp_V;
     } 
   }
   else{//start:generalspin
#pragma omp parallel for default(none) shared(list_Diagonal) firstprivate(i_max, dtmp_V, isite1, isite2, isigma1, isigma2, X) private(j, num1)
     for(j = 1;j <= i_max; j++){
       num1=BitCheckGeneral (j-1, isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
       if(num1 != 0){
         num1=BitCheckGeneral (j-1, isite2, isigma2, X->Def.SiteToBit, X->Def.Tpow);
         list_Diagonal[j] += dtmp_V*num1;
       }
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
 * @brief Update the vector by the general two-body diagonal interaction, \f$ H_{i\sigma_1 j\sigma_2} n_ {i\sigma_1}n_{j\sigma_2}\f$.\n
 * (Using in Time Evolution mode).
 * @param isite1 [in] a site number \f$i \f$
 * @param isite2 [in] a site number \f$j \f$
 * @param isigma1 [in] a spin index at \f$i \f$ site.
 * @param isigma2 [in] a spin index at \f$j \f$ site.
 * @param dtmp_V [in] A value of general two-body diagonal interaction \f$ H_{i\sigma_1 j\sigma_2} \f$
 * @param X [in] Define list to get the operator information.
 * @param tmp_v0 [in,out] Result vector
 * @param tmp_v1 [in] Input produced vector
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 2.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalTEInterAll(
        long unsigned int isite1,
        long unsigned int isite2,
        long unsigned int isigma1,
        long unsigned int isigma2,
        double dtmp_V,
        struct BindStruct *X,
        double complex *tmp_v0,
        double complex *tmp_v1
){
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
  double complex dam_pr=0.0;


  /*
  Forse isite1 <= isite2
  */
  if (isite2 < isite1) {
    j = isite2;
    isite2 = isite1;
    isite1 = j;
    j = isigma2;
    isigma2 = isigma1;
    isigma1 = j;
  }
  /*
  When isite1 & site2 are in the inter process regino
  */
  if (isite1 > X->Def.Nsite) {

    switch (X->Def.iCalcModel) {

      case HubbardGC:
      case KondoGC:
      case Hubbard:
      case Kondo:
        is1_spin = X->Def.Tpow[2 * isite1 - 2 + isigma1];
        is2_spin = X->Def.Tpow[2 * isite2 - 2 + isigma2];
        num1 = 0;
        ibit1_spin = (unsigned long int)myrank&is1_spin;
        num1 += ibit1_spin / is1_spin;
        num2 = 0;
        ibit2_spin = (unsigned long int)myrank&is2_spin;
        num2 += ibit2_spin / is2_spin;
        break;/*case HubbardGC, KondoGC, Hubbard, Kondo:*/

      case SpinGC:
      case Spin:
        if (X->Def.iFlgGeneralSpin == FALSE) {
          is1_up = X->Def.Tpow[isite1 - 1];
          is2_up = X->Def.Tpow[isite2 - 1];
          num1 = child_SpinGC_CisAis((unsigned long int) myrank + 1, X, is1_up, isigma1);
          num2 = child_SpinGC_CisAis((unsigned long int) myrank + 1, X, is2_up, isigma2);
        }/*if (X->Def.iFlgGeneralSpin == FALSE)*/
        else {//start:generalspin
          num1 = BitCheckGeneral((unsigned long int) myrank, isite1, isigma1,
                                 X->Def.SiteToBit, X->Def.Tpow);
          num2 = BitCheckGeneral((unsigned long int) myrank, isite2, isigma2,
                                 X->Def.SiteToBit, X->Def.Tpow);
        }
        break;/*case SpinGC, Spin:*/

      default:
        fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
        return -1;
    }/*if (isite1 > X->Def.Nsite)*/

    if (num1 * num2 != 0) {
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1) \
firstprivate(i_max, dtmp_V) private(j)
      for (j = 1; j <= i_max; j++) {
        tmp_v0[j] += dtmp_V * tmp_v1[j];
        dam_pr += dtmp_V * conj(tmp_v1[j]) * tmp_v1[j];
      }
    }
    dam_pr=SumMPI_dc(dam_pr);
    X->Large.prdct += dam_pr;
    return 0;

  }/*if (isite1 > X->Def.Nsite)*/
  else if (isite2 > X->Def.Nsite) {

    switch (X->Def.iCalcModel) {

      case HubbardGC:

        is1_spin = X->Def.Tpow[2 * isite1 - 2 + isigma1];
        is2_spin = X->Def.Tpow[2 * isite2 - 2 + isigma2];

        num2 = 0;
        ibit2_spin = (unsigned long int)myrank&is2_spin;
        num2 += ibit2_spin / is2_spin;
        if(num2 !=0) {
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1)\
                     firstprivate(i_max, dtmp_V, is1_spin) private(num1, ibit1_spin, j)
          for (j = 1; j <= i_max; j++) {
            num1 = 0;
            ibit1_spin = (j - 1) & is1_spin;
            num1 += ibit1_spin / is1_spin;
            tmp_v0[j] += dtmp_V * num1 * tmp_v1[j];
            dam_pr += dtmp_V * num1 * conj(tmp_v1[j]) * tmp_v1[j];
          }
        }
        break;/*case HubbardGC:*/

      case KondoGC:
      case Hubbard:
      case Kondo:

        is1_spin = X->Def.Tpow[2 * isite1 - 2 + isigma1];
        is2_spin = X->Def.Tpow[2 * isite2 - 2 + isigma2];

        num2 = 0;
        ibit2_spin = (unsigned long int)myrank&is2_spin;
        num2 += ibit2_spin / is2_spin;
        if(num2 !=0) {
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1, list_1)\
                     firstprivate(i_max, dtmp_V, is1_spin) private(num1, ibit1_spin, j)
          for (j = 1; j <= i_max; j++) {
            num1 = 0;
            ibit1_spin = list_1[j] & is1_spin;
            num1 += ibit1_spin / is1_spin;
            tmp_v0[j] += dtmp_V *num1*tmp_v1[j];
            dam_pr += dtmp_V * num1*conj(tmp_v1[j]) * tmp_v1[j];
          }
        }
        break;/*case KondoGC, Hubbard, Kondo:*/

      case SpinGC:

        if (X->Def.iFlgGeneralSpin == FALSE) {
          is1_up = X->Def.Tpow[isite1 - 1];
          is2_up = X->Def.Tpow[isite2 - 1];
          num2 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is2_up, isigma2);

          if(num2 !=0) {
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1)\
                     firstprivate(i_max, dtmp_V, is1_up, isigma1, X) private(num1, j)
            for (j = 1; j <= i_max; j++) {
              num1 = child_SpinGC_CisAis(j, X, is1_up, isigma1);
              tmp_v0[j] += dtmp_V * num1 * tmp_v1[j];
              dam_pr += dtmp_V * num1*conj(tmp_v1[j]) * tmp_v1[j];
            }
          }
        }/* if (X->Def.iFlgGeneralSpin == FALSE)*/
        else {//start:generalspin
          num2 = BitCheckGeneral((unsigned long int)myrank, isite2, isigma2,
                                 X->Def.SiteToBit, X->Def.Tpow);
          if (num2 != 0) {
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1) \
firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
            for (j = 1; j <= i_max; j++) {
              num1 = BitCheckGeneral(j - 1, isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
              tmp_v0[j] += dtmp_V * num1 * tmp_v1[j];
              dam_pr += dtmp_V * num1*conj(tmp_v1[j]) * tmp_v1[j];
            }
          }
        }/* if (X->Def.iFlgGeneralSpin == TRUE)*/

        break;/*case SpinGC:*/

      case Spin:

        if (X->Def.iFlgGeneralSpin == FALSE) {
          is1_up = X->Def.Tpow[isite1 - 1];
          is2_up = X->Def.Tpow[isite2 - 1];
          num2 = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is2_up, isigma2);

          if(num2 !=0) {
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1) \
firstprivate(i_max, dtmp_V, is1_up, isigma1, X, num2) private(j, num1)
            for (j = 1; j <= i_max; j++) {
              num1 = child_Spin_CisAis(j, X, is1_up, isigma1);
              tmp_v0[j] += dtmp_V * num1 * tmp_v1[j];
              dam_pr += dtmp_V * num1*conj(tmp_v1[j]) * tmp_v1[j];
            }
          }
        }/* if (X->Def.iFlgGeneralSpin == FALSE)*/
        else /* if (X->Def.iFlgGeneralSpin == TRUE)*/{
          num2 = BitCheckGeneral((unsigned long int)myrank, isite2, isigma2, \
          X->Def.SiteToBit, X->Def.Tpow);
          if (num2 != 0) {
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1, list_1)\
firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
            for (j = 1; j <= i_max; j++) {
              num1 = BitCheckGeneral(list_1[j], isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
              tmp_v0[j] += dtmp_V * num1 * tmp_v1[j];
              dam_pr += dtmp_V * num1*conj(tmp_v1[j]) * tmp_v1[j];
            }
          }
        } /* if (X->Def.iFlgGeneralSpin == TRUE)*/

        break;/*case Spin:*/

      default:
        fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
        return -1;

    }/*switch (X->Def.iCalcModel)*/
    dam_pr=SumMPI_dc(dam_pr);
    X->Large.prdct += dam_pr;
    return 0;
  }/*else if (isite2 > X->Def.Nsite)*/

  switch (X->Def.iCalcModel){
    case HubbardGC: //list_1[j] -> j-1
      is1_spin   = X->Def.Tpow[2*isite1-2+isigma1];
      is2_spin   = X->Def.Tpow[2*isite2-2+isigma2];
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, is1_spin, is2_spin) private(num1, ibit1_spin, num2, ibit2_spin)
      for(j = 1;j <= i_max;j++){
        num1=0;
        num2=0;
        ibit1_spin=(j-1)&is1_spin;
        num1+=ibit1_spin/is1_spin;
        ibit2_spin=(j-1)&is2_spin;
        num2+=ibit2_spin/is2_spin;
        tmp_v0[j] += dtmp_V * num1*num2*tmp_v1[j];
        dam_pr += dtmp_V * num1*num2*conj(tmp_v1[j]) * tmp_v1[j];
      }
      break;
    case KondoGC:
    case Hubbard:
    case Kondo:
      is1_spin  = X->Def.Tpow[2*isite1-2+isigma1];
      is2_spin = X->Def.Tpow[2*isite2-2+isigma2];

#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1, list_1) firstprivate(i_max, dtmp_V, is1_spin, is2_spin) private(num1, ibit1_spin, num2, ibit2_spin)
      for(j = 1;j <= i_max;j++){
        num1=0;
        num2=0;
        ibit1_spin=list_1[j]&is1_spin;
        num1+=ibit1_spin/is1_spin;

        ibit2_spin=list_1[j]&is2_spin;
        num2+=ibit2_spin/is2_spin;
        tmp_v0[j] += dtmp_V * num1*num2*tmp_v1[j];
        dam_pr += dtmp_V * num1*num2*conj(tmp_v1[j]) * tmp_v1[j];
      }
      break;

    case Spin:
      if(X->Def.iFlgGeneralSpin==FALSE){
        is1_up   = X->Def.Tpow[isite1-1];
        is2_up   = X->Def.Tpow[isite2-1];
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1)  firstprivate(i_max, dtmp_V, is1_up, is2_up, isigma1, isigma2, X) private(j, num1, num2)
        for(j = 1;j <= i_max; j++){
          num1=child_Spin_CisAis(j, X, is1_up, isigma1);
          num2=child_Spin_CisAis(j, X, is2_up, isigma2);
          tmp_v0[j] += dtmp_V * num1*num2*tmp_v1[j];
          dam_pr += dtmp_V * num1*num2*conj(tmp_v1[j]) * tmp_v1[j];
        }
      }
      else{
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1, list_1) firstprivate(i_max, dtmp_V, isite1, isite2, isigma1, isigma2, X) private(j, num1)
        for(j = 1;j <= i_max; j++){
          num1=BitCheckGeneral (list_1[j], isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
          if(num1 != 0){
            num1=BitCheckGeneral (list_1[j], isite2, isigma2, X->Def.SiteToBit, X->Def.Tpow);
            tmp_v0[j] += dtmp_V *num1*tmp_v1[j];
            dam_pr += dtmp_V * num1*conj(tmp_v1[j]) * tmp_v1[j];
          }
        }

      }
      break;

    case SpinGC:
      if(X->Def.iFlgGeneralSpin==FALSE){
        is1_up   = X->Def.Tpow[isite1-1];
        is2_up   = X->Def.Tpow[isite2-1];
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, is1_up, is2_up, isigma1, isigma2, X) private(j, num1, num2)
        for(j = 1;j <= i_max; j++){
          num1=child_SpinGC_CisAis(j, X, is1_up, isigma1);
          num2=child_SpinGC_CisAis(j, X, is2_up, isigma2);
          tmp_v0[j] += dtmp_V * num1*num2*tmp_v1[j];
          dam_pr += dtmp_V * num1*num2*conj(tmp_v1[j]) * tmp_v1[j];
        }
      }
      else{//start:generalspin
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, isite1, isite2, isigma1, isigma2, X) private(j, num1)
        for(j = 1;j <= i_max; j++){
          num1=BitCheckGeneral (j-1, isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
          if(num1 != 0){
            num1=BitCheckGeneral (j-1, isite2, isigma2, X->Def.SiteToBit, X->Def.Tpow);
            tmp_v0[j] += dtmp_V *num1*tmp_v1[j];
            dam_pr += dtmp_V * num1*conj(tmp_v1[j]) * tmp_v1[j];
          }
        }
      }
      break;

    default:
      fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
      return -1;
  }
  dam_pr=SumMPI_dc(dam_pr);
  X->Large.prdct += dam_pr;
  return 0;
}

/**
 *
 *
 * @brief Update the vector by the chemical potential \f$ \mu_{i \sigma_1} n_ {i \sigma_1} \f$ \n
 * generated by the commutation relation in terms of the general two-body interaction, \n
 * \f$ c_ {i \sigma_1} a_{j\sigma_2}c_ {j \sigma_2}a_ {i \sigma_1} = c_ {i \sigma_1}a_ {i \sigma_1}-c_ {i \sigma_1} a_ {i \sigma_1} c_ {j \sigma_2}a_{j\sigma_2}\f$ .
 * (Using in Time Evolution mode).
 * @param isite1 [in] a site number
 * @param spin [in] a spin number
 * @param dtmp_V [in] A value of coulombintra interaction \f$ \mu_{i \sigma_1} \f$
 * @param X [in] Define list to get dimension number
 * @param tmp_v0 [in,out] Result vector
 * @param tmp_v1 [in] Input produced vector
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 2.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int SetDiagonalTEChemi(
        long unsigned int isite1,
        long unsigned int spin,
        double dtmp_V,
        struct BindStruct *X,
        double complex *tmp_v0,
        double complex *tmp_v1
){
  long unsigned int is1_up;
  long unsigned int num1;
  long unsigned int isigma1 =spin;
  long unsigned int is1,ibit1;

  long unsigned int j;
  long unsigned int i_max=X->Check.idim_max;
  double complex dam_pr=0;

  /*
    When isite1 is in the inter process region
  */
  if (isite1 > X->Def.Nsite){

    switch (X->Def.iCalcModel) {

      case HubbardGC:
      case KondoGC:
      case Hubbard:
      case Kondo:

        if (spin == 0) {
          is1 = X->Def.Tpow[2 * isite1 - 2];
        }
        else {
          is1 = X->Def.Tpow[2 * isite1 - 1];
        }
        ibit1 = (unsigned long int)myrank & is1;
        num1 = ibit1 / is1;
        break;/*case HubbardGC, case KondoGC, Hubbard, Kondo:*/

      case SpinGC:
      case Spin:

        if (X->Def.iFlgGeneralSpin == FALSE) {
          is1_up = X->Def.Tpow[isite1 - 1];
          num1 = (((unsigned long int)myrank& is1_up) / is1_up) ^ (1 - spin);
        } /*if (X->Def.iFlgGeneralSpin == FALSE)*/
        else /*if (X->Def.iFlgGeneralSpin == TRUE)*/ {
          num1 = BitCheckGeneral((unsigned long int)myrank,
                                 isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
        }/*if (X->Def.iFlgGeneralSpin == TRUE)*/
        break;/*case SpinGC, Spin:*/

      default:
        fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
        return -1;

    } /*switch (X->Def.iCalcModel)*/
    if (num1 != 0) {
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1)  \
firstprivate(i_max, dtmp_V) private(j)
      for (j = 1; j <= i_max; j++){
        tmp_v0[j] += dtmp_V * tmp_v1[j];
        dam_pr += dtmp_V * conj(tmp_v1[j]) * tmp_v1[j];
      }
    }/*if (num1 != 0)*/
    dam_pr=SumMPI_dc(dam_pr);
    X->Large.prdct += dam_pr;
    return 0;

  }/*if (isite1 >= X->Def.Nsite*/

  switch (X->Def.iCalcModel){
    case HubbardGC:
      if(spin==0){
        is1   = X->Def.Tpow[2*isite1-2];
      }else{
        is1 = X->Def.Tpow[2*isite1-1];
      }

      #pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, is1) private(num1, ibit1)
      for(j = 1;j <= i_max;j++){
        ibit1 = (j-1)&is1;
        num1  = ibit1/is1;
        tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
        dam_pr += dtmp_V * num1*conj(tmp_v1[j]) * tmp_v1[j];
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

#pragma omp parallel for default(none) reduction(+:dam_pr) shared(list_1, tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, is1) private(num1, ibit1)
      for(j = 1;j <= i_max;j++){

        ibit1 = list_1[j]&is1;
        num1  = ibit1/is1;
        tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
        dam_pr += dtmp_V * num1*conj(tmp_v1[j]) * tmp_v1[j];
      }
      break;

    case SpinGC:
      if(X->Def.iFlgGeneralSpin==FALSE){
        is1_up   = X->Def.Tpow[isite1-1];
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(list_1, tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, is1_up, spin) private(num1)
        for(j = 1;j <= i_max;j++){
          num1=(((j-1)& is1_up)/is1_up)^(1-spin);
          tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
          dam_pr += dtmp_V * num1*conj(tmp_v1[j]) * tmp_v1[j];
        }
      }
      else{
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
        for(j = 1;j <= i_max; j++){
          num1=BitCheckGeneral (j-1, isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
          if(num1 != 0){
            tmp_v0[j] += dtmp_V * tmp_v1[j];
            dam_pr += dtmp_V * conj(tmp_v1[j]) * tmp_v1[j];
          }
        }
      }
      break;

    case Spin:
      if(X->Def.iFlgGeneralSpin==FALSE){
        is1_up   = X->Def.Tpow[isite1-1];
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(list_1, tmp_v0, tmp_v1) firstprivate(i_max, dtmp_V, is1_up, spin) private(num1)
        for(j = 1;j <= i_max;j++){
          num1=((list_1[j]& is1_up)/is1_up)^(1-spin);
          tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
          dam_pr += dtmp_V * num1*conj(tmp_v1[j]) * tmp_v1[j];
        }
      }
      else{
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1, list_1) firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
        for(j = 1;j <= i_max; j++){
          num1=BitCheckGeneral (list_1[j], isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
          if(num1 != 0){
            tmp_v0[j] += dtmp_V * tmp_v1[j];
            dam_pr += dtmp_V * conj(tmp_v1[j]) * tmp_v1[j];

          }
        }
      }

      break;
    default:
      fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
      return -1;
  }
  dam_pr=SumMPI_dc(dam_pr);
  X->Large.prdct += dam_pr;
  return 0;
}

/**
 *
 * @brief Update the vector by the general one-body diagonal interaction, \f$ \mu_{i\sigma_1} n_ {i\sigma_1}\f$.\n
 * (Using in Time Evolution mode).
 * @param isite1 [in] a site number \f$i \f$
 * @param dtmp_V [in] A value of general one-body diagonal interaction \f$ \mu_{i\sigma_1} \f$
 * @param spin [in] a spin index at \f$i \f$ site.
 * @param X [in] Define list to get the operator information.
 * @param tmp_v0 [in,out] Result vector
 * @param tmp_v1 [in] Input produced vector
 * @retval -1 fail to calculate the diagonal component.
 * @retval  0 succeed to calculate the diagonal component.
 *
 * @version 2.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */

int SetDiagonalTETransfer
        (
                long unsigned int isite1,
                double dtmp_V,
                long unsigned int spin,
                struct BindStruct *X,
                double complex *tmp_v0,
                double complex *tmp_v1
        ){
  long unsigned int is1_up;
  long unsigned int ibit1_up;
  long unsigned int num1;
  long unsigned int isigma1 =spin;
  long unsigned int is1,ibit1;
  double complex dam_pr = 0.0;

  long unsigned int j;
  long unsigned int i_max=X->Check.idim_max;

  /*
    When isite1 is in the inter process region
  */
  if (isite1 > X->Def.Nsite){

    switch (X->Def.iCalcModel) {

      case HubbardGC:
      case KondoGC:
      case Hubbard:
      case Kondo:
        if (spin == 0) {
          is1 = X->Def.Tpow[2 * isite1 - 2];
        }
        else {
          is1 = X->Def.Tpow[2 * isite1 - 1];
        }
        ibit1 = (unsigned long int)myrank & is1;
        num1 = ibit1 / is1;
        break;/*case HubbardGC, case KondoGC, Hubbard, Kondo:*/

      case SpinGC:
      case Spin:
        if (X->Def.iFlgGeneralSpin == FALSE) {
          is1_up = X->Def.Tpow[isite1 - 1];
          num1 = (((unsigned long int)myrank& is1_up) / is1_up) ^ (1 - spin);
        } /*if (X->Def.iFlgGeneralSpin == FALSE)*/
        else /*if (X->Def.iFlgGeneralSpin == TRUE)*/ {
          num1 = BitCheckGeneral((unsigned long int)myrank,
                                 isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
        }/*if (X->Def.iFlgGeneralSpin == TRUE)*/
        break;/*case SpinGC, Spin:*/

      default:
        fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
        return -1;

    } /*switch (X->Def.iCalcModel)*/

    if(num1 !=0) {
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1)\
                     firstprivate(i_max, dtmp_V) private(j)
      for (j = 1; j <= i_max; j++) {
        tmp_v0[j] += dtmp_V * tmp_v1[j];
        dam_pr += dtmp_V * conj(tmp_v1[j]) * tmp_v1[j];
      }
    }
  }/*if (isite1 >= X->Def.Nsite*/
  else {//(isite1 < X->Def.Nsite)
    switch (X->Def.iCalcModel) {
      case HubbardGC:
        if (spin == 0) {
          is1 = X->Def.Tpow[2 * isite1 - 2];
        } else {
          is1 = X->Def.Tpow[2 * isite1 - 1];
        }
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(list_1, tmp_v0, tmp_v1) \
        firstprivate(i_max, dtmp_V, is1) private(num1, ibit1)
        for (j = 1; j <= i_max; j++) {
          ibit1 = (j - 1) & is1;
          num1 = ibit1 / is1;
          tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
          dam_pr += dtmp_V * num1*conj(tmp_v1[j]) * tmp_v1[j];
        }
        break;

      case KondoGC:
      case Hubbard:
      case Kondo:
        if (spin == 0) {
          is1 = X->Def.Tpow[2 * isite1 - 2];
        } else {
          is1 = X->Def.Tpow[2 * isite1 - 1];
        }
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(list_1, tmp_v0, tmp_v1) \
        firstprivate(i_max, dtmp_V, is1) private(num1, ibit1)
        for (j = 1; j <= i_max; j++) {
          ibit1 = list_1[j] & is1;
          num1 = ibit1 / is1;
          tmp_v0[j] += dtmp_V * num1*tmp_v1[j];
          dam_pr += dtmp_V * num1*conj(tmp_v1[j]) * tmp_v1[j];
        }
        break;

      case SpinGC:
        if (X->Def.iFlgGeneralSpin == FALSE) {
          is1_up = X->Def.Tpow[isite1 - 1];
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(list_1, tmp_v0, tmp_v1) \
          firstprivate(i_max, dtmp_V, is1_up, spin) private(num1, ibit1_up)
          for (j = 1; j <= i_max; j++) {
            ibit1_up = (((j - 1) & is1_up) / is1_up) ^ (1 - spin);
            tmp_v0[j] += dtmp_V * ibit1_up*tmp_v1[j];
            dam_pr += dtmp_V * ibit1_up*conj(tmp_v1[j]) * tmp_v1[j];
          }
        } else {
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1) \
          firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
          for (j = 1; j <= i_max; j++) {
            num1 = BitCheckGeneral(j - 1, isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
            if (num1 != 0) {
              tmp_v0[j] += dtmp_V *tmp_v1[j];
              dam_pr += dtmp_V *conj(tmp_v1[j]) * tmp_v1[j];
            }
          }
        }
        break;

      case Spin:
        if (X->Def.iFlgGeneralSpin == FALSE) {
          is1_up = X->Def.Tpow[isite1 - 1];
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(list_1, tmp_v0, tmp_v1)\
           firstprivate(i_max, dtmp_V, is1_up, spin) private(num1, ibit1_up)
          for (j = 1; j <= i_max; j++) {
            ibit1_up = ((list_1[j] & is1_up) / is1_up) ^ (1 - spin);
            tmp_v0[j] += dtmp_V * ibit1_up * tmp_v1[j];
            dam_pr += dtmp_V * ibit1_up * conj(tmp_v1[j]) * tmp_v1[j];
          }
        } else {
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(list_1, tmp_v0, tmp_v1)\
          firstprivate(i_max, dtmp_V, isite1, isigma1, X) private(j, num1)
          for (j = 1; j <= i_max; j++) {
            num1 = BitCheckGeneral(list_1[j], isite1, isigma1, X->Def.SiteToBit, X->Def.Tpow);
            tmp_v0[j] += dtmp_V * num1 * tmp_v1[j];
            dam_pr += dtmp_V * num1 * conj(tmp_v1[j]) * tmp_v1[j];
          }
        }
        break;

      default:
        fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
        return -1;
    }
  }
  dam_pr=SumMPI_dc(dam_pr);
  X->Large.prdct += dam_pr;
  return 0;
}
