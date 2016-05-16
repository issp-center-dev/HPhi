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
#include "Common.h"
#include "mltply.h"
#include "bitcalc.h"
#include "CalcSpectrum.h"
#include "CalcSpectrumByLanczos.h"
#include "SingleEx.h"
#include "wrapperMPI.h"
#include "mltplyMPI.h"

/**
 * @file   CalcSpectrum.c
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  File for givinvg functions of calculating spectrum 
 * 
 * 
 */

/** 
 * @brief A main function to calculate spectrum 
 * 
 * @param[in,out] X CalcStruct list for getting and pushing calculation information 
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 *
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Youhei Yamaji (The University of Tokyo)
 * 
 */
int CalcSpectrum(		 
		 struct EDMainCalStruct *X
				 )
{
  char sdt[D_FileNameMax];
  double diff_ene,var;
  double complex cdnorm;
  unsigned long int i;
  unsigned long int i_max=0;
  FILE *fp;
  double dnorm;
  
  fp = fopen(sdt, "rb");
  //set omega
  if(SetOmega(&(X->Bind.Def)) != TRUE){
    fprintf(stderr, "Error: Fail to set Omega.\n");
    fclose(fp);
    exitMPI(-1);
  }
  else{
    if(X->Bind.Def.iFlgSpecOmegaIm == FALSE){
      X->Bind.Def.dOmegaIm = (X->Bind.Def.dOmegaMax - X->Bind.Def.dOmegaMin)/(double) X->Bind.Def.iNOmega;
    }
  }

  
  //input eigen vector
  fprintf(stdoutMPI, "An Eigenvector is inputted in CalcSpectrum.\n");
  sprintf(sdt, cFileNameInputEigen, X->Bind.Def.CDataFileHead, X->Bind.Def.k_exct-1, myrank);

  if(fp==NULL){
    fprintf(stderr, "Error: A file of Inputvector does not exist.\n");
    fclose(fp);
    exitMPI(-1);
  }
  fread(&i_max, sizeof(long int), 1, fp);
  if(i_max != X->Bind.Check.idim_max){
    fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
    fclose(fp);
    exitMPI(-1);
  }
  fread(v1, sizeof(complex double), X->Bind.Check.idim_max+1, fp);
  fclose(fp);
  //mltply Operator
  fprintf(stdoutMPI, "Starting mltply operators in CalcSpectrum.\n");
  GetExcitedState( &(X->Bind), v0, v1);  

  //calculate norm
  fprintf(stdoutMPI, "Calculationg norm in CalcSpectrum.\n");
  dnorm = NormMPI_dc(i_max, v0);
  
  //normalize vector
  fprintf(stdoutMPI, "Normalizing the wave function in CalcSpectrum.\n");
#pragma omp parallel for default(none) private(i) shared(v1, v0) firstprivate(i_max, dnorm)
  for(i=1;i<=i_max;i++){
    v1[i] = v0[i]/dnorm;
  }
  fprintf(stdoutMPI, "The wave function normalized in CalcSpectrum.\n");

  int CalcSpecByLanczos=0;
  int iCalcSpecType=CalcSpecByLanczos;
  int iret=TRUE;
  switch (iCalcSpecType){
  case 0:
    iret = CalcSpectrumByLanczos(X, v1, dnorm);
    if(iret != TRUE){
      //Error Message will be added.
      return FALSE;
    }    
    break;    
    // case CalcSpecByShiftedKlyrov will be added
  default:
    break;
  }
  fprintf(stdoutMPI, "End of CalcSpectrumBy* in CalcSpectrum.\n");

  return TRUE;
}

int GetExcitedState
(
 struct BindStruct *X,
 double complex *tmp_v0,/**< [out] Result v0 = H v1*/
 double complex *tmp_v1 /**< [in] v0 = H v1*/
 )
{
  if(!GetSingleExcitedState(X,tmp_v0, tmp_v1)==TRUE){
    return FALSE;
  }

  if(!GetPairExcitedState(X,tmp_v0, tmp_v1)==TRUE){
    return FALSE;
  }
  
  return TRUE;
}

int GetSingleExcitedState
(
 struct BindStruct *X,
 double complex *tmp_v0, /**< [out] Result v0 = H v1*/
  double complex *tmp_v1 /**< [in] v0 = H v1*/
 ){

  long int idim_max;
  long unsigned int i,j;
  long unsigned int org_isite,ispin,itype;
  long unsigned int is1_spin;
  double complex tmpphi;
  double complex tmp_dam_pr;
  long unsigned int *tmp_off;

  idim_max = X->Check.idim_max;

  //tmp_v0

  switch(X->Def.iCalcModel){
  case HubbardGC:
    // SingleEx  
    fprintf(stdoutMPI, "SingleOperation in GetSingleExcitedState Re= %lf ; Im= %lf.\n",
    creal(X->Def.ParaSingleExcitationOperator[0]),cimag(X->Def.ParaSingleExcitationOperator[0]));
    // X->Def.NSingleExcitationOperator 
    // X->Def.SingleExcitationOperator[0][0]
    // X->Def.ParaSingleExcitationOperator[0]
    // clear all elements of tmp_v0 to zero  
    for(i=0;i<X->Def.NSingleExcitationOperator;i++){
      org_isite = X->Def.SingleExcitationOperator[i][0];
      ispin     = X->Def.SingleExcitationOperator[i][1];
      itype     = X->Def.SingleExcitationOperator[i][2];
      tmpphi    = X->Def.ParaSingleExcitationOperator[i];
      if(itype == 1){
        if( org_isite >= X->Def.Nsite){
          tmp_dam_pr = X_GC_Cis_MPI(org_isite,ispin,tmpphi,tmp_v0,tmp_v1,idim_max,v1buf,X->Def.Tpow);
        }
        else{
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1, X)	\
  firstprivate(idim_max, tmpphi, org_isite, ispin) private(j, is1_spin, tmp_dam_pr, tmp_off)
          for(j=1;j<=idim_max;j++){
            is1_spin = X->Def.Tpow[2*org_isite+ispin];
            tmp_dam_pr = GC_Cis(j,tmp_v0,tmp_v1,is1_spin,tmpphi,tmp_off); 
          }
        }
      }
      else if(itype == 0){
        if( org_isite >= X->Def.Nsite){
          tmp_dam_pr = X_GC_Ajt_MPI(org_isite,ispin,tmpphi,tmp_v0,tmp_v1,idim_max,v1buf,X->Def.Tpow);
        }
        else{
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1, X)	\
  firstprivate(idim_max, tmpphi, org_isite, ispin) private(j, is1_spin, tmp_dam_pr, tmp_off)
          for(j=1;j<=idim_max;j++){
            is1_spin = X->Def.Tpow[2*org_isite+ispin];
            tmp_dam_pr = GC_Ajt(j,tmp_v0,tmp_v1,is1_spin,tmpphi,tmp_off); 
          }
        }
      }
    }

    break;
  
  case KondoGC:
  case Hubbard:
  case Kondo:
  case Spin:
  case SpinGC:
    return FALSE;
    break;

  default:
    return FALSE;
  }

  return TRUE;
}


int GetPairExcitedState
(
 struct BindStruct *X,
 double complex *tmp_v0, /**< [out] Result v0 = H v1*/
 double complex *tmp_v1 /**< [in] v0 = H v1*/
 )
{

  FILE *fp;
  char sdt[D_FileNameMax];

  long unsigned int i,j;
  long unsigned int irght,ilft,ihfbit;
  long unsigned int isite1;
  long unsigned int org_isite1,org_isite2,org_sigma1,org_sigma2;
  long unsigned int tmp_off=0;
  double complex tmp_trans=0;
  long int i_max;
  int tmp_sgn, num1;
  int idx=0;
  int ihermite=0;
  long int ibit1, ibit;
  long unsigned int is1_up, is;
  //For TPQ
  int step=0;
  int rand_i=0;

  i_max = X->Check.idim_max;
  if(GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit)!=0){
    return -1;
  }
  X->Large.i_max    = i_max;
  X->Large.irght    = irght;
  X->Large.ilft     = ilft;
  X->Large.ihfbit   = ihfbit;
  X->Large.mode     = M_MLTPLY;

  switch(X->Def.iCalcModel){
  case HubbardGC:

    for(i=0;i<X->Def.NPairExcitationOperator;i++){
      org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
      org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
      org_sigma1 = X->Def.PairExcitationOperator[i][1];
      org_sigma2 = X->Def.PairExcitationOperator[i][3];
      tmp_trans = X->Def.ParaPairExcitationOperator[i];
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
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1)	\
  firstprivate(i_max, tmp_trans) private(j)
	    for (j = 1; j <= i_max; j++) tmp_v0[j] += tmp_trans*tmp_v1[j];
	  }
	}
	else{
	  X_GC_child_general_hopp_MPIdouble(org_isite1-1, org_sigma1, org_isite2-1, org_sigma2, -tmp_trans, X, tmp_v0, tmp_v1);
	}
      }
      else if (org_isite2  > X->Def.Nsite || org_isite1  > X->Def.Nsite){
	if(org_isite1<org_isite2){
	  X_GC_child_general_hopp_MPIsingle(org_isite1-1, org_sigma1, org_isite2-1, org_sigma2, -tmp_trans, X, tmp_v0, tmp_v1);
	}
	else{
	  X_GC_child_general_hopp_MPIsingle(org_isite2-1, org_sigma2, org_isite1-1, org_sigma1, -conj(tmp_trans), X, tmp_v0, tmp_v1);
	}
      }
      else{
	if(child_general_hopp_GetInfo( X,org_isite1,org_isite2,org_sigma1,org_sigma2)!=0){
	  return -1;
	}
	GC_child_general_hopp(tmp_v0, tmp_v1, X, tmp_trans);
      }
    }
    break;

  case KondoGC:
  case Hubbard:
  case Kondo:
    for(i=0;i<X->Def.NPairExcitationOperator;i++){
      org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
      org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
      org_sigma1 = X->Def.PairExcitationOperator[i][1];
      org_sigma2 = X->Def.PairExcitationOperator[i][3];
      tmp_trans = X->Def.ParaPairExcitationOperator[i];

      if(X->Def.iFlgSzConserved ==TRUE){
	if(org_sigma1 != org_sigma2){
	  continue;
	}
      }
      if (org_isite1  > X->Def.Nsite &&
	  org_isite2  > X->Def.Nsite) {
	if(org_isite1==org_isite2 && org_sigma1==org_sigma2){//diagonal

	  is   = X->Def.Tpow[2 * org_isite1 - 2+org_sigma1];
	  ibit = (unsigned long int)myrank & is;
	  if (ibit == is) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1)	\
  firstprivate(i_max, tmp_trans) private(j)
	    for (j = 1; j <= i_max; j++) tmp_v0[j] += tmp_trans*tmp_v1[j];
	  }

	}
	else{
	  X_child_general_hopp_MPIdouble(org_isite1-1, org_sigma1, org_isite2-1, org_sigma2, -tmp_trans, X, tmp_v0, tmp_v1);
	}
      }
      else if (org_isite2  > X->Def.Nsite || org_isite1  > X->Def.Nsite){
	if(org_isite1 < org_isite2){
	  X_child_general_hopp_MPIsingle(org_isite1-1, org_sigma1,org_isite2-1, org_sigma2, -tmp_trans, X, tmp_v0, tmp_v1);
	}
	else{
	  X_child_general_hopp_MPIsingle(org_isite2-1, org_sigma2, org_isite1-1, org_sigma1, -conj(tmp_trans), X, tmp_v0, tmp_v1);
	}
      }
      else{
	if(child_general_hopp_GetInfo( X,org_isite1,org_isite2,org_sigma1,org_sigma2)!=0){
	  return -1;
	}
	if(org_isite1==org_isite2 && org_sigma1==org_sigma2){

	  is   = X->Def.Tpow[2 * org_isite1 - 2 + org_sigma1];

#pragma omp parallel for default(none) shared(list_1, tmp_v0, tmp_v1) firstprivate(i_max, is, tmp_trans) private(num1, ibit)
	  for(j = 1;j <= i_max;j++){
	    ibit = list_1[j]&is;
	    num1  = ibit/is;
	    tmp_v0[j] += tmp_trans*num1*tmp_v1[j];
	  }
	}
	else{
	  child_general_hopp(tmp_v0, tmp_v1,X,tmp_trans);
	}
      }
    }
    break;

  case Spin: // for the Sz-conserved spin system

    if(X->Def.iFlgGeneralSpin==FALSE){
      for(i=0;i<X->Def.NPairExcitationOperator;i++) {
	org_isite1 = X->Def.PairExcitationOperator[i][0] + 1;
	org_isite2 = X->Def.PairExcitationOperator[i][2] + 1;
	org_sigma1 = X->Def.PairExcitationOperator[i][1];
	org_sigma2 = X->Def.PairExcitationOperator[i][3];
	tmp_trans = X->Def.ParaPairExcitationOperator[i];
	if (org_sigma1 == org_sigma2) {
	  if (org_isite1 == org_isite2) {
	    if (org_isite1 > X->Def.Nsite) {
	      is1_up = X->Def.Tpow[org_isite1 - 1];
	      ibit1 = X_SpinGC_CisAis((unsigned long int) myrank + 1, X, is1_up, org_sigma1);
	      if (ibit1 != 0) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1)	\
  firstprivate(i_max, tmp_trans) private(j)
		for (j = 1; j <= i_max; j++) tmp_v0[j] += tmp_trans * tmp_v1[j];
	      }
	    }// org_isite1 > X->Def.Nsite
	    else {
	      isite1 = X->Def.Tpow[org_isite1 - 1];
#pragma omp parallel for default(none) private(j) firstprivate(i_max, isite1, org_sigma1, X, tmp_trans) shared(tmp_v0, tmp_v1)
	      for (j = 1; j <= i_max; j++) {
		tmp_v0[j] += X_Spin_CisAis(j, X, isite1, org_sigma1) * tmp_v1[j] * tmp_trans;
	      }
	    }
	  } else {
	    // for the canonical case
	  }

	}
      }
    }//FlgGeneralSpin=FALSE
    else{
      for(i=0;i<X->Def.NPairExcitationOperator;i++){
	org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
	org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
	org_sigma1 = X->Def.PairExcitationOperator[i][1];
	org_sigma2 = X->Def.PairExcitationOperator[i][3];
	tmp_trans = X->Def.ParaPairExcitationOperator[i];
	if(org_isite1 == org_isite2){
	  if(org_isite1 >X->Def.Nsite){
	    if(org_sigma1==org_sigma2){
	      // longitudinal magnetic field
	      num1 = BitCheckGeneral((unsigned long int)myrank,
				     org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
	      if (num1 != 0) {
#pragma omp parallel for default(none) private(j) firstprivate(i_max, tmp_trans) shared(tmp_v0,tmp_v1)
		for(j=1;j<=i_max;j++){
		  tmp_v0[j]+= tmp_trans*tmp_v1[j];
		}
	      }
	    }
	  }
	  else {//org_isite1 <= X->Def.Nsite
	    if(org_sigma1==org_sigma2){
	      // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, X, tmp_trans) shared(tmp_v0,tmp_v1, list_1)
	      for(j=1;j<=i_max;j++){
		num1 = BitCheckGeneral(list_1[j], org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
		tmp_v0[j]+= tmp_trans*tmp_v1[j]*num1;
	      }
	    }
	  }
	}else{
	  // hopping is not allowed in localized spin system
	}//org_isite1 != org_isite2
      }
    }//general spin

    break;

  case SpinGC:

    if(X->Def.iFlgGeneralSpin==FALSE){
      for(i=0;i<X->Def.NPairExcitationOperator;i++){
	org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
	org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
	org_sigma1 = X->Def.PairExcitationOperator[i][1];
	org_sigma2 = X->Def.PairExcitationOperator[i][3];
	tmp_trans = X->Def.ParaPairExcitationOperator[i];

	if(org_isite1 == org_isite2){
	  if(org_isite1 > X->Def.Nsite){
	    if(org_sigma1==org_sigma2){  // longitudinal magnetic field
	      X_GC_child_CisAis_spin_MPIdouble(org_isite1-1, org_sigma1, tmp_trans, X, tmp_v0, tmp_v1);
	    }
	    else{  // transverse magnetic field
	      X_GC_child_CisAit_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, tmp_trans, X, tmp_v0, tmp_v1);
	    }
	  }else{
	    isite1 = X->Def.Tpow[org_isite1-1];

	    if(org_sigma1==org_sigma2){
	      // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, tmp_sgn) firstprivate(i_max, isite1, org_sigma1, X,tmp_trans) shared(tmp_v0, tmp_v1)
	      for(j=1;j<=i_max;j++){
		tmp_v0[j] += X_SpinGC_CisAis(j, X, isite1, org_sigma1)*tmp_v1[j]*tmp_trans;
	      }
	    }else{
	      // transverse magnetic field
#pragma omp parallel for default(none) private(j, tmp_sgn, tmp_off) firstprivate(i_max, isite1, org_sigma2, X, tmp_trans) shared(tmp_v0, tmp_v1)
	      for(j=1;j<=i_max;j++){
		tmp_sgn  =  X_SpinGC_CisAit(j,X, isite1,org_sigma2,&tmp_off);
		if(tmp_sgn !=0){
		  tmp_v0[j]+= tmp_sgn*tmp_v1[j]*tmp_trans;
		}
	      }
	    }
	  }
	}else{
	  // hopping is not allowed in localized spin system
	}
      }

    }//FlgGeneralSpin=FALSE
    else{
      for(i=0;i<X->Def.NPairExcitationOperator;i++){
	org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
	org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
	org_sigma1 = X->Def.PairExcitationOperator[i][1];
	org_sigma2 = X->Def.PairExcitationOperator[i][3];
	tmp_trans = X->Def.ParaPairExcitationOperator[i];
	if(org_isite1 == org_isite2){
	  if(org_isite1 > X->Def.Nsite){
	    if(org_sigma1==org_sigma2){
	      // longitudinal magnetic field
	      X_GC_child_CisAis_GeneralSpin_MPIdouble(org_isite1-1, org_sigma1, tmp_trans, X, tmp_v0, tmp_v1);
	    }else{
	      // transverse magnetic field
	      X_GC_child_CisAit_GeneralSpin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, tmp_trans, X, tmp_v0, tmp_v1);
	    }
	  }
	  else{//org_isite1 <= X->Def.Nsite
	    if(org_sigma1==org_sigma2){
	      // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, X, tmp_trans) shared(tmp_v0, tmp_v1)
	      for(j=1;j<=i_max;j++){
		num1 = BitCheckGeneral(j-1, org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
		tmp_v0[j] += tmp_trans* tmp_v1[j]*num1;
	      }
	    }else{
	      // transverse magnetic field
#pragma omp parallel for default(none) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, org_sigma2, X,tmp_off, tmp_trans) shared(tmp_v0, tmp_v1)
	      for(j=1;j<=i_max;j++){
		num1 = GetOffCompGeneralSpin(j-1, org_isite1, org_sigma2, org_sigma1, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
		if(num1 !=0){
		  tmp_v0[j]  +=  tmp_trans*tmp_v1[j]*num1;
		}
	      }
	    }
	  }
	}else{
	  // hopping is not allowed in localized spin system
	}
      }
    }
    break;

  default:
    return FALSE;
  }

  return TRUE;
}

int SetOmega
(
 struct DefineList *X
){
  FILE *fp;
  char sdt[D_FileNameMax],ctmp[256];
  double domegaMax;
  double domegaMin;
  int istp;
  double E1, E2, E3, E4, Emax;

  if(X->iFlgSpecOmegaMax == TRUE && X->iFlgSpecOmegaMin == TRUE){
    return TRUE;
  }
  else{
    sprintf(sdt, cFileNameLanczosStep, X->CDataFileHead);
    fp=fopenMPI(sdt, "r");
    if(fp == NULL){
      fprintf(stdoutMPI, "Error: xx_Lanczos_Step.dat does not exist.\n");
      return FALSE;
    }      
    while(fgetsMPI(ctmp, 256, fp) != NULL){
      sscanf(ctmp, "%d %lf %lf %lf %lf %lf\n",
	     &istp,
	     &E1,
	     &E2,
	     &E3,
	     &E4,
	     &Emax);
    }
    if(istp < 4){
      fprintf(stdoutMPI, "Error: Lanczos step must be greater than 4 for using spectrum calculation.\n");
      return FALSE;
    }  
    //Read Lanczos_Step
    if(X->iFlgSpecOmegaMax == FALSE){
      X->dOmegaMax= Emax*(double)X->NsiteMPI;
    }
    if(X->iFlgSpecOmegaMax == FALSE){
      X->dOmegaMax= E1;
    }
  }
  
  return TRUE;
}
