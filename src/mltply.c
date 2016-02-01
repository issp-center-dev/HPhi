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

//Define Mode for mltply
// complex version
#ifdef MPI
#include "mpi.h"
#endif
#include <bitcalc.h>
#include "mfmemory.h"
#include "xsetmem.h"
#include "mltply.h"
#include "mltplyMPI.h"
#include "wrapperMPI.h"

/**
 *
 *
 * @param X
 * @param tmp_v0
 * @param tmp_v1
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int mltply(struct BindStruct *X, double complex *tmp_v0,double complex *tmp_v1) {

  long unsigned int j;
  long unsigned int i;
  long unsigned int off = 0;
  long unsigned int tmp_off = 0;
  long unsigned int tmp_off2 = 0;
  long unsigned int is1_spin = 0;
  long unsigned int irght=0;
  long unsigned int ilft=0;
  long unsigned int ihfbit=0;
  long unsigned int isite1, isite2, sigma1, sigma2;
  long unsigned int isite3, isite4, sigma3, sigma4;
  long unsigned int ibitsite1, ibitsite2, ibitsite3, ibitsite4;

  double complex dam_pr;
  double complex tmp_trans;
  long int tmp_sgn;
  double num1 = 0;
  /*[s] For InterAll */
  double complex tmp_V;
  double complex dmv=0;
  /*[e] For InterAll */

  /* SpinGCBoost */
  int flagBoost;
  char *filename = "inputBoost";
  FILE *fp;
  double complex *tmp_v2, *tmp_v3;

  long unsigned int i_max;
  int ihermite=0;
  int idx=0;
  i_max = X->Check.idim_max;
  X->Large.prdct = 0.0;
  dam_pr = 0.0;

  /* SpinGCBoost */
  flagBoost=0;
  c_malloc1(tmp_v2, i_max+1);
  c_malloc1(tmp_v3, i_max+1);

  if(i_max!=0){
    if (X->Def.iFlgGeneralSpin == FALSE) {
      if (GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit) != 0) {
	return -1;
      }
    }
    else{
      if(X->Def.iCalcModel==Spin){
	if (GetSplitBitForGeneralSpin(X->Def.Nsite, &ihfbit, X->Def.SiteToBit) != 0) {
	  return -1;
	}
      }
    }
  }
  else{
    irght=0;
    ilft=0;
    ihfbit=0;
  }
  X->Large.i_max = i_max;
  X->Large.irght = irght;
  X->Large.ilft = ilft;
  X->Large.ihfbit = ihfbit;
  X->Large.mode = M_MLTPLY;

#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max) shared(tmp_v0, tmp_v1, list_Diagonal)
  for (j = 1; j <= i_max; j++) {
    tmp_v0[j] += (list_Diagonal[j]) * tmp_v1[j];
    dam_pr += (list_Diagonal[j]) * conj(tmp_v1[j]) * tmp_v1[j];
  }
  X->Large.prdct += dam_pr;
  
  switch (X->Def.iCalcModel) {
    case HubbardGC:
      //Transfer
      for (i = 0; i < X->Def.EDNTransfer; i += 2) {

        if (X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite &&
          X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite) {
          GC_child_general_hopp_MPIdouble(i, X, tmp_v0, tmp_v1);
        }
        else if (X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite){
          GC_child_general_hopp_MPIsingle(i, X, tmp_v0, tmp_v1);
        }
        else if (X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite) {
          GC_child_general_hopp_MPIsingle(i+1, X, tmp_v0, tmp_v1);
        }
        else {
          for (ihermite = 0; ihermite<2; ihermite++) {
            idx = i + ihermite;
            isite1 = X->Def.EDGeneralTransfer[idx][0] + 1;
            isite2 = X->Def.EDGeneralTransfer[idx][2] + 1;
            sigma1 = X->Def.EDGeneralTransfer[idx][1];
            sigma2 = X->Def.EDGeneralTransfer[idx][3];
            if (child_general_hopp_GetInfo(X, isite1, isite2, sigma1, sigma2) != 0) {
              return -1;
            }
            tmp_trans = -X->Def.EDParaGeneralTransfer[idx];
            dam_pr = GC_child_general_hopp(tmp_v0, tmp_v1, X, tmp_trans);
            X->Large.prdct += dam_pr;
          }
        }
      }

      for (i = 0; i < X->Def.NInterAll_OffDiagonal; i+=2) {
	  isite1 = X->Def.InterAll_OffDiagonal[i][0] + 1;
	  isite2 = X->Def.InterAll_OffDiagonal[i][2] + 1;
	  isite3 = X->Def.InterAll_OffDiagonal[i][4] + 1;
	  isite4 = X->Def.InterAll_OffDiagonal[i][6] + 1;
	  sigma1 = X->Def.InterAll_OffDiagonal[i][1];
	  sigma2 = X->Def.InterAll_OffDiagonal[i][3];
	  sigma3 = X->Def.InterAll_OffDiagonal[i][5];
	  sigma4 = X->Def.InterAll_OffDiagonal[i][7];
	  tmp_V = X->Def.ParaInterAll_OffDiagonal[i];

	if(CheckPE(isite1-1, X)==TRUE || CheckPE(isite2-1, X)==TRUE ||
	   CheckPE(isite3-1, X)==TRUE || CheckPE(isite4-1, X)==TRUE){
	  ibitsite1 = X->Def.OrgTpow[2*isite1-2+sigma1] ;
	  ibitsite2 = X->Def.OrgTpow[2*isite2-2+sigma2] ;
	  ibitsite3 = X->Def.OrgTpow[2*isite3-2+sigma3] ;
	  ibitsite4 = X->Def.OrgTpow[2*isite4-2+sigma4] ;
	  if(ibitsite1 == ibitsite2 && ibitsite3 == ibitsite4){
	    
	    dam_pr = X_GC_child_CisAisCjtAjt_Hubbard_MPI(isite1-1, sigma1, 
							 isite3-1, sigma3, 
							 tmp_V, X, tmp_v0, tmp_v1);
	  }
	  else if(ibitsite1 == ibitsite2 && ibitsite3 != ibitsite4){
	    
	    dam_pr = X_GC_child_CisAisCjtAku_Hubbard_MPI(isite1-1, sigma1, 
							 isite3-1, sigma3, isite4-1, sigma4,
							 tmp_V, X, tmp_v0, tmp_v1);
	    
	  }
	  else if(ibitsite1 != ibitsite2 && ibitsite3 == ibitsite4){
	    
	    dam_pr = X_GC_child_CisAjtCkuAku_Hubbard_MPI(isite1-1, sigma1, isite2-1, sigma2,
							 isite3-1, sigma3, 
							 tmp_V, X, tmp_v0, tmp_v1);
	    
	  }
	  else if(ibitsite1 != ibitsite2 && ibitsite3 != ibitsite4){
	    dam_pr = X_GC_child_CisAjtCkuAlv_Hubbard_MPI(isite1-1, sigma1, isite2-1, sigma2,
							 isite3-1, sigma3, isite4-1, sigma4,
							 tmp_V, X, tmp_v0, tmp_v1);
	  }
      }//InterPE
      else{
	dam_pr=0.0;
	for(ihermite=0; ihermite<2; ihermite++){
	  idx=i+ihermite;
	  isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
	  isite2 = X->Def.InterAll_OffDiagonal[idx][2] + 1;
	  isite3 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
	  isite4 = X->Def.InterAll_OffDiagonal[idx][6] + 1;
	  sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
	  sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
	  sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
	  sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
	  tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];
	  
	  child_general_int_GetInfo(
				    i,
				    X,
				    isite1,
				    isite2,
				    isite3,
				    isite4,
				    sigma1,
				    sigma2,
				    sigma3,
				    sigma4,
				    tmp_V
				    );
	  dam_pr += GC_child_general_int(tmp_v0, tmp_v1, X);
	}
      }
	X->Large.prdct += dam_pr;
      }
      
      //Pair hopping
      for (i = 0; i < X->Def.NPairHopping; i += 2) {
        if (X->Def.PairHopping[i][0] + 1 > X->Def.Nsite &&
          X->Def.PairHopping[i][1] + 1 > X->Def.Nsite) {
          fprintf(stderr, "In mltply Pairhop\n Sorry This Interaction has not be supported in MPI yet.\n");
          exitMPI(-1);
        }
        else if (X->Def.PairHopping[i][1] + 1 > X->Def.Nsite) {
          fprintf(stderr, "In mltply Pairhop\n Sorry This Interaction has not be supported in MPI yet.\n");
          exitMPI(-1);
        }
        else if (X->Def.PairHopping[i][0] + 1 > X->Def.Nsite) {
          fprintf(stderr, "In mltply Pairhop\n Sorry This Interaction has not be supported in MPI yet.\n");
          exitMPI(-1);
        }
        else {
          for (ihermite = 0; ihermite<2; ihermite++) {
            idx = i + ihermite;
            child_pairhopp_GetInfo(idx, X);
            dam_pr = GC_child_pairhopp(tmp_v0, tmp_v1, X);
            X->Large.prdct += dam_pr;
          }
        }
      }/*for (i = 0; i < X->Def.NPairHopping; i += 2)*/
      
      //Exchange
      for (i = 0; i < X->Def.NExchangeCoupling; i += 2) {
        if(X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite &&
          X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite){
          fprintf(stderr, "In mltply ExchangeCoupl\n Sorry This Interaction has not be supported in MPI yet.\n");
          exitMPI(-1);
        }
        else if(X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite){
          fprintf(stderr, "In mltply ExchangeCoupl\n Sorry This Interaction has not be supported in MPI yet.\n");
          exitMPI(-1);
        }
        else if (X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite) {
          fprintf(stderr, "In mltply ExchangeCoupl\n Sorry This Interaction has not be supported in MPI yet.\n");
          exitMPI(-1);
        }
        else {
          for (ihermite = 0; ihermite<2; ihermite++) {
            idx = i + ihermite;
            child_exchange_GetInfo(idx, X);
            dam_pr = GC_child_exchange(tmp_v0, tmp_v1, X);
            X->Large.prdct += dam_pr;
          }
        }
      }/*for (i = 0; i < X->Def.NExchangeCoupling; i+=2)*/
      break;
      
  case KondoGC:
  case Hubbard:
  case Kondo:

      //Transfer
      for (i = 0; i < X->Def.EDNTransfer; i+=2) {
        if (X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite &&
          X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite) {
          child_general_hopp_MPIdouble(i, X, tmp_v0, tmp_v1);
        }
        else if (X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite) {
          child_general_hopp_MPIsingle(i, X, tmp_v0, tmp_v1);
        }
        else if (X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite) {
          child_general_hopp_MPIsingle(i + 1, X, tmp_v0, tmp_v1);
        }
        else {
          for (ihermite = 0; ihermite<2; ihermite++) {
            idx = i + ihermite;
            isite1 = X->Def.EDGeneralTransfer[idx][0] + 1;
            isite2 = X->Def.EDGeneralTransfer[idx][2] + 1;
            sigma1 = X->Def.EDGeneralTransfer[idx][1];
            sigma2 = X->Def.EDGeneralTransfer[idx][3];
            if (child_general_hopp_GetInfo(X, isite1, isite2, sigma1, sigma2) != 0) {
              return -1;
            }
            tmp_trans = -X->Def.EDParaGeneralTransfer[idx];
            X->Large.tmp_trans = tmp_trans;
            dam_pr = child_general_hopp(tmp_v0, tmp_v1, X, tmp_trans);
            X->Large.prdct += dam_pr;
          }
        }
      }
      
          //InterAll
      for (i = 0; i < X->Def.NInterAll_OffDiagonal; i+=2) {
	
	isite1 = X->Def.InterAll_OffDiagonal[i][0] + 1;
	isite2 = X->Def.InterAll_OffDiagonal[i][2] + 1;
	isite3 = X->Def.InterAll_OffDiagonal[i][4] + 1;
	isite4 = X->Def.InterAll_OffDiagonal[i][6] + 1;
	sigma1 = X->Def.InterAll_OffDiagonal[i][1];
	sigma2 = X->Def.InterAll_OffDiagonal[i][3];
	sigma3 = X->Def.InterAll_OffDiagonal[i][5];
	sigma4 = X->Def.InterAll_OffDiagonal[i][7];
	tmp_V = X->Def.ParaInterAll_OffDiagonal[i];

	dam_pr =0.0;
	if(CheckPE(isite1-1, X)==TRUE || CheckPE(isite2-1, X)==TRUE ||
	   CheckPE(isite3-1, X)==TRUE || CheckPE(isite4-1, X)==TRUE){
	  ibitsite1 = X->Def.OrgTpow[2*isite1-2+sigma1] ;
	  ibitsite2 = X->Def.OrgTpow[2*isite2-2+sigma2] ;
	  ibitsite3 = X->Def.OrgTpow[2*isite3-2+sigma3] ;
	  ibitsite4 = X->Def.OrgTpow[2*isite4-2+sigma4] ;
	  if(ibitsite1 == ibitsite2 && ibitsite3 == ibitsite4){	    
	    dam_pr += X_child_CisAisCjtAjt_Hubbard_MPI(isite1-1, sigma1, 
						      isite3-1, sigma3, 
						      tmp_V, X, tmp_v0, tmp_v1);
	  }
	  else if(ibitsite1 == ibitsite2 && ibitsite3 != ibitsite4){
	    dam_pr += X_child_CisAisCjtAku_Hubbard_MPI(isite1-1, sigma1, 
						      isite3-1, sigma3, isite4-1, sigma4,
						      tmp_V, X, tmp_v0, tmp_v1);
	  }
	  else if(ibitsite1 != ibitsite2 && ibitsite3 == ibitsite4){	
	    dam_pr += X_child_CisAjtCkuAku_Hubbard_MPI(isite1-1, sigma1, isite2-1, sigma2,
						      isite3-1, sigma3, 
						      tmp_V, X, tmp_v0, tmp_v1);
	  }
	  else if(ibitsite1 != ibitsite2 && ibitsite3 != ibitsite4){
	    dam_pr += X_child_CisAjtCkuAlv_Hubbard_MPI(isite1-1, sigma1, isite2-1, sigma2,
						      isite3-1, sigma3, isite4-1, sigma4,
						      tmp_V, X, tmp_v0, tmp_v1);
	  }	 
	}
	else{
	  for(ihermite=0; ihermite<2; ihermite++){
	    idx=i+ihermite;
	    isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
	    isite2 = X->Def.InterAll_OffDiagonal[idx][2] + 1;
	    isite3 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
	    isite4 = X->Def.InterAll_OffDiagonal[idx][6] + 1;
	    sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
	    sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
	    sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
	    sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
	    tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];

	    child_general_int_GetInfo(
				      i,
				      X,
				      isite1,
				      isite2,
				      isite3,
				      isite4,
				      sigma1,
				      sigma2,
				      sigma3,
				      sigma4,
				      tmp_V
				      );

	    dam_pr += child_general_int(tmp_v0, tmp_v1, X);
	    
	  }
	}
	X->Large.prdct += dam_pr;
      }
	//Pair hopping
      for (i = 0; i < X->Def.NPairHopping; i += 2) {
        if (X->Def.PairHopping[i][0] + 1 > X->Def.Nsite &&
          X->Def.PairHopping[i][1] + 1 > X->Def.Nsite) {
          fprintf(stderr, "In mltply Pairhop\n Sorry This Interaction has not be supported in MPI yet.\n");
          exitMPI(-1);
        }
        else if(X->Def.PairHopping[i][1] + 1 > X->Def.Nsite){
          fprintf(stderr, "In mltply Pairhop\n Sorry This Interaction has not be supported in MPI yet.\n");
          exitMPI(-1);
        }
        else if(X->Def.PairHopping[i][0] + 1 > X->Def.Nsite){
          fprintf(stderr, "In mltply Pairhop\n Sorry This Interaction has not be supported in MPI yet.\n");
          exitMPI(-1);
        }
        else {
          for (ihermite = 0; ihermite<2; ihermite++) {
            idx = i + ihermite;
            child_pairhopp_GetInfo(idx, X);
            dam_pr = child_pairhopp(tmp_v0, tmp_v1, X);
            X->Large.prdct += dam_pr;
          }/*for (ihermite = 0; ihermite<2; ihermite++)*/
        }
      }/*for (i = 0; i < X->Def.NPairHopping; i += 2)*/

	//Exchange
      for (i = 0; i < X->Def.NExchangeCoupling; i += 2) {
        if (X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite &&
          X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite) {
          fprintf(stderr, "In mltply ExchangeCoupl\n Sorry This Interaction has not be supported in MPI yet.\n");
          exitMPI(-1);
        }
        else if (X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite) {
          fprintf(stderr, "In mltply ExchangeCoupl\n Sorry This Interaction has not be supported in MPI yet.\n");
          exitMPI(-1);
        }
        else if (X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite) {
          fprintf(stderr, "In mltply ExchangeCoupl\n Sorry This Interaction has not be supported in MPI yet.\n");
          exitMPI(-1);
        }
        else {
          for (ihermite = 0; ihermite<2; ihermite++) {
            idx = i + ihermite;
            child_exchange_GetInfo(idx, X);
            dam_pr = child_exchange(tmp_v0, tmp_v1, X);
            X->Large.prdct += dam_pr;
          }/*for (ihermite = 0; ihermite<2; ihermite++)*/
        }
      }/*for (i = 0; i < X->Def.NExchangeCoupling; i += 2)*/
  
      break;
      
    case Spin:
      if (X->Def.iFlgGeneralSpin == FALSE) {
	//Transfer absorbed in Diagonal term.
	//InterAll
	for (i = 0; i < X->Def.NInterAll_OffDiagonal; i+=2) {

          if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite &&
            X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
            child_general_int_spin_MPIdouble(i, X, tmp_v0, tmp_v1);
          }
          else if (X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
            child_general_int_spin_MPIsingle(i, X, tmp_v0, tmp_v1);
          }
          else if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite) {
            child_general_int_spin_MPIsingle(i + 1, X, tmp_v0, tmp_v1);
          }
          else {
            for (ihermite = 0; ihermite<2; ihermite++) {
              idx = i + ihermite;
              isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
              isite2 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
              sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
              sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
              sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
              sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
              tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];
              child_general_int_spin_GetInfo(X, isite1, isite2, sigma1, sigma2, sigma3, sigma4, tmp_V);
              dam_pr = child_general_int_spin(tmp_v0, tmp_v1, X);
              X->Large.prdct += dam_pr;
            }
          }
	}
      	
	//Exchange	
        for (i = 0; i < X->Def.NExchangeCoupling; i += 2) {
          if (X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite &&
            X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite) {
            fprintf(stderr, "In mltply ExchangeCoupl\n Sorry This Interaction has not be supported in MPI yet.\n");
            exitMPI(-1);
          }
          else if (X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite) {
            fprintf(stderr, "In mltply ExchangeCoupl\n Sorry This Interaction has not be supported in MPI yet.\n");
            exitMPI(-1);
          }
          else if (X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite) {
            fprintf(stderr, "In mltply ExchangeCoupl\n Sorry This Interaction has not be supported in MPI yet.\n");
            exitMPI(-1);
          }
          else {
            for (ihermite = 0; ihermite<2; ihermite++) {
              idx = i + ihermite;
              child_exchange_spin_GetInfo(idx, X);
              dam_pr = child_exchange_spin(tmp_v0, tmp_v1, X);
              X->Large.prdct += dam_pr;
            }
          }
	}/*for (i = 0; i < X->Def.NExchangeCoupling; i += 2)*/
      }
      else{
	//Transfer absorbed in Diagonal term.
	//InterAll
	ihfbit =X->Check.sdim;
        for (i = 0; i < X->Def.NInterAll_OffDiagonal; i += 2) {

          if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite &&
            X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
            child_general_int_GeneralSpin_MPIdouble(i, X, tmp_v0, tmp_v1);
          }
          else if (X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
            child_general_int_GeneralSpin_MPIsingle(i, X, tmp_v0, tmp_v1);
          }
          else if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite) {
            child_general_int_GeneralSpin_MPIsingle(i + 1, X, tmp_v0, tmp_v1);
          }
          else {
            for (ihermite = 0; ihermite < 2; ihermite++) {
              idx = i + ihermite;
              isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
              isite2 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
              sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
              sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
              sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
              sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
              tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];
              dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) \
private(j, tmp_sgn, dmv, off, tmp_off, tmp_off2) \
firstprivate(i_max, isite1, isite2, sigma1, sigma2, sigma3, sigma4, X, tmp_V, ihfbit) \
shared(tmp_v0, tmp_v1, list_1, list_2_1, list_2_2)
              for (j = 1; j <= i_max; j++) {
                tmp_sgn = GetOffCompGeneralSpin(list_1[j], isite2, sigma4, sigma3, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
                if (tmp_sgn == TRUE) {
                  tmp_sgn = GetOffCompGeneralSpin(tmp_off, isite1, sigma2, sigma1, &tmp_off2, X->Def.SiteToBit, X->Def.Tpow);
                  if (tmp_sgn == TRUE) {
                    ConvertToList1GeneralSpin(tmp_off2, ihfbit, &off);
                    dmv = tmp_v1[j] * tmp_V;
                    if (X->Large.mode == M_MLTPLY) { // for multply
                      tmp_v0[off] += dmv;
                    }
                    dam_pr += conj(tmp_v1[off]) * dmv;
                  }
                }
              }
              X->Large.prdct += dam_pr;
            }
          }
        }
      }	
      break;

/* SpinGCBoost */
    case SpinGC:

    if((fp = fopen(filename, "r")) == NULL){
      fprintf(stderr, "\n\n ###Boost### failed to open a file %s\n\n", filename);
      exit(EXIT_FAILURE);
    }
    if(myrank==0){
      fscanf(fp, "%d", &flagBoost);
#ifdef MPI
      MPI_Bcast(&flagBoost, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
      fclose(fp);
    }
    if(myrank==0){printf("\n\n###Boost### SpinGC Boost mode flagBoost %d \n\n", flagBoost);}

    if(flagBoost == 0){
     
      if (X->Def.iFlgGeneralSpin == FALSE) {	
        for (i = 0; i < X->Def.EDNTransfer; i+=2 ) {
	  if(X->Def.EDGeneralTransfer[i][0]+1 > X->Def.Nsite){
	    dam_pr=0;
	    if(X->Def.EDGeneralTransfer[i][1]==X->Def.EDGeneralTransfer[i][3]){
	      fprintf(stderr, "Transverse_OffDiagonal component is illegal.\n");
	    }
	    else{
	      dam_pr += X_GC_child_CisAit_spin_MPIdouble(X->Def.EDGeneralTransfer[i][0], X->Def.EDGeneralTransfer[i][1], X->Def.EDGeneralTransfer[i][3], -X->Def.EDParaGeneralTransfer[i], X, tmp_v0, tmp_v1);
	    }
	  }
	  else{
	    dam_pr=0;
	    for(ihermite=0; ihermite<2; ihermite++){
	      idx=i+ihermite;
	      isite1 = X->Def.EDGeneralTransfer[idx][0] + 1;
	      isite2 = X->Def.EDGeneralTransfer[idx][2] + 1;
	      sigma1 = X->Def.EDGeneralTransfer[idx][1];
	      sigma2 = X->Def.EDGeneralTransfer[idx][3];
	      tmp_trans = -X->Def.EDParaGeneralTransfer[idx];
	      if (child_general_hopp_GetInfo(X, isite1, isite2, sigma1, sigma2) != 0) {
		return -1;
	      }
	      
	      if(sigma1==sigma2){
		fprintf(stderr, "Transverse_OffDiagonal component is illegal.\n");
	      }
	      else{
		// longitudinal magnetic field (considerd in diagonalcalc.c)
		// transverse magnetic field
		is1_spin = X->Def.Tpow[isite1 - 1];
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn) firstprivate(i_max, is1_spin, sigma2, X,off, tmp_trans) shared(tmp_v0, tmp_v1)
		for (j = 1; j <= i_max; j++) {
		  tmp_sgn = X_SpinGC_CisAit(j, X, is1_spin, sigma2, &off);
		  if(tmp_sgn !=0){
		    tmp_v0[off+1] += tmp_v1[j]*tmp_trans;
		    dam_pr += tmp_trans * conj(tmp_v1[off + 1]) * tmp_v1[j];
		  }
		}
	      }//sigma1 != sigma2
	    }
	  }
	  X->Large.prdct += dam_pr;
	}
      

	//InterAll	
        for (i = 0; i < X->Def.NInterAll_OffDiagonal; i+=2) {
          if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite &&
            X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
            GC_child_general_int_spin_MPIdouble(i, X, tmp_v0, tmp_v1);
          }
          else if (X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
            GC_child_general_int_spin_MPIsingle(i, X, tmp_v0, tmp_v1);
          }
          else if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite) {
            GC_child_general_int_spin_MPIsingle(i + 1, X, tmp_v0, tmp_v1);
          }
          else {
            for (ihermite = 0; ihermite < 2; ihermite++) {
              idx = i + ihermite;
              isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
              isite2 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
              sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
              sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
              sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
              sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
              tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];
              child_general_int_spin_GetInfo(X, isite1, isite2, sigma1, sigma2, sigma3, sigma4, tmp_V);
              dam_pr = GC_child_general_int_spin(tmp_v0, tmp_v1, X);
              X->Large.prdct += dam_pr;
            }
          }
	}
	
        //Exchange
        for (i = 0; i < X->Def.NExchangeCoupling; i += 2) {
          if (X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite &&
            X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite) {
            fprintf(stderr, "In mltply ExchangeCoupl\n Sorry This Interaction has not be supported in MPI yet.\n");
            exitMPI(-1);
          }
          else if (X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite) {
            fprintf(stderr, "In mltply ExchangeCoupl\n Sorry This Interaction has not be supported in MPI yet.\n");
            exitMPI(-1);
          }
          else if (X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite) {
            fprintf(stderr, "In mltply ExchangeCoupl\n Sorry This Interaction has not be supported in MPI yet.\n");
            exitMPI(-1);
          }
          else {
            for (ihermite = 0; ihermite<2; ihermite++) {
              idx = i + ihermite;
              child_exchange_spin_GetInfo(idx, X);
              dam_pr = GC_child_exchange_spin(tmp_v0, tmp_v1, X);
              X->Large.prdct += dam_pr;
            }/*for (ihermite = 0; ihermite<2; ihermite++)*/
          }
	}/*for (i = 0; i < X->Def.NExchangeCoupling; i += 2)*/

        //PairLift
        for (i = 0; i < X->Def.NPairLiftCoupling; i += 2) {
          if (X->Def.PairLiftCoupling[i][0] + 1 > X->Def.Nsite &&
            X->Def.PairLiftCoupling[i][1] + 1 > X->Def.Nsite) {
            fprintf(stderr, "In mltply PairLift\n Sorry This Interaction has not be supported in MPI yet.\n");
            exitMPI(-1);
          }
          else if (X->Def.PairLiftCoupling[i][1] + 1 > X->Def.Nsite) {
            fprintf(stderr, "In mltply PairLift\n Sorry This Interaction has not be supported in MPI yet.\n");
            exitMPI(-1);
          }
          else if (X->Def.PairLiftCoupling[i][0] + 1 > X->Def.Nsite) {
            fprintf(stderr, "In mltply PairLift\n Sorry This Interaction has not be supported in MPI yet.\n");
            exitMPI(-1);
          }
          else {
            for (ihermite = 0; ihermite<2; ihermite++) {
              idx = i + ihermite;
              child_pairlift_spin_GetInfo(idx, X);
              dam_pr = GC_child_pairlift_spin(tmp_v0, tmp_v1, X);
              X->Large.prdct += dam_pr;
            }/*for (ihermite = 0; ihermite<2; ihermite++)*/
          }
	}/*for (i = 0; i < X->Def.NPairLiftCoupling; i += 2)*/
      }
      else {//For General spin
        for (i = 0; i < X->Def.EDNTransfer; i += 2) {
          if (X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite &&
            X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite) {
            fprintf(stderr, "In mltply GSpin+TransMag\n Sorry This Interaction has not be supported in MPI yet.\n");
            exitMPI(-1);
          }
          else if(X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite){
            fprintf(stderr, "In mltply GSpin+TransMag\n Sorry This Interaction has not be supported in MPI yet.\n");
            exitMPI(-1);
          }
          else if(X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite){
            fprintf(stderr, "In mltply GSpin+TransMag\n Sorry This Interaction has not be supported in MPI yet.\n");
            exitMPI(-1);
          }
          else {
            for (ihermite = 0; ihermite<2; ihermite++) {
              idx = i + ihermite;
              isite1 = X->Def.EDGeneralTransfer[idx][0] + 1;
              isite2 = X->Def.EDGeneralTransfer[idx][2] + 1;
              sigma1 = X->Def.EDGeneralTransfer[idx][1];
              sigma2 = X->Def.EDGeneralTransfer[idx][3];
              tmp_trans = -X->Def.EDParaGeneralTransfer[idx];
              if (isite1 == isite2) {
                if (isite1 > X->Def.Nsite) {
                  if (sigma1 != sigma2) {
                    dam_pr = X_GC_child_CisAit_GeneralSpin_MPIdouble(isite1 - 1, sigma1, sigma2, tmp_trans, X, tmp_v0, tmp_v1);
                  }
                }
                else {
                  if (sigma1 != sigma2) {//sigma1 != sigma2
                                         // transverse magnetic field
                    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, num1) firstprivate(i_max, isite1, sigma1, sigma2, X, off, tmp_trans) shared(tmp_v0, tmp_v1)
                    for (j = 1; j <= i_max; j++) {
                      num1 = GetOffCompGeneralSpin(j - 1, isite1, sigma2, sigma1, &off, X->Def.SiteToBit, X->Def.Tpow);
                      if (num1 != 0) { // for multply
                        tmp_v0[off + 1] += tmp_v1[j] * tmp_trans;
                        dam_pr += conj(tmp_v1[off + 1]) * tmp_v1[j] * tmp_trans;
                      }
                    }
                  }
                }
              }
              else {
                // hopping is not allowed in localized spin system
                return -1;
              }
              X->Large.prdct += dam_pr;
            }/*for (ihermite = 0; ihermite<2; ihermite++)*/
          }
	}
      
        //InterAll        
        for (i = 0; i< X->Def.NInterAll_OffDiagonal; i += 2) {

          if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite &&
            X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
            GC_child_general_int_GeneralSpin_MPIdouble(i, X, tmp_v0, tmp_v1);
          }
          else if (X->Def.InterAll_OffDiagonal[i][4] + 1 > X->Def.Nsite) {
            GC_child_general_int_GeneralSpin_MPIsingle(i, X, tmp_v0, tmp_v1);
          }
          else if (X->Def.InterAll_OffDiagonal[i][0] + 1 > X->Def.Nsite) {
            GC_child_general_int_GeneralSpin_MPIsingle(i + 1, X, tmp_v0, tmp_v1);
          }
          else {
            for (ihermite = 0; ihermite < 2; ihermite++) {
              idx = i + ihermite;
              isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
              isite2 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
              sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
              sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
              sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
              sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
              tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];

              dam_pr = 0.0;
              if (sigma1 == sigma2) {
                if (sigma3 == sigma4) {
                  fprintf(stderr, "InterAll_OffDiagonal component is illegal.\n");
                  return -1;
                }
                else {
                  //sigma3=sigma4 term is considerd as a diagonal term.
#pragma omp parallel for default(none) reduction(+:dam_pr) \
private(j, tmp_sgn, dmv, off) \
firstprivate(i_max, isite1, isite2, sigma1, sigma3, sigma4, X, tmp_V) \
shared(tmp_v0, tmp_v1)
                  for (j = 1; j <= i_max; j++) {
                    tmp_sgn = GetOffCompGeneralSpin(j - 1, isite2, sigma4, sigma3, &off, X->Def.SiteToBit, X->Def.Tpow);
                    if (tmp_sgn == TRUE) {
                      tmp_sgn = BitCheckGeneral(off, isite1, sigma1, X->Def.SiteToBit, X->Def.Tpow);
                      if (tmp_sgn == TRUE) {
                        dmv = tmp_v1[j] * tmp_V;
                        if (X->Large.mode == M_MLTPLY) { // for multply
                          tmp_v0[off + 1] += dmv;
                        }
                        dam_pr += conj(tmp_v1[off + 1]) * dmv;
                      }
                    }
                  }
                }
              }
              else if (sigma3 == sigma4) {
                //sigma1=sigma2 term is considerd as a diagonal term.
#pragma omp parallel for default(none) reduction(+:dam_pr) \
private(j, tmp_sgn, dmv, off, tmp_off) \
firstprivate(i_max, isite1, isite2, sigma1, sigma2, sigma3, sigma4, X, tmp_V) \
shared(tmp_v0, tmp_v1)
                for (j = 1; j <= i_max; j++) {
                  tmp_sgn = BitCheckGeneral(j - 1, isite2, sigma3, X->Def.SiteToBit, X->Def.Tpow);
                  if (tmp_sgn == TRUE) {
                    tmp_sgn = GetOffCompGeneralSpin(j - 1, isite1, sigma2, sigma1, &off, X->Def.SiteToBit, X->Def.Tpow);
                    if (tmp_sgn == TRUE) {
                      dmv = tmp_v1[j] * tmp_V;
                      if (X->Large.mode == M_MLTPLY) { // for multply
                        tmp_v0[off + 1] += dmv;
                      }
                      dam_pr += conj(tmp_v1[off + 1]) * dmv;
                    }
                  }
                }
              }
              else {
#pragma omp parallel for default(none) reduction(+:dam_pr) \
private(j, tmp_sgn, dmv, off, tmp_off) \
firstprivate(i_max, isite1, isite2, sigma1, sigma2, sigma3, sigma4, X, tmp_V) \
shared(tmp_v0, tmp_v1)
                for (j = 1; j <= i_max; j++) {
                  tmp_sgn = GetOffCompGeneralSpin(j - 1, isite2, sigma4, sigma3, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
                  if (tmp_sgn == TRUE) {
                    tmp_sgn = GetOffCompGeneralSpin(tmp_off, isite1, sigma2, sigma1, &off, X->Def.SiteToBit, X->Def.Tpow);
                    if (tmp_sgn == TRUE) {
                      dmv = tmp_v1[j] * tmp_V;
                      if (X->Large.mode == M_MLTPLY) { // for multply
                        tmp_v0[off + 1] += dmv;
                      }
                      dam_pr += conj(tmp_v1[off + 1]) * dmv;
                    }
                  }
                }
              }
              X->Large.prdct += dam_pr;
            }
          }
        }
      }  //end:generalspin
  }
  else{  
    if(myrank==0){printf("\n\n###Boost### SpinGC Boost mode start \n\n");}
     
    child_general_int_spin_MPIBoost(X, tmp_v0, tmp_v1, tmp_v2, tmp_v3);
    
    if(myrank==0){printf("\n\n###Boost### SpinGC Boost mode step \n\n");}
  }
  break;
/* SpinGCBoost */
      
  default:
    return -1;
  }
  
  X->Large.prdct = SumMPI_dc(X->Large.prdct);
  //  fprintf(stdoutMPI, "debug: prdct=%lf, %lf\n",creal( X->Large.prdct), cimag( X->Large.prdct ) );
  //FinalizeMPI();
  //exit(0);

  /* SpinGCBoost */
  c_free1(tmp_v2, i_max+1);  
  c_free1(tmp_v3, i_max+1);  

  return 0;
}


/******************************************************************************/
//[s] child functions
/******************************************************************************/
/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_pairhopp
          (
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X
          ) {
    long int j;
    long unsigned int i_max = X->Large.i_max;
    long unsigned int off = 0;
    double complex dam_pr = 0.0;

#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max, X,off) private(j) shared(tmp_v0, tmp_v1)
    for (j = 1; j <= i_max; j++) {
      dam_pr += child_pairhopp_element(j, tmp_v0, tmp_v1, X, &off);
    }

    return dam_pr;
  }


/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_exchange
          (
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X
          ) {
    long int j;
    long unsigned int i_max = X->Large.i_max;
    long unsigned int off = 0;
    double complex dam_pr = 0;

#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max, X,off) private(j) shared(tmp_v0, tmp_v1)
    for (j = 1; j <= i_max; j++) {
      dam_pr += child_exchange_element(j, tmp_v0, tmp_v1, X, &off);
    }
    return dam_pr;
  }

/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_exchange_spin
          (
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X
          ) {
    long unsigned int j;
    long unsigned int i_max = X->Large.i_max;
    long unsigned int off = 0;
    double complex dam_pr = 0;

#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max, X,off) private(j) shared(tmp_v0, tmp_v1)
    for (j = 1; j <= i_max; j++) {
      dam_pr += child_exchange_spin_element(j, tmp_v0, tmp_v1, X, &off);
    }
    return dam_pr;
  }

/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_pairlift_spin
          (
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X
          ) {
    long int j;
    long unsigned int i_max = X->Large.i_max;
    long unsigned int off = 0;
    double complex dam_pr = 0;

#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max, X,off) private(j) shared(tmp_v0, tmp_v1)
    for (j = 1; j <= i_max; j++) {
      dam_pr += child_pairlift_spin_element(j, tmp_v0, tmp_v1, X, &off);
    }
    return dam_pr;
  }

/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param trans
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_general_hopp
          (
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  double complex trans
          ) {

    long unsigned int j, isite1, isite2, Asum, Adiff;
    long unsigned int i_max = X->Large.i_max;

    isite1 = X->Large.is1_spin;
    isite2 = X->Large.is2_spin;
    Asum = X->Large.isA_spin;
    Adiff = X->Large.A_spin;

    double complex dam_pr = 0;
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max,X,Asum,Adiff,isite1,isite2,trans) private(j) shared(tmp_v0, tmp_v1)
    for (j = 1; j <= i_max; j++) {
      dam_pr += CisAjt(j, tmp_v0, tmp_v1, X, isite1, isite2, Asum, Adiff, trans) * trans;
    }
    return dam_pr;
  }

/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param trans
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_general_hopp
          (
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  double complex trans
          ) {

    long unsigned int j, isite1, isite2, Asum, Adiff;
    long unsigned int tmp_off;
    long unsigned int i_max = X->Large.i_max;

    isite1 = X->Large.is1_spin;
    isite2 = X->Large.is2_spin;
    Asum = X->Large.isA_spin;
    Adiff = X->Large.A_spin;

    double complex dam_pr = 0;
    if(isite1==isite2 ){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isite1, trans) shared(tmp_v0, tmp_v1)
      for(j=1;j<=i_max;j++){
	dam_pr += GC_CisAis(j, tmp_v0, tmp_v1, X, isite1, trans) * trans;
	 
      }
    }
    else{
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max,X,Asum,Adiff,isite1,isite2,trans) private(j,tmp_off) shared(tmp_v0, tmp_v1)
      for (j = 1; j <= i_max; j++) {
	dam_pr += GC_CisAjt(j, tmp_v0, tmp_v1, X, isite1, isite2, Asum, Adiff, trans, &tmp_off) * trans;
      }
    }
    return dam_pr;
  }

/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_general_int(double complex *tmp_v0, double complex *tmp_v1, struct BindStruct *X) {
    double complex dam_pr, tmp_V;
    long unsigned int j, i_max;
    long unsigned int isite1, isite2, isite3, isite4;
    long unsigned int Asum, Bsum, Adiff, Bdiff;
    long unsigned int tmp_off = 0;
    long unsigned int tmp_off_2 = 0;

    //note: this site is labeled for not only site but site with spin.
    i_max = X->Large.i_max;
    isite1 = X->Large.is1_spin;
    isite2 = X->Large.is2_spin;
    Asum = X->Large.isA_spin;
    Adiff = X->Large.A_spin;

    isite3 = X->Large.is3_spin;
    isite4 = X->Large.is4_spin;
    Bsum = X->Large.isB_spin;
    Bdiff = X->Large.B_spin;

    tmp_V = X->Large.tmp_V;
    dam_pr = 0.0;

    if (isite1 == isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_off) firstprivate(i_max,X,isite1,isite3,tmp_V) shared(tmp_v0, tmp_v1)
      for (j = 1; j <= i_max; j++) {
        dam_pr += child_CisAisCisAis_element(j, isite1, isite3, tmp_V, tmp_v0, tmp_v1, X, &tmp_off);
      }
    } else if (isite1 == isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_off) firstprivate(i_max,X,isite1,isite4,isite3, Bsum, Bdiff, tmp_V) shared(tmp_v0, tmp_v1)
      for (j = 1; j <= i_max; j++) {
        dam_pr += child_CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, tmp_V, tmp_v0, tmp_v1, X,
                                             &tmp_off);
      }
    } else if (isite1 != isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j,tmp_off) firstprivate(i_max,X,isite1,isite2,isite3,Asum,Adiff,tmp_V) shared(tmp_v0, tmp_v1)
      for (j = 1; j <= i_max; j++) {
        dam_pr += child_CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, tmp_V, tmp_v0, tmp_v1, X,
                                             &tmp_off);
      }
    } else if (isite1 != isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_off_2) firstprivate(i_max,X,isite1,isite2,isite3,isite4,Asum,Bsum,Adiff,Bdiff, tmp_V) shared(tmp_v0, tmp_v1)
      for (j = 1; j <= i_max; j++) {
        dam_pr += child_CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, tmp_V, tmp_v0,
                                             tmp_v1, X, &tmp_off_2);

      }
    }
    return dam_pr;
  }

/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_general_int(double complex *tmp_v0, double complex *tmp_v1, struct BindStruct *X) {
    double complex dam_pr, tmp_V;
    long unsigned int j, i_max;
    long unsigned int isite1, isite2, isite3, isite4;
    long unsigned int Asum, Bsum, Adiff, Bdiff;
    long unsigned int tmp_off = 0;
    long unsigned int tmp_off_2 = 0;

    i_max = X->Large.i_max;
    isite1 = X->Large.is1_spin;
    isite2 = X->Large.is2_spin;
    Asum = X->Large.isA_spin;
    Adiff = X->Large.A_spin;

    isite3 = X->Large.is3_spin;
    isite4 = X->Large.is4_spin;
    Bsum = X->Large.isB_spin;
    Bdiff = X->Large.B_spin;

    tmp_V = X->Large.tmp_V;
    dam_pr = 0.0;

    if (isite1 == isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) shared(tmp_v0, tmp_v1)
      for (j = 1; j <= i_max; j++) {
        dam_pr += GC_child_CisAisCisAis_element(j, isite1, isite3, tmp_V, tmp_v0, tmp_v1, X, &tmp_off);
      }
    } else if (isite1 == isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isite1,isite2,isite4,isite3,Asum,Bsum,Adiff,Bdiff,tmp_off,tmp_off_2,tmp_V) shared(tmp_v0, tmp_v1)
      for (j = 1; j <= i_max; j++) {
        dam_pr += GC_child_CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, tmp_V, tmp_v0, tmp_v1, X,
                                                &tmp_off);
      }
    } else if (isite1 != isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isite1,isite2,isite3,Asum,Adiff,tmp_off,tmp_V) shared(tmp_v0, tmp_v1)
      for (j = 1; j <= i_max; j++) {
        dam_pr += GC_child_CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, tmp_V, tmp_v0, tmp_v1, X,
                                                &tmp_off);
      }
    } else if (isite1 != isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isite1,isite2,isite3,isite4,Asum,Bsum,Adiff,Bdiff, tmp_off_2,tmp_V) shared(tmp_v0, tmp_v1)
      for (j = 1; j <= i_max; j++) {
        dam_pr += GC_child_CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, tmp_V,
                                                tmp_v0, tmp_v1, X, &tmp_off_2);
      }
    }

    return dam_pr;
  }

/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_pairhopp
          (
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X
          ) {
    long int j;
    long unsigned int i_max = X->Large.i_max;
    long unsigned int off = 0;
    double complex dam_pr = 0.0;

#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max,X,off) private(j) shared(tmp_v0, tmp_v1)
    for (j = 1; j <= i_max; j++) {
      dam_pr += GC_child_pairhopp_element(j, tmp_v0, tmp_v1, X, &off);
    }

    return dam_pr;
  }

/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_exchange
          (
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X
          ) {
    long int j;
    long unsigned int i_max = X->Large.i_max;
    long unsigned int off = 0;
    double complex dam_pr = 0.0;

#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max, X,off) private(j) shared(tmp_v0, tmp_v1)
    for (j = 1; j <= i_max; j++) {
      dam_pr += GC_child_exchange_element(j, tmp_v0, tmp_v1, X, &off);
    }
    return dam_pr;
  }

/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_exchange_spin
          (
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X
          ) {
    long unsigned int j;
    long unsigned int i_max = X->Large.i_max;
    long unsigned int off = 0;
    double complex dam_pr = 0;

#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max, X,off) private(j) shared(tmp_v0, tmp_v1)
    for (j = 1; j <= i_max; j++) {
      dam_pr += GC_child_exchange_spin_element(j, tmp_v0, tmp_v1, X, &off);
    }
    return dam_pr;
  }

/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_pairlift_spin
          (
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X
          ) {
    long int j;
    long unsigned int i_max = X->Large.i_max;
    long unsigned int off = 0;
    double complex dam_pr = 0;

#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max, X,off) private(j) shared(tmp_v0, tmp_v1)
    for (j = 1; j <= i_max; j++) {
      dam_pr += GC_child_pairlift_spin_element(j, tmp_v0, tmp_v1, X, &off);
    }
    return dam_pr;
  }

/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_general_int_spin(double complex *tmp_v0, double complex *tmp_v1, struct BindStruct *X) {
    double complex dam_pr, tmp_V, dmv;
    long unsigned int j, i_max;
    long unsigned int org_sigma2, org_sigma4;
    long unsigned int isA_up, isB_up;
    long unsigned int tmp_off = 0;
    int tmp_sgn;

    i_max = X->Large.i_max;
    org_sigma2 = X->Large.is2_spin;
    org_sigma4 = X->Large.is4_spin;
    tmp_V = X->Large.tmp_V;
    isA_up = X->Large.is1_up;
    isB_up = X->Large.is2_up;
    dam_pr = 0.0;

#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, dmv) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(tmp_v1, tmp_v0)
    for (j = 1; j <= i_max; j++) {
      tmp_sgn = X_child_exchange_spin_element(j, X, isA_up, isB_up, org_sigma2, org_sigma4, &tmp_off);
      if (tmp_sgn != 0) {
        dmv = tmp_v1[j] * tmp_sgn * tmp_V;
        tmp_v0[tmp_off] += dmv;
        dam_pr += conj(tmp_v1[tmp_off]) * dmv;
      }
    }

    return dam_pr;
  }

/**
 *
 *
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_general_int_spin(double complex *tmp_v0, double complex *tmp_v1, struct BindStruct *X) {
    double complex dam_pr, tmp_V;
    long unsigned int j, i_max;
    long unsigned int org_isite1, org_isite2;
    long unsigned int org_sigma1, org_sigma2, org_sigma3, org_sigma4;
    long unsigned int isA_up, isB_up;
    long unsigned int tmp_off = 0;

    i_max = X->Large.i_max;
    org_isite1 = X->Large.isite1;
    org_isite2 = X->Large.isite2;
    org_sigma1 = X->Large.is1_spin;
    org_sigma2 = X->Large.is2_spin;
    org_sigma3 = X->Large.is3_spin;
    org_sigma4 = X->Large.is4_spin;
    tmp_V = X->Large.tmp_V;
    dam_pr = 0.0;
    isA_up = X->Def.Tpow[org_isite1 - 1];
    isB_up = X->Def.Tpow[org_isite2 - 1];

    if (org_sigma1 == org_sigma2 && org_sigma3 == org_sigma4) { //diagonal
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off, tmp_V) shared(tmp_v0, tmp_v1)
      for (j = 1; j <= i_max; j++) {
        dam_pr += GC_child_CisAisCisAis_spin_element(j, isA_up, isB_up, org_sigma2, org_sigma4, tmp_V, tmp_v0, tmp_v1,
                                                     X);
      }
    }
    else if (org_sigma1 == org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(tmp_v0, tmp_v1)
      for (j = 1; j <= i_max; j++) {
        dam_pr += GC_child_CisAisCitAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, tmp_v0, tmp_v1,
                                                     X, &tmp_off);
      }
    } else if (org_sigma1 != org_sigma2 && org_sigma3 == org_sigma4) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(tmp_v0, tmp_v1)
      for (j = 1; j <= i_max; j++) {
        dam_pr += GC_child_CisAitCiuAiu_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, tmp_v0, tmp_v1,
                                                     X, &tmp_off);
      }
    } else if (org_sigma1 != org_sigma2 && org_sigma3 != org_sigma4) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max,X,isA_up,isB_up,org_sigma2,org_sigma4,tmp_off,tmp_V) shared(tmp_v0, tmp_v1)
      for (j = 1; j <= i_max; j++) {
        dam_pr += GC_child_CisAitCiuAiv_spin_element(j, org_sigma2, org_sigma4, isA_up, isB_up, tmp_V, tmp_v0, tmp_v1,
                                                     X, &tmp_off);
      }
    }

    return dam_pr;
  }
/******************************************************************************/
//[e] child functions
/******************************************************************************/


/******************************************************************************/
//[s] GetInfo functions
/******************************************************************************/
/**
 *
 *
 * @param X
 * @param isite1
 * @param isite2
 * @param sigma1
 * @param sigma2
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int child_general_hopp_GetInfo
          (
                  struct BindStruct *X,
                  unsigned long int isite1,
                  unsigned long int isite2,
                  unsigned long int sigma1,
                  unsigned long int sigma2
          ) {
    X->Large.is1_spin = X->Def.Tpow[2 * isite1 - 2 + sigma1];
    X->Large.is2_spin = X->Def.Tpow[2 * isite2 - 2 + sigma2];

    if (isite1 > isite2) {
      X->Large.A_spin = (X->Def.Tpow[2 * isite1 - 2 + sigma1] - X->Def.Tpow[2 * isite2 - 1 + sigma2]);
    } else if (isite1 < isite2) {
      X->Large.A_spin = (X->Def.Tpow[2 * isite2 - 2 + sigma2] - X->Def.Tpow[2 * isite1 - 1 + sigma1]);
    }
    else {
      if (sigma1 > sigma2) {
        X->Large.A_spin = (X->Def.Tpow[2 * isite1 - 2 + sigma1] - X->Def.Tpow[2 * isite2 - 1 + sigma2]);
      }
      else {
        X->Large.A_spin = (X->Def.Tpow[2 * isite2 - 2 + sigma2] - X->Def.Tpow[2 * isite1 - 1 + sigma1]);
      }
    }
    X->Large.isA_spin = X->Large.is1_spin + X->Large.is2_spin;
    return 0;
  }

/**
 *
 *
 * @param iInterAll
 * @param X
 * @param isite1
 * @param isite2
 * @param isite3
 * @param isite4
 * @param sigma1
 * @param sigma2
 * @param sigma3
 * @param sigma4
 * @param tmp_V
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int child_general_int_GetInfo
          (
                  const int iInterAll,
                  struct BindStruct *X,
                  long unsigned int isite1,
                  long unsigned int isite2,
                  long unsigned int isite3,
                  long unsigned int isite4,
                  long unsigned int sigma1,
                  long unsigned int sigma2,
                  long unsigned int sigma3,
                  long unsigned int sigma4,
                  double complex tmp_V
          ) {
    //int isite1, isite2, isite3, isite4;
    //int sigma1, sigma2, sigma3, sigma4;
    //double complex tmp_V;
    long unsigned int is1_spin, is2_spin, is3_spin, is4_spin;
    long unsigned int A_spin, B_spin;
    long unsigned int isA_spin, isB_spin;

    is1_spin = X->Def.Tpow[2 * isite1 - 2 + sigma1];
    is2_spin = X->Def.Tpow[2 * isite2 - 2 + sigma2];
    if (isite1 > isite2) {
      A_spin = (X->Def.Tpow[2 * isite1 - 2 + sigma1] - X->Def.Tpow[2 * isite2 - 1 + sigma2]);
    } else if (isite2 > isite1) {
      A_spin = (X->Def.Tpow[2 * isite2 - 2 + sigma2] - X->Def.Tpow[2 * isite1 - 1 + sigma1]);
    }
    else {//isite1=isite2
      if (sigma1 > sigma2) {
        A_spin = (X->Def.Tpow[2 * isite1 - 2 + sigma1] - X->Def.Tpow[2 * isite2 - 1 + sigma2]);
      }
      else {
        A_spin = (X->Def.Tpow[2 * isite2 - 2 + sigma2] - X->Def.Tpow[2 * isite1 - 1 + sigma1]);
      }
    }

    is3_spin = X->Def.Tpow[2 * isite3 - 2 + sigma3];
    is4_spin = X->Def.Tpow[2 * isite4 - 2 + sigma4];
    if (isite3 > isite4) {
      B_spin = (X->Def.Tpow[2 * isite3 - 2 + sigma3] - X->Def.Tpow[2 * isite4 - 1 + sigma4]);
    } else if (isite3 < isite4) {
      B_spin = (X->Def.Tpow[2 * isite4 - 2 + sigma4] - X->Def.Tpow[2 * isite3 - 1 + sigma3]);
    }
    else {//isite3=isite4
      if (sigma3 > sigma4) {
        B_spin = (X->Def.Tpow[2 * isite3 - 2 + sigma3] - X->Def.Tpow[2 * isite4 - 1 + sigma4]);
      }
      else {
        B_spin = (X->Def.Tpow[2 * isite4 - 2 + sigma4] - X->Def.Tpow[2 * isite3 - 1 + sigma3]);
      }
    }

    isA_spin = is1_spin + is2_spin;
    isB_spin = is3_spin + is4_spin;

    X->Large.is1_spin = is1_spin;
    X->Large.is2_spin = is2_spin;
    X->Large.is3_spin = is3_spin;
    X->Large.is4_spin = is4_spin;
    X->Large.isA_spin = isA_spin;
    X->Large.isB_spin = isB_spin;
    X->Large.A_spin = A_spin;
    X->Large.B_spin = B_spin;
    X->Large.tmp_V = tmp_V;
    X->Large.isite1 = isite1;
    X->Large.isite2 = isite2;
    X->Large.isite3 = isite3;
    X->Large.isite4 = isite4;

    return 0;
  }

/**
 *
 *
 * @param iPairHopp
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int child_pairhopp_GetInfo
          (
                  const int iPairHopp,
                  struct BindStruct *X
          ) {
    int isite1 = X->Def.PairHopping[iPairHopp][0] + 1;
    int isite2 = X->Def.PairHopping[iPairHopp][1] + 1;
    X->Large.tmp_J = X->Def.ParaPairHopping[iPairHopp];

    X->Large.is1_up = X->Def.Tpow[2 * isite1 - 2];
    X->Large.is1_down = X->Def.Tpow[2 * isite1 - 1];
    X->Large.is2_up = X->Def.Tpow[2 * isite2 - 2];
    X->Large.is2_down = X->Def.Tpow[2 * isite2 - 1];

    return 0;
  }

/**
 *
 *
 * @param iExchange
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int child_exchange_GetInfo
          (
                  const int iExchange,
                  struct BindStruct *X
          ) {
    int isite1 = X->Def.ExchangeCoupling[iExchange][0] + 1;
    int isite2 = X->Def.ExchangeCoupling[iExchange][1] + 1;
    X->Large.tmp_J = -X->Def.ParaExchangeCoupling[iExchange];

    X->Large.is1_up = X->Def.Tpow[2 * isite1 - 2];
    X->Large.is1_down = X->Def.Tpow[2 * isite1 - 1];
    X->Large.is2_up = X->Def.Tpow[2 * isite2 - 2];
    X->Large.is2_down = X->Def.Tpow[2 * isite2 - 1];

    return 0;
  }

/**
 *
 *
 * @param iExchange
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int child_exchange_spin_GetInfo
          (
                  const int iExchange,
                  struct BindStruct *X
          ) {
    int isite1 = X->Def.ExchangeCoupling[iExchange][0] + 1;
    int isite2 = X->Def.ExchangeCoupling[iExchange][1] + 1;
    X->Large.tmp_J = X->Def.ParaExchangeCoupling[iExchange];
    X->Large.is1_up = X->Def.Tpow[isite1 - 1];
    X->Large.is2_up = X->Def.Tpow[isite2 - 1];
    X->Large.isA_spin = X->Large.is1_up + X->Large.is2_up;
    return 0;
  }

/**
 *
 *
 * @param iPairLift
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int child_pairlift_spin_GetInfo
          (
                  const int iPairLift,
                  struct BindStruct *X
          ) {
    int isite1 = X->Def.PairLiftCoupling[iPairLift][0] + 1;
    int isite2 = X->Def.PairLiftCoupling[iPairLift][1] + 1;
    X->Large.tmp_J = X->Def.ParaPairLiftCoupling[iPairLift];
    X->Large.is1_up = X->Def.Tpow[isite1 - 1];
    X->Large.is2_up = X->Def.Tpow[isite2 - 1];
    X->Large.isA_spin = X->Large.is1_up + X->Large.is2_up;
    return 0;
  }

/**
 *
 *
 * @param X
 * @param isite1
 * @param isite2
 * @param sigma1
 * @param sigma2
 * @param sigma3
 * @param sigma4
 * @param tmp_V
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int child_general_int_spin_GetInfo
          (
                  struct BindStruct *X,
                  long unsigned int isite1,
                  long unsigned int isite2,
                  long unsigned int sigma1,
                  long unsigned int sigma2,
                  long unsigned int sigma3,
                  long unsigned int sigma4,
                  double complex tmp_V
          ) {
    X->Large.tmp_V = tmp_V;
    X->Large.isite1 = isite1;
    X->Large.isite2 = isite2;
    X->Large.is1_up = X->Def.Tpow[isite1 - 1];
    X->Large.is2_up = X->Def.Tpow[isite2 - 1];
    X->Large.is1_spin = sigma1;
    X->Large.is2_spin = sigma2;
    X->Large.is3_spin = sigma3;
    X->Large.is4_spin = sigma4;
    return 0;
  }

/******************************************************************************/
//[e] GetInfo functions
/******************************************************************************/


/******************************************************************************/
//[s] core routines
/******************************************************************************/
/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param is1_spin
 * @param tmp_trans
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_CisAis(
          long unsigned int j,
          double complex *tmp_v0,
          double complex *tmp_v1,
          struct BindStruct *X,
          long unsigned int is1_spin,
          double complex tmp_trans
  ) {

    long unsigned int A_ibit_tmp;
    long unsigned int list_1_j;
    double complex dmv;
    double complex dam_pr;

    // off = j-1;

    list_1_j = j - 1;
    A_ibit_tmp = (list_1_j & is1_spin) / is1_spin;
    dmv = tmp_v1[j] * A_ibit_tmp;
    if (X->Large.mode == M_MLTPLY) { // for multply
      tmp_v0[j] += dmv * tmp_trans;
    }
    dam_pr = dmv * conj(tmp_v1[j]);
    return dam_pr;
  }

/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param is1_spin
 * @param tmp_trans
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex CisAis(
          long unsigned int j,
          double complex *tmp_v0,
          double complex *tmp_v1,
          struct BindStruct *X,
          long unsigned int is1_spin,
          double complex tmp_trans
  ) {

    long unsigned int A_ibit_tmp;
    double complex dmv;
    double complex dam_pr;

    // off = j
    A_ibit_tmp = (list_1[j] & is1_spin) / is1_spin;
    dmv = tmp_v1[j] * A_ibit_tmp;
    if (X->Large.mode == M_MLTPLY) { // for multply
      tmp_v0[j] += dmv * tmp_trans;
    }
    dam_pr = dmv * conj(tmp_v1[j]);
    return dam_pr;

  }

/**
 *
 *
 * @param list_1_j
 * @param X
 * @param is1_spin
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int X_CisAis(
          long unsigned int list_1_j,
          struct BindStruct *X,
          long unsigned int is1_spin
  ) {

    int A_ibit_tmp;

    // off = j
    A_ibit_tmp = (list_1_j & is1_spin) / is1_spin;
    return A_ibit_tmp;

  }

/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param is1_spin
 * @param is2_spin
 * @param sum_spin
 * @param diff_spin
 * @param tmp_V
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex CisAjt(
          long unsigned int j,
          double complex *tmp_v0,
          double complex *tmp_v1,
          struct BindStruct *X,
          long unsigned int is1_spin,
          long unsigned int is2_spin,
          long unsigned int sum_spin,
          long unsigned int diff_spin,
          double complex tmp_V
  ) {
    long unsigned int ibit_tmp_1, ibit_tmp_2;
    long unsigned int bit, iexchg, off;
    int sgn;
    double complex dmv, dam_pr;

    ibit_tmp_1 = (list_1[j] & is1_spin);
    ibit_tmp_2 = (list_1[j] & is2_spin);
    if (ibit_tmp_1 == 0 && ibit_tmp_2 != 0) {
      bit = list_1[j] & diff_spin;
      SgnBit(bit, &sgn); // Fermion sign
      iexchg = list_1[j] ^ sum_spin;
      GetOffComp(list_2_1, list_2_2, iexchg, X->Large.irght, X->Large.ilft, X->Large.ihfbit, &off);

      dmv = sgn * tmp_v1[j];
      if (X->Large.mode == M_MLTPLY) { // for multply
        tmp_v0[off] += tmp_V * dmv;
      }
      dam_pr = dmv * conj(tmp_v1[off]);
      return dam_pr;
    } else {
      return 0;
    }
  }

/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param is1_spin
 * @param is2_spin
 * @param sum_spin
 * @param diff_spin
 * @param tmp_V
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_CisAjt(
          long unsigned int j,
          double complex *tmp_v0,
          double complex *tmp_v1,
          struct BindStruct *X,
          long unsigned int is1_spin,
          long unsigned int is2_spin,
          long unsigned int sum_spin,
          long unsigned int diff_spin,
          double complex tmp_V,
          long unsigned int *tmp_off
  ) {

    long unsigned int list_1_j, list_1_off;
    long unsigned int ibit_tmp_1, ibit_tmp_2;
    long unsigned int bit;
    int sgn;
    double complex dmv, dam_pr;

    list_1_j = j - 1;
    ibit_tmp_1 = (list_1_j & is1_spin);
    ibit_tmp_2 = (list_1_j & is2_spin);
    *tmp_off = 0;

    if (ibit_tmp_1 == 0 && ibit_tmp_2 != 0) {
      bit = list_1_j & diff_spin;
      SgnBit(bit, &sgn); // Fermion sign
      list_1_off = list_1_j ^ sum_spin;
      *tmp_off = list_1_off;
      dmv = sgn * tmp_v1[j];
      if (X->Large.mode == M_MLTPLY) { // for multply
        tmp_v0[list_1_off + 1] += dmv * tmp_V;
      }
      dam_pr = dmv * conj(tmp_v1[list_1_off + 1]);
      return dam_pr;
    } else {
      return 0;
    }
  }

/**
 *
 *
 * @param list_1_j
 * @param X
 * @param is1_spin
 * @param is2_spin
 * @param sum_spin
 * @param diff_spin
 * @param iexchg
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int X_CisAjt(
          long unsigned int list_1_j,
          struct BindStruct *X,
          long unsigned int is1_spin,
          long unsigned int is2_spin,
          long unsigned int sum_spin,
          long unsigned int diff_spin,
          long unsigned int *tmp_off
  ) {
    long unsigned int off;
    int sgn = 1;

    sgn = X_GC_CisAjt(list_1_j, X, is1_spin, is2_spin, sum_spin, diff_spin, tmp_off);
    if(sgn !=0){
      GetOffComp(list_2_1, list_2_2, *tmp_off, X->Large.irght, X->Large.ilft, X->Large.ihfbit, &off);
      *tmp_off =off;
      return sgn;
    }
    else {
      *tmp_off = 1;
      return 0;
    }
  }

/**
 *
 *
 * @param list_1_j
 * @param X
 * @param is1_spin
 * @param is2_spin
 * @param sum_spin
 * @param diff_spin
 * @param iexchg
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int X_GC_CisAjt(
          long unsigned int list_1_j,
          struct BindStruct *X,
          long unsigned int is1_spin,
          long unsigned int is2_spin,
          long unsigned int sum_spin,
          long unsigned int diff_spin,
          long unsigned int *tmp_off
  ) {
    long unsigned int ibit_tmp_1, ibit_tmp_2;
    long unsigned int bit, off;
    int sgn = 1;

    ibit_tmp_1 = (list_1_j & is1_spin);
    ibit_tmp_2 = (list_1_j & is2_spin);
    if (ibit_tmp_1 == 0 && ibit_tmp_2 != 0) {
      bit = list_1_j & diff_spin;
      SgnBit(bit, &sgn); // Fermion sign
      off = list_1_j ^ sum_spin;
      *tmp_off = off;
      return sgn; // pm 1
    } else {
      *tmp_off = 1;
      return 0;
    }
  }

// for Spin
/**
 *
 *
 * @param j
 * @param X
 * @param is1_spin
 * @param sigma1
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int X_Spin_CisAis(
          long unsigned int j,
          struct BindStruct *X,
          long unsigned int is1_spin,
          long unsigned int sigma1
  ) {
    int A_ibit_tmp;
    // off = j
    A_ibit_tmp = ((list_1[j] & is1_spin) / is1_spin) ^ (1 - sigma1);
    return A_ibit_tmp;
  }
//

/**
 *
 *
 * @param j
 * @param X
 * @param is1_spin
 * @param is2_spin
 * @param sigma1
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int X_Spin_CisAjs
          (
                  long unsigned int j,
                  struct BindStruct *X,
                  long unsigned int is1_spin,
                  long unsigned int is2_spin,
                  long unsigned int sigma1,
                  long unsigned int *tmp_off
          ) {
    int ibit_tmp;
    int jbit_tmp;
    long unsigned int iexchg;
    long unsigned int irght = X->Large.irght;
    long unsigned int ilft = X->Large.ilft;
    long unsigned int ihfbit = X->Large.ihfbit;
    long unsigned int is_up;
    is_up = is1_spin + is2_spin;

    //bit_tmp =0 -> down spin, 1 -> up spin
    ibit_tmp = ((list_1[j] & is1_spin) / is1_spin) ^ (1 - sigma1);//create: 1=(0^1, 1^0) is OK
    jbit_tmp = ((list_1[j] & is2_spin) / is2_spin) ^ (1 - sigma1);//anihilate: 0=(0^0, 1^1) is OK
    //1-sigma1 = 0 -> down spin, 1->up spin
    if (ibit_tmp != 0 && jbit_tmp == 0) {
      iexchg = list_1[j] ^ is_up;
      GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, tmp_off);
      return 1;
    }
    *tmp_off = 1;
    return 0;
  }

// for SpinGC
/**
 *
 *
 * @param j
 * @param X
 * @param is1_spin
 * @param sigma1
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int X_SpinGC_CisAis(
          long unsigned int j,
          struct BindStruct *X,
          long unsigned int is1_spin,
          long unsigned int sigma1
  ) {
    int A_ibit_tmp;
    long unsigned int list_1_j;
    // off = j
    list_1_j = j - 1;
    A_ibit_tmp = ((list_1_j & is1_spin) / is1_spin) ^ (1 - sigma1);
    return A_ibit_tmp;
  }

/**
 *
 *
 * @param j
 * @param X
 * @param is1_spin
 * @param sigma2
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int X_SpinGC_CisAit(
          long unsigned int j,
          struct BindStruct *X,
          long unsigned int is1_spin,
          long unsigned int sigma2,
          long unsigned int *tmp_off
  ) {
    long unsigned int list_1_j, ibit_tmp_1;

    list_1_j = j - 1;
    
    ibit_tmp_1 = list_1_j & is1_spin;
    if (ibit_tmp_1 == 0 && sigma2 == 0) {    // down -> up
      *tmp_off = list_1_j + is1_spin;
      return 1;
    } else if (ibit_tmp_1 != 0 && sigma2 == 1) { // up -> down
      *tmp_off = list_1_j - is1_spin;
      return 1;
    } else {
      *tmp_off = 1;
      return 0;
      }        
  }
/******************************************************************************/
//[e] core routines
/******************************************************************************/

/******************************************************************************/
//[s] child element functions
/******************************************************************************/
/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_exchange_element
          (
                  const long unsigned int j,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    long unsigned int off;
    long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
    double complex dmv;
    long unsigned int iexchg;
    long unsigned int is1_up = X->Large.is1_up;
    long unsigned int is2_up = X->Large.is2_up;
    long unsigned int is1_down = X->Large.is1_down;
    long unsigned int is2_down = X->Large.is2_down;
    long unsigned int irght = X->Large.irght;
    long unsigned int ilft = X->Large.ilft;
    long unsigned int ihfbit = X->Large.ihfbit;
    double complex tmp_J = X->Large.tmp_J;
    int mode = X->Large.mode;
    double complex dam_pr = 0;

    ibit1_up = list_1[j] & is1_up;
    ibit2_up = list_1[j] & is2_up;
    ibit1_down = list_1[j] & is1_down;
    ibit2_down = list_1[j] & is2_down;

    if (ibit1_up == 0 && ibit1_down != 0 && ibit2_up != 0 && ibit2_down == 0) {

      iexchg = list_1[j] - (is1_down + is2_up);
      iexchg += (is1_up + is2_down);
      GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
      *tmp_off = off;
      dmv = tmp_J * tmp_v1[j];
      if (mode == M_MLTPLY) {
        tmp_v0[off] += dmv;
      }
      dam_pr += dmv * conj(tmp_v1[off]);
    } else if (ibit1_up != 0 && ibit1_down == 0 && ibit2_up == 0 && ibit2_down != 0) {
      iexchg = list_1[j] - (is1_up + is2_down);
      iexchg += (is1_down + is2_up);
      GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
      *tmp_off = off;
      dmv = tmp_J * tmp_v1[j];
      if (mode == M_MLTPLY) {
        tmp_v0[off] += dmv;
      }
      dam_pr += dmv * conj(tmp_v1[off]);
    }

    return dam_pr;

  }

/**
 *
 *
 * @param j
 * @param X
 * @param isA_up
 * @param isB_up
 * @param sigmaA
 * @param sigmaB
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  int X_child_exchange_spin_element
          (
                  const long unsigned int j,
                  struct BindStruct *X,
                  const long unsigned int isA_up,
                  const long unsigned int isB_up,
                  const long unsigned int sigmaA,
                  const long unsigned int sigmaB,
                  long unsigned int *tmp_off
          ) {
    long unsigned int iexchg, off;
    long unsigned int irght = X->Large.irght;
    long unsigned int ilft = X->Large.ilft;
    long unsigned int ihfbit = X->Large.ihfbit;
    long unsigned int ibit_tmp_A, ibit_tmp_B;

    ibit_tmp_A = ((list_1[j] & isA_up) / isA_up);
    ibit_tmp_B = ((list_1[j] & isB_up) / isB_up);
    if (ibit_tmp_A == sigmaA && ibit_tmp_B == sigmaB) {
      iexchg = list_1[j] ^ (isA_up + isB_up);
      GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
      *tmp_off = off;
      return 1;
    } else {
      *tmp_off = 1; // just tentative
      return 0;
    }
  }


/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_exchange_spin_element
          (
                  const long unsigned int j,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    long unsigned int off;
    double complex dmv;
    long unsigned int iexchg;
    long unsigned int is_up = X->Large.isA_spin;
    long unsigned int irght = X->Large.irght;
    long unsigned int ilft = X->Large.ilft;
    long unsigned int ihfbit = X->Large.ihfbit;
    double complex tmp_J = X->Large.tmp_J;
    int mode = X->Large.mode;
    double complex dam_pr = 0;
    long unsigned int ibit_tmp;

    ibit_tmp = (list_1[j] & is_up);
    if (ibit_tmp == 0 || ibit_tmp == is_up) {
      return dam_pr;
    } else {
      iexchg = list_1[j] ^ is_up;
      GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
      *tmp_off = off;
      dmv = tmp_J * tmp_v1[j];
      if (mode == M_MLTPLY) {
        tmp_v0[off] += dmv;
      }
      dam_pr += dmv * conj(tmp_v1[off]);
      return dam_pr;
    }

  }

/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_exchange_spin_element
          (
                  const long unsigned int j,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    double complex dmv;
    long unsigned int is_up = X->Large.isA_spin;
    double complex tmp_J = X->Large.tmp_J;
    int mode = X->Large.mode;
    long unsigned int list_1_j, list_1_off;

    double complex dam_pr = 0;
    list_1_j = j - 1;

    long unsigned int ibit_tmp;
    ibit_tmp = (list_1_j & is_up);
    if (ibit_tmp == 0 || ibit_tmp == is_up) {
      return dam_pr;
    } else {
      list_1_off = list_1_j ^ is_up;
      *tmp_off = list_1_off;
      dmv = tmp_J * tmp_v1[j];
      if (mode == M_MLTPLY) {
        tmp_v0[list_1_off + 1] += dmv;
      }
      dam_pr += dmv * conj(tmp_v1[list_1_off + 1]);
      return dam_pr;
    }
  }

/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_pairlift_spin_element
          (
                  const long unsigned int j,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    long unsigned int off;
    double complex dmv;
    long unsigned int iexchg;
    long unsigned int is1_up = X->Large.is1_up;
    long unsigned int is2_up = X->Large.is2_up;
    long unsigned int is_up = X->Large.isA_spin;
    long unsigned int irght = X->Large.irght;
    long unsigned int ilft = X->Large.ilft;
    long unsigned int ihfbit = X->Large.ihfbit;
    double complex tmp_J = X->Large.tmp_J;
    int mode = X->Large.mode;
    double complex dam_pr = 0;

    long unsigned int ibit_tmp;
    ibit_tmp = ((list_1[j] & is1_up) / is1_up) ^ ((list_1[j] & is2_up) / is2_up);
    if (ibit_tmp == 0) {
      iexchg = list_1[j] ^ is_up; //Change: ++ -> -- or -- -> ++
      GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
      dmv = tmp_J * tmp_v1[j] * ibit_tmp;
      *tmp_off = off;
      if (mode == M_MLTPLY) {
        tmp_v0[off] += dmv;
      }
      dam_pr += dmv * conj(tmp_v1[off]);
    }
    return dam_pr;

  }

/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_pairlift_spin_element
          (
                  const long unsigned int j,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    double complex dmv;
    long unsigned int is_up = X->Large.isA_spin;
    double complex tmp_J = X->Large.tmp_J;
    int mode = X->Large.mode;
    double complex dam_pr = 0;
    long unsigned int list_1_off;
    long unsigned int list_1_j = j - 1;
    long unsigned int ibit_tmp;
    //ibit_tmp = ((list_1_j & is1_up) / is1_up) ^ ((list_1_j & is2_up) / is2_up);
    ibit_tmp = (list_1_j & is_up);
    if (ibit_tmp == 0|| ibit_tmp==is_up) {
      list_1_off = list_1_j ^ is_up; //Change: ++ -> -- or -- -> ++
      *tmp_off = list_1_off;
      dmv = tmp_J * tmp_v1[j] ;//* ibit_tmp;
      if (mode == M_MLTPLY) {
        tmp_v0[list_1_off + 1] += dmv;
      }
      dam_pr += dmv * conj(tmp_v1[list_1_off + 1]);
      return dam_pr;
    }else{
      return dam_pr;
    } 
  }

/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_pairhopp_element
          (
                  const long unsigned int j,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    long unsigned int off;
    long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
    double complex dmv;
    long unsigned int iexchg;
    long unsigned int is1_up = X->Large.is1_up;
    long unsigned int is2_up = X->Large.is2_up;
    long unsigned int is1_down = X->Large.is1_down;
    long unsigned int is2_down = X->Large.is2_down;
    long unsigned int irght = X->Large.irght;
    long unsigned int ilft = X->Large.ilft;
    long unsigned int ihfbit = X->Large.ihfbit;
    double complex tmp_J = X->Large.tmp_J;
    int mode = X->Large.mode;
    double complex dam_pr = 0;

    ibit1_up = list_1[j] & is1_up;
    ibit2_up = list_1[j] & is2_up;
    ibit1_down = list_1[j] & is1_down;
    ibit2_down = list_1[j] & is2_down;

    if (ibit1_up == 0 && ibit1_down == 0 && ibit2_up != 0 && ibit2_down != 0) {
      iexchg = list_1[j] - (is2_up + is2_down);
      iexchg += (is1_up + is1_down);
      GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
      *tmp_off = off;
      dmv = tmp_J * tmp_v1[j];
      if (mode == M_MLTPLY) {
        tmp_v0[off] += dmv;
      }
      dam_pr += dmv * conj(tmp_v1[off]);
    }
    return dam_pr;
  }

/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_exchange_element
          (
                  const long unsigned int j,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
    double complex dmv;
    long unsigned int iexchg;
    long unsigned int is1_up = X->Large.is1_up;
    long unsigned int is2_up = X->Large.is2_up;
    long unsigned int is1_down = X->Large.is1_down;
    long unsigned int is2_down = X->Large.is2_down;
    long unsigned int list_1_j, list_1_off;
    double complex tmp_J = X->Large.tmp_J;
    int mode = X->Large.mode;
    double complex dam_pr = 0;

    list_1_j = j - 1;
    ibit1_up = list_1_j & is1_up;
    ibit2_up = list_1_j & is2_up;
    ibit1_down = list_1_j & is1_down;
    ibit2_down = list_1_j & is2_down;

    if (ibit1_up == 0 && ibit1_down != 0 && ibit2_up != 0 && ibit2_down == 0) {

      iexchg = list_1_j - (is1_down + is2_up);
      iexchg += (is1_up + is2_down);
      list_1_off = iexchg;
      *tmp_off = list_1_off;

      dmv = tmp_J * tmp_v1[j];
      if (mode == M_MLTPLY) {
        tmp_v0[list_1_off + 1] += dmv;
      }
      dam_pr += dmv * conj(tmp_v1[list_1_off + 1]);
    } else if (ibit1_up != 0 && ibit1_down == 0 && ibit2_up == 0 && ibit2_down != 0) {
      iexchg = list_1_j - (is1_up + is2_down);
      iexchg += (is1_down + is2_up);
      list_1_off = iexchg;
      *tmp_off = list_1_off;

      dmv = tmp_J * tmp_v1[j];
      if (mode == M_MLTPLY) {
        tmp_v0[list_1_off + 1] += dmv;
      }
      dam_pr += dmv * conj(tmp_v1[list_1_off + 1]);
    }

    return dam_pr;
  }

/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_pairhopp_element
          (
                  const long unsigned int j,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
    double complex dmv;
    long unsigned int iexchg;
    long unsigned int is1_up = X->Large.is1_up;
    long unsigned int is2_up = X->Large.is2_up;
    long unsigned int is1_down = X->Large.is1_down;
    long unsigned int is2_down = X->Large.is2_down;
    long unsigned int list_1_j, list_1_off;
    double complex tmp_J = X->Large.tmp_J;
    int mode = X->Large.mode;

    double complex dam_pr = 0 + 0 * I;
    list_1_j = j - 1;

    ibit1_up = list_1_j & is1_up;

    ibit2_up = list_1_j & is2_up;

    ibit1_down = list_1_j & is1_down;

    ibit2_down = list_1_j & is2_down;

    if (ibit1_up == 0 && ibit1_down == 0 && ibit2_up != 0 && ibit2_down != 0) {
      iexchg = list_1_j - (is2_up + is2_down);
      iexchg += (is1_up + is1_down);
      list_1_off = iexchg;
      *tmp_off = list_1_off;
      dmv = tmp_J * tmp_v1[j];
      if (mode == M_MLTPLY) {
        tmp_v0[list_1_off + 1] += dmv;
      }
      dam_pr += dmv * conj(tmp_v1[list_1_off + 1]);
    }

    return dam_pr;
  }

/**
 *
 *
 * @param j
 * @param isite1
 * @param isite3
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_CisAisCisAis_element
          (
                  long unsigned int j,
                  long unsigned int isite1,
                  long unsigned int isite3,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    int tmp_sgn;
    double complex dmv;
    double complex dam_pr = 0 + 0 * I;
    tmp_sgn = X_CisAis(list_1[j], X, isite3);
    tmp_sgn *= X_CisAis(list_1[j], X, isite1);
    dmv = tmp_V * tmp_v1[j] * tmp_sgn;
    if (X->Large.mode == M_MLTPLY) { // for multply
      tmp_v0[j] += dmv;
    }
    dam_pr = conj(tmp_v1[j]) * dmv;
    return dam_pr;
  }

/**
 *
 *
 * @param j
 * @param isite1
 * @param isite3
 * @param isite4
 * @param Bsum
 * @param Bdiff
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_CisAisCjtAku_element
          (
                  long unsigned int j,
                  long unsigned int isite1,
                  long unsigned int isite3,
                  long unsigned int isite4,
                  long unsigned int Bsum,
                  long unsigned int Bdiff,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    int tmp_sgn;
    double complex dmv;
    double complex dam_pr = 0 + 0 * I;
    tmp_sgn = X_CisAjt(list_1[j], X, isite3, isite4, Bsum, Bdiff, tmp_off);
    if (tmp_sgn != 0) {
      tmp_sgn *= X_CisAis(list_1[*tmp_off], X, isite1);
      if (tmp_sgn != 0) {
        dmv = tmp_V * tmp_v1[j] * tmp_sgn;
        if (X->Large.mode == M_MLTPLY) { // for multply
          tmp_v0[*tmp_off] += dmv;
        }
        dam_pr = conj(tmp_v1[*tmp_off]) * dmv;
      }
    }
    return dam_pr;
  }

/**
 *
 *
 * @param j
 * @param isite1
 * @param isite2
 * @param isite3
 * @param Asum
 * @param Adiff
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_CisAjtCkuAku_element
          (
                  long unsigned int j,
                  long unsigned int isite1,
                  long unsigned int isite2,
                  long unsigned int isite3,
                  long unsigned int Asum,
                  long unsigned int Adiff,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    int tmp_sgn;
    double complex dmv;
    double complex dam_pr;
    dam_pr = 0;
    tmp_sgn = X_CisAis(list_1[j], X, isite3);
    if (tmp_sgn != 0) {
      tmp_sgn *= X_CisAjt(list_1[j], X, isite1, isite2, Asum, Adiff, tmp_off);
      if (tmp_sgn != 0) {
        dmv = tmp_V * tmp_v1[j] * tmp_sgn;
        if (X->Large.mode == M_MLTPLY) { // for multply
          tmp_v0[*tmp_off] += dmv;
        }
        dam_pr = conj(tmp_v1[*tmp_off]) * dmv;
      }
    }
    return dam_pr;
  }

/**
 *
 *
 * @param j
 * @param isite1
 * @param isite2
 * @param isite3
 * @param isite4
 * @param Asum
 * @param Adiff
 * @param Bsum
 * @param Bdiff
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off_2
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_CisAjtCkuAlv_element
          (
                  long unsigned int j,
                  long unsigned int isite1,
                  long unsigned int isite2,
                  long unsigned int isite3,
                  long unsigned int isite4,
                  long unsigned int Asum,
                  long unsigned int Adiff,
                  long unsigned int Bsum,
                  long unsigned int Bdiff,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off_2
          ) {
    int tmp_sgn;
    long unsigned int tmp_off_1;

    double complex dmv;
    double complex dam_pr = 0;
    tmp_sgn = X_GC_CisAjt(list_1[j], X, isite3, isite4, Bsum, Bdiff,  &tmp_off_1);

    if (tmp_sgn != 0) {
      tmp_sgn *= X_CisAjt(tmp_off_1, X, isite1, isite2, Asum, Adiff, tmp_off_2);       
      if (tmp_sgn != 0) {
        dmv = tmp_V * tmp_v1[j] * tmp_sgn;
        if (X->Large.mode == M_MLTPLY) { // for multply
          tmp_v0[*tmp_off_2] += dmv;
        }
        dam_pr = conj(tmp_v1[*tmp_off_2]) * dmv;
      }
    }
    return dam_pr;
  }

//[s]Spin
/**
 *
 *
 * @param j
 * @param isA_up
 * @param isB_up
 * @param org_sigma2
 * @param org_sigma4
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_CisAisCisAis_spin_element
          (
                  long unsigned int j,
                  long unsigned int isA_up,
                  long unsigned int isB_up,
                  long unsigned int org_sigma2,
                  long unsigned int org_sigma4,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X
          ) {
    int tmp_sgn;
    double complex dmv;
    double complex dam_pr = 0;

    tmp_sgn = X_Spin_CisAis(j, X, isB_up, org_sigma4);
    tmp_sgn *= X_Spin_CisAis(j, X, isA_up, org_sigma2);
    dmv = tmp_v1[j] * tmp_sgn * tmp_V;
    if (X->Large.mode == M_MLTPLY) { // for multply
      tmp_v0[j] += dmv;
    }
    dam_pr = conj(tmp_v1[j]) * dmv;
    return dam_pr;

  }

/**
 *
 *
 * @param j
 * @param isA_up
 * @param isB_up
 * @param org_sigma2
 * @param org_sigma4
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_CisAjsCjtAit_spin_element
          (
                  long unsigned int j,
                  long unsigned int isA_up,
                  long unsigned int isB_up,
                  long unsigned int org_sigma2,
                  long unsigned int org_sigma4,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    int tmp_sgn;
    long unsigned int tmp_off_1;
    double complex dmv = 0;
    double complex dam_pr = 0;
    tmp_sgn = X_Spin_CisAjs(j, X, isB_up, isA_up, org_sigma4, &tmp_off_1);
    if (tmp_sgn != 0) {
      tmp_sgn *= X_Spin_CisAjs(tmp_off_1, X, isA_up, isB_up, org_sigma2, tmp_off);
      if (tmp_sgn != 0) {
        dmv = tmp_v1[j] * tmp_V;
        if (X->Large.mode == M_MLTPLY) { // for multply
          tmp_v0[*tmp_off] += dmv;
        }
        dam_pr = conj(tmp_v1[*tmp_off]) * dmv;
      }
    }
    return dam_pr;
  }
//[e]Spin

//[s]GC Spin
/**
 *
 *
 * @param j
 * @param isA_up
 * @param isB_up
 * @param org_sigma2
 * @param org_sigma4
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_CisAisCisAis_spin_element
          (
                  long unsigned int j,
                  long unsigned int isA_up,
                  long unsigned int isB_up,
                  long unsigned int org_sigma2,
                  long unsigned int org_sigma4,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X
          ) {
    int tmp_sgn;
    double complex dmv = 0;
    double complex dam_pr = 0;

    tmp_sgn = X_SpinGC_CisAis(j, X, isB_up, org_sigma4);
    tmp_sgn *= X_SpinGC_CisAis(j, X, isA_up, org_sigma2);
    if (tmp_sgn != 0) {
      dmv = tmp_v1[j] * tmp_sgn * tmp_V;
      if (X->Large.mode == M_MLTPLY) { // for multply
        tmp_v0[j] += dmv;
      }
      dam_pr = conj(tmp_v1[j]) * dmv;
    }
    return dam_pr;

  }

/**
 *
 *
 * @param j
 * @param org_sigma2
 * @param org_sigma4
 * @param isA_up
 * @param isB_up
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_CisAisCitAiu_spin_element
          (
                  long unsigned int j,
                  long unsigned int org_sigma2,
                  long unsigned int org_sigma4,
                  long unsigned int isA_up,
                  long unsigned int isB_up,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    int tmp_sgn;
    double complex dmv;
    double complex dam_pr = 0 + 0 * I;
    tmp_sgn = X_SpinGC_CisAit(j, X, isB_up, org_sigma4, tmp_off);
    if (tmp_sgn != 0) {
      tmp_sgn *= X_SpinGC_CisAis((*tmp_off + 1), X, isA_up, org_sigma2);
      if (tmp_sgn != 0) {
        dmv = tmp_v1[j] * tmp_sgn * tmp_V;
        if (X->Large.mode == M_MLTPLY) { // for multply
          tmp_v0[*tmp_off + 1] += dmv;
        }
        dam_pr = conj(tmp_v1[*tmp_off + 1]) * dmv;
      }
    }
    return dam_pr;
  }

/**
 *
 *
 * @param j
 * @param org_sigma2
 * @param org_sigma4
 * @param isA_up
 * @param isB_up
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_CisAitCiuAiu_spin_element
          (
                  long unsigned int j,
                  long unsigned int org_sigma2,
                  long unsigned int org_sigma4,
                  long unsigned int isA_up,
                  long unsigned int isB_up,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    int tmp_sgn;
    double complex dmv;
    double complex dam_pr = 0 + 0 * I;
    tmp_sgn = X_SpinGC_CisAis(j, X, isB_up, org_sigma4);
    if (tmp_sgn != 0) {
      tmp_sgn *= X_SpinGC_CisAit(j, X, isA_up, org_sigma2, tmp_off);
      if (tmp_sgn != 0) {
        dmv = tmp_v1[j] * tmp_sgn * tmp_V;
        if (X->Large.mode == M_MLTPLY) { // for multply
          tmp_v0[*tmp_off + 1] += dmv;
        }
        dam_pr = conj(tmp_v1[*tmp_off + 1]) * dmv;
      }
    }
    return dam_pr;
  }

/**
 *
 *
 * @param j
 * @param org_sigma2
 * @param org_sigma4
 * @param isA_up
 * @param isB_up
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off_2
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_CisAitCiuAiv_spin_element
          (
                  long unsigned int j,
                  long unsigned int org_sigma2,
                  long unsigned int org_sigma4,
                  long unsigned int isA_up,
                  long unsigned int isB_up,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off_2
          ) {
    int tmp_sgn;
    long unsigned int tmp_off_1;
    double complex dmv;
    double complex dam_pr = 0 + 0 * I;
    tmp_sgn = X_SpinGC_CisAit(j, X, isB_up, org_sigma4, &tmp_off_1);
    if (tmp_sgn != 0) {
      tmp_sgn *= X_SpinGC_CisAit((tmp_off_1 + 1), X, isA_up, org_sigma2, tmp_off_2);
      if (tmp_sgn != 0) {
        dmv = tmp_v1[j] * tmp_sgn * tmp_V;
        if (X->Large.mode == M_MLTPLY) { // for multply
          tmp_v0[*tmp_off_2 + 1] += dmv;
        }
        dam_pr = conj(tmp_v1[*tmp_off_2 + 1]) * dmv;
      }
    }
    return dam_pr;
  }
//[e]GC Spin

//[s] Grand Canonical
/**
 *
 *
 * @param j
 * @param isite1
 * @param isite3
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_CisAisCisAis_element
          (
                  long unsigned int j,
                  long unsigned int isite1,
                  long unsigned int isite3,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    int tmp_sgn;
    double complex dmv;
    double complex dam_pr = 0;
    tmp_sgn = X_CisAis(j - 1, X, isite3);
    tmp_sgn *= X_CisAis(j - 1, X, isite1);
    if (tmp_sgn != 0) {
      dmv = tmp_V * tmp_v1[j] * tmp_sgn;
      if (X->Large.mode == M_MLTPLY) { // for multply
        tmp_v0[j] += dmv;
      }
      dam_pr = conj(tmp_v1[j]) * dmv;
    }
    return dam_pr;
  }

/**
 *
 *
 * @param j
 * @param isite1
 * @param isite3
 * @param isite4
 * @param Bsum
 * @param Bdiff
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_CisAisCjtAku_element
          (
                  long unsigned int j,
                  long unsigned int isite1,
                  long unsigned int isite3,
                  long unsigned int isite4,
                  long unsigned int Bsum,
                  long unsigned int Bdiff,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    int tmp_sgn;
    double complex dmv;
    double complex dam_pr = 0 + 0 * I;
    tmp_sgn = X_GC_CisAjt((j - 1), X, isite3, isite4, Bsum, Bdiff, tmp_off);
    if (tmp_sgn != 0) {
      tmp_sgn *= X_CisAis(*tmp_off, X, isite1);
      if (tmp_sgn != 0) {
        dmv = tmp_V * tmp_v1[j] * tmp_sgn;
        if (X->Large.mode == M_MLTPLY) { // for multply
          tmp_v0[*tmp_off + 1] += dmv;
        }
        dam_pr = conj(tmp_v1[*tmp_off + 1]) * dmv;
      }
    }
    return dam_pr;
  }

/**
 *
 *
 * @param j
 * @param isite1
 * @param isite2
 * @param isite3
 * @param Asum
 * @param Adiff
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_CisAjtCkuAku_element
          (
                  long unsigned int j,
                  long unsigned int isite1,
                  long unsigned int isite2,
                  long unsigned int isite3,
                  long unsigned int Asum,
                  long unsigned int Adiff,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off
          ) {
    int tmp_sgn;
    double complex dmv;
    double complex dam_pr = 0 + 0 * I;
    tmp_sgn = X_CisAis((j - 1), X, isite3);
    if (tmp_sgn != 0) {
      tmp_sgn *= X_GC_CisAjt((j - 1), X, isite1, isite2, Asum, Adiff, tmp_off);
      if (tmp_sgn != 0) {
        dmv = tmp_V * tmp_v1[j] * tmp_sgn;
        if (X->Large.mode == M_MLTPLY) { // for multply
          tmp_v0[*tmp_off + 1] += dmv;
        }
        dam_pr = conj(tmp_v1[*tmp_off + 1]) * dmv;
      }
    }
    return dam_pr;
  }

/**
 *
 *
 * @param j
 * @param isite1
 * @param isite2
 * @param isite3
 * @param isite4
 * @param Asum
 * @param Adiff
 * @param Bsum
 * @param Bdiff
 * @param tmp_V
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param tmp_off_2
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex GC_child_CisAjtCkuAlv_element
          (
                  long unsigned int j,
                  long unsigned int isite1,
                  long unsigned int isite2,
                  long unsigned int isite3,
                  long unsigned int isite4,
                  long unsigned int Asum,
                  long unsigned int Adiff,
                  long unsigned int Bsum,
                  long unsigned int Bdiff,
                  double complex tmp_V,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X,
                  long unsigned int *tmp_off_2
          ) {
    int tmp_sgn;
    long unsigned int tmp_off_1;
    double complex dmv;
    double complex dam_pr = 0 + 0 * I;

    tmp_sgn = X_GC_CisAjt((j - 1), X, isite3, isite4, Bsum, Bdiff, &tmp_off_1);
    if (tmp_sgn != 0) {
      tmp_sgn *= X_GC_CisAjt(tmp_off_1, X, isite1, isite2, Asum, Adiff, tmp_off_2);
      if (tmp_sgn != 0) {
        dmv = tmp_V * tmp_v1[j] * tmp_sgn;
        if (X->Large.mode == M_MLTPLY) { // for multply
          tmp_v0[*tmp_off_2 + 1] += dmv;
        }
        dam_pr = conj(tmp_v1[*tmp_off_2 + 1]) * dmv;
      }
    }
    return dam_pr;
  }
//[e] Grand Canonical

//[s] for debug
/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 *
 * @return
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
  double complex child_general_hopp_element
          (
                  const long unsigned int j,
                  double complex *tmp_v0,
                  double complex *tmp_v1,
                  struct BindStruct *X
          ) {
    long unsigned int bit;
    int sgn;
    long unsigned int iexchg;
    long unsigned int off;
    double complex dam_pr, dmv;
    long unsigned int ibit_1, ibit_2;
    long unsigned int is1 = X->Large.is1_spin;
    long unsigned int is2 = X->Large.is2_spin;
    long unsigned int is = X->Large.isA_spin;
    long unsigned int b_sigma = X->Large.A_spin;
    long unsigned int irght = X->Large.irght;
    long unsigned int ilft = X->Large.ilft;
    long unsigned int ihfbit = X->Large.ihfbit;
    double complex trans = X->Large.tmp_trans;


    ibit_1 = list_1[j] & is1;
    ibit_2 = list_1[j] & is2;
    dam_pr = 0;
    if (ibit_1 == 0 && ibit_2 != 0) {
      bit = list_1[j] & b_sigma;
      SgnBit(bit, &sgn); // Fermion sign
      //SplitBit(bit,irght, ilft, ihfbit, &b_sigma_rght, &b_sigma_lft);
      //sgn=list_3[b_sigma_rght]*list_3[b_sigma_lft];

      iexchg = list_1[j] ^ is;
      GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
      dmv = (sgn) * tmp_v1[off] * trans;
      if (X->Large.mode == M_MLTPLY) { // for multply
        tmp_v0[j] += dmv;
      }
      dam_pr += dmv * conj(tmp_v1[j]);
    }
    return dam_pr;
  }
//[e] for debug

/******************************************************************************/
//[e] child element functions
/******************************************************************************/
