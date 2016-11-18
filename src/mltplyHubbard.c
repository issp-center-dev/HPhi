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


//Define Mode for mltply
// complex version
#include <bitcalc.h>
#include "mfmemory.h"
#include "xsetmem.h"
#include "mltply.h"
#include "mltplyHubbard.h"
#include "mltplyMPI.h"
#include "wrapperMPI.h"
#include "CalcTime.h"
/**

/**
 *
 *
 * @param X
 * @param tmp_v0
 * @param tmp_v1
 *
 * @return
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int mltplyHubbard(struct BindStruct *X, double complex *tmp_v0,double complex *tmp_v1){
  
  long unsigned int i;
  long unsigned int isite1, isite2, sigma1, sigma2;
  long unsigned int isite3, isite4, sigma3, sigma4;
  long unsigned int ibitsite1, ibitsite2, ibitsite3, ibitsite4;

  double complex dam_pr;
  double complex tmp_trans;
  /*[s] For InterAll */
  double complex tmp_V;
  /*[e] For InterAll */

  int ihermite=0;
  int idx=0;

  StartTimer(300);
  StartTimer(310);
  //Transfer
  for (i = 0; i < X->Def.EDNTransfer; i+=2) {
    if (X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite &&
        X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite) {
      StartTimer(311);
      child_general_hopp_MPIdouble(i, X, tmp_v0, tmp_v1);
      StopTimer(311);
    }
    else if (X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite) {
      StartTimer(312);
      child_general_hopp_MPIsingle(i, X, tmp_v0, tmp_v1);
      StopTimer(312);
    }
    else if (X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite) {
      StartTimer(312);
      child_general_hopp_MPIsingle(i + 1, X, tmp_v0, tmp_v1);
      StopTimer(312);
    }
    else {
      StartTimer(313);
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
      StopTimer(313);
    }
  }
  StopTimer(310);
  StartTimer(320);
      
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
      StartTimer(321);
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
      StopTimer(321);
	}
	else{
      StartTimer(322);
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
      StopTimer(322); 
	}
	X->Large.prdct += dam_pr;
  }
  StopTimer(320);
  StartTimer(330);
      
  //Pair hopping
  for (i = 0; i < X->Def.NPairHopping; i +=2) {
    sigma1=0;
    sigma2=1;
    dam_pr = 0.0;
        
    if (
	    X->Def.PairHopping[i][0] + 1 > X->Def.Nsite ||
	    X->Def.PairHopping[i][1] + 1 > X->Def.Nsite
	    )
	  {
        StartTimer(331);
	    dam_pr = X_child_CisAjtCkuAlv_Hubbard_MPI
          (
           X->Def.PairHopping[i][0], sigma1,
           X->Def.PairHopping[i][1], sigma1,
           X->Def.PairHopping[i][0], sigma2,
           X->Def.PairHopping[i][1], sigma2,
           X->Def.ParaPairHopping[i], X, tmp_v0, tmp_v1
           );
        StopTimer(331);
	  }
    else {
      StartTimer(332);
      for (ihermite = 0; ihermite<2; ihermite++) {
        idx = i + ihermite;
        child_pairhopp_GetInfo(idx, X);
        dam_pr += child_pairhopp(tmp_v0, tmp_v1, X);            
      }/*for (ihermite = 0; ihermite<2; ihermite++)*/
      StopTimer(332);
    }
	X->Large.prdct += dam_pr;
  }/*for (i = 0; i < X->Def.NPairHopping; i += 2)*/
  StopTimer(330);
      
  StartTimer(340);
  //Exchange
  for (i = 0; i < X->Def.NExchangeCoupling; i ++) {
    sigma1=0; sigma2=1;
    dam_pr=0.0;
    if (X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite ||
        X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite) {
      StartTimer(341);
      dam_pr = X_child_CisAjtCkuAlv_Hubbard_MPI
        (
         X->Def.ExchangeCoupling[i][0], sigma1,
         X->Def.ExchangeCoupling[i][1], sigma1,
         X->Def.ExchangeCoupling[i][1], sigma2,
         X->Def.ExchangeCoupling[i][0], sigma2,
         X->Def.ParaExchangeCoupling[i], X, tmp_v0, tmp_v1
         );
      StopTimer(341);
    }
    else {
      StartTimer(342);
      child_exchange_GetInfo(i, X);
      dam_pr = child_exchange(tmp_v0, tmp_v1, X);
      StopTimer(342);
    }
    X->Large.prdct += dam_pr;
  }/*for (i = 0; i < X->Def.NExchangeCoupling; i ++)*/
  StopTimer(340);
  StopTimer(300);

  
  return 0;
}

int mltplyHubbardGC(struct BindStruct *X, double complex *tmp_v0,double complex *tmp_v1){

  long unsigned int i;
  long unsigned int isite1, isite2, sigma1, sigma2;
  long unsigned int isite3, isite4, sigma3, sigma4;
  long unsigned int ibitsite1, ibitsite2, ibitsite3, ibitsite4;

  double complex dam_pr;
  double complex tmp_trans;
  /*[s] For InterAll */
  double complex tmp_V;
  /*[e] For InterAll */

  int ihermite=0;
  int idx=0;

  //Transfer
  StartTimer(200);
  StartTimer(210);
  for (i = 0; i < X->Def.EDNTransfer; i += 2) {

    if (X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite &&
        X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite) {
      StartTimer(211);
      GC_child_general_hopp_MPIdouble(i, X, tmp_v0, tmp_v1);
      StopTimer(211);
    }
    else if (X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite){
      StartTimer(212);
      GC_child_general_hopp_MPIsingle(i, X, tmp_v0, tmp_v1);
      StopTimer(212);
    }
    else if (X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite) {
      StartTimer(212);
      GC_child_general_hopp_MPIsingle(i+1, X, tmp_v0, tmp_v1);
      StopTimer(212);
    }
    else {
      StartTimer(213);
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
      StopTimer(213);
    }
  }
  StopTimer(210);

  StartTimer(220);
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
      StartTimer(221);
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
      StopTimer(221);
    }//InterPE
    else{
      StartTimer(222);
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
      StopTimer(222);
    }
	X->Large.prdct += dam_pr;
  }

  StopTimer(220);
  StartTimer(230);
  //Pair hopping
  for (i = 0; i < X->Def.NPairHopping; i +=2) {
	sigma1=0;
	sigma2=1;
	dam_pr=0.0;
    if (X->Def.PairHopping[i][0] + 1 > X->Def.Nsite ||
	    X->Def.PairHopping[i][1] + 1 > X->Def.Nsite) {
      StartTimer(231);
	  dam_pr = X_GC_child_CisAjtCkuAlv_Hubbard_MPI
	    (
	     X->Def.PairHopping[i][0], sigma1,
	     X->Def.PairHopping[i][1], sigma1,
	     X->Def.PairHopping[i][0], sigma2,
	     X->Def.PairHopping[i][1], sigma2,
	     X->Def.ParaPairHopping[i], X, tmp_v0, tmp_v1
	     );
      StopTimer(231);
    }
    else {
      StartTimer(232);
      for (ihermite = 0; ihermite<2; ihermite++) {
        idx = i + ihermite;
	    child_pairhopp_GetInfo(idx, X);
	    dam_pr += GC_child_pairhopp(tmp_v0, tmp_v1, X);
	  }
      StopTimer(232);
	}
	X->Large.prdct += dam_pr;
  }/*for (i = 0; i < X->Def.NPairHopping; i += 2)*/
      
  StopTimer(230);
  StartTimer(240);
  //Exchange
  for (i = 0; i < X->Def.NExchangeCoupling; i++) {
    sigma1=0; sigma2=1;
    dam_pr=0.0;
    if(X->Def.ExchangeCoupling[i][0] + 1 > X->Def.Nsite ||
       X->Def.ExchangeCoupling[i][1] + 1 > X->Def.Nsite){
      StartTimer(241);
      dam_pr = X_GC_child_CisAjtCkuAlv_Hubbard_MPI
	    (
	     X->Def.ExchangeCoupling[i][0], sigma1,
	     X->Def.ExchangeCoupling[i][1], sigma1,
	     X->Def.ExchangeCoupling[i][1], sigma2,
	     X->Def.ExchangeCoupling[i][0], sigma2,
	     X->Def.ParaExchangeCoupling[i], X, tmp_v0, tmp_v1
	     );
      StopTimer(241);
    }
    else {
      StartTimer(242);
	  child_exchange_GetInfo(i, X);
	  dam_pr = GC_child_exchange(tmp_v0, tmp_v1, X);
      StopTimer(242);
    }
	X->Large.prdct += dam_pr;
  }/*for (i = 0; i < X->Def.NExchangeCoupling; i++)*/
  StopTimer(240);
  StopTimer(200);
  return 0;
}

