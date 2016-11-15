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

#include "bitcalc.h"
#include "mltply.h"
#include "expec_energy_flct.h"
#include "wrapperMPI.h"
#include "CalcTime.h"

/** 
 * 
 * 
 * @param X 
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return 
 */
int expec_energy_flct(struct BindStruct *X){

  long unsigned int i,j;
  long unsigned int irght,ilft,ihfbit;
  long unsigned int isite1;
  long unsigned int is1_up,is1_down;
  long unsigned int is1;
  double complex dam_pr,dam_pr1;

  long int num1_up, num1_down;
  long unsigned int ibit1;
  double tmp_num_up, tmp_num_down;
  double D,tmp_D,tmp_D2;
  double N,tmp_N,tmp_N2;
  double Sz,tmp_Sz, tmp_Sz2;
  double tmp_v02;  
  long unsigned int i_max,tmp_list_1;
  
  switch(X->Def.iCalcType){
  case Lanczos:
    fprintf(stdoutMPI, "%s", cLogExpecEnergyStart);
    TimeKeeper(X, cFileNameTimeKeep, cExpecStart, "a");

    break;
  case TPQCalc:
#ifdef _DEBUG
    fprintf(stdoutMPI, "%s", cLogExpecEnergyStart);
    TimeKeeperWithStep(X, cFileNameTimeKeep, cTPQExpecStart, "a", step_i);
#endif
    break;
  case FullDiag:
    break;
  default:
    return -1;
    //break;
  }

  i_max=X->Check.idim_max;      
  if(GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit)!=0){
    return -1;
  }

  X->Large.i_max    = i_max;
  X->Large.irght    = irght;
  X->Large.ilft     = ilft;
  X->Large.ihfbit   = ihfbit;
  X->Large.mode     = M_ENERGY;
  X->Phys.energy=0.0;    
  dam_pr=0.0;
        
  // tentative doublon
  tmp_D        = 0.0;
  tmp_D2       = 0.0;
  tmp_N        = 0.0;
  tmp_N2       = 0.0;
  tmp_Sz       = 0.0;
  tmp_Sz2      = 0.0;
  tmp_num_up   = 0.0;
  tmp_num_down = 0.0;


  int nCalcFlct;
  if(X->Def.iCalcType == Lanczos){
    nCalcFlct=4301;
  }
  else if (X->Def.iCalcType == TPQCalc){
    nCalcFlct=3201;
  }
  else{//For FullDiag
    nCalcFlct=5301;
  }
  StartTimer(nCalcFlct);
  
  switch(X->Def.iCalcModel){
  case HubbardGC:
#pragma omp parallel for reduction(+:tmp_D,tmp_D2,tmp_N,tmp_N2,tmp_Sz,tmp_Sz2, tmp_num_up, tmp_num_down) default(none) shared(v0) \
  firstprivate(i_max, num1_up, num1_down,X,myrank) private(j, tmp_v02,D,N,Sz,isite1,is1_up,is1_down,is1,ibit1)
  for(j = 1; j <= i_max; j++){
    tmp_v02 = conj(v0[j])*v0[j];
    D       = 0.0;
    N       = 0.0;
    Sz      = 0.0;
    for(isite1=1;isite1<=X->Def.NsiteMPI;isite1++){
      if(isite1 > X->Def.Nsite){
        is1_up    = X->Def.Tpow[2 * isite1 - 2];
        is1_down  = X->Def.Tpow[2 * isite1 - 1];
        is1       = is1_up+is1_down;
        ibit1     = (unsigned long int)myrank & is1;
        num1_up   = (ibit1&is1_up) / is1_up;
        num1_down = (ibit1&is1_down) / is1_down;

        D        += num1_up*num1_down;
        N        += num1_up+num1_down;
        Sz       += num1_up-num1_down;
      }else{
        is1_up    = X->Def.Tpow[2*isite1-2];
        is1_down  = X->Def.Tpow[2*isite1-1];
        is1       = is1_up+is1_down;
        ibit1     = (j-1)&is1;
        num1_up   = ((j-1)&is1_up)/is1_up;
        num1_down = ((j-1)&is1_down)/is1_down;

        D        += num1_up*num1_down;
        N        += num1_up+num1_down;
        Sz       += num1_up-num1_down;
      }
    }
    tmp_D   += tmp_v02*D;
    tmp_D2  += tmp_v02*D*D;
    tmp_N   += tmp_v02*N;
    tmp_N2  += tmp_v02*N*N;
    tmp_Sz  += tmp_v02*Sz;
    tmp_Sz2 += tmp_v02*Sz*Sz;
  } 
  break;
  case KondoGC:
  case Hubbard:
  case Kondo:
#pragma omp parallel for reduction(+:tmp_D,tmp_D2,tmp_N,tmp_N2,tmp_Sz,tmp_Sz2) default(none) shared(v0,list_1) \
  firstprivate(i_max, num1_up, num1_down,X,myrank) private(j, tmp_v02,D,N,Sz,isite1,is1_up,is1_down,is1,ibit1,tmp_list_1)
  for(j = 1; j <= i_max; j++){
    tmp_v02 = conj(v0[j])*v0[j];
    D       = 0.0;
    N       = 0.0;
    Sz      = 0.0;
    tmp_list_1 = list_1[j];
    for(isite1=1;isite1<=X->Def.NsiteMPI;isite1++){
      //printf("DEBUG: j=%d %d %d\n",j,isite1,myrank);
      if(isite1 > X->Def.Nsite){
        is1_up    = X->Def.Tpow[2 * isite1 - 2];
        is1_down  = X->Def.Tpow[2 * isite1 - 1];
        is1       = is1_up+is1_down;
        ibit1     = (unsigned long int)myrank & is1;
        num1_up   = (ibit1&is1_up) / is1_up;
        num1_down = (ibit1&is1_down) / is1_down;

        D        += num1_up*num1_down;
        N        += num1_up+num1_down;
        Sz       += num1_up-num1_down;
      }else{
        is1_up    = X->Def.Tpow[2*isite1-2];
        is1_down  = X->Def.Tpow[2*isite1-1];
        is1       = is1_up+is1_down;
        //ibit1     = tmp_list_1&is1;
        num1_up   = (tmp_list_1&is1_up)/is1_up;
        num1_down = (tmp_list_1&is1_down)/is1_down;
        D        += num1_up*num1_down;
        N        += num1_up+num1_down;
        Sz       += num1_up-num1_down;
      } 
    }
    tmp_D   += tmp_v02*D;
    tmp_D2  += tmp_v02*D*D;
    tmp_N   += tmp_v02*N;
    tmp_N2  += tmp_v02*N*N;
    tmp_Sz  += tmp_v02*Sz;
    tmp_Sz2 += tmp_v02*Sz*Sz;
  }
  break;
  
  case SpinGC:
  if(X->Def.iFlgGeneralSpin == FALSE) {
#pragma omp parallel for reduction(+:tmp_Sz,tmp_Sz2)default(none) shared(v0)   \
  firstprivate(i_max,X,myrank) private(j,Sz, is1_up,ibit1,isite1,tmp_v02)
    for(j = 1; j <= i_max; j++){ 
      tmp_v02  = conj(v0[j])*v0[j];
      Sz       = 0.0;
      for(isite1=1;isite1<=X->Def.NsiteMPI;isite1++){
        if(isite1 > X->Def.Nsite){
          is1_up = X->Def.Tpow[isite1 - 1];
          ibit1  = (unsigned long int)myrank& is1_up;
          if(ibit1==is1_up){
            Sz += 1.0; 
          }else{
            Sz += -1.0; 
          }
        }else{
          is1_up=X->Def.Tpow[isite1-1];
          ibit1=(j-1)&is1_up;
          if(ibit1==is1_up){
            Sz += 1.0; 
          }else{
            Sz += -1.0; 
          }
        }
      }
      tmp_Sz   += Sz*tmp_v02;
      tmp_Sz2  += Sz*Sz*tmp_v02;
    }
  }
  else{//for generalspin
    for(j = 1; j <= i_max; j++){ 
      tmp_v02  = conj(v0[j])*v0[j];
      Sz       = 0.0;
      for(isite1=1;isite1<=X->Def.NsiteMPI;isite1++){
        //prefactor 0.5 is added later.
        if(isite1 > X->Def.Nsite){
          Sz += GetLocal2Sz(isite1, myrank, X->Def.SiteToBit, X->Def.Tpow);
        }else{
          Sz += GetLocal2Sz(isite1, j-1, X->Def.SiteToBit, X->Def.Tpow);
        }
      }
      tmp_Sz   += Sz*tmp_v02;
      tmp_Sz2  += Sz*Sz*tmp_v02;
    }
  }
  break;/*case SpinGC*/
  /* SpinGCBoost */
  case Spin:
    break;
  default:
    return -1;
  }

  tmp_D        = SumMPI_d(tmp_D);
  tmp_D2       = SumMPI_d(tmp_D2);
  tmp_N        = SumMPI_d(tmp_N);
  tmp_N2       = SumMPI_d(tmp_N2);
  tmp_Sz       = SumMPI_d(tmp_Sz);
  tmp_Sz2      = SumMPI_d(tmp_Sz2);
//  tmp_num_up   = SumMPI_d(tmp_num_up);
//  tmp_num_down = SumMPI_d(tmp_num_down);
  
  switch(X->Def.iCalcModel){
  case HubbardGC:
  case KondoGC:
  case Hubbard:
  case Kondo:
    X->Phys.doublon   = tmp_D;
    X->Phys.doublon2  = tmp_D2;
    X->Phys.num       = tmp_N;
    X->Phys.num2      = tmp_N2;
    X->Phys.Sz        = tmp_Sz*0.5;
    X->Phys.Sz2       = tmp_Sz2*0.25;
    X->Phys.num_up    = 0.5*(tmp_N+tmp_Sz);
    X->Phys.num_down  = 0.5*(tmp_N-tmp_Sz);
    break;      
      
  case SpinGC:
    X->Phys.doublon   = 0.0;
    X->Phys.doublon2  = 0.0;
    X->Phys.num       = X->Def.NsiteMPI;
    X->Phys.num2      = X->Def.NsiteMPI*X->Def.NsiteMPI;
    X->Phys.Sz        = tmp_Sz*0.5;
    X->Phys.Sz2       = tmp_Sz2*0.25;
    X->Phys.num_up    = 0.5*(X->Def.NsiteMPI+tmp_Sz);
    X->Phys.num_down  = 0.5*(X->Def.NsiteMPI-tmp_Sz);
    break;
  case Spin:
    X->Phys.doublon   = 0.0;
    X->Phys.doublon2  = 0.0;
    X->Phys.num_up    = X->Def.Nup;
    X->Phys.num_down  = X->Def.Ndown;    
    X->Phys.num       = (X->Def.Nup+X->Def.Ndown);
    X->Phys.num2      = (X->Def.Nup+X->Def.Ndown)*(X->Def.Nup+X->Def.Ndown);
    X->Phys.Sz        = 0.5*(X->Def.Total2SzMPI);
    X->Phys.Sz2       = 0.25*pow((X->Def.Total2SzMPI),2);
    break;
  default:
    return -1;
  }
  StopTimer(nCalcFlct);       

#pragma omp parallel for default(none) private(i) shared(v1,v0) firstprivate(i_max)
  for(i = 1; i <= i_max; i++){
    v1[i]=v0[i];
    v0[i]=0.0+0.0*I;
  }


  int nCalcExpec;
  if(X->Def.iCalcType == Lanczos){
    nCalcExpec=4302;
  }
  else if (X->Def.iCalcType == TPQCalc){
    nCalcExpec=3202;
  }
  else{//For FullDiag
    nCalcExpec=5302;
  }
  StartTimer(nCalcExpec);
  mltply(X, v0, v1); // v0+=H*v1
  StopTimer(nCalcExpec);
/* switch -> SpinGCBoost */

  dam_pr=0.0;
  dam_pr1=0.0;
  #pragma omp parallel for default(none) reduction(+:dam_pr, dam_pr1) private(j) shared(v0, v1)firstprivate(i_max) 
  for(j=1;j<=i_max;j++){
    dam_pr   += conj(v1[j])*v0[j]; // E   = <v1|H|v1>=<v1|v0>
    dam_pr1  += conj(v0[j])*v0[j]; // E^2 = <v1|H*H|v1>=<v0|v0>
    //v0[j]=v1[j]; v1-> orginal v0=H*v1
  }  
  dam_pr = SumMPI_dc(dam_pr);
  dam_pr1 = SumMPI_dc(dam_pr1);
  //  fprintf(stdoutMPI, "Debug: ene=%lf, var=%lf\n", creal(dam_pr), creal(dam_pr1));
  
  X->Phys.energy = dam_pr;
  X->Phys.var    = dam_pr1;

  switch(X->Def.iCalcType){
  case Lanczos:
    fprintf(stdoutMPI, "%s", cLogExpecEnergyEnd);
    TimeKeeper(X, cFileNameTimeKeep, cExpecEnd, "a");
    break;
  case TPQCalc:
#ifdef _DEBUG
      fprintf(stdoutMPI, "%s", cLogExpecEnergyEnd);
      TimeKeeperWithStep(X, cFileNameTimeKeep, cTPQExpecEnd, "a", step_i);
#endif
    break;
      default:
    break;
  }
  

  return 0;
}
