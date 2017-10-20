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
#include "mltplyCommon.h"
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
  double complex dam_pr,dam_pr1;
  long unsigned int i_max;

#if defined(_ORG)
  long unsigned int ibit1, ibit2;
  long unsigned int ibit_up,ibit_down,ibit_D;
  double D,tmp_D,tmp_D2;
  double N,tmp_N,tmp_N2;
  double Sz,Sz2,tmp_Sz,tmp_Sz2;
  double tmp_v02;
  long unsigned int tmp_list_1;
  unsigned int l_ibit1,u_ibit1,i_32;

  i_32   = (unsigned int)(pow(2,32)-1);
#endif

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
  case CG:
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
#ifdef _ORG
  tmp_D        = 0.0;
  tmp_D2       = 0.0;
  tmp_N        = 0.0;
  tmp_N2       = 0.0;
  tmp_Sz       = 0.0;
  tmp_Sz2      = 0.0;
#endif

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
      expec_energy_flct_HubbardGC(X);
#ifdef _ORG
//[s] for bit count
  is1_up_a   = 0;
  is1_up_b   = 0;
  is1_down_a = 0;
  is1_down_b = 0;
  for(isite1=1;isite1<=X->Def.NsiteMPI;isite1++){
    if(isite1 > X->Def.Nsite){
      is1_up_a   += X->Def.Tpow[2*isite1 - 2];
      is1_down_a += X->Def.Tpow[2*isite1 - 1];
    }else{
      is1_up_b   += X->Def.Tpow[2*isite1 - 2];
      is1_down_b += X->Def.Tpow[2*isite1 - 1];
    }
  }
//[e]
#pragma omp parallel for reduction(+:tmp_D,tmp_D2,tmp_N,tmp_N2,tmp_Sz,tmp_Sz2) default(none) shared(v0,list_1) \
  firstprivate(i_max, X,myrank,is1_up_a,is1_down_a,is1_up_b,is1_down_b,i_32) \
  private(j, tmp_v02,D,N,Sz,isite1,ibit1,tmp_list_1,bit_up,bit_down,bit_D,u_ibit1,l_ibit1,ibit_up,ibit_down,ibit_D)
  for(j = 1; j <= i_max; j++){
    tmp_v02 = conj(v0[j])*v0[j];
    bit_up     = 0;
    bit_down   = 0;
    bit_D      = 0;
// isite1 > X->Def.Nsite
    ibit_up    = (unsigned long int) myrank & is1_up_a;
    u_ibit1    = ibit_up >> 32;
    l_ibit1    = ibit_up & i_32;
    bit_up    += pop(u_ibit1);
    bit_up    += pop(l_ibit1);

    ibit_down  = (unsigned long int) myrank & is1_down_a;
    u_ibit1    = ibit_down >> 32;
    l_ibit1    = ibit_down & i_32;
    bit_down  += pop(u_ibit1);
    bit_down  += pop(l_ibit1);

    ibit_D     = (ibit_up) & (ibit_down >> 1);
    u_ibit1    = ibit_D >> 32;
    l_ibit1    = ibit_D & i_32;
    bit_D     += pop(u_ibit1);
    bit_D     += pop(l_ibit1);

// isite1 <= X->Def.Nsite
    ibit_up    = (unsigned long int) (j-1) & is1_up_b;
    u_ibit1    = ibit_up >> 32;
    l_ibit1    = ibit_up & i_32;
    bit_up    += pop(u_ibit1);
    bit_up    += pop(l_ibit1);

    ibit_down  = (unsigned long int) (j-1) & is1_down_b;
    u_ibit1    = ibit_down >> 32;
    l_ibit1    = ibit_down & i_32;
    bit_down  += pop(u_ibit1);
    bit_down  += pop(l_ibit1);

    ibit_D     = (ibit_up) & (ibit_down >> 1);
    u_ibit1    = ibit_D >> 32;
    l_ibit1    = ibit_D & i_32;
    bit_D     += pop(u_ibit1);
    bit_D     += pop(l_ibit1);

    D          =  bit_D;
    N          =  bit_up+bit_down;
    Sz         =  bit_up-bit_down;


/*
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
        ibit1     =  (unsigned long int)(j-1)&is1;
        num1_up   = ((j-1)&is1_up)/is1_up;
        num1_down = ((j-1)&is1_down)/is1_down;

        D        += num1_up*num1_down;
        N        += num1_up+num1_down;
        Sz       += num1_up-num1_down;
      }
    }
*/
    tmp_D   += tmp_v02*D;
    tmp_D2  += tmp_v02*D*D;
    tmp_N   += tmp_v02*N;
    tmp_N2  += tmp_v02*N*N;
    tmp_Sz  += tmp_v02*Sz;
    tmp_Sz2 += tmp_v02*Sz*Sz;
  }
#endif
  break;
  case KondoGC:
  case Hubbard:
  case Kondo:
      expec_energy_flct_Hubbard(X);
#ifdef _ORG
//[s] for bit count
  is1_up_a   = 0;
  is1_up_b   = 0;
  is1_down_a = 0;
  is1_down_b = 0;
  for(isite1=1;isite1<=X->Def.NsiteMPI;isite1++){
    if(isite1 > X->Def.Nsite){
      is1_up_a   += X->Def.Tpow[2*isite1 - 2];
      is1_down_a += X->Def.Tpow[2*isite1 - 1];
    }else{
      is1_up_b   += X->Def.Tpow[2*isite1 - 2];
      is1_down_b += X->Def.Tpow[2*isite1 - 1];
    }
  }
//[e]
#pragma omp parallel for reduction(+:tmp_D,tmp_D2,tmp_N,tmp_N2,tmp_Sz,tmp_Sz2) default(none) shared(v0,list_1) \
  firstprivate(i_max, X,myrank,is1_up_a,is1_down_a,is1_up_b,is1_down_b,i_32) \
  private(j, tmp_v02,D,N,Sz,isite1,ibit1,tmp_list_1,bit_up,bit_down,bit_D,u_ibit1,l_ibit1,ibit_up,ibit_down,ibit_D)
  for(j = 1; j <= i_max; j++){
    tmp_v02 = conj(v0[j])*v0[j];
    bit_up     = 0;
    bit_down   = 0;
    bit_D      = 0;
    tmp_list_1 = list_1[j];
// isite1 > X->Def.Nsite
    ibit_up    = (unsigned long int) myrank & is1_up_a;
    u_ibit1    = ibit_up >> 32;
    l_ibit1    = ibit_up & i_32;
    bit_up    += pop(u_ibit1);
    bit_up    += pop(l_ibit1);

    ibit_down  = (unsigned long int) myrank & is1_down_a;
    u_ibit1    = ibit_down >> 32;
    l_ibit1    = ibit_down & i_32;
    bit_down  += pop(u_ibit1);
    bit_down  += pop(l_ibit1);

    ibit_D     = (ibit_up) & (ibit_down >> 1);
    u_ibit1    = ibit_D >> 32;
    l_ibit1    = ibit_D & i_32;
    bit_D     += pop(u_ibit1);
    bit_D     += pop(l_ibit1);

// isite1 <= X->Def.Nsite
    ibit_up    = (unsigned long int) tmp_list_1 & is1_up_b;
    u_ibit1    = ibit_up >> 32;
    l_ibit1    = ibit_up & i_32;
    bit_up    += pop(u_ibit1);
    bit_up    += pop(l_ibit1);

    ibit_down  = (unsigned long int) tmp_list_1 & is1_down_b;
    u_ibit1    = ibit_down >> 32;
    l_ibit1    = ibit_down & i_32;
    bit_down  += pop(u_ibit1);
    bit_down  += pop(l_ibit1);

    ibit_D     = (ibit_up) & (ibit_down >> 1);
    u_ibit1    = ibit_D >> 32;
    l_ibit1    = ibit_D & i_32;
    bit_D     += pop(u_ibit1);
    bit_D     += pop(l_ibit1);

    D          =  bit_D;
    N          =  bit_up+bit_down;
    Sz         =  bit_up-bit_down;

/*
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
        //ibit1     = tmp_list_1&is1;
        num1_up   = (tmp_list_1&is1_up)/is1_up;
        num1_down = (tmp_list_1&is1_down)/is1_down;
        D        += num1_up*num1_down;
        N        += num1_up+num1_down;
        Sz       += num1_up-num1_down;
      } 
    }
*/
    tmp_D   += tmp_v02*D;
    tmp_D2  += tmp_v02*D*D;
    tmp_N   += tmp_v02*N;
    tmp_N2  += tmp_v02*N*N;
    tmp_Sz  += tmp_v02*Sz;
    tmp_Sz2 += tmp_v02*Sz*Sz;
  }
#endif
  break;
  
  case SpinGC:
  if(X->Def.iFlgGeneralSpin == FALSE) {
      expec_energy_flct_HalfSpinGC(X);
#ifdef _ORG
//[s] for bit count
    is1_up_a = 0;
    is1_up_b = 0;
    for(isite1=1;isite1<=X->Def.NsiteMPI;isite1++){
      if(isite1 > X->Def.Nsite){
        is1_up_a += X->Def.Tpow[isite1 - 1];
      }else{
        is1_up_b += X->Def.Tpow[isite1 - 1];
      }
    }
//[e]
#pragma omp parallel for reduction(+:tmp_Sz,tmp_Sz2)default(none) shared(v0)   \
  firstprivate(i_max,X,myrank,i_32,is1_up_a,is1_up_b) private(j,Sz,ibit1,isite1,tmp_v02,u_ibit1,l_ibit1)
    for(j = 1; j <= i_max; j++){ 
      tmp_v02  = conj(v0[j])*v0[j];
      Sz       = 0.0;

// isite1 > X->Def.Nsite
      ibit1   = (unsigned long int) myrank & is1_up_a;
      u_ibit1 = ibit1 >> 32;
      l_ibit1 = ibit1 & i_32;
      Sz      += pop(u_ibit1);
      Sz      += pop(l_ibit1);
// isite1 <= X->Def.Nsite
      ibit1   = (unsigned long int) (j-1)&is1_up_b;
      u_ibit1 = ibit1 >> 32;
      l_ibit1 = ibit1 & i_32;
      Sz     += pop(u_ibit1);
      Sz     += pop(l_ibit1);
      Sz      = 2*Sz-X->Def.NsiteMPI;

/*
      for(isite1=1;isite1<=X->Def.NsiteMPI;isite1++){
        if(isite1 > X->Def.Nsite){
          is1_up = X->Def.Tpow[isite1 - 1];
          ibit1  = (unsigned long int)myrank& is1_up;
          if(ibit1==is1_up){
            Sz = 1.0; 
          }else{
            Sz = -1.0; 
          }
        }else{
          is1_up=X->Def.Tpow[isite1-1];
          ibit1=(j-1)&is1_up;
          if(ibit1==is1_up){
            Sz = 1.0; 
          }else{
            Sz = -1.0; 
          }
        }
     }

*/
      tmp_Sz   += Sz*tmp_v02;
      tmp_Sz2  += Sz*Sz*tmp_v02;
    }
#endif
  }
  else{//for generalspin
      expec_energy_flct_GeneralSpinGC(X);
#ifdef _ORG
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
#endif
  }
  break;/*case SpinGC*/
  /* SpinGCBoost */
  case Spin:
    break;
  default:
    return -1;
  }

#ifdef _ORG
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
#endif

  StopTimer(nCalcFlct);       

#pragma omp parallel for default(none) private(i) shared(v1,v0) firstprivate(i_max)
  for(i = 1; i <= i_max; i++){
    v1[i]=v0[i];
    v0[i]=0.0;
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

int expec_energy_flct_HubbardGC(struct BindStruct *X) {
    long unsigned int j;
    long unsigned int isite1;
    long unsigned int is1_up_a, is1_up_b;
    long unsigned int is1_down_a, is1_down_b;
    int bit_up, bit_down, bit_D;
    long unsigned int ibit_up, ibit_down, ibit_D;
    double D, tmp_D, tmp_D2;
    double N, tmp_N, tmp_N2;
    double Sz, tmp_Sz, tmp_Sz2;
    double tmp_v02;
    long unsigned int i_max;
    unsigned int l_ibit1, u_ibit1, i_32;
    i_max=X->Check.idim_max;

    i_32 = 0xFFFFFFFF; //2^32 - 1
    // tentative doublon
    tmp_D        = 0.0;
    tmp_D2       = 0.0;
    tmp_N        = 0.0;
    tmp_N2       = 0.0;
    tmp_Sz       = 0.0;
    tmp_Sz2      = 0.0;

//[s] for bit count
    is1_up_a = 0;
    is1_up_b = 0;
    is1_down_a = 0;
    is1_down_b = 0;
    for (isite1 = 1; isite1 <= X->Def.NsiteMPI; isite1++) {
        if (isite1 > X->Def.Nsite) {
            is1_up_a += X->Def.Tpow[2 * isite1 - 2];
            is1_down_a += X->Def.Tpow[2 * isite1 - 1];
        } else {
            is1_up_b += X->Def.Tpow[2 * isite1 - 2];
            is1_down_b += X->Def.Tpow[2 * isite1 - 1];
        }
    }
//[e]
#pragma omp parallel for reduction(+:tmp_D,tmp_D2,tmp_N,tmp_N2,tmp_Sz,tmp_Sz2) default(none) shared(v0,list_1) \
  firstprivate(i_max, X,myrank,is1_up_a,is1_down_a,is1_up_b,is1_down_b,i_32) \
  private(j, tmp_v02,D,N,Sz,isite1,bit_up,bit_down,bit_D,u_ibit1,l_ibit1,ibit_up,ibit_down,ibit_D)
    for (j = 1; j <= i_max; j++) {
        tmp_v02 = conj(v0[j]) * v0[j];
        bit_up = 0;
        bit_down = 0;
        bit_D = 0;
// isite1 > X->Def.Nsite
        ibit_up = (unsigned long int) myrank & is1_up_a;
        u_ibit1 = ibit_up >> 32;
        l_ibit1 = ibit_up & i_32;
        bit_up += pop(u_ibit1);
        bit_up += pop(l_ibit1);

        ibit_down = (unsigned long int) myrank & is1_down_a;
        u_ibit1 = ibit_down >> 32;
        l_ibit1 = ibit_down & i_32;
        bit_down += pop(u_ibit1);
        bit_down += pop(l_ibit1);

        ibit_D = (ibit_up) & (ibit_down >> 1);
        u_ibit1 = ibit_D >> 32;
        l_ibit1 = ibit_D & i_32;
        bit_D += pop(u_ibit1);
        bit_D += pop(l_ibit1);

// isite1 <= X->Def.Nsite
        ibit_up = (unsigned long int) (j - 1) & is1_up_b;
        u_ibit1 = ibit_up >> 32;
        l_ibit1 = ibit_up & i_32;
        bit_up += pop(u_ibit1);
        bit_up += pop(l_ibit1);

        ibit_down = (unsigned long int) (j - 1) & is1_down_b;
        u_ibit1 = ibit_down >> 32;
        l_ibit1 = ibit_down & i_32;
        bit_down += pop(u_ibit1);
        bit_down += pop(l_ibit1);

        ibit_D = (ibit_up) & (ibit_down >> 1);
        u_ibit1 = ibit_D >> 32;
        l_ibit1 = ibit_D & i_32;
        bit_D += pop(u_ibit1);
        bit_D += pop(l_ibit1);

        D = bit_D;
        N = bit_up + bit_down;
        Sz = bit_up - bit_down;

        tmp_D += tmp_v02 * D;
        tmp_D2 += tmp_v02 * D * D;
        tmp_N += tmp_v02 * N;
        tmp_N2 += tmp_v02 * N * N;
        tmp_Sz += tmp_v02 * Sz;
        tmp_Sz2 += tmp_v02 * Sz * Sz;
    }
    tmp_D        = SumMPI_d(tmp_D);
    tmp_D2       = SumMPI_d(tmp_D2);
    tmp_N        = SumMPI_d(tmp_N);
    tmp_N2       = SumMPI_d(tmp_N2);
    tmp_Sz       = SumMPI_d(tmp_Sz);
    tmp_Sz2      = SumMPI_d(tmp_Sz2);

    X->Phys.doublon   = tmp_D;
    X->Phys.doublon2  = tmp_D2;
    X->Phys.num       = tmp_N;
    X->Phys.num2      = tmp_N2;
    X->Phys.Sz        = tmp_Sz*0.5;
    X->Phys.Sz2       = tmp_Sz2*0.25;
    X->Phys.num_up    = 0.5*(tmp_N+tmp_Sz);
    X->Phys.num_down  = 0.5*(tmp_N-tmp_Sz);

    return 0;
}

int expec_energy_flct_Hubbard(struct BindStruct *X){
    long unsigned int j;
    long unsigned int isite1;
    long unsigned int is1_up_a,is1_up_b;
    long unsigned int is1_down_a,is1_down_b;
    int bit_up,bit_down,bit_D;

    long unsigned int ibit_up,ibit_down,ibit_D;
    double D,tmp_D,tmp_D2;
    double N,tmp_N,tmp_N2;
    double Sz,tmp_Sz, tmp_Sz2;
    double tmp_v02;
    long unsigned int i_max,tmp_list_1;
    unsigned int l_ibit1,u_ibit1,i_32;
    i_max=X->Check.idim_max;

    i_32   = (unsigned int)(pow(2,32)-1);

    tmp_D        = 0.0;
    tmp_D2       = 0.0;
    tmp_N        = 0.0;
    tmp_N2       = 0.0;
    tmp_Sz       = 0.0;
    tmp_Sz2      = 0.0;

    //[s] for bit count
    is1_up_a   = 0;
    is1_up_b   = 0;
    is1_down_a = 0;
    is1_down_b = 0;
    for(isite1=1;isite1<=X->Def.NsiteMPI;isite1++){
        if(isite1 > X->Def.Nsite){
            is1_up_a   += X->Def.Tpow[2*isite1 - 2];
            is1_down_a += X->Def.Tpow[2*isite1 - 1];
        }else{
            is1_up_b   += X->Def.Tpow[2*isite1 - 2];
            is1_down_b += X->Def.Tpow[2*isite1 - 1];
        }
    }
//[e]
#pragma omp parallel for reduction(+:tmp_D,tmp_D2,tmp_N,tmp_N2,tmp_Sz,tmp_Sz2) default(none) shared(v0,list_1) \
  firstprivate(i_max, X,myrank,is1_up_a,is1_down_a,is1_up_b,is1_down_b,i_32) \
  private(j, tmp_v02,D,N,Sz,isite1,tmp_list_1,bit_up,bit_down,bit_D,u_ibit1,l_ibit1,ibit_up,ibit_down,ibit_D)
    for(j = 1; j <= i_max; j++) {
        tmp_v02 = conj(v0[j]) * v0[j];
        bit_up = 0;
        bit_down = 0;
        bit_D = 0;
        tmp_list_1 = list_1[j];
// isite1 > X->Def.Nsite
        ibit_up = (unsigned long int) myrank & is1_up_a;
        u_ibit1 = ibit_up >> 32;
        l_ibit1 = ibit_up & i_32;
        bit_up += pop(u_ibit1);
        bit_up += pop(l_ibit1);

        ibit_down = (unsigned long int) myrank & is1_down_a;
        u_ibit1 = ibit_down >> 32;
        l_ibit1 = ibit_down & i_32;
        bit_down += pop(u_ibit1);
        bit_down += pop(l_ibit1);

        ibit_D = (ibit_up) & (ibit_down >> 1);
        u_ibit1 = ibit_D >> 32;
        l_ibit1 = ibit_D & i_32;
        bit_D += pop(u_ibit1);
        bit_D += pop(l_ibit1);

// isite1 <= X->Def.Nsite
        ibit_up = (unsigned long int) tmp_list_1 & is1_up_b;
        u_ibit1 = ibit_up >> 32;
        l_ibit1 = ibit_up & i_32;
        bit_up += pop(u_ibit1);
        bit_up += pop(l_ibit1);

        ibit_down = (unsigned long int) tmp_list_1 & is1_down_b;
        u_ibit1 = ibit_down >> 32;
        l_ibit1 = ibit_down & i_32;
        bit_down += pop(u_ibit1);
        bit_down += pop(l_ibit1);

        ibit_D = (ibit_up) & (ibit_down >> 1);
        u_ibit1 = ibit_D >> 32;
        l_ibit1 = ibit_D & i_32;
        bit_D += pop(u_ibit1);
        bit_D += pop(l_ibit1);

        D = bit_D;
        N = bit_up + bit_down;
        Sz = bit_up - bit_down;

        tmp_D += tmp_v02 * D;
        tmp_D2 += tmp_v02 * D * D;
        tmp_N += tmp_v02 * N;
        tmp_N2 += tmp_v02 * N * N;
        tmp_Sz += tmp_v02 * Sz;
        tmp_Sz2 += tmp_v02 * Sz * Sz;
    }


    tmp_D        = SumMPI_d(tmp_D);
    tmp_D2       = SumMPI_d(tmp_D2);
    tmp_N        = SumMPI_d(tmp_N);
    tmp_N2       = SumMPI_d(tmp_N2);
    tmp_Sz       = SumMPI_d(tmp_Sz);
    tmp_Sz2      = SumMPI_d(tmp_Sz2);

    X->Phys.doublon   = tmp_D;
    X->Phys.doublon2  = tmp_D2;
    X->Phys.num       = tmp_N;
    X->Phys.num2      = tmp_N2;
    X->Phys.Sz        = tmp_Sz*0.5;
    X->Phys.Sz2       = tmp_Sz2*0.25;
    X->Phys.num_up    = 0.5*(tmp_N+tmp_Sz);
    X->Phys.num_down  = 0.5*(tmp_N-tmp_Sz);
    return 0;
}

int expec_energy_flct_HalfSpinGC(struct BindStruct *X){
    long unsigned int j;
    long unsigned int isite1;
    long unsigned int is1_up_a,is1_up_b;

    long unsigned int ibit1;
    double Sz,tmp_Sz, tmp_Sz2;
    double tmp_v02;
    long unsigned int i_max;
    unsigned int l_ibit1,u_ibit1,i_32;
    i_max=X->Check.idim_max;

    i_32 = 0xFFFFFFFF; //2^32 - 1

    // tentative doublon
    tmp_Sz       = 0.0;
    tmp_Sz2      = 0.0;

//[s] for bit count
    is1_up_a = 0;
    is1_up_b = 0;
    for(isite1=1;isite1<=X->Def.NsiteMPI;isite1++){
        if(isite1 > X->Def.Nsite){
            is1_up_a += X->Def.Tpow[isite1 - 1];
        }else{
            is1_up_b += X->Def.Tpow[isite1 - 1];
        }
    }
//[e]
#pragma omp parallel for reduction(+:tmp_Sz,tmp_Sz2)default(none) shared(v0)   \
  firstprivate(i_max,X,myrank,i_32,is1_up_a,is1_up_b) private(j,Sz,ibit1,isite1,tmp_v02,u_ibit1,l_ibit1)
    for(j = 1; j <= i_max; j++){
        tmp_v02  = conj(v0[j])*v0[j];
        Sz       = 0.0;

// isite1 > X->Def.Nsite
        ibit1   = (unsigned long int) myrank & is1_up_a;
        u_ibit1 = ibit1 >> 32;
        l_ibit1 = ibit1 & i_32;
        Sz      += pop(u_ibit1);
        Sz      += pop(l_ibit1);
// isite1 <= X->Def.Nsite
        ibit1   = (unsigned long int) (j-1)&is1_up_b;
        u_ibit1 = ibit1 >> 32;
        l_ibit1 = ibit1 & i_32;
        Sz     += pop(u_ibit1);
        Sz     += pop(l_ibit1);
        Sz      = 2*Sz-X->Def.NsiteMPI;

        tmp_Sz   += Sz*tmp_v02;
        tmp_Sz2  += Sz*Sz*tmp_v02;
    }
    tmp_Sz       = SumMPI_d(tmp_Sz);
    tmp_Sz2      = SumMPI_d(tmp_Sz2);

    X->Phys.doublon   = 0.0;
    X->Phys.doublon2  = 0.0;
    X->Phys.num       = X->Def.NsiteMPI;
    X->Phys.num2      = X->Def.NsiteMPI*X->Def.NsiteMPI;
    X->Phys.Sz        = tmp_Sz*0.5;
    X->Phys.Sz2       = tmp_Sz2*0.25;
    X->Phys.num_up    = 0.5*(X->Def.NsiteMPI+tmp_Sz);
    X->Phys.num_down  = 0.5*(X->Def.NsiteMPI-tmp_Sz);

    return 0;
}

int expec_energy_flct_GeneralSpinGC(struct BindStruct *X){
    long unsigned int j;
    long unsigned int isite1;

    double Sz,tmp_Sz, tmp_Sz2;
    double tmp_v02;
    long unsigned int i_max;
    i_max=X->Check.idim_max;

    // tentative doublon
    tmp_Sz       = 0.0;
    tmp_Sz2      = 0.0;

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

    tmp_Sz       = SumMPI_d(tmp_Sz);
    tmp_Sz2      = SumMPI_d(tmp_Sz2);

    X->Phys.doublon   = 0.0;
    X->Phys.doublon2  = 0.0;
    X->Phys.num       = X->Def.NsiteMPI;
    X->Phys.num2      = X->Def.NsiteMPI*X->Def.NsiteMPI;
    X->Phys.Sz        = tmp_Sz*0.5;
    X->Phys.Sz2       = tmp_Sz2*0.25;
    X->Phys.num_up    = 0.5*(X->Def.NsiteMPI+tmp_Sz);
    X->Phys.num_down  = 0.5*(X->Def.NsiteMPI-tmp_Sz);

    return 0;
}
