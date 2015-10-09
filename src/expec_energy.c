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
#include "expec_energy.h"

/** 
 * 
 * 
 * @param X 
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return 
 */
int expec_energy(struct BindStruct *X){

  char sdt[D_FileNameMax];

  long unsigned int i,j;
  long unsigned int irght,ilft,ihfbit;
  long unsigned int isite1;
  long unsigned int is1_up,is1_down;
  long unsigned int is1;
  double complex dam_pr,dam_pr1;
    
  long unsigned int num1_up, num1_down;
  long unsigned int ibit1;
  double tmp_num_up, tmp_num_down, tmp_doublon;
  double tmp_v02;  
  
  long unsigned int i_max;

  switch(X->Def.iCalcType){
  case Lanczos:
    printf("%s", cLogExpecEnergyStart);
    break;
  case TPQCalc:
    TimeKeeperWithStep(X, cFileNameTimeKeep, cExpecStart, "a", step_i);
    break;
  case FullDiag:
    break;
  default:
    return -1;
    break;
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
  switch(X->Def.iCalcModel){
  case HubbardGC:
    tmp_doublon=0.0;
    tmp_num_up=0.0;
    tmp_num_down=0.0;
    for(isite1=1;isite1<=X->Def.Nsite;isite1++){
#pragma omp parallel for reduction(+:tmp_doublon, tmp_num_up, tmp_num_down) default(none) private(j, is1_up, is1_down, is1, ibit1,num1_up,num1_down, tmp_v02) shared(list_1, v0) firstprivate(i_max, X, isite1) 
      for(j=1;j<=i_max;j++){  
		is1_up=X->Def.Tpow[2*isite1-2];
		is1_down=X->Def.Tpow[2*isite1-1];
		is1=is1_up+is1_down;
		ibit1=(j-1)&is1;
		tmp_v02 = conj(v0[j])*v0[j];
		if(ibit1==is1){
		  tmp_doublon  += tmp_v02;
		}
		num1_up   = ((j-1)&is1_up)/is1_up;
		num1_down = ((j-1)&is1_down)/is1_down;

		tmp_num_up   += tmp_v02*num1_up;
		tmp_num_down += tmp_v02*num1_down;
      }
    }
    X->Phys.doublon  = tmp_doublon;
    X->Phys.num_up   = tmp_num_up;
    X->Phys.num_down = tmp_num_down;    
    X->Phys.num      = tmp_num_up+tmp_num_down;
    break;
  case KondoGC:
  case Hubbard:
  case Kondo:
    X->Phys.num=X->Def.Nup+X->Def.Ndown;//canonical
    tmp_doublon=0.0;
    tmp_num_up=0;
    tmp_num_down=0;
    for(isite1=1;isite1<=X->Def.Nsite;isite1++){
#pragma omp parallel for reduction(+:tmp_doublon, tmp_num_up, tmp_num_down) default(none) private(j, is1_up, is1_down, is1, ibit1, tmp_v02, num1_up, num1_down) shared(list_1,v0) firstprivate(i_max, X, isite1) 
      for(j=1;j<=i_max;j++){  
		is1_up=X->Def.Tpow[2*isite1-2];
		is1_down=X->Def.Tpow[2*isite1-1];
		is1=is1_up+is1_down;
		ibit1=list_1[j]&is1;
		tmp_v02 = conj(v0[j])*v0[j];
		if(ibit1==is1){
		  tmp_doublon  += tmp_v02;
		}
		num1_up   = (list_1[j]&is1_up)/is1_up;
		num1_down = (list_1[j]&is1_down)/is1_down;
		tmp_num_up += tmp_v02*num1_up;
		tmp_num_down += tmp_v02*num1_down;
      }
    }
    X->Phys.num_up   = tmp_num_up;
    X->Phys.num_down = tmp_num_down;    
    X->Phys.num      = tmp_num_up+tmp_num_down;
    X->Phys.doublon  = tmp_doublon;
    break;

  case SpinGC:
    tmp_num_up=0.0;
    tmp_num_down=0.0;
    for(isite1=1;isite1<=X->Def.Nsite;isite1++){
      is1_up=X->Def.Tpow[isite1-1];
#pragma omp parallel for reduction(+: tmp_num_up, tmp_num_down) default(none) private(j,  ibit1,num1_up,num1_down, tmp_v02) shared(list_1, v0) firstprivate(i_max, X, isite1, is1_up) 
      for(j=1;j<=i_max;j++){
	ibit1=(j-1)&is1_up;	
	tmp_v02 = conj(v0[j])*v0[j];
	if(ibit1==is1_up){
	  tmp_num_up  += tmp_v02;
	}
	else{
	  tmp_num_down +=tmp_v02;
	}
      }
    }
    X->Phys.doublon  = 0.0;
    X->Phys.num_up   = tmp_num_up;
    X->Phys.num_down = tmp_num_down;    
    X->Phys.num      = tmp_num_up+tmp_num_down;    
    break;
    
  case Spin:
    X->Phys.num_up   = X->Def.Nup;
    X->Phys.num_down = X->Def.Ndown;    
    X->Phys.num      = X->Def.Nup+X->Def.Ndown;//canonical
    X->Phys.doublon  = 0.0;// spin
    break;
    
  default:
    return -1;
  }
  
#pragma omp parallel for default(none) private(i) shared(v1,v0) firstprivate(i_max)
  for(i = 1; i <= i_max; i++){
    v1[i]=v0[i];
    v0[i]=0.0+0.0*I;
  }

  mltply(X, v0, v1); // v0+=H*v1

  dam_pr=0.0;
  dam_pr1=0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr, dam_pr1) private(j) shared(v0, v1)firstprivate(i_max) 
  for(j=1;j<=i_max;j++){
    dam_pr   += conj(v1[j])*v0[j]; // E   = <v1|H|v1>=<v1|v0>
    dam_pr1  += conj(v0[j])*v0[j]; // E^2 = <v1|H*H|v1>=<v0|v0>
    //v0[j]=v1[j]; v1-> orginal v0=H*v1
  }  
  X->Phys.energy = dam_pr;
  X->Phys.var    = dam_pr1;

  switch(X->Def.iCalcType){
  case Lanczos:
    if(X->Def.St==0){
      sprintf(sdt, cFileNameEnergy_Lanczos, X->Def.CDataFileHead);
    }else if(X->Def.St==1){
      sprintf(sdt, cFileNameEnergy_CG, X->Def.CDataFileHead);
    }
	
    FILE *fp;    
    if(childfopen(sdt, "w", &fp)!=0){
      return -1;
    }  
    fprintf(fp,"Energy  %.10lf \n",X->Phys.energy);
    fprintf(fp,"Doublon  %.10lf \n",X->Phys.doublon);
    fclose(fp);
	
    break;    
  case TPQCalc:
    TimeKeeperWithStep(X, cFileNameTimeKeep, cExpecEnd, "a", step_i);
    break;
  case FullDiag:
    break;
  default:
    return -1;
    break;
  }
  
  if(X->Def.iCalcType==Lanczos){
    if(X->Def.St==0){
      TimeKeeper(X, cFileNameTimeKeep, cLanczosExpecEnergyEnd, "a");
      printf("%s", cLogLanczosExpecEnergyEnd);
    }else if(X->Def.St==1){
      TimeKeeper(X, cFileNameTimeKeep, cCGExpecEnergyEnd, "a");
      printf("%s", cLogCGExpecEnergyEnd);
    }
    printf(cLogEnergy, X->Phys.energy);
  }
  
  return 0;
}
