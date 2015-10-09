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
#include "expec_totalspin.h"

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
int expec_totalspin
(
 struct BindStruct *X,
 double complex *vec
 )
{
  switch(X->Def.iCalcModel){
  case Spin:
     totalspin_Spin(X,vec);
     break;
  case SpinGC:
     totalspin_SpinGC(X,vec);
     break;
   case Hubbard: 
   case Kondo:
     totalspin_Hubbard(X,vec);
     break;
  case HubbardGC: 
  case KondoGC:
     totalspin_HubbardGC(X,vec);
     break;
  //default:
     //X->Phys.s2=0.0;   
  }
  return 0;
}

void totalspin_Hubbard(struct BindStruct *X,double complex *vec){ 
  long unsigned int j;
  long unsigned int irght,ilft,ihfbit;
  long unsigned int isite1,isite2;
  long unsigned int is1_up,is2_up,is1_down,is2_down;
  long unsigned int iexchg, off;
  int num1_up,num2_up;
  int num1_down,num2_down; 
  long unsigned int ibit1_up,ibit2_up,ibit1_down,ibit2_down; 
  double complex spn_z;
  double complex spn;
  long unsigned i_max;    
  i_max=X->Check.idim_max;

  GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit);
  spn=0.0;
  for(isite1=1;isite1<=X->Def.Nsite;isite1++){
    for(isite2=1;isite2<=X->Def.Nsite;isite2++){
      is1_up=X->Def.Tpow[2*isite1-2];
      is1_down=X->Def.Tpow[2*isite1-1];
      is2_up=X->Def.Tpow[2*isite2-2];
      is2_down=X->Def.Tpow[2*isite2-1];


#pragma omp parallel for reduction(+:spn) default(none) firstprivate(i_max, is1_up, is1_down, is2_up, is2_down, irght, ilft, ihfbit, isite1, isite2) private(ibit1_up, num1_up, ibit2_up, num2_up, ibit1_down, num1_down, ibit2_down, num2_down, spn_z, iexchg, off) shared(vec, list_1, list_2_1, list_2_2)
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
	spn+= conj(vec[j])*vec[j]*spn_z/4.0;
	if(isite1==isite2){
	  spn+=conj(vec[j])*vec[j]*(num1_up+num1_down-2*num1_up*num1_down)/2.0;
	}else{
	  if(ibit1_up!=0 && ibit1_down==0 && ibit2_up==0 &&ibit2_down!=0 ){
	    iexchg= list_1[j]-(is1_up+is2_down);
	    iexchg+=(is2_up+is1_down);
	    GetOffComp(list_2_1, list_2_2, iexchg,  irght, ilft, ihfbit, &off);                    
	    spn+=conj(vec[j])*vec[off]/2.0;
	  }else if(ibit1_up==0 && ibit1_down!=0 && ibit2_up!=0 && ibit2_down==0){
	    iexchg= list_1[j]-(is1_down+is2_up);
	    iexchg+=(is2_down+is1_up);
	    GetOffComp(list_2_1, list_2_2, iexchg,  irght, ilft, ihfbit, &off);
	    spn+=conj(vec[j])*vec[off]/2.0;
	  }
	}
      }
    }
  }
  X->Phys.s2=creal(spn);
}

void totalspin_HubbardGC(struct BindStruct *X,double complex *vec){ 
  long unsigned int j;
  long unsigned int isite1,isite2;
  long unsigned int is1_up,is2_up,is1_down,is2_down;
  long unsigned int iexchg, off;
  int num1_up,num2_up;
  int num1_down,num2_down; 
  long unsigned int ibit1_up,ibit2_up,ibit1_down,ibit2_down,list_1_j; 
  double complex spn_z;
  double complex spn;
  long unsigned int i_max;
    
  i_max=X->Check.idim_max;

  spn=0.0;
  for(isite1=1;isite1<=X->Def.Nsite;isite1++){
    for(isite2=1;isite2<=X->Def.Nsite;isite2++){
      is1_up=X->Def.Tpow[2*isite1-2];
      is1_down=X->Def.Tpow[2*isite1-1];
      is2_up=X->Def.Tpow[2*isite2-2];
      is2_down=X->Def.Tpow[2*isite2-1];


#pragma omp parallel for reduction(+:spn) default(none) firstprivate(i_max, is1_up, is1_down, is2_up, is2_down, isite1, isite2) private(list_1_j, ibit1_up, num1_up, ibit2_up, num2_up, ibit1_down, num1_down, ibit2_down, num2_down, spn_z, iexchg, off) shared(vec)
      for(j=1;j<=i_max;j++){                    
        list_1_j    = j-1;
	ibit1_up= list_1_j&is1_up;
	num1_up=ibit1_up/is1_up;            
	ibit2_up= list_1_j&is2_up;
	num2_up=ibit2_up/is2_up;
            
	ibit1_down= list_1_j&is1_down;
	num1_down=ibit1_down/is1_down;
	ibit2_down= list_1_j&is2_down;
	num2_down=ibit2_down/is2_down;
            
	spn_z=(num1_up-num1_down)*(num2_up-num2_down);
	spn+= conj(vec[j])*vec[j]*spn_z/4.0;
	if(isite1==isite2){
	  spn+=conj(vec[j])*vec[j]*(num1_up+num1_down-2*num1_up*num1_down)/2.0;
	}else{
	  if(ibit1_up!=0 && ibit1_down==0 && ibit2_up==0 &&ibit2_down!=0 ){
	    iexchg= list_1_j-(is1_up+is2_down);
	    iexchg+=(is2_up+is1_down);
            off    = iexchg+1;
	    spn+=conj(vec[j])*vec[off]/2.0;
	  }else if(ibit1_up==0 && ibit1_down!=0 && ibit2_up!=0 && ibit2_down==0){
	    iexchg= list_1_j-(is1_down+is2_up);
	    iexchg+=(is2_down+is1_up);
            off    = iexchg+1;
	    spn+=conj(vec[j])*vec[off]/2.0;
	  }
	}
      }
    }
  }
  X->Phys.s2=creal(spn);
}



void totalspin_Spin(struct BindStruct *X,double complex *vec){ 

  long unsigned int j;
  long unsigned int irght,ilft,ihfbit;
  long unsigned int isite1,isite2;
  long unsigned int is1_up,is2_up;
  long unsigned int iexchg, off;
  int num1_up,num2_up;
  int num1_down,num2_down; 
  long unsigned int ibit1_up,ibit2_up,ibit_tmp,is_up; 
  double complex spn_z;
  double complex spn;
  long unsigned int i_max;
    
  i_max=X->Check.idim_max;
  GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit);
  spn=0.0;
  for(isite1=1;isite1<=X->Def.Nsite;isite1++){
    for(isite2=1;isite2<=X->Def.Nsite;isite2++){
      is1_up      = X->Def.Tpow[isite1-1];
      is2_up      = X->Def.Tpow[isite2-1];
      is_up       = is1_up+is2_up;
      
#pragma omp parallel for reduction(+: spn) default(none) firstprivate(i_max, is_up, is1_up, is2_up, irght, ilft, ihfbit, isite1, isite2) private(ibit1_up, num1_up, ibit2_up, num2_up, num1_down, num2_down, spn_z, iexchg, off, ibit_tmp) shared(list_1, list_2_1, list_2_2, vec)
      for(j=1;j<=i_max;j++){                    
	ibit1_up  = list_1[j]&is1_up;
	num1_up   = ibit1_up/is1_up;            
        num1_down = 1-num1_up;
	ibit2_up  = list_1[j]&is2_up;
	num2_up   = ibit2_up/is2_up;
        num2_down = 1-num2_up;
        
	spn_z  = (num1_up-num1_down)*(num2_up-num2_down);
	spn   += conj(vec[j])*vec[j]*spn_z/4.0;
            
	if(isite1==isite2){
	  spn      += conj(vec[j])*vec[j]/2.0;
	}else{
          ibit_tmp  = (num1_up) ^ (num2_up);
          if(ibit_tmp!=0){
            iexchg  = list_1[j] ^ (is_up);
            GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
	    spn    += conj(vec[j])*vec[off]/2.0;
	  }
	}
      }
    }
  }
  X->Phys.s2=creal(spn);
}


void totalspin_SpinGC(struct BindStruct *X,double complex *vec){ 

  long unsigned int j;
  long unsigned int irght,ilft,ihfbit;
  long unsigned int isite1,isite2;
  long unsigned int is1_up,is2_up;
  long unsigned int iexchg, off;
  int num1_up,num2_up;
  int num1_down,num2_down; 
  long unsigned int ibit1_up,ibit2_up,ibit_tmp,is_up; 
  double complex spn_z;
  double complex spn;
  long unsigned int list_1_j;
  long unsigned int i_max;
    
  i_max=X->Check.idim_max;
  GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit);
    
  spn=0.0;
  for(isite1=1;isite1<=X->Def.Nsite;isite1++){
    for(isite2=1;isite2<=X->Def.Nsite;isite2++){
            
      is1_up      = X->Def.Tpow[isite1-1];
      is2_up      = X->Def.Tpow[isite2-1];
      is_up       = is1_up+is2_up;

#pragma omp parallel for reduction(+: spn) default(none) firstprivate(i_max, is_up, is1_up, is2_up, irght, ilft, ihfbit, isite1, isite2) private(list_1_j, ibit1_up, num1_up, ibit2_up, num2_up, num1_down, num2_down, spn_z, iexchg, off, ibit_tmp) shared(vec)

      for(j=1;j<=i_max;j++){
	list_1_j=j-1;
	ibit1_up  = list_1_j&is1_up;
	num1_up   = ibit1_up/is1_up;            
        num1_down = 1-num1_up;
	ibit2_up  = list_1_j&is2_up;
	num2_up   = ibit2_up/is2_up;
        num2_down = 1-num2_up;
            
	spn_z  = (num1_up-num1_down)*(num2_up-num2_down);
	spn   += conj(vec[j])*vec[j]*spn_z/4.0;
            
	if(isite1==isite2){
	  spn      += conj(vec[j])*vec[j]/2.0;
	}else{
          ibit_tmp  = (num1_up) ^ (num2_up);
          if(ibit_tmp!=0){
            iexchg  = list_1_j ^ (is_up);
            off    = iexchg+1;
	    spn    += conj(vec[j])*vec[off]/2.0;
	  }
	}
      }
    }
  }

  X->Phys.s2=creal(spn);
}
