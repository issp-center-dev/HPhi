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
#include "sgn.h"

/** 
 * 
 * 
 * @param X 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void sgn(struct BindStruct *X){
  long unsigned int i,j,num,div;
  printf("%s", cProStartCalcSgn);
  for(i=0;i< X->Check.sdim;i++){
    num=0;
    for(j=0;j<=X->Def.Nsite-1;j+=1){
      div=i & X->Def.Tpow[j];
      div=div/X->Def.Tpow[j];
      num+=div;       
      //printf("ADEBUG: i=%d j=%d div=%d num=%d: Tpow[%d]=%ld\n",i,j,div,num,j,X->Def.Tpow[j]);  
    }
    
    if(num%2==1){
      list_3[i]=-1;
    }else{
      list_3[i]=1;
    }  
    //printf("CDEBUG: i=%ld list_3=%d: sdim=%ld\n",i,list_3[i],X->Check.sdim);  
  }
  printf("%s", cProEndCalcSgn);
}    
    
