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
#include "output.h"
#include "FileIO.h"

/** 
 * 
 * 
 * @param X 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int output(struct BindStruct *X){
  
  FILE *fp;
  char sdt[D_FileNameMax];
  long int i,i_max;
  i_max=X->Check.idim_max;
    
  if(X->Def.iCalcType==FullDiag){
    double tmp_N,tmp_Sz;
    tmp_N=X->Phys.num_up+X->Phys.num_down;
    tmp_Sz=X->Phys.sz;
    switch(X->Def.iCalcModel){
    case Spin:
    case Hubbard:
    case Kondo:	
      sprintf(sdt,cFileNamePhys_FullDiag, X->Def.CDataFileHead, X->Def.Nup,X->Def.Ndown);
      break;
    case SpinGC:
    case HubbardGC:
    case KondoGC:	
      sprintf(sdt,cFileNamePhys_FullDiag_GC, X->Def.CDataFileHead);
      break;
    default:
      break;
    }
    if(childfopenMPI(sdt,"w",&fp)!=0){
      return -1;
    }
    fprintf(fp,"  <H>         <N>        <Sz>       <S2>       <D> \n");
    for(i=0;i<i_max;i++){
      fprintf(fp," %10lf %10lf %10lf %10lf %10lf\n",X->Phys.all_energy[i],tmp_N,tmp_Sz,X->Phys.all_s2[i],X->Phys.all_doublon[i]);
    }
    fclose(fp);
  }
  
  return 0;
}
