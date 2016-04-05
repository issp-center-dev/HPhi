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
#include "output.h"
#include "FileIO.h"
#include "wrapperMPI.h"
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
    double tmp_N;
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

     if(X->Def.iCalcModel==Spin || X->Def.iCalcModel==SpinGC){
	tmp_N =X->Def.Nsite;
      }
      else{
	tmp_N  = X->Phys.num_up + X->Phys.num_down;
      }
    
    fprintf(fp,"  <H>         <N>        <Sz>       <S2>       <D> \n");
    for(i=0;i<i_max;i++){
      fprintf(fp," %10lf %10lf %10lf %10lf %10lf\n",X->Phys.all_energy[i],tmp_N, X->Phys.all_sz[i],X->Phys.all_s2[i],X->Phys.all_doublon[i]);
    }
    fclose(fp);
  }
  
  return 0;
}

int outputHam(struct BindStruct *X){
  //Output Ham
  long int i=0;
  long int j=0;
  long int imax = X->Check.idim_max;
  long int ihermite=0;

  FILE *fp;
  char sdt[D_FileNameMax];
  sprintf(sdt,cFileNamePhys_FullDiag_Ham, X->Def.CDataFileHead);
  if(childfopenMPI(sdt,"w",&fp)!=0){
    return -1;
  }

  fprintf(fp, "#Ham#MatrixMarket matrix coordinate complex hermitian\n");
  fprintf(fp, "#Ham# put cout %d %d \n",imax,imax);
  ihermite=0;
  for (i=1; i<=imax; i++){
    for (j=1; j<=i; j++){
      if(cabs(Ham[i][j])>1.0e-13)
        fprintf(fp, "#Ham# %d %d %lf %lf\n",i,j,creal(Ham[i][j]),cimag(Ham[i][j]));
        ihermite += 1;
    }
  }
  fprintf(fp, "#Ham# count %d\n",ihermite);
  fclose(fp);
  return 0;
}