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
#include "output_list.h"
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
int output_list(struct BindStruct *X){
  
  FILE *fp;
  char sdt[D_FileNameMax];
  int i,i_max;
  
  fprintf(stdoutMPI, "%s", cProStartOutputList);
  i_max=X->Check.idim_max;
  switch(X->Def.iCalcModel){
  case HubbardGC:
  case Hubbard:
  case Spin:
  case SpinGC:
  sprintf(sdt, cFileNameListModel, X->Def.Nsite,X->Def.Nup,X->Def.Ndown);
  break;
  case Kondo:
  case KondoGC:
    sprintf(sdt, "ListForKondo_Ns%d_Ncond%d", X->Def.Nsite,X->Def.Ne);
    break;
  default:
    return -1;
  }
  if(childfopenMPI(sdt,"w",&fp)!=0){
    return -1;
  }
  for(i=1;i<=i_max;i++){
    fprintf(fp," %lu \n",list_1[i]);
  }
  fclose(fp);
  
  fprintf(stdoutMPI, "%s", cProEndOutputList);
  return 0;
}
