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
/*-------------------------------------------------------------*/

#include "Common.h"
#include "readdef.h"
#include "StdFace_main.h"
#include "wrapperMPI.h"

/** 
 * @brief Main program for HPhi
 * 
 * @param argc 
 * @param argv 
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @return 
 */
int main(int argc, char* argv[]){

  int mode=0;
  char cFileListName[D_FileNameMax];
  if(JudgeDefType(argc, argv, &mode)!=0){
    return -1;
  }  

  //MakeDirectory for output
  struct stat tmpst;
  if(stat(cParentOutputFolder,&tmpst)!=0){
    if(mkdir(cParentOutputFolder, 0777)!=0){
      fprintf(stdoutMPI, "%s", cErrOutput);
      return -1;
    }
  }

  strcpy(cFileListName, argv[2]);
  
  if(mode==STANDARD_MODE || mode == STANDARD_DRY_MODE){
    StdFace_main(argv[2]);
    strcpy(cFileListName, "namelist.def");
    if (mode == STANDARD_DRY_MODE){
      fprintf(stdoutMPI, "Dry run is Finished. \n\n");
      return 0;
    }
  }

  setmem_HEAD(&X.Bind);
  //  if(ReadDefFileNInt(argv[1], &(X.Bind.Def))!=0){
  if(ReadDefFileNInt(cFileListName, &(X.Bind.Def))!=0){
    fprintf(stdoutMPI, "%s", cErrDefFile);
    return (-1);
  }
  if(X.Bind.Def.nvec < X.Bind.Def.k_exct){
    fprintf(stdoutMPI, "%s", cErrnvec);
    fprintf(stdoutMPI, cErrnvecShow, X.Bind.Def.nvec, X.Bind.Def.k_exct);
    return (-1);
  }	  
  fprintf(stdoutMPI, "Definition files are correct.\n");
  
  
  /*ALLOCATE-------------------------------------------*/
  setmem_def(&X.Bind);
  /*-----------------------------------------------------*/

  /*Read Def files.*/
  if(ReadDefFileIdxPara(&(X.Bind.Def))!=0){
    fprintf(stdoutMPI, "%s", cErrIndices);
    return (-1);
  }
  else{
    if(!check(&(X.Bind))==0){
      return (-1);
    }
  }
  
  /*LARGE VECTORS ARE ALLOCATED*/
  if(!setmem_large(&X.Bind)==0){
    fprintf(stdoutMPI, cErrLargeMem, iErrCodeMem);
    return (-1);
  }
  /*Set convergence Factor*/
  SetConvergenceFactor(&(X.Bind.Def));
  /*---------------------------*/
  HPhiTrans(&(X.Bind));
  
  switch (X.Bind.Def.iCalcModel){
  case HubbardGC:
  case Hubbard:
  case Kondo:
  case KondoGC:
    sgn(&(X.Bind));
    break;
 
  case Spin:     // not having sign
  case SpinGC:     
    break;
    
  default :
    break;
  }
  
  if(!sz(&(X.Bind))==0){
    return -1;
  }

  if(X.Bind.Def.WRITE==1){
    output_list(&(X.Bind));
    return 0;
  }

  diagonalcalc(&(X.Bind));
  
  //Start Calculation
  switch (X.Bind.Def.iCalcType){
  case Lanczos:
    if(!CalcByLanczos(&X)==0){
      return -1;
    }
    break;
  case FullDiag:
    if(!CalcByFullDiag(&X)==0){
      return -1;
    }
    break;
  case TPQCalc:
    if(!CalcBySSM(NumAve, ExpecInterval, &X)==0){
      return -1;
    }
    break;
  default:
    return -1;
  }  
  
  return 0;
}
