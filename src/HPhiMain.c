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

#include <sz.h>
#include <HPhiTrans.h>
#include <output_list.h>
#include <diagonalcalc.h>
#include <CalcByLanczos.h>
#include <CalcByFullDiag.h>
#include <CalcByTPQ.h>
#include <CalcSpectrum.h>
#include <check.h>
#include "Common.h"
#include "readdef.h"
#include "StdFace_main.h"
#include "wrapperMPI.h"

/*!
@mainpage

<H2>Introduction</H2>
A numerical solver package for a wide range of quantum lattice models including Hubbard-type itinerant electron hamiltonians, quantum spin models, and Kondo-type hamiltonians for itinerant electrons coupled with quantum spins. The Lanczos algorithm for finding ground states and newly developed Lanczos-based algorithm for finite-temperature properties of these models are implemented for parallel computing. A broad spectrum of users including experimental researchers is cordially welcome.
<HR>
<H2>Developers</H2>
Youhei Yamaji (Quantum-Phase Electronics Center, The University of Tokyo)\n
Takahiro Misawa (Department of Applied Physics, The University of Tokyo)\n
Synge Todo (Department of Physics, The University of Tokyo)\n
Kazuyoshi Yoshimi (Institute for Solid State Physics, The University of Tokyo)\n
Mitsuaki Kawamura (Institute for Solid State Physics, The University of Tokyo)\n
Naoki Kawashima (Institute for Solid State Physics, The University of Tokyo)
<HR>
<H2>Methods</H2>
Lanczos algorithm, thermal pure quantum state, full diagonalization
<HR>
<H2>Target models</H2>
Hubbard model, Heisenberg model, Kondo lattice model, Kitaev model, Kitaev-Heisenberg model, multi-orbital Hubbard model
<HR>
<H2>Link</H2>
https://github.com/QLMS/HPhi
<HR>
<H2>Download</H2>
https://github.com/QLMS/HPhi/releases
<HR>
<H2>Forum</H2>
http://ma.cms-initiative.jp/ja/community/materiapps-messageboard/e5hes9
<HR>
<H2>licence</H2>
<B>GNU GPL version 3</B>\n
This software is developed under the support of "Project for advancement of software usability in materials science" by The Institute for Solid State Physics, The University of Tokyo.\n
*/

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
    FinalizeMPI();
    return TRUE;
  }  

  if (mode == STANDARD_DRY_MODE) {
    myrank = 0;
    nproc = 1;
    stdoutMPI = stdout;
  }
  else InitializeMPI(argc, argv);

  //MakeDirectory for output
  struct stat tmpst;
  if (myrank == 0) {
    if (stat(cParentOutputFolder, &tmpst) != 0) {
      if (mkdir(cParentOutputFolder, 0777) != 0) {
        fprintf(stdoutMPI, "%s", cErrOutput);
        exitMPI(-1);
      }
    }
  }/*if (myrank == 0)*/

  strcpy(cFileListName, argv[2]);
  
  if(mode==STANDARD_MODE || mode == STANDARD_DRY_MODE){
    if (myrank == 0) StdFace_main(argv[2]);
    strcpy(cFileListName, "namelist.def");
    if (mode == STANDARD_DRY_MODE){
      fprintf(stdout, "Dry run is Finished. \n\n");
      return FALSE;
    }
  }

  setmem_HEAD(&X.Bind);
  if(ReadDefFileNInt(cFileListName, &(X.Bind.Def), &(X.Bind.Boost))!=0){
    fprintf(stderr, "%s", cErrDefFile);
    FinalizeMPI();
    return FALSE;
  }
  if (X.Bind.Def.nvec < X.Bind.Def.k_exct){
    fprintf(stdoutMPI, "%s", cErrnvec);
    fprintf(stdoutMPI, cErrnvecShow, X.Bind.Def.nvec, X.Bind.Def.k_exct);
    FinalizeMPI();
    return FALSE;
  }	  
  fprintf(stdoutMPI,  cProFinishDefFiles);
  
  /*ALLOCATE-------------------------------------------*/
  setmem_def(&X.Bind, &X.Bind.Boost);
  /*-----------------------------------------------------*/

  /*Read Def files.*/
  if(ReadDefFileIdxPara(&(X.Bind.Def), &(X.Bind.Boost))!=0){
    fprintf(stdoutMPI, "%s", cErrIndices);
    FinalizeMPI();
    return FALSE;
  }
  
  fprintf(stdoutMPI, cProFinishDefCheck);
  if(check(&(X.Bind))==MPIFALSE){
    FinalizeMPI();
    return FALSE;
  }
  
  
  /*LARGE VECTORS ARE ALLOCATED*/
  if(!setmem_large(&X.Bind)==0){
    fprintf(stdoutMPI, cErrLargeMem, iErrCodeMem);
    exitMPI(-1);
  }

  /*Set convergence Factor*/
  SetConvergenceFactor(&(X.Bind.Def));

  /*---------------------------*/
  if(!HPhiTrans(&(X.Bind))==0){
    exitMPI(-1);
  }

  if(!sz(&(X.Bind))==0){
    exitMPI(-1);
  }

  if(X.Bind.Def.WRITE==1){
    output_list(&(X.Bind));
    FinalizeMPI();
    return FALSE;
  }

  diagonalcalc(&(X.Bind));
  
  //Start Calculation
  switch (X.Bind.Def.iCalcType){
  case Lanczos:    
    if(!CalcByLanczos(&X)==TRUE){
      FinalizeMPI();
      return FALSE;
    }    
    break;
  case FullDiag:
    if(!CalcByFullDiag(&X)==TRUE){
      FinalizeMPI();
      return FALSE;
    }
    break;
    
  case TPQCalc:
    if(!CalcByTPQ(NumAve, ExpecInterval, &X)==TRUE){
      FinalizeMPI();
      return FALSE;
    }
    break;

    case Spectrum:
      if(!CalcSpectrum(&X) == TRUE){
        FinalizeMPI();
        return FALSE;
      }
      break;

    default:
      FinalizeMPI();
      return FALSE;
  }  

  FinalizeMPI();
  return TRUE;
}
