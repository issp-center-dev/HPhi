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
#include "CalcByFullDiag.h"
#include "wrapperMPI.h"
#include "CalcTime.h"
/** 
 * 
 * 
 * @param X 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @return 
 */
int CalcByFullDiag(
		   struct EDMainCalStruct *X
		   )
{
  fprintf(stdoutMPI, "%s", cLogFullDiag_SetHam_Start);
  StartTimer(5100);
  makeHam(&(X->Bind));
  StopTimer(5100);
  fprintf(stdoutMPI, "%s", cLogFullDiag_SetHam_End);

  if(X->Bind.Def.iOutputHam == TRUE){
    fprintf(stdoutMPI, "%s", cLogFullDiag_OutputHam_Start);
    StartTimer(5500);
    outputHam(&(X->Bind));
    StopTimer(5500);
    fprintf(stdoutMPI, "%s", cLogFullDiag_OutputHam_End);
    return TRUE;
  }
  fprintf(stdoutMPI, "%s", cLogFullDiag_Start);
  StartTimer(5200);
  lapack_diag(&(X->Bind));
  StopTimer(5200);
  fprintf(stdoutMPI, "%s", cLogFullDiag_End);

  X->Bind.Def.St=0;
  fprintf(stdoutMPI, "%s", cLogFullDiag_ExpecValue_Start);
  StartTimer(5300);
  phys(&(X->Bind), X->Bind.Check.idim_max);
  StopTimer(5300);
  fprintf(stdoutMPI, "%s", cLogFullDiag_ExpecValue_End);
  StartTimer(5400);  
  output(&(X->Bind));
  StopTimer(5400);  
  fprintf(stdoutMPI, "%s", cLogFinish);
  return TRUE;
}
