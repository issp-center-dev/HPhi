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
#include "input.h"
#include "wrapperMPI.h"
#include "CalcTime.h"

 /// \brief Parent function for FullDiag mode
 /// \param X [in,out] Struct to get information about file header names, dimension of hirbert space, calc type, physical quantities.
 /// \retval TRUE(=1) normally finished.
 /// \retval FALSE(=0) abnormally finished.

int CalcByFullDiag(
       struct EDMainCalStruct *X
       )
{
  int iret=0;
  fprintf(stdoutMPI, "%s", cLogFullDiag_SetHam_Start);
  StartTimer(5100);
  if(X->Bind.Def.iInputHam==FALSE){
    makeHam(&(X->Bind));
  }
  else if(X->Bind.Def.iInputHam==TRUE){
    fprintf(stdoutMPI, "%s", cLogFullDiag_InputHam_Start);
    inputHam(&(X->Bind));
    fprintf(stdoutMPI, "%s", cLogFullDiag_InputHam_End);
  }
  StopTimer(5100);
  fprintf(stdoutMPI, "%s", cLogFullDiag_SetHam_End);
  if(iret != 0) return FALSE;


  if(X->Bind.Def.iOutputHam == TRUE){
    fprintf(stdoutMPI, "%s", cLogFullDiag_OutputHam_Start);
    StartTimer(5500);
    iret=outputHam(&(X->Bind));
    StopTimer(5500);
    fprintf(stdoutMPI, "%s", cLogFullDiag_OutputHam_End);
    if(iret != 0) return FALSE;
    return TRUE;
  }

  fprintf(stdoutMPI, "%s", cLogFullDiag_Start);
  StartTimer(5200);
  iret=lapack_diag(&(X->Bind));
  StopTimer(5200);
  fprintf(stdoutMPI, "%s", cLogFullDiag_End);
  if(iret != 0) return FALSE;

  X->Bind.Def.St=0;
  fprintf(stdoutMPI, "%s", cLogFullDiag_ExpecValue_Start);
  StartTimer(5300);
  phys(&(X->Bind), X->Bind.Check.idim_max);
  StopTimer(5300);
  fprintf(stdoutMPI, "%s", cLogFullDiag_ExpecValue_End);

  StartTimer(5400);
  iret=output(&(X->Bind));
  StopTimer(5400);  
  fprintf(stdoutMPI, "%s", cLogFinish);
  if(iret != 0) return FALSE;

  return TRUE;
}
