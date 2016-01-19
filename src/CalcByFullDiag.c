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
#include "CalcByFullDiag.h"
#include "wrapperMPI.h"

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
  fprintf(stdoutMPI, cLogFullDiag_SetHam_Start);
  makeHam(&(X->Bind));
  fprintf(stdoutMPI, cLogFullDiag_SetHam_End);

  fprintf(stdoutMPI,cLogFullDiag_Start);
  lapack_diag(&(X->Bind));
  fprintf(stdoutMPI,cLogFullDiag_End);

  X->Bind.Def.St=0;
  fprintf(stdoutMPI, cLogFullDiag_ExpecValue_Start);
  phys(&(X->Bind));
  fprintf(stdoutMPI, cLogFullDiag_ExpecValue_End);
  output(&(X->Bind));
  fprintf(stdoutMPI, cLogFinish);
  return 0;
}
