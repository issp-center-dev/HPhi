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
  printf("Start Setting Hamiltonian.\n");
  makeHam(&(X->Bind));
  printf("Finish Setting Hamiltonian.\n");
  lapack_diag(&(X->Bind));
  printf("Diagonalization is done.\n");  
  X->Bind.Def.St=0;
  printf("Start: Calc Expected value.\n");
  phys(&(X->Bind));
  printf("End: Calc Expected value.\n");
  output(&(X->Bind));
  printf("Finish. \n");
  return 0;
}
