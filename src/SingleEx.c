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
/**@file
@brief Functions to compute singly excited state
*/
#include "bitcalc.h"
#include "SingleEx.h"
#include "SingleExHubbard.h"
/**
@brief Calculation of single excited state
Target System: Hubbard, Kondo
@returns TRUE: Normally finished
@returns FALSE: Abnormally finished
@author Kazuyoshi Yoshimi
@version 1.2
*/
int GetSingleExcitedState(
  struct BindStruct *X,//!<define list to get and put information of calcuation
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
) {
  int iret = 0;
  //tmp_v0
  if (X->Def.NSingleExcitationOperator == 0) return TRUE;

  switch (X->Def.iCalcModel) {
  case HubbardGC:
    iret = GetSingleExcitedStateHubbardGC(X, tmp_v0, tmp_v1);
    break;

  case KondoGC:
  case Hubbard:
  case Kondo:
    iret = GetSingleExcitedStateHubbard(X, tmp_v0, tmp_v1);
    break;

  case Spin:
  case SpinGC:
    iret = FALSE;
    break;

  default:
    iret = FALSE;
    break;
  }/*switch (X->Def.iCalcModel)*/
  return iret;
}/*int GetSingleExcitedState*/
