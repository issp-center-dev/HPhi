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
  unsigned int NSingleExcitationOperator,//!<[in] number of single excitation operators
  int **SingleExcitationOperator,//!<[in] [n][3] = {site, spin, type}
  double complex *ParaSingleExcitationOperator,//!<[in] coefficient of each operator
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
) {
  int iret = 0;
  //tmp_v0
  /* Both the canonical and grand-canonical Hubbard leaves take the operator set
     explicitly, so the ket (X->Def) and the bra (SingleExcitationBra /
     single_ex_bra_*.def) sets are built independently. Spin/SpinGC single
     excitation is rejected below. */
  if (NSingleExcitationOperator == 0) return TRUE;

  switch (X->Def.iCalcModel) {
  case HubbardGC:
    iret = GetSingleExcitedStateHubbardGC(X, NSingleExcitationOperator,
                                          SingleExcitationOperator,
                                          ParaSingleExcitationOperator,
                                          tmp_v0, tmp_v1);
    break;

  case Hubbard:
  case Kondo:
  case KondoGC:
  case tJ:
  case tJGC:
    iret = GetSingleExcitedStateHubbard(X, NSingleExcitationOperator,
                                        SingleExcitationOperator,
                                        ParaSingleExcitationOperator,
                                        tmp_v0, tmp_v1);
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
