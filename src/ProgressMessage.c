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

#include "ProgressMessage.h"

const char* cProFinishAlloc ="LARGE ALLOCATE FINISH !\n";
const char* cProEDNTrans ="EDTrans EDNTransfer=%d \n";
const char* cProEDNChemi ="EDTrans EDNChemi=%d \n";
const char* cProStartCalcSgn = "Start: Calc sgn. \n";
const char* cProEndCalcSgn = "End  : Calc sgn. \n";
const char* cProStartCalcSz = "Start: Calc Sz. \n";
const char* cProEndCalcSz = "End  : Calc Sz. \n";
const char* cProStartOutputList = "Start: output list. \n";
const char* cProEndOutputList = "End  : output list. \n";
const char* cProStartCalcDiag = "Start: calc diagaonal components of Hamiltonian. \n";
const char* cProEndCalcDiag = "End  : calc diagaonal components of Hamiltonian. \n";
