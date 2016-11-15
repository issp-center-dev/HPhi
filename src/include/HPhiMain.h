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
#pragma once
#include "Common.h"
/*[s] For TPQCalc */
#include "dSFMT.h"
/*[e] For TPQCalc */
#include "../FileIO.c"
#include "../bitcalc.c"
#include "../output.c"
#include "../time.c"

#include "../matrixlapack.c"
#include "../mfmemory.c"
#include "../readdef.c"
#include "../xsetmem.c"
#include "../log.c"

#include "../HPhiTrans.c"
#include "../sgn.c"
#include "../sz.c"
#include "../check.c"
#include "../output_list.c"
#include "../diagonalcalc.c"
#include "../mltply.c"
#include "../expec_cisajs.c"
#include "../expec_cisajscktaltdc.c"
#include "../expec_totalspin.c"
#include "../PowerLanczos.c"
#include "../CalcByLanczos.c"

/* [s] For ExaDiag */
#include "../bisec.c"
#include "../vec12.c"
//#include "mltply_CG.c"
#include "../CG_EigenVector.c"
#include "../Lanczos_EigenValue.c"
#include "../Lanczos_EigenVector.c"
#include "../expec_energy.c"
/* [e] For ExaDiag */

/* [s] For TPQCalc */
#include "../FirstMultiply.c"
#include "../Multiply.c"
#include "../CalcByTPQ.c"
/* [e] For TPQCalc */

/* [s] For AllDiag */
#include "../makeHam.c"
#include "../lapack_diag.c"
#include "../phys.c"
#include "../CalcByFullDiag.c"
/* [e] For AllDiag */

