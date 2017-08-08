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

int diagonalcalc
(
 struct BindStruct *X
 );

int SetDiagonalCoulombIntra
(
 long unsigned int isite1,
 double dtmp_V,
 struct BindStruct *X
 );

int SetDiagonalChemi
(
 long unsigned int isite1,
 double dtmp_V,
 long unsigned int spin,
 struct BindStruct *X
 );

int SetDiagonalCoulombInter
(
 long unsigned int isite1,
 long unsigned int isite2,
 double dtmp_V,
 struct BindStruct *X
 );

int SetDiagonalHund
(
 long unsigned int isite1,
 long unsigned int isite2,
 double dtmp_V,
 struct BindStruct *X
 );

int SetDiagonalInterAll
(
 long unsigned int isite1,
 long unsigned int isite2,
 long unsigned int isigma1,
 long unsigned int isigma2,
 double dtmp_V,
 struct BindStruct *X
 );

int diagonalcalcForTE(
        const int _istep,
        struct BindStruct *X,
        double complex *tmp_v0,
        double complex *tmp_v1
);

