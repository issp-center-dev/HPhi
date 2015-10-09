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
#pragma once
#include "Common.h"

int expec_cisajs(
		 struct BindStruct *X,
                 double complex *vec
		 );

double complex child_Green_0(
		     long int j,
		     long int is1,
		     struct BindStruct *X,
		     double complex *vec
		     );

double complex child_Green_0_Spin
(
 long int j,
 long int is1,
 int sigma1,
 struct BindStruct *X,
 double complex *vec
 );

double complex child_Green_1(
		     long int j,
                     long int is,
                     long int is1,
                     long int is2,
                     long int b_sigma,
                     long int irght,
                     long int ilft,
                     long int ihfbit,
		     struct BindStruct *X,
		     double complex *vec
		     );

double complex child_Green_1_Spin(
		     long int j,
                     long int is,
                     long int is1,
                     long int is2,
		     int sigma1,
		     int sigma2,
                     long int irght,
                     long int ilft,
                     long int ihfbit,
		     struct BindStruct *X,
		     double complex *vec
		     );
