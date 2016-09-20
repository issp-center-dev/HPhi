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

static unsigned long int mfint[7];/*for malloc*/

void setmem_HEAD
(
 struct BindStruct *X
 );

void setmem_def
(
 struct BindStruct *X,
 struct BoostList *xBoost
);

int setmem_large
(
 struct BindStruct *X
);

/*
int setmem_Dir_large
(
 struct BindStruct *X,
 unsigned long int *list_1_,
 unsigned long int *list_1_buf_,
 unsigned long int *list_2_1_,
 unsigned long int *list_2_2_,
 unsigned long int *list_jb_,
 unsigned long int *list_2_1_Sz_,
 unsigned long int *list_2_2_Sz_,
 unsigned long int *list_Diagonal_,
 double complex *v0_,
 double complex *v1_,
 double complex *v1buf_,
 double complex *vg_,
 double *alpha_,// -> common used
 double *beta_,//  -> common used
 double complex *vec_
);
*/

void setmem_IntAll_Diagonal
(
 struct DefineList *X
);
