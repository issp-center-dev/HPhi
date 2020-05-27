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

double complex child_Cor_1(long unsigned int j,
                   long unsigned int is_A,
                   long unsigned int is1,
                   long unsigned int is2,
                   long unsigned int is3,
                   long unsigned int b_sigma_A,
                   long unsigned int irght,
                   long unsigned int ilft,
                   long unsigned int ihfbit,
                   struct BindStruct *X,
                   double complex *vec
                   );

double complex child_Cor_2(long unsigned int j,
                   long unsigned int is_B,
                   long unsigned int is1,
                   long unsigned int is3,
                   long unsigned int is4,
                   long unsigned int b_sigma_B,
                   long unsigned int irght,
                   long unsigned int ilft,
                   long unsigned int ihfbit,
                   struct BindStruct *X,
                   double complex *vec
                   );

double complex child_Cor_3(long unsigned int j,
                   long unsigned int is_A,
                   long unsigned int is_B,
                   long unsigned int is1,
                   long unsigned int is2,
                   long unsigned int is3,
                   long unsigned int is4,
                   long unsigned int b_sigma_A,
                   long unsigned int b_sigma_B,
                   long unsigned int irght,
                   long unsigned int ilft,
                   long unsigned int ihfbit,
                   struct BindStruct *X,
                   double complex *vec
                   );


int expec_cisajscktaltdc(
                         struct BindStruct *X,
                         double complex *vec
                         );

void expec_cisajscktaltdc_alldiag_spin(
                                  struct BindStruct *X,
                                  double complex *vec
                                  );

int expec_Threebody_SpinGCHalf(struct BindStruct *X,double complex *vec, FILE **_fp);
int expec_Fourbody_SpinGCHalf(struct BindStruct *X,double complex *vec, FILE **_fp);
