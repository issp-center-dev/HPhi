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

int child_omp_sz(
		 long unsigned int ib, 
		 long unsigned int ihfbit,
		 int N2,
		 struct BindStruct *X
		 );


int child_omp_sz_Kondo(
		      long unsigned int ib, 
		      long unsigned int ihfbit,
		      int N, 
		      struct BindStruct *X
		      );



int child_omp_sz_KondoGC(
		      long unsigned int ib, 
		      long unsigned int ihfbit,
		      int N, 
		      struct BindStruct *X
		      );



int child_omp_sz_spin(
		      long unsigned int ib, 
		      long unsigned int ihfbit,
		      int N, 
		      struct BindStruct *X
		      );

int child_omp_sz_GeneralSpin(
		      long unsigned int ib, 
		      long unsigned int ihfbit,
		      int N,
		      struct BindStruct *X
			     );

long int Binomial(
	     int n,
	     int k,
	     long int **comb,
	     int Nsite
	     );

int sz(
       struct BindStruct *X
       );

int Read_sz
(
 struct BindStruct *X,
 const long unsigned int irght,
 const long unsigned int ilft,
 const long unsigned int ihfbit,
 long unsigned int *i_max
 );
