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

int omp_sz(
    long unsigned int ib, 
    long unsigned int ihfbit,
    struct BindStruct *X,
    long unsigned int *list_1_,
    long unsigned int *list_2_1_,
    long unsigned int *list_2_2_,
    long unsigned int *list_jb_
);

int omp_tJ(
    long unsigned int ib, 
    long unsigned int ihfbit,
    struct BindStruct *X,
    long unsigned int *list_1_,
    long unsigned int *list_2_1_,
    long unsigned int *list_2_2_,
    long unsigned int *list_jb_
);

int omp_sz_hacker(
    long unsigned int ib, 
    long unsigned int ihfbit,
    struct BindStruct *X,
    long unsigned int *list_1_,
    long unsigned int *list_2_1_,
    long unsigned int *list_2_2_,
    long unsigned int *list_jb_
);

int omp_sz_Kondo(
    long unsigned int ib, 
    long unsigned int ihfbit,
    struct BindStruct *X,
    long unsigned int *list_1_,
    long unsigned int *list_2_1_,
    long unsigned int *list_2_2_,
    long unsigned int *list_jb_              
);

int omp_sz_KondoNConserved(
    long unsigned int ib, 
    long unsigned int ihfbit,
    struct BindStruct *X,
    long unsigned int *list_1_,
    long unsigned int *list_2_1_,
    long unsigned int *list_2_2_,
    long unsigned int *list_jb_              
);

int omp_sz_KondoGC(
    long unsigned int ib, 
    long unsigned int ihfbit,
    struct BindStruct *X,
    long unsigned int *list_1_,
    long unsigned int *list_2_1_,
    long unsigned int *list_2_2_,
    long unsigned int *list_jb_
);

int omp_sz_spin(
    long unsigned int ib, 
    long unsigned int ihfbit,
    unsigned int N, 
    struct BindStruct *X,
    long unsigned int *list_1_,
    long unsigned int *list_2_1_,
    long unsigned int *list_2_2_,
    long unsigned int *list_jb_
);

int omp_sz_spin_hacker(
    long unsigned int ib, 
    long unsigned int ihfbit,
    unsigned int N, 
    struct BindStruct *X,
    long unsigned int *list_1_,
    long unsigned int *list_2_1_,
    long unsigned int *list_2_2_,
    long unsigned int *list_jb_      
);

int omp_sz_Kondo_hacker(
    long unsigned int ib,
    long unsigned int ihfbit,
    struct BindStruct *X,
    long unsigned int *list_1_,
    long unsigned int *list_2_1_,
    long unsigned int *list_2_2_,
    long unsigned int *list_jb_
);


int omp_sz_GeneralSpin(
    long unsigned int ib, 
    long unsigned int ihfbit,
    struct BindStruct *X,
    long unsigned int *list_1_,
    long unsigned int *list_2_1_,
    long unsigned int *list_2_2_,
    long int *list_2_1_Sz_,
    long int *list_2_2_Sz_,
    long unsigned int *list_jb_
);

long int Binomial(
    int n,
    int k,
    long int **comb,
    int Nsite
);

int sz(
    struct BindStruct *X,
    long unsigned int *list_1_,
    long unsigned int *list_2_1_,
    long unsigned int *list_2_2_
);

int Read_sz(
    struct BindStruct *X,
    const long unsigned int irght,
    const long unsigned int ilft,
    const long unsigned int ihfbit,
    long unsigned int *i_max
);

/*[s] func. for refactoring */
int count_localized_spins(struct BindStruct *X);
void calculate_jb_GeneralSpin(struct BindStruct *X, long unsigned int *list_jb, long int *list_2_1_Sz,long int *list_2_2_Sz, long unsigned int ihfbit,long unsigned int ilftdim,unsigned int N);
void calculate_jb_Spin_m1(struct BindStruct *X, long unsigned int *list_jb, long unsigned int *list_1_, long unsigned int *list_2_1,long unsigned int *list_2_2,\
long unsigned int ihfbit,long unsigned int irght,long unsigned int ilft,long unsigned int ibpatn, unsigned int N);
void calculate_jb_Spin_Old(struct BindStruct *X, long unsigned int *list_jb, long unsigned int ihfbit,unsigned int N);
void calculate_jb_Spin_Hacker(struct BindStruct *X, long unsigned int *list_jb, long unsigned int ihfbit,unsigned int N);
void calculate_jb_Kondo(struct BindStruct *X, long unsigned int *list_jb, long unsigned int ihfbit);
void calculate_jb_KondoGC(struct BindStruct *X, int num_loc, long unsigned int *list_jb, long unsigned int ihfbit);
void calculate_jb_KondoNConserved(struct BindStruct *X, long unsigned int *list_jb, long unsigned int ihfbit);
void calculate_jb_Hubbard(struct BindStruct *X,long unsigned int *list_jb, long unsigned int ihfbit, unsigned int N2);
void calculate_jb_Hubbard_Hacker(struct BindStruct *X,long unsigned int *list_jb, long unsigned int ihfbit, unsigned int N2);
void calculate_jb_HubbardNCoserved(struct BindStruct *X,long unsigned int *list_jb, long unsigned int ihfbit, unsigned int N2);
void calculate_jb_HubbardNCoserved_Hacker(struct BindStruct *X,long unsigned int *list_jb, long unsigned int ihfbit, unsigned int N2);
void calculate_jb_tJ(struct BindStruct *X,long unsigned int *list_jb, long unsigned int ihfbit, unsigned int N2);
void calculate_jb_tJNConserved(struct BindStruct *X,long unsigned int *list_jb, long unsigned int ihfbit, unsigned int N2);
void calculate_jb_tJGC(struct BindStruct *X,long unsigned int *list_jb, long unsigned int ihfbit, unsigned int N2);
/*[e] func. for refactoring */
