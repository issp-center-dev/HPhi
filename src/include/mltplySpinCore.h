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

#ifndef HPHI_MLTPLYSPINCORE_H
#define HPHI_MLTPLYSPINCORE_H

#include "Common.h"

void child_exchange_spin_element
(
 long unsigned int j,
 double complex       *tmp_v0,
 double complex       *tmp_v1,
 struct BindStruct *X,
  long unsigned int *tmp_off
 );

void GC_child_pairlift_spin_element
(
 long unsigned int j,
 double complex       *tmp_v0,
 double complex       *tmp_v1,
 struct BindStruct *X, 
 long unsigned int *tmp_off
 );

void GC_child_exchange_spin_element
(
 long unsigned int j,
 double complex       *tmp_v0,
 double complex       *tmp_v1,
 struct BindStruct *X,  
  long unsigned int *tmp_off 
 );

int X_child_exchange_spin_element
(
 long unsigned int j,
 struct BindStruct *X,
 long unsigned int isA_up,
 long unsigned int isB_up,
 long unsigned int sigmaA,
 long unsigned int sigmaB,
 long unsigned int *tmp_off
);

//[s]Spin
void child_CisAisCisAis_spin_element
        (
                long unsigned int j,
                long unsigned int isA_up,
                long unsigned int isB_up,
                long unsigned int org_sigma2,
                long unsigned int org_sigma4,
                double complex tmp_V,
                int nstate, double complex **tmp_v0,
                double complex **tmp_v1,
                struct BindStruct *X
        );

void child_CisAisCitAiu_spin_element
        (
                long unsigned int j,
                long unsigned int org_sigma2,
                long unsigned int org_sigma4,
                long unsigned int isA_up,
                long unsigned int isB_up,
                double complex tmp_V,
                int nstate, double complex **tmp_v0,
                double complex **tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

void child_CisAitCiuAiu_spin_element
        (
                long unsigned int j,
                long unsigned int org_sigma2,
                long unsigned int org_sigma4,
                long unsigned int isA_up,
                long unsigned int isB_up,
                double complex tmp_V,
                int nstate, double complex **tmp_v0,
                double complex **tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

void child_CisAitCiuAiv_spin_element
        (
                long unsigned int j,
                long unsigned int org_sigma2,
                long unsigned int org_sigma4,
                long unsigned int isA_up,
                long unsigned int isB_up,
                double complex tmp_V,
                int nstate, double complex **tmp_v0,
                double complex **tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off_2
        );
//[e]Spin

//[s]GC Spin
void GC_child_CisAisCisAis_spin_element
        (
                long unsigned int j,
                long unsigned int isA_up,
                long unsigned int isB_up,
                long unsigned int org_sigma2,
                long unsigned int org_sigma4,
                double complex tmp_V,
                int nstate, double complex **tmp_v0,
                double complex **tmp_v1,
                struct BindStruct *X
        );

void GC_child_CisAisCitAiu_spin_element
        (
                long unsigned int j,
                long unsigned int org_sigma2,
                long unsigned int org_sigma4,
                long unsigned int isA_up,
                long unsigned int isB_up,
                double complex tmp_V,
                int nstate, double complex **tmp_v0,
                double complex **tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

void GC_child_CisAitCiuAiu_spin_element
        (
                long unsigned int j,
                long unsigned int org_sigma2,
                long unsigned int org_sigma4,
                long unsigned int isA_up,
                long unsigned int isB_up,
                double complex tmp_V,
                int nstate, double complex **tmp_v0,
                double complex **tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

void GC_child_CisAitCiuAiv_spin_element
        (
                long unsigned int j,
                long unsigned int org_sigma2,
                long unsigned int org_sigma4,
                long unsigned int isA_up,
                long unsigned int isB_up,
                double complex tmp_V,
                int nstate, double complex **tmp_v0,
                double complex **tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off_2
        );
//[e]GC Spin

int child_general_int_spin_GetInfo
(
 struct BindStruct *X,
 long unsigned int isite1,
 long unsigned int isite2,
 long unsigned int sigma1,
 long unsigned int sigma2,
 long unsigned int sigma3,
 long unsigned int sigma4,
 double complex tmp_V
 );


int child_exchange_spin_GetInfo
(
 int iExchange,
 struct BindStruct *X 
 );

int child_pairlift_spin_GetInfo
(
 int iPairLift,
 struct BindStruct *X 
 );

int X_SpinGC_CisAit(
long unsigned int j,
struct BindStruct *X,
long unsigned int is1_spin,
long unsigned int sigma2,
long unsigned int *tmp_off
);

int X_Spin_CisAit(
        long unsigned int j,
        struct BindStruct *X,
        long unsigned int is1_spin,
        long unsigned int sigma2,
        long unsigned int *list_1_Org_,
        long unsigned int *list_2_1_,
        long unsigned int *list_2_2_,
        long unsigned int *tmp_off
);

int X_Spin_CisAis(
long unsigned int j,
struct BindStruct *X,
long unsigned int is1_spin,
long unsigned int sigma1
);

int X_SpinGC_CisAis(
long unsigned int j,
struct BindStruct *X,
long unsigned int is1_spin,
long unsigned int sigma1
);

#endif
