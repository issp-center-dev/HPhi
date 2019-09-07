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

#ifndef HPHI_MLTPLYHUBBARDCORE_H
#define HPHI_MLTPLYHUBBARDCORE_H

#include "Common.h"

double complex pairhopp_element
(
 long unsigned int j,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex GC_exchange_element
(
 long unsigned int j,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex GC_pairhopp_element
(
 long unsigned int j,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex exchange_element
(
 long unsigned int j,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex CisAisCisAis_element
(
 long unsigned int j,
 long unsigned int isite1,
 long unsigned int isite3,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex CisAisCjtAku_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite3,
                long unsigned int isite4,
                long unsigned int Bsum,
                long unsigned int Bdiff,
                double complex tmp_V,
                double complex *tmp_v0,
                double complex *tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

double complex CisAjtCkuAku_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite2,
                long unsigned int isite3,
                long unsigned int Asum,
                long unsigned int Adiff,
                double complex tmp_V,
                double complex *tmp_v0,
                double complex *tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

double complex CisAjtCkuAlv_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite2,
                long unsigned int isite3,
                long unsigned int isite4,
                long unsigned int Asum,
                long unsigned int Adiff,
                long unsigned int Bsum,
                long unsigned int Bdiff,
                double complex tmp_V,
                double complex *tmp_v0,
                double complex *tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off_2
        );
//[s]Grand canonical
double complex GC_CisAisCisAis_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite3,
                double complex tmp_V,
                double complex *tmp_v0,
                double complex *tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

double complex GC_CisAisCjtAku_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite3,
                long unsigned int isite4,
                long unsigned int Bsum,
                long unsigned int Bdiff,
                double complex tmp_V,
                double complex *tmp_v0,
                double complex *tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

double complex GC_CisAjtCkuAku_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite2,
                long unsigned int isite3,
                long unsigned int Asum,
                long unsigned int Adiff,
                double complex tmp_V,
                double complex *tmp_v0,
                double complex *tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off
        );

double complex GC_CisAjtCkuAlv_element
        (
                long unsigned int j,
                long unsigned int isite1,
                long unsigned int isite2,
                long unsigned int isite3,
                long unsigned int isite4,
                long unsigned int Asum,
                long unsigned int Adiff,
                long unsigned int Bsum,
                long unsigned int Bdiff,
                double complex tmp_V,
                double complex *tmp_v0,
                double complex *tmp_v1,
                struct BindStruct *X,
                long unsigned int *tmp_off_2
        );
//[e]Grand canonical

double complex GC_CisAis
(
 long unsigned int j,
 double  complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int is1_spin,
 double complex tmp_trans
);

double complex GC_AisCis(
        long unsigned int j,
        double complex *tmp_v0,
        double complex *tmp_v1,
        struct BindStruct *X,
        long unsigned int is1_spin,
        double complex tmp_trans
);

int child_CisAis
(
 long unsigned int list_1_j,
 struct BindStruct *X,
 long unsigned int is1_spin
 );

int child_CisAjt
(
 long unsigned int list_1_j,
 struct BindStruct *X,
 long unsigned int is1_spin,
 long unsigned int is2_spin,
 long unsigned int sum_spin,
 long unsigned int diff_spin,
 long unsigned int *tmp_off
 );


int child_GC_CisAjt
(
 long unsigned int list_1_j,
 struct BindStruct *X,
 long unsigned int is1_spin,
 long unsigned int is2_spin,
 long unsigned int sum_spin,
 long unsigned int diff_spin,
 long unsigned int *tmp_off
 );


double complex CisAjt
(
 long unsigned int j,
 double  complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int is1_spin,
 long unsigned int is2_spin,
 long unsigned int sum_spin,
 long unsigned int diff_spin,
 double complex tmp_V
 );


double complex GC_CisAjt
(
 long unsigned int j,
 double  complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int is1_spin,
 long unsigned int is2_spin,
 long unsigned int sum_spin,
 long unsigned int diff_spin,
 double complex tmp_V,
 long unsigned int *tmp_off
);


int general_hopp_GetInfo
(
 struct BindStruct *X,
 unsigned long int isite1,
 unsigned long int isite2,
 unsigned long int sigma1,
 unsigned long int sigma2
 );

int general_int_GetInfo
(
 int iInterAll,
 struct BindStruct *X,
 long unsigned int isite1,
 long unsigned int isite2,
 long unsigned int isite3,
 long unsigned int isite4,
 long unsigned int sigma1,
 long unsigned int sigma2,
 long unsigned int sigma3,
 long unsigned int sigma4,
 double complex tmp_V
 );


int pairhopp_GetInfo
(
 int iPairHopp,
 struct BindStruct *X 
  );

int exchange_GetInfo
(
 int iExchange,
 struct BindStruct *X 
 );


double complex GC_Ajt
(
 long unsigned int j,
 double complex *tmp_v0,
 double complex *tmp_v1,
 long unsigned int is1_spin,
 double complex tmp_V,
 long unsigned int *tmp_off
 );

double complex GC_Cis
(
 long unsigned int j,
 double complex *tmp_v0,
 double complex *tmp_v1,
 long unsigned int is1_spin,
 double complex tmp_V,
 long unsigned int *tmp_off
 );



double complex GC_Ajt
(
 long unsigned int j,
 double complex *tmp_v0,
 double complex *tmp_v1,
 long unsigned int is1_spin,
 double complex tmp_V,
 long unsigned int *tmp_off
 );

int child_Cis
(
 long unsigned int j,
 long unsigned int is1_spin,
 long unsigned int *tmp_off,
 long unsigned int *list_1_org,
 long unsigned int *list_2_1_target,
 long unsigned int *list_2_2_target,
 long unsigned int _irght,
 long unsigned int _ilft,
 long unsigned int _ihfbit
 );



double complex child_Ajt
(
 long unsigned int j,
 long unsigned int is1_spin,
 long unsigned int *tmp_off,
 long unsigned int *list_1_org,
 long unsigned int *list_2_1_target,
 long unsigned int *list_2_2_target,
 long unsigned int _irght,
 long unsigned int _ilft,
 long unsigned int _ihfbit
 );

#endif
