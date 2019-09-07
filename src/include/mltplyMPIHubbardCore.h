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
#include <complex.h>
#include "struct.h"

int CheckPE
(
 int isite,
 struct BindStruct *X
 );

int CheckBit_Cis
(
 long unsigned int is1_spin,
 long unsigned int orgbit,
 long unsigned int *offbit
 );

int CheckBit_Ajt
(
 long unsigned int is1_spin,
 long unsigned int orgbit,
 long unsigned int *offbit
 );

int CheckBit_InterAllPE
(
 int isite1,
 int isigma1,
 int isite2,
 int isigma2,
 int isite3,
 int isigma3,
 int isite4,
 int isigma4,
 struct BindStruct *X,
 unsigned long int orgbit,
 unsigned long int *offbit
 );

int CheckBit_PairPE
(
 int isite1,
 int isigma1,
 int isite3,
 int isigma3,
 struct BindStruct *X,
 unsigned long int orgbit
 );

int GetSgnInterAll
(
 unsigned long int isite1,
 unsigned long int isite2,
 unsigned long int isite3,
 unsigned long int isite4,
 int *Fsgn,
 struct BindStruct *X,
 unsigned long int orgbit,
 unsigned long int *offbit
 );

double complex child_GC_CisAisCjtAjt_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 );

double complex child_GC_CisAjtCkuAlv_Hubbard_MPI
(
 int isite1,
 int isigma1,
 int isite2,
 int isigma2,
 int isite3,
 int isigma3,
 int isite4,
 int isigma4,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 );

double complex child_GC_CisAjtCkuAku_Hubbard_MPI
(
 int isite1,
 int isigma1,
 int isite2,
 int isigma2,
 int isite3,
 int isigma3,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 );

double complex child_GC_CisAisCjtAku_Hubbard_MPI
(
 int isite1,
 int isigma1,
 int isite3,
 int isigma3,
 int isite4,
 int isigma4,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 );

double complex child_GC_CisAis_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
);

double complex child_GC_CisAjt_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite2,
 int org_ispin2,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
);

double complex child_CisAisCjtAjt_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 );

double complex child_CisAjtCkuAlv_Hubbard_MPI
(
 int isite1,
 int isigma1,
 int isite2,
 int isigma2,
 int isite3,
 int isigma3,
 int isite4,
 int isigma4,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 );

double complex child_CisAjtCkuAku_Hubbard_MPI
(
 int isite1,
 int isigma1,
 int isite2,
 int isigma2,
 int isite3,
 int isigma3,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 );

double complex child_CisAisCjtAku_Hubbard_MPI
(
 int isite1,
 int isigma1,
 int isite3,
 int isigma3,
 int isite4,
 int isigma4,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 );

double complex child_CisAis_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
);

double complex child_CisAjt_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_isite2,
 int org_ispin2,
 double complex tmp_trans,
 struct BindStruct *X,
 double complex *tmp_v0, 
 double complex *tmp_v1,
 double complex *v1buf,
 long unsigned int *list_1_org,
 long unsigned int *list_1buf_org,
 long unsigned int *list_2_1_target,
 long unsigned int *list_2_2_target
 );

double complex child_CisAjt_MPIsingle
(
 int org_isite1,
 int org_ispin1,
 int org_isite2,
 int org_ispin2,
 double complex tmp_trans,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1,
 double complex *v1buf,
 long unsigned int *list_1_org,
 long unsigned int *list_1buf_org,
 long unsigned int *list_2_1_target,
 long unsigned int *list_2_2_target
 );


double complex child_GC_Cis_MPI
(
 int org_isite,
 int org_ispin,
 double complex tmp_trans,
 double complex *tmp_v0,
 double complex *tmp_v1,
 unsigned long int idim_max,
 double complex *tmp_v1buf,
 unsigned long int *Tpow
 );

double complex child_GC_Ajt_MPI
(
 int org_isite,
 int org_ispin,
 double complex tmp_trans,
 double complex *tmp_v0,
 double complex *tmp_v1,
 unsigned long int idim_max,
 double complex *tmp_v1buf,
 long unsigned int *Tpow
 );

double complex child_Cis_MPI
(
 int org_isite,
 unsigned int org_ispin,
 double complex tmp_trans,
 double complex *tmp_v0,
 double complex *tmp_v1,
 double complex *tmp_v1buf,
 unsigned long int idim_max,
 long unsigned int *Tpow,
 long unsigned int *list_1_org,
 long unsigned int *list_1buf_org,
 long unsigned int *list_2_1_target,
 long unsigned int *list_2_2_target,
 long unsigned int _irght,
 long unsigned int _ilft,
 long unsigned int _ihfbit
 );

double complex child_Ajt_MPI
(
 int org_isite,
 unsigned int org_ispin,
 double complex tmp_trans,
 double complex *tmp_v0,
 double complex *tmp_v1,
 double complex *tmp_v1buf,
 unsigned long int idim_max,
 long unsigned int *Tpow,
 long unsigned int *list_1_org,
 long unsigned int *list_1buf_org,
 long unsigned int *list_2_1_target,
 long unsigned int *list_2_2_target,
 long unsigned int _irght,
 long unsigned int _ilft,
 long unsigned int _ihfbit
 );
