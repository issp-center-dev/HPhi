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

//Define Mode for mltply
// complex version

#pragma once
#include <complex.h>
#include "struct.h"

void GC_child_general_hopp_MPIdouble(unsigned long int itrans, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1);

double complex X_GC_child_general_hopp_MPIdouble(int org_isite1, int org_ispin1, int org_isite2, int org_ispin2, double complex tmp_trans,
  struct BindStruct *X ,double complex *tmp_v0, double complex *tmp_v1);

void GC_child_general_hopp_MPIsingle(unsigned long int itrans, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1);

double complex X_GC_child_general_hopp_MPIsingle(int org_isite1, int org_ispin1, int org_isite2, int org_ispin2, double complex tmp_trans,
  struct BindStruct *X ,double complex *tmp_v0, double complex *tmp_v1);


void child_general_hopp_MPIdouble(unsigned long int itrans, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1);

double complex X_child_general_hopp_MPIdouble(int org_isite1, int org_ispin1, int org_isite2, int org_ispin2, double complex tmp_trans,
  struct BindStruct *X ,double complex *tmp_v0, double complex *tmp_v1);

void child_general_hopp_MPIsingle(unsigned long int itrans, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1);

double complex X_child_general_hopp_MPIsingle(int org_isite1, int org_ispin1, int org_isite2, int org_ispin2, double complex tmp_trans,
  struct BindStruct *X ,double complex *tmp_v0, double complex *tmp_v1);


void child_general_int_spin_MPIdouble(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1);

double complex X_child_general_int_spin_MPIdouble(int org_isite1, int org_ispin1, int org_ispin2, int org_isite3, int org_ispin3, int org_ispin4, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1);


double complex X_child_general_int_spin_TotalS_MPIdouble(
						  int org_isite1,
						  int org_isite3,
						  struct BindStruct *X,
						  double complex *tmp_v0,
						  double complex *tmp_v1
							 );

void child_general_int_spin_MPIsingle(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1);

double complex X_child_general_int_spin_MPIsingle(int org_isite1, int org_ispin1, int org_ispin2, int org_isite3, int org_ispin3, int org_ispin4, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1);


void GC_child_general_int_spin_MPIdouble(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1);
  
void GC_child_general_int_spin_MPIsingle(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1);


void GC_child_general_int_GeneralSpin_MPIdouble(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1);

void GC_child_general_int_GeneralSpin_MPIsingle(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1);

void child_general_int_GeneralSpin_MPIdouble(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1);

void child_general_int_GeneralSpin_MPIsingle(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1);

double complex X_GC_child_CisAitCiuAiv_spin_MPIdouble(
int org_isite1, int org_ispin1, int org_ispin2, int org_isite3, int org_ispin3, int org_ispin4, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1);

double complex X_GC_child_CisAisCjuAjv_spin_MPIdouble( int org_isite1, int org_ispin1, int org_isite3, int org_ispin3, int org_ispin4, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1);

double complex X_GC_child_CisAitCjuAju_spin_MPIdouble(int org_isite1, int org_ispin1, int org_ispin2, int org_isite3, int org_ispin3, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1);

double complex X_GC_child_CisAisCjuAju_spin_MPIdouble(int org_isite1, int org_ispin1, int org_isite3, int org_ispin3, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1);

double complex X_GC_child_CisAitCiuAiv_spin_MPIsingle(int org_isite1, int org_ispin1, int org_ispin2, int org_isite3, int org_ispin3, int org_ispin4, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1);

double complex X_GC_child_CisAisCjuAjv_spin_MPIsingle(int org_isite1, int org_ispin1, int org_isite3, int org_ispin3, int org_ispin4, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1);

double complex X_GC_child_CisAitCjuAju_spin_MPIsingle(int org_isite1, int org_ispin1, int org_ispin2, int org_isite3, int org_ispin3, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1);

double complex X_GC_child_CisAisCjuAju_spin_MPIsingle(int org_isite1, int org_ispin1, int org_isite3, int org_ispin3, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1);

double complex X_GC_child_CisAisCjuAju_spin_MPIsingle(int org_isite1, int org_ispin1, int org_isite3, int org_ispin3, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1);

double complex X_GC_child_CisAit_spin_MPIdouble( int org_isite1,  int org_ispin1, int org_ispin2, double complex tmp_trans, struct BindStruct *X , double complex *tmp_v0,  double complex *tmp_v1);

double complex X_GC_child_CisAis_spin_MPIdouble(int org_isite1, int org_ispin1, double complex tmp_trans, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1);


