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

//Define Mode for mltply
// complex version

#pragma once
#include <complex.h>
#include "struct.h"
#include "mltplyCommon.h"

void child_general_hopp_Spinless_MPIdouble(unsigned long int itrans, struct BindStruct *X,
                                           double complex *tmp_v0, double complex *tmp_v1);

double complex X_child_general_hopp_Spinless_MPIdouble(int org_isite1,
                                                       int org_isite2,
                                                       double complex tmp_trans,
                                                       struct BindStruct *X,
                                                       double complex *tmp_v0,
                                                       double complex *tmp_v1);

void child_general_hopp_Spinless_MPIsingle(unsigned long int itrans, struct BindStruct *X,
                                           double complex *tmp_v0, double complex *tmp_v1);

double complex X_child_general_hopp_Spinless_MPIsingle(int org_isite1,
                                                       int org_isite2,
                                                       double complex tmp_trans,
                                                       struct BindStruct *X,
                                                       double complex *tmp_v0,
                                                       double complex *tmp_v1);

void child_general_hopp_Spinless_MPIsingle_per_site(unsigned long int itrans, struct BindStruct *X,
                                                    double complex *tmp_v0, double complex *tmp_v1);

double complex X_child_general_hopp_Spinless_MPIsingle_per_site(int org_isite1,
                                                                int org_isite2,
                                                                double complex tmp_trans,
                                                                struct BindStruct *X,
                                                                double complex *tmp_v0,
                                                                double complex *tmp_v1);

// SpinlessFermionGC (Grand Canonical) MPI functions
void child_GC_general_hopp_SpinlessFermion_MPIdouble(unsigned long int itrans, struct BindStruct *X,
                                                      double complex *tmp_v0, double complex *tmp_v1);

double complex X_child_GC_general_hopp_SpinlessFermion_MPIdouble(int org_isite1,
                                                                  int org_isite2,
                                                                  double complex tmp_trans,
                                                                  struct BindStruct *X,
                                                                  double complex *tmp_v0,
                                                                  double complex *tmp_v1);

void child_GC_general_hopp_SpinlessFermion_MPIsingle(unsigned long int itrans, struct BindStruct *X,
                                                      double complex *tmp_v0, double complex *tmp_v1);

double complex X_child_GC_general_hopp_SpinlessFermion_MPIsingle(int org_isite1,
                                                                  int org_isite2,
                                                                  double complex tmp_trans,
                                                                  struct BindStruct *X,
                                                                  double complex *tmp_v0,
                                                                  double complex *tmp_v1);

void child_GC_general_hopp_SpinlessFermion_MPIsingle_per_site(unsigned long int org_isite2, struct BindStruct *X,
                                                               double complex *tmp_v0, double complex *tmp_v1);

// Two-body Green's function MPI functions for SpinlessFermionGC
double complex X_CisAjtCkuAlv_SpinlessFermion_MPI(int org_isite1, int org_isite2, int org_isite3, int org_isite4,
                                                   struct BindStruct *X, double complex *vec);

double complex X_GC_CisAjtCkuAlv_SpinlessFermion_MPI(int org_isite1, int org_isite2, int org_isite3, int org_isite4,
                                                      struct BindStruct *X, double complex *vec);
