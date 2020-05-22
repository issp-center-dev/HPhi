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

#ifndef HPHI_MLTPLYSPIN_H
#define HPHI_MLTPLYSPIN_H

#include "Common.h"

int mltplySpin(struct BindStruct *X, double complex *tmp_v0,double complex *tmp_v1);

int mltplyHalfSpin(struct BindStruct *X, double complex *tmp_v0,double complex *tmp_v1);

int mltplyGeneralSpin(struct BindStruct *X, double complex *tmp_v0,double complex *tmp_v1);

void mltplyHalfSpinGC_mini(
  struct BindStruct *X,//!<[inout]
  int site_i,
  int spin_i,
  int site_j,
  int spin_j,
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1//!<[in] Input producted vector
);

int mltplySpinGC(struct BindStruct *X, double complex *tmp_v0,double complex *tmp_v1);

int mltplyHalfSpinGC(struct BindStruct *X, double complex *tmp_v0,double complex *tmp_v1);

int mltplyGeneralSpinGC(struct BindStruct *X, double complex *tmp_v0,double complex *tmp_v1);

int mltplySpinGCBoost(struct BindStruct *X, double complex *tmp_v0,double complex *tmp_v1);

double complex GC_general_int_spin
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex general_int_spin
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );


double complex GC_exchange_spin
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex exchange_spin
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex GC_pairlift_spin
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

#endif
