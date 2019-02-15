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

#ifndef HPHI_MLTPLYHUBBARD_H
#define HPHI_MLTPLYHUBBARD_H

#include "Common.h"

int mltplyHubbard(struct BindStruct *X, int nstate, double complex **tmp_v0,double complex **tmp_v1);

int mltplyHubbardGC(struct BindStruct *X, int nstate, double complex **tmp_v0,double complex **tmp_v1);

double complex GC_child_general_hopp
(
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1,
 struct BindStruct *X,
 double complex trans
 );


double complex GC_child_general_int(
                         int nstate, double complex **tmp_v0,
                         double complex **tmp_v1,
                         struct BindStruct *X
                         );


double complex child_general_int
(
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1,
 struct BindStruct *X
 );


double complex child_general_hopp
(
 double complex       *tmp_v0,
 double complex       *tmp_v1,
 struct BindStruct *X,
 double complex trans
 );

double complex child_exchange
(
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1,
 struct BindStruct *X
 );

double complex child_pairhopp
(
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1,
 struct BindStruct *X
 );

double complex GC_child_exchange
(
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1,
 struct BindStruct *X
 );

double complex GC_child_pairlift
(
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1,
 struct BindStruct *X
 );

double complex GC_child_pairhopp
(
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1,
 struct BindStruct *X
 );

#endif
