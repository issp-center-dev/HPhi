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

void child_GC_CisAisCjuAjv_GeneralSpin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 int org_ispin4,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAitCjuAju_GeneralSpin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAitCjuAjv_GeneralSpin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 int org_ispin4,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

//general spin - single 
void child_GC_CisAisCjuAjv_GeneralSpin_MPIsingle
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 int org_ispin4,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAitCjuAju_GeneralSpin_MPIsingle
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAitCjuAjv_GeneralSpin_MPIsingle
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 int org_ispin4,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAit_GeneralSpin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 double complex tmp_trans,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAis_GeneralSpin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_AisCis_GeneralSpin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAisCjuAju_GeneralSpin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 double complex tmp_trans,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAisCjuAju_GeneralSpin_MPIsingle
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 double complex tmp_trans,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_CisAit_GeneralSpin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 double complex tmp_trans,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1,
 unsigned long int idim_max
 );


void child_GC_CisAitCiuAiv_spin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 int org_ispin4,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAisCjuAjv_spin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 int org_ispin4,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAitCjuAju_spin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAisCjuAju_spin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAitCiuAiv_spin_MPIsingle
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 int org_ispin4,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAisCjuAjv_spin_MPIsingle
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 int org_ispin4,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAitCjuAju_spin_MPIsingle
(
 int org_isite1,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAisCjuAju_spin_MPIsingle
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAisCjuAju_spin_MPIsingle
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAit_spin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 double complex tmp_trans,
 struct BindStruct *X ,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_CisAis_spin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 double complex tmp_trans,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_GC_AisCis_spin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 double complex tmp_trans,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_CisAit_spin_MPIdouble
(
 int org_isite1,
 int org_ispin2,
 double complex tmp_trans,
 struct BindStruct *X /**< [inout]*/,
 int nstate,
 double complex **tmp_v0 /**< [out] Result v0 = H v1*/, 
 double complex **tmp_v1, /**< [in] v0 = H v1*/
 unsigned long int idim_max
 );

void child_CisAisCjuAju_GeneralSpin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_CisAitCjuAjv_GeneralSpin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 int org_ispin4,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

//general spin - single 
void child_CisAisCjuAju_GeneralSpin_MPIsingle
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void child_CisAitCjuAjv_GeneralSpin_MPIsingle
(
 int org_isite1,
 int org_ispin1,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 int org_ispin4,
 double complex tmp_J,
 struct BindStruct *X,
 int nstate, double complex **tmp_v0,
 double complex **tmp_v1
 );

void GC_CisAisCjuAjv_spin_MPIdouble
(
 unsigned long int i_int /**< [in] Interaction ID*/,
 struct BindStruct *X /**< [inout]*/,
 int nstate, double complex **tmp_v0 /**< [out] Result v0 = H v1*/,
 double complex **tmp_v1 /**< [in] v0 = H v1*/
 );

void GC_CisAitCjuAju_spin_MPIdouble
(
 unsigned long int i_int /**< [in] Interaction ID*/,
 struct BindStruct *X /**< [inout]*/,
 int nstate, double complex **tmp_v0 /**< [out] Result v0 = H v1*/,
 double complex **tmp_v1 /**< [in] v0 = H v1*/
 );

void GC_CisAitCiuAiv_spin_MPIdouble
(
 unsigned long int i_int /**< [in] Interaction ID*/,
 struct BindStruct *X /**< [inout]*/,
 int nstate, double complex **tmp_v0 /**< [out] Result v0 = H v1*/,
 double complex **tmp_v1 /**< [in] v0 = H v1*/
 );

void GC_CisAisCjuAjv_spin_MPIsingle
(
 unsigned long int i_int /**< [in] Interaction ID*/,
 struct BindStruct *X /**< [inout]*/,
 int nstate, double complex **tmp_v0 /**< [out] Result v0 = H v1*/,
 double complex **tmp_v1 /**< [in] v0 = H v1*/
 );

void GC_CisAitCjuAju_spin_MPIsingle
(
 unsigned long int i_int /**< [in] Interaction ID*/,
 struct BindStruct *X /**< [inout]*/,
 int nstate, double complex **tmp_v0 /**< [out] Result v0 = H v1*/,
 double complex **tmp_v1 /**< [in] v0 = H v1*/
 );

void GC_CisAitCiuAiv_spin_MPIsingle
(
 unsigned long int i_int /**< [in] Interaction ID*/,
 struct BindStruct *X /**< [inout]*/,
 int nstate, double complex **tmp_v0 /**< [out] Result v0 = H v1*/,
 double complex **tmp_v1 /**< [in] v0 = H v1*/
 );
