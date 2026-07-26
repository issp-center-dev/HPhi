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
/**
 * @file rearray_interactions.h
 *
 * @brief Canonical declaration of Rearray_Interactions() -- the two-body
 * Green's-function operator reordering shared by the Mode-1 expec dispatch
 * (src/expec_cisajscktaltdc.c), the ExpecMode-2 mapping extraction
 * (src/expec_trace.c), and the trace-map unit test
 * (test/unit/expec_trace_map_check.c, which links the real definition).
 *
 * The definition lives in src/rearray_interactions.c (moved verbatim out of
 * src/expec_cisajscktaltdc.c in phase 3b Task 2 so the unit test can link the
 * REAL function instead of a drift-prone copy). It reads only
 * X->Def.{CisAjtCkuAlvDC,TBody,FBody,SBody} and performs no MPI.
 */
#pragma once
#include <complex.h>
#include "struct.h"

int Rearray_Interactions(
    int i,
    long unsigned int *org_isite1,
    long unsigned int *org_isite2,
    long unsigned int *org_isite3,
    long unsigned int *org_isite4,
    long unsigned int *org_sigma1,
    long unsigned int *org_sigma2,
    long unsigned int *org_sigma3,
    long unsigned int *org_sigma4,
    double complex *tmp_V,
    struct BindStruct *X,
    int type);
