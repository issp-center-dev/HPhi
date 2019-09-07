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
/*-------------------------------------------------------------
 *[ver.2009.05.25]
 *-------------------------------------------------------------
 * Copyright (C) 2007-2009 Daisuke Tahara. All rights reserved.
 * Copyright (C) 2009-     Takahiro Misawa. All rights reserved.
 * Some functions are added by TM.
 *-------------------------------------------------------------*/
/*=================================================================================================*/

#ifndef MATRIchild_LAPACK_MAGMA_
#define MATRIchild_LAPACK_MAGMA_


#include <stdio.h>
#include <stdlib.h>

#include "magma_v2.h"
#include "magma_operators.h"
#include <complex.h>

int diag_magma_cmp(int xNsize, double complex **A,
                   double complex *r,double complex **vec, int ngpu);

int diag_magma_real(int xNsize, double **A,
                    double *r, double **vec, int ngpu);

#endif
