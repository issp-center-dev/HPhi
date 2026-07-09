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
#ifndef HPHI_MATRIXLAPACK_ELPA_H
#define HPHI_MATRIXLAPACK_ELPA_H

#ifdef _ELPA
#include <complex.h>

/* ELPA block size used for all 2D block-cyclic descriptors on the ELPA
   path (design doc section 3). */
#define ELPA_NBLK 64

int diag_elpa_cmp(int xNsize, double complex *A_distr,
                  double complex *Z_distr, double *w,
                  int local_nrows, int local_ncols,
                  int myrow, int mycol, int ngpu);
#endif /* _ELPA */

#endif /* HPHI_MATRIXLAPACK_ELPA_H */
