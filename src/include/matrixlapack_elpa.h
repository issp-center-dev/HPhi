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

/* Preferred ELPA block size for the 2D block-cyclic descriptors on the ELPA
   path (design doc section 3). */
#define ELPA_NBLK 64

/* Effective block size for an n x n matrix on an nprow x npcol grid.
   ELPA rejects setups where a process row/column owns no block
   (ELPA_ERROR_SETUP), so the preferred ELPA_NBLK is capped such that
   ceil(n/nblk) >= max(nprow, npcol). Must be used consistently for the
   descriptors, numroc_, and the value passed to diag_elpa_cmp. */
static inline long int ElpaBlockSize(long int n, int nprow, int npcol) {
  long int maxdim = (nprow > npcol) ? nprow : npcol;
  long int nblk = ELPA_NBLK;
  if (maxdim > 0 && nblk * maxdim > n) {
    nblk = n / maxdim;
    if (nblk < 1) nblk = 1;
  }
  return nblk;
}

int diag_elpa_cmp(int xNsize, double complex *A_distr,
                  double complex *Z_distr, double *w,
                  int local_nrows, int local_ncols,
                  int myrow, int mycol, int nblk, int ngpu);
#endif /* _ELPA */

#endif /* HPHI_MATRIXLAPACK_ELPA_H */
