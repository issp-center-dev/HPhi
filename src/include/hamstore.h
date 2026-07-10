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

#ifndef HPHI_HAMSTORE_H
#define HPHI_HAMSTORE_H

#include <assert.h>
#include "global.h"

/**
 * Storage abstraction for the dense FullDiag Hamiltonian
 * (design doc section 3, phase 2).
 *
 * Replicated mode (iHamPanelActive == 0): writes go to the global
 * Ham[i][j] (1-based, as before).
 * Distributed-panel mode (Solver 3, nproc > 1): each rank stores only
 * its owned column range [HamColBegin, HamColEnd] (1-based, inclusive)
 * in the contiguous column-major panel Ham_local; element (i, j) maps
 * to Ham_local[(j - HamColBegin) * HamPanelLd + (i - 1)].
 *
 * Every Hamiltonian write in the generation code MUST go through
 * AddHamElem (or be guarded by HAM_OWNED_COL); a direct Ham[i][j]
 * write silently corrupts nothing in panel mode (Ham is not allocated)
 * but crashes on NULL — the debug assert below catches ownership bugs
 * before that.
 */

static inline int HamOwnedCol(long int jcol) {
  return !iHamPanelActive || (jcol >= HamColBegin && jcol <= HamColEnd);
}
#define HAM_OWNED_COL(jcol) HamOwnedCol((long int)(jcol))

#define AddHamElem(irow, jcol, val)                                    \
  do {                                                                 \
    long int hs_i_ = (long int)(irow);                                 \
    long int hs_j_ = (long int)(jcol);                                 \
    if (iHamPanelActive) {                                             \
      assert(hs_j_ >= HamColBegin && hs_j_ <= HamColEnd);              \
      Ham_local[(hs_j_ - HamColBegin) * HamPanelLd + (hs_i_ - 1)]      \
        += (val);                                                      \
    } else {                                                           \
      Ham[hs_i_][hs_j_] += (val);                                      \
    }                                                                  \
  } while (0)

#endif /* HPHI_HAMSTORE_H */
