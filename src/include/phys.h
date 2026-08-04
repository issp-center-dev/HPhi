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
#pragma once
#include "Common.h"

/**
 * @brief Write one FullDiag eigenvector in the existing eigenvec.dat format.
 *
 * FullDiag eigenvectors are global (not split by the MPI site decomposition),
 * so every state is written as ``*_eigenvec_<state>_rank_0.dat`` even when
 * the distributed eigensolver or state-parallel expectation-value path is
 * used. The input array is 0-based and has @p X->Check.idim_max elements;
 * the writer prepends the unused element required by the existing binary
 * format.
 *
 * @param[in] X calculation parameters
 * @param[in] state 0-based eigenstate index
 * @param[in] eigenvector 0-based global eigenvector
 * @return 0 on success, -1 on file-name, open, write, or close failure
 */
int FullDiagOutputEigenvector(const struct BindStruct *X,
                              unsigned long int state,
                              const double complex *eigenvector);

void phys(struct BindStruct *X, unsigned long int neig);
