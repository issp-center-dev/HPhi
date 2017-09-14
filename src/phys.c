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
#include "phys.h"
#include "expec_energy_flct.h"
#include "expec_totalspin.h"
#include "expec_cisajs.h"
#include "expec_cisajscktaltdc.h"
#include "wrapperMPI.h"

/**
 * @file   phys.c
 * @version 0.1, 0.2
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  File for giving a parent function to calculate physical quantities  by full diagonalization method 
 * 
 * 
 */



/** 
 * 
 * @brief A main function to calculate physical quantities by full diagonalization method.
 * @param[in,out] X CalcStruct list for getting and pushing calculation information 
 * @param neig number of eigenvalues
 * @version 0.2
 * @details add output process of calculation results for general spin
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void phys(struct BindStruct *X, //!<[inout]
          unsigned long int neig //!<[in]
) {

  long unsigned int i, j, i_max;
  double tmp_N;

  i_max = X->Check.idim_max;
  for (i = 0; i < neig; i++) {
    for (j = 0; j < i_max; j++) {
      v0[j + 1] = L_vec[i][j];
    }
    X->Phys.eigen_num = i;
    if (expec_energy_flct(X) != 0) {
      fprintf(stderr, "Error: calc expec_energy.\n");
      exitMPI(-1);
    }
    if (expec_cisajs(X, v1) != 0) {
      fprintf(stderr, "Error: calc OneBodyG.\n");
      exitMPI(-1);
    }
    if (expec_cisajscktaltdc(X, v1) != 0) {
      fprintf(stderr, "Error: calc TwoBodyG.\n");
      exitMPI(-1);
    }
    if (X->Def.iCalcType == FullDiag) {
      if (expec_totalspin(X, v1) != 0) {
        fprintf(stderr, "Error: calc TotalSpin.\n");
        exitMPI(-1);
      }
    }

    if (X->Def.iCalcModel == Spin || X->Def.iCalcModel == SpinGC) {
      tmp_N = X->Def.NsiteMPI;
    } else {
      tmp_N = X->Phys.num_up + X->Phys.num_down;
    }

    if (X->Def.iCalcType == FullDiag)
      fprintf(stdoutMPI, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf S2=%10lf Doublon=%10lf \n", i, X->Phys.energy, tmp_N,
              X->Phys.Sz, X->Phys.s2, X->Phys.doublon);
    else if (X->Def.iCalcType == CG)
      fprintf(stdoutMPI, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf Doublon=%10lf \n", i, X->Phys.energy, tmp_N,
              X->Phys.Sz, X->Phys.doublon);
    X->Phys.all_energy[i] = X->Phys.energy;
    X->Phys.all_doublon[i] = X->Phys.doublon;
    X->Phys.all_sz[i] = X->Phys.Sz;
    X->Phys.all_s2[i] = X->Phys.s2;
    X->Phys.all_num_up[i] = X->Phys.num_up;
    X->Phys.all_num_down[i] = X->Phys.num_down;
  }
}
