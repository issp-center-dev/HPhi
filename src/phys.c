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
#ifdef _SCALAPACK
#include "matrixscalapack.h"
#endif

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
  long unsigned int i;
  double tmp_N;
#ifdef _SCALAPACK
  double complex *vec_tmp;
  int ictxt, ierr, rank;
  long unsigned int j, i_max;

  i_max = X->Check.idim_max;

  if(use_scalapack){
  fprintf(stdoutMPI, "In scalapack fulldiag, total spin is not calculated !\n");
  vec_tmp = malloc(i_max*sizeof(double complex));
  }
  for (i = 0; i < neig; i++) {
    if (use_scalapack) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      GetEigenVector(i, i_max, Z_vec, descZ_vec, vec_tmp);
      if (rank == 0) {
        for (j = 0; j < i_max; j++) {
          v0[j + 1][i] = vec_tmp[j];
        }
      }
      else {
        for (j = 0; j < i_max; j++) {
          v0[j + 1][i] = 0.0;
        }
      }
    }
    else {
      if (X->Def.iCalcType == FullDiag) {
        for (j = 0; j < i_max; j++) {
          v0[j + 1][i] = v1[j][i];
        }
      }
    }
  }/*for (i = 0; i < neig; i++)*/
#endif

  if (expec_energy_flct(X, neig, v0, v1) != 0) {
    fprintf(stderr, "Error: calc expec_energy.\n");
    exitMPI(-1);
  }
  if (expec_cisajs(X, neig, v0, v1) != 0) {
    fprintf(stderr, "Error: calc OneBodyG.\n");
    exitMPI(-1);
  }
  if (expec_cisajscktaltdc(X, neig, v0, v1) != 0) {
    fprintf(stderr, "Error: calc TwoBodyG.\n");
    exitMPI(-1);
  }
    
  if (X->Def.iCalcType == FullDiag) {
    if (expec_totalspin(X, neig, v1) != 0) {
      fprintf(stderr, "Error: calc TotalSpin.\n");
      exitMPI(-1);
    }
  }

  for (i = 0; i < neig; i++) {
    if (X->Def.iCalcModel == Spin || X->Def.iCalcModel == SpinGC) {
      tmp_N = X->Def.NsiteMPI;
    }
    else {
      tmp_N = X->Phys.num_up[i] + X->Phys.num_down[i];
    }
    if (X->Def.iCalcType == FullDiag) {
      fprintf(stdoutMPI, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf S2=%10lf Doublon=%10lf \n",
        i, X->Phys.energy[i], tmp_N, X->Phys.Sz[i], X->Phys.s2[i], X->Phys.doublon[i]);
    }
    else if (X->Def.iCalcType == CG)
      fprintf(stdoutMPI, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf Doublon=%10lf \n",
        i, X->Phys.energy[i], tmp_N, X->Phys.Sz[i], X->Phys.doublon[i]);
  }
#ifdef _SCALAPACK
  if(use_scalapack) free(vec_tmp);
#endif  
}
