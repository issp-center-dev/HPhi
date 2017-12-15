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
#include "lapack_diag.h"
#include "matrixlapack.h"
#include "FileIO.h"
#ifdef _MAGMA
#include "matrixlapack_magma.h"
#endif

/** 
 * 
 * @brief performing full diagonalization using lapack
 * @param[in,out] X 
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return 
 */
int lapack_diag(
struct BindStruct *X//!<[inout]
) {

  FILE *fp;
  char sdt[D_FileNameMax] = "";
  int i, j, i_max, xMsize;

  i_max = X->Check.idim_max;

  for (i = 0; i < i_max; i++) {
    for (j = 0; j < i_max; j++) {
      Ham[i][j] = Ham[i + 1][j + 1];
    }
  }
  xMsize = i_max;
  if (X->Def.iNGPU == 0) {
    ZHEEVall(xMsize, Ham, v0, L_vec);
  } else {
#ifdef _MAGMA
    if(diag_magma_cmp(xMsize, Ham, v0, L_vec, X->Def.iNGPU) != 0) {
    return -1;
  }
#else
    fprintf(stdoutMPI, "Warning: MAGMA is not used in this calculation.");
    ZHEEVall(xMsize, Ham, v0, L_vec);
#endif
  }

  strcpy(sdt, cFileNameEigenvalue_Lanczos);
  if (childfopenMPI(sdt, "w", &fp) != 0) {
    return -1;
  }
  for (i = 0; i < i_max; i++) {
    fprintf(fp, " %d %.10lf \n", i, creal(v0[i]));
  }
  fclose(fp);


  return 0;
}
