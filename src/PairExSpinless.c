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
/*-------------------------------------------------------------*/
#include "PairExSpinless.h"
#include "mltplySpinless.h"
#include "mltplyMPISpinlessFermion.h"
#include "mltply.h"
#ifdef MPI
#include "common/setmemory.h"
#endif

///
/// \brief Calculating the pair excited state for spinless fermion system
/// \param X [in,out] define list to get and put information of calculation
/// \param tmp_v0 [out] Result v0 = H v1
/// \param tmp_v1 [in] v0 = H v1
/// \returns 0: Normally finished
/// \returns -1: Abnormally finished
/// \author Kazuyoshi Yoshimi
int GetPairExcitedState_SpinlessFermion(
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1
) {
  long int i;
  int isite1, isite2;
  complex double dam_pr, tmp_trans;
  //Transfer
  for (i = 0; i < X->Def.NPairExcitationOperator; i++) {
    tmp_trans = X->Def.ParaPairExcitationOperator[i];

    // Diagonal operator c^+_i c_i = n_i (or 1 - n_i for itype 0): it is diagonal
    // in the occupation basis, so it needs no MPI communication regardless of
    // whether site i is held by this process. Apply it directly; the hopping
    // kernels return 0 for creation and annihilation on the same site.
    if (X->Def.PairExcitationOperator[i][0] == X->Def.PairExcitationOperator[i][2]) {
      int isite0 = X->Def.PairExcitationOperator[i][0]; // 0-based site
      int itype = X->Def.PairExcitationOperator[i][4];
      double complex coef = (itype == 0) ? -tmp_trans : tmp_trans;
      long unsigned int jj, i_max = X->Large.i_max;
      if (isite0 + 1 > X->Def.Nsite) {
        // Inter-process site: occupation is fixed by this process rank.
        int occ = (myrank >> (isite0 - X->Def.Nsite)) & 1;
        double complex num1 = (itype == 0) ? (double)(1 - occ) : (double)occ;
        for (jj = 1; jj <= i_max; jj++) tmp_v0[jj] += coef * num1 * tmp_v1[jj];
      } else {
        // Local site: occupation read from the basis bit pattern.
        long unsigned int is = X->Def.Tpow[isite0];
        int isGC = (X->Def.iCalcModel == SpinlessFermionGC);
        for (jj = 1; jj <= i_max; jj++) {
          long unsigned int bitrep = isGC ? (jj - 1) : list_1[jj];
          long unsigned int ibit = (bitrep & is) / is;
          double complex num1 = (itype == 0) ? (double)(1 - ibit) : (double)ibit;
          tmp_v0[jj] += coef * num1 * tmp_v1[jj];
        }
      }
      continue;
    }

    if (X->Def.PairExcitationOperator[i][0] + 1 > X->Def.Nsite &&
        X->Def.PairExcitationOperator[i][2] + 1 > X->Def.Nsite) {
      dam_pr = X_child_general_hopp_Spinless_MPIdouble(X->Def.PairExcitationOperator[i][0],
                                                       X->Def.PairExcitationOperator[i][2],
                                                       tmp_trans, X,
                                                       tmp_v0,
                                                       tmp_v1);
    } else if (X->Def.PairExcitationOperator[i][2] + 1 > X->Def.Nsite) {
      dam_pr = X_child_general_hopp_Spinless_MPIsingle(X->Def.PairExcitationOperator[i][0],
                                                       X->Def.PairExcitationOperator[i][2],
                                                       tmp_trans, X,
                                                       tmp_v0,
                                                       tmp_v1);
    } else if (X->Def.PairExcitationOperator[i][0] + 1 > X->Def.Nsite) {
      dam_pr = X_child_general_hopp_Spinless_MPIsingle(X->Def.PairExcitationOperator[i][2],
                                                       X->Def.PairExcitationOperator[i][0],
                                                       conj(tmp_trans), X,
                                                       tmp_v0,
                                                       tmp_v1);
    } else {
      isite1 = X->Def.PairExcitationOperator[i][0] + 1;
      isite2 = X->Def.PairExcitationOperator[i][2] + 1;
      if (child_general_hopp_GetInfo_Spinless(X, isite1, isite2) != 0) {
        return -1;
      }
      dam_pr = child_general_hopp_Spinless(tmp_v0, tmp_v1, X, tmp_trans);
    }
    X->Large.prdct += dam_pr;
  }
  return 0;
}
