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

// Define Mode for mltply
// complex version
#include <bitcalc.h>
#include "mltply.h"
#include "mltplySpin.h"
#include "mltplySpinSym.h"
#include "mltplyHubbard.h"
#include "mltplySpinless.h"
#include "wrapperMPI.h"
#include "CalcTime.h"
#include "mltplyCommon.h"
#include "diagonalcalc.h"
#include "symmetry_basis.h"

/**
 * @file   mltply.c
 *
 * @brief  Main entry point for Hamiltonian-vector multiplication H|v1> -> |v0>
 *
 * This is the central routine called by Lanczos, LOBPCG, and other iterative
 * solvers. It dispatches to model-specific implementations:
 * - Hubbard/HubbardGC: mltplyHubbard(), mltplyHubbardGC()
 * - Spin/SpinGC: mltplySpin(), mltplySpinGC()
 * - SpinlessFermion/SpinlessFermionGC: mltplySpinlessFermion()
 * - Kondo: mltplyKondo(), mltplyKondoGC()
 *
 * The computation proceeds in two phases:
 * 1. Diagonal terms: list_Diagonal[j] * v1[j] (on-site energies, interactions)
 * 2. Off-diagonal terms: Hopping, exchange, pair-hopping (model-specific)
 *
 * @version 0.2 Added general spin support
 * @version 0.1
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */


/**
 * @brief Compute H|v1> and accumulate into |v0>
 *
 * Main Hamiltonian application: tmp_v0 += H * tmp_v1
 *
 * Processing order:
 * 1. Initialize bit masks (irght, ilft, ihfbit) for split-index scheme
 * 2. Apply diagonal terms: tmp_v0[j] += list_Diagonal[j] * tmp_v1[j]
 * 3. Dispatch to model-specific off-diagonal routine
 *
 * Mode branching (X->Large.mode):
 * - M_MLTPLY: Full H|v> computation (used by Lanczos/LOBPCG)
 *   - Updates tmp_v0 with H*tmp_v1
 *   - Returns energy expectation <v1|H|v1> in X->Large.prdct
 * - M_CORR: Correlation function mode
 *   - Only computes expectation values, no tmp_v0 update
 *   - Used by expec_cisajs, expec_cisajscktaltdc
 *
 * Timer IDs:
 * - 1: Total mltply time
 * - 100: Diagonal term
 * - 200-600: Off-diagonal terms (model-dependent)
 *
 * @param X Struct containing Hamiltonian and model parameters [in]
 * @param tmp_v0 Output vector: updated as v0 += H*v1 [in,out]
 * @param tmp_v1 Input vector [in]
 *
 * @return 0 on success, -1 on error
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int mltply(struct BindStruct *X, double complex *tmp_v0,double complex *tmp_v1) {
  long unsigned int j=0;
  long unsigned int irght=0;
  long unsigned int ilft=0;
  long unsigned int ihfbit=0;

  double complex dam_pr;

  long unsigned int i_max;

  StartTimer(1);
  i_max = X->Check.idim_max;
  X->Large.prdct = 0.0;
  dam_pr = 0.0;

  if (X->Def.iFlgGeneralSpin == FALSE) {
    if (GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit) != 0) {
      StopTimer(1);
      return -1;
    }
  }
  else if (i_max != 0 && X->Def.iCalcModel == Spin) {
    if (GetSplitBitForGeneralSpin(X->Def.Nsite, &ihfbit, X->Def.SiteToBit) != 0) {
      StopTimer(1);
      return -1;
    }
  }
  X->Large.i_max = i_max;
  X->Large.irght = irght;
  X->Large.ilft = ilft;
  X->Large.ihfbit = ihfbit;
  X->Large.mode = M_MLTPLY;

  if (X->Def.iFlgSymmetryBasis == TRUE) {
    if (X->Sym == NULL || X->Sym->enabled != TRUE) {
      fprintf(stdoutMPI, "Error: symmetry basis is requested but not built.\n");
      StopTimer(1);
      return -1;
    }
    if (X->Def.iCalcModel != Spin) {
      fprintf(stdoutMPI, "Error: symmetry basis mltply supports only Spin in v1.\n");
      StopTimer(1);
      return -1;
    }
    if (mltplySpinSym(X, tmp_v0, tmp_v1) != 0) {
      StopTimer(1);
      return -1;
    }
    X->Large.prdct = SumMPI_dc(X->Large.prdct);
    StopTimer(1);
    return 0;
  }

  StartTimer(100);
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max) shared(tmp_v0, tmp_v1, list_Diagonal)
  for (j = 1; j <= i_max; j++) {
    tmp_v0[j] += (list_Diagonal[j]) * tmp_v1[j];
    dam_pr += (list_Diagonal[j]) * conj(tmp_v1[j]) * tmp_v1[j];
  }
  X->Large.prdct += dam_pr;
  StopTimer(100);
  if (X->Def.iCalcType == TimeEvolution) diagonalcalcForTE(step_i, X, tmp_v0, tmp_v1);
  
  switch (X->Def.iCalcModel) {
  case HubbardGC:
    mltplyHubbardGC(X, tmp_v0, tmp_v1);
    break;
      
  case Hubbard:
  case tJ:
  case tJGC:
  case Kondo:
  case KondoGC:
    mltplyHubbard(X, tmp_v0, tmp_v1);
    break;
      
  case Spin:
    mltplySpin(X, tmp_v0, tmp_v1);
    break;
      
  case SpinGC:
    mltplySpinGC(X, tmp_v0, tmp_v1);
    break;

  case SpinlessFermion:
  case SpinlessFermionGC:
    mltplySpinlessFermion(X, tmp_v0, tmp_v1);
    break;

  default:
    StopTimer(1);
    return -1;
  }
  
  X->Large.prdct = SumMPI_dc(X->Large.prdct);  
  StopTimer(1);
  return 0;
}
