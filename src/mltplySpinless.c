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
#include "mltplySpinless.h"
#include "mltplyMPISpinlessFermion.h"
#include "mltplyMPIBatched.h"
#include "wrapperMPI.h"
#include "CalcTime.h"
#include "mltplyCommon.h"
#include "DefCommon.h"

#ifdef MPI
// Static storage for batched transfers (initialized once per calculation)
static MPIBatchedTransfers batched_transfers_GC = {0, NULL, 0};
static int batched_transfers_GC_initialized = 0;
static MPIBatchedTransfers batched_transfers_canonical = {0, NULL, 0};
static int batched_transfers_canonical_initialized = 0;
#endif

double complex child_general_hopp_Spinless
    (
        double complex *tmp_v0,
        double complex *tmp_v1,
        struct BindStruct *X,
        double complex trans
    );

int child_general_int_GetInfo_Spinless
    (
        const int iInterAll,
        struct BindStruct *X,
        long unsigned int isite1,
        long unsigned int isite2,
        long unsigned int isite3,
        long unsigned int isite4,
        double complex tmp_V
    );

double complex CisAjt_Hermite(
    long unsigned int j,
    double complex *tmp_v0,
    double complex *tmp_v1,
    struct BindStruct *X,
    long unsigned int is1_spin,
    long unsigned int is2_spin,
    long unsigned int sum_spin,
    long unsigned int diff_spin,
    double complex *tmp_V
);

/**
 * @file   mltply.c
 *
 * @brief  Multiplying the wavefunction by the Hamiltonian. @f$ H v_1@f$.
 *
 * @version 0.2
 * @details add function to treat the case of generalspin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */


/**
 * @brief Matrix-vector multiplication for SpinlessFermion/SpinlessFermionGC models
 *
 * Computes tmp_v0 += H * tmp_v1, where H contains hopping (transfer) and
 * interaction (CoulombInter, etc.) terms for spinless fermion systems.
 *
 * Mode branching (X->Large.mode):
 * - M_MLTPLY: Full H|psi> calculation with batched MPI optimization
 *   - Uses InitializeMPIBatchedTransfers to group transfers by MPI partner
 *   - Calls batched functions to reduce MPI_Sendrecv overhead
 * - M_CORR/M_CALCSPEC: Unbatched per-transfer processing
 *   - Used when computing correlation functions or spectrum
 *
 * Transfer classification (site numbering is 1-based internally):
 * - MPIdouble: Both sites > Nsite (both inter-process)
 *   - Timer 611, requires full MPI exchange
 * - MPIsingle: One site > Nsite (one local, one inter-process)
 *   - Timer 612, partial MPI exchange
 * - Local: Both sites <= Nsite
 *   - Timer 613, no MPI communication
 *
 * Note: Transfers are stored in pairs (i, i+1) for Hermitian conjugates,
 * so the loop increments by 2.
 *
 * @param X Struct containing Hamiltonian parameters and calculation mode [in]
 * @param tmp_v0 Output vector: updated as v0 += H*v1 [in,out]
 * @param tmp_v1 Input vector [in]
 *
 * @return 0 on success, -1 on error
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int mltplySpinlessFermion(struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1) {
  long unsigned int i, j = 0;
  long unsigned int isite1, isite2;

  double complex dam_pr;
  double complex tmp_trans;

  int ihermite = 0;
  int idx = 0;
  int isGC = (X->Def.iCalcModel == SpinlessFermionGC);
  dam_pr = 0.0;

  //Transfer
  if (X->Large.mode == M_MLTPLY && MPIBatchingEnabled()) {
    StartTimer(600);
    StartTimer(610);
    for (i = 0; i < X->Def.EDNTransfer; i += 2) {
      if (X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite &&
          X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite) {
        StartTimer(611);
        if (isGC) {
          child_GC_general_hopp_SpinlessFermion_MPIdouble(i, X, tmp_v0, tmp_v1);
        } else {
          child_general_hopp_Spinless_MPIdouble(i, X, tmp_v0, tmp_v1);
        }
        StopTimer(611);
      } else {
        continue;
      }
    }

#ifdef MPI
    // Use batched MPI communication for SpinlessFermionGC
    if (isGC) {
      StartTimer(612);
      // Initialize batched transfers on first call
      if (!batched_transfers_GC_initialized) {
        if (InitializeMPIBatchedTransfers_SpinlessFermionGC(X, &batched_transfers_GC) != 0) {
          fprintf(stderr, "Error: Failed to initialize batched MPI transfers\n");
          return -1;
        }
        batched_transfers_GC_initialized = 1;
      }

      // Process all groups with batched communication
      for (int g = 0; g < batched_transfers_GC.num_groups; g++) {
        dam_pr = X_child_GC_general_hopp_SpinlessFermion_MPIsingle_batched(
            &batched_transfers_GC.groups[g], X, tmp_v0, tmp_v1);
        X->Large.prdct += dam_pr;
      }
      StopTimer(612);
    } else {
      // Canonical mode: use batched MPI processing
      StartTimer(612);
      if (!batched_transfers_canonical_initialized) {
        if (InitializeMPIBatchedTransfers_SpinlessFermion(X, &batched_transfers_canonical) != 0) {
          fprintf(stderr, "Error: Failed to initialize batched MPI transfers (canonical)\n");
          return -1;
        }
        batched_transfers_canonical_initialized = 1;
      }

      for (int g = 0; g < batched_transfers_canonical.num_groups; g++) {
        dam_pr = X_child_general_hopp_Spinless_MPIsingle_batched(
            &batched_transfers_canonical.groups[g], X, tmp_v0, tmp_v1);
        X->Large.prdct += dam_pr;
      }
      StopTimer(612);
    }
#else
    // Non-MPI: skip inter-process transfers
    (void)i;
#endif

    for (i = 0; i < X->Def.EDNTransfer; i += 2) {
      if (X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite ||
          X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite) {
        continue;
      } else {
        StartTimer(613);
        isite1 = X->Def.EDGeneralTransfer[i][0] + 1;
        isite2 = X->Def.EDGeneralTransfer[i][2] + 1;
        if (child_general_hopp_GetInfo_Spinless(X, isite1, isite2) != 0) {
          return -1;
        }
        tmp_trans = -X->Def.EDParaGeneralTransfer[i];
        dam_pr = child_general_hopp_Spinless(tmp_v0, tmp_v1, X, tmp_trans);
        X->Large.prdct += dam_pr;
        StopTimer(613);
      }
    }
  } else {
    for (i = 0; i < X->Def.EDNTransfer; i += 2) {
      if (X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite &&
          X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite) {
        StartTimer(611);
        if (isGC) {
          child_GC_general_hopp_SpinlessFermion_MPIdouble(i, X, tmp_v0, tmp_v1);
        } else {
          child_general_hopp_Spinless_MPIdouble(i, X, tmp_v0, tmp_v1);
        }
        StopTimer(611);
      } else if (X->Def.EDGeneralTransfer[i][2] + 1 > X->Def.Nsite) {
        StartTimer(612);
        if (isGC) {
          child_GC_general_hopp_SpinlessFermion_MPIsingle(i, X, tmp_v0, tmp_v1);
        } else {
          child_general_hopp_Spinless_MPIsingle(i, X, tmp_v0, tmp_v1);
        }
        StopTimer(612);
      } else if (X->Def.EDGeneralTransfer[i][0] + 1 > X->Def.Nsite) {
        StartTimer(612);
        if (isGC) {
          child_GC_general_hopp_SpinlessFermion_MPIsingle(i + 1, X, tmp_v0, tmp_v1);
        } else {
          child_general_hopp_Spinless_MPIsingle(i + 1, X, tmp_v0, tmp_v1);
        }
        StopTimer(612);
      } else {
        StartTimer(613);
        isite1 = X->Def.EDGeneralTransfer[i][0] + 1;
        isite2 = X->Def.EDGeneralTransfer[i][2] + 1;
        if (child_general_hopp_GetInfo_Spinless(X, isite1, isite2) != 0) {
          return -1;
        }
        tmp_trans = -X->Def.EDParaGeneralTransfer[i];
        dam_pr = child_general_hopp_Spinless(tmp_v0, tmp_v1, X, tmp_trans);
        X->Large.prdct += dam_pr;
        StopTimer(613);
      }
    }
  }
  StopTimer(610);
  StopTimer(600);

  return 0;
}

/**
 * @brief Apply local hopping term c†_i c_j to wavefunction for spinless fermions
 *
 * Computes contribution from a single hopping term -t * c†_i c_j to H|psi>.
 * Both sites i and j must be local (not inter-process).
 *
 * The hopping is applied using CisAjt_Hermite which handles:
 * - Occupation check (site j must be occupied, site i must be empty)
 * - Fermion sign from anticommutation (via SgnBit on intermediate bits)
 * - Hermitian conjugate contribution
 *
 * Prerequisites:
 * - X->Large.is1_spin, is2_spin, isA_spin, A_spin must be set via
 *   child_general_hopp_GetInfo_Spinless before calling this function
 *
 * Parallel safety:
 * - OpenMP parallel for with reduction on dam_pr
 * - Each thread processes independent j indices
 * - tmp_v0 updates are to distinct indices (no race condition)
 *
 * @param tmp_v0 Output vector (updated atomically per index) [in,out]
 * @param tmp_v1 Input vector [in]
 * @param X Struct with site masks in X->Large [in]
 * @param trans Hopping coefficient -t [in]
 *
 * @return Diagonal contribution <v1|H|v1> for this hopping term
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex child_general_hopp_Spinless
    (
        double complex *tmp_v0,
        double complex *tmp_v1,
        struct BindStruct *X,
        double complex trans
    ) {

  long unsigned int j, isite1, isite2, Asum, Adiff;
  long unsigned int i_max = X->Large.i_max;

  isite1 = X->Large.is1_spin;
  isite2 = X->Large.is2_spin;
  Asum = X->Large.isA_spin;
  Adiff = X->Large.A_spin;

  double complex dam_pr = 0;
  double complex X_dam_pr = 0;
  double complex tmp_trans;
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max, X, Asum, Adiff, isite1, isite2, trans) private(j, X_dam_pr, tmp_trans) shared(tmp_v0, tmp_v1)
  for (j = 1; j <= i_max; j++) {
    tmp_trans = trans;
    X_dam_pr = CisAjt_Hermite(j, tmp_v0, tmp_v1, X, isite1, isite2, Asum, Adiff, &tmp_trans);
    dam_pr += X_dam_pr * tmp_trans;
  }
  return dam_pr;
}

/**
 *
 *
 * @param X
 * @param isite1
 * @param isite2
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int child_general_hopp_GetInfo_Spinless
    (
        struct BindStruct *X,
        unsigned long int isite1,
        unsigned long int isite2
    ) {
  X->Large.is1_spin = X->Def.Tpow[isite1 - 1];
  X->Large.is2_spin = X->Def.Tpow[isite2 - 1];
  X->Large.isA_spin = X->Large.is1_spin + X->Large.is2_spin;

  if (isite1 > isite2) {
    X->Large.A_spin = (X->Def.Tpow[isite1 - 1] - X->Def.Tpow[isite2]);
  } else if (isite1 < isite2) {
    X->Large.A_spin = (X->Def.Tpow[isite2 - 1] - X->Def.Tpow[isite1]);
  }
  return 0;
}

/**
 *
 *
 * @param iInterAll
 * @param X
 * @param isite1
 * @param isite2
 * @param isite3
 * @param isite4
 * @param tmp_V
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int child_general_int_GetInfo_Spinless
    (
        const int iInterAll,
        struct BindStruct *X,
        long unsigned int isite1,
        long unsigned int isite2,
        long unsigned int isite3,
        long unsigned int isite4,
        double complex tmp_V
    ) {
  //int isite1, isite2, isite3, isite4;
  //double complex tmp_V;
  long unsigned int is1_spin, is2_spin, is3_spin, is4_spin;
  long unsigned int A_spin, B_spin;
  long unsigned int isA_spin, isB_spin;

  is1_spin = X->Def.Tpow[isite1 - 1];
  is2_spin = X->Def.Tpow[isite2 - 1];
  if (isite1 > isite2) {
    A_spin = (X->Def.Tpow[isite1 - 1] - X->Def.Tpow[isite2]);
  } else if (isite2 > isite1) {
    A_spin = (X->Def.Tpow[isite2 - 1] - X->Def.Tpow[isite1]);
  } else {
    A_spin = 0;
  }

  is3_spin = X->Def.Tpow[isite3 - 1];
  is4_spin = X->Def.Tpow[isite4 - 1];
  if (isite3 > isite4) {
    B_spin = (X->Def.Tpow[isite3 - 1] - X->Def.Tpow[isite4]);
  } else if (isite3 < isite4) {
    B_spin = (X->Def.Tpow[isite4 - 1] - X->Def.Tpow[isite3]);
  } else {//isite3=isite4
    B_spin = 0;
  }

  isA_spin = is1_spin + is2_spin;
  isB_spin = is3_spin + is4_spin;

  X->Large.is1_spin = is1_spin;
  X->Large.is2_spin = is2_spin;
  X->Large.is3_spin = is3_spin;
  X->Large.is4_spin = is4_spin;
  X->Large.isA_spin = isA_spin;
  X->Large.isB_spin = isB_spin;
  X->Large.A_spin = A_spin;
  X->Large.B_spin = B_spin;
  X->Large.tmp_V = tmp_V;
  X->Large.isite1 = isite1;
  X->Large.isite2 = isite2;
  X->Large.isite3 = isite3;
  X->Large.isite4 = isite4;

  return 0;
}

/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param is1_spin
 * @param is2_spin
 * @param sum_spin
 * @param diff_spin
 * @param tmp_V
 *
 * @return
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex CisAjt_Hermite(
    long unsigned int j,
    double complex *tmp_v0,
    double complex *tmp_v1,
    struct BindStruct *X,
    long unsigned int is1_spin,
    long unsigned int is2_spin,
    long unsigned int sum_spin,
    long unsigned int diff_spin,
    double complex *tmp_V
) {
  long unsigned int ibit_tmp_1, ibit_tmp_2;
  long unsigned int bit, iexchg, off;
  long unsigned int org_bit;
  int sgn;
  double complex dmv, dam_pr;
  dam_pr = 0.0;

  // For GC mode, use (j-1) as the bit representation; for canonical, use list_1[j]
  if (X->Def.iCalcModel == SpinlessFermionGC) {
    org_bit = j - 1;
  } else {
    org_bit = list_1[j];
  }

  ibit_tmp_1 = org_bit & sum_spin;
  if (ibit_tmp_1 == 0 || ibit_tmp_1 == sum_spin) {
    return (0.0);
  }

  ibit_tmp_1 = (org_bit & is1_spin);
  ibit_tmp_2 = (org_bit & is2_spin);
  if (ibit_tmp_1 != 0 && ibit_tmp_2 == 0) {
    *tmp_V = conj(*tmp_V);
  } else if (!(ibit_tmp_1 == 0 && ibit_tmp_2 != 0)) {
    return (0.0);
  }
  bit = org_bit & diff_spin;
  SgnBit(bit, &sgn);
  iexchg = org_bit ^ sum_spin;

  // For GC mode, calculate off directly; for canonical, use GetOffComp
  if (X->Def.iCalcModel == SpinlessFermionGC) {
    off = iexchg + 1;
  } else {
    if (GetOffComp(list_2_1, list_2_2, iexchg, X->Large.irght, X->Large.ilft, X->Large.ihfbit, &off) != TRUE) {
      return (0.0);
    }
  }

  dmv = sgn * tmp_v1[j];
  if (X->Large.mode == M_MLTPLY) {
    tmp_v0[off] += *tmp_V * dmv;
  }
  dam_pr = dmv * conj(tmp_v1[off]);

  return dam_pr;
}
