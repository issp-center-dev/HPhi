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

//Define Mode for mltply
// complex version
// Note: MPI support for SpinlessFermion is not yet implemented in this version.
// The MPI functions are placeholder stubs.

#ifdef MPI
#include "mpi.h"
#endif
#include "Common.h"
#include "mltply.h"
#include "bitcalc.h"
#include "wrapperMPI.h"
#include "mltplyMPISpinlessFermion.h"

/**
 *
 * Hopping term in SpinlessFermion + Canonical ensemble
 * When both site1 and site2 are in the inter process region.
 * (MPI version - placeholder for future implementation)
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void child_general_hopp_Spinless_MPIdouble(
    unsigned long int itrans,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1) {
#ifdef MPI
  double complex dam_pr;
  dam_pr = X_child_general_hopp_Spinless_MPIdouble(X->Def.EDGeneralTransfer[itrans][0],
                                                   X->Def.EDGeneralTransfer[itrans][2],
                                                   X->Def.EDParaGeneralTransfer[itrans],
                                                   X,
                                                   tmp_v0,
                                                   tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}

/**
 *
 * Hopping term in SpinlessFermion + Canonical ensemble
 * When both site1 and site2 are in the inter process region.
 * (MPI version - placeholder for future implementation)
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_child_general_hopp_Spinless_MPIdouble(
    int org_isite1,
    int org_isite2,
    double complex tmp_trans,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1) {
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin, bitdiff, Fsgn;
  unsigned long int idim_max_buf, j, ioff;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;

  mask1 = (int) X->Def.Tpow[org_isite1];
  mask2 = (int) X->Def.Tpow[org_isite2];

  if (mask2 > mask1)
    bitdiff = mask2 - mask1 * 2;
  else
    bitdiff = mask1 - mask2 * 2;
  origin = myrank ^ (mask1 + mask2);

  state1 = origin & mask1;
  state2 = origin & mask2;

  SgnBit((unsigned long int) (origin & bitdiff), &Fsgn);

  if (state1 == 0 && state2 == mask2) {
    trans = -(double) Fsgn * tmp_trans;
  } else if (state1 == mask1 && state2 == 0) {
    trans = -(double) Fsgn * conj(tmp_trans);
    if (X->Large.mode == M_CORR) {
      trans = 0;
    }
  } else
    return 0.0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0)
    exitMPI(-1);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                      list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0)
    exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0)
    exitMPI(-1);

  dam_pr = 0.0;
  if (X->Large.mode == M_MLTPLY) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff) \
  firstprivate(idim_max_buf, trans, X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {
      if (GetOffComp(list_2_1, list_2_2, list_1buf[j],
                     X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == FALSE) {
        ioff = 0;
        continue;
      }
      dmv = trans * v1buf[j];
      tmp_v0[ioff] += dmv;
      dam_pr += conj(tmp_v1[ioff]) * dmv;
    }
  } else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff) \
  firstprivate(idim_max_buf, trans, X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {
      if (GetOffComp(list_2_1, list_2_2, list_1buf[j],
                     X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == FALSE) {
        ioff = 0;
        continue;
      }
      dmv = trans * v1buf[j];
      dam_pr += conj(tmp_v1[ioff]) * dmv;
    }
  }
  return (dam_pr);
#else
  return 0.0;
#endif
}

/**
 *
 * Hopping term in SpinlessFermion + Canonical ensemble
 * When only site2 is in the inter process region.
 * (MPI version - placeholder for future implementation)
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void child_general_hopp_Spinless_MPIsingle(
    unsigned long int itrans,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1) {
#ifdef MPI
  double complex dam_pr;
  dam_pr = X_child_general_hopp_Spinless_MPIsingle(X->Def.EDGeneralTransfer[itrans][0],
                                                   X->Def.EDGeneralTransfer[itrans][2],
                                                   X->Def.EDParaGeneralTransfer[itrans],
                                                   X,
                                                   tmp_v0,
                                                   tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}

/**
 *
 * Hopping term in SpinlessFermion + Canonical ensemble
 * When only site2 is in the inter process region.
 * (MPI version - placeholder for future implementation)
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_child_general_hopp_Spinless_MPIsingle(
    int org_isite1,
    int org_isite2,
    double complex tmp_trans,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1) {
#ifdef MPI
  unsigned long int mask1, state1;
  int mask2, state2, ierr, origin, bit2diff, Fsgn;
  unsigned long int idim_max_buf, j, state1check, bit1diff, ioff, jreal;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;

  mask2 = (int) X->Def.Tpow[org_isite2];
  bit2diff = mask2 - 1;
  origin = myrank ^ mask2;

  state2 = origin & mask2;

  SgnBit((unsigned long int) (origin & bit2diff), &Fsgn);

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0)
    exitMPI(-1);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                      list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0)
    exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0)
    exitMPI(-1);

  mask1 = X->Def.Tpow[org_isite1];
  if (state2 == mask2) {
    trans = -(double) Fsgn * tmp_trans;
    state1check = 0;
  } else if (state2 == 0) {
    state1check = mask1;
    trans = -(double) Fsgn * conj(tmp_trans);
    if (X->Large.mode == M_CORR) {
      trans = 0;
    }
  } else
    return 0.0;

  bit1diff = X->Def.Tpow[X->Def.Nsite - 1] * 2 - mask1 * 2;

  dam_pr = 0.0;
  if (X->Large.mode == M_MLTPLY) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff, jreal, state1) \
  firstprivate(idim_max_buf, trans, X, mask1, state1check, bit1diff, myrank) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0, list_1)
    for (j = 1; j <= idim_max_buf; j++) {
      Fsgn = 1;
      jreal = list_1buf[j];
      state1 = jreal & mask1;
      if (state1 == state1check) {
        SgnBit((unsigned long int) (jreal & bit1diff), &Fsgn);
        if (GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
                       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == FALSE) {
          ioff = 0;
          continue;
        }
        dmv = (double) Fsgn * trans * v1buf[j];
        tmp_v0[ioff] += dmv;
        dam_pr += conj(tmp_v1[ioff]) * dmv;
      }
    }
  } else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff, jreal, state1) \
  firstprivate(idim_max_buf, trans, X, mask1, state1check, bit1diff) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {
      jreal = list_1buf[j];
      state1 = jreal & mask1;
      if (state1 == state1check) {
        SgnBit(jreal & bit1diff, &Fsgn);
        if (GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
                       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == FALSE) {
          ioff = 0;
          continue;
        }
        dmv = (double) Fsgn * trans * v1buf[j];
        dam_pr += conj(tmp_v1[ioff]) * dmv;
      }
    }
  }
  return (dam_pr);
#else
  return 0.0;
#endif
}

/**
 *
 * Hopping term in SpinlessFermion + Canonical ensemble
 * When only site2 is in the inter process region.
 * (Per-site version - iterates over all transfer terms for this site)
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void child_general_hopp_Spinless_MPIsingle_per_site(
    unsigned long int org_isite2,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1) {
#ifdef MPI
  // Iterate over all transfer terms and find ones involving org_isite2
  unsigned long int i;
  double complex dam_pr;

  for (i = 0; i < X->Def.EDNTransfer; i += 2) {
    // Check if this transfer involves org_isite2 as the inter-process site
    if (X->Def.EDGeneralTransfer[i][2] == (int)org_isite2 &&
        X->Def.EDGeneralTransfer[i][0] + 1 <= X->Def.Nsite) {
      dam_pr = X_child_general_hopp_Spinless_MPIsingle(
          X->Def.EDGeneralTransfer[i][0],
          X->Def.EDGeneralTransfer[i][2],
          X->Def.EDParaGeneralTransfer[i],
          X,
          tmp_v0,
          tmp_v1);
      X->Large.prdct += dam_pr;
    }
    if (X->Def.EDGeneralTransfer[i][0] == (int)org_isite2 &&
        X->Def.EDGeneralTransfer[i][2] + 1 <= X->Def.Nsite) {
      dam_pr = X_child_general_hopp_Spinless_MPIsingle(
          X->Def.EDGeneralTransfer[i][2],
          X->Def.EDGeneralTransfer[i][0],
          conj(X->Def.EDParaGeneralTransfer[i]),
          X,
          tmp_v0,
          tmp_v1);
      X->Large.prdct += dam_pr;
    }
  }
#endif
}

/**
 * Placeholder for X_child_general_hopp_Spinless_MPIsingle_per_site
 * (Not yet implemented - requires EDGeneralTransferSingle structure)
 */
double complex X_child_general_hopp_Spinless_MPIsingle_per_site(
    int org_isite1,
    int org_isite2,
    double complex tmp_trans,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1) {
  // Placeholder - will be implemented with EDGeneralTransferSingle support
  return 0.0;
}

/*******************************************************************************
 * SpinlessFermionGC (Grand Canonical) MPI functions
 ******************************************************************************/

/**
 * @brief Hopping term in SpinlessFermionGC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void child_GC_general_hopp_SpinlessFermion_MPIdouble(
    unsigned long int itrans,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1) {
#ifdef MPI
  double complex dam_pr;
  dam_pr = X_child_GC_general_hopp_SpinlessFermion_MPIdouble(
      X->Def.EDGeneralTransfer[itrans][0],
      X->Def.EDGeneralTransfer[itrans][2],
      X->Def.EDParaGeneralTransfer[itrans],
      X,
      tmp_v0,
      tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}

/**
 * @brief Hopping term in SpinlessFermionGC
 * When both site1 and site2 are in the inter process region.
 * (Core function)
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_child_GC_general_hopp_SpinlessFermion_MPIdouble(
    int org_isite1,
    int org_isite2,
    double complex tmp_trans,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1) {
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin, bitdiff, Fsgn;
  unsigned long int idim_max_buf, j;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;

  mask1 = (int) X->Def.Tpow[org_isite1];
  mask2 = (int) X->Def.Tpow[org_isite2];

  if (mask2 > mask1)
    bitdiff = mask2 - mask1 * 2;
  else
    bitdiff = mask1 - mask2 * 2;
  origin = myrank ^ (mask1 + mask2);

  state1 = origin & mask1;
  state2 = origin & mask2;

  SgnBit((unsigned long int) (origin & bitdiff), &Fsgn);

  if (state1 == 0 && state2 == mask2) {
    trans = -(double) Fsgn * tmp_trans;
  } else if (state1 == mask1 && state2 == 0) {
    trans = -(double) Fsgn * conj(tmp_trans);
    if (X->Large.mode == M_CORR || X->Large.mode == M_CALCSPEC) {
      trans = 0;
    }
  } else
    return 0.0;

  // For GC mode, we only need to exchange vector (no list_1)
  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0)
    exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0)
    exitMPI(-1);

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, trans, X) shared(v1buf, tmp_v1, tmp_v0)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {
        dmv = trans * v1buf[j];
        tmp_v0[j] += dmv;
        dam_pr += conj(tmp_v1[j]) * dmv;
      }
    } else {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {
        dmv = trans * v1buf[j];
        dam_pr += conj(tmp_v1[j]) * dmv;
      }
    }
  }
  return (dam_pr);
#else
  return 0.0;
#endif
}

/**
 * @brief Hopping term in SpinlessFermionGC
 * When only site2 is in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void child_GC_general_hopp_SpinlessFermion_MPIsingle(
    unsigned long int itrans,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1) {
#ifdef MPI
  double complex dam_pr;
  dam_pr = X_child_GC_general_hopp_SpinlessFermion_MPIsingle(
      X->Def.EDGeneralTransfer[itrans][0],
      X->Def.EDGeneralTransfer[itrans][2],
      X->Def.EDParaGeneralTransfer[itrans],
      X,
      tmp_v0,
      tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}

/**
 * @brief Hopping term in SpinlessFermionGC
 * When only site2 is in the inter process region.
 * (Core function)
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_child_GC_general_hopp_SpinlessFermion_MPIsingle(
    int org_isite1,
    int org_isite2,
    double complex tmp_trans,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1) {
#ifdef MPI
  unsigned long int mask1, state1;
  int mask2, state2, ierr, origin, bit2diff, Fsgn;
  unsigned long int idim_max_buf, j, state1check, bit1diff, ioff, jreal;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;

  mask2 = (int) X->Def.Tpow[org_isite2];
  bit2diff = mask2 - 1;
  origin = myrank ^ mask2;

  state2 = origin & mask2;

  SgnBit((unsigned long int) (origin & bit2diff), &Fsgn);

  // For GC mode, we only need to exchange vector (no list_1)
  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0)
    exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0)
    exitMPI(-1);

  mask1 = X->Def.Tpow[org_isite1];
  if (state2 == mask2) {
    trans = -(double) Fsgn * tmp_trans;
    state1check = 0;
  } else if (state2 == 0) {
    state1check = mask1;
    trans = -(double) Fsgn * conj(tmp_trans);
    if (X->Large.mode == M_CORR || X->Large.mode == M_CALCSPEC) {
      trans = 0;
    }
  } else
    return 0.0;

  bit1diff = X->Def.Tpow[X->Def.Nsite - 1] * 2 - mask1 * 2;

  dam_pr = 0.0;
  if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff, jreal, state1) \
  firstprivate(idim_max_buf, trans, X, mask1, state1check, bit1diff, myrank) shared(v1buf, tmp_v1, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {
      Fsgn = 1;
      // For GC mode, jreal = j - 1 (direct bit representation)
      jreal = j - 1;
      state1 = jreal & mask1;
      if (state1 == state1check) {
        SgnBit((unsigned long int) (jreal & bit1diff), &Fsgn);
        // For GC mode, ioff is calculated directly
        ioff = (jreal ^ mask1) + 1;
        dmv = (double) Fsgn * trans * v1buf[j];
        tmp_v0[ioff] += dmv;
        dam_pr += conj(tmp_v1[ioff]) * dmv;
      }
    }
  } else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff, jreal, state1) \
  firstprivate(idim_max_buf, trans, X, mask1, state1check, bit1diff) shared(v1buf, tmp_v1, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {
      Fsgn = 1;
      jreal = j - 1;
      state1 = jreal & mask1;
      if (state1 == state1check) {
        SgnBit(jreal & bit1diff, &Fsgn);
        ioff = (jreal ^ mask1) + 1;
        dmv = (double) Fsgn * trans * v1buf[j];
        dam_pr += conj(tmp_v1[ioff]) * dmv;
      }
    }
  }
  return (dam_pr);
#else
  return 0.0;
#endif
}

/**
 * @brief Hopping term in SpinlessFermionGC
 * Per-site version for single MPI process site
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void child_GC_general_hopp_SpinlessFermion_MPIsingle_per_site(
    unsigned long int org_isite2,
    struct BindStruct *X,
    double complex *tmp_v0,
    double complex *tmp_v1) {
#ifdef MPI
  // For SpinlessFermionGC, iterate over all transfer terms for this site
  // Note: X_child_GC_...MPIsingle internally handles both directions
  // (forward and hermitian conjugate), so we only need to call it once
  // per transfer with the inter-process site as org_isite2
  unsigned long int i;
  double complex dam_pr;

  for (i = 0; i < X->Def.EDNTransfer; i += 2) {
    // Case 1: destination site (creation) is the inter-process site
    if (X->Def.EDGeneralTransfer[i][2] == (int)org_isite2 &&
        X->Def.EDGeneralTransfer[i][0] + 1 <= X->Def.Nsite) {
      dam_pr = X_child_GC_general_hopp_SpinlessFermion_MPIsingle(
          X->Def.EDGeneralTransfer[i][0],
          X->Def.EDGeneralTransfer[i][2],
          X->Def.EDParaGeneralTransfer[i],
          X,
          tmp_v0,
          tmp_v1);
      X->Large.prdct += dam_pr;
    }
    // Case 2: source site (annihilation) is the inter-process site
    // Use entry i+1 which has the sites swapped (hermitian conjugate pair)
    else if (X->Def.EDGeneralTransfer[i][0] == (int)org_isite2 &&
             X->Def.EDGeneralTransfer[i][2] + 1 <= X->Def.Nsite) {
      // Entry i+1 has sites swapped: [i+1][0] = [i][2], [i+1][2] = [i][0]
      // So [i+1] has the inter-process site as destination
      dam_pr = X_child_GC_general_hopp_SpinlessFermion_MPIsingle(
          X->Def.EDGeneralTransfer[i + 1][0],
          X->Def.EDGeneralTransfer[i + 1][2],
          X->Def.EDParaGeneralTransfer[i + 1],
          X,
          tmp_v0,
          tmp_v1);
      X->Large.prdct += dam_pr;
    }
  }
#endif
}

/**
 * @brief Compute two-body Green's function <c†_i c_j c†_k c_l> for SpinlessFermionGC with MPI
 *
 * This function handles the case where any of the sites i,j,k,l may be inter-process.
 * It computes the expectation value by iterating over all local states and determining
 * the contribution from each state.
 *
 * For SpinlessFermionGC with MPI:
 * - Sites 0 to Nsite-1 are local (stored in local index j-1)
 * - Sites Nsite to NsiteMPI-1 are inter-process (stored in myrank bits)
 * - Tpow[site] for site < Nsite gives 2^site (mask for local bits)
 * - Tpow[site] for site >= Nsite gives 2^(site-Nsite) (mask for myrank bits)
 *
 * @param org_isite1 Site i (creation operator c†_i)
 * @param org_isite2 Site j (annihilation operator c_j)
 * @param org_isite3 Site k (creation operator c†_k)
 * @param org_isite4 Site l (annihilation operator c_l)
 * @param X BindStruct with calculation parameters
 * @param vec Wavefunction vector
 * @return Expectation value contribution from this process
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_GC_CisAjtCkuAlv_SpinlessFermion_MPI(
    int org_isite1,
    int org_isite2,
    int org_isite3,
    int org_isite4,
    struct BindStruct *X,
    double complex *vec) {
#ifdef MPI
  double complex dam_pr = 0.0;
  unsigned long int j;
  unsigned long int i_max = X->Check.idim_max;
  int Nsite = X->Def.Nsite;

  // Bit masks for each site
  unsigned long int is1 = X->Def.Tpow[org_isite1];  // c†_i (create at i)
  unsigned long int is2 = X->Def.Tpow[org_isite2];  // c_j (annihilate at j)
  unsigned long int is3 = X->Def.Tpow[org_isite3];  // c†_k (create at k)
  unsigned long int is4 = X->Def.Tpow[org_isite4];  // c_l (annihilate at l)

  // Determine which sites are inter-process (site index >= Nsite)
  int site1_interPE = (org_isite1 >= Nsite) ? 1 : 0;
  int site2_interPE = (org_isite2 >= Nsite) ? 1 : 0;
  int site3_interPE = (org_isite3 >= Nsite) ? 1 : 0;
  int site4_interPE = (org_isite4 >= Nsite) ? 1 : 0;

  // Compute the rank flip mask: XOR of inter-process site masks that flip
  // The operator c†_i c_j c†_k c_l flips all four sites (assuming distinct sites)
  unsigned long int rank_flip_mask = 0;
  if (site1_interPE) rank_flip_mask ^= is1;
  if (site2_interPE) rank_flip_mask ^= is2;
  if (site3_interPE) rank_flip_mask ^= is3;
  if (site4_interPE) rank_flip_mask ^= is4;

  // The origin process (for MPI communication) is determined by which rank bits flip
  int origin = myrank ^ (int)rank_flip_mask;

  // Check if inter-process site conditions are satisfied for this rank
  // For the operator c†_i c_j c†_k c_l applied to initial state:
  // - site l must be occupied (to annihilate)
  // - site k must be empty (to create) - checked after c_l
  // - site j must be occupied (to annihilate) - checked after c†_k c_l
  // - site i must be empty (to create) - checked after c_j c†_k c_l
  //
  // For inter-process sites, we need to check these conditions on myrank
  // and track how the bits change through the operations
  //
  // IMPORTANT: Do NOT return early here because we need all processes to
  // participate in MPI communication to avoid deadlock.

  // Track the rank bits through the operator sequence
  int tmp_rank = myrank;
  int rank_sgn = 1;
  int tmp_sgn;
  int rank_valid = 1;  // Flag to track if inter-process conditions are met

  // Apply c_l (if inter-process)
  if (site4_interPE) {
    if ((tmp_rank & (int)is4) == 0) {
      rank_valid = 0;  // site l must be occupied, but it's not
    } else {
      tmp_rank ^= (int)is4;
      SgnBit((unsigned long int)(tmp_rank & (is4 - 1)), &tmp_sgn);
      rank_sgn *= tmp_sgn;
    }
  }

  // Apply c†_k (if inter-process)
  if (site3_interPE && rank_valid) {
    if ((tmp_rank & (int)is3) != 0) {
      rank_valid = 0;  // site k must be empty, but it's occupied
    } else {
      tmp_rank ^= (int)is3;
      SgnBit((unsigned long int)(tmp_rank & (is3 - 1)), &tmp_sgn);
      rank_sgn *= tmp_sgn;
    }
  }

  // Apply c_j (if inter-process)
  if (site2_interPE && rank_valid) {
    if ((tmp_rank & (int)is2) == 0) {
      rank_valid = 0;  // site j must be occupied, but it's not
    } else {
      tmp_rank ^= (int)is2;
      SgnBit((unsigned long int)(tmp_rank & (is2 - 1)), &tmp_sgn);
      rank_sgn *= tmp_sgn;
    }
  }

  // Apply c†_i (if inter-process)
  if (site1_interPE && rank_valid) {
    if ((tmp_rank & (int)is1) != 0) {
      rank_valid = 0;  // site i must be empty, but it's occupied
    } else {
      tmp_rank ^= (int)is1;
      SgnBit((unsigned long int)(tmp_rank & (is1 - 1)), &tmp_sgn);
      rank_sgn *= tmp_sgn;
    }
  }

  // If origin == myrank, final states are local (no MPI exchange needed)
  if (origin == myrank) {
    // No MPI communication needed - can safely skip if rank_valid is false
    if (!rank_valid) return 0.0;

    // Iterate over local states and compute contributions
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(vec) \
  firstprivate(i_max, is1, is2, is3, is4, site1_interPE, site2_interPE, site3_interPE, site4_interPE, rank_sgn, Nsite) private(j)
    for (j = 1; j <= i_max; j++) {
      unsigned long int local_bit = j - 1;
      unsigned long int tmp_bit = local_bit;
      int sgn = rank_sgn, tmp_sgn_local;

      // Apply c_l (if local)
      if (!site4_interPE) {
        if ((tmp_bit & is4) == 0) continue;  // site l must be occupied
        tmp_bit ^= is4;
        SgnBit(tmp_bit & (is4 - 1), &tmp_sgn_local);
        sgn *= tmp_sgn_local;
      }

      // Apply c†_k (if local)
      if (!site3_interPE) {
        if ((tmp_bit & is3) != 0) continue;  // site k must be empty
        tmp_bit ^= is3;
        SgnBit(tmp_bit & (is3 - 1), &tmp_sgn_local);
        sgn *= tmp_sgn_local;
      }

      // Apply c_j (if local)
      if (!site2_interPE) {
        if ((tmp_bit & is2) == 0) continue;  // site j must be occupied
        tmp_bit ^= is2;
        SgnBit(tmp_bit & (is2 - 1), &tmp_sgn_local);
        sgn *= tmp_sgn_local;
      }

      // Apply c†_i (if local)
      if (!site1_interPE) {
        if ((tmp_bit & is1) != 0) continue;  // site i must be empty
        tmp_bit ^= is1;
        SgnBit(tmp_bit & (is1 - 1), &tmp_sgn_local);
        sgn *= tmp_sgn_local;
      }

      // Final local state index
      unsigned long int off_local = tmp_bit;
      dam_pr += sgn * conj(vec[off_local + 1]) * vec[j];
    }
  } else {
    // Need MPI communication - ALL processes must participate to avoid deadlock
    unsigned long int idim_max_buf;
    MPI_Status statusMPI;
    int ierr;

    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(vec, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    // Skip computation if rank conditions are not met (but we still did the communication)
    if (!rank_valid) return 0.0;

    // Iterate over local states and compute contributions
    // The final state's local bits are computed, and we get the bra vector from v1buf
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(vec, v1buf) \
  firstprivate(i_max, idim_max_buf, is1, is2, is3, is4, site1_interPE, site2_interPE, site3_interPE, site4_interPE, rank_sgn, Nsite) private(j)
    for (j = 1; j <= i_max; j++) {
      unsigned long int local_bit = j - 1;
      unsigned long int tmp_bit = local_bit;
      int sgn = rank_sgn, tmp_sgn_local;

      // Apply c_l (if local)
      if (!site4_interPE) {
        if ((tmp_bit & is4) == 0) continue;
        tmp_bit ^= is4;
        SgnBit(tmp_bit & (is4 - 1), &tmp_sgn_local);
        sgn *= tmp_sgn_local;
      }

      // Apply c†_k (if local)
      if (!site3_interPE) {
        if ((tmp_bit & is3) != 0) continue;
        tmp_bit ^= is3;
        SgnBit(tmp_bit & (is3 - 1), &tmp_sgn_local);
        sgn *= tmp_sgn_local;
      }

      // Apply c_j (if local)
      if (!site2_interPE) {
        if ((tmp_bit & is2) == 0) continue;
        tmp_bit ^= is2;
        SgnBit(tmp_bit & (is2 - 1), &tmp_sgn_local);
        sgn *= tmp_sgn_local;
      }

      // Apply c†_i (if local)
      if (!site1_interPE) {
        if ((tmp_bit & is1) != 0) continue;
        tmp_bit ^= is1;
        SgnBit(tmp_bit & (is1 - 1), &tmp_sgn_local);
        sgn *= tmp_sgn_local;
      }

      // Final local state index (on origin process)
      unsigned long int off_local = tmp_bit;
      if (off_local + 1 <= idim_max_buf) {
        // vec[j] is ψ[n] on this process
        // v1buf[off_local + 1] is ψ[m] from origin process
        dam_pr += sgn * conj(v1buf[off_local + 1]) * vec[j];
      }
    }
  }

  return dam_pr;
#else
  return 0.0;
#endif
}
