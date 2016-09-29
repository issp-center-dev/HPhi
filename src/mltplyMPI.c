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

#ifdef MPI
#include "mpi.h"
#endif
#include "Common.h"
#include "mltply.h"
#include "bitcalc.h"
#include "wrapperMPI.h"
#include "mltplyMPI.h"


/**
 *
 * Hopping term in Hubbard + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void GC_child_general_hopp_MPIdouble
(
 unsigned long int itrans /**< [in] Transfer ID*/,
 struct BindStruct *X /**< [inout]*/,
 double complex *tmp_v0 /**< [out] Result v0 = H v1*/, 
 double complex *tmp_v1 /**< [in] v0 = H v1*/)
{

#ifdef MPI
  double complex dam_pr=0;
  dam_pr=X_GC_child_general_hopp_MPIdouble(
				    X->Def.EDGeneralTransfer[itrans][0],
				    X->Def.EDGeneralTransfer[itrans][1],
				    X->Def.EDGeneralTransfer[itrans][2],
				    X->Def.EDGeneralTransfer[itrans][3],
				    X->Def.EDParaGeneralTransfer[itrans],
				    X,
				    tmp_v0,
				    tmp_v1
				    );
  X->Large.prdct += dam_pr;
#endif
}/*void GC_child_general_hopp_MPIdouble*/

/**
 *
 * Hopping term in Hubbard + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
double complex X_GC_child_general_hopp_MPIdouble(
				       int org_isite1,
				       int org_ispin1,
				       int org_isite2,
				       int org_ispin2,
				       double complex tmp_trans,
				       struct BindStruct *X /**< [inout]*/,
				       double complex *tmp_v0 /**< [out] Result v0 = H v1*/, 
				       double complex *tmp_v1 /**< [in] v0 = H v1*/) {
#ifdef MPI
    int mask1, mask2, state1, state2, ierr, origin, bitdiff, Fsgn;
    unsigned long int idim_max_buf, j;
    MPI_Status statusMPI;
    double complex trans, dmv, dam_pr;

    mask1 = (int) X->Def.Tpow[2 * org_isite1 + org_ispin1];
    mask2 = (int) X->Def.Tpow[2 * org_isite2 + org_ispin2];
    if (mask2 > mask1) bitdiff = mask2 - mask1 * 2;
    else bitdiff = mask1 - mask2 * 2;
    origin = myrank ^ (mask1 + mask2);

    state1 = origin & mask1;
    state2 = origin & mask2;

    SgnBit((unsigned long int) (origin & bitdiff), &Fsgn); // Fermion sign

    if (state1 == 0 && state2 == mask2) {
        trans = -(double) Fsgn * tmp_trans;
    }
    else if (state1 == mask1 && state2 == 0) {
        trans = -(double) Fsgn * conj(tmp_trans);
        if (X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) {
            trans = 0;
        }
    }
    else return 0;

    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) {
        exitMPI(-1);
    }
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) {
        exitMPI(-1);
    }
    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) firstprivate(idim_max_buf, trans, X) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {
            dmv = trans * v1buf[j];
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) firstprivate(idim_max_buf, trans, X) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {
            dmv = trans * v1buf[j];
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    return (dam_pr);

#endif
}/*void GC_child_general_hopp_MPIdouble*/

/**
 *
 * Hopping term in Hubbard + Canonical
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
double complex X_child_CisAjt_MPIdouble(
				       int org_isite1,
				       int org_ispin1,
				       int org_isite2,
				       int org_ispin2,
				       double complex tmp_trans,
				       struct BindStruct *X /**< [inout]*/,
				       double complex *tmp_v0 /**< [out] Result v0 = H v1*/, 
				       double complex *tmp_v1 /**< [in] v0 = H v1*/,
                       long unsigned int *list_1_org,
                       long unsigned int *list_1buf_org,
                       long unsigned int *list_2_1_target,
                       long unsigned int *list_2_2_target
) {
#ifdef MPI
    int mask1, mask2, state1, state2, ierr, origin, bitdiff, Fsgn;
    unsigned long int idim_max_buf, j, ioff;
    MPI_Status statusMPI;
    double complex trans, dmv, dam_pr;

    mask1 = (int) X->Def.Tpow[2 * org_isite1 + org_ispin1];
    mask2 = (int) X->Def.Tpow[2 * org_isite2 + org_ispin2];
    if (mask2 > mask1) bitdiff = mask2 - mask1 * 2;
    else bitdiff = mask1 - mask2 * 2;
    origin = myrank ^ (mask1 + mask2);

    state1 = origin & mask1;
    state2 = origin & mask2;

    SgnBit((unsigned long int) (origin & bitdiff), &Fsgn); // Fermion sign

    if (state1 == 0 && state2 == mask2) {
        trans = -(double) Fsgn * tmp_trans;
    }
    else if (state1 == mask1 && state2 == 0) {
        trans = -(double) Fsgn * conj(tmp_trans);
        if (X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) {
            trans = 0;
        }
    }
    else return 0;

    ierr = MPI_Sendrecv(&X->Check.idim_maxOrg, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) {
        exitMPI(-1);
    }

    ierr = MPI_Sendrecv(list_1_org, X->Check.idim_maxOrg + 1, MPI_UNSIGNED_LONG, origin, 0,
                        list_1buf_org, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_maxOrg + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) {
        exitMPI(-1);
    }
    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) private(j, dmv, ioff)                 \
  firstprivate(idim_max_buf, trans, X, list_2_1_target, list_2_2_target, list_1buf_org) shared(v1buf, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++){
            dmv = trans * v1buf[j];
            GetOffComp(list_2_1_target, list_2_2_target, list_1buf_org[j],
                       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
            tmp_v0[ioff] += dmv;
        }
    }
    else {
      return 0.0;
    }
    return 1.0;
#endif
}/*void child_CisAjt_MPIdouble*/


/**
  *
  * Hopping term in Hubbard + GC
  * When only site2 is in the inter process region.
  *
  * @author Mitsuaki Kawamura (The University of Tokyo)
  */
void GC_child_general_hopp_MPIsingle(
  unsigned long int itrans /**< [in] Transfer ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  double complex dam_pr=0;
  dam_pr=X_GC_child_general_hopp_MPIsingle(
				    X->Def.EDGeneralTransfer[itrans][0],
				    X->Def.EDGeneralTransfer[itrans][1],
				    X->Def.EDGeneralTransfer[itrans][2],
				    X->Def.EDGeneralTransfer[itrans][3],
				    X->Def.EDParaGeneralTransfer[itrans],
				    X,
				    tmp_v0,
				    tmp_v1
				    );
  X->Large.prdct += dam_pr;
#endif
}/*void GC_child_general_hopp_MPIsingle*/

/**
  *
  * Hopping term in Hubbard + GC
  * When only site2 is in the inter process region.
  *
  * @author Mitsuaki Kawamura (The University of Tokyo)
  * @author Kazuyoshi Yoshimi (The University of Tokyo)
  */
double complex X_GC_child_general_hopp_MPIsingle(
				       int org_isite1,
				       int org_ispin1,
				       int org_isite2,
				       int org_ispin2,
				       double complex tmp_trans,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/) {
#ifdef MPI
    int mask2, state1, state2, ierr, origin, bit2diff, Fsgn;
    unsigned long int idim_max_buf, j, mask1, state1check, bit1diff, ioff;
    MPI_Status statusMPI;
    double complex trans, dmv, dam_pr;
    /*
     Prepare index in the inter PE
    */
    mask2 = (int) X->Def.Tpow[2 * org_isite2 + org_ispin2];
    bit2diff = mask2 - 1;
    origin = myrank ^ mask2;
    state2 = origin & mask2;

    SgnBit((unsigned long int) (origin & bit2diff), &Fsgn); // Fermion sign

    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    /*
     Index in the intra PE
    */
    mask1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];

    if (state2 == mask2) {
        trans = -(double) Fsgn * tmp_trans;
        state1check = 0;
    }
    else if (state2 == 0) {
        state1check = mask1;
        trans = -(double) Fsgn * conj(tmp_trans);
        if (X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) {
            trans = 0;
        }

    }
    else return 0;

    bit1diff = X->Def.Tpow[2 * X->Def.Nsite - 1] * 2 - mask1 * 2;

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, state1, Fsgn, ioff) \
  firstprivate(idim_max_buf, trans, X, mask1, state1check, bit1diff) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 0; j < idim_max_buf; j++) {

            state1 = j & mask1;

            if (state1 == state1check) {

                SgnBit(j & bit1diff, &Fsgn);
                ioff = j ^ mask1;

                dmv = (double) Fsgn * trans * v1buf[j + 1];
                tmp_v0[ioff + 1] += dmv;
                dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
            }
        }

    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, state1, Fsgn, ioff) \
  firstprivate(idim_max_buf, trans, X, mask1, state1check, bit1diff) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 0; j < idim_max_buf; j++) {

            state1 = j & mask1;

            if (state1 == state1check) {

                SgnBit(j & bit1diff, &Fsgn);
                ioff = j ^ mask1;

                dmv = (double) Fsgn * trans * v1buf[j + 1];
                dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
            }
        }
    }
    return (dam_pr);
#endif
}/*void GC_child_general_hopp_MPIsingle*/

/**
 *
 * Hopping term in Hubbard (Kondo) + Canonical ensemble
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void child_general_hopp_MPIdouble(
  unsigned long int itrans /**< [in] Transfer ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  double complex dam_pr;
  dam_pr =X_child_general_hopp_MPIdouble( X->Def.EDGeneralTransfer[itrans][0],
				  X->Def.EDGeneralTransfer[itrans][1],
				  X->Def.EDGeneralTransfer[itrans][2],
				  X->Def.EDGeneralTransfer[itrans][3],
				  X->Def.EDParaGeneralTransfer[itrans],
				  X,
				  tmp_v0,
				  tmp_v1);
  X->Large.prdct += dam_pr;

#endif
}/*void child_general_hopp_MPIdouble*/

/**
 *
 * Hopping term in Hubbard (Kondo) + Canonical ensemble
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
double complex X_child_general_hopp_MPIdouble(
				    int org_isite1,
				    int org_ispin1,
				    int org_isite2,
				    int org_ispin2,
				    double complex tmp_trans,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/) {
#ifdef MPI
    int mask1, mask2, state1, state2, ierr, origin, bitdiff, Fsgn;
    unsigned long int idim_max_buf, j, ioff;
    MPI_Status statusMPI;
    double complex trans, dmv, dam_pr;

    mask1 = (int) X->Def.Tpow[2 * org_isite1 + org_ispin1];
    mask2 = (int) X->Def.Tpow[2 * org_isite2 + org_ispin2];

    if (mask2 > mask1) bitdiff = mask2 - mask1 * 2;
    else bitdiff = mask1 - mask2 * 2;
    origin = myrank ^ (mask1 + mask2);

    state1 = origin & mask1;
    state2 = origin & mask2;

    SgnBit((unsigned long int) (origin & bitdiff), &Fsgn); // Fermion sign

    if (state1 == 0 && state2 == mask2) {
        trans = -(double) Fsgn * tmp_trans;
    }
    else if (state1 == mask1 && state2 == 0) {
        trans = -(double) Fsgn * conj(tmp_trans);
        if (X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) {
            trans = 0;
        }
    }
    else return 0;

    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                        list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff) \
  firstprivate(idim_max_buf, trans, X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {
            GetOffComp(list_2_1, list_2_2, list_1buf[j],
                       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
            dmv = trans * v1buf[j];
            tmp_v0[ioff] += dmv;
            dam_pr += conj(tmp_v1[ioff]) * dmv;
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff) \
  firstprivate(idim_max_buf, trans, X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {
            GetOffComp(list_2_1, list_2_2, list_1buf[j],
                       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
            dmv = trans * v1buf[j];
            dam_pr += conj(tmp_v1[ioff]) * dmv;
        }
    }
    return (dam_pr);

#endif
}/*void child_general_hopp_MPIdouble*/


/**
 *
 * Hopping term in Hubbard (Kondo) + Canonical ensemble
 * When only site2 is in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void child_general_hopp_MPIsingle(
  unsigned long int itrans /**< [in] Transfer ID*/, 
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  double complex dam_pr;
  dam_pr =X_child_general_hopp_MPIsingle( X->Def.EDGeneralTransfer[itrans][0],
				  X->Def.EDGeneralTransfer[itrans][1],
				  X->Def.EDGeneralTransfer[itrans][2],
				  X->Def.EDGeneralTransfer[itrans][3],
				  X->Def.EDParaGeneralTransfer[itrans],
				  X,
				  tmp_v0,
				  tmp_v1);
  X->Large.prdct += dam_pr;

#endif
}/*void child_general_hopp_MPIsingle*/

/**
 *
 * Hopping term in Hubbard (Kondo) + Canonical ensemble
 * When only site2 is in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
double complex X_child_general_hopp_MPIsingle(
				    int org_isite1,
				    int org_ispin1,
				    int org_isite2,
				    int org_ispin2,
				    double complex tmp_trans,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/) {
#ifdef MPI
  int mask2, state2, ierr, origin, bit2diff, Fsgn;
  unsigned long int mask1, state1, idim_max_buf, j, state1check, bit1diff, ioff, jreal;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[2 * org_isite2+org_ispin2];
  bit2diff = mask2 - 1;
  origin = myrank ^ mask2;

    state2 = origin & mask2;

    SgnBit((unsigned long int) (origin & bit2diff), &Fsgn); // Fermion sign

    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                        list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    /*
    Index in the intra PE
    */
    mask1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
    if (state2 == mask2) {
        trans = -(double) Fsgn * tmp_trans;
        state1check = 0;
    }
    else if (state2 == 0) {
        state1check = mask1;
        trans = -(double) Fsgn * conj(tmp_trans);
        if (X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) {
            trans = 0;
        }
    }
    else return 0;

    bit1diff = X->Def.Tpow[2 * X->Def.Nsite - 1] * 2 - mask1 * 2;

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff, jreal, state1) \
  firstprivate(idim_max_buf, trans, X, mask1, state1check, bit1diff) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {

            jreal = list_1buf[j];
            state1 = jreal & mask1;

            if (state1 == state1check) {

                SgnBit(jreal & bit1diff, &Fsgn);
                GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
                           X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);

                dmv = (double) Fsgn * trans * v1buf[j];
                tmp_v0[ioff] += dmv;
                dam_pr += conj(tmp_v1[ioff]) * dmv;
            }
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff, jreal, state1) \
  firstprivate(idim_max_buf, trans, X, mask1, state1check, bit1diff) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {

            jreal = list_1buf[j];
            state1 = jreal & mask1;

            if (state1 == state1check) {

                SgnBit(jreal & bit1diff, &Fsgn);
                GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
                           X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);

                dmv = (double) Fsgn * trans * v1buf[j];
                dam_pr += conj(tmp_v1[ioff]) * dmv;
            }
        }
    }
    return (dam_pr);
#endif
}/*double complex child_general_hopp_MPIsingle*/

/**
 *
 * Exchange term in Spin model
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void child_general_int_spin_MPIdouble(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  
  double complex dam_pr = 0;
  dam_pr=X_child_general_int_spin_MPIdouble
    (
     (int) X->Def.InterAll_OffDiagonal[i_int][0], (int)X->Def.InterAll_OffDiagonal[i_int][1], (int)X->Def.InterAll_OffDiagonal[i_int][3],
     (int) X->Def.InterAll_OffDiagonal[i_int][4], (int)X->Def.InterAll_OffDiagonal[i_int][5], (int)X->Def.InterAll_OffDiagonal[i_int][7],
     X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  
  X->Large.prdct += dam_pr;

#endif
}/*void child_general_int_spin_MPIdouble*/

/**
 *
 * Exchange term in Spin model
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
double complex X_child_general_int_spin_MPIdouble(
						  int org_isite1,
						  int org_ispin1,
						  int org_ispin2,
						  int org_isite3,
						  int org_ispin3,
						  int org_ispin4,
						  double complex tmp_J,
						  struct BindStruct *X,
						  double complex *tmp_v0,
						  double complex *tmp_v1
						  ) {
#ifdef MPI
    int mask1, mask2, state1, state2, ierr, origin;
    unsigned long int idim_max_buf, j, ioff;
    MPI_Status statusMPI;
    double complex Jint, dmv, dam_pr;

    mask1 = (int) X->Def.Tpow[org_isite1];
    mask2 = (int) X->Def.Tpow[org_isite3];
    origin = myrank ^ (mask1 + mask2);

    state1 = (origin & mask1) / mask1;
    state2 = (origin & mask2) / mask2;

    if (state1 == org_ispin2 && state2 == org_ispin4) {
        Jint = tmp_J;
    }
    else if (state1 == org_ispin1 && state2 == org_ispin3) {
        Jint = conj(tmp_J);
        if (X->Large.mode == M_CORR || X->Large.mode ==M_CALCSPEC) {
            Jint = 0;
        }
    }
    else return 0;

    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                        list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff) \
  firstprivate(idim_max_buf, Jint, X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {
            GetOffComp(list_2_1, list_2_2, list_1buf[j],
                       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
            dmv = Jint * v1buf[j];
            tmp_v0[ioff] += dmv;
            dam_pr += conj(tmp_v1[ioff]) * dmv;
        }

    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff) \
  firstprivate(idim_max_buf, Jint, X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {
            GetOffComp(list_2_1, list_2_2, list_1buf[j],
                       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
            dmv = Jint * v1buf[j];
            dam_pr += conj(tmp_v1[ioff]) * dmv;
        }
    }
    return dam_pr;
#endif
}/*double complex X_child_general_int_spin_MPIdouble*/


/**
 *
 * Exchange term in Spin model
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
double complex X_child_general_int_spin_TotalS_MPIdouble(
						  int org_isite1,
						  int org_isite3,
						  struct BindStruct *X,
						  double complex *tmp_v0,
						  double complex *tmp_v1
						  )
{
#ifdef MPI
  int mask1, mask2, num1_up, num2_up, ierr, origin;
  unsigned long int idim_max_buf, j, ioff, ibit_tmp;
  MPI_Status statusMPI;
  double complex dmv, dam_pr;

  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  if(mask1 == mask2){
     origin = myrank ^ mask1;
  }
  else{
    origin = myrank ^ (mask1 + mask2);
  }
  num1_up = (origin & mask1) / mask1;
  num2_up = (origin & mask2) / mask2;

  ibit_tmp=(num1_up)^(num2_up);
  if(ibit_tmp ==0) return 0;
  
  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
    list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
  #pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff) \
    firstprivate(idim_max_buf,  X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
  for (j = 1; j <= idim_max_buf; j++) {    
    GetOffComp(list_2_1, list_2_2, list_1buf[j],
	       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
    dmv = 0.5 * v1buf[j];
    dam_pr += conj(tmp_v1[ioff]) * dmv;
  }
  return dam_pr;  
#endif
}/*double complex X_child_general_int_spin_MPIdouble*/


/**
 *
 * Exchange term in Spin model
 * When only site2 is in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void child_general_int_spin_MPIsingle(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  double complex dam_pr = 0;
  
  dam_pr = X_child_general_int_spin_MPIsingle(
  (int)X->Def.InterAll_OffDiagonal[i_int][0], (int)X->Def.InterAll_OffDiagonal[i_int][1], (int)X->Def.InterAll_OffDiagonal[i_int][3],
  (int)X->Def.InterAll_OffDiagonal[i_int][4], (int)X->Def.InterAll_OffDiagonal[i_int][5], (int)X->Def.InterAll_OffDiagonal[i_int][7],
  X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  
  X->Large.prdct += dam_pr;
  
#endif
}/*void child_general_int_spin_MPIsingle*/


double complex X_child_general_int_spin_MPIsingle(
						  int org_isite1,
						  int org_ispin1,
						  int org_ispin2,
						  int org_isite3,
						  int org_ispin3,
						  int org_ispin4,
						  double complex tmp_J,
						  struct BindStruct *X,
						  double complex *tmp_v0,
						  double complex *tmp_v1
						  ) {
#ifdef MPI
    int mask2, state2, ierr, origin;
    unsigned long int mask1, idim_max_buf, j, ioff, state1, jreal, state1check;
    MPI_Status statusMPI;
    double complex Jint, dmv, dam_pr;
    /*
    Prepare index in the inter PE
    */
    mask2 = (int) X->Def.Tpow[org_isite3];
    origin = myrank ^ mask2;
    state2 = (origin & mask2) / mask2;

    if (state2 == org_ispin4) {
        state1check = (unsigned long int) org_ispin2;
        Jint = tmp_J;
    }
    else if (state2 == org_ispin3) {
        state1check = (unsigned long int) org_ispin1;
        Jint = conj(tmp_J);
        if (X->Large.mode == M_CORR || X->Large.mode ==M_CALCSPEC) {
            Jint = 0;
        }
    }
    else return 0;

    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                        list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    /*
    Index in the intra PE
    */
    mask1 = X->Def.Tpow[org_isite1];

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff, jreal, state1) \
  firstprivate(idim_max_buf, Jint, X, mask1, state1check, org_isite1) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {

            jreal = list_1buf[j];

            state1 = (jreal & mask1) / mask1;
            if (state1 == state1check) {
                GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
                           X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);

                dmv = Jint * v1buf[j];
                tmp_v0[ioff] += dmv;
                dam_pr += conj(tmp_v1[ioff]) * dmv;
            }
        }

    }
    else if (X->Large.mode == M_TOTALS) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff, jreal, state1) \
  firstprivate(idim_max_buf, Jint, X, mask1, state1check, org_isite1) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {

            jreal = list_1buf[j];

            state1 = (jreal & mask1) / mask1;
            if (state1 == state1check) {
                GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
                           X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);

                dmv = Jint * v1buf[j];
                dmv = 0.5 * v1buf[j];
                dam_pr += conj(tmp_v1[ioff]) * dmv;
            }
        }
    }
    else{
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff, jreal, state1) \
  firstprivate(idim_max_buf, Jint, X, mask1, state1check, org_isite1) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {

            jreal = list_1buf[j];

            state1 = (jreal & mask1) / mask1;
            if (state1 == state1check) {
                GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
                           X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
                dmv = Jint * v1buf[j];
                dam_pr += conj(tmp_v1[ioff]) * dmv;
            }
        }
    }
    return dam_pr;

#endif
}

/**
 *
 * Exchange and Pairlifting term in Spin model + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void GC_child_CisAitCiuAiv_spin_MPIdouble(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  double complex dam_pr;  
  dam_pr =  X_GC_child_CisAitCiuAiv_spin_MPIdouble( X->Def.InterAll_OffDiagonal[i_int][0],  X->Def.InterAll_OffDiagonal[i_int][1],  X->Def.InterAll_OffDiagonal[i_int][3],  X->Def.InterAll_OffDiagonal[i_int][4],  X->Def.InterAll_OffDiagonal[i_int][5],  X->Def.InterAll_OffDiagonal[i_int][7],X->Def.ParaInterAll_OffDiagonal[i_int],X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}/*void GC_child_CisAitCiuAiv_spin_MPIdouble*/


/**
 *
 * Exchange and Pairlifting term in Spin model + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
double complex X_GC_child_CisAitCiuAiv_spin_MPIdouble(
  int org_isite1, int org_ispin1, int org_ispin2,
  int org_isite3, int org_ispin3, int org_ispin4,
  double complex tmp_J, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1) {
#ifdef MPI
    int mask1, mask2, state1, state2, ierr, origin;
    unsigned long int idim_max_buf, j;
    MPI_Status statusMPI;
    double complex Jint, dmv, dam_pr;

    mask1 = (int) X->Def.Tpow[org_isite1];
    mask2 = (int) X->Def.Tpow[org_isite3];
    if (org_isite1 != org_isite3) {
        origin = myrank ^ (mask1 + mask2);
    }
    else {
        if (org_ispin1 == org_ispin4 && org_ispin2 == org_ispin3) { //CisAitCitAis=CisAis
            dam_pr = X_GC_child_CisAis_spin_MPIdouble(org_isite1, org_ispin1, tmp_J, X, tmp_v0, tmp_v1);
            return (dam_pr);
        }
        else { //CisAitCisAit=0
            return 0.0;
        }
    }


    state1 = (origin & mask1) / mask1;
    state2 = (origin & mask2) / mask2;

    if (state1 == org_ispin2 && state2 == org_ispin4) {
        Jint = tmp_J;
    }
    else if (state1 == org_ispin1 && state2 == org_ispin3) {
        Jint = conj(tmp_J);
        if (X->Large.mode == M_CORR || X->Large.mode ==M_CALCSPEC) {
            Jint = 0;
        }
    }
    else {
        return 0;
    }

    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {
            dmv = Jint * v1buf[j];
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {
            dmv = Jint * v1buf[j];
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    return dam_pr;

#endif
}/*void GC_child_CisAitCiuAiv_spin_MPIdouble*/

/**
 *
 * Wrapper for calculating CisAisCjuAjv term in Spin model + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void GC_child_CisAisCjuAjv_spin_MPIdouble(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI

  double complex dam_pr;
  dam_pr = X_GC_child_CisAisCjuAjv_spin_MPIdouble( X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],
						   X->Def.InterAll_OffDiagonal[i_int][4], X->Def.InterAll_OffDiagonal[i_int][5], X->Def.InterAll_OffDiagonal[i_int][7],
						   X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0,  tmp_v1);
  X->Large.prdct += dam_pr;

#endif
}/*void GC_child_CisAitCiuAiv_spin_MPIdouble*/


/**
 *
 * CisAisCjuAjv term in Spin model + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_GC_child_CisAisCjuAjv_spin_MPIdouble(
					    int org_isite1,
					    int org_ispin1,
					    int org_isite3,
					    int org_ispin3,
					    int org_ispin4,
					    double complex tmp_J,
					    struct BindStruct *X,
					    double complex *tmp_v0,
					    double complex *tmp_v1) {
#ifdef MPI
    int mask1, mask2, state2, ierr;
    long int origin, num1;
    unsigned long int idim_max_buf, j;
    MPI_Status statusMPI;
    double complex Jint, dmv, dam_pr;

    if (org_isite1 == org_isite3 && org_ispin1 == org_ispin4) {//CisAisCitAis
        return 0.0;
    }

    mask1 = (int) X->Def.Tpow[org_isite1];
    mask2 = (int) X->Def.Tpow[org_isite3];
    origin = myrank ^ mask2;
    state2 = (origin & mask2) / mask2;
    num1 = X_SpinGC_CisAis((unsigned long int) myrank + 1, X, mask1, org_ispin1);
    if (num1 != 0 && state2 == org_ispin4) {
        Jint = tmp_J;
    }
    else if (X_SpinGC_CisAis(origin + 1, X, mask1, org_ispin1) == TRUE && state2 == org_ispin3) {
        Jint = conj(tmp_J);
        if (X->Large.mode == M_CORR || X->Large.mode ==M_CALCSPEC) Jint = 0;
    }
    else {
        return 0.0;
    }

    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {
            dmv = Jint * v1buf[j];
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {
            dmv = Jint * v1buf[j];
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    return (dam_pr);
#endif
}/*double complex X_GC_child_CisAisCjuAjv_spin_MPIdouble*/

/**
 *
 * Wrapper for calculating CisAitCjuAju term in Spin model + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void GC_child_CisAitCjuAju_spin_MPIdouble(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI

  double complex dam_pr;
  dam_pr = X_GC_child_CisAitCjuAju_spin_MPIdouble( X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],X->Def.InterAll_OffDiagonal[i_int][3],
						   X->Def.InterAll_OffDiagonal[i_int][4], X->Def.InterAll_OffDiagonal[i_int][5],
						   X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0,  tmp_v1);
  X->Large.prdct += dam_pr;

#endif
}/*void GC_child_CisAitCiuAiv_spin_MPIdouble*/


/**
 *
 * CisAisCjuAjv term in Spin model + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_GC_child_CisAitCjuAju_spin_MPIdouble(
					    int org_isite1,
					    int org_ispin1,
					    int org_ispin2,
					    int org_isite3,
					    int org_ispin3,
					    double complex tmp_J,
					    struct BindStruct *X,
					    double complex *tmp_v0,
					    double complex *tmp_v1) {
#ifdef MPI
    int mask1, mask2, state1, ierr, num1;
    long int origin;
    unsigned long int idim_max_buf, j;
    MPI_Status statusMPI;
    double complex Jint, dmv, dam_pr;

    if (org_isite1 == org_isite3 && org_ispin1 == org_ispin3) {//cisaitcisais
        return 0.0;
    }

    mask1 = (int) X->Def.Tpow[org_isite1];
    origin = myrank ^ mask1;
    state1 = (origin & mask1) / mask1;
    mask2 = (int) X->Def.Tpow[org_isite3];
    num1 = X_SpinGC_CisAis(origin + 1, X, mask2, org_ispin3);
    if (state1 == org_ispin2 && num1 != 0) {
        Jint = tmp_J;
    }
    else if (state1 == org_ispin1) {
        num1 = X_SpinGC_CisAis((unsigned long int) myrank + 1, X, mask2, org_ispin3);
        if (num1 != 0) {
            Jint = conj(tmp_J);
            if (X->Large.mode == M_CORR || X->Large.mode ==M_CALCSPEC) {
                Jint = 0;
            }
        }
        else {
            return 0.0;
        }
    }

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {
            dmv = Jint * v1buf[j];
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= idim_max_buf; j++) {
            dmv = Jint * v1buf[j];
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    return (dam_pr);
#endif
}/*double complex X_GC_child_CisAisCjuAjv_spin_MPIdouble*/

/**
 *
 * CisAisCjuAjv term in Spin model + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_GC_child_CisAisCjuAju_spin_MPIdouble(
					    int org_isite1,
					    int org_ispin1,
					    int org_isite3,
					    int org_ispin3,
					    double complex tmp_J,
					    struct BindStruct *X,
					    double complex *tmp_v0,
					    double complex *tmp_v1)
{
#ifdef MPI
  long unsigned int mask1, mask2, num1,num2;
  unsigned long int  j;
//  MPI_Status statusMPI;
  double complex dmv, dam_pr;
  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  num1 =  X_SpinGC_CisAis((unsigned long int)myrank + 1, X, mask1, org_ispin1);
  num2 = X_SpinGC_CisAis((unsigned long int)myrank + 1, X, mask2, org_ispin3);
  
  dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
    firstprivate(tmp_J, X, num1, num2) shared(tmp_v1, tmp_v0)
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = num1*num2*tmp_v1[j]*tmp_J;
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
    firstprivate(tmp_J, X, num1, num2) shared(tmp_v1, tmp_v0)
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = num1 * num2 * tmp_v1[j] * tmp_J;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
  return(dam_pr);
#endif
}/*double complex X_GC_child_CisAisCjuAju_spin_MPIdouble*/

/**
 *
 * CisAisCjuAjv term in Spin model + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_GC_child_CisAisCjuAju_spin_MPIsingle(
					    int org_isite1,
					    int org_ispin1,
					    int org_isite3,
					    int org_ispin3,
					    double complex tmp_J,
					    struct BindStruct *X,
					    double complex *tmp_v0,
					    double complex *tmp_v1) {
#ifdef MPI
    long unsigned int mask1, mask2, num1, num2;
    unsigned long int j;
//  MPI_Status statusMPI;
    double complex Jint, dmv, dam_pr;
    Jint = tmp_J;
    mask1 = (int) X->Def.Tpow[org_isite1];
    mask2 = (int) X->Def.Tpow[org_isite3];
    num2 = X_SpinGC_CisAis((unsigned long int) myrank + 1, X, mask2, org_ispin3);

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, num1) \
  firstprivate(Jint, X, num2, mask1, org_ispin1) shared(tmp_v1, tmp_v0)
        for (j = 1; j <= X->Check.idim_max; j++) {
            num1 = X_SpinGC_CisAis(j, X, mask1, org_ispin1);
            dmv = Jint * num1 * num2 * tmp_v1[j];
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }

    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, num1) \
  firstprivate(Jint, X, num2, mask1, org_ispin1) shared(tmp_v1, tmp_v0)
        for (j = 1; j <= X->Check.idim_max; j++) {
            num1 = X_SpinGC_CisAis(j, X, mask1, org_ispin1);
            dmv = Jint * num1 * num2 * tmp_v1[j];
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    return (dam_pr);
#endif
}/*double complex X_GC_child_CisAisCjuAju_spin_MPIdouble*/


/**
 *
 * General interaction term in the Spin model + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void GC_child_general_int_spin_MPIdouble(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
      X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
    GC_child_CisAisCjuAjv_spin_MPIdouble(i_int, X, tmp_v0, tmp_v1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
    GC_child_CisAitCjuAju_spin_MPIdouble(i_int, X, tmp_v0, tmp_v1);
  }
  else {
      GC_child_CisAitCiuAiv_spin_MPIdouble(i_int, X, tmp_v0, tmp_v1);
  }
}/*void GC_child_general_int_spin_MPIdouble*/


/**
 *
 * Exchange and Pairlifting term in Spin model + GC
 * When only site2 is in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void GC_child_CisAitCiuAiv_spin_MPIsingle(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  double complex dam_pr;  
  dam_pr =X_GC_child_CisAitCiuAiv_spin_MPIsingle(X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1], X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4], X->Def.InterAll_OffDiagonal[i_int][5], X->Def.InterAll_OffDiagonal[i_int][7], X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;

#endif
}/*void GC_child_CisAitCiuAiv_spin_MPIsingle*/

/**
 *
 * Exchange and Pairlifting term in Spin model + GC
 * When only site2 is in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
double complex X_GC_child_CisAitCiuAiv_spin_MPIsingle(
					    int org_isite1, int org_ispin1, int org_ispin2,
					    int org_isite3, int org_ispin3, int org_ispin4,
					    double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1) {
#ifdef MPI
    int mask2, state2, ierr, origin;
    unsigned long int mask1, idim_max_buf, j, ioff, state1, state1check;
    MPI_Status statusMPI;
    double complex Jint, dmv, dam_pr;
    /*
    Prepare index in the inter PE
    */
    mask2 = (int) X->Def.Tpow[org_isite3];
    origin = myrank ^ mask2;
    state2 = (origin & mask2) / mask2;

    if (state2 == org_ispin4) {
        state1check = (unsigned long int) org_ispin2;
        Jint = tmp_J;
    }
    else if (state2 == org_ispin3) {
        state1check = (unsigned long int) org_ispin1;
        Jint = conj(tmp_J);
        if (X->Large.mode == M_CORR || X->Large.mode ==M_CALCSPEC) {
            Jint = 0;
        }
    }
    else return 0.0;

    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    /*
    Index in the intra PE
    */
    mask1 = X->Def.Tpow[org_isite1];

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, state1, ioff) \
    firstprivate(idim_max_buf, Jint, X, state1check, mask1) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 0; j < idim_max_buf; j++) {
            state1 = X_SpinGC_CisAit(j + 1, X, mask1, state1check, &ioff);
            if (state1 != 0) {
                dmv = Jint * v1buf[j + 1];
                tmp_v0[ioff + 1] += dmv;
                dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
            }
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, state1, ioff) \
    firstprivate(idim_max_buf, Jint, X, state1check, mask1) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 0; j < idim_max_buf; j++) {
            state1 = X_SpinGC_CisAit(j + 1, X, mask1, state1check, &ioff);
            if (state1 != 0) {
                dmv = Jint * v1buf[j + 1];
                dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
            }
        }
    }

    return (dam_pr);

#endif
}/*void GC_child_CisAitCiuAiv_spin_MPIsingle*/


/**
 *
 * Wrapper for CisAisCjuAjv term in Spin model + GC
 * When only site2 is in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void GC_child_CisAisCjuAjv_spin_MPIsingle(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  double complex dam_pr;  
  dam_pr =X_GC_child_CisAisCjuAjv_spin_MPIsingle(X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1], X->Def.InterAll_OffDiagonal[i_int][4], X->Def.InterAll_OffDiagonal[i_int][5], X->Def.InterAll_OffDiagonal[i_int][7], X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;

#endif
}/*void GC_child_CisAisCjuAjv_spin_MPIsingle*/

/**
 *
 * CisAisCjuAjv term in Spin model + GC
 * When only site2 is in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_GC_child_CisAisCjuAjv_spin_MPIsingle( int org_isite1, int org_ispin1,  int org_isite3, int org_ispin3, int org_ispin4, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1) {
#ifdef MPI
    int mask2, state2, ierr, origin;
    unsigned long int mask1, idim_max_buf, j, state1, state1check;
    MPI_Status statusMPI;
    double complex Jint, dmv, dam_pr;
    /*
    Prepare index in the inter PE
    */
    mask2 = (int) X->Def.Tpow[org_isite3];
    origin = myrank ^ mask2;
    state2 = (origin & mask2) / mask2;
    if (state2 == org_ispin4) {
        state1check = (unsigned long int) org_ispin1;
        Jint = tmp_J;
    }
    else if (state2 == org_ispin3) {
        state1check = (unsigned long int) org_ispin1;
        Jint = conj(tmp_J);
        if (X->Large.mode == M_CORR || X->Large.mode ==M_CALCSPEC) {
            Jint = 0;
        }
    }
    else return 0.0;

    ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    /*
    Index in the intra PE
    */
    mask1 = X->Def.Tpow[org_isite1];

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, state1) \
    firstprivate(idim_max_buf, Jint, X, state1check, mask1) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 0; j < idim_max_buf; j++) {
            state1 = (j & mask1) / mask1;
            if (state1 == state1check) {
                dmv = Jint * v1buf[j + 1];
                tmp_v0[j + 1] += dmv;
                dam_pr += conj(tmp_v1[j + 1]) * dmv;
            }
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, state1) \
    firstprivate(idim_max_buf, Jint, X, state1check, mask1) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 0; j < idim_max_buf; j++) {
            state1 = (j & mask1) / mask1;
            if (state1 == state1check) {
                dmv = Jint * v1buf[j + 1];
                dam_pr += conj(tmp_v1[j + 1]) * dmv;
            }
        }
    }
    return (dam_pr);

#endif
}/*void GC_child_CisAitCiuAiv_spin_MPIsingle*/

/**
 *
 * Wrapper for CisAisCjuAjv term in Spin model + GC
 * When only site2 is in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void GC_child_CisAitCjuAju_spin_MPIsingle(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  double complex dam_pr;  
  dam_pr =X_GC_child_CisAitCjuAju_spin_MPIsingle(X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1], X->Def.InterAll_OffDiagonal[i_int][3], X->Def.InterAll_OffDiagonal[i_int][4], X->Def.InterAll_OffDiagonal[i_int][5], X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;

#endif
}/*void GC_child_CisAisCjuAjv_spin_MPIsingle*/

/**
 *
 * CisAisCjuAjv term in Spin model + GC
 * When only site2 is in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_GC_child_CisAitCjuAju_spin_MPIsingle( int org_isite1, int org_ispin1, int org_ispin2,  int org_isite3, int org_ispin3, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1) {
#ifdef MPI
    int mask2, state2;
    unsigned long int mask1, j, ioff, state1, state1check;
    //MPI_Status statusMPI;
    double complex Jint, dmv, dam_pr;
    /*
    Prepare index in the inter PE
    */
    mask2 = (int) X->Def.Tpow[org_isite3];
    state2 = (myrank & mask2) / mask2;

    if (state2 == org_ispin3) {
        state1check = org_ispin2;
        Jint = tmp_J;
    }
    else {
        return 0.0;
    }

    mask1 = (int) X->Def.Tpow[org_isite1];

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, state1, ioff) \
  firstprivate(Jint, X, state1check, mask1) shared( tmp_v1, tmp_v0)
        for (j = 0; j < X->Check.idim_max; j++) {

            state1 = (j & mask1) / mask1;
            ioff = j ^ mask1;
            if (state1 == state1check) {
                dmv = Jint * tmp_v1[j + 1];
            }
            else {
                dmv = conj(Jint) * tmp_v1[j + 1];
            }
            tmp_v0[ioff + 1] += dmv;
            dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
        }
    }
    else if (X->Large.mode == M_CORR ) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, state1, ioff) \
  firstprivate(Jint, X, state1check, mask1) shared( tmp_v1, tmp_v0)
        for (j = 0; j < X->Check.idim_max; j++) {

            state1 = (j & mask1) / mask1;
            ioff = j ^ mask1;
            if (state1 == state1check) {
                dmv = Jint * tmp_v1[j + 1];
            }
            else {
                dmv = 0.0;
            }
            dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
        }

    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, state1, ioff) \
  firstprivate(Jint, X, state1check, mask1) shared( tmp_v1, tmp_v0)
        for (j = 0; j < X->Check.idim_max; j++) {

            state1 = (j & mask1) / mask1;
            ioff = j ^ mask1;
            if (state1 == state1check) {
                dmv = Jint * tmp_v1[j + 1];
            }
            else {
                dmv = conj(Jint) * tmp_v1[j + 1];
            }
            dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
        }
    }
    return (dam_pr);

#endif
}/*void GC_child_CisAitCiuAiv_spin_MPIsingle*/


/**
 *
 * General interaction term in the Spin model + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void GC_child_general_int_spin_MPIsingle(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
      X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
        GC_child_CisAisCjuAjv_spin_MPIsingle(i_int, X, tmp_v0, tmp_v1);

  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
     GC_child_CisAitCjuAju_spin_MPIsingle(i_int, X, tmp_v0, tmp_v1);  
  }
  else {
    GC_child_CisAitCiuAiv_spin_MPIsingle(i_int, X, tmp_v0, tmp_v1);
  }
}/*void GC_child_general_int_spin_MPIsingle*/

/**
 *
 * General interaction term in the Spin model + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void GC_child_general_int_GeneralSpin_MPIdouble(
    unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
    double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
    double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  double complex dam_pr;
 // MPI_Status statusMPI;

  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
     dam_pr = X_GC_child_CisAisCjuAjv_GeneralSpin_MPIdouble(X->Def.InterAll_OffDiagonal[i_int][0],
							X->Def.InterAll_OffDiagonal[i_int][1],
							X->Def.InterAll_OffDiagonal[i_int][4],
							X->Def.InterAll_OffDiagonal[i_int][5],X->Def.InterAll_OffDiagonal[i_int][7],
							X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
    dam_pr = X_GC_child_CisAitCjuAju_GeneralSpin_MPIsingle(X->Def.InterAll_OffDiagonal[i_int][0],
							X->Def.InterAll_OffDiagonal[i_int][1],X->Def.InterAll_OffDiagonal[i_int][3],
							X->Def.InterAll_OffDiagonal[i_int][4],
							X->Def.InterAll_OffDiagonal[i_int][5],
							X->Def.ParaInterAll_OffDiagonal[i_int],X, tmp_v0, tmp_v1);
  }
  else {
    dam_pr = X_GC_child_CisAitCjuAjv_GeneralSpin_MPIdouble(X->Def.InterAll_OffDiagonal[i_int][0],
							X->Def.InterAll_OffDiagonal[i_int][1],X->Def.InterAll_OffDiagonal[i_int][3],
							X->Def.InterAll_OffDiagonal[i_int][4],
							X->Def.InterAll_OffDiagonal[i_int][5],X->Def.InterAll_OffDiagonal[i_int][7],
							X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  }
  X->Large.prdct += dam_pr;
#endif
}/*void GC_child_general_int_spin_MPIdouble*/


double complex X_GC_child_CisAisCjuAjv_GeneralSpin_MPIdouble(
							  int org_isite1,
							  int org_ispin1,
							  int org_isite3,
							  int org_ispin3,
							  int org_ispin4,
							  double complex tmp_J,
							  struct BindStruct *X,
							  double complex *tmp_v0,
							  double complex *tmp_v1
							  ) {
#ifdef MPI
    unsigned long int off, j;
    int origin, ierr;
    double complex tmp_V, dmv, dam_pr;
    MPI_Status statusMPI;
    int ihermite = TRUE;
    if (org_isite1 == org_isite3 && org_ispin1 == org_ispin4) {//cisaisciuais=0 && cisaiucisais=0
        return 0.0;
    }

    if (GetOffCompGeneralSpin((unsigned long int) myrank, org_isite3 + 1, org_ispin4, org_ispin3,
                              &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
        if (BitCheckGeneral(off, org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
            tmp_V = tmp_J;
        }
        else {
            ihermite = FALSE;
        }
    }
    else {
        ihermite = FALSE;
    }

    if (ihermite == FALSE) {
        if (BitCheckGeneral((unsigned long int) myrank, org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow) ==
            TRUE &&
            GetOffCompGeneralSpin((unsigned long int) myrank, org_isite3 + 1, org_ispin3, org_ispin4, &off,
                                  X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
            tmp_V = conj(tmp_J);
            if (X->Large.mode == M_CORR || X->Large.mode ==M_CALCSPEC) tmp_V = 0.0;
        }
        else {
            return 0.0;
        }
    }

    origin = (int) off;
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1, v1buf)
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = v1buf[j] * tmp_V;
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1, v1buf) 
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = v1buf[j] * tmp_V;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    return dam_pr;
#endif
}

double complex X_GC_child_CisAitCjuAju_GeneralSpin_MPIdouble(
							  int org_isite1,
							  int org_ispin1,
							  int org_ispin2,
							  int org_isite3,
							  int org_ispin3,
							  double complex tmp_J,
							  struct BindStruct *X,
							  double complex *tmp_v0,
							  double complex *tmp_v1
							  ) {
#ifdef MPI
    unsigned long int j, off;
    int origin, ierr;
    double complex tmp_V, dmv, dam_pr;
    MPI_Status statusMPI;

    if (org_isite1 == org_isite3 && org_ispin1 == org_ispin3) {//cisaitcisais=0 && cisaiscitais=0
        return 0.0;
    }

    if (BitCheckGeneral((unsigned long int) myrank, org_isite3 + 1, org_ispin3, X->Def.SiteToBit, X->Def.Tpow) == TRUE
        && GetOffCompGeneralSpin((unsigned long int) myrank, org_isite1 + 1, org_ispin2, org_ispin1,
                                 &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
        tmp_V = tmp_J;
    }
    else if (GetOffCompGeneralSpin((unsigned long int) myrank, org_isite1 + 1, org_ispin1, org_ispin2, &off,
                                   X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
        if (BitCheckGeneral(off, org_isite3 + 1, org_ispin3, X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
            tmp_V = conj(tmp_J);
            if (X->Large.mode == M_CORR || X->Large.mode ==M_CALCSPEC) tmp_V = 0.0;
        }
        else return 0.0;
    }
    else return 0.0;

    origin = (int) off;

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel  for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1, v1buf)
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = v1buf[j] * tmp_V;
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }

    }
    else {
#pragma omp parallel  for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1, v1buf) 
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = v1buf[j] * tmp_V;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    return dam_pr;
#endif
}

double complex X_GC_child_CisAitCjuAjv_GeneralSpin_MPIdouble(
							  int org_isite1,
							  int org_ispin1,
							  int org_ispin2,
							  int org_isite3,
							  int org_ispin3,
							  int org_ispin4,
							  double complex tmp_J,
							  struct BindStruct *X,
							  double complex *tmp_v0,
							  double complex *tmp_v1
							  ) {
#ifdef MPI
    unsigned long int tmp_off, off, j;
    int origin, ierr, ihermite;
    double complex tmp_V, dmv, dam_pr;
    MPI_Status statusMPI;

    ihermite = TRUE;

    if (org_isite1 == org_isite3 && org_ispin1 == org_ispin4 &&
        org_ispin2 == org_ispin3) { //cisaitcitais=cisais && cisaitcitais =cisais
        dam_pr = X_GC_child_CisAis_GeneralSpin_MPIdouble(org_isite1, org_ispin1, tmp_J, X, tmp_v0, tmp_v1);
        return (dam_pr);
    }
    //cisaitcisait
    if (GetOffCompGeneralSpin((unsigned long int) myrank, org_isite1 + 1, org_ispin1, org_ispin2,
                              &tmp_off, X->Def.SiteToBit, X->Def.Tpow) == TRUE) {

        if (GetOffCompGeneralSpin(tmp_off, org_isite3 + 1, org_ispin3, org_ispin4,
                                  &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE) {

            tmp_V = tmp_J;
        }
        else ihermite = FALSE;
    }
    else {
        ihermite = FALSE;
    }

    if (ihermite == FALSE) {
        if (GetOffCompGeneralSpin((unsigned long int) myrank, org_isite3 + 1, org_ispin4, org_ispin3, &tmp_off,
                                  X->Def.SiteToBit, X->Def.Tpow) == TRUE) {

            if (GetOffCompGeneralSpin(tmp_off, org_isite1 + 1, org_ispin2, org_ispin1, &off, X->Def.SiteToBit,
                                      X->Def.Tpow) == TRUE) {
                tmp_V = conj(tmp_J);
                if (X->Large.mode == M_CORR || X->Large.mode ==M_CALCSPEC) tmp_V = 0.0;
            }
            else return 0.0;
        }
        else return 0.0;
    }

    origin = (int) off;

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1, v1buf)
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = v1buf[j] * tmp_V;
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }

    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1, v1buf)
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = v1buf[j] * tmp_V;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    return dam_pr;
#endif
}

double complex X_GC_child_CisAisCjuAju_GeneralSpin_MPIdouble(
							  int org_isite1,
							  int org_ispin1,
							  int org_isite3,
							  int org_ispin3,
							  double complex tmp_J,
							  struct BindStruct *X,
							  double complex *tmp_v0,
							  double complex *tmp_v1
							  ) {
#ifdef MPI
    unsigned long int j, num1;
    double complex tmp_V, dmv, dam_pr;
    //MPI_Status statusMPI;

    num1 = BitCheckGeneral((unsigned long int) myrank, org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);

    if (num1 == TRUE) {
        num1 = BitCheckGeneral((unsigned long int) myrank, org_isite3 + 1, org_ispin3, X->Def.SiteToBit, X->Def.Tpow);
        if (num1 == TRUE) {
            tmp_V = tmp_J;
        }
        else return 0.0;
    }
    else return 0.0;

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1)
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = tmp_v1[j] * tmp_V;
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }

    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1) 
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = tmp_v1[j] * tmp_V;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    return dam_pr;
#endif
}

double complex X_GC_child_CisAit_GeneralSpin_MPIdouble(
						       int org_isite1,
						       int org_ispin1,
						       int org_ispin2,
						       double complex tmp_trans,
						       struct BindStruct *X,
						       double complex *tmp_v0,
						       double complex *tmp_v1) {
#ifdef MPI
    unsigned long int off, j;
    int origin, ierr;
    double complex tmp_V, dmv, dam_pr;
    MPI_Status statusMPI;

    if (GetOffCompGeneralSpin((unsigned long int) myrank, org_isite1 + 1, org_ispin1, org_ispin2,
                              &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
        tmp_V = tmp_trans;
    }
    else if (GetOffCompGeneralSpin((unsigned long int) myrank,
                                   org_isite1 + 1, org_ispin2, org_ispin1, &off,
                                   X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
        tmp_V = conj(tmp_trans);
        if (X->Large.mode == M_CORR || X->Large.mode ==M_CALCSPEC) tmp_V = 0.0;
    }
    else return 0.0;

    origin = (int) off;

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode ==M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1, v1buf)
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = v1buf[j] * tmp_V;
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }

    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1, v1buf) 
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = v1buf[j] * tmp_V;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    return dam_pr;
#endif
}

double complex X_GC_child_CisAis_GeneralSpin_MPIdouble(
						       int org_isite1,
						       int org_ispin1,
						       double complex tmp_trans,
						       struct BindStruct *X,
						       double complex *tmp_v0,
						       double complex *tmp_v1) {
#ifdef MPI
    unsigned long int j, num1;
    double complex tmp_V, dmv, dam_pr;
    //MPI_Status statusMPI;

    num1 = BitCheckGeneral((unsigned long int) myrank,
                           org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);
    if (num1 != 0) {
        tmp_V = tmp_trans;
    }
    else return 0.0;

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1)
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = tmp_v1[j] * tmp_V;
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1) 
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = tmp_v1[j] * tmp_V;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    return dam_pr;
#endif
}

double complex X_GC_child_AisCis_GeneralSpin_MPIdouble(
        int org_isite1,
        int org_ispin1,
        double complex tmp_trans,
        struct BindStruct *X,
        double complex *tmp_v0,
        double complex *tmp_v1) {
#ifdef MPI
    unsigned long int j, num1;
    double complex tmp_V, dmv, dam_pr;
    //MPI_Status statusMPI;

    num1 = BitCheckGeneral((unsigned long int) myrank,
                           org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);
    if (num1 == 0) {
        tmp_V = tmp_trans;
    }
    else return 0.0;

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1)
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = tmp_v1[j] * tmp_V;
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1)
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = tmp_v1[j] * tmp_V;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    return dam_pr;
#endif
}

 /**
 *
 * General interaction term in the Spin model + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void GC_child_general_int_GeneralSpin_MPIsingle(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  double complex dam_pr;
  //MPI_Status statusMPI;

  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
    dam_pr = X_GC_child_CisAisCjuAjv_GeneralSpin_MPIsingle(X->Def.InterAll_OffDiagonal[i_int][0],
							X->Def.InterAll_OffDiagonal[i_int][1],
							X->Def.InterAll_OffDiagonal[i_int][4],
							X->Def.InterAll_OffDiagonal[i_int][5],X->Def.InterAll_OffDiagonal[i_int][7],
							X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
    dam_pr = X_GC_child_CisAitCjuAju_GeneralSpin_MPIsingle(X->Def.InterAll_OffDiagonal[i_int][0],
							X->Def.InterAll_OffDiagonal[i_int][1],X->Def.InterAll_OffDiagonal[i_int][3],
							X->Def.InterAll_OffDiagonal[i_int][4],
							X->Def.InterAll_OffDiagonal[i_int][5],
							X->Def.ParaInterAll_OffDiagonal[i_int],X, tmp_v0, tmp_v1);
  }
  else {
    dam_pr = X_GC_child_CisAitCjuAjv_GeneralSpin_MPIsingle(X->Def.InterAll_OffDiagonal[i_int][0],
							X->Def.InterAll_OffDiagonal[i_int][1],X->Def.InterAll_OffDiagonal[i_int][3],
							X->Def.InterAll_OffDiagonal[i_int][4],
							X->Def.InterAll_OffDiagonal[i_int][5],X->Def.InterAll_OffDiagonal[i_int][7],
							X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  }

  X->Large.prdct += dam_pr;
#endif
}/*void GC_child_general_int_spin_MPIsingle*/
 
double complex X_GC_child_CisAisCjuAjv_GeneralSpin_MPIsingle(
							  int org_isite1,
							  int org_ispin1,
							  int org_isite3,
							  int org_ispin3,
							  int org_ispin4,
							  double complex tmp_J,
							  struct BindStruct *X,
							  double complex *tmp_v0,
							  double complex *tmp_v1
							  )
{
#ifdef MPI
  unsigned long int off, j, num1;
  int origin, ierr, isite, IniSpin;
  double complex tmp_V, dmv, dam_pr;
  MPI_Status statusMPI;


  if (GetOffCompGeneralSpin((unsigned long int)myrank,
			    org_isite3+ 1, org_ispin3, org_ispin4, &off,
			    X->Def.SiteToBit, X->Def.Tpow) == TRUE)
    {
      tmp_V = tmp_J;
      isite = org_isite1 + 1;
      IniSpin = org_ispin1;
    }
  else if (GetOffCompGeneralSpin((unsigned long int)myrank,
				 org_isite3+ 1, org_ispin4, org_ispin3, &off,
				 X->Def.SiteToBit, X->Def.Tpow) == TRUE)
    {
      tmp_V = conj(tmp_J);
      if(X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) tmp_V=0.0;
      isite = org_isite1 + 1;
      IniSpin = org_ispin1;
    }
  else return 0.0;
  
  origin = (int)off;
  
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
		      v1buf, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
		      MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC){
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, isite, IniSpin) private(j, dmv, num1) shared (tmp_v0, tmp_v1, v1buf)
        for (j = 1; j <= X->Check.idim_max; j++) {
            num1 = BitCheckGeneral(j-1, isite, IniSpin, X->Def.SiteToBit, X->Def.Tpow);
            if (num1 !=0)
            {
                dmv = v1buf[j] * tmp_V;
                tmp_v0[j] += dmv;
                dam_pr += conj(tmp_v1[j]) * dmv;
            }
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, isite, IniSpin) private(j, dmv, num1) shared (tmp_v0, tmp_v1, v1buf) 
        for (j = 1; j <= X->Check.idim_max; j++) {
            num1 = BitCheckGeneral(j - 1, isite, IniSpin, X->Def.SiteToBit, X->Def.Tpow);
            if (num1 != 0) {
                dmv = v1buf[j] * tmp_V;
                dam_pr += conj(tmp_v1[j]) * dmv;
            }
        }
    }
  return dam_pr;
#endif
}

double complex X_GC_child_CisAitCjuAju_GeneralSpin_MPIsingle(
							  int org_isite1,
							  int org_ispin1,
							  int org_ispin2,
							  int org_isite3,
							  int org_ispin3,
							  double complex tmp_J,
							  struct BindStruct *X,
							  double complex *tmp_v0,
							  double complex *tmp_v1
							  )
{
#ifdef MPI
  unsigned long int num1, j, off;
  int isite, IniSpin, FinSpin;
  double complex tmp_V, dmv, dam_pr;
  //MPI_Status statusMPI;

  num1 = BitCheckGeneral((unsigned long int)myrank, 
			 org_isite3+1, org_ispin3, X->Def.SiteToBit, X->Def.Tpow);
  if(num1 != 0){
    tmp_V = tmp_J;
    isite = org_isite1 + 1;
    IniSpin = org_ispin2;
    FinSpin = org_ispin1;
  }
  else return 0.0;

  dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC){
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, isite, IniSpin, FinSpin) private(j, dmv, num1, off) shared (tmp_v0, tmp_v1, v1buf)
        for (j = 1; j <= X->Check.idim_max; j++) {
            if (GetOffCompGeneralSpin(j - 1, isite, IniSpin, FinSpin, &off,
                                      X->Def.SiteToBit, X->Def.Tpow) == TRUE)
            {
                dmv = tmp_v1[j] * tmp_V;
                tmp_v0[off + 1] += dmv;
                dam_pr += conj(tmp_v1[off + 1]) * dmv;
            }
        }

    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, isite, IniSpin, FinSpin) private(j, dmv, num1, off) shared (tmp_v0, tmp_v1, v1buf)   
        for (j = 1; j <= X->Check.idim_max; j++) {
            if (GetOffCompGeneralSpin(j - 1, isite, IniSpin, FinSpin, &off,
                                      X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
                dmv = tmp_v1[j] * tmp_V;
                dam_pr += conj(tmp_v1[off + 1]) * dmv;
            }
        }
    }
  return dam_pr;
#endif
}

double complex X_GC_child_CisAitCjuAjv_GeneralSpin_MPIsingle(
							  int org_isite1,
							  int org_ispin1,
							  int org_ispin2,
							  int org_isite3,
							  int org_ispin3,
							  int org_ispin4,
							  double complex tmp_J,
							  struct BindStruct *X,
							  double complex *tmp_v0,
							  double complex *tmp_v1
							  )
{
#ifdef MPI
  unsigned long int off, j;
  int origin, ierr, isite, IniSpin, FinSpin;
  double complex tmp_V, dmv, dam_pr;
  MPI_Status statusMPI;

    if (GetOffCompGeneralSpin((unsigned long int)myrank,
      org_isite3+ 1, org_ispin3, org_ispin4, &off,
      X->Def.SiteToBit, X->Def.Tpow) == TRUE)
    {
      tmp_V = tmp_J;
      isite = org_isite1 + 1;
      IniSpin = org_ispin2;
      FinSpin = org_ispin1;
    }
    else if (GetOffCompGeneralSpin((unsigned long int)myrank,
				   org_isite3+ 1, org_ispin4, org_ispin3, &off,
				   X->Def.SiteToBit, X->Def.Tpow) == TRUE)
    {
      tmp_V = conj(tmp_J);
      if(X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) tmp_V=0.0;
      isite = org_isite1 + 1;
      IniSpin = org_ispin1;
      FinSpin = org_ispin2;
    }
    else return 0.0;

    origin = (int)off;

    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
      v1buf, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
      MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC){
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, isite, IniSpin, FinSpin) private(j, dmv, off) shared (tmp_v0, tmp_v1, v1buf)
        for (j = 1; j <= X->Check.idim_max; j++) {

            if (GetOffCompGeneralSpin(j - 1, isite, IniSpin, FinSpin, &off,
                                      X->Def.SiteToBit, X->Def.Tpow) == TRUE)
            {
                dmv = v1buf[j] * tmp_V;
                tmp_v0[off + 1] += dmv;
                dam_pr += conj(tmp_v1[off + 1]) * dmv;
            }
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, isite, IniSpin, FinSpin) private(j, dmv, off) shared (tmp_v0, tmp_v1, v1buf) 
        for (j = 1; j <= X->Check.idim_max; j++) {

            if (GetOffCompGeneralSpin(j - 1, isite, IniSpin, FinSpin, &off,
                                      X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
                dmv = v1buf[j] * tmp_V;
                dam_pr += conj(tmp_v1[off + 1]) * dmv;
            }
        }
    }
    return dam_pr;
#endif
}


double complex X_GC_child_CisAisCjuAju_GeneralSpin_MPIsingle(
							  int org_isite1,
							  int org_ispin1,
							  int org_isite3,
							  int org_ispin3,
							  double complex tmp_J,
							  struct BindStruct *X,
							  double complex *tmp_v0,
							  double complex *tmp_v1
							  )
{
#ifdef MPI
  unsigned long int j, num1;
  double complex tmp_V, dmv, dam_pr;
  //MPI_Status statusMPI;

  num1 = BitCheckGeneral((unsigned long int)myrank, org_isite3+1, org_ispin3, X->Def.SiteToBit, X->Def.Tpow);
  if(num1 != FALSE){
      tmp_V = tmp_J;
  }
  else return 0.0;
  
  dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, org_isite1, org_ispin1) private(j, dmv, num1) shared (tmp_v0, tmp_v1)
        for (j = 1; j <= X->Check.idim_max; j++) {
            num1 = BitCheckGeneral(j - 1, org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);

            dmv = tmp_v1[j] * tmp_V * num1;
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, org_isite1, org_ispin1) private(j, dmv, num1) shared (tmp_v0, tmp_v1) 
        for (j = 1; j <= X->Check.idim_max; j++) {
            num1 = BitCheckGeneral(j - 1, org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);
            dmv = tmp_v1[j] * tmp_V * num1;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
  return dam_pr;
#endif
}


 /**
 *
 * General interaction term in the Spin model + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void child_general_int_GeneralSpin_MPIdouble(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
  double complex dam_pr;
  dam_pr = X_child_CisAitCjuAjv_GeneralSpin_MPIdouble(X->Def.InterAll_OffDiagonal[i_int][0],
						      X->Def.InterAll_OffDiagonal[i_int][1],X->Def.InterAll_OffDiagonal[i_int][3],
						      X->Def.InterAll_OffDiagonal[i_int][4],
						      X->Def.InterAll_OffDiagonal[i_int][5],X->Def.InterAll_OffDiagonal[i_int][7],
						      X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;

}/*void GC_child_general_int_spin_MPIdouble*/


double complex X_child_CisAitCjuAjv_GeneralSpin_MPIdouble(
							  int org_isite1,
							  int org_ispin1,
							  int org_ispin2,
							  int org_isite3,
							  int org_ispin3,
							  int org_ispin4,
							  double complex tmp_J,
							  struct BindStruct *X,
							  double complex *tmp_v0,
							  double complex *tmp_v1
							  )
{
#ifdef MPI
  unsigned long int tmp_off, off, j, idim_max_buf;
  int origin, ierr;
  double complex tmp_V, dmv, dam_pr;
  MPI_Status statusMPI;
  int ihermite=TRUE;

  if (GetOffCompGeneralSpin((unsigned long int)myrank, org_isite1 + 1, org_ispin1, org_ispin2, &tmp_off, X->Def.SiteToBit, X->Def.Tpow) == TRUE)
  {
    if (GetOffCompGeneralSpin(tmp_off, org_isite3 + 1, org_ispin3, org_ispin4, &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE)
    {
      tmp_V = tmp_J;
    }
    else{
      ihermite =FALSE;
    }
  }
  else{
    ihermite=FALSE;
  }
  
  if(ihermite==FALSE){
    if(GetOffCompGeneralSpin((unsigned long int)myrank, org_isite3 + 1, org_ispin4, org_ispin3, &tmp_off, X->Def.SiteToBit, X->Def.Tpow) == TRUE)
      {
	if (GetOffCompGeneralSpin(tmp_off, org_isite1 + 1, org_ispin2, org_ispin1, &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE)
	  {
	    tmp_V = conj(tmp_J);
	    if(X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC){
	      tmp_V=0.0;
	    }
	  }
	else return 0.0;
      }
    else return 0.0;
  }
  
  
  origin = (int)off;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
    list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, idim_max_buf) private(j, dmv, off) shared (tmp_v0, tmp_v1, list_1buf, v1buf)
        for (j = 1; j <= idim_max_buf; j++) {
            ConvertToList1GeneralSpin(list_1buf[j], X->Check.sdim, &off);
            dmv = v1buf[j] * tmp_V;
            tmp_v0[off] += dmv;
            dam_pr += conj(tmp_v1[off]) * dmv;
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, idim_max_buf) private(j, dmv, off) shared (tmp_v0, tmp_v1, list_1buf, v1buf) 
        for (j = 1; j <= idim_max_buf; j++) {
            ConvertToList1GeneralSpin(list_1buf[j], X->Check.sdim, &off);
            dmv = v1buf[j] * tmp_V;
            dam_pr += conj(tmp_v1[off]) * dmv;
        }
    }
  return dam_pr;
#endif
}

 
double complex X_child_CisAisCjuAju_GeneralSpin_MPIdouble(
							  int org_isite1,
							  int org_ispin1,
							  int org_isite3,
							  int org_ispin3,
							  double complex tmp_J,
							  struct BindStruct *X,
							  double complex *tmp_v0,
							  double complex *tmp_v1
							  ) {
#ifdef MPI
    unsigned long int j, num1;
    double complex tmp_V, dmv, dam_pr;

    if (org_isite1 == org_isite3 && org_ispin1 == org_ispin3) {
        num1 = BitCheckGeneral((unsigned long int) myrank, org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);
        if (num1 != FALSE) {
            tmp_V = tmp_J;
        }
        else {
            return 0.0;
        }
    }
    else {
        num1 = BitCheckGeneral((unsigned long int) myrank, org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);
        if (num1 != FALSE) {
            num1 = BitCheckGeneral((unsigned long int) myrank, org_isite3 + 1, org_ispin3, X->Def.SiteToBit,
                                   X->Def.Tpow);
            if (num1 != FALSE) {
                tmp_V = tmp_J;
            }
            else {
                return 0.0;
            }
        }
        else {
            return 0.0;
        }
    }

    dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1)
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = tmp_v1[j] * tmp_V;
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V) private(j, dmv) shared (tmp_v0, tmp_v1) 
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = tmp_v1[j] * tmp_V;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    return dam_pr;
#endif
}

 double complex X_child_CisAisCjuAju_GeneralSpin_MPIsingle(
							  int org_isite1,
							  int org_ispin1,
							  int org_isite3,
							  int org_ispin3,
							  double complex tmp_J,
							  struct BindStruct *X,
							  double complex *tmp_v0,
							  double complex *tmp_v1
							  )
 {
#ifdef MPI
     unsigned long int j, num1;
     double complex tmp_V, dmv, dam_pr;
     //MPI_Status statusMPI;

     num1 = BitCheckGeneral((unsigned long int) myrank, org_isite3 + 1, org_ispin3, X->Def.SiteToBit, X->Def.Tpow);
     if (num1 != FALSE) {
         tmp_V = tmp_J;
     }
     else return 0.0;

     dam_pr = 0.0;
     if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, org_isite1, org_ispin1) private(j, dmv, num1) shared (tmp_v0, tmp_v1, list_1)
         for (j = 1; j <= X->Check.idim_max; j++) {
             num1 = BitCheckGeneral(list_1[j], org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);

             dmv = tmp_v1[j] * tmp_V * num1;
             tmp_v0[j] += dmv;
             dam_pr += conj(tmp_v1[j]) * dmv;
         }
     }
     else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, org_isite1, org_ispin1) private(j, dmv, num1) shared (tmp_v0, tmp_v1, list_1)
         for (j = 1; j <= X->Check.idim_max; j++) {
             num1 = BitCheckGeneral(list_1[j], org_isite1 + 1, org_ispin1, X->Def.SiteToBit, X->Def.Tpow);

             dmv = tmp_v1[j] * tmp_V * num1;
             dam_pr += conj(tmp_v1[j]) * dmv;
         }
     }

     return dam_pr;
#endif
 }
 
/**
  *
  * General interaction term in the Spin model + GC
  * When both site1 and site2 are in the inter process region.
  *
  * @author Mitsuaki Kawamura (The University of Tokyo)
  */
void child_general_int_GeneralSpin_MPIsingle(
  unsigned long int i_int /**< [in] Interaction ID*/,
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
  double complex dam_pr;

  dam_pr = X_child_CisAitCjuAjv_GeneralSpin_MPIsingle(X->Def.InterAll_OffDiagonal[i_int][0],
						      X->Def.InterAll_OffDiagonal[i_int][1],X->Def.InterAll_OffDiagonal[i_int][3],
						      X->Def.InterAll_OffDiagonal[i_int][4],
						      X->Def.InterAll_OffDiagonal[i_int][5],X->Def.InterAll_OffDiagonal[i_int][7],
						      X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);

  X->Large.prdct += dam_pr;

}/*void GC_child_general_int_spin_MPIsingle*/

double complex X_child_CisAitCjuAjv_GeneralSpin_MPIsingle(
							  int org_isite1,
							  int org_ispin1,
							  int org_ispin2,
							  int org_isite3,
							  int org_ispin3,
							  int org_ispin4,
							  double complex tmp_J,
							  struct BindStruct *X,
							  double complex *tmp_v0,
							  double complex *tmp_v1
							  )
{
#ifdef MPI
  unsigned long int tmp_off, off, j, idim_max_buf;
  int origin, ierr, isite, IniSpin, FinSpin;
  double complex tmp_V, dmv, dam_pr;
  MPI_Status statusMPI;
  
  if (GetOffCompGeneralSpin((unsigned long int)myrank,
    org_isite3 + 1, org_ispin3, org_ispin4, &off,
    X->Def.SiteToBit, X->Def.Tpow) == TRUE)
  {
    tmp_V = tmp_J;
    isite = org_isite1 + 1;
    IniSpin = org_ispin2;
    FinSpin = org_ispin1;
  }
  else if (GetOffCompGeneralSpin((unsigned long int)myrank,
    org_isite3 + 1, org_ispin4, org_ispin3, &off, X->Def.SiteToBit, X->Def.Tpow) == TRUE)
  {
    tmp_V = conj(tmp_J);
    isite = org_isite1 + 1;
    IniSpin = org_ispin1;
    FinSpin = org_ispin2;
  }
  else return 0.0;

  origin = (int)off;
  
  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
    list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, idim_max_buf, IniSpin, FinSpin, isite) private(j, dmv, off, tmp_off) shared (tmp_v0, tmp_v1, list_1buf, v1buf)
        for (j = 1; j <= idim_max_buf; j++) {

            if (GetOffCompGeneralSpin(list_1buf[j], isite, IniSpin, FinSpin, &tmp_off,
                                      X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
                ConvertToList1GeneralSpin(tmp_off, X->Check.sdim, &off);
                dmv = v1buf[j] * tmp_V;
                tmp_v0[off] += dmv;
                dam_pr += conj(tmp_v1[off]) * dmv;
            }
        }
    }
else {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(X, tmp_V, idim_max_buf, IniSpin, FinSpin, isite) private(j, dmv, off, tmp_off) shared (tmp_v0, tmp_v1, list_1buf, v1buf)
        for (j = 1; j <= idim_max_buf; j++) {

            if (GetOffCompGeneralSpin(list_1buf[j], isite, IniSpin, FinSpin, &tmp_off,
                                      X->Def.SiteToBit, X->Def.Tpow) == TRUE) {
                ConvertToList1GeneralSpin(tmp_off, X->Check.sdim, &off);
                dmv = v1buf[j] * tmp_V;
                dam_pr += conj(tmp_v1[off]) * dmv;
            }
        }
    }
  return dam_pr;
#endif
}

/**
 *
 * Hopping term in Spin + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_GC_child_CisAit_spin_MPIdouble(
				       int org_isite1,
				       int org_ispin1,
				       int org_ispin2,
				       double complex tmp_trans,
				       struct BindStruct *X /**< [inout]*/,
				       double complex *tmp_v0 /**< [out] Result v0 = H v1*/, 
				       double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  int mask1, state1, ierr, origin;
  unsigned long int idim_max_buf, j;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;
  
  mask1 = (int)X->Def.Tpow[org_isite1];
  origin = myrank ^ mask1;
  state1 = (origin & mask1)/mask1;

  //fprintf(stdout, "Debug: myrank=%d, origin=%d, state1=%d\n", myrank, origin, state1);

  if(state1 ==  org_ispin2){
    trans = tmp_trans;
  }
  else if(state1 == org_ispin1) {
    trans = conj(tmp_trans);
    if(X->Large.mode == M_CORR|| X->Large.mode ==M_CALCSPEC){
      trans = 0.0;
    }
  }
  else{
    return 0.0;
  }

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
		      &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
		      v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);

    dam_pr = 0.0;
    if(X->Large.mode == M_MLTPLY|| X->Large.mode ==M_CALCSPEC){
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) firstprivate(idim_max_buf, trans, X) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= X->Check.idim_max ; j++) {
            dmv = trans * v1buf[j];
            tmp_v0[j] += dmv;
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
    else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) firstprivate(idim_max_buf, trans, X) shared(v1buf, tmp_v1, tmp_v0)
        for (j = 1; j <= X->Check.idim_max; j++) {
            dmv = trans * v1buf[j];
            dam_pr += conj(tmp_v1[j]) * dmv;
        }
    }
  return (dam_pr);
  
#endif
}/*double complex  X_GC_child_CisAit_spin_MPIdouble*/

/**
 *
 * Hopping term in Spin + Canonical for CalcSpectrum
 * When both site1 and site2 are in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_child_CisAit_spin_MPIdouble(
				       int org_isite1,
				       int org_ispin2,
				       double complex tmp_trans,
				       struct BindStruct *X /**< [inout]*/,
				       double complex *tmp_v0 /**< [out] Result v0 = H v1*/, 
				       double complex *tmp_v1, /**< [in] v0 = H v1*/
                       double complex *tmp_v1buf,
                       unsigned long int idim_max,
                       long unsigned int *Tpow,
                       long unsigned int *list_1_org,
                       long unsigned int *list_1buf_org,
                       long unsigned int *list_2_1_target,
                       long unsigned int *list_2_2_target,
                       const long unsigned int _irght,
                       const long unsigned int _ilft,
                       const long unsigned int _ihfbit
)
{
#ifdef MPI
  int mask1, state1, ierr, origin;
  unsigned long int idim_max_buf, j;
  unsigned long int tmp_off;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;
  
  mask1 = (int)X->Def.Tpow[org_isite1];
  origin = myrank ^ mask1;
  state1 = (origin & mask1)/mask1;

  if(state1 ==  org_ispin2){
    trans = tmp_trans;
  }
  else{
    trans =0.0;
  }

  //  fprintf(stdout, "Debug: myrank=%d, origin=%d, trans=%lf\n", myrank, origin, trans);
  
  ierr = MPI_Sendrecv(&idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
		      &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(list_1_org, idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                        list_1buf_org, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
		      v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if(ierr != 0) exitMPI(-1);
    
    dam_pr = 0.0;
    if(X->Large.mode == M_MLTPLY|| X->Large.mode ==M_CALCSPEC){
#pragma omp parallel for default(none) private(j, dmv, tmp_off) \
  firstprivate(idim_max_buf, trans, X, list_1buf_org,                   \
               list_2_1_target, list_2_2_target) shared(v1buf, tmp_v0)
      for (j = 1; j <= idim_max_buf ; j++) {
        GetOffComp(list_2_1_target, list_2_2_target, list_1buf_org[j], X->Large.irght, X->Large.ilft, X->Large.ihfbit, &tmp_off);
        dmv = trans * v1buf[j];
        tmp_v0[tmp_off] += dmv;
      }
    }
    else {
      tmp_off=0;
      return 0;
    }
    return 1;
  
#endif
}/*double complex  X_child_CisAit_spin_MPIdouble*/

/**
 *
 * Hopping term in Spin + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_GC_child_CisAis_spin_MPIdouble
(
 int org_isite1,
 int org_ispin1,
 double complex tmp_trans,
 struct BindStruct *X /**< [inout]*/,
 double complex *tmp_v0 /**< [out] Result v0 = H v1*/, 
 double complex *tmp_v1 /**< [in] v0 = H v1*/
 )
{
#ifdef MPI
  long unsigned int j;
  int mask1;
  int ibit1;
  double complex dam_pr;
  mask1 = (int)X->Def.Tpow[org_isite1];
  ibit1 = (((unsigned long int)myrank& mask1)/mask1)^(1-org_ispin1);

  dam_pr =0.0;
  if(ibit1 != 0) {
      if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
#pragma omp parallel for reduction(+:dam_pr)default(none) shared(tmp_v1, tmp_v0) \
  firstprivate(X, tmp_trans) private(j)
          for (j = 1; j <= X->Check.idim_max; j++) {
              tmp_v0[j] += tmp_v1[j] * tmp_trans;
              dam_pr += tmp_trans * conj(tmp_v1[j]) * tmp_v1[j];
          }
      }
      else {
#pragma omp parallel for reduction(+:dam_pr)default(none) shared(tmp_v1, tmp_v0) \
  firstprivate(X, tmp_trans) private(j)
          for (j = 1; j <= X->Check.idim_max; j++) {
              dam_pr += tmp_trans * conj(tmp_v1[j]) * tmp_v1[j];
          }
      }
  }

  return dam_pr;
#endif
}

/**
 *
 * Hopping term in Spin + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_GC_child_AisCis_spin_MPIdouble
        (
                int org_isite1,
                int org_ispin1,
                double complex tmp_trans,
                struct BindStruct *X /**< [inout]*/,
                double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
                double complex *tmp_v1 /**< [in] v0 = H v1*/
        )
{
#ifdef MPI
    long unsigned int j;
    int mask1;
    int ibit1;
    double complex dam_pr;
    mask1 = (int)X->Def.Tpow[org_isite1];
    ibit1 = (((unsigned long int)myrank& mask1)/mask1)^(1-org_ispin1);

    dam_pr =0.0;
    if(ibit1 == 0) {
        if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
#pragma omp parallel for reduction(+:dam_pr)default(none) shared(tmp_v1, tmp_v0) \
  firstprivate(X, tmp_trans) private(j)
            for (j = 1; j <= X->Check.idim_max; j++) {
                tmp_v0[j] += tmp_v1[j] * tmp_trans;
                dam_pr += tmp_trans * conj(tmp_v1[j]) * tmp_v1[j];
            }
        }
        else {
#pragma omp parallel for reduction(+:dam_pr)default(none) shared(tmp_v1, tmp_v0) \
  firstprivate(X, tmp_trans) private(j)
            for (j = 1; j <= X->Check.idim_max; j++) {
                dam_pr += tmp_trans * conj(tmp_v1[j]) * tmp_v1[j];
            }
        }
    }

    return dam_pr;
#endif
}


int CheckPE(
	    int org_isite,
	    struct BindStruct *X
	    )
{
  if(org_isite+1 > X->Def.Nsite){
    return TRUE;
  }
  else{
    return FALSE;
  }
}

int CheckBit_Cis(
		 long unsigned int is1_spin,
		 long unsigned int orgbit,
		 long unsigned int *offbit
){
  long unsigned int ibit_tmp;
  ibit_tmp = orgbit & is1_spin;
  if(ibit_tmp == 0){
    *offbit = orgbit+is1_spin;
    return TRUE;
  }
  *offbit=0;
  return FALSE;
}

int CheckBit_Ajt(
		 long unsigned int is1_spin,
		 long unsigned int orgbit,
		 long unsigned int *offbit
){
  long unsigned int ibit_tmp;
  ibit_tmp = orgbit & is1_spin;
  if(ibit_tmp != 0){
    *offbit = orgbit-is1_spin;
    return TRUE;
  }
  *offbit=0;
  return FALSE;
}

int CheckBit_InterAllPE(
			int org_isite1,
			int org_isigma1,
			int org_isite2,
			int org_isigma2,
			int org_isite3,
			int org_isigma3,
			int org_isite4,
			int org_isigma4,
			struct BindStruct *X,
			long unsigned int orgbit,
			long unsigned int *offbit
			)
{
  long unsigned int tmp_ispin;
  long unsigned int tmp_org, tmp_off;
  int iflgBitExist = TRUE;
  tmp_org=orgbit;
  
  if(CheckPE(org_isite1, X)==TRUE){
    tmp_ispin = X->Def.Tpow[2*org_isite1+org_isigma1];
    if(CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE){
      iflgBitExist=FALSE;
    }
    tmp_org = tmp_off;
  }

  if(CheckPE(org_isite2, X)==TRUE){
    tmp_ispin = X->Def.Tpow[2*org_isite2+org_isigma2];
    if(CheckBit_Cis(tmp_ispin, tmp_org, &tmp_off) != TRUE){
      iflgBitExist=FALSE;
    }
    tmp_org = tmp_off;
  }
  
  if(CheckPE(org_isite3, X)==TRUE){
    tmp_ispin = X->Def.Tpow[2*org_isite3+org_isigma3];
    if(CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE){
      iflgBitExist=FALSE;
    }
    tmp_org = tmp_off;
  }

  if(CheckPE(org_isite4, X)==TRUE){
    tmp_ispin = X->Def.Tpow[2*org_isite4+org_isigma4];
    if(CheckBit_Cis(tmp_ispin, tmp_org, &tmp_off) != TRUE){
      iflgBitExist=FALSE;
    }
    tmp_org = tmp_off;
  }

  if(iflgBitExist != TRUE){
    *offbit=0;
    return FALSE;
  }
  
  *offbit=tmp_org;
  return TRUE;
}

int CheckBit_PairPE(
		    int org_isite1,
		    int org_isigma1,
		    int org_isite3,
		    int org_isigma3,
		    struct BindStruct *X,
		    long unsigned int orgbit
		    )
{
  long unsigned int tmp_ispin;
  long unsigned int tmp_org, tmp_off;
  int iflgBitExist = TRUE;
  tmp_org=orgbit;
  
  if(CheckPE(org_isite1, X)==TRUE){
    tmp_ispin = X->Def.Tpow[2*org_isite1+org_isigma1];
    if(CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE){
      iflgBitExist=FALSE;
    }
  }
  
  if(CheckPE(org_isite3, X)==TRUE){
    tmp_ispin = X->Def.Tpow[2*org_isite3+org_isigma3];
    if(CheckBit_Ajt(tmp_ispin, tmp_org, &tmp_off) != TRUE){
      iflgBitExist=FALSE;
    }
  }
  
  if(iflgBitExist != TRUE){
    return FALSE;
  }

  return TRUE;
}

int GetSgnInterAll(
		   int isite1,
		   int isite2,
		   int isite3,
		   int isite4,
		   int *Fsgn,
		   struct BindStruct *X,
		   unsigned long int orgbit,
		   unsigned long int *offbit
		   )
{
  long unsigned int diffA;
  long unsigned int tmp_off;
  long unsigned int tmp_ispin1, tmp_ispin2;
  int tmp_sgn=0;

  tmp_ispin1=isite1;
  tmp_ispin2=isite2;
 
  if(tmp_ispin1 == tmp_ispin2){
    if( (orgbit & tmp_ispin1) == 0){
      *offbit =0;
      *Fsgn = tmp_sgn;
      return FALSE;
    }
    tmp_sgn=1;
    tmp_off =orgbit;
  }
  else{
    
    if(tmp_ispin2 > tmp_ispin1) diffA = tmp_ispin2 - tmp_ispin1*2;
    else diffA = tmp_ispin1-tmp_ispin2*2;  
    
    tmp_sgn=X_GC_CisAjt(orgbit, X, tmp_ispin1, tmp_ispin2, tmp_ispin1+tmp_ispin2, diffA, &tmp_off);
    
    if(tmp_sgn ==0){
      *offbit =0;
      *Fsgn = 0;
      return FALSE;
    }
  }

  tmp_ispin1 = isite3;
  tmp_ispin2 = isite4;
  if(tmp_ispin1 == tmp_ispin2){
    if( (tmp_off & tmp_ispin1) == 0){
      *offbit =0;
      *Fsgn = 0;
      return FALSE;
    }
    *offbit=tmp_off;
  }
  else{
    if(tmp_ispin2 > tmp_ispin1) diffA = tmp_ispin2 - tmp_ispin1*2;
    else diffA = tmp_ispin1-tmp_ispin2*2;  
    
    tmp_sgn *=X_GC_CisAjt(tmp_off, X, tmp_ispin1, tmp_ispin2, tmp_ispin1+tmp_ispin2, diffA, offbit);
    
    if(tmp_sgn ==0){
      *offbit =0;
      *Fsgn = 0;
      return FALSE;
    }
  }
  
  *Fsgn =tmp_sgn;
  *offbit = *offbit%X->Def.OrgTpow[2*X->Def.Nsite];
  
  // exitMPI(-1);
  return TRUE;
}

double complex X_GC_child_CisAisCjtAjt_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 )
{
#ifdef MPI
  double complex dam_pr=0.0;
  int iCheck;
  unsigned long int tmp_ispin1;
  unsigned long int i_max = X->Check.idim_max;
  unsigned long int tmp_off, j;
  double complex dmv;
//  MPI_Status statusMPI;

  iCheck=CheckBit_PairPE(org_isite1, org_ispin1, org_isite3, org_ispin3, X, (long unsigned int) myrank);
  if(iCheck != TRUE){
    return 0.0;
  }
  if(org_isite1+1 > X->Def.Nsite && org_isite3+1 > X->Def.Nsite) {
      if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X) private(dmv, j, tmp_off)
          for (j = 1; j <= i_max; j++) {
              dmv = tmp_v1[j] * tmp_V;
              tmp_v0[j] += dmv;
              dam_pr += conj(tmp_v1[j]) * dmv;
          }
      }
      else {
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X) private(dmv, j, tmp_off)
          for (j = 1; j <= i_max; j++) {
              dmv = tmp_v1[j] * tmp_V;
              dam_pr += conj(tmp_v1[j]) * dmv;
          }
      }
  }
  else if (org_isite1+1 > X->Def.Nsite || org_isite3+1 > X->Def.Nsite) {
      if (org_isite1 > org_isite3) {
          tmp_ispin1 = X->Def.Tpow[2 * org_isite3 + org_ispin3];
      }
      else {
          tmp_ispin1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
      }
      if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) { // for multply
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X, tmp_ispin1) private(dmv, j, tmp_off)
          for (j = 1; j <= i_max; j++) {
              if (CheckBit_Ajt(tmp_ispin1, j - 1, &tmp_off) == TRUE) {
                  dmv = tmp_v1[j] * tmp_V;
                  tmp_v0[j] += dmv;
                  dam_pr += conj(tmp_v1[j]) * dmv;

              }
          }
      }
      else {
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X, tmp_ispin1) private(dmv, j, tmp_off)
          for (j = 1; j <= i_max; j++) {
              if (CheckBit_Ajt(tmp_ispin1, j - 1, &tmp_off) == TRUE) {
                  dmv = tmp_v1[j] * tmp_V;
                  dam_pr += conj(tmp_v1[j]) * dmv;
              }
          }
      }
  }
  return dam_pr;
#endif
}

double complex X_GC_child_CisAjtCkuAku_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite2,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 ) {
#ifdef MPI
    double complex dam_pr = 0.0;
    unsigned long int i_max = X->Check.idim_max;
    unsigned long int idim_max_buf;
    int iCheck, ierr, Fsgn;
    unsigned long int isite1, isite2, isite3;
    unsigned long int tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4;
    unsigned long int j, Asum, Adiff;
    double complex dmv;
    unsigned long int origin, tmp_off;
    unsigned long int org_rankbit;
    MPI_Status statusMPI;

    iCheck = CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite3,
                                 org_ispin3, X, (long unsigned int) myrank, &origin);
    isite1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
    isite2 = X->Def.Tpow[2 * org_isite2 + org_ispin2];
    isite3 = X->Def.Tpow[2 * org_isite3 + org_ispin3];

    if (iCheck == TRUE) {
        tmp_isite1 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
        tmp_isite2 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
        tmp_isite3 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
        tmp_isite4 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
        Asum = tmp_isite1 + tmp_isite2;
        if (tmp_isite2 > tmp_isite1) Adiff = tmp_isite2 - tmp_isite1 * 2;
        else Adiff = tmp_isite1 - tmp_isite2 * 2;
    }
    else {
        iCheck = CheckBit_InterAllPE(org_isite3, org_ispin3, org_isite3, org_ispin3, org_isite2, org_ispin2, org_isite1,
                                     org_ispin1, X, (long unsigned int) myrank, &origin);
        if (iCheck == TRUE) {
            tmp_V = conj(tmp_V);
            tmp_isite4 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
            tmp_isite3 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
            tmp_isite2 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
            tmp_isite1 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
            Asum = tmp_isite3 + tmp_isite4;
            if (tmp_isite4 > tmp_isite3) Adiff = tmp_isite4 - tmp_isite3 * 2;
            else Adiff = tmp_isite3 - tmp_isite4 * 2;
            if (X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) {
                tmp_V = 0;
            }
        }
        else {
            return 0.0;
        }
    }

    if (myrank == origin) {// only k is in PE

        if (CheckBit_Ajt(isite3, myrank, &tmp_off) == FALSE) return 0.0;

#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max,X,Asum,Adiff,isite1,isite2, tmp_V) private(j,tmp_off) shared(tmp_v0, tmp_v1)
        for (j = 1; j <= i_max; j++) {
            dam_pr += GC_CisAjt(j, tmp_v0, tmp_v1, X, isite2, isite1, Asum, Adiff, tmp_V, &tmp_off);
        }

        if (X->Large.mode != M_CORR) {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max,X,Asum,Adiff,isite1,isite2, tmp_V) private(j,tmp_off) shared(tmp_v0, tmp_v1)
            for (j = 1; j <= i_max; j++) {
                dam_pr += GC_CisAjt(j, tmp_v0, tmp_v1, X, isite1, isite2, Asum, Adiff, tmp_V, &tmp_off);
            }
        }
        return dam_pr;
    }//myrank =origin
    else {
        ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                            &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
        if (ierr != 0) exitMPI(-1);
        ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                            v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
        if (ierr != 0) exitMPI(-1);

        if (org_isite1 + 1 > X->Def.Nsite && org_isite2 + 1 > X->Def.Nsite) {
            if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
            else Adiff = isite1 - isite2 * 2;
            SgnBit(((long unsigned int) myrank & Adiff), &Fsgn);
            tmp_V *= Fsgn;

            if (org_isite3 + 1 > X->Def.Nsite) {
                if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) firstprivate(idim_max_buf, tmp_V, X) shared(v1buf, tmp_v1, tmp_v0)
                    for (j = 1; j <= idim_max_buf; j++) {
                        dmv = tmp_V * v1buf[j];
                        tmp_v0[j] += dmv;
                        dam_pr += conj(tmp_v1[j]) * dmv;
                    }


                }
                else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) firstprivate(idim_max_buf, tmp_V, X) shared(v1buf, tmp_v1, tmp_v0)
                    for (j = 1; j <= idim_max_buf; j++) {
                        dmv = tmp_V * v1buf[j];
                        dam_pr += conj(tmp_v1[j]) * dmv;
                    }
                }
            }
            else { //org_isite3 <= X->Def.Nsite
                if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, tmp_off) firstprivate(idim_max_buf, tmp_V, X, isite3) shared(v1buf, tmp_v1, tmp_v0)
                    for (j = 1; j <= idim_max_buf; j++) {
                        if (CheckBit_Ajt(isite3, j - 1, &tmp_off) == TRUE) {
                            dmv = tmp_V * v1buf[j];
                            tmp_v0[j] += dmv;
                            dam_pr += conj(tmp_v1[j]) * dmv;
                        }
                    }
                }
                else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, tmp_off) firstprivate(idim_max_buf, tmp_V, X, isite3) shared(v1buf, tmp_v1, tmp_v0)
                    for (j = 1; j <= idim_max_buf; j++) {
                        if (CheckBit_Ajt(isite3, j - 1, &tmp_off) == TRUE) {
                            dmv = tmp_V * v1buf[j];
                            dam_pr += conj(tmp_v1[j]) * dmv;
                        }
                    }
                }
            }
        }
        else {
            org_rankbit = X->Def.OrgTpow[2 * X->Def.Nsite] * origin;
            if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {

#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, tmp_off, Fsgn) firstprivate(idim_max_buf, tmp_V, X, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit) shared(v1buf, tmp_v1, tmp_v0)
                for (j = 1; j <= idim_max_buf; j++) {
                    /*
                    if(GetSgnInterAll(tmp_isite3, tmp_isite4, tmp_isite1, tmp_isite2, &Fsgn, X, (j-1)+org_rankbit, &tmp_off)==TRUE){
                    */
                    if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, X, (j - 1) + org_rankbit,
                                       &tmp_off) == TRUE) {
                        dmv = tmp_V * v1buf[j] * Fsgn;
                        tmp_v0[tmp_off + 1] += dmv;
                        dam_pr += conj(tmp_v1[tmp_off + 1]) * dmv;

                    }
                }
            }
            else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, tmp_off, Fsgn) firstprivate(idim_max_buf, tmp_V, X, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit) shared(v1buf, tmp_v1, tmp_v0)
                for (j = 1; j <= idim_max_buf; j++) {
                    /*
                    if(GetSgnInterAll(tmp_isite3, tmp_isite4, tmp_isite1, tmp_isite2, &Fsgn, X, (j-1)+org_rankbit, &tmp_off)==TRUE){
                    */
                    if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, X, (j - 1) + org_rankbit,
                                       &tmp_off) == TRUE) {
                        dmv = tmp_V * v1buf[j] * Fsgn;
                        dam_pr += conj(tmp_v1[tmp_off + 1]) * dmv;

                    }
                }
            }
        }
    }
    return dam_pr;
#endif
}


double complex X_GC_child_CisAisCjtAku_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 int org_isite4,
 int org_ispin4,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 )
{
#ifdef MPI
  double complex dam_pr=0;
  dam_pr=X_GC_child_CisAjtCkuAku_Hubbard_MPI
    (
     org_isite4, org_ispin4, org_isite3, org_ispin3,
     org_isite1, org_ispin1, conj(tmp_V), X, tmp_v0, tmp_v1
     );

  return conj(dam_pr);
#endif
}


double complex X_GC_child_CisAjtCkuAlv_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite2,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 int org_isite4,
 int org_ispin4,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 ) {
#ifdef MPI
    double complex dam_pr = 0;
    unsigned long int i_max = X->Check.idim_max;
    unsigned long int idim_max_buf;
    int iCheck, ierr, Fsgn;
    unsigned long int isite1, isite2, isite3, isite4;
    unsigned long int tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4;
    unsigned long int j, Adiff, Bdiff;
    double complex dmv;
    unsigned long int origin, tmp_off, tmp_off2;
    unsigned long int org_rankbit;
    int iFlgHermite = FALSE;
    MPI_Status statusMPI;

    iCheck = CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2,
                                 org_isite3, org_ispin3, org_isite4, org_ispin4,
                                 X, (long unsigned int) myrank, &origin);
    isite1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
    isite2 = X->Def.Tpow[2 * org_isite2 + org_ispin2];
    isite3 = X->Def.Tpow[2 * org_isite3 + org_ispin3];
    isite4 = X->Def.Tpow[2 * org_isite4 + org_ispin4];

    if (iCheck == TRUE) {
        tmp_isite1 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
        tmp_isite2 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
        tmp_isite3 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
        tmp_isite4 = X->Def.OrgTpow[2 * org_isite4 + org_ispin4];
    }
    else {
        iCheck = CheckBit_InterAllPE(org_isite4, org_ispin4, org_isite3, org_ispin3,
                                     org_isite2, org_ispin2, org_isite1, org_ispin1,
                                     X, (long unsigned int) myrank, &origin);
        if (iCheck == TRUE) {
            tmp_V = conj(tmp_V);
            tmp_isite4 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
            tmp_isite3 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
            tmp_isite2 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
            tmp_isite1 = X->Def.OrgTpow[2 * org_isite4 + org_ispin4];
            iFlgHermite = TRUE;
            if (X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) {
                tmp_V = 0;
            }
        }
        else {
            return 0.0;
        }
    }

    if (myrank == origin) {
        if (isite1 == isite4 && isite2 == isite3) { // CisAjvCjvAis =Cis(1-njv)Ais=nis-nisnjv
            //calc nis
            dam_pr = X_GC_child_CisAis_Hubbard_MPI(org_isite1, org_ispin1, tmp_V, X, tmp_v0, tmp_v1);
            //calc -nisniv
            dam_pr -= X_GC_child_CisAisCjtAjt_Hubbard_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3, tmp_V, X,
                                                          tmp_v0, tmp_v1);
        }
        else if (isite2 == isite3) { // CisAjvCjvAku= Cis(1-njv)Aku=-CisAkunjv+CisAku: j is in PE
            //calc CisAku
            if (isite4 > isite1) Adiff = isite4 - isite1 * 2;
            else Adiff = isite1 - isite4 * 2;

#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_off) firstprivate(i_max, tmp_V, X, isite1, isite4, Adiff) shared(tmp_v1, tmp_v0)
            for (j = 1; j <= i_max; j++) {
                dam_pr += GC_CisAjt(j - 1, tmp_v0, tmp_v1, X, isite1, isite4, (isite1 + isite4), Adiff, tmp_V,
                                    &tmp_off);
            }
            //calc -CisAku njv
            dam_pr -= X_GC_child_CisAjtCkuAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite4, org_ispin4, org_isite2,
                                                          org_ispin2, tmp_V, X, tmp_v0, tmp_v1);
            if (X->Large.mode != M_CORR) { //for hermite
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_off) firstprivate(i_max, tmp_V, X, isite1, isite4, Adiff) shared(tmp_v1, tmp_v0)
                for (j = 1; j <= i_max; j++) {
                    dam_pr += GC_CisAjt(j - 1, tmp_v0, tmp_v1, X, isite4, isite1, (isite1 + isite4), Adiff, tmp_V,
                                        &tmp_off);
                }
                //calc -njvCkuAis
                dam_pr -= X_GC_child_CisAisCjtAku_Hubbard_MPI(org_isite2, org_ispin2, org_isite4, org_ispin4,
                                                              org_isite1, org_ispin1, tmp_V, X, tmp_v0, tmp_v1);

            }
        }
        else {// CisAjtCkuAis = -CisAisCkuAjt: i is in PE

            dam_pr = -X_GC_child_CisAisCjtAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3, org_isite2,
                                                          org_ispin2, tmp_V, X, tmp_v0, tmp_v1);
            if (X->Large.mode != M_CORR) { //for hermite

                dam_pr += -X_GC_child_CisAisCjtAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite2, org_ispin2,
                                                               org_isite3, org_ispin3, tmp_V, X, tmp_v0, tmp_v1);
            }
        }
        return dam_pr;
    }//myrank =origin
    else {

        ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                            &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
        if (ierr != 0) exitMPI(-1);
        ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                            v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
        if (ierr != 0) exitMPI(-1);

        if (org_isite1 + 1 > X->Def.Nsite && org_isite2 + 1 > X->Def.Nsite
            && org_isite3 + 1 > X->Def.Nsite && org_isite4 + 1 > X->Def.Nsite) {

            if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
            else Adiff = isite1 - isite2 * 2;
            if (isite4 > isite3) Bdiff = isite4 - isite3 * 2;
            else Bdiff = isite3 - isite4 * 2;

            if (iFlgHermite == FALSE) {
                Fsgn = X_GC_CisAjt((long unsigned int) myrank, X, isite2, isite1, (isite1 + isite2), Adiff, &tmp_off2);
                Fsgn *= X_GC_CisAjt(tmp_off2, X, isite4, isite3, (isite3 + isite4), Bdiff, &tmp_off);
                tmp_V *= Fsgn;
            }
            else {
                Fsgn = X_GC_CisAjt((long unsigned int) myrank, X, isite3, isite4, (isite3 + isite4), Bdiff, &tmp_off2);
                Fsgn *= X_GC_CisAjt(tmp_off2, X, isite1, isite2, (isite1 + isite2), Adiff, &tmp_off);
                tmp_V *= Fsgn;
            }
            if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {

#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) firstprivate(idim_max_buf, tmp_V, X) shared(v1buf, tmp_v1, tmp_v0)
                for (j = 1; j <= idim_max_buf; j++) {
                    dmv = tmp_V * v1buf[j];
                    tmp_v0[j] += dmv;
                    dam_pr += conj(tmp_v1[j]) * dmv;
                }
            }
            else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) firstprivate(idim_max_buf, tmp_V, X) shared(v1buf, tmp_v1, tmp_v0)
                for (j = 1; j <= idim_max_buf; j++) {
                    dmv = tmp_V * v1buf[j];
                    dam_pr += conj(tmp_v1[j]) * dmv;
                }
            }
        }
        else {
            org_rankbit = X->Def.OrgTpow[2 * X->Def.Nsite] * origin;
            if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, tmp_off, Fsgn) firstprivate(idim_max_buf, tmp_V, X, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit) shared(v1buf, tmp_v1, tmp_v0)
                for (j = 1; j <= idim_max_buf; j++) {

                    if (GetSgnInterAll(tmp_isite3, tmp_isite4, tmp_isite1, tmp_isite2, &Fsgn, X, (j - 1) + org_rankbit,
                                       &tmp_off) == TRUE) {
                        dmv = tmp_V * v1buf[j] * Fsgn;
                        tmp_v0[tmp_off + 1] += dmv;
                        dam_pr += conj(tmp_v1[tmp_off + 1]) * dmv;
                    }
                }
            }
            else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, tmp_off, Fsgn) firstprivate(idim_max_buf, tmp_V, X, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit) shared(v1buf, tmp_v1, tmp_v0)
                for (j = 1; j <= idim_max_buf; j++) {

                    if (GetSgnInterAll(tmp_isite3, tmp_isite4, tmp_isite1, tmp_isite2, &Fsgn, X, (j - 1) + org_rankbit,
                                       &tmp_off) == TRUE) {
                        dmv = tmp_V * v1buf[j] * Fsgn;
                        dam_pr += conj(tmp_v1[tmp_off + 1]) * dmv;
                    }
                }
            }
        }
    }

    return dam_pr;
#endif
}

double complex X_GC_child_CisAis_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 ) {
#ifdef MPI
    double complex dam_pr = 0.0;
    unsigned long int i_max = X->Check.idim_max;
    unsigned long int j, isite1, tmp_off;
    double complex dmv;
//  MPI_Status statusMPI;

    isite1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
    if (org_isite1 + 1 > X->Def.Nsite) {
        if (CheckBit_Ajt(isite1, (unsigned long int) myrank, &tmp_off) == FALSE) return 0.0;

        if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) { // for multply

#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X) private(dmv, j, tmp_off)
            for (j = 1; j <= i_max; j++) {
                dmv = tmp_v1[j] * tmp_V;
                tmp_v0[j] += dmv;
                dam_pr += conj(tmp_v1[j]) * dmv;
            }
        }
        else {
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X) private(dmv, j, tmp_off)
            for (j = 1; j <= i_max; j++) {
                dmv = tmp_v1[j] * tmp_V;
                dam_pr += conj(tmp_v1[j]) * dmv;
            }
        }
    }
    else {
        if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) { // for multply

#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X, isite1) private(dmv, j, tmp_off)
            for (j = 1; j <= i_max; j++) {
                if (CheckBit_Ajt(isite1, j - 1, &tmp_off) == TRUE) {
                    dmv = tmp_v1[j] * tmp_V;
                    tmp_v0[j] += dmv;
                    dam_pr += conj(tmp_v1[j]) * dmv;
                }
            }
        }
        else {
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X, isite1) private(dmv, j, tmp_off)
            for (j = 1; j <= i_max; j++) {
                if (CheckBit_Ajt(isite1, j - 1, &tmp_off) == TRUE) {
                    dmv = tmp_v1[j] * tmp_V;
                    dam_pr += conj(tmp_v1[j]) * dmv;
                }
            }

        }
    }

    return dam_pr;
#endif
}

double complex X_GC_child_CisAjt_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite2,
 int org_ispin2,
 double complex tmp_trans,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 ){
  #ifdef MPI
  double complex dam_pr=0.0;
//  MPI_Status statusMPI;

  if(org_isite1 + 1 > X->Def.Nsite && org_isite2+1 > X->Def.Nsite){
    dam_pr = X_GC_child_general_hopp_MPIdouble(org_isite1, org_ispin1, org_isite2, org_ispin2, tmp_trans, X, tmp_v0, tmp_v1);    
  }
  else if(org_isite1 +1 > X->Def.Nsite || org_isite2+1>X->Def.Nsite){
    dam_pr = X_GC_child_general_hopp_MPIsingle(org_isite1, org_ispin1, org_isite2, org_ispin2, tmp_trans, X, tmp_v0, tmp_v1);    
  }
  else{
    //error message will be added.
    exitMPI(-1);
  }  
  return dam_pr;
  #endif
}

double complex X_child_CisAisCjtAjt_Hubbard_MPI
(int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1) {
#ifdef MPI
    double complex dam_pr = 0.0;
    int iCheck;
    unsigned long int tmp_ispin1;
    unsigned long int i_max = X->Check.idim_max;
    unsigned long int tmp_off, j;
    double complex dmv;
//  MPI_Status statusMPI;

    iCheck = CheckBit_PairPE(org_isite1, org_ispin1, org_isite3, org_ispin3, X, (long unsigned int) myrank);
    if (iCheck != TRUE) {
        return 0.0;
    }
    if (org_isite1 + 1 > X->Def.Nsite && org_isite3 + 1 > X->Def.Nsite) {
        if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) { // for multply
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X) private(dmv, j, tmp_off)
            for (j = 1; j <= i_max; j++) {
                dmv = tmp_v1[j] * tmp_V;
                tmp_v0[j] += dmv;
                dam_pr += conj(tmp_v1[j]) * dmv;
            }

        }
        else {
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X) private(dmv, j, tmp_off)
            for (j = 1; j <= i_max; j++) {
                dmv = tmp_v1[j] * tmp_V;
                dam_pr += conj(tmp_v1[j]) * dmv;
            }
        }
    }
    else if (org_isite1 + 1 > X->Def.Nsite || org_isite3 + 1 > X->Def.Nsite) {
        if (org_isite1 > org_isite3) {
            tmp_ispin1 = X->Def.Tpow[2 * org_isite3 + org_ispin3];
        }
        else {
            tmp_ispin1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
        }

        if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1, list_1) \
  firstprivate(i_max, tmp_V, X, tmp_ispin1) private(dmv, j, tmp_off)
            for (j = 1; j <= i_max; j++) {
                if (CheckBit_Ajt(tmp_ispin1, list_1[j], &tmp_off) == TRUE) {
                    dmv = tmp_v1[j] * tmp_V;
                    tmp_v0[j] += dmv;
                    dam_pr += conj(tmp_v1[j]) * dmv;
                }
            }
        }
        else {
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1, list_1) \
  firstprivate(i_max, tmp_V, X, tmp_ispin1) private(dmv, j, tmp_off)
            for (j = 1; j <= i_max; j++) {
                if (CheckBit_Ajt(tmp_ispin1, list_1[j], &tmp_off) == TRUE) {
                    dmv = tmp_v1[j] * tmp_V;
                    dam_pr += conj(tmp_v1[j]) * dmv;
                }
            }
        }
    }
    return dam_pr;
#endif
}


double complex X_child_CisAjtCkuAlv_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite2,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 int org_isite4,
 int org_ispin4,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 ) {
#ifdef MPI
    double complex dam_pr = 0;
    unsigned long int i_max = X->Check.idim_max;
    unsigned long int idim_max_buf;
    int iCheck, ierr, Fsgn;
    unsigned long int isite1, isite2, isite3, isite4;
    unsigned long int tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4;
    unsigned long int j, Adiff, Bdiff;
    double complex dmv;
    unsigned long int origin, tmp_off, tmp_off2;
    unsigned long int org_rankbit, ioff;
    int iFlgHermite = FALSE;
    MPI_Status statusMPI;

    iCheck = CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2,
                                 org_isite3, org_ispin3, org_isite4, org_ispin4,
                                 X, (long unsigned int) myrank, &origin);
    isite1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
    isite2 = X->Def.Tpow[2 * org_isite2 + org_ispin2];
    isite3 = X->Def.Tpow[2 * org_isite3 + org_ispin3];
    isite4 = X->Def.Tpow[2 * org_isite4 + org_ispin4];

    //printf("debug: myrank=%d, isite1=%d, ispin1=%d, isite2=%d, ispin2=%d, isite3=%d, ispin3=%d, isite4=%d, ispin4=%d\n", myrank, org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite4, org_ispin4);
    if (iCheck == TRUE) {
        tmp_isite1 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
        tmp_isite2 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
        tmp_isite3 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
        tmp_isite4 = X->Def.OrgTpow[2 * org_isite4 + org_ispin4];
    }
    else {
        iCheck = CheckBit_InterAllPE(org_isite4, org_ispin4, org_isite3, org_ispin3,
                                     org_isite2, org_ispin2, org_isite1, org_ispin1,
                                     X, (long unsigned int) myrank, &origin);
        if (iCheck == TRUE) {
            tmp_V = conj(tmp_V);
            tmp_isite4 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
            tmp_isite3 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
            tmp_isite2 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
            tmp_isite1 = X->Def.OrgTpow[2 * org_isite4 + org_ispin4];
            iFlgHermite = TRUE;
            if (X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) {
                tmp_V = 0;
            }
        }
        else {
            return 0.0;
        }
    }

    if (myrank == origin) {
        if (isite1 == isite4 && isite2 == isite3) { // CisAjvCjvAis =Cis(1-njv)Ais=nis-nisnjv
            //calc nis
            dam_pr = X_child_CisAis_Hubbard_MPI(org_isite1, org_ispin1, tmp_V, X, tmp_v0, tmp_v1);
            //calc -nisniv
            dam_pr -= X_child_CisAisCjtAjt_Hubbard_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3, tmp_V, X, tmp_v0,
                                                       tmp_v1);
        }
        else if (isite2 == isite3) { // CisAjvCjvAku= Cis(1-njv)Aku=-CisAkunjv+CisAku: j is in PE
            if (isite4 > isite1) Adiff = isite4 - isite1 * 2;
            else Adiff = isite1 - isite4 * 2;

            //calc CisAku
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_off) firstprivate(i_max, tmp_V, X, isite1, isite4, Adiff) shared(tmp_v1, tmp_v0, list_1)
            for (j = 1; j <= i_max; j++) {
                dam_pr += CisAjt(j, tmp_v0, tmp_v1, X, isite1, isite4, (isite1 + isite4), Adiff, tmp_V);
            }
            //calc -CisAku njv
            dam_pr -= X_child_CisAjtCkuAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite4, org_ispin4, org_isite2,
                                                       org_ispin2, tmp_V, X, tmp_v0, tmp_v1);

            if (X->Large.mode != M_CORR) {  //for hermite
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_off) firstprivate(i_max, tmp_V, X, isite1, isite4, Adiff) shared(tmp_v1, tmp_v0)
                for (j = 1; j <= i_max; j++) {
                    dam_pr += CisAjt(j, tmp_v0, tmp_v1, X, isite4, isite1, (isite1 + isite4), Adiff, tmp_V);
                }
                //calc -njvCkuAis
                dam_pr -= X_child_CisAisCjtAku_Hubbard_MPI(org_isite2, org_ispin2, org_isite4, org_ispin4, org_isite1,
                                                           org_ispin1, tmp_V, X, tmp_v0, tmp_v1);
            }
        }
        else {// CisAjtCkuAis = -CisAisCkuAjt: i is in PE
            dam_pr = -X_child_CisAisCjtAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite3, org_ispin3, org_isite2,
                                                       org_ispin2, tmp_V, X, tmp_v0, tmp_v1);

            if (X->Large.mode != M_CORR) {  //for hermite: CisAkuCjtAis=-CisAisCjtAku
                dam_pr = -X_child_CisAisCjtAku_Hubbard_MPI(org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3,
                                                           org_ispin3, tmp_V, X, tmp_v0, tmp_v1);
            }
        }

        return dam_pr;
    }//myrank =origin
    else {
        //printf("debug: myrank=%d, origin=%d\n", myrank, origin);
        ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                            &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
        if (ierr != 0) exitMPI(-1);

        ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                            list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
        if (ierr != 0) exitMPI(-1);

        ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                            v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
        if (ierr != 0) exitMPI(-1);

        if (org_isite1 + 1 > X->Def.Nsite && org_isite2 + 1 > X->Def.Nsite
            && org_isite3 + 1 > X->Def.Nsite && org_isite4 + 1 > X->Def.Nsite) {

            if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
            else Adiff = isite1 - isite2 * 2;
            if (isite4 > isite3) Bdiff = isite4 - isite3 * 2;
            else Bdiff = isite3 - isite4 * 2;

            if (iFlgHermite == FALSE) {
                Fsgn = X_GC_CisAjt((long unsigned int) myrank, X, isite2, isite1, (isite1 + isite2), Adiff, &tmp_off2);
                Fsgn *= X_GC_CisAjt(tmp_off2, X, isite4, isite3, (isite3 + isite4), Bdiff, &tmp_off);
                tmp_V *= Fsgn;
            }
            else {
                Fsgn = X_GC_CisAjt((long unsigned int) myrank, X, isite3, isite4, (isite3 + isite4), Bdiff, &tmp_off2);
                Fsgn *= X_GC_CisAjt(tmp_off2, X, isite1, isite2, (isite1 + isite2), Adiff, &tmp_off);
                tmp_V *= Fsgn;
            }
            dam_pr = 0;
            if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff) firstprivate(idim_max_buf, tmp_V, X) shared(v1buf, tmp_v1, tmp_v0, list_2_1, list_2_2, list_1buf)
                for (j = 1; j <= idim_max_buf; j++) {
                    if (GetOffComp(list_2_1, list_2_2, list_1buf[j],
                                   X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == TRUE) {
                        dmv = tmp_V * v1buf[j];
                        tmp_v0[ioff] += dmv;
                        dam_pr += conj(tmp_v1[ioff]) * dmv;
                    }
                }
            }
            else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff) firstprivate(idim_max_buf, tmp_V, X) shared(v1buf, tmp_v1, tmp_v0, list_2_1, list_2_2, list_1buf)
                for (j = 1; j <= idim_max_buf; j++) {
                    if (GetOffComp(list_2_1, list_2_2, list_1buf[j],
                                   X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == TRUE) {
                        dmv = tmp_V * v1buf[j];
                        dam_pr += conj(tmp_v1[ioff]) * dmv;
                    }
                }
            }
        }//org_isite1+1 > X->Def.Nsite && org_isite2+1 > X->Def.Nsite
            // && org_isite3+1 > X->Def.Nsite && org_isite4+1 > X->Def.Nsite
        else {

            org_rankbit = X->Def.OrgTpow[2 * X->Def.Nsite] * origin;
            dam_pr = 0;
            if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {

#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, tmp_off, Fsgn, ioff) firstprivate(idim_max_buf, tmp_V, X, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit, org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite4, org_ispin4) shared(v1buf, tmp_v1, tmp_v0, list_1buf, list_2_1, list_2_2)
                for (j = 1; j <= idim_max_buf; j++) {
                    if (GetSgnInterAll(tmp_isite3, tmp_isite4, tmp_isite1, tmp_isite2, &Fsgn, X,
                                       list_1buf[j] + org_rankbit,
                                       &tmp_off) == TRUE) {

                        if (GetOffComp(list_2_1, list_2_2, tmp_off,
                                       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == TRUE) {
                            dmv = tmp_V * v1buf[j] * Fsgn;
                            tmp_v0[ioff] += dmv;
                            dam_pr += conj(tmp_v1[ioff]) * dmv;

                        }
                    }
                }
            }
            else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, tmp_off, Fsgn, ioff) firstprivate(idim_max_buf, tmp_V, X, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit, org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite4, org_ispin4) shared(v1buf, tmp_v1, tmp_v0, list_1buf, list_2_1, list_2_2)
                for (j = 1; j <= idim_max_buf; j++) {
                    if (GetSgnInterAll(tmp_isite3, tmp_isite4, tmp_isite1, tmp_isite2, &Fsgn, X,
                                       list_1buf[j] + org_rankbit, &tmp_off) == TRUE) {

                        if (GetOffComp(list_2_1, list_2_2, tmp_off,
                                       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff) == TRUE) {
                            dmv = tmp_V * v1buf[j] * Fsgn;
                            dam_pr += conj(tmp_v1[ioff]) * dmv;

                        }
                    }
                }
            }
        }
    }

    return dam_pr;
#endif
}

double complex X_child_CisAjtCkuAku_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite2,
 int org_ispin2,
 int org_isite3,
 int org_ispin3,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 ) {
#ifdef MPI
    double complex dam_pr = 0.0;
    unsigned long int i_max = X->Check.idim_max;
    unsigned long int idim_max_buf, ioff;
    int iCheck, ierr, Fsgn;
    unsigned long int isite1, isite2, isite3;
    unsigned long int tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4;
    unsigned long int j, Asum, Adiff;
    double complex dmv;
    unsigned long int origin, tmp_off;
    unsigned long int org_rankbit;
    MPI_Status statusMPI;

    iCheck = CheckBit_InterAllPE(org_isite1, org_ispin1, org_isite2, org_ispin2, org_isite3, org_ispin3, org_isite3,
                                 org_ispin3, X, (long unsigned int) myrank, &origin);
    isite1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
    isite2 = X->Def.Tpow[2 * org_isite2 + org_ispin2];
    isite3 = X->Def.Tpow[2 * org_isite3 + org_ispin3];

    if (iCheck == TRUE) {

        tmp_isite1 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
        tmp_isite2 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
        tmp_isite3 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
        tmp_isite4 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
        Asum = tmp_isite1 + tmp_isite2;
        if (tmp_isite2 > tmp_isite1) Adiff = tmp_isite2 - tmp_isite1 * 2;
        else Adiff = tmp_isite1 - tmp_isite2 * 2;
    }
    else {
        iCheck = CheckBit_InterAllPE(org_isite3, org_ispin3, org_isite3, org_ispin3, org_isite2, org_ispin2, org_isite1,
                                     org_ispin1, X, (long unsigned int) myrank, &origin);
        if (iCheck == TRUE) {
            tmp_V = conj(tmp_V);
            tmp_isite4 = X->Def.OrgTpow[2 * org_isite1 + org_ispin1];
            tmp_isite3 = X->Def.OrgTpow[2 * org_isite2 + org_ispin2];
            tmp_isite2 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
            tmp_isite1 = X->Def.OrgTpow[2 * org_isite3 + org_ispin3];
            Asum = tmp_isite3 + tmp_isite4;
            if (tmp_isite4 > tmp_isite3) Adiff = tmp_isite4 - tmp_isite3 * 2;
            else Adiff = tmp_isite3 - tmp_isite4 * 2;
            if (X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) {
                tmp_V = 0;
            }
        }
        else {
            return 0.0;
        }
    }

    if (myrank == origin) {// only k is in PE
        //for hermite
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max, Asum, Adiff, isite1, isite2, tmp_V, X) private(j) shared(tmp_v0, tmp_v1)
        for (j = 1; j <= i_max; j++) {
            dam_pr += CisAjt(j, tmp_v0, tmp_v1, X, isite1, isite2, Asum, Adiff, tmp_V);
        }

        if (X->Large.mode != M_CORR) {
#pragma omp parallel for default(none) reduction(+:dam_pr) firstprivate(i_max,X,Asum,Adiff,isite1,isite2, tmp_V) private(j) shared(tmp_v0, tmp_v1)
            for (j = 1; j <= i_max; j++) {
                dam_pr += CisAjt(j, tmp_v0, tmp_v1, X, isite2, isite1, Asum, Adiff, tmp_V);
            }
        }
        return dam_pr;
    }//myrank =origin
    else {
        ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                            &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
        if (ierr != 0) exitMPI(-1);
        ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                            list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
        if (ierr != 0) exitMPI(-1);

        ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                            v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
        if (ierr != 0) exitMPI(-1);

        if (org_isite1 + 1 > X->Def.Nsite && org_isite2 + 1 > X->Def.Nsite) {
            if (isite2 > isite1) Adiff = isite2 - isite1 * 2;
            else Adiff = isite1 - isite2 * 2;
            SgnBit(((long unsigned int) myrank & Adiff), &Fsgn);
            tmp_V *= Fsgn;

            if (org_isite3 + 1 > X->Def.Nsite) {
                if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, ioff, dmv) firstprivate(idim_max_buf, tmp_V, X) shared(v1buf, tmp_v1, tmp_v0, list_1buf, list_2_1, list_2_2)
                    for (j = 1; j <= idim_max_buf; j++) {
                        dmv = tmp_V * v1buf[j];
                        GetOffComp(list_2_1, list_2_2, list_1buf[j],
                                   X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
                        tmp_v0[ioff] += dmv;
                        dam_pr += conj(tmp_v1[ioff]) * dmv;
                    }
                }
                else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, ioff, dmv) firstprivate(idim_max_buf, tmp_V, X) shared(v1buf, tmp_v1, tmp_v0, list_1buf, list_2_1, list_2_2)
                    for (j = 1; j <= idim_max_buf; j++) {
                        dmv = tmp_V * v1buf[j];
                        GetOffComp(list_2_1, list_2_2, list_1buf[j],
                                   X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
                        dam_pr += conj(tmp_v1[ioff]) * dmv;
                    }
                }
            }
            else { //org_isite3 <= X->Def.Nsite
                if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, ioff, dmv, tmp_off) firstprivate(idim_max_buf, tmp_V, X, isite3) shared(v1buf, tmp_v1, tmp_v0,list_1buf, list_2_1, list_2_2)
                    for (j = 1; j <= idim_max_buf; j++) {
                        if (CheckBit_Ajt(isite3, list_1buf[j], &tmp_off) == TRUE) {
                            dmv = tmp_V * v1buf[j];
                            GetOffComp(list_2_1, list_2_2, list_1buf[j],
                                       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
                            tmp_v0[ioff] += dmv;
                            dam_pr += conj(tmp_v1[ioff]) * dmv;
                        }
                    }

                }
                else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, ioff, dmv, tmp_off) firstprivate(idim_max_buf, tmp_V, X, isite3) shared(v1buf, tmp_v1, tmp_v0,list_1buf, list_2_1, list_2_2)
                    for (j = 1; j <= idim_max_buf; j++) {
                        if (CheckBit_Ajt(isite3, list_1buf[j], &tmp_off) == TRUE) {
                            dmv = tmp_V * v1buf[j];
                            GetOffComp(list_2_1, list_2_2, list_1buf[j],
                                       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
                            dam_pr += conj(tmp_v1[ioff]) * dmv;
                        }
                    }
                }
            }
        }
        else {
            org_rankbit = X->Def.OrgTpow[2 * X->Def.Nsite] * origin;
            if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff, tmp_off, Fsgn) firstprivate(idim_max_buf, tmp_V, X, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit) shared(v1buf, tmp_v1, tmp_v0, list_1buf, list_2_1, list_2_2)
                for (j = 1; j <= idim_max_buf; j++) {
                    /*
                    if(GetSgnInterAll(tmp_isite3, tmp_isite4, tmp_isite1, tmp_isite2, &Fsgn, X, list_1buf[j]+org_rankbit, &tmp_off)==TRUE){
                    */
                    if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, X,
                                       list_1buf[j] + org_rankbit, &tmp_off) == TRUE) {
                        dmv = tmp_V * v1buf[j] * Fsgn;
                        GetOffComp(list_2_1, list_2_2, tmp_off,
                                   X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
                        tmp_v0[ioff] += dmv;
                        dam_pr += conj(tmp_v1[ioff]) * dmv;

                    }
                }
            }
            else {
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff, tmp_off, Fsgn) firstprivate(idim_max_buf, tmp_V, X, tmp_isite1, tmp_isite2, tmp_isite3, tmp_isite4, org_rankbit) shared(v1buf, tmp_v1, tmp_v0, list_1buf, list_2_1, list_2_2)
                for (j = 1; j <= idim_max_buf; j++) {
                    /*
                    if(GetSgnInterAll(tmp_isite3, tmp_isite4, tmp_isite1, tmp_isite2, &Fsgn, X, list_1buf[j]+org_rankbit, &tmp_off)==TRUE){
                    */
                    if (GetSgnInterAll(tmp_isite4, tmp_isite3, tmp_isite2, tmp_isite1, &Fsgn, X,
                                       list_1buf[j] + org_rankbit, &tmp_off) == TRUE) {
                        dmv = tmp_V * v1buf[j] * Fsgn;
                        GetOffComp(list_2_1, list_2_2, tmp_off,
                                   X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
                        dam_pr += conj(tmp_v1[ioff]) * dmv;

                    }
                }
            }
        }
    }
    return dam_pr;
#endif
}
  

double complex X_child_CisAisCjtAku_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 int org_isite3,
 int org_ispin3,
 int org_isite4,
 int org_ispin4,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 ){
#ifdef MPI
  double complex dam_pr=0;
  
  dam_pr=X_child_CisAjtCkuAku_Hubbard_MPI
    (
     org_isite4, org_ispin4, org_isite3, org_ispin3,
     org_isite1, org_ispin1, conj(tmp_V), X, tmp_v0, tmp_v1
     );
  
  return conj(dam_pr);
  
#endif
}

double complex X_child_CisAis_Hubbard_MPI
(
 int org_isite1,
 int org_ispin1,
 double complex tmp_V,
 struct BindStruct *X,
 double complex *tmp_v0,
 double complex *tmp_v1
 ) {
#ifdef MPI
    double complex dam_pr = 0.0;
    unsigned long int i_max = X->Check.idim_max;
    unsigned long int j, isite1, tmp_off;
    double complex dmv;
//  MPI_Status statusMPI;

    isite1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];
    if (org_isite1 + 1 > X->Def.Nsite) {
        if (CheckBit_Ajt(isite1, (unsigned long int) myrank, &tmp_off) == FALSE) {
            return 0.0;
        }
        if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) { // for multply
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X) private(dmv, j, tmp_off)
            for (j = 1; j <= i_max; j++) {
                dmv = tmp_v1[j] * tmp_V;
                tmp_v0[j] += dmv;
                dam_pr += conj(tmp_v1[j]) * dmv;
            }
        }
        else {
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1) \
  firstprivate(i_max, tmp_V, X) private(dmv, j, tmp_off)
            for (j = 1; j <= i_max; j++) {
                dmv = tmp_v1[j] * tmp_V;
                dam_pr += conj(tmp_v1[j]) * dmv;
            }
        }
    }
    else {
        if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) { // for multply
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1, list_1) \
  firstprivate(i_max, tmp_V, X, isite1) private(dmv, j, tmp_off)
            for (j = 1; j <= i_max; j++) {
                if (X_CisAis(list_1[j], X, isite1) != 0) {
                    dmv = tmp_v1[j] * tmp_V;
                    tmp_v0[j] += dmv;
                    dam_pr += conj(tmp_v1[j]) * dmv;
                }
            }

        }
        else {
#pragma omp parallel for reduction(+:dam_pr) default(none) shared(tmp_v0, tmp_v1, list_1) \
  firstprivate(i_max, tmp_V, X, isite1) private(dmv, j, tmp_off)
            for (j = 1; j <= i_max; j++) {
                if (X_CisAis(list_1[j], X, isite1) != 0) {
                    dmv = tmp_v1[j] * tmp_V;
                    dam_pr += conj(tmp_v1[j]) * dmv;
                }
            }
        }
    }
    return dam_pr;
#endif
}
