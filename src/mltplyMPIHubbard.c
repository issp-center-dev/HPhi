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
/**@file
@brief Functions for Hubbard Hamiltonian + MPI
*/
#ifdef MPI
#include "mpi.h"
#endif
#include "Common.h"
#include "bitcalc.h"
#include "wrapperMPI.h"
#include "mltplyCommon.h"
#include "mltplyMPIHubbard.h"
/**
@brief Hopping term in Hubbard + GC
When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_general_hopp_MPIdouble
(
 unsigned long int itrans,//!<[in] Transfer ID
 struct BindStruct *X,//!<[inout]
 double complex *tmp_v0,//!<[out] Result v0 = H v1
 double complex *tmp_v1//!<[in] v0 = H v1
){
#ifdef MPI
  double complex dam_pr = 0;
  dam_pr = child_GC_general_hopp_MPIdouble(
    X->Def.EDGeneralTransfer[itrans][0], X->Def.EDGeneralTransfer[itrans][1],
    X->Def.EDGeneralTransfer[itrans][2], X->Def.EDGeneralTransfer[itrans][3],
    X->Def.EDParaGeneralTransfer[itrans], X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}/*void GC_general_hopp_MPIdouble*/
/**
@brief Hopping term in Hubbard + GC
When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
@return fragment of @f$\langle v_1|{\hat H}|v_1\rangle@f$
*/
double complex child_GC_general_hopp_MPIdouble(
  int org_isite1,//!<[in] @f$i_1@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  int org_ispin1,//!<[in] @f$\sigma_1@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  int org_isite2,//!<[in] @f$i_2@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  int org_ispin2,//!<[in] @f$\sigma_2@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  double complex tmp_trans,//!<[in] Transfer @f$t@f$
  struct BindStruct *X,//!< [inout]
  double complex *tmp_v0,//!< [out] Result v0 = H v1
  double complex *tmp_v1 //!< [in] v0 = H v1
) {
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin, bitdiff, Fsgn;
  unsigned long int idim_max_buf, j;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;

  mask1 = (int)X->Def.Tpow[2 * org_isite1 + org_ispin1];
  mask2 = (int)X->Def.Tpow[2 * org_isite2 + org_ispin2];
  if (mask2 > mask1) bitdiff = mask2 - mask1 * 2;
  else bitdiff = mask1 - mask2 * 2;
  origin = myrank ^ (mask1 + mask2);

  state1 = origin & mask1;
  state2 = origin & mask2;

  SgnBit((unsigned long int) (origin & bitdiff), &Fsgn); // Fermion sign

  if (state1 == 0 && state2 == mask2) {
    trans = -(double)Fsgn * tmp_trans;
  }/*if (state1 == 0 && state2 == mask2)*/
  else if (state1 == mask1 && state2 == 0) {
    trans = -(double)Fsgn * conj(tmp_trans);
    if (X->Large.mode == M_CORR || X->Large.mode == M_CALCSPEC) trans = 0.0;
  }/*if (state1 == mask1 && state2 == 0)*/
  else return 0.0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0, 
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,       idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, 
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

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
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }/*if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC)*/
    else {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {
        dmv = trans * v1buf[j];
        dam_pr += conj(tmp_v1[j]) * dmv;
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }
  }/*End of parallel region*/
  return (dam_pr);
#else
  return 0.0;
#endif
}/*void GC_general_hopp_MPIdouble*/
/**
@brief Hopping term in Hubbard + MPI
When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
@return fragment of @f$\langle v_1|{\hat H}|v_1\rangle@f$
*/
double complex child_CisAjt_MPIdouble(
  int org_isite1,//!<[in] @f$i_1@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  int org_ispin1,//!<[in] @f$\sigma_1@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  int org_isite2,//!<[in] @f$i_2@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  int org_ispin2,//!<[in] @f$\sigma_2@f$ of @f$c_{i_1 \sigma_1}^\dagger c_{i_2 \sigma_2}@f$
  double complex tmp_trans,//!<[in] Transfer @f$t@f$
  struct BindStruct *X,//!< [inout]
  double complex *tmp_v0,//!< [out] Result v0 = H v1
  double complex *tmp_v1,//!< [in] v0 = H v1
  double complex *v1buf,//!<[in]
  long unsigned int *list_1_org,//!<[in]
  long unsigned int *list_1buf_org,//!<[in]
  long unsigned int *list_2_1_target,//!<[in]
  long unsigned int *list_2_2_target//!<[in]
) {
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin, bitdiff, Fsgn;
  unsigned long int idim_max_buf, j, ioff;
  MPI_Status statusMPI;
  double complex trans, dmv;

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
  }/*if (state1 == 0 && state2 == mask2)*/
  else if (state1 == mask1 && state2 == 0) {
    trans = -(double) Fsgn * conj(tmp_trans);
    if (X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) {
      trans = 0;
    }
  }/*if (state1 == mask1 && state2 == 0)*/
  else return 0;

  ierr = MPI_Sendrecv(&X->Check.idim_maxOrg, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf,         1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  ierr = MPI_Sendrecv(list_1_org, X->Check.idim_maxOrg + 1, MPI_UNSIGNED_LONG, origin, 0,
                      list_1buf_org,      idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_maxOrg + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,          idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  
  if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) private(j, dmv, ioff) \
  firstprivate(idim_max_buf, trans, X, list_2_1_target, list_2_2_target, list_1buf_org) \
  shared(v1buf, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++){
      dmv = trans * v1buf[j];
      GetOffComp(list_2_1_target, list_2_2_target, list_1buf_org[j],
                 X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
      tmp_v0[ioff] += dmv;
    }/*for (j = 1; j <= idim_max_buf; j++)*/
  }/*if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC)*/
  else return 0.0;
  
  return 1.0;
#else
  return 0.0;
#endif
}/*void CisAjt_MPIdouble*/
/**
@brief Hopping term in Hubbard + GC
 When only site2 is in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void GC_general_hopp_MPIsingle(
  unsigned long int itrans,//!<[in] Transfer ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
#ifdef MPI
  double complex dam_pr=0;
  dam_pr=child_GC_general_hopp_MPIsingle(
    X->Def.EDGeneralTransfer[itrans][0], X->Def.EDGeneralTransfer[itrans][1],
    X->Def.EDGeneralTransfer[itrans][2], X->Def.EDGeneralTransfer[itrans][3],
    X->Def.EDParaGeneralTransfer[itrans], X, tmp_v0, tmp_v1       );
  X->Large.prdct += dam_pr;
#endif
}/*void GC_general_hopp_MPIsingle*/
/**
@brief Hopping term in Hubbard + GC
 When only site2 is in the inter process region.
@return @f$\langle v_1|{\hat H}_{\rm this}|v_1\rangle@f$
@author Mitsuaki Kawamura (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex child_GC_general_hopp_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  double complex tmp_trans,//!<[in] Hopping integral
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
) {
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
                      &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,       idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  /*
    Index in the intra PE
  */
  mask1 = X->Def.Tpow[2 * org_isite1 + org_ispin1];

  if (state2 == mask2) {
    trans = -(double) Fsgn * tmp_trans;
    state1check = 0;
  }/*if (state2 == mask2)*/
  else if (state2 == 0) {
    state1check = mask1;
    trans = -(double) Fsgn * conj(tmp_trans);
    if (X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) trans = 0;
  }/*if (state2 != mask2)*/
  else return 0;

  bit1diff = X->Def.Tpow[2 * X->Def.Nsite - 1] * 2 - mask1 * 2;

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) private(j, dmv, state1, Fsgn, ioff) \
  firstprivate(idim_max_buf, trans, X, mask1, state1check, bit1diff) shared(v1buf, tmp_v1, tmp_v0)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 0; j < idim_max_buf; j++) {

        state1 = j & mask1;

        if (state1 == state1check) {

          SgnBit(j & bit1diff, &Fsgn);
          ioff = j ^ mask1;

          dmv = (double)Fsgn * trans * v1buf[j + 1];
          tmp_v0[ioff + 1] += dmv;
          dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
        }/*if (state1 == state1check)*/
      }/*for (j = 0; j < idim_max_buf; j++)*/
    }/*if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC)*/
    else {
#pragma omp for
      for (j = 0; j < idim_max_buf; j++) {

        state1 = j & mask1;

        if (state1 == state1check) {

          SgnBit(j & bit1diff, &Fsgn);
          ioff = j ^ mask1;

          dmv = (double)Fsgn * trans * v1buf[j + 1];
          dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
        }/*for (j = 0; j < idim_max_buf; j++)*/
      }/*for (j = 0; j < idim_max_buf; j++)*/
    }/*if (! (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC))*/
  }/*End of parallel region*/
  return (dam_pr);
#else
  return 0.0;
#endif
}/*void GC_general_hopp_MPIsingle*/
/**
@brief Hopping term in Hubbard (Kondo) + Canonical ensemble
 When both site1 and site2 are in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void general_hopp_MPIdouble(
  unsigned long int itrans,//!<[in] Transfer ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
#ifdef MPI
  double complex dam_pr;
  dam_pr =child_general_hopp_MPIdouble( 
    X->Def.EDGeneralTransfer[itrans][0], X->Def.EDGeneralTransfer[itrans][1],
    X->Def.EDGeneralTransfer[itrans][2], X->Def.EDGeneralTransfer[itrans][3],
    X->Def.EDParaGeneralTransfer[itrans], X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}/*void general_hopp_MPIdouble*/
/**
@brief Hopping term in Hubbard (Kondo) + Canonical ensemble
 When both site1 and site2 are in the inter process region.
@return @f$\langle v_1|{\hat H}_{\rm this}|v_1\rangle@f$
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double complex child_general_hopp_MPIdouble(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  double complex tmp_trans,//!<[in] Hopping integral
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
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
    if (X->Large.mode == M_CORR|| X->Large.mode == M_CALCSPEC) trans = 0;
  }
  else return 0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                      list_1buf,   idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,       idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, 
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);

  dam_pr = 0.0;
#pragma omp parallel default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff) \
  firstprivate(idim_max_buf, trans, X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {
        GetOffComp(list_2_1, list_2_2, list_1buf[j],
          X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
        dmv = trans * v1buf[j];
        tmp_v0[ioff] += dmv;
        dam_pr += conj(tmp_v1[ioff]) * dmv;
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }/*if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC)*/
    else {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {
        GetOffComp(list_2_1, list_2_2, list_1buf[j],
          X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
        dmv = trans * v1buf[j];
        dam_pr += conj(tmp_v1[ioff]) * dmv;
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }/*for (j = 1; j <= idim_max_buf; j++)*/
  }/*End of parallel region*/
  return (dam_pr);
#else
  return 0.0;
#endif
}/*void general_hopp_MPIdouble*/
/**
@brief Hopping term in Hubbard (Kondo) + Canonical ensemble
 When only site2 is in the inter process region.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void general_hopp_MPIsingle(
  unsigned long int itrans,//!<[in] Transfer ID
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
){
#ifdef MPI
  double complex dam_pr;
  dam_pr =child_general_hopp_MPIsingle(
    X->Def.EDGeneralTransfer[itrans][0], X->Def.EDGeneralTransfer[itrans][1],
    X->Def.EDGeneralTransfer[itrans][2], X->Def.EDGeneralTransfer[itrans][3],
    X->Def.EDParaGeneralTransfer[itrans], X, tmp_v0, tmp_v1);
  X->Large.prdct += dam_pr;
#endif
}/*void general_hopp_MPIsingle*/
/**
@brief Hopping term in Hubbard (Kondo) + Canonical ensemble
 When only site2 is in the inter process region.
@return @f$\langle v_1|{\hat H}_{\rm this}|v_1\rangle@f$
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double complex child_general_hopp_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  double complex tmp_trans,//!<[in] Hopping integral
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1//!<[in] v0 = H v1
) {
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
                      &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                      list_1buf,   idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,       idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, 
                      MPI_COMM_WORLD, &statusMPI);
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
#pragma omp parallel default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff, jreal, state1) \
  firstprivate(idim_max_buf, trans, X, mask1, state1check, bit1diff, myrank) shared(list_1, list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
  {
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {

        jreal = list_1buf[j];
        state1 = jreal & mask1;

        if (state1 == state1check) {
          SgnBit(jreal & bit1diff, &Fsgn);
          GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
            X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);

          dmv = (double)Fsgn * trans * v1buf[j];
          tmp_v0[ioff] += dmv;
          dam_pr += conj(tmp_v1[ioff]) * dmv;
        }/*if (state1 == state1check)*/
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }/*if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC)*/
    else {
#pragma omp for
      for (j = 1; j <= idim_max_buf; j++) {

        jreal = list_1buf[j];
        state1 = jreal & mask1;

        if (state1 == state1check) {
          SgnBit(jreal & bit1diff, &Fsgn);
          GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
            X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
/*
          if(X->Large.mode==M_CORR){
            printf("DEBUG: myrank=%d, org=%d, bit=%d, iexchg=%d, ioff=%d, list_1=%d\n", myrank, jreal, state1, jreal ^ mask1, ioff, list_1[ioff]);
          }
*/
          dmv = (double)Fsgn * trans * v1buf[j];
          dam_pr += conj(tmp_v1[ioff]) * dmv;
        }/*if (state1 == state1check)*/
      }/*for (j = 1; j <= idim_max_buf; j++)*/
    }/*if (! (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC))*/
  }/*End of parallel region*/
  return (dam_pr);
#else
  return 0.0;
#endif
}/*double complex general_hopp_MPIsingle*/
/**
@brief Hopping term in Hubbard (Kondo) + Canonical ensemble
  When only site2 is in the inter process region.
@return @f$\langle v_1|{\hat H}_{\rm this}|v_1\rangle@f$
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double complex child_CisAjt_MPIsingle(
  int org_isite1,//!<[in] Site 1
  int org_ispin1,//!<[in] Spin 1
  int org_isite2,//!<[in] Site 2
  int org_ispin2,//!<[in] Spin 2
  double complex tmp_trans,//!<[in] Hopping integral
  struct BindStruct *X,//!<[inout]
  double complex *tmp_v0,//!<[out] Result v0 = H v1
  double complex *tmp_v1,//!<[in] v0 = H v1
  double complex *v1buf,//!<[in] Buffer for sendrecv of wavefunction
  long unsigned int *list_1_org,//!<[in] Similler to ::list_1
  long unsigned int *list_1buf_org,//!<[in] Similler to ::list_1buf
  long unsigned int *list_2_1_target,//!<[in] ???
  long unsigned int *list_2_2_target//!<[in] ???
){
#ifdef MPI
  int mask2, state2, ierr, origin, bit2diff, Fsgn;
  unsigned long int mask1, state1, idim_max_buf, j, state1check, bit1diff, ioff, jreal;
  MPI_Status statusMPI;
  double complex trans, dmv;
  /*
    Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[2 * org_isite2+org_ispin2];
  bit2diff = mask2 - 1;
  origin = myrank ^ mask2;

  state2 = origin & mask2;

  SgnBit((unsigned long int) (origin & bit2diff), &Fsgn); // Fermion sign

  ierr = MPI_Sendrecv(&X->Check.idim_maxOrg, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &idim_max_buf,         1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(list_1_org, X->Check.idim_maxOrg + 1, MPI_UNSIGNED_LONG, origin, 0,
                      list_1buf_org,      idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(-1);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_maxOrg + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      v1buf,          idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
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
  }
  else return 0;

  bit1diff = X->Def.Tpow[2 * X->Def.Nsite - 1] * 2 - mask1 * 2;

  if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC) {
#pragma omp parallel for default(none) private(j, dmv, Fsgn, ioff, jreal, state1) \
  firstprivate(idim_max_buf, trans, X, mask1, state1check, bit1diff,list_2_1_target, list_2_2_target, list_1buf_org, list_1) shared(v1buf, tmp_v0)
    for (j = 1; j <= idim_max_buf; j++) {
      jreal = list_1buf_org[j];
      state1 = jreal & mask1;
      if (state1 == state1check) {
        SgnBit(jreal & bit1diff, &Fsgn);
        GetOffComp(list_2_1_target, list_2_2_target, jreal ^ mask1,
                   X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
        if(ioff !=0){
          dmv = (double) Fsgn * trans * v1buf[j];
          tmp_v0[ioff] += dmv;
        }/*if(ioff !=0)*/
      }/*if (state1 == state1check)*/
    }/*for (j = 1; j <= idim_max_buf; j++)*/
  }/*if (X->Large.mode == M_MLTPLY|| X->Large.mode == M_CALCSPEC)*/
  else return 0;

  return 1;
#else
  return 0.0;
#endif
}/*double complex general_hopp_MPIsingle*/
