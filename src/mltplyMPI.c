/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima */

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
void GC_child_general_hopp_MPIdouble(
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
				       double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
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

  SgnBit((unsigned long int)(origin & bitdiff), &Fsgn); // Fermion sign

  if(state1 == 0 && state2 == mask2){
    trans = - (double)Fsgn * tmp_trans;
  }
  else if(state1 == mask1 && state2 == 0) {
    trans = - (double)Fsgn * conj(tmp_trans);
    if(X->Large.mode == M_CORR){
      trans = 0;
    }
  }
  else return 0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

  
  dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) firstprivate(idim_max_buf, trans, X) shared(v1buf, tmp_v1, tmp_v0)
  for (j = 1; j <= idim_max_buf; j++) {
    dmv = trans * v1buf[j];
    if (X->Large.mode == M_MLTPLY) tmp_v0[j] += dmv;
    dam_pr += conj(tmp_v1[j]) * dmv;
  }
  return (dam_pr);
  
#endif
}/*void GC_child_general_hopp_MPIdouble*/


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
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  int mask2, state1, state2, ierr, origin, bit2diff, Fsgn;
  unsigned long int idim_max_buf, j, mask1, state1check, bit1diff, ioff;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;
  /*
   Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[2 * org_isite2+org_ispin2];
  bit2diff = mask2 - 1;
  origin = myrank ^ mask2;
  state2 = origin & mask2;

  SgnBit((unsigned long int)(origin & bit2diff), &Fsgn); // Fermion sign

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

  /*
   Index in the intra PE
  */
  mask1 = X->Def.Tpow[2 * org_isite1+ org_ispin1];
  
  if (state2 == mask2) {
    trans = - (double)Fsgn * tmp_trans;
    state1check = 0;
  }
  else if (state2 == 0) {
    state1check = mask1;
    trans = - (double)Fsgn * conj(tmp_trans);
    if(X->Large.mode == M_CORR){
      trans = 0;
    }

  }
  else return 0;

  bit1diff = X->Def.Tpow[2 * X->Def.Nsite - 1] * 2 - mask1 * 2;

  dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, state1, Fsgn, ioff) \
  firstprivate(idim_max_buf, trans, X, mask1, state1check, bit1diff) shared(v1buf, tmp_v1, tmp_v0)
  for (j = 0; j < idim_max_buf; j++) {

    state1 = j & mask1;

    if (state1 == state1check) {

      SgnBit(j & bit1diff, &Fsgn);
      ioff = j ^ mask1;

      dmv = (double)Fsgn * trans * v1buf[j + 1];
      if (X->Large.mode == M_MLTPLY) tmp_v0[ioff + 1] += dmv;
      dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
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
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin, bitdiff, Fsgn;
  unsigned long int idim_max_buf, j, ioff;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;

  mask1 = (int)X->Def.Tpow[2 * org_isite1+org_ispin1];
  mask2 = (int)X->Def.Tpow[2 * org_isite2+org_ispin2];
  if(mask2 > mask1) bitdiff = mask2 - mask1 * 2;
  else bitdiff = mask1 - mask2 * 2;
  origin = myrank ^ (mask1 + mask2);

  state1 = origin & mask1;
  state2 = origin & mask2;

   SgnBit((unsigned long int)(origin & bitdiff), &Fsgn); // Fermion sign

  if (state1 == 0 && state2 == mask2) {
    trans = - (double)Fsgn * tmp_trans;
  }
  else if (state1 == mask1 && state2 == 0) {
    trans = - (double)Fsgn * conj(tmp_trans);
    if(X->Large.mode == M_CORR){
      trans = 0;
    }
  }
  else return 0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
    list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

  dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff) \
  firstprivate(idim_max_buf, trans, X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
  for (j = 1; j <= idim_max_buf; j++) {
    GetOffComp(list_2_1, list_2_2, list_1buf[j], 
      X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
    dmv = trans * v1buf[j];
    if (X->Large.mode == M_MLTPLY) tmp_v0[ioff] += dmv;
    dam_pr += conj(tmp_v1[ioff]) * dmv;
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
  double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin, bit2diff, Fsgn;
  unsigned long int idim_max_buf, j, state1check, bit1diff, ioff, jreal;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[2 * org_isite2+org_ispin2];
  bit2diff = mask2 - 1;
  origin = myrank ^ mask2;

  state2 = origin & mask2;

  SgnBit((unsigned long int)(origin & bit2diff), &Fsgn); // Fermion sign

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
    list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

  /*
  Index in the intra PE
  */
  mask1 = X->Def.Tpow[2 * org_isite1+ org_ispin1];
  if (state2 == mask2) {
    trans = - (double)Fsgn * tmp_trans;
    state1check = 0;
  }
  else if (state2 == 0) {
    state1check = mask1;
    trans = - (double)Fsgn * conj(tmp_trans);
    if(X->Large.mode == M_CORR){
      trans = 0;
    }
  }
  else return 0;

  bit1diff = X->Def.Tpow[2 * X->Def.Nsite - 1] * 2 - mask1 * 2;

  dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, Fsgn, ioff, jreal, state1) \
  firstprivate(idim_max_buf, trans, X, mask1, state1check, bit1diff) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
  for (j = 1; j <= idim_max_buf; j++) {

    jreal = list_1buf[j];
    state1 = jreal & mask1;

    if (state1 == state1check) {

      SgnBit(jreal & bit1diff, &Fsgn);
      GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
        X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);

      dmv = (double)Fsgn * trans * v1buf[j];
      if (X->Large.mode == M_MLTPLY) tmp_v0[ioff] += dmv;
      dam_pr += conj(tmp_v1[ioff]) * dmv;
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
  dam_pr=X_child_general_int_spin_MPIdouble(
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
						  )
{
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin;
  unsigned long int idim_max_buf, j, ioff;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;

  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  origin = myrank ^ (mask1 + mask2);

  state1 = (origin & mask1) / mask1;
  state2 = (origin & mask2) / mask2;

  if (state1 == org_ispin2 &&  state2 == org_ispin4) {
    Jint = tmp_J;
  }
  else if (state1 == org_ispin1 && state2 == org_ispin3) {
    Jint = conj(tmp_J);
    if(X->Large.mode == M_CORR){
      Jint = 0;
    }
  }
  else return 0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
    list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

  dam_pr = 0.0;
  #pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff) \
  firstprivate(idim_max_buf, Jint, X) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
  for (j = 1; j <= idim_max_buf; j++) {
    GetOffComp(list_2_1, list_2_2, list_1buf[j],
	       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);
    dmv = Jint * v1buf[j];
    if (X->Large.mode == M_MLTPLY) tmp_v0[ioff] += dmv;
    dam_pr += conj(tmp_v1[ioff]) * dmv;
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
  int mask1, mask2, num1_up, num2_up, num1_down, num2_down, ierr, origin;
  unsigned long int idim_max_buf, j, ioff, ibit_tmp;
  MPI_Status statusMPI;
  double complex dmv, dam_pr;

  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  origin = myrank ^ (mask1 + mask2);

  num1_up = (origin & mask1) / mask1;
  num1_down = 1- num1_up;
  num2_up = (origin & mask2) / mask2;
  num2_down = 1- num2_up;

  ibit_tmp=(num1_up)^(num2_up);
  if(ibit_tmp ==0) return 0;
  
  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
    list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

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
						  )
{
#ifdef MPI
  int mask2, state2, ierr, origin;
  int num1_up, num1_down, num2_up,num2_down;
  unsigned long int is1_up, ibit1_up;
  unsigned long int mask1, idim_max_buf, j, ioff, state1, jreal, state1check;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr, spn_z;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[org_isite3];
  origin = myrank ^ mask2;
  state2 = (origin & mask2) / mask2;

  if (state2 == org_ispin4) {
    state1check = (unsigned long int)org_ispin2;
    Jint = tmp_J;
  }
  else if (state2 == org_ispin3) {
    state1check = (unsigned long int) org_ispin1;
    Jint = conj(tmp_J);
    if(X->Large.mode == M_CORR){
      Jint = 0;
    }   
  }
  else return 0;

  if(X->Large.mode==M_TOTALS){
    num2_up= X_SpinGC_CisAis((unsigned long int)myrank + 1, X, mask2, 0);
    num2_down = 1-num2_up;
  }
  
  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
    list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
  /*
  Index in the intra PE
  */
  mask1 = X->Def.Tpow[org_isite1];

  dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, ioff, jreal, state1, num1_up, num1_down, is1_up, ibit1_up, spn_z) \
  firstprivate(idim_max_buf, Jint, X, mask1, state1check, num2_up, num2_down, org_isite1) shared(list_2_1, list_2_2, list_1buf, v1buf, tmp_v1, tmp_v0)
  for (j = 1; j <= idim_max_buf; j++) {

    jreal = list_1buf[j];

    state1 = (jreal & mask1) / mask1;
    if (state1 == state1check) {
      GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
		 X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);

      dmv = Jint * v1buf[j];
      if (X->Large.mode == M_MLTPLY) tmp_v0[ioff] += dmv;
      else if(X->Large.mode==M_TOTALS){
	dmv=0.5*v1buf[j];
      }
      dam_pr += conj(tmp_v1[ioff]) * dmv;
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
  double complex *tmp_v0, double complex *tmp_v1
						      )
{
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin;
  unsigned long int idim_max_buf, j;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;

  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  origin = myrank ^ (mask1 + mask2);

  state1 = (origin & mask1) / mask1;
  state2 = (origin & mask2) / mask2;

  if (state1 == org_ispin2 && state2 == org_ispin4) {
    Jint = tmp_J;
  }
  else if (state1 == org_ispin1 && state2 == org_ispin3) {
    Jint = conj(tmp_J);
    if(X->Large.mode == M_CORR){
      Jint = 0;
    }
  }
  else return 0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

  dam_pr = 0.0;
  #pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
  for (j = 1; j <= idim_max_buf; j++) {
    dmv = Jint * v1buf[j];
    if (X->Large.mode == M_MLTPLY) tmp_v0[j] += dmv;
    dam_pr += conj(tmp_v1[j]) * dmv;
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
					    double complex *tmp_v1)
{
#ifdef MPI
  int mask1, mask2, state1, state2, ierr;
  long int origin;
  unsigned long int idim_max_buf, j;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr,  tmp_off;
  int tmp_sgn;
  
  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  origin = myrank ^ mask2;
  state2 = (origin & mask2) / mask2;
  if(state2 == org_ispin4){
    Jint = tmp_J;
  }
  else if(state2 == org_ispin3){
    Jint = conj(tmp_J);
    if(X->Large.mode == M_CORR){
      Jint = 0;
    }
  }
  else return 0.0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

  dam_pr = 0.0;
  #pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
  for (j = 1; j <= idim_max_buf; j++) {
    dmv = Jint * v1buf[j];
    if (X->Large.mode == M_MLTPLY) tmp_v0[j] += dmv;
    dam_pr += conj(tmp_v1[j]) * dmv;
  }
  return(dam_pr);
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
  dam_pr = X_GC_child_CisAitCjuAju_spin_MPIdouble( X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1],X->Def.InterAll_OffDiagonal[i_int][2],
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
					    double complex *tmp_v1)
{
#ifdef MPI
  int mask1, mask2, state1, state2, ierr;
  long int origin;
  unsigned long int idim_max_buf, j;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr,  tmp_off;
  
  mask1 = (int)X->Def.Tpow[org_isite1];
  mask2 = (int)X->Def.Tpow[org_isite3];
  origin = myrank ^ mask1;
  state1 = (origin & mask1)/mask1;

  if(state1 == org_ispin2){
    Jint = tmp_J;
  }
  else if(state1 == org_ispin1){
    Jint = conj(tmp_J);
     if(X->Large.mode == M_CORR){
      Jint = 0;
    }
  }
  else return 0.0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

  dam_pr = 0.0;
  #pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, Jint, X) shared(v1buf, tmp_v1, tmp_v0)
  for (j = 1; j <= idim_max_buf; j++) {
    dmv = Jint * v1buf[j];
    if (X->Large.mode == M_MLTPLY) tmp_v0[j] += dmv;
    dam_pr += conj(tmp_v1[j]) * dmv;
  }
  return(dam_pr);
#endif
}/*double complex X_GC_child_CisAisCjuAjv_spin_MPIdouble*/

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
					    double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1)
{
#ifdef MPI
  int mask2, state2, ierr, origin;
  unsigned long int mask1, idim_max_buf, j, ioff, state1, state1check;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[org_isite3];
  origin = myrank ^ mask2;
  state2 = (origin & mask2) / mask2;

  if (state2 == org_ispin4) {
    state1check = (unsigned long int)org_ispin2;
    Jint = tmp_J;
  }
  else if (state2 == org_ispin3) {
    state1check = (unsigned long int)org_ispin1;
    Jint = conj(tmp_J);
    if(X->Large.mode == M_CORR){
      Jint = 0;
    }
  }
  else return 0.0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
  /*
  Index in the intra PE
  */
  mask1 = X->Def.Tpow[org_isite1];

  dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, state1, ioff) \
    firstprivate(idim_max_buf, Jint, X, state1check, mask1) shared(v1buf, tmp_v1, tmp_v0)
  for (j = 0; j < idim_max_buf; j++) {

    state1 = (j & mask1) / mask1;
    if (state1 == state1check) {

      ioff = j ^ mask1;

      dmv = Jint * v1buf[j + 1];
      if (X->Large.mode == M_MLTPLY) tmp_v0[ioff + 1] += dmv;
      dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
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
double complex X_GC_child_CisAisCjuAjv_spin_MPIsingle( int org_isite1, int org_ispin1,  int org_isite3, int org_ispin3, int org_ispin4, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1)
{
#ifdef MPI
  int mask2, state2, ierr, origin;
  unsigned long int mask1, idim_max_buf, j, ioff, state1, state1check;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[org_isite3];
  origin = myrank ^ mask2;
  state2 = (origin & mask2) / mask2;
  if (state2 == org_ispin4) {
    state1check = (unsigned long int) org_ispin1;
    Jint = tmp_J;
  }
  else if (state2 == org_ispin3) {
    state1check = (unsigned long int)org_ispin1;
    Jint = conj(tmp_J);
    if(X->Large.mode == M_CORR){
      Jint = 0;
    }
  }
  else return 0.0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
  /*
  Index in the intra PE
  */
  mask1 = X->Def.Tpow[org_isite1];

  dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, state1, ioff) \
    firstprivate(idim_max_buf, Jint, X, state1check, mask1) shared(v1buf, tmp_v1, tmp_v0)
  for (j = 0; j < idim_max_buf; j++) {
    state1 = (j & mask1) / mask1;
    if (state1 == state1check) {
      dmv = Jint * v1buf[j + 1];
      if (X->Large.mode == M_MLTPLY) tmp_v0[j + 1] += dmv;
      dam_pr += conj(tmp_v1[j + 1]) * dmv;
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
  dam_pr =X_GC_child_CisAitCjuAju_spin_MPIsingle(X->Def.InterAll_OffDiagonal[i_int][0], X->Def.InterAll_OffDiagonal[i_int][1], X->Def.InterAll_OffDiagonal[i_int][2], X->Def.InterAll_OffDiagonal[i_int][4], X->Def.InterAll_OffDiagonal[i_int][5], X->Def.ParaInterAll_OffDiagonal[i_int], X, tmp_v0, tmp_v1);
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
double complex X_GC_child_CisAitCjuAju_spin_MPIsingle( int org_isite1, int org_ispin1, int org_ispin2,  int org_isite3, int org_ispin3, double complex tmp_J, struct BindStruct *X, double complex *tmp_v0, double complex *tmp_v1)
{
#ifdef MPI
  int mask2, state2, ierr, origin;
  unsigned long int mask1, idim_max_buf, j, ioff, state1, state1check;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[org_isite3];
  state2 = (origin & mask2) / mask2;

  if (state2 == org_ispin3) {
    state1check = org_ispin2;
    Jint = tmp_J;
  }
  else return 0.0;
    
  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
  /*
  Index in the intra PE
  */
  mask1 = (int)X->Def.Tpow[org_isite1];

  dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv, state1, ioff) \
  firstprivate(idim_max_buf, Jint, X, state1check, mask1) shared(v1buf, tmp_v1, tmp_v0)
  for (j = 0; j < idim_max_buf; j++) {
 
    state1 = (j & mask1) / mask1;
    if (state1 == state1check) {
      ioff = j ^ mask1;
      dmv = Jint * v1buf[j + 1];
      if (X->Large.mode == M_MLTPLY) tmp_v0[ioff + 1] += dmv;
      dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
    }
    else{
      ioff = j ^ mask1;
      dmv = conj(Jint) * v1buf[j + 1];
      if(X->Large.mode == M_CORR) dmv=0.0;
      if (X->Large.mode == M_MLTPLY) tmp_v0[ioff + 1] += dmv;
      
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
  int mask1, mask2, state1, state2, ierr, origin, bitdiff, Fsgn;
  unsigned long int idim_max_buf, j;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;
  
  mask1 = (int)X->Def.Tpow[org_isite1];
  origin = myrank ^ mask1;
  state1 = (origin & mask1)/mask1;
  
  if(state1 ==  org_ispin1){
    trans = tmp_trans;
  }
  else if(state1 == org_ispin2) {
    trans = conj(tmp_trans);
    if(X->Large.mode == M_CORR){
      trans = 0;
    }
  }
  else return 0.0;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

  dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) firstprivate(idim_max_buf, trans, X) shared(v1buf, tmp_v1, tmp_v0)
  for (j = 1; j <= idim_max_buf; j++) {
    dmv = trans * v1buf[j];
    if (X->Large.mode == M_MLTPLY) tmp_v0[j] += dmv;
    dam_pr += conj(tmp_v1[j]) * dmv;
  }
  return (dam_pr);
  
#endif
}/*double complex  X_GC_child_CisAit_spin_MPIdouble*/

/**
 *
 * Hopping term in Spin + GC
 * When both site1 and site2 are in the inter process region.
 *
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
double complex X_GC_child_CisAis_spin_MPIdouble(
				       int org_isite1,
				       int org_ispin1,
				       double complex tmp_trans,
				       struct BindStruct *X /**< [inout]*/,
				       double complex *tmp_v0 /**< [out] Result v0 = H v1*/, 
				       double complex *tmp_v1 /**< [in] v0 = H v1*/)
{
#ifdef MPI
  long unsigned int j;
  int mask1, state1, state2, ierr, origin, bitdiff, Fsgn;
  int ibit1;
  double complex trans, dam_pr;
  mask1 = (int)X->Def.Tpow[org_isite1];
  ibit1 = ((unsigned long int)myrank& mask1)^(1-org_ispin1);

  dam_pr =0.0;
  if(ibit1 != 0){
#pragma omp parallel for reduction(+:dam_pr)default(none) shared(tmp_v1, tmp_v0) \
  firstprivate(X, tmp_trans) private(j)
    for (j = 1; j <= X->Check.idim_max; j++){
      if (X->Large.mode == M_MLTPLY) { // for multply
	tmp_v0[j] += tmp_v1[j]*tmp_trans;
      }
      dam_pr += tmp_trans*conj(tmp_v1[j])*tmp_v1[j];
    }
  }
  return dam_pr;
#endif
}
