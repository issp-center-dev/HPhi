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

void GC_child_general_hopp_MPIdouble(unsigned long int itrans, struct BindStruct *X, 
  double complex *tmp_v0, double complex *tmp_v1)
{
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin, bitdiff, Fsgn;
  unsigned long int idim_max_buf, j;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;

  mask1 = (int)X->Def.Tpow[2 * X->Def.EDGeneralTransfer[itrans][0] 
                            + X->Def.EDGeneralTransfer[itrans][1]];
  mask2 = (int)X->Def.Tpow[2 * X->Def.EDGeneralTransfer[itrans][2]
                            + X->Def.EDGeneralTransfer[itrans][3]];
  if (mask2 > mask1) bitdiff = mask2 - mask1 * 2;
  else bitdiff = mask1 - mask2 * 2;
  origin = myrank ^ (mask1 + mask2);

  state1 = origin & mask1;
  state2 = origin & mask2;

  SgnBit((unsigned long int)(origin & bitdiff), &Fsgn); // Fermion sign

  if(state1 == 0 && state2 == mask2){
    trans = - (double)Fsgn * X->Def.EDParaGeneralTransfer[itrans];
  }
  else if(state1 == mask1 && state2 == 0) {
    trans = - (double)Fsgn * conj(X->Def.EDParaGeneralTransfer[itrans]);
  }
  else return;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

  dam_pr = 0.0;
  for (j = 1; j <= idim_max_buf; j++) {
    dmv = trans * v1buf[j];
    if (X->Large.mode == M_MLTPLY) tmp_v0[j] += dmv;
    dam_pr += conj(tmp_v1[j]) * dmv;
  }

  X->Large.prdct += dam_pr;
  
#endif
}

void GC_child_general_hopp_MPIsingle(unsigned long int itrans, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1)
{
#ifdef MPI
  int mask2, state1, state2, ierr, origin, bit2diff, Fsgn;
  unsigned long int idim_max_buf, j, mask1, state1check, bit1diff, ioff;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;
  /*
   Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[2 * X->Def.EDGeneralTransfer[itrans][2]
    + X->Def.EDGeneralTransfer[itrans][3]];
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
  mask1 = X->Def.Tpow[2 * X->Def.EDGeneralTransfer[itrans][0]
    + X->Def.EDGeneralTransfer[itrans][1]];
  
  if (state2 == mask2) {
    trans = - (double)Fsgn * X->Def.EDParaGeneralTransfer[itrans];
    state1check = 0;
  }
  else if (state2 == 0) {
    trans = - (double)Fsgn * conj(X->Def.EDParaGeneralTransfer[itrans]);
    state1check = mask1;
  }
  else return;

  bit1diff = X->Def.Tpow[2 * X->Def.Nsite - 1] * 2 - mask1 * 2;

  dam_pr = 0.0;
  for (j = 0; j < idim_max_buf; j++) {

    state1 = j & mask1;

    if (state1 == state1check) {

      SgnBit(j & bit1diff, &Fsgn);
      ioff = j ^ mask1;

      dmv = (double)Fsgn * trans * v1buf[j + 1];
      if (X->Large.mode == M_MLTPLY) tmp_v0[ioff + 1] += dmv;
      dam_pr += conj(tmp_v1[j + 1]) * dmv;
    }
  }

  X->Large.prdct += dam_pr;

#endif
}

void child_general_hopp_MPIdouble(unsigned long int itrans, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1)
{
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin, bitdiff, Fsgn;
  unsigned long int idim_max_buf, j, ioff;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;

  mask1 = (int)X->Def.Tpow[2 * X->Def.EDGeneralTransfer[itrans][0]
    + X->Def.EDGeneralTransfer[itrans][1]];
  mask2 = (int)X->Def.Tpow[2 * X->Def.EDGeneralTransfer[itrans][2]
    + X->Def.EDGeneralTransfer[itrans][3]];
  if(mask2 > mask1) bitdiff = mask2 - mask1 * 2;
  else bitdiff = mask1 - mask2 * 2;
  origin = myrank ^ (mask1 + mask2);

  state1 = origin & mask1;
  state2 = origin & mask2;

   SgnBit((unsigned long int)(origin & bitdiff), &Fsgn); // Fermion sign

  if (state1 == 0 && state2 == mask2) {
    trans = - (double)Fsgn * X->Def.EDParaGeneralTransfer[itrans];
  }
  else if (state1 == mask1 && state2 == 0) {
    trans = - (double)Fsgn * conj(X->Def.EDParaGeneralTransfer[itrans]);
  }
  else return;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
    list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

  dam_pr = 0.0;
  for (j = 1; j <= idim_max_buf; j++) {

    GetOffComp(list_2_1, list_2_2, list_1buf[j], 
      X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);

    dmv = trans * v1buf[j];
    if (X->Large.mode == M_MLTPLY) tmp_v0[ioff] += dmv;
    dam_pr += conj(tmp_v1[ioff]) * dmv;
  }

  X->Large.prdct += dam_pr;

#endif
}

void child_general_hopp_MPIsingle(unsigned long int itrans, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1)
{
#ifdef MPI
  int mask2, state1, state2, ierr, origin, bit2diff, Fsgn;
  unsigned long int idim_max_buf, j, mask1, state1check, bit1diff, ioff, jreal;
  MPI_Status statusMPI;
  double complex trans, dmv, dam_pr;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[2 * X->Def.EDGeneralTransfer[itrans][2]
    + X->Def.EDGeneralTransfer[itrans][3]];
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
  mask1 = X->Def.Tpow[2 * X->Def.EDGeneralTransfer[itrans][0]
    + X->Def.EDGeneralTransfer[itrans][1]];

  if (state2 == mask2) {
    trans = - (double)Fsgn * X->Def.EDParaGeneralTransfer[itrans];
    state1check = 0;
  }
  else if (state2 == 0) {
    trans = - (double)Fsgn * conj(X->Def.EDParaGeneralTransfer[itrans]);
    state1check = mask1;
  }
  else return;

  bit1diff = X->Def.Tpow[2 * X->Def.Nsite - 1] * 2 - mask1 * 2;

  dam_pr = 0.0;
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

  X->Large.prdct += dam_pr;

#endif
}

void child_general_int_spin_MPIdouble(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1)
{
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin;
  unsigned long int idim_max_buf, j, ioff;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;

  mask1 = (int)X->Def.Tpow[X->Def.InterAll_OffDiagonal[i_int][0]];
  mask2 = (int)X->Def.Tpow[X->Def.InterAll_OffDiagonal[i_int][4]];
  origin = myrank ^ (mask1 + mask2);

  state1 = (origin & mask1) / mask1;
  state2 = (origin & mask2) / mask2;

  if (state1 == X->Def.InterAll_OffDiagonal[i_int][3] && 
    state2 == X->Def.InterAll_OffDiagonal[i_int][7]) {
    Jint = X->Def.ParaInterAll_OffDiagonal[i_int];
  }
  else if (state1 == X->Def.InterAll_OffDiagonal[i_int][1] && 
    state2 == X->Def.InterAll_OffDiagonal[i_int][5]) {
    Jint = conj(X->Def.ParaInterAll_OffDiagonal[i_int]);
  }
  else return;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
    list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

  dam_pr = 0.0;
  for (j = 1; j <= idim_max_buf; j++) {

    GetOffComp(list_2_1, list_2_2, list_1buf[j],
      X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);

    dmv = Jint * v1buf[j];
    if (X->Large.mode == M_MLTPLY) tmp_v0[ioff] += dmv;
    dam_pr += conj(tmp_v1[ioff]) * dmv;
  }

  X->Large.prdct += dam_pr;

#endif
}

void child_general_int_spin_MPIsingle(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1)
{
#ifdef MPI
  int mask2, state2, ierr, origin;
  unsigned long int mask1, idim_max_buf, j, ioff, state1, jreal, state1check;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[X->Def.InterAll_OffDiagonal[i_int][4]];
  origin = myrank ^ mask2;
  state2 = (origin & mask2) / mask2;

  if (state2 == X->Def.InterAll_OffDiagonal[i_int][7]) {
    state1check = (unsigned long int)X->Def.InterAll_OffDiagonal[i_int][3];
    Jint = X->Def.ParaInterAll_OffDiagonal[i_int];
  }
  else if (state2 == X->Def.InterAll_OffDiagonal[i_int][5]) {
    state1check = (unsigned long int)X->Def.InterAll_OffDiagonal[i_int][1];
    Jint = conj(X->Def.ParaInterAll_OffDiagonal[i_int]);
  }
  else return;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
    list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
  /*
  Index in the intra PE
  */
  mask1 = X->Def.Tpow[X->Def.InterAll_OffDiagonal[i_int][0]];

  dam_pr = 0.0;
  for (j = 1; j <= idim_max_buf; j++) {

    jreal = list_1buf[j];

    state1 = (jreal & mask1) / mask1;
    if (state1 == state1check) {
      GetOffComp(list_2_1, list_2_2, jreal ^ mask1,
        X->Large.irght, X->Large.ilft, X->Large.ihfbit, &ioff);

      dmv = Jint * v1buf[j];
      if (X->Large.mode == M_MLTPLY) tmp_v0[ioff] += dmv;
      dam_pr += conj(tmp_v1[ioff]) * dmv;
    }
  }

  X->Large.prdct += dam_pr;

#endif
}

void GC_child_CisAitCiuAiv_spin_MPIdouble(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1)
{
#ifdef MPI
  int mask1, mask2, state1, state2, ierr, origin;
  unsigned long int idim_max_buf, j;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;

  mask1 = (int)X->Def.Tpow[X->Def.InterAll_OffDiagonal[i_int][0]];
  mask2 = (int)X->Def.Tpow[X->Def.InterAll_OffDiagonal[i_int][4]];
  origin = myrank ^ (mask1 + mask2);

  state1 = (origin & mask1) / mask1;
  state2 = (origin & mask2) / mask2;

  if (state1 == X->Def.InterAll_OffDiagonal[i_int][3] &&
    state2 == X->Def.InterAll_OffDiagonal[i_int][7]) {
    Jint = X->Def.ParaInterAll_OffDiagonal[i_int];
  }
  else if (state1 == X->Def.InterAll_OffDiagonal[i_int][1] &&
    state2 == X->Def.InterAll_OffDiagonal[i_int][5]) {
    Jint = conj(X->Def.ParaInterAll_OffDiagonal[i_int]);
  }
  else return;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);

  dam_pr = 0.0;
  for (j = 1; j <= idim_max_buf; j++) {
    dmv = Jint * v1buf[j];
    if (X->Large.mode == M_MLTPLY) tmp_v0[j] += dmv;
    dam_pr += conj(tmp_v1[j]) * dmv;
  }

  X->Large.prdct += dam_pr;

#endif
}

void GC_child_general_int_spin_MPIdouble(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1)
{
  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
      X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) { //diagonal
    fprintf(stderr, "\nThis interaction has not been supported yet.\n");
    exitMPI(-1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
    fprintf(stderr, "\nThis interaction has not been supported yet.\n");
    exitMPI(-1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
    fprintf(stderr, "\nThis interaction has not been supported yet.\n");
    exitMPI(-1);
  }
  else {
    GC_child_CisAitCiuAiv_spin_MPIdouble(i_int, X, tmp_v0, tmp_v1);
  }
}

void GC_child_CisAitCiuAiv_spin_MPIsingle(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1)
{
#ifdef MPI
  int mask2, state2, ierr, origin;
  unsigned long int mask1, idim_max_buf, j, ioff, state1, state1check;
  MPI_Status statusMPI;
  double complex Jint, dmv, dam_pr;
  /*
  Prepare index in the inter PE
  */
  mask2 = (int)X->Def.Tpow[X->Def.InterAll_OffDiagonal[i_int][4]];
  origin = myrank ^ mask2;
  state2 = (origin & mask2) / mask2;

  if (state2 == X->Def.InterAll_OffDiagonal[i_int][7]) {
    state1check = (unsigned long int)X->Def.InterAll_OffDiagonal[i_int][3];
    Jint = X->Def.ParaInterAll_OffDiagonal[i_int];
  }
  else if (state2 == X->Def.InterAll_OffDiagonal[i_int][5]) {
    state1check = (unsigned long int)X->Def.InterAll_OffDiagonal[i_int][1];
    Jint = conj(X->Def.ParaInterAll_OffDiagonal[i_int]);
  }
  else return;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
    v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
  /*
  Index in the intra PE
  */
  mask1 = X->Def.Tpow[X->Def.InterAll_OffDiagonal[i_int][0]];

  dam_pr = 0.0;
  for (j = 0; j < idim_max_buf; j++) {

    state1 = (j & mask1) / mask1;
    if (state1 == state1check) {

      ioff = j ^ mask1;

      dmv = Jint * v1buf[j + 1];
      if (X->Large.mode == M_MLTPLY) tmp_v0[ioff + 1] += dmv;
      dam_pr += conj(tmp_v1[ioff + 1]) * dmv;
    }
  }

  X->Large.prdct += dam_pr;

#endif
}

void GC_child_general_int_spin_MPIsingle(unsigned long int i_int, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1)
{
  if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
      X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) { //diagonal
    fprintf(stderr, "\nThis interaction has not been supported yet.\n");
    exitMPI(-1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] == X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] != X->Def.InterAll_OffDiagonal[i_int][7]) {
    fprintf(stderr, "\nThis interaction has not been supported yet.\n");
    exitMPI(-1);
  }
  else if (X->Def.InterAll_OffDiagonal[i_int][1] != X->Def.InterAll_OffDiagonal[i_int][3] &&
           X->Def.InterAll_OffDiagonal[i_int][5] == X->Def.InterAll_OffDiagonal[i_int][7]) {
    fprintf(stderr, "\nThis interaction has not been supported yet.\n");
    exitMPI(-1);
  }
  else {
    GC_child_CisAitCiuAiv_spin_MPIsingle(i_int, X, tmp_v0, tmp_v1);
  }
}


