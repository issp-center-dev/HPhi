/*
HPhi  -  Quantum Lattice Model Simulator
Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "StdFace_vals.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "StdFace_ModelUtil.h"

/**
 *
 * Setup a Hamiltonian for the Hubbard model on a Chain lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void FermionHubbard_ChainLattice(
struct StdIntList *StdI,
  int nelec /**< [in] The number of electrons */, 
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int ispin;
  int iL, iL2;
  int ktrans, kintr;
  double trans0, LargeValue0;
  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_NotUsed_i("W", StdI->W);
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  /**/
  StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
  StdFace_PrintVal_d("U", &StdI->U, 0.0);
  StdFace_PrintVal_d("t", &StdI->t, 1.0);
  StdFace_PrintVal_d("V", &StdI->V, 0.0);
  StdFace_PrintVal_d("t'", &StdI->tp, 0.0);
  StdFace_PrintVal_d("V'", &StdI->Vp, 0.0);
  /**/
  StdFace_NotUsed_i("2S", StdI->S2);
  StdFace_NotUsed_d("tpp", StdI->tpp);
  StdFace_NotUsed_d("t0", StdI->t0);
  StdFace_NotUsed_d("t1", StdI->t1);
  StdFace_NotUsed_d("t2", StdI->t2);
  StdFace_NotUsed_d("Vpp", StdI->Vpp);
  StdFace_NotUsed_d("V0", StdI->V0);
  StdFace_NotUsed_d("V1", StdI->V1);
  StdFace_NotUsed_d("V2", StdI->V2);
  /*
  Local Spin
  */
  StdI->nsite = StdI->L;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = 0;
  /*
  Transfer
  */
  StdI->ntrans = StdI->nsite * 2 * 5;
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double *)malloc(sizeof(double) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    StdI->transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  StdI->ntrans = 0;
  for (iL = 0; iL < StdI->L; iL++){
    isite = iL;
    for (ispin = 0; ispin < 2; ispin++){
      for (iL2 = -2; iL2 <= 2; iL2++){

        if (iL2 == 0) trans0 = StdI->mu;
        else if(abs(iL2) == 1) trans0 = StdI->t;
        else trans0 = StdI->tp;

        jsite = (iL + iL2 + 2 * StdI->L) % StdI->L;
        StdFace_trans(StdI, trans0, isite, ispin, jsite, ispin);
      }
    }
  }
  /*
  Interaction
  */
  StdI->nintr = StdI->nsite * (1 + 4 * 2);
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++){
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  StdI->nintr = 0;
  for (iL = 0; iL < StdI->L; iL++){

    isite = iL;

    StdFace_intr(StdI, StdI->U, isite, 0, isite, 0, isite, 1, isite, 1);

    iL2 = (iL + 1) % StdI->L;
    jsite = iL2;
    StdFace_Coulomb(StdI, StdI->V, isite, jsite);

    iL2 = (iL + 2) % StdI->L;
    jsite = iL2;
    StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
  }
  /*
  Set mTPQ parameter
  */
  if (lGC == 0){
    LargeValue0 = fabs(StdI->mu) * (double)nelec / (double)(StdI->L * StdI->W)
      + 2.0 * 2.0 * fabs(StdI->t) + 2.0 * 2.0 * fabs(StdI->tp)
      + fabs(StdI->U) + 2.0 * 2.0 * fabs(StdI->V) + 2.0 * 2.0 * fabs(StdI->Vp);
  }
  else{
    LargeValue0 = fabs(StdI->mu) * 2.0 + 2.0 * 2.0 * fabs(StdI->t) + 2.0 * 2.0 * fabs(StdI->tp)
        +fabs(StdI->U) + 2.0 * 2.0 * fabs(StdI->V) + 2.0 * 2.0 * fabs(StdI->Vp);
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
}

/**
  *
  * Setup a Hamiltonian for the generalized Heisenberg model on a Chain lattice
  *
  * @author Mitsuaki Kawamura (The University of Tokyo)
  */
void Spin_ChainLattice(
struct StdIntList *StdI,
  int Sz2 /**< [in] 2 * Total Sz */, 
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int iL, iL2;
  int ktrans, kintr;
  double LargeValue0, S;

  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_NotUsed_i("W", StdI->W);
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  /**/
  StdFace_PrintVal_i("2S", &StdI->S2, 1);
  StdFace_PrintVal_d("h", &StdI->h, 0.0);
  StdFace_PrintVal_d("Gamma", &StdI->Gamma, 0.0);
  StdFace_PrintVal_d("D", &StdI->D, 0.0);
  StdFace_PrintVal_d("J", &StdI->J, 1.0);
  StdFace_PrintVal_d("Jz", &StdI->Jz, StdI->J);
  StdFace_PrintVal_d("Jxy", &StdI->Jxy, StdI->J);
  StdFace_PrintVal_d("Jx", &StdI->Jx, StdI->Jxy);
  StdFace_PrintVal_d("Jy", &StdI->Jy, StdI->Jxy);
  StdFace_PrintVal_d("J'", &StdI->Jp, 0.0);
  StdFace_PrintVal_d("Jz'", &StdI->Jzp, StdI->Jp);
  StdFace_PrintVal_d("Jxy'", &StdI->Jxyp, StdI->Jp);
  StdFace_PrintVal_d("Jx'", &StdI->Jxp, StdI->Jxyp);
  StdFace_PrintVal_d("Jy'", &StdI->Jyp, StdI->Jxyp);
  StdI->Jxy = 0.5 * (StdI->Jx + StdI->Jy);
  StdI->Jxyp = 0.5 * (StdI->Jxp + StdI->Jyp);
  /**/
  StdFace_NotUsed_d("J0", StdI->J0);
  StdFace_NotUsed_d("J1", StdI->J1);
  StdFace_NotUsed_d("J2", StdI->J2);
  StdFace_NotUsed_d("Jz0", StdI->Jz0);
  StdFace_NotUsed_d("Jz1", StdI->Jz1);
  StdFace_NotUsed_d("Jxy0", StdI->Jxy0);
  StdFace_NotUsed_d("Jxy1", StdI->Jxy1);
  StdFace_NotUsed_d("K", StdI->K);
  /*
  Local Spin
  */
  StdI->nsite = StdI->L;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = StdI->S2;
  /*
  Transfer
  */
  StdI->ntrans = StdI->L * (StdI->S2 + 1 + 2 * StdI->S2);
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double *)malloc(sizeof(double) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    StdI->transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  StdI->ntrans = 0;
  for (isite = 0; isite < StdI->nsite; isite++){
    StdFace_MagLong(StdI, StdI->S2, -StdI->h, isite);
    StdFace_MagTrans(StdI, StdI->S2, -StdI->Gamma, isite);
  }
  /*
  Interaction
  */
  StdI->nintr = StdI->L * ((StdI->S2 + 1) * (StdI->S2 + 1) * 3 
    + 2 * StdI->S2 * StdI->S2 * 4);
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++){
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  StdI->nintr = 0;
  for (iL = 0; iL < StdI->L; iL++){

    isite = iL;
    StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->D, isite, isite);

    iL2 = (iL + 1) % StdI->L;
    jsite = iL2;
    StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->Jz, isite, jsite);
    StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->Jxy, isite, jsite);
    StdFace_Kitaev(StdI, StdI->S2, StdI->S2, 0.5 * (StdI->Jx - StdI->Jy), isite, jsite);

    iL2 = (iL + 2) % StdI->L;
    jsite = iL2;
    StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->Jzp, isite, jsite);
    StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->Jxyp, isite, jsite);
    StdFace_Kitaev(StdI, StdI->S2, StdI->S2, 0.5 * (StdI->Jxp - StdI->Jyp), isite, jsite);
  }
  /*
  Set mTPQ parameter
  */
  S = (double)StdI->S2 * 0.5;
  if (lGC == 0){
    LargeValue0 = (double)Sz2 / (double)(2 * StdI->nsite) * fabs(StdI->h) + S * fabs(StdI->D) + S * S * fabs(StdI->Gamma)
      + 2.0 / 2.0 * S * S * (fabs(StdI->Jx) + fabs(StdI->Jy) + fabs(StdI->Jz))
      + 2.0 / 2.0 * S * S * (fabs(StdI->Jxp) + fabs(StdI->Jyp) + fabs(StdI->Jzp));
  }
  else{
    LargeValue0 = S * fabs(StdI->h) + S * fabs(StdI->D) + S * S * fabs(StdI->Gamma)
      + 2.0 / 2.0 * S * S * (fabs(StdI->Jx) + fabs(StdI->Jy) + fabs(StdI->Jz))
      + 2.0 / 2.0 * S * S * (fabs(StdI->Jxp) + fabs(StdI->Jyp) + fabs(StdI->Jzp));
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
}

/**
*
* Setup a Hamiltonian for the generalized Heisenberg model on a Chain lattice
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void Spin_ChainLattice_Boost(
struct StdIntList *StdI,
  int Sz2 /**< [in] 2 * Total Sz */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int iL, iL2;
  int ktrans, kintr;
  int j, iW;
  double LargeValue0, S;
  FILE *fp;

  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_NotUsed_i("W", StdI->W);
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  /**/
  StdFace_PrintVal_i("2S", &StdI->S2, 1);
  StdFace_PrintVal_d("h", &StdI->h, 0.0);
  StdFace_PrintVal_d("Gamma", &StdI->Gamma, 0.0);
  StdFace_PrintVal_d("D", &StdI->D, 0.0);
  StdFace_PrintVal_d("J", &StdI->J, 1.0);
  StdFace_PrintVal_d("Jz", &StdI->Jz, StdI->J);
  StdFace_PrintVal_d("Jxy", &StdI->Jxy, StdI->J);
  StdFace_PrintVal_d("Jx", &StdI->Jx, StdI->Jxy);
  StdFace_PrintVal_d("Jy", &StdI->Jy, StdI->Jxy);
  StdFace_PrintVal_d("J'", &StdI->Jp, 0.0);
  StdFace_PrintVal_d("Jz'", &StdI->Jzp, StdI->Jp);
  StdFace_PrintVal_d("Jxy'", &StdI->Jxyp, StdI->Jp);
  StdFace_PrintVal_d("Jx'", &StdI->Jxp, StdI->Jxyp);
  StdFace_PrintVal_d("Jy'", &StdI->Jyp, StdI->Jxyp);
  /**/
  StdFace_NotUsed_d("J0", StdI->J0);
  StdFace_NotUsed_d("J1", StdI->J1);
  StdFace_NotUsed_d("J2", StdI->J2);
  StdFace_NotUsed_d("Jz0", StdI->Jz0);
  StdFace_NotUsed_d("Jz1", StdI->Jz1);
  StdFace_NotUsed_d("Jxy0", StdI->Jxy0);
  StdFace_NotUsed_d("Jxy1", StdI->Jxy1);
  StdFace_NotUsed_d("K", StdI->K);
  /*
  Local Spin
  */
  StdI->nsite = StdI->L;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = StdI->S2;
  /*
  Transfer
  */
  StdI->ntrans = 1;
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double *)malloc(sizeof(double) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++) {
    StdI->transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }
  StdI->ntrans = 0;
  /*
  Interaction
  */
  StdI->nintr = 1;
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++) {
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  StdI->nintr = 0;
  /*
  Set mTPQ parameter
  */
  S = (double)StdI->S2 * 0.5;
  if (lGC == 0) {
    LargeValue0 = (double)Sz2 / (double)(2 * StdI->nsite) * fabs(StdI->h) + S * fabs(StdI->D) + S * S * fabs(StdI->Gamma)
      + 2.0 / 2.0 * S * S * (fabs(StdI->Jx) + fabs(StdI->Jy) + fabs(StdI->Jz))
      + 2.0 / 2.0 * S * S * (fabs(StdI->Jxp) + fabs(StdI->Jyp) + fabs(StdI->Jzp));
  }
  else {
    LargeValue0 = S * fabs(StdI->h) + S * fabs(StdI->D) + S * S * fabs(StdI->Gamma)
      + 2.0 / 2.0 * S * S * (fabs(StdI->Jx) + fabs(StdI->Jy) + fabs(StdI->Jz))
      + 2.0 / 2.0 * S * S * (fabs(StdI->Jxp) + fabs(StdI->Jyp) + fabs(StdI->Jzp));
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
  /*
  Magnetic field
  */
  fp = fopen("boost.def", "w");
  fprintf(fp, "# Magnetic field\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    -0.5 * StdI->Gamma, 0.0, 0.0, 0.0, -0.5 *StdI->h, 0.0);
  /*
  Interaction
  */
  fprintf(fp, "%d  # Number of type of J\n", 2);
  fprintf(fp, "# J 1\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->Jx, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * StdI->Jy, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * StdI->Jz, 0.0);
  fprintf(fp, "# J 2\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->Jxp, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * StdI->Jyp, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * StdI->Jzp, 0.0);
  /*
  Topology
  */
  if (StdI->S2 != 1) {
    fprintf(stderr, "\n ERROR! S2 must be 1 in Boost. \n\n");
    exit(-1);
  }
  StdI->ishift_nspin = 4;
  StdI->W = 4;
  if(StdI->L % 4 != 0){
    fprintf(stderr, "\n ERROR! L % 4 != 0 \n\n");
    exit(-1);
  }
  StdI->L = StdI->L / 4;
  if (StdI->L < 2) {
    fprintf(stderr, "\n ERROR! L < 8 \n\n");
    exit(-1);
  }
  StdI->num_pivot = 1;
/**/
  fprintf(fp, "# W0  R0  StdI->num_pivot  StdI->ishift_nspin\n");
  fprintf(fp, "%d %d %d %d\n", StdI->W, StdI->L, StdI->num_pivot, StdI->ishift_nspin);

  StdI->list_6spin_star = (int **)malloc(sizeof(int*) * StdI->num_pivot);
  for (j = 0; j < StdI->num_pivot; j++) {
    StdI->list_6spin_star[j] = (int *)malloc(sizeof(int) * 7);
  }

  StdI->list_6spin_star[0][0] = 8; // num of J
  StdI->list_6spin_star[0][1] = 1;
  StdI->list_6spin_star[0][2] = 1;
  StdI->list_6spin_star[0][3] = 1;
  StdI->list_6spin_star[0][4] = 1;
  StdI->list_6spin_star[0][5] = 1;
  StdI->list_6spin_star[0][6] = 1; // flag

  fprintf(fp, "# StdI->list_6spin_star\n");
  for (j = 0; j < StdI->num_pivot; j++) {
    fprintf(fp, "# pivot %d\n", j);
    for (iW = 0; iW < 7; iW++) {
      fprintf(fp, "%d ", StdI->list_6spin_star[j][iW]);
    }
    fprintf(fp, "\n");
  }

  StdI->list_6spin_pair = (int ***)malloc(sizeof(int**) * StdI->num_pivot);
  for (j = 0; j < StdI->num_pivot; j++) {
    StdI->list_6spin_pair[j] = (int **)malloc(sizeof(int*) * 7);
    for (iW = 0; iW < 7; iW++) {
      StdI->list_6spin_pair[j][iW] = (int *)malloc(sizeof(int) * StdI->list_6spin_star[j][0]);
    }
  }

  StdI->list_6spin_pair[0][0][0] = 0; 
  StdI->list_6spin_pair[0][1][0] = 1;
  StdI->list_6spin_pair[0][2][0] = 2;
  StdI->list_6spin_pair[0][3][0] = 3;
  StdI->list_6spin_pair[0][4][0] = 4;
  StdI->list_6spin_pair[0][5][0] = 5; 
  StdI->list_6spin_pair[0][6][0] = 1; // type of J
  StdI->list_6spin_pair[0][0][1] = 1;
  StdI->list_6spin_pair[0][1][1] = 2;
  StdI->list_6spin_pair[0][2][1] = 0;
  StdI->list_6spin_pair[0][3][1] = 3;
  StdI->list_6spin_pair[0][4][1] = 4;
  StdI->list_6spin_pair[0][5][1] = 5;
  StdI->list_6spin_pair[0][6][1] = 1; // type of J
  StdI->list_6spin_pair[0][0][2] = 2;
  StdI->list_6spin_pair[0][1][2] = 3;
  StdI->list_6spin_pair[0][2][2] = 0;
  StdI->list_6spin_pair[0][3][2] = 1;
  StdI->list_6spin_pair[0][4][2] = 4;
  StdI->list_6spin_pair[0][5][2] = 5;
  StdI->list_6spin_pair[0][6][2] = 1; // type of J
  StdI->list_6spin_pair[0][0][3] = 3;
  StdI->list_6spin_pair[0][1][3] = 4;
  StdI->list_6spin_pair[0][2][3] = 0;
  StdI->list_6spin_pair[0][3][3] = 1;
  StdI->list_6spin_pair[0][4][3] = 2;
  StdI->list_6spin_pair[0][5][3] = 5;
  StdI->list_6spin_pair[0][6][3] = 1; // type of J
  StdI->list_6spin_pair[0][0][4] = 0;
  StdI->list_6spin_pair[0][1][4] = 2;
  StdI->list_6spin_pair[0][2][4] = 1;
  StdI->list_6spin_pair[0][3][4] = 3;
  StdI->list_6spin_pair[0][4][4] = 4;
  StdI->list_6spin_pair[0][5][4] = 5;
  StdI->list_6spin_pair[0][6][4] = 2; // type of J
  StdI->list_6spin_pair[0][0][5] = 1;
  StdI->list_6spin_pair[0][1][5] = 3;
  StdI->list_6spin_pair[0][2][5] = 0;
  StdI->list_6spin_pair[0][3][5] = 2;
  StdI->list_6spin_pair[0][4][5] = 4;
  StdI->list_6spin_pair[0][5][5] = 5;
  StdI->list_6spin_pair[0][6][5] = 2; // type of J
  StdI->list_6spin_pair[0][0][6] = 2;
  StdI->list_6spin_pair[0][1][6] = 4;
  StdI->list_6spin_pair[0][2][6] = 0;
  StdI->list_6spin_pair[0][3][6] = 1;
  StdI->list_6spin_pair[0][4][6] = 3;
  StdI->list_6spin_pair[0][5][6] = 5;
  StdI->list_6spin_pair[0][6][6] = 2; // type of J
  StdI->list_6spin_pair[0][0][7] = 3;
  StdI->list_6spin_pair[0][1][7] = 5;
  StdI->list_6spin_pair[0][2][7] = 0;
  StdI->list_6spin_pair[0][3][7] = 1;
  StdI->list_6spin_pair[0][4][7] = 2;
  StdI->list_6spin_pair[0][5][7] = 4;
  StdI->list_6spin_pair[0][6][7] = 2; // type of J

  fprintf(fp, "# StdI->list_6spin_pair\n");
  for (j = 0; j < StdI->num_pivot; j++) {
    fprintf(fp, "# pivot %d\n", j);
    for (iL = 0; iL < StdI->list_6spin_star[j][0]; iL++) {
      for (iW = 0; iW < 7; iW++) {
        fprintf(fp, "%d ", StdI->list_6spin_pair[j][iW][iL]);
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);

  for (j = 0; j < StdI->num_pivot; j++) {
    free(StdI->list_6spin_star[j]);
  }
  free(StdI->list_6spin_star);

  for (j = 0; j < StdI->num_pivot; j++) {
    for (iW = 0; iW < 7; iW++) {
      free(StdI->list_6spin_pair[j][iW]);
    }
    free(StdI->list_6spin_pair[j]);
  }
  free(StdI->list_6spin_pair);

}

/**
 *
 * Setup a Hamiltonian for the Kondo lattice model on a Chain lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void KondoLattice_ChainLattice(
struct StdIntList *StdI,
  int nelec /**< [in] The number of valence electrons */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int ispin;
  int iL;
  int ktrans, kintr;
  double LargeValue0, S;

  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_NotUsed_i("W", StdI->W);
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  /**/
  StdFace_PrintVal_i("2S", &StdI->S2, 1);
  StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
  StdFace_PrintVal_d("t", &StdI->t, 1.0);
  StdFace_PrintVal_d("J", &StdI->J, 0.0);
  /**/
  StdFace_NotUsed_d("U", StdI->U);
  StdFace_NotUsed_d("t'", StdI->tp);
  StdFace_NotUsed_d("tpp", StdI->tpp);
  StdFace_NotUsed_d("t0", StdI->t0);
  StdFace_NotUsed_d("t1", StdI->t1);
  StdFace_NotUsed_d("t2", StdI->t2);
  StdFace_NotUsed_d("V'", StdI->Vp);
  StdFace_NotUsed_d("Vpp", StdI->Vpp);
  StdFace_NotUsed_d("V0", StdI->V0);
  StdFace_NotUsed_d("V1", StdI->V1);
  StdFace_NotUsed_d("V2", StdI->V2);
  /*
  Local Spin
  */
  StdI->nsite = StdI->L * 2;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  for (iL = 0; iL < StdI->L; iL++){
    StdI->locspinflag[iL] = StdI->S2;
    StdI->locspinflag[iL + StdI->L] = 0;
  }
  /*
  Transfer
  */
  StdI->ntrans = StdI->L * 2 * 3;
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double *)malloc(sizeof(double) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    StdI->transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  StdI->ntrans = 0;
  for (iL = 0; iL < StdI->L; iL++){
    isite = StdI->L + iL;
    for (ispin = 0; ispin < 2; ispin++){

      StdFace_trans(StdI, StdI->mu, isite, ispin, isite, ispin);

      jsite = StdI->L + (iL + 1) % StdI->L;
      StdFace_trans(StdI, StdI->t, isite, ispin, jsite, ispin);
      StdFace_trans(StdI, StdI->t, jsite, ispin, isite, ispin);
    }
  }
  /*
  Interaction
  */
  StdI->nintr = StdI->L * ((StdI->S2 + 1) * (1 + 1) + 2 * StdI->S2 * 1);
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++){
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  StdI->nintr = 0;
  for (iL = 0; iL < StdI->L; iL++){

    isite = iL + StdI->L;
    jsite = iL;

    StdFace_exchange(StdI, 1, StdI->S2, StdI->J, isite, jsite);
    StdFace_SzSz(StdI, 1, StdI->S2, StdI->J, isite, jsite);
  }
  /*
  Set mTPQ parameter
  */
  S = (double)StdI->S2 * 0.5;
  if (lGC == 0){
    LargeValue0 = fabs(StdI->mu) * (double)nelec / (double)(StdI->L * StdI->W)
      + 2.0 * 2.0 * fabs(StdI->t) + 0.5 * S * fabs(StdI->J);
  }
  else{
    LargeValue0 = fabs(StdI->mu) * 2.0 + 2.0 * 2.0 * fabs(StdI->t) + 0.5 * S * fabs(StdI->J);
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
}
