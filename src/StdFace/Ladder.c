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
#include "StdFace_ModelUtil.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**
 *
 * Setup a Hamiltonian for the Hubbard model on a square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void FermionHubbard_Ladder(
struct StdIntList *StdI,
  int nelec /**< [in] The number of electrons */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int ispin;
  int iL, iW;
  int ktrans, kintr;
  double LargeValue0;
  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_RequiredVal_i("W", StdI->W);
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  /**/
  StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
  StdFace_PrintVal_d("U", &StdI->U, 0.0);
  StdFace_PrintVal_d("t1", &StdI->t1, 1.0);
  StdFace_PrintVal_d("t1'", &StdI->t1p, 1.0);
  StdFace_PrintVal_d("t0", &StdI->t0, 1.0);
  StdFace_PrintVal_d("t2", &StdI->t2, 1.0);
  StdFace_PrintVal_d("t2", &StdI->t2p, 1.0);
  StdFace_PrintVal_d("V0", &StdI->V0, 0.0);
  StdFace_PrintVal_d("V1", &StdI->V1, 0.0);
  StdFace_PrintVal_d("V1'", &StdI->V1p, 0.0);
  StdFace_PrintVal_d("V2", &StdI->V2, 0.0);
  /**/
  StdFace_NotUsed_i("2S", StdI->S2);
  StdFace_NotUsed_d("t", StdI->t);
  StdFace_NotUsed_d("t'", StdI->tp);
  StdFace_NotUsed_d("t''", StdI->tpp);
  StdFace_NotUsed_d("V", StdI->V);
  StdFace_NotUsed_d("V'", StdI->Vp);
  StdFace_NotUsed_d("V''", StdI->Vpp);
  /*
  Local Spin
  */
  StdI->nsite = StdI->L * StdI->W;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = 0;
  /*
  Transfer
  */
  StdI->ntrans = 2 * (5 * StdI->L*StdI->W + 6* StdI->L*(StdI->W - 1));
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double *)malloc(sizeof(double) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    StdI->transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  StdI->ntrans = 0;
  for (iL = 0; iL < StdI->L; iL++) {
    for (iW = 0; iW < StdI->W; iW++){
      isite = iW + iL * StdI->W;
      for (ispin = 0; ispin < 2; ispin++){

        StdFace_trans(StdI, StdI->mu, isite, ispin, isite, ispin);

        jsite = iW + ((iL + 1 + 2 * StdI->L) % StdI->L) * StdI->W;
        StdFace_trans(StdI, StdI->t1, isite, ispin, jsite, ispin);
        StdFace_trans(StdI, StdI->t1, jsite, ispin, isite, ispin);

        jsite = iW + ((iL + 2 + 2 * StdI->L) % StdI->L) * StdI->W;
        StdFace_trans(StdI, StdI->t1p, isite, ispin, jsite, ispin);
        StdFace_trans(StdI, StdI->t1p, jsite, ispin, isite, ispin);

        if (iW < StdI->W - 1){
          jsite = (iW + 1) + iL * StdI->W;
          StdFace_trans(StdI, StdI->t0, isite, ispin, jsite, ispin);
          StdFace_trans(StdI, StdI->t0, jsite, ispin, isite, ispin);

          jsite = (iW + 1) + ((iL + 1 + 2 * StdI->L) % StdI->L) * StdI->W;
          StdFace_trans(StdI, StdI->t2, isite, ispin, jsite, ispin);
          StdFace_trans(StdI, StdI->t2, jsite, ispin, isite, ispin);

          jsite = (iW + 1) + ((iL - 1 + 2 * StdI->L) % StdI->L) * StdI->W;
          StdFace_trans(StdI, StdI->t2p, isite, ispin, jsite, ispin);
          StdFace_trans(StdI, StdI->t2p, jsite, ispin, isite, ispin);
        }

      }
    }
  }
  /*
  Interaction
  */
  StdI->nintr = StdI->L*StdI->W + 4 * 
    (StdI->L*(StdI->W - 1) + StdI->L*StdI->W + StdI->L*StdI->W + 2 * StdI->L*(StdI->W - 1));
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++){
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  StdI->nintr = 0;
  for (iW = 0; iW <StdI->W; iW++){
    for (iL = 0; iL < StdI->L; iL++){

      isite = iW + iL * StdI->W;

      StdFace_intr(StdI, StdI->U, isite, 0, isite, 0, isite, 1, isite, 1);

      jsite = iW + ((iL + 1 + 2 * StdI->L) % StdI->L) * StdI->W;
      StdFace_Coulomb(StdI, StdI->V1, isite, jsite);

      jsite = iW + ((iL + 2 + 2 * StdI->L) % StdI->L) * StdI->W;
      StdFace_Coulomb(StdI, StdI->V1p, isite, jsite);

      if (iW < StdI->W - 1){
        jsite = (iW + 1) + iL * StdI->W;
        StdFace_Coulomb(StdI, StdI->V0, isite, jsite);

        jsite = (iW + 1) + ((iL + 1 + 2 * StdI->L) % StdI->L) * StdI->W;
        StdFace_Coulomb(StdI, StdI->V2, isite, jsite);

        jsite = (iW + 1) + ((iL - 1 + 2 * StdI->L) % StdI->L) * StdI->W;
        StdFace_Coulomb(StdI, StdI->V2p, isite, jsite);
      }
    }
  }
  /*
  Set mTPQ parameter
  */
  if (lGC == 0){
    LargeValue0 =  fabs(StdI->mu) * (double)nelec / (double)(StdI->L * StdI->W)
      + 2.0 * 2.0 * fabs(StdI->t1) + 2.0 * 2.0 * fabs(StdI->t1p)
      + 2.0 * 2.0 * fabs(StdI->t0) + 2.0 * 2.0 * fabs(StdI->t2) + 2.0 * 2.0 * fabs(StdI->t2p)
      + fabs(StdI->U)
      + 4.0 * fabs(StdI->V1) + 4.0 * fabs(StdI->V1p) + 4.0 * fabs(StdI->V0)
      + 2.0 * 2.0 * fabs(StdI->V2) + 2.0 * 2.0 * fabs(StdI->V2p);
  }
  else{
    LargeValue0 = fabs(StdI->mu) * 2.0
      + 2.0 * 2.0 * fabs(StdI->t1) + 2.0 * 2.0 * fabs(StdI->t1p)
      + 2.0 * 2.0 * fabs(StdI->t0) + 2.0 * 2.0 * fabs(StdI->t2) + 2.0 * 2.0 * fabs(StdI->t2p)
      + fabs(StdI->U)
      + 4.0 * fabs(StdI->V1) + 4.0 * fabs(StdI->V1p) + 4.0 * fabs(StdI->V0)
      + 2.0 * 2.0 * fabs(StdI->V2) + 2.0 * 2.0 * fabs(StdI->V2p);
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
}

/**
 *
 * Setup a Hamiltonian for the generalized Heisenberg model on a square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void Spin_Ladder(
struct StdIntList *StdI,
  int Sz2 /**< [in] 2 * Total Sz */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int iL, iW;
  int ktrans, kintr;
  double LargeValue0, S;

  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_RequiredVal_i("W", StdI->W);
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  /**/
  StdFace_PrintVal_i("2S", &StdI->S2, 1);
  StdFace_PrintVal_d("h", &StdI->h, 0.0);
  StdFace_PrintVal_d("Gamma", &StdI->Gamma, 0.0);
  StdFace_PrintVal_d("D", &StdI->D, 0.0);
  StdFace_PrintVal_d("J1", &StdI->J1, 1.0);
  StdFace_PrintVal_d("J1'", &StdI->J1p, 1.0);
  StdFace_PrintVal_d("J0", &StdI->J0, 1.0);
  StdFace_PrintVal_d("J2", &StdI->J2, 1.0);
  StdFace_PrintVal_d("J2'", &StdI->J2p, 1.0);
  /**/
  StdFace_NotUsed_d("J", StdI->J);
  StdFace_NotUsed_d("J'", StdI->Jp);
  StdFace_NotUsed_d("J''", StdI->Jpp);
  StdFace_NotUsed_d("Jxy", StdI->Jxy);
  StdFace_NotUsed_d("Jx", StdI->Jx);
  StdFace_NotUsed_d("Jy", StdI->Jy);
  StdFace_NotUsed_d("Jz", StdI->Jz);
  StdFace_NotUsed_d("Jz0", StdI->Jz0);
  StdFace_NotUsed_d("Jz1", StdI->Jz1);
  StdFace_NotUsed_d("Jxy0", StdI->Jxy0);
  StdFace_NotUsed_d("Jxy1", StdI->Jxy1);
  StdFace_NotUsed_d("K", StdI->K);
  /*
  Local Spin
  */
  StdI->nsite = StdI->L * StdI->W;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = StdI->S2;
  /*
  Transfer
  */
  StdI->ntrans = StdI->L * StdI->W * (StdI->S2 + 1 + 2 * StdI->S2);
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
  StdI->nintr = (3* StdI->L*StdI->W + 3* StdI->L*(StdI->W - 1)) 
    * (StdI->S2 + 1) * (StdI->S2 + 1) 
    + (2* StdI->L*StdI->W + 3* StdI->L*(StdI->W - 1))*StdI->S2*StdI->S2*2;
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++){
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  StdI->nintr = 0;

  for (iL = 0; iL < StdI->L; iL++) {
    for (iW = 0; iW < StdI->W; iW++){

      isite = iW + iL * StdI->W;

      StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->D, isite, isite);
      
      jsite = iW + ((iL + 1 + 2 * StdI->L) % StdI->L) * StdI->W;
      StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->J1, isite, jsite);
      StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->J1, isite, jsite);
      
      jsite = iW + ((iL + 2 + 2 * StdI->L) % StdI->L) * StdI->W;
      StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->J1p, isite, jsite);
      StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->J1p, isite, jsite);
      
      
      if (iW < StdI->W - 1){
        jsite = (iW + 1) + iL * StdI->W;
        StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->J0, isite, jsite);
        StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->J0, isite, jsite);

        jsite = (iW + 1) + ((iL + 1 + 2 * StdI->L) % StdI->L) * StdI->W;
        StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->J2, isite, jsite);
        StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->J2, isite, jsite);

        jsite = (iW + 1) + ((iL - 1 + 2 * StdI->L) % StdI->L) * StdI->W;
        StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->J2p, isite, jsite);
        StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->J2p, isite, jsite);
      }
    }
  }
  /*
  Set mTPQ parameter
  */
  S = (double)StdI->S2 * 0.5;
  if (lGC == 0){
    LargeValue0 = (double)Sz2 / (double)(2 * StdI->nsite) * fabs(StdI->h) 
      + S * fabs(StdI->D) + 0.5 * S * fabs(StdI->Gamma)
    + S*S*fabs(StdI->J1) + S*S*fabs(StdI->J1p) + S*S*fabs(StdI->J0) 
      + S*S*fabs(StdI->J2) + S*S*fabs(StdI->J2p);
  }
  else{
    LargeValue0 = S * fabs(StdI->h) + S * fabs(StdI->D) + 0.5 * S * fabs(StdI->Gamma)
      + S*S*fabs(StdI->J1) + S*S*fabs(StdI->J1p)
      + S*S*fabs(StdI->J0) + S*S*fabs(StdI->J2) + S*S*fabs(StdI->J2p);
   }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
}

/**
*
* Setup a Hamiltonian for the generalized Heisenberg model on a square lattice
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void Spin_Ladder_Boost(
struct StdIntList *StdI,
  int Sz2 /**< [in] 2 * Total Sz */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int iL, iW;
  int ktrans, kintr;
  double LargeValue0, S;
  FILE *fp;
  int j;

  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_RequiredVal_i("W", StdI->W);
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  /**/
  StdFace_PrintVal_i("2S", &StdI->S2, 1);
  StdFace_PrintVal_d("h", &StdI->h, 0.0);
  StdFace_PrintVal_d("Gamma", &StdI->Gamma, 0.0);
  StdFace_PrintVal_d("D", &StdI->D, 0.0);
  StdFace_PrintVal_d("J1", &StdI->J1, 1.0);
  StdFace_PrintVal_d("J1'", &StdI->J1p, 1.0);
  StdFace_PrintVal_d("J0", &StdI->J0, 1.0);
  StdFace_PrintVal_d("J2", &StdI->J2, 1.0);
  StdFace_PrintVal_d("J2'", &StdI->J2p, 1.0);
  /**/
  StdFace_NotUsed_d("J", StdI->J);
  StdFace_NotUsed_d("J'", StdI->Jp);
  StdFace_NotUsed_d("J''", StdI->Jpp);
  StdFace_NotUsed_d("Jxy", StdI->Jxy);
  StdFace_NotUsed_d("Jx", StdI->Jx);
  StdFace_NotUsed_d("Jy", StdI->Jy);
  StdFace_NotUsed_d("Jz", StdI->Jz);
  StdFace_NotUsed_d("Jz0", StdI->Jz0);
  StdFace_NotUsed_d("Jz1", StdI->Jz1);
  StdFace_NotUsed_d("Jxy0", StdI->Jxy0);
  StdFace_NotUsed_d("Jxy1", StdI->Jxy1);
  StdFace_NotUsed_d("K", StdI->K);
  /*
  Local Spin
  */
  StdI->nsite = StdI->L * StdI->W;
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
    LargeValue0 = (double)Sz2 / (double)(2 * StdI->nsite) * fabs(StdI->h)
      + S * fabs(StdI->D) + 0.5 * S * fabs(StdI->Gamma)
      + S*S*fabs(StdI->J1) + S*S*fabs(StdI->J1p) + S*S*fabs(StdI->J0) 
      + S*S*fabs(StdI->J2) + S*S*fabs(StdI->J2p);
  }
  else {
    LargeValue0 = S * fabs(StdI->h) + S * fabs(StdI->D) + 0.5 * S * fabs(StdI->Gamma)
      + S*S*fabs(StdI->J1) + S*S*fabs(StdI->J1p)
      + S*S*fabs(StdI->J0) + S*S*fabs(StdI->J2) + S*S*fabs(StdI->J2p);
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
  /*
  Magnetic field
  */
  fp = fopen("boost.def", "w");
  fprintf(fp, "# Magnetic field\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    -0.5 * StdI->Gamma, 0.0, 0.0, 0.0, -0.5 * StdI->h, 0.0);
  /*
  Interaction
  */
  fprintf(fp, "%d  # Number of type of J\n", 5);
  fprintf(fp, "# J 1 (inter chain, vertical)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J0, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * StdI->J0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * StdI->J0, 0.0);
  fprintf(fp, "# J 2 (Nearest neighbor, along chain)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * StdI->J1, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * StdI->J1, 0.0);
  fprintf(fp, "# J 3 (Second nearest neighbor, along chain)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1p, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * StdI->J1p, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * StdI->J1p, 0.0);
  fprintf(fp, "# J 4 (inter chain, diagonal1)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * StdI->J2, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * StdI->J2, 0.0);
  fprintf(fp, "# J 5 (inter chain, diagonal2)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2p, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 *StdI->J2p, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * StdI->J2p, 0.0);
  /*
  Topology
  */
  if (StdI->S2 != 1) {
    fprintf(stderr, "\n ERROR! S2 must be 1 in Boost. \n\n");
    exit(-1);
  }
  StdI->ishift_nspin = 2;
  if (StdI->W != 2) {
    fprintf(stderr, "\n ERROR! W != 2 \n\n");
    exit(-1);
  }
  if (StdI->L < 3) {
    fprintf(stderr, "\n ERROR! L < 3 \n\n");
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

  StdI->list_6spin_star[0][0] = 7; // num of J
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
  StdI->list_6spin_pair[0][0][1] = 0;
  StdI->list_6spin_pair[0][1][1] = 2;
  StdI->list_6spin_pair[0][2][1] = 1;
  StdI->list_6spin_pair[0][3][1] = 3;
  StdI->list_6spin_pair[0][4][1] = 4;
  StdI->list_6spin_pair[0][5][1] = 5;
  StdI->list_6spin_pair[0][6][1] = 2; // type of J
  StdI->list_6spin_pair[0][0][2] = 1;
  StdI->list_6spin_pair[0][1][2] = 3;
  StdI->list_6spin_pair[0][2][2] = 0;
  StdI->list_6spin_pair[0][3][2] = 2;
  StdI->list_6spin_pair[0][4][2] = 4;
  StdI->list_6spin_pair[0][5][2] = 5;
  StdI->list_6spin_pair[0][6][2] = 2; // type of J
  StdI->list_6spin_pair[0][0][3] = 0;
  StdI->list_6spin_pair[0][1][3] = 4;
  StdI->list_6spin_pair[0][2][3] = 1;
  StdI->list_6spin_pair[0][3][3] = 2;
  StdI->list_6spin_pair[0][4][3] = 3;
  StdI->list_6spin_pair[0][5][3] = 5;
  StdI->list_6spin_pair[0][6][3] = 3; // type of J
  StdI->list_6spin_pair[0][0][4] = 1;
  StdI->list_6spin_pair[0][1][4] = 5;
  StdI->list_6spin_pair[0][2][4] = 0;
  StdI->list_6spin_pair[0][3][4] = 2;
  StdI->list_6spin_pair[0][4][4] = 3;
  StdI->list_6spin_pair[0][5][4] = 4;
  StdI->list_6spin_pair[0][6][4] = 3; // type of J
  StdI->list_6spin_pair[0][0][5] = 0;
  StdI->list_6spin_pair[0][1][5] = 3;
  StdI->list_6spin_pair[0][2][5] = 1;
  StdI->list_6spin_pair[0][3][5] = 2;
  StdI->list_6spin_pair[0][4][5] = 4;
  StdI->list_6spin_pair[0][5][5] = 5;
  StdI->list_6spin_pair[0][6][5] = 4; // type of J
  StdI->list_6spin_pair[0][0][6] = 1;
  StdI->list_6spin_pair[0][1][6] = 2;
  StdI->list_6spin_pair[0][2][6] = 0;
  StdI->list_6spin_pair[0][3][6] = 3;
  StdI->list_6spin_pair[0][4][6] = 4;
  StdI->list_6spin_pair[0][5][6] = 5;
  StdI->list_6spin_pair[0][6][6] = 5; // type of J

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
 * Setup a Hamiltonian for the Kondo lattice model on a square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void KondoLattice_Ladder(
struct StdIntList *StdI,
  int nelec /**< [in] The number of valence electrons */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int ispin;
  int iL, iW;
  int ktrans, kintr;
  double LargeValue0, S;
  /**/
  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_RequiredVal_i("W", StdI->W);
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  /**/
  StdFace_PrintVal_i("2S", &StdI->S2, 1);
  StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
  StdFace_PrintVal_d("t0", &StdI->t0, 1.0);
  StdFace_PrintVal_d("t1", &StdI->t1, 1.0);
  StdFace_PrintVal_d("t1'", &StdI->t1p, 1.0);
  StdFace_PrintVal_d("t2", &StdI->t2, 1.0);
  StdFace_PrintVal_d("t2'", &StdI->t2p, 1.0);
  StdFace_PrintVal_d("J", &StdI->J, 0.0);
  /**/
  StdFace_NotUsed_d("U", StdI->U);
  StdFace_NotUsed_d("t", StdI->t);
  StdFace_NotUsed_d("t'", StdI->tp);
  StdFace_NotUsed_d("tpp", StdI->tpp);
  StdFace_NotUsed_d("V'", StdI->Vp);
  StdFace_NotUsed_d("Vpp", StdI->Vpp);
  StdFace_NotUsed_d("V0", StdI->V0);
  StdFace_NotUsed_d("V1", StdI->V1);
  StdFace_NotUsed_d("V2", StdI->V2);
  /*
  Local Spin
  */
  StdI->nsite = 2 * StdI->L * StdI->W;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  for (iL = 0; iL < StdI->L * StdI->W; iL++){
    StdI->locspinflag[iL] = StdI->S2;
    StdI->locspinflag[iL + StdI->L * StdI->W] = 0;
  }
  /*
  Transfer
  */
  StdI->ntrans = 2 * (5 * StdI->L*StdI->W + 6 * StdI->L*(StdI->W - 1));
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double *)malloc(sizeof(double) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    StdI->transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  StdI->ntrans = 0;
  for (iL = 0; iL <StdI->L; iL++) {
    for (iW = 0; iW < StdI->W; iW++){
      isite = 2 * StdI->L * StdI->W + iW + iL * StdI->W;
      for (ispin = 0; ispin < 2; ispin++){

        StdFace_trans(StdI, StdI->mu, isite, ispin, isite, ispin);

        jsite = StdI->L * StdI->W + iW + ((iL + 1 + 2 * StdI->L) % StdI->L) * StdI->W;
        StdFace_trans(StdI, StdI->t1, isite, ispin, jsite, ispin);
        StdFace_trans(StdI, StdI->t1, jsite, ispin, isite, ispin);

        jsite = StdI->L * StdI->W + iW + ((iL + 2 + 2 * StdI->L) % StdI->L) * StdI->W;
        StdFace_trans(StdI, StdI->t1p, isite, ispin, jsite, ispin);
        StdFace_trans(StdI, StdI->t1p, jsite, ispin, isite, ispin);

        if (iW < StdI->W - 1){
          jsite = StdI->L * StdI->W + (iW + 1) + iL * StdI->W;
          StdFace_trans(StdI, StdI->t0, isite, ispin, jsite, ispin);
          StdFace_trans(StdI, StdI->t0, jsite, ispin, isite, ispin);

          jsite = StdI->L * StdI->W + (iW + 1) 
            + ((iL + 1 + 2 * StdI->L) % StdI->L) * StdI->W;
          StdFace_trans(StdI, StdI->t2, isite, ispin, jsite, ispin);
          StdFace_trans(StdI, StdI->t2, jsite, ispin, isite, ispin);

          jsite = StdI->L * StdI->W
            + (iW + 1) + ((iL - 1 + 2 * StdI->L) % StdI->L) * StdI->W;
          StdFace_trans(StdI, StdI->t2p, isite, ispin, jsite, ispin);
          StdFace_trans(StdI, StdI->t2p, jsite, ispin, isite, ispin);
        }

      }
    }
  }
  /*
  Interaction
  */
  StdI->nintr = StdI->L * StdI->W * ((StdI->S2 + 1) * (1 + 1) + 2 * StdI->S2 * 1);
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++){
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  StdI->nintr = 0;
  for (iL = 0; iL < StdI->L; iL++) {
    for (iW = 0; iW < StdI->W; iW++) {

      isite = iL + iW * StdI->L + StdI->L * StdI->W;
      jsite = iL + iW * StdI->L;

      StdFace_exchange(StdI, 1, StdI->S2, StdI->J, isite, jsite);
      StdFace_SzSz(StdI, 1, StdI->S2, StdI->J, isite, jsite);
    }
  }
  /*
  Set mTPQ parameter
  */
  S = (double)StdI->S2 * 0.5;
  if (lGC == 0){
    LargeValue0 = fabs(StdI->mu) * (double)nelec / (double)(StdI->L * StdI->W)
      + 2.0 * 2.0 * fabs(StdI->t1) + 2.0 * 2.0 * fabs(StdI->t1p)
      + 2.0 * 2.0 * fabs(StdI->t0) + 2.0 * 2.0 * fabs(StdI->t2) + 2.0 * 2.0 * fabs(StdI->t2p)
      + 0.5 * S * fabs(StdI->J);
  }
  else{
    LargeValue0 = fabs(StdI->mu) * 2.0
      + 2.0 * 2.0 * fabs(StdI->t1) + 2.0 * 2.0 * fabs(StdI->t1p)
      + 2.0 * 2.0 * fabs(StdI->t0) + 2.0 * 2.0 * fabs(StdI->t2) + 2.0 * 2.0 * fabs(StdI->t2p)
      + 0.5 * S * fabs(StdI->J);
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
}
