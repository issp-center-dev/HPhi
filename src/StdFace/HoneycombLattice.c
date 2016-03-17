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
#include "../include/mfmemory.h"

/**
 *
 * Setup a Hamiltonian for the Hubbard model on a Honeycomb lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void FermionHubbard_HoneycombLattice(
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
  StdFace_PrintVal_d("t", &StdI->t, 1.0);
  StdFace_PrintVal_d("t0", &StdI->t0, StdI->t);
  StdFace_PrintVal_d("t1", &StdI->t1, StdI->t);
  StdFace_PrintVal_d("t2", &StdI->t2, StdI->t);
  StdFace_PrintVal_d("V", &StdI->V, 0.0);
  StdFace_PrintVal_d("V0", &StdI->V0, StdI->V);
  StdFace_PrintVal_d("V1", &StdI->V1, StdI->V);
  StdFace_PrintVal_d("V2", &StdI->V2, StdI->V);
  /**/
  StdFace_NotUsed_i("2S", StdI->S2);
  StdFace_NotUsed_d("t'", StdI->tp);
  StdFace_NotUsed_d("tpp", StdI->tpp);
  StdFace_NotUsed_d("V'", StdI->Vp);
  StdFace_NotUsed_d("Vpp", StdI->Vpp);
  /*
  Local Spin
  */
  StdI->nsite = StdI->L * StdI->W * 2;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = 0;
  /*
  Transfer
  */
  StdI->ntrans = StdI->L * StdI->W * 2 * 2 * 4;
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double *)malloc(sizeof(double) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    StdI->transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  StdI->ntrans = 0;
  for (iW = 0; iW < StdI->W; iW++){
    for (iL = 0; iL < StdI->L; iL++){
      for (ispin = 0; ispin < 2; ispin++){
        isite = (iL + iW * StdI->L) * 2;
        StdFace_trans(StdI, StdI->mu, isite, ispin, isite, ispin);
        jsite = (iL + iW * StdI->L) * 2 + 1;
        StdFace_trans(StdI, StdI->t0, isite, ispin, jsite, ispin);
        jsite = (((iL - 1 + 2 * StdI->L) % StdI->L) + ((iW + 0 + 2 * StdI->W) % StdI->W) * StdI->L) * 2 + 1;
        StdFace_trans(StdI, StdI->t1, isite, ispin, jsite, ispin);
        jsite = (((iL + 0 + 2 * StdI->L) % StdI->L) + ((iW - 1 + 2 * StdI->W) % StdI->W) * StdI->L) * 2 + 1;
        StdFace_trans(StdI, StdI->t2, isite, ispin, jsite, ispin);

        isite = (iL + iW * StdI->L) * 2 + 1;
        StdFace_trans(StdI, StdI->mu, isite, ispin, isite, ispin);
        jsite = (iL + iW * StdI->L) * 2;
        StdFace_trans(StdI, StdI->t0, isite, ispin, jsite, ispin);
        jsite = (((iL + 1 + 2 * StdI->L) % StdI->L) + ((iW + 0 + 2 * StdI->W) % StdI->W) * StdI->L) * 2;
        StdFace_trans(StdI, StdI->t1, isite, ispin, jsite, ispin);
        jsite = (((iL + 0 + 2 * StdI->L) % StdI->L) + ((iW + 1 + 2 * StdI->W) % StdI->W) * StdI->L) * 2;
        StdFace_trans(StdI, StdI->t2, isite, ispin, jsite, ispin);
      }
    }
  }
  /*
  Interaction
  */
  StdI->nintr = StdI->L * StdI->W * (2 + 4 * 3);
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++){
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  StdI->nintr = 0;
  for (iW = 0; iW < StdI->W; iW++){
    for (iL = 0; iL < StdI->L; iL++){
      isite = (iL + iW * StdI->L) * 2;
      StdFace_intr(StdI, StdI->U, isite, 0, isite, 0, isite, 1, isite, 1);
      isite = (iL + iW * StdI->L) * 2 + 1;
      StdFace_intr(StdI, StdI->U, isite, 0, isite, 0, isite, 1, isite, 1);

      jsite = (iL + iW * StdI->L) * 2;
      StdFace_Coulomb(StdI, StdI->V0, isite, jsite);
      jsite = (((iL + 1 + 2 * StdI->L) % StdI->L) + ((iW + 0 + 2 * StdI->W) % StdI->W) * StdI->L) * 2;
      StdFace_Coulomb(StdI, StdI->V1, isite, jsite);
      jsite = (((iL + 0 + 2 * StdI->L) % StdI->L) + ((iW + 1 + 2 * StdI->W) % StdI->W) * StdI->L) * 2;
      StdFace_Coulomb(StdI, StdI->V2, isite, jsite);
    }
  }
  /*
  Set mTPQ parameter
  */
  if (lGC == 0){
    LargeValue0 =  fabs(StdI->mu) * (double)nelec / (double)(2 * StdI->L * StdI->W)
      + 2.0 * 3.0 * fabs(StdI->t)
      + fabs(StdI->U) + 2.0 * 3.0 * fabs(StdI->V);
  }
  else{
    LargeValue0 = fabs(StdI->mu) * 2.0 + 2.0 * 3.0 * fabs(StdI->t)
      + fabs(StdI->U) + 2.0 * 3.0 * fabs(StdI->V);
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
}

/**
 *
 * Setup a Hamiltonian for the generalized Heisenberg model on a Heisenberg lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void Spin_HoneycombLattice(
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
  StdFace_PrintVal_d("J", &StdI->J, 0.0);
  StdFace_PrintVal_d("J0", &StdI->J0, StdI->J);
  StdFace_PrintVal_d("J1", &StdI->J1, StdI->J);
  StdFace_PrintVal_d("J2", &StdI->J2, StdI->J);
  StdFace_PrintVal_d("Jz", &StdI->Jz, StdI->J);
  StdFace_PrintVal_d("Jz0", &StdI->Jz0, StdI->Jz);
  StdFace_PrintVal_d("Jz1", &StdI->Jz1, StdI->Jz);
  StdFace_PrintVal_d("Jz2", &StdI->Jz2, StdI->Jz);
  StdFace_PrintVal_d("Jxy", &StdI->Jxy, StdI->J);
  StdFace_PrintVal_d("Jxy0", &StdI->Jxy0, StdI->Jxy);
  StdFace_PrintVal_d("Jxy1", &StdI->Jxy1, StdI->Jxy);
  StdFace_PrintVal_d("Jxy2", &StdI->Jxy2, StdI->Jxy);
  StdFace_PrintVal_d("Jx", &StdI->Jx, StdI->Jxy);
  StdFace_PrintVal_d("Jx0", &StdI->Jx0, StdI->Jx);
  StdFace_PrintVal_d("Jx1", &StdI->Jx1, StdI->Jx);
  StdFace_PrintVal_d("Jx2", &StdI->Jx2, StdI->Jx);
  StdFace_PrintVal_d("Jy", &StdI->Jy, StdI->Jxy);
  StdFace_PrintVal_d("Jy0", &StdI->Jy0, StdI->Jy);
  StdFace_PrintVal_d("Jy1", &StdI->Jy1, StdI->Jy);
  StdFace_PrintVal_d("Jy2", &StdI->Jy2, StdI->Jy);
  StdI->Jxy = 0.5 * (StdI->Jx + StdI->Jy);
  StdI->Jxy0 = 0.5 * (StdI->Jx0 + StdI->Jy0);
  StdI->Jxy1 = 0.5 * (StdI->Jx1 + StdI->Jy1);
  StdI->Jxy2 = 0.5 * (StdI->Jx2 + StdI->Jy2);
  /**/
  StdFace_NotUsed_d("J'", StdI->Jp);
  StdFace_NotUsed_d("Jz'", StdI->Jzp);
  StdFace_NotUsed_d("Jxy'", StdI->Jxyp);
  StdFace_NotUsed_d("Jx'", StdI->Jxp);
  StdFace_NotUsed_d("Jy'", StdI->Jyp);
  StdFace_NotUsed_d("K", StdI->K);
  /*
  Local Spin
  */
  StdI->nsite = StdI->L * StdI->W * 2;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = StdI->S2;
  /*
  Transfer
  */
  StdI->ntrans = StdI->L * StdI->W * 2 * (StdI->S2 + 1 + 2 * StdI->S2);
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
  StdI->nintr = StdI->L * StdI->W * ((StdI->S2 + 1) * (StdI->S2 + 1) * 5 + 2 * StdI->S2 * StdI->S2 * 6);
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++){
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  StdI->nintr = 0;
  for (iW = 0; iW < StdI->W; iW++){
    for (iL = 0; iL < StdI->L; iL++){

      isite = (iL + iW * StdI->L) * 2;
      StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->D, isite, isite);
      isite = (iL + iW * StdI->L) * 2 + 1;
      StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->D, isite, isite);

      jsite = (iL + iW * StdI->L) * 2;
      StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->Jz0, isite, jsite);
      StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->Jxy0, isite, jsite);
      StdFace_Kitaev(StdI, StdI->S2, StdI->S2, 0.5 * (StdI->Jx0 - StdI->Jy0), isite, jsite);
      jsite = (((iL + 1 + 2 * StdI->L) % StdI->L) + ((iW + 0 + 2 * StdI->W) % StdI->W) * StdI->L) * 2;
      StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->Jz1, isite, jsite);
      StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->Jxy1, isite, jsite);
      StdFace_Kitaev(StdI, StdI->S2, StdI->S2, 0.5 * (StdI->Jx1 - StdI->Jy1), isite, jsite);
      jsite = (((iL + 0 + 2 * StdI->L) % StdI->L) + ((iW + 1 + 2 * StdI->W) % StdI->W) * StdI->L) * 2;
      StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->Jz2, isite, jsite);
      StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->Jxy2, isite, jsite);
      StdFace_Kitaev(StdI, StdI->S2, StdI->S2, 0.5 * (StdI->Jx2 - StdI->Jy2), isite, jsite);
    }
  }
  /*
  Set mTPQ parameter
  */
  S = (double)StdI->S2 * 0.5;
  if (lGC == 0){
    LargeValue0 = (double)Sz2 / (double)(2 * StdI->nsite) * fabs(StdI->h) + S * fabs(StdI->D) + S * S * fabs(StdI->Gamma)
      + 1.0 / 2.0 * S * S * (fabs(StdI->Jx0) + fabs(StdI->Jy0) + fabs(StdI->Jz0))
      + 1.0 / 2.0 * S * S * (fabs(StdI->Jx1) + fabs(StdI->Jy1) + fabs(StdI->Jz1))
      + 1.0 / 2.0 * S * S * (fabs(StdI->Jx2) + fabs(StdI->Jy2) + fabs(StdI->Jz2));
  }
  else{
    LargeValue0 = S * fabs(StdI->h) + S * fabs(StdI->D) + S * S * fabs(StdI->Gamma)
      + 1.0 / 2.0 * S * S * (fabs(StdI->Jx0) + fabs(StdI->Jy0) + fabs(StdI->Jz0))
      + 1.0 / 2.0 * S * S * (fabs(StdI->Jx1) + fabs(StdI->Jy1) + fabs(StdI->Jz1))
      + 1.0 / 2.0 * S * S * (fabs(StdI->Jx2) + fabs(StdI->Jy2) + fabs(StdI->Jz2));
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
}

/**
*
* Setup a Hamiltonian for the generalized Heisenberg model on a Heisenberg lattice
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void Spin_HoneycombLattice_Boost(
struct StdIntList *StdI,
  int Sz2 /**< [in] 2 * Total Sz */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, ipivot;
  int ktrans, kintr;
  double LargeValue0, S;
  FILE *fp;

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
  StdFace_PrintVal_d("J", &StdI->J, 0.0);
  StdFace_PrintVal_d("J0", &StdI->J0, StdI->J);
  StdFace_PrintVal_d("J1", &StdI->J1, StdI->J);
  StdFace_PrintVal_d("J2", &StdI->J2, StdI->J);
  StdFace_PrintVal_d("Jz", &StdI->Jz, StdI->J);
  StdFace_PrintVal_d("Jz0", &StdI->Jz0, StdI->Jz);
  StdFace_PrintVal_d("Jz1", &StdI->Jz1, StdI->Jz);
  StdFace_PrintVal_d("Jz2", &StdI->Jz2, StdI->Jz);
  StdFace_PrintVal_d("Jxy", &StdI->Jxy, StdI->J);
  StdFace_PrintVal_d("Jxy0", &StdI->Jxy0, StdI->Jxy);
  StdFace_PrintVal_d("Jxy1", &StdI->Jxy1, StdI->Jxy);
  StdFace_PrintVal_d("Jxy2", &StdI->Jxy2, StdI->Jxy);
  StdFace_PrintVal_d("Jx", &StdI->Jx, StdI->Jxy);
  StdFace_PrintVal_d("Jx0", &StdI->Jx0, StdI->Jx);
  StdFace_PrintVal_d("Jx1", &StdI->Jx1, StdI->Jx);
  StdFace_PrintVal_d("Jx2", &StdI->Jx2, StdI->Jx);
  StdFace_PrintVal_d("Jy", &StdI->Jy, StdI->Jxy);
  StdFace_PrintVal_d("Jy0", &StdI->Jy0, StdI->Jy);
  StdFace_PrintVal_d("Jy1", &StdI->Jy1, StdI->Jy);
  StdFace_PrintVal_d("Jy2", &StdI->Jy2, StdI->Jy);
  StdI->Jxy = 0.5 * (StdI->Jx + StdI->Jy);
  StdI->Jxy0 = 0.5 * (StdI->Jx0 + StdI->Jy0);
  StdI->Jxy1 = 0.5 * (StdI->Jx1 + StdI->Jy1);
  StdI->Jxy2 = 0.5 * (StdI->Jx2 + StdI->Jy2);
  /**/
  StdFace_NotUsed_d("J'", StdI->Jp);
  StdFace_NotUsed_d("Jz'", StdI->Jzp);
  StdFace_NotUsed_d("Jxy'", StdI->Jxyp);
  StdFace_NotUsed_d("Jx'", StdI->Jxp);
  StdFace_NotUsed_d("Jy'", StdI->Jyp);
  StdFace_NotUsed_d("K", StdI->K);
  /*
  Local Spin
  */
  StdI->nsite = StdI->L * StdI->W * 2;
  StdI->S2 = 1;
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
      + 1.0 / 2.0 * S * S * (fabs(StdI->Jx0) + fabs(StdI->Jy0) + fabs(StdI->Jz0))
      + 1.0 / 2.0 * S * S * (fabs(StdI->Jx1) + fabs(StdI->Jy1) + fabs(StdI->Jz1))
      + 1.0 / 2.0 * S * S * (fabs(StdI->Jx2) + fabs(StdI->Jy2) + fabs(StdI->Jz2));
  }
  else {
    LargeValue0 = S * fabs(StdI->h) + S * fabs(StdI->D) + S * S * fabs(StdI->Gamma)
      + 1.0 / 2.0 * S * S * (fabs(StdI->Jx0) + fabs(StdI->Jy0) + fabs(StdI->Jz0))
      + 1.0 / 2.0 * S * S * (fabs(StdI->Jx1) + fabs(StdI->Jy1) + fabs(StdI->Jz1))
      + 1.0 / 2.0 * S * S * (fabs(StdI->Jx2) + fabs(StdI->Jy2) + fabs(StdI->Jz2));
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
  fprintf(fp, "%d  # Number of type of J\n", 3);
  fprintf(fp, "# J 1\n"); 
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->Jx0, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * StdI->Jy0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * StdI->Jz0, 0.0);
  fprintf(fp, "# J 2\n"); 
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->Jx1, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * StdI->Jy1, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * StdI->Jz1, 0.0);
  fprintf(fp, "# J 3\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->Jx2, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * StdI->Jy2, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * StdI->Jz2, 0.0);
  /*
   Topology
  */
  if (StdI->S2 != 1) {
    fprintf(stderr, "\n ERROR! S2 must be 1 in Boost. \n\n");
    exit(-1);
  }
  StdI->ishift_nspin = 3;
  if (StdI->L < 2) {
    fprintf(stderr, "\n ERROR! L < 2 \n\n");
    exit(-1);
  }
  if (StdI->W % StdI->ishift_nspin != 0) {
    fprintf(stderr, "\n ERROR! W %% %d != 0 \n\n", StdI->ishift_nspin);
    exit(-1);
  }
  StdI->num_pivot = 2;
  if (StdI->W != 3) {
    fprintf(stderr, "DEBUG: W != 3\n");
    exit(-1);
  }
  StdI->W = 6;
  fprintf(fp, "# W0  R0  StdI->num_pivot  StdI->ishift_nspin\n");
  fprintf(fp, "%d %d %d %d\n", StdI->W, StdI->L, StdI->num_pivot, StdI->ishift_nspin);

  StdI->list_6spin_star = (int **)malloc(sizeof(int*) * StdI->num_pivot);
  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    StdI->list_6spin_star[ipivot] = (int *)malloc(sizeof(int) * 7);
  }

  StdI->list_6spin_star[0][0] = 5; // num of J
  StdI->list_6spin_star[0][1] = 1;
  StdI->list_6spin_star[0][2] = 1;
  StdI->list_6spin_star[0][3] = 1;
  StdI->list_6spin_star[0][4] = 2;
  StdI->list_6spin_star[0][5] = 1;
  StdI->list_6spin_star[0][6] = 1; // flag

  StdI->list_6spin_star[1][0] = 4; //(0,2+2*j)=4 ! num of J
  StdI->list_6spin_star[1][1] = 1; //(1,2+2*j)=1
  StdI->list_6spin_star[1][2] = 1; //(2,2+2*j)=1
  StdI->list_6spin_star[1][3] = 1; //(3,2+2*j)=1
  StdI->list_6spin_star[1][4] = 2; //(4,2+2*j)=2
  StdI->list_6spin_star[1][5] = 2; //(5,2+2*j)=2
  StdI->list_6spin_star[1][6] = 1; //(6,2+2*j)=1 ! flag

  fprintf(fp, "# StdI->list_6spin_star\n");
  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    fprintf(fp, "# pivot %d\n", ipivot);
    for (isite = 0; isite < 7; isite++) {
      fprintf(fp, "%d ", StdI->list_6spin_star[ipivot][isite]);
    }
    fprintf(fp, "\n");
  }

  StdI->list_6spin_pair = (int ***)malloc(sizeof(int**) * StdI->num_pivot);
  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    StdI->list_6spin_pair[ipivot] = (int **)malloc(sizeof(int*) * 7);
    for (isite = 0; isite < 7; isite++) {
      StdI->list_6spin_pair[ipivot][isite] = (int *)malloc(sizeof(int) * StdI->list_6spin_star[ipivot][0]);
    }
  }

  StdI->list_6spin_pair[0][0][0] = 0; //(1,1,1+2*j)=0 
  StdI->list_6spin_pair[0][1][0] = 1; //(2,1,1+2*j)=1
  StdI->list_6spin_pair[0][2][0] = 2; //(3,1,1+2*j)=2
  StdI->list_6spin_pair[0][3][0] = 3; //(4,1,1+2*j)=3
  StdI->list_6spin_pair[0][4][0] = 4; //(5,1,1+2*j)=4
  StdI->list_6spin_pair[0][5][0] = 5; //(6,1,1+2*j)=5
  StdI->list_6spin_pair[0][6][0] = 3; //(7,1,1+2*j)=3 ! type of J
  StdI->list_6spin_pair[0][0][1] = 1; //(1,2,1+2*j)=1 
  StdI->list_6spin_pair[0][1][1] = 2; //(2,2,1+2*j)=2
  StdI->list_6spin_pair[0][2][1] = 0; //(3,2,1+2*j)=0
  StdI->list_6spin_pair[0][3][1] = 3; //(4,2,1+2*j)=3
  StdI->list_6spin_pair[0][4][1] = 4; //(5,2,1+2*j)=4
  StdI->list_6spin_pair[0][5][1] = 5; //(6,2,1+2*j)=5
  StdI->list_6spin_pair[0][6][1] = 1; //(7,2,1+2*j)=1 ! type of J
  StdI->list_6spin_pair[0][0][2] = 2; //(1,3,1+2*j)=2 
  StdI->list_6spin_pair[0][1][2] = 3; //(2,3,1+2*j)=3
  StdI->list_6spin_pair[0][2][2] = 0; //(3,3,1+2*j)=0
  StdI->list_6spin_pair[0][3][2] = 1; //(4,3,1+2*j)=1
  StdI->list_6spin_pair[0][4][2] = 4; //(5,3,1+2*j)=4
  StdI->list_6spin_pair[0][5][2] = 5; //(6,3,1+2*j)=5
  StdI->list_6spin_pair[0][6][2] = 3; //(7,3,1+2*j)=3 ! type of J
  StdI->list_6spin_pair[0][0][3] = 0; //(1,4,1+2*j)=0 
  StdI->list_6spin_pair[0][1][3] = 4; //(2,4,1+2*j)=4
  StdI->list_6spin_pair[0][2][3] = 1; //(3,4,1+2*j)=1
  StdI->list_6spin_pair[0][3][3] = 2; //(4,4,1+2*j)=2
  StdI->list_6spin_pair[0][4][3] = 3; //(5,4,1+2*j)=3
  StdI->list_6spin_pair[0][5][3] = 5; //(6,4,1+2*j)=5
  StdI->list_6spin_pair[0][6][3] = 1; //(7,4,1+2*j)=1 ! type of J
  StdI->list_6spin_pair[0][0][4] = 1; //(1,5,1+2*j)=1 
  StdI->list_6spin_pair[0][1][4] = 5; //(2,5,1+2*j)=5
  StdI->list_6spin_pair[0][2][4] = 0; //(3,5,1+2*j)=0
  StdI->list_6spin_pair[0][3][4] = 2; //(4,5,1+2*j)=2
  StdI->list_6spin_pair[0][4][4] = 3; //(5,5,1+2*j)=3
  StdI->list_6spin_pair[0][5][4] = 4; //(6,5,1+2*j)=4
  StdI->list_6spin_pair[0][6][4] = 2; //(7,5,1+2*j)=2 ! type of J

  StdI->list_6spin_pair[1][0][0] = 0; //(1,1,2+2*j)=0 
  StdI->list_6spin_pair[1][1][0] = 1; //(2,1,2+2*j)=1
  StdI->list_6spin_pair[1][2][0] = 2; //(3,1,2+2*j)=2
  StdI->list_6spin_pair[1][3][0] = 3; //(4,1,2+2*j)=3
  StdI->list_6spin_pair[1][4][0] = 4; //(5,1,2+2*j)=4
  StdI->list_6spin_pair[1][5][0] = 5; //(6,1,2+2*j)=5
  StdI->list_6spin_pair[1][6][0] = 1; //(7,1,2+2*j)=1 ! type of J
  StdI->list_6spin_pair[1][0][1] = 1; //(1,2,2+2*j)=1 
  StdI->list_6spin_pair[1][1][1] = 2; //(2,2,2+2*j)=2
  StdI->list_6spin_pair[1][2][1] = 0; //(3,2,2+2*j)=0
  StdI->list_6spin_pair[1][3][1] = 3; //(4,2,2+2*j)=3
  StdI->list_6spin_pair[1][4][1] = 4; //(5,2,2+2*j)=4
  StdI->list_6spin_pair[1][5][1] = 5; //(6,2,2+2*j)=5
  StdI->list_6spin_pair[1][6][1] = 3; //(7,2,2+2*j)=3 ! type of J
  StdI->list_6spin_pair[1][0][2] = 0; //(1,3,2+2*j)=0 
  StdI->list_6spin_pair[1][1][2] = 4; //(2,3,2+2*j)=4
  StdI->list_6spin_pair[1][2][2] = 1; //(3,3,2+2*j)=1
  StdI->list_6spin_pair[1][3][2] = 2; //(4,3,2+2*j)=2
  StdI->list_6spin_pair[1][4][2] = 3; //(5,3,2+2*j)=3
  StdI->list_6spin_pair[1][5][2] = 5; //(6,3,2+2*j)=5
  StdI->list_6spin_pair[1][6][2] = 2; //(7,3,2+2*j)=2 ! type of J
  StdI->list_6spin_pair[1][0][3] = 2; //(1,4,2+2*j)=2 
  StdI->list_6spin_pair[1][1][3] = 5; //(2,4,2+2*j)=5
  StdI->list_6spin_pair[1][2][3] = 0; //(3,4,2+2*j)=0
  StdI->list_6spin_pair[1][3][3] = 1; //(4,4,2+2*j)=1
  StdI->list_6spin_pair[1][4][3] = 3; //(5,4,2+2*j)=3
  StdI->list_6spin_pair[1][5][3] = 4; //(6,4,2+2*j)=4
  StdI->list_6spin_pair[1][6][3] = 2; //(7,4,2+2*j)=2 ! type of J

  fprintf(fp, "# StdI->list_6spin_pair\n");
  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    fprintf(fp, "# pivot %d\n", ipivot);
    for (kintr = 0; kintr < StdI->list_6spin_star[ipivot][0]; kintr++) {
      for (isite = 0; isite < 7; isite++) {
        fprintf(fp, "%d ", StdI->list_6spin_pair[ipivot][isite][kintr]);
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);

  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    free(StdI->list_6spin_star[ipivot]);
  }
  free(StdI->list_6spin_star);

  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    for (isite = 0; isite < 7; isite++) {
      free(StdI->list_6spin_pair[ipivot][isite]);
    }
    free(StdI->list_6spin_pair[ipivot]);
  }
  free(StdI->list_6spin_pair);

}

/**
 *
 * Setup a Hamiltonian for the Kondo lattice model on a Honeycomb lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void KondoLattice_HoneycombLattice(
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
  StdI->nsite = 4 * StdI->L * StdI->W;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  for (iL = 0; iL < 2 * StdI->L * StdI->W; iL++){
    StdI->locspinflag[iL] = StdI->S2;
    StdI->locspinflag[2 * StdI->L * StdI->W + iL] = 0;
   }
  /*
  Transfer
  */
  StdI->ntrans = StdI->L * StdI->W * 2 * 2 * 4;
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double *)malloc(sizeof(double) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    StdI->transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  StdI->ntrans = 0;
  for (iW = 0; iW < StdI->W; iW++){
    for (iL = 0; iL < StdI->L; iL++){
      for (ispin = 0; ispin < 2; ispin++){
        isite = 2 * StdI->L * StdI->W + (iL + iW * StdI->L) * 2;
        StdFace_trans(StdI, StdI->mu, isite, ispin, isite, ispin);
        jsite = 2 * StdI->L * StdI->W + (iL + iW * StdI->L) * 2 + 1;
        StdFace_trans(StdI, StdI->t0, isite, ispin, jsite, ispin);
        jsite = 2 * StdI->L * StdI->W + (((iL - 1 + 2 * StdI->L) % StdI->L) + ((iW + 0 + 2 * StdI->W) % StdI->W) * StdI->L) * 2 + 1;
        StdFace_trans(StdI, StdI->t1, isite, ispin, jsite, ispin);
        jsite = 2 * StdI->L * StdI->W + (((iL + 0 + 2 * StdI->L) % StdI->L) + ((iW - 1 + 2 * StdI->W) % StdI->W) * StdI->L) * 2 + 1;
        StdFace_trans(StdI, StdI->t2, isite, ispin, jsite, ispin);

        isite = 2 * StdI->L * StdI->W + (iL + iW * StdI->L) * 2 + 1;
        StdFace_trans(StdI, StdI->mu, isite, ispin, isite, ispin);
        jsite = 2 * StdI->L * StdI->W + (iL + iW * StdI->L) * 2;
        StdFace_trans(StdI, StdI->t0, isite, ispin, jsite, ispin);
        jsite = 2 * StdI->L * StdI->W + (((iL + 1 + 2 * StdI->L) % StdI->L) + ((iW + 0 + 2 * StdI->W) % StdI->W) * StdI->L) * 2;
        StdFace_trans(StdI, StdI->t1, isite, ispin, jsite, ispin);
        jsite = 2 * StdI->L * StdI->W + (((iL + 0 + 2 * StdI->L) % StdI->L) + ((iW + 1 + 2 * StdI->W) % StdI->W) * StdI->L) * 2;
        StdFace_trans(StdI, StdI->t2, isite, ispin, jsite, ispin);
      }
    }
  }
  /*
  Interaction
  */
  StdI->nintr = StdI->L * StdI->W * 2 * ((StdI->S2 + 1) * (1 + 1) + 2 * StdI->S2 * 1);
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++){
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  StdI->nintr = 0;
  for (iW = 0; iW < StdI->W; iW++){
    for (iL = 0; iL < StdI->L; iL++){
      isite = 2 * (iL + iW * StdI->L) + 2 * StdI->L * StdI->W;
      jsite = 2 * (iL + iW * StdI->L);
      StdFace_exchange(StdI, 1, StdI->S2, StdI->J, isite, jsite);
      StdFace_SzSz(StdI, 1, StdI->S2, StdI->J, isite, jsite);

      isite = 2 * (iL + iW * StdI->L) + 1 + 2 * StdI->L * StdI->W;
      jsite = 2 * (iL + iW * StdI->L) + 1;
      StdFace_exchange(StdI, 1, StdI->S2, StdI->J, isite, jsite);
      StdFace_SzSz(StdI, 1, StdI->S2, StdI->J, isite, jsite);
    }
  }
  /*
  Set mTPQ parameter
  */
  S = (double)StdI->S2 * 0.5;
  if (lGC == 0){
    LargeValue0 = fabs(StdI->mu) * (double)nelec / (double)(2 * StdI->L * StdI->W) + 2.0 * 3.0 * fabs(StdI->t) + 0.5 * S * fabs(StdI->J);
  }
  else{
    LargeValue0 = fabs(StdI->mu) * 2.0 + 2.0 * 3.0 * fabs(StdI->t) + 0.5 * S * fabs(StdI->J);
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
}
