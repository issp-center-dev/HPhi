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
void FermionHubbard_SquareLattice(
struct StdIntList *StdI,
  int nelec /**< [in] The number of electrons */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int ispin;
  int iL, iW, iL2, iW2;
  int ktrans, kintr;
  double trans0, LargeValue0;
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
  StdFace_PrintVal_d("V", &StdI->V, 0.0);
  StdFace_PrintVal_d("t'", &StdI->tp, 0.0);
  StdFace_PrintVal_d("V'", &StdI->Vp, 0.0);
  /**/
  StdFace_NotUsed_i("2S", StdI->S2);
  StdFace_NotUsed_d("t''", StdI->tpp);
  StdFace_NotUsed_d("t0", StdI->t0);
  StdFace_NotUsed_d("t1", StdI->t1);
  StdFace_NotUsed_d("t2", StdI->t2);
  StdFace_NotUsed_d("V''", StdI->Vpp);
  StdFace_NotUsed_d("V0", StdI->V0);
  StdFace_NotUsed_d("V1", StdI->V1);
  StdFace_NotUsed_d("V2", StdI->V2);
  /*
  Local Spin
  */
  StdI->nsite = StdI->L * StdI->W;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = 0;
  /*
  Transfer
  */
  StdI->ntrans = StdI->L * StdI->W * 2 * 9;
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double *)malloc(sizeof(double) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    StdI->transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  StdI->ntrans = 0;
  for (iW = 0; iW < StdI->W; iW++) {
    for (iL = 0; iL < StdI->L; iL++){
      isite = iL + iW * StdI->L;
      for (ispin = 0; ispin < 2; ispin++){
        for (iW2 = -1; iW2 <= 1; iW2++){
          for (iL2 = -1; iL2 <= 1; iL2++){

            if (iW2 == 0 && iL2 == 0) trans0 = StdI->mu;
            else if (iW2 == 0 || iL2 == 0) trans0 = StdI->t;
            else trans0 = StdI->tp;

            jsite = ((iL + iL2 + 2 * StdI->L) % StdI->L) + ((iW + iW2 + 2 * StdI->W) % StdI->W) * StdI->L;
            StdFace_trans(StdI, trans0, isite, ispin, jsite, ispin);
          }
        }
      }
    }
  }
  /*
  Interaction
  */
  StdI->nintr = StdI->L * StdI->W * (1 + 4 * 4);
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++){
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  StdI->nintr = 0;
  for (iW = 0; iW < StdI->W; iW++) {

    iW2 = (iW + 1) % StdI->W;

    for (iL = 0; iL < StdI->L; iL++){

      iL2 = (iL + 1) % StdI->L;

      isite = iL + iW * StdI->L;
      StdFace_intr(StdI, StdI->U, isite, 0, isite, 0, isite, 1, isite, 1);

      isite = iL + iW * StdI->L;
      jsite = iL + iW2 * StdI->L;
      StdFace_Coulomb(StdI, StdI->V, isite, jsite);

      jsite = iL2 + iW * StdI->L;
      StdFace_Coulomb(StdI, StdI->V, isite, jsite);

      jsite = iL2 + iW2 * StdI->L;
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);

      isite = iL + iW2 * StdI->L;
      jsite = iL2 + iW * StdI->L;
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
     }
  }
  /*
  Set mTPQ parameter
  */
  if (lGC == 0){
    LargeValue0 =  fabs(StdI->mu) * (double)nelec / (double)(StdI->L * StdI->W)
      + 2.0 * 4.0 * fabs(StdI->t) + 2.0 * 4.0 * fabs(StdI->tp)
      + fabs(StdI->U) + 2.0 * 4.0 * fabs(StdI->V) + 2.0 * 4.0 * fabs(StdI->Vp);
  }
  else{
    LargeValue0 = fabs(StdI->mu) * 2.0 + 2.0 * 4.0 * fabs(StdI->t) + 2.0 * 4.0 * fabs(StdI->tp)
      + fabs(StdI->U) + 2.0 * 4.0 * fabs(StdI->V) + 2.0 * 4.0 * fabs(StdI->Vp);
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
}

/**
 *
 * Setup a Hamiltonian for the generalized Heisenberg model on a square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void Spin_SquareLattice(
struct StdIntList *StdI,
  int Sz2 /**< [in] 2 * Total Sz */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int iL, iW, iW2, iL2;
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
  StdI->nsite = StdI->L * StdI->W;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = StdI->S2;
  /*
  Transfer
  */
  StdI->ntrans = StdI->L * StdI->W * (StdI->S2 + 1 + 2* StdI->S2);
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
  StdI->nintr = StdI->L * StdI->W * 
    ((StdI->S2 + 1) * (StdI->S2 + 1) * 5 + 2 * StdI->S2 *StdI->S2 * 8);
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++){
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  StdI->nintr = 0;
  for (iW = 0; iW < StdI->W; iW++){

    iW2 = (iW + 1) % StdI->W;

    for (iL = 0; iL <StdI->L; iL++){

      iL2 = (iL + 1) % StdI->L;

      isite = iL + iW * StdI->L;
      StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->D, isite, isite);

      jsite = iL + iW2 * StdI->L;
      StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->Jz, isite, jsite);
      StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->Jxy, isite, jsite);
      StdFace_Kitaev(StdI, StdI->S2, StdI->S2, 0.5 * (StdI->Jx - StdI->Jy), isite, jsite);

      jsite = iL2 + iW * StdI->L;
      StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->Jz, isite, jsite);
      StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->Jxy, isite, jsite);
      StdFace_Kitaev(StdI, StdI->S2, StdI->S2, 0.5 * (StdI->Jx - StdI->Jy), isite, jsite);

      jsite = iL2 + iW2 * StdI->L;
      StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->Jzp, isite, jsite);
      StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->Jxyp, isite, jsite);
      StdFace_Kitaev(StdI, StdI->S2, StdI->S2, 0.5 * (StdI->Jxp - StdI->Jyp), isite, jsite);

      isite = iL + iW2 * StdI->L;
      jsite = iL2 + iW * StdI->L;
      StdFace_SzSz(StdI, StdI->S2, StdI->S2, StdI->Jzp, isite, jsite);
      StdFace_exchange(StdI, StdI->S2, StdI->S2, StdI->Jxyp, isite, jsite);
      StdFace_Kitaev(StdI, StdI->S2, StdI->S2, 0.5 * (StdI->Jxp - StdI->Jyp), isite, jsite);
    }
  }
  /*
  Set mTPQ parameter
  */
  S = (double)StdI->S2 * 0.5;
  if (lGC == 0){
    LargeValue0 = (double)Sz2 / (double)(2 * StdI->nsite) * fabs(StdI->h) + S * fabs(StdI->D) + 0.5 * S * fabs(StdI->Gamma)
      + 4.0 / 2.0 * S * S * (fabs(StdI->Jx) + fabs(StdI->Jy) + fabs(StdI->Jz))
      + 4.0 / 2.0 * S * S * (fabs(StdI->Jxp) + fabs(StdI->Jyp) + fabs(StdI->Jzp));
  }
  else{
    LargeValue0 = S * fabs(StdI->h) + S * fabs(StdI->D) + 0.5 * S * fabs(StdI->Gamma)
      + 4.0 / 2.0 * S * S * (fabs(StdI->Jx) + fabs(StdI->Jy) + fabs(StdI->Jz))
      + 4.0 / 2.0 * S * S * (fabs(StdI->Jxp) + fabs(StdI->Jyp) + fabs(StdI->Jzp));
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
}

/**
 *
 * Setup a Hamiltonian for the Kondo lattice model on a square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void KondoLattice_SquareLattice(
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
  StdI->nsite = 2 * StdI->L * StdI->W;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  for (iL = 0; iL <StdI->L * StdI->W; iL++){
    StdI->locspinflag[iL] = StdI->S2;
    StdI->locspinflag[iL + StdI->L * StdI->W] = 0;
  }
  /*
  Transfer
  */
  StdI->ntrans = StdI->L * StdI->W * 2 * 5;
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double *)malloc(sizeof(double) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    StdI->transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  StdI->ntrans = 0;
  for (iW = 0; iW < StdI->W; iW++){
    for (iL = 0; iL < StdI->L; iL++){
      isite = StdI->L * StdI->W + iL + iW * StdI->L;
      for (ispin = 0; ispin < 2; ispin++){

        StdFace_trans(StdI, StdI->mu, isite, ispin, isite, ispin);

        jsite = StdI->L * StdI->W + ((iL + 1) % StdI->L) + (iW * StdI->L);
        StdFace_trans(StdI, StdI->t, isite, ispin, jsite, ispin);
        StdFace_trans(StdI, StdI->t, jsite, ispin, isite, ispin);

        jsite = StdI->L * StdI->W + iL + (((iW + 1) % StdI->W) * StdI->L);
        StdFace_trans(StdI, StdI->t, isite, ispin, jsite, ispin);
        StdFace_trans(StdI, StdI->t, jsite, ispin, isite, ispin);
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
  for (iW = 0; iW < StdI->W; iW++){
    for (iL = 0; iL < StdI->L; iL++){

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
    LargeValue0 = fabs(StdI->mu) * (double)nelec / (double)(StdI->L * StdI->W) + 2.0 * 4.0 * fabs(StdI->t) + 0.5 * S * fabs(StdI->J);
  }
  else{
    LargeValue0 = fabs(StdI->mu) * 2.0 + 2.0 * 4.0 * fabs(StdI->t) + 0.5 * S * fabs(StdI->J);
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
}
