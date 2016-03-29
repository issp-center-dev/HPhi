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
#include <complex.h>
#include "../include/wrapperMPI.h"
#include <string.h>

/**
 *
 * Setup a Hamiltonian for the generalized Heisenberg model on a square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_Ladder(struct StdIntList *StdI, char *model)
{
  int isite, jsite, ksite;
  int iL, iW, ispin;
  int ktrans, kintr;

  StdI->NsiteUC = 1;
  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_RequiredVal_i("W", StdI->W);
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  /**/
  if (strcmp(StdI->model, "spin") == 0 ) {
    StdFace_PrintVal_i("2S", &StdI->S2, 1);
    StdFace_PrintVal_d("h", &StdI->h, 0.0);
    StdFace_PrintVal_d("Gamma", &StdI->Gamma, 0.0);
    StdFace_PrintVal_d("D", &StdI->D[2][2], 0.0);
    StdFace_PrintVal_d("J0", &StdI->J0All, 1.0);
    StdFace_PrintVal_d("J0x", &StdI->J0[0][0], StdI->J0All);
    StdFace_PrintVal_d("J0y", &StdI->J0[1][1], StdI->J0All);
    StdFace_PrintVal_d("J0z", &StdI->J0[2][2], StdI->J0All);
    StdFace_PrintVal_d("J0xy", &StdI->J0[0][1], 0.0);
    StdFace_PrintVal_d("J0xz", &StdI->J0[0][2], 0.0);
    StdFace_PrintVal_d("J0yz", &StdI->J0[1][2], 0.0);
    StdFace_PrintVal_d("J0yx", &StdI->J0[1][0], StdI->J0[0][1]);
    StdFace_PrintVal_d("J0zx", &StdI->J0[2][0], StdI->J0[0][2]);
    StdFace_PrintVal_d("J0zy", &StdI->J0[2][1], StdI->J0[1][2]);
    StdFace_PrintVal_d("J1", &StdI->J1All, 0.0);
    StdFace_PrintVal_d("J1x", &StdI->J1[0][0], StdI->J1All);
    StdFace_PrintVal_d("J1y", &StdI->J1[1][1], StdI->J1All);
    StdFace_PrintVal_d("J1z", &StdI->J1[2][2], StdI->J1All);
    StdFace_PrintVal_d("J1xy", &StdI->J1[0][1], 0.0);
    StdFace_PrintVal_d("J1xz", &StdI->J1[0][2], 0.0);
    StdFace_PrintVal_d("J1yz", &StdI->J1[1][2], 0.0);
    StdFace_PrintVal_d("J1yx", &StdI->J1[1][0], StdI->J1[0][1]);
    StdFace_PrintVal_d("J1zx", &StdI->J1[2][0], StdI->J1[0][2]);
    StdFace_PrintVal_d("J1zy", &StdI->J1[2][1], StdI->J1[1][2]);   
    StdFace_PrintVal_d("J1'", &StdI->J1pAll, 0.0);
    StdFace_PrintVal_d("J1'x", &StdI->J1p[0][0], StdI->J1pAll);
    StdFace_PrintVal_d("J1'y", &StdI->J1p[1][1], StdI->J1pAll);
    StdFace_PrintVal_d("J1'z", &StdI->J1p[2][2], StdI->J1pAll);
    StdFace_PrintVal_d("J1'xy", &StdI->J1p[0][1], 0.0);
    StdFace_PrintVal_d("J1'xz", &StdI->J1p[0][2], 0.0);
    StdFace_PrintVal_d("J1'yz", &StdI->J1p[1][2], 0.0);
    StdFace_PrintVal_d("J1'yx", &StdI->J1p[1][0], StdI->J1p[0][1]);
    StdFace_PrintVal_d("J1'zx", &StdI->J1p[2][0], StdI->J1p[0][2]);
    StdFace_PrintVal_d("J1'zy", &StdI->J1p[2][1], StdI->J1p[1][2]);
    StdFace_PrintVal_d("J2", &StdI->J2All, 0.0);
    StdFace_PrintVal_d("J2x", &StdI->J2[0][0], StdI->J2All);
    StdFace_PrintVal_d("J2y", &StdI->J2[1][1], StdI->J2All);
    StdFace_PrintVal_d("J2z", &StdI->J2[2][2], StdI->J2All);
    StdFace_PrintVal_d("J2xy", &StdI->J2[0][1], 0.0);
    StdFace_PrintVal_d("J2xz", &StdI->J2[0][2], 0.0);
    StdFace_PrintVal_d("J2yz", &StdI->J2[1][2], 0.0);
    StdFace_PrintVal_d("J2yx", &StdI->J2[1][0], StdI->J2[0][1]);
    StdFace_PrintVal_d("J2zx", &StdI->J2[2][0], StdI->J2[0][2]);
    StdFace_PrintVal_d("J2zy", &StdI->J2[2][1], StdI->J2[1][2]);    
    StdFace_PrintVal_d("J2'", &StdI->J2pAll, 0.0);
    StdFace_PrintVal_d("J2'x", &StdI->J2p[0][0], StdI->J2pAll);
    StdFace_PrintVal_d("J2'y", &StdI->J2p[1][1], StdI->J2pAll);
    StdFace_PrintVal_d("J2'z", &StdI->J2p[2][2], StdI->J2pAll);
    StdFace_PrintVal_d("J2'xy", &StdI->J2p[0][1], 0.0);
    StdFace_PrintVal_d("J2'xz", &StdI->J2p[0][2], 0.0);
    StdFace_PrintVal_d("J2'yz", &StdI->J2p[1][2], 0.0);
    StdFace_PrintVal_d("J2'yx", &StdI->J2p[1][0], StdI->J2p[0][1]);
    StdFace_PrintVal_d("J2'zx", &StdI->J2p[2][0], StdI->J2p[0][2]);
    StdFace_PrintVal_d("J2'zy", &StdI->J2p[2][1], StdI->J2p[1][2]);
    /**/
    StdFace_NotUsed_d("J", StdI->JAll);
    StdFace_NotUsed_d("J'", StdI->JpAll);
    StdFace_NotUsed_d("K", StdI->K);
  }/*if (strcmp(StdI->model, "spin") == 0 )*/
  else {
    StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
    StdFace_PrintVal_d("U", &StdI->U, 0.0);
    StdFace_PrintVal_c("t1", &StdI->t1, 0.0);
    StdFace_PrintVal_c("t1'", &StdI->t1p, 0.0);
    StdFace_PrintVal_c("t0", &StdI->t0, 1.0);
    StdFace_PrintVal_c("t2", &StdI->t2, 0.0);
    StdFace_PrintVal_c("t2", &StdI->t2p, 0.0);
    StdFace_PrintVal_d("V0", &StdI->V0, 0.0);
    StdFace_PrintVal_d("V1", &StdI->V1, 0.0);
    StdFace_PrintVal_d("V1'", &StdI->V1p, 0.0);
    StdFace_PrintVal_d("V2", &StdI->V2, 0.0);
    StdFace_PrintVal_d("V2'", &StdI->V2p, 0.0);
    /**/
    StdFace_NotUsed_c("t", StdI->t);
    StdFace_NotUsed_c("t'", StdI->tp);
    StdFace_NotUsed_d("V", StdI->V);
    StdFace_NotUsed_d("V'", StdI->Vp);

    if (strcmp(StdI->model, "hubbard") == 0 ) {
      StdFace_NotUsed_i("2S", StdI->S2);
      StdFace_NotUsed_d("J", StdI->JAll);
    }
    else {
      StdFace_PrintVal_i("2S", &StdI->S2, 1);
      StdFace_PrintVal_d("J", &StdI->JAll, 1.0);
      StdFace_PrintVal_d("Jx", &StdI->J[0][0], StdI->JAll);
      StdFace_PrintVal_d("Jy", &StdI->J[1][1], StdI->JAll);
      StdFace_PrintVal_d("Jz", &StdI->J[2][2], StdI->JAll);
      StdFace_PrintVal_d("Jxy", &StdI->J[0][1], 0.0);
      StdFace_PrintVal_d("Jxz", &StdI->J[0][2], 0.0);
      StdFace_PrintVal_d("Jyz", &StdI->J[1][2], 0.0);
      StdFace_PrintVal_d("Jyx", &StdI->J[1][0], StdI->J[0][1]);
      StdFace_PrintVal_d("Jzx", &StdI->J[2][0], StdI->J[0][2]);
      StdFace_PrintVal_d("Jzy", &StdI->J[2][1], StdI->J[1][2]);
    }
  }/*if (model != "spin")*/
  /*
  Local Spin
  */
  StdI->nsite = StdI->L * StdI->W;
  if (strcmp(StdI->model, "kondo") == 0 ) StdI->nsite *= 2;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  /**/
  if (strcmp(StdI->model, "spin") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = StdI->S2;
  else if (strcmp(StdI->model, "hubbard") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = 0;
  else if (strcmp(StdI->model, "kondo") == 0 )
    for (isite = 0; isite < StdI->nsite / 2; isite++) {
      StdI->locspinflag[isite] = StdI->S2;
      StdI->locspinflag[isite + StdI->nsite / 2] = 0;
    }
  /*
  The number of Transfer & Interaction
  */
  if (strcmp(StdI->model, "spin") == 0 ) {
    StdI->ntrans = StdI->L * StdI->W * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
    StdI->nintr = StdI->L * StdI->W * (1/*D*/ + 1/*J1*/ + 1/*J1'*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1)
      + StdI->L * (StdI->W - 1) * (1/*J0*/ + 1/*J2*/ + 1/*J2'*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }/*if (strcmp(StdI->model, "spin") == 0 )*/
  else {
    StdI->ntrans = StdI->L*StdI->W * (1/*mu*/ + 2/*t1*/ + 2/*t1'*/)
      + StdI->L*(StdI->W - 1) * (2/*t0*/ + 2/*t2*/ + 2/*t2'*/);
    StdI->nintr = StdI->L*StdI->W * 1/*U*/
      + StdI->L*StdI->W * 4 * (1/*V1*/ + 1/*V1'*/)
      + StdI->L*(StdI->W - 1) * 4 * (1/*V0*/ + 1/*V2*/ + 1/*V2'*/);

    if (strcmp(StdI->model, "kondo") == 0 ) StdI->nintr += StdI->nsite / 2 * (3 * 1 + 1) * (3 * StdI->S2 + 1);
  }
  /**/
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double complex *)malloc(sizeof(double complex) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    StdI->transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }
  /**/
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double complex *)malloc(sizeof(double complex) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++) {
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  /*
   Set Transfer & Interaction
  */
  StdI->ntrans = 0;
  StdI->nintr = 0;
  for (iL = 0; iL < StdI->L; iL++) {
    for (iW = 0; iW < StdI->W; iW++) {

      isite = iW + iL * StdI->W;
      if (strcmp(StdI->model, "kondo") == 0 ) isite += StdI->L * StdI->W;
      /*
       Local term
      */
      if (strcmp(StdI->model, "spin") == 0 ) {
        StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, isite);
        StdFace_GeneralJ(StdI, StdI->D, StdI->S2, StdI->S2, isite, jsite);
      }/*if (strcmp(StdI->model, "spin") == 0 )*/
      else {
        StdFace_Hopping(StdI, StdI->mu, isite, isite);
        StdFace_intr(StdI, StdI->U, isite, 0, isite, 0, isite, 1, isite, 1);
        /**/
        if (strcmp(StdI->model, "kondo") == 0 ) {
          jsite = iW + iL * StdI->W;
          StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite, jsite);
        }/*if (strcmp(StdI->model, "kondo") == 0 )*/
      }/*if (model != "spin")*/
      /*
       Nearest neighbor along the ladder
      */
      jsite = iW + ((iL + 1) % StdI->L) * StdI->W;
      if (strcmp(StdI->model, "kondo") == 0 ) jsite += StdI->L * StdI->W;
      /**/
      if (strcmp(StdI->model, "spin") == 0 ) {
        StdFace_GeneralJ(StdI, StdI->J1, StdI->S2, StdI->S2, isite, jsite);
      }/*if (strcmp(StdI->model, "spin") == 0 )*/
      else {
        StdFace_Hopping(StdI, StdI->t1, isite, jsite);
        StdFace_Coulomb(StdI, StdI->V1, isite, jsite);
      }/*if (model != "spin")*/
      /*
       Second nearest neighbor along the ladder
      */
      jsite = iW + ((iL + 2) % StdI->L) * StdI->W;
      if (strcmp(StdI->model, "kondo") == 0 ) jsite += StdI->L * StdI->W;
      if (strcmp(StdI->model, "spin") == 0 ) {
        StdFace_GeneralJ(StdI, StdI->J1p, StdI->S2, StdI->S2, isite, jsite);
      }/*if (strcmp(StdI->model, "spin") == 0 )*/
      else {
        StdFace_Hopping(StdI, StdI->t1p, isite, jsite);
        StdFace_Coulomb(StdI, StdI->V1p, isite, jsite);
      }/*if (model != "spin")*/
      /*
      Across rung
      */
      if (iW < StdI->W - 1) {
        /*
         Vertical
        */
        jsite = (iW + 1) + iL * StdI->W;
        if (strcmp(StdI->model, "kondo") == 0 ) jsite += StdI->L * StdI->W;
        /**/
        if (strcmp(StdI->model, "spin") == 0 ) {
          StdFace_GeneralJ(StdI, StdI->J0, StdI->S2, StdI->S2, isite, jsite);
        }/*if (strcmp(StdI->model, "spin") == 0 )*/
        else {
          StdFace_Hopping(StdI, StdI->t0, isite, jsite);
          StdFace_Coulomb(StdI, StdI->V0, isite, jsite);
        }/*if (model != "spin")*/
        /*
         Diagonal 1
        */
        jsite = (iW + 1) + ((iL + 1) % StdI->L) * StdI->W;
        if (strcmp(StdI->model, "kondo") == 0 ) jsite += StdI->L * StdI->W;
        /**/
        if (strcmp(StdI->model, "spin") == 0 ) {
          StdFace_GeneralJ(StdI, StdI->J2, StdI->S2, StdI->S2, isite, jsite);
        }/*if (strcmp(StdI->model, "spin") == 0 )*/
        else {
          StdFace_Hopping(StdI, StdI->t2, isite, jsite);
          StdFace_Coulomb(StdI, StdI->V2, isite, jsite);
        }/*if (model != "spin")*/
        /*
         Diagonal 2
        */
        /**/
        if (strcmp(StdI->model, "spin") == 0 ) {
        }
        ksite = (iW + 1) + iL * StdI->W;
        if (strcmp(StdI->model, "kondo") == 0 ) ksite += StdI->L * StdI->W;
        jsite = iW + ((iL + 1) % StdI->L) * StdI->W;
        if (strcmp(StdI->model, "kondo") == 0 ) jsite += StdI->L * StdI->W;
        if (strcmp(StdI->model, "spin") == 0 ) {
          StdFace_GeneralJ(StdI, StdI->J2p, StdI->S2, StdI->S2, isite, jsite);
        }/*if (strcmp(StdI->model, "spin") == 0 )*/
        else {
          StdFace_Hopping(StdI, StdI->t2p, isite, jsite);
          StdFace_Coulomb(StdI, StdI->V2p, isite, jsite);
        }/*if (model != "spin")*/

      }/*if (iW < StdI->W - 1)*/

    }/*for (iW = 0; iW < StdI->W; iW++)*/
  }/*for (iL = 0; iL < StdI->L; iL++)*/

}

/**
*
* Setup a Hamiltonian for the generalized Heisenberg model on a square lattice
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_Ladder_Boost(struct StdIntList *StdI)
{
  int isite, ipivot;
  int ktrans, kintr;
  double LargeValue0, S;
  FILE *fp;

  StdI->NsiteUC = 1;
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
  StdFace_PrintVal_d("J0", &StdI->J0All, 1.0);
  StdFace_PrintVal_d("J0x", &StdI->J0[0][0], StdI->J0All);
  StdFace_PrintVal_d("J0y", &StdI->J0[1][1], StdI->J0All);
  StdFace_PrintVal_d("J0z", &StdI->J0[2][2], StdI->J0All);
  StdFace_PrintVal_d("J0xy", &StdI->J0[0][1], 0.0);
  StdFace_PrintVal_d("J0xz", &StdI->J0[0][2], 0.0);
  StdFace_PrintVal_d("J0yz", &StdI->J0[1][2], 0.0);
  StdFace_PrintVal_d("J0yx", &StdI->J0[1][0], StdI->J0[0][1]);
  StdFace_PrintVal_d("J0zx", &StdI->J0[2][0], StdI->J0[0][2]);
  StdFace_PrintVal_d("J0zy", &StdI->J0[2][1], StdI->J0[1][2]);
  StdFace_PrintVal_d("J1", &StdI->J1All, 0.0);
  StdFace_PrintVal_d("J1x", &StdI->J1[0][0], StdI->J1All);
  StdFace_PrintVal_d("J1y", &StdI->J1[1][1], StdI->J1All);
  StdFace_PrintVal_d("J1z", &StdI->J1[2][2], StdI->J1All);
  StdFace_PrintVal_d("J1xy", &StdI->J1[0][1], 0.0);
  StdFace_PrintVal_d("J1xz", &StdI->J1[0][2], 0.0);
  StdFace_PrintVal_d("J1yz", &StdI->J1[1][2], 0.0);
  StdFace_PrintVal_d("J1yx", &StdI->J1[1][0], StdI->J1[0][1]);
  StdFace_PrintVal_d("J1zx", &StdI->J1[2][0], StdI->J1[0][2]);
  StdFace_PrintVal_d("J1zy", &StdI->J1[2][1], StdI->J1[1][2]);
  StdFace_PrintVal_d("J1'", &StdI->J1pAll, 0.0);
  StdFace_PrintVal_d("J1'x", &StdI->J1p[0][0], StdI->J1pAll);
  StdFace_PrintVal_d("J1'y", &StdI->J1p[1][1], StdI->J1pAll);
  StdFace_PrintVal_d("J1'z", &StdI->J1p[2][2], StdI->J1pAll);
  StdFace_PrintVal_d("J1'xy", &StdI->J1p[0][1], 0.0);
  StdFace_PrintVal_d("J1'xz", &StdI->J1p[0][2], 0.0);
  StdFace_PrintVal_d("J1'yz", &StdI->J1p[1][2], 0.0);
  StdFace_PrintVal_d("J1'yx", &StdI->J1p[1][0], StdI->J1p[0][1]);
  StdFace_PrintVal_d("J1'zx", &StdI->J1p[2][0], StdI->J1p[0][2]);
  StdFace_PrintVal_d("J1'zy", &StdI->J1p[2][1], StdI->J1p[1][2]);
  StdFace_PrintVal_d("J2", &StdI->J2All, 0.0);
  StdFace_PrintVal_d("J2x", &StdI->J2[0][0], StdI->J2All);
  StdFace_PrintVal_d("J2y", &StdI->J2[1][1], StdI->J2All);
  StdFace_PrintVal_d("J2z", &StdI->J2[2][2], StdI->J2All);
  StdFace_PrintVal_d("J2xy", &StdI->J2[0][1], 0.0);
  StdFace_PrintVal_d("J2xz", &StdI->J2[0][2], 0.0);
  StdFace_PrintVal_d("J2yz", &StdI->J2[1][2], 0.0);
  StdFace_PrintVal_d("J2yx", &StdI->J2[1][0], StdI->J2[0][1]);
  StdFace_PrintVal_d("J2zx", &StdI->J2[2][0], StdI->J2[0][2]);
  StdFace_PrintVal_d("J2zy", &StdI->J2[2][1], StdI->J2[1][2]);
  StdFace_PrintVal_d("J2'", &StdI->J2pAll, 0.0);
  StdFace_PrintVal_d("J2'x", &StdI->J2p[0][0], StdI->J2pAll);
  StdFace_PrintVal_d("J2'y", &StdI->J2p[1][1], StdI->J2pAll);
  StdFace_PrintVal_d("J2'z", &StdI->J2p[2][2], StdI->J2pAll);
  StdFace_PrintVal_d("J2'xy", &StdI->J2p[0][1], 0.0);
  StdFace_PrintVal_d("J2'xz", &StdI->J2p[0][2], 0.0);
  StdFace_PrintVal_d("J2'yz", &StdI->J2p[1][2], 0.0);
  StdFace_PrintVal_d("J2'yx", &StdI->J2p[1][0], StdI->J2p[0][1]);
  StdFace_PrintVal_d("J2'zx", &StdI->J2p[2][0], StdI->J2p[0][2]);
  StdFace_PrintVal_d("J2'zy", &StdI->J2p[2][1], StdI->J2p[1][2]);
  /**/
  StdFace_NotUsed_d("D", StdI->D[2][2]);
  StdFace_NotUsed_d("J", StdI->JAll);
  StdFace_NotUsed_d("J'", StdI->JpAll);
  StdFace_NotUsed_d("K", StdI->K);
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
    0.25 * StdI->J0[0][0], 0.0, 0.25 * StdI->J0[0][1], 0.0, 0.25 * StdI->J0[0][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J0[0][1], 0.0, 0.25 * StdI->J0[1][1], 0.0, 0.25 * StdI->J0[1][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J0[0][2], 0.0, 0.25 * StdI->J0[1][2], 0.0, 0.25 * StdI->J0[2][2], 0.0);
  fprintf(fp, "# J 2 (Nearest neighbor, along chain)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1[0][0], 0.0, 0.25 * StdI->J1[0][1], 0.0, 0.25 * StdI->J1[0][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1[0][1], 0.0, 0.25 * StdI->J1[1][1], 0.0, 0.25 * StdI->J1[1][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1[0][2], 0.0, 0.25 * StdI->J1[1][2], 0.0, 0.25 * StdI->J1[2][2], 0.0);
  fprintf(fp, "# J 3 (Second nearest neighbor, along chain)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1p[0][0], 0.0, 0.25 * StdI->J1p[0][1], 0.0, 0.25 * StdI->J1p[0][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1p[0][1], 0.0, 0.25 * StdI->J1p[1][1], 0.0, 0.25 * StdI->J1p[1][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1p[0][2], 0.0, 0.25 * StdI->J1p[1][2], 0.0, 0.25 * StdI->J1p[2][2], 0.0);
  fprintf(fp, "# J 4 (inter chain, diagonal1)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2[0][0], 0.0, 0.25 * StdI->J2[0][1], 0.0, 0.25 * StdI->J2[0][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2[0][1], 0.0, 0.25 * StdI->J2[1][1], 0.0, 0.25 * StdI->J2[1][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2[0][2], 0.0, 0.25 * StdI->J2[1][2], 0.0, 0.25 * StdI->J2[2][2], 0.0);
  fprintf(fp, "# J 5 (inter chain, diagonal2)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2p[0][0], 0.0, 0.25 * StdI->J2p[0][1], 0.0, 0.25 * StdI->J2p[0][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2p[0][1], 0.0, 0.25 * StdI->J2p[1][1], 0.0, 0.25 * StdI->J2p[1][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2p[0][2], 0.0, 0.25 * StdI->J2p[1][2], 0.0, 0.25 * StdI->J2p[2][2], 0.0);
  /*
  Topology
  */
  if (StdI->S2 != 1) {
    fprintf(stderr, "\n ERROR! S2 must be 1 in Boost. \n\n");
    exitMPI(-1);
  }
  StdI->ishift_nspin = 2;
  if (StdI->W != 2) {
    fprintf(stderr, "\n ERROR! W != 2 \n\n");
    exitMPI(-1);
  }
  if (StdI->L % 2 != 0) {
    fprintf(stderr, "\n ERROR! L % 2 != 0 \n\n");
    exitMPI(-1);
  }
  if (StdI->L < 4) {
    fprintf(stderr, "\n ERROR! L < 4 \n\n");
    exitMPI(-1);
  }
  StdI->W = StdI->L;
  StdI->L = 2;
  StdI->num_pivot = StdI->W / 2;
  /**/
  fprintf(fp, "# W0  R0  StdI->num_pivot  StdI->ishift_nspin\n");
  fprintf(fp, "%d %d %d %d\n", StdI->W, StdI->L, StdI->num_pivot, StdI->ishift_nspin);

  StdI->list_6spin_star = (int **)malloc(sizeof(int*) * StdI->num_pivot);
  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    StdI->list_6spin_star[ipivot] = (int *)malloc(sizeof(int) * 7);
  }

  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    StdI->list_6spin_star[ipivot][0] = 7; // num of J
    StdI->list_6spin_star[ipivot][1] = 1;
    StdI->list_6spin_star[ipivot][2] = 1;
    StdI->list_6spin_star[ipivot][3] = 1;
    StdI->list_6spin_star[ipivot][4] = 1;
    StdI->list_6spin_star[ipivot][5] = 1;
    StdI->list_6spin_star[ipivot][6] = 1; // flag
  }

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

  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    StdI->list_6spin_pair[ipivot][0][0] = 0;
    StdI->list_6spin_pair[ipivot][1][0] = 1;
    StdI->list_6spin_pair[ipivot][2][0] = 2;
    StdI->list_6spin_pair[ipivot][3][0] = 3;
    StdI->list_6spin_pair[ipivot][4][0] = 4;
    StdI->list_6spin_pair[ipivot][5][0] = 5;
    StdI->list_6spin_pair[ipivot][6][0] = 1; // type of J
    StdI->list_6spin_pair[ipivot][0][1] = 0;
    StdI->list_6spin_pair[ipivot][1][1] = 2;
    StdI->list_6spin_pair[ipivot][2][1] = 1;
    StdI->list_6spin_pair[ipivot][3][1] = 3;
    StdI->list_6spin_pair[ipivot][4][1] = 4;
    StdI->list_6spin_pair[ipivot][5][1] = 5;
    StdI->list_6spin_pair[ipivot][6][1] = 2; // type of J
    StdI->list_6spin_pair[ipivot][0][2] = 1;
    StdI->list_6spin_pair[ipivot][1][2] = 3;
    StdI->list_6spin_pair[ipivot][2][2] = 0;
    StdI->list_6spin_pair[ipivot][3][2] = 2;
    StdI->list_6spin_pair[ipivot][4][2] = 4;
    StdI->list_6spin_pair[ipivot][5][2] = 5;
    StdI->list_6spin_pair[ipivot][6][2] = 2; // type of J
    StdI->list_6spin_pair[ipivot][0][3] = 0;
    StdI->list_6spin_pair[ipivot][1][3] = 4;
    StdI->list_6spin_pair[ipivot][2][3] = 1;
    StdI->list_6spin_pair[ipivot][3][3] = 2;
    StdI->list_6spin_pair[ipivot][4][3] = 3;
    StdI->list_6spin_pair[ipivot][5][3] = 5;
    StdI->list_6spin_pair[ipivot][6][3] = 3; // type of J
    StdI->list_6spin_pair[ipivot][0][4] = 1;
    StdI->list_6spin_pair[ipivot][1][4] = 5;
    StdI->list_6spin_pair[ipivot][2][4] = 0;
    StdI->list_6spin_pair[ipivot][3][4] = 2;
    StdI->list_6spin_pair[ipivot][4][4] = 3;
    StdI->list_6spin_pair[ipivot][5][4] = 4;
    StdI->list_6spin_pair[ipivot][6][4] = 3; // type of J
    StdI->list_6spin_pair[ipivot][0][5] = 0;
    StdI->list_6spin_pair[ipivot][1][5] = 3;
    StdI->list_6spin_pair[ipivot][2][5] = 1;
    StdI->list_6spin_pair[ipivot][3][5] = 2;
    StdI->list_6spin_pair[ipivot][4][5] = 4;
    StdI->list_6spin_pair[ipivot][5][5] = 5;
    StdI->list_6spin_pair[ipivot][6][5] = 4; // type of J
    StdI->list_6spin_pair[ipivot][0][6] = 1;
    StdI->list_6spin_pair[ipivot][1][6] = 2;
    StdI->list_6spin_pair[ipivot][2][6] = 0;
    StdI->list_6spin_pair[ipivot][3][6] = 3;
    StdI->list_6spin_pair[ipivot][4][6] = 4;
    StdI->list_6spin_pair[ipivot][5][6] = 5;
    StdI->list_6spin_pair[ipivot][6][6] = 5; // type of J
  }

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
