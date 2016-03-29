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
#include <complex.h>
#include "../include/wrapperMPI.h"
#include <string.h>

/**
 *
 * Setup a Hamiltonian for the Hubbard model on a Honeycomb lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_Honeycomb(struct StdIntList *StdI, char *model)
{
  int isite, jsite, kCell, ispin;
  int iL, iW, iLp1;
  int ktrans, kintr;
  FILE *fp;
  /*
  Initialize Cell
  */
  fp = fopen("lattice.gp", "w");
  /**/
  StdI->NsiteUC = 2;
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_RequiredVal_i("W", StdI->W);
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  StdFace_PrintVal_d("a0", &StdI->a0, StdI->a);
  StdFace_PrintVal_d("a1", &StdI->a1, StdI->a);
  StdFace_PrintVal_d("Wx", &StdI->Wx, StdI->a0);
  StdFace_PrintVal_d("Wy", &StdI->Wy, 0.0);
  StdFace_PrintVal_d("Lx", &StdI->Lx, StdI->a1 * 0.5);
  StdFace_PrintVal_d("Ly", &StdI->Ly, StdI->a1 * 0.5 * sqrt(3.0));
  StdFace_PrintVal_i("a0W", &StdI->a0W, StdI->W);
  StdFace_PrintVal_i("a0L", &StdI->a0L, 0);
  StdFace_PrintVal_i("a1W", &StdI->a1W, 0);
  StdFace_PrintVal_i("a1L", &StdI->a1L, StdI->L);
  /**/
  StdFace_InitSite2D(StdI, fp);
  StdI->tau[0][0] = 0.0; StdI->tau[0][1] = 0.0;
  StdI->tau[1][0] = 1.0 / 3.0; StdI->tau[1][1] = 1.0 / 3.0;
  /**/
  if (strcmp(StdI->model, "spin") == 0 ) {
    StdFace_PrintVal_i("2S", &StdI->S2, 1);
    StdFace_PrintVal_d("h", &StdI->h, 0.0);
    StdFace_PrintVal_d("Gamma", &StdI->Gamma, 0.0);
    StdFace_PrintVal_d("D", &StdI->D[2][2], 0.0);
    /**/
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
    /**/
    StdFace_PrintVal_dd("J0yx", &StdI->J0[1][0], StdI->J0[0][1], StdI->J[1][0]);
    StdFace_PrintVal_dd("J0zx", &StdI->J0[2][0], StdI->J0[0][2], StdI->J[2][0]);
    StdFace_PrintVal_dd("J0zy", &StdI->J0[2][1], StdI->J0[1][2], StdI->J[2][1]);
    StdFace_PrintVal_dd("J0x", &StdI->J0[0][0], StdI->J0All, StdI->J[0][0]);
    StdFace_PrintVal_dd("J0y", &StdI->J0[1][1], StdI->J0All, StdI->J[1][1]);
    StdFace_PrintVal_dd("J0z", &StdI->J0[2][2], StdI->J0All, StdI->J[2][2]);
    StdFace_PrintVal_d("J0xy", &StdI->J0[0][1], StdI->J[0][1]);
    StdFace_PrintVal_d("J0xz", &StdI->J0[0][2], StdI->J[0][2]);
    StdFace_PrintVal_d("J0yz", &StdI->J0[1][2], StdI->J[1][2]);
    /**/
    StdFace_PrintVal_dd("J1yx", &StdI->J1[1][0], StdI->J1[0][1], StdI->J[1][0]);
    StdFace_PrintVal_dd("J1zx", &StdI->J1[2][0], StdI->J1[0][2], StdI->J[2][0]);
    StdFace_PrintVal_dd("J1zy", &StdI->J1[2][1], StdI->J1[1][2], StdI->J[2][1]);
    StdFace_PrintVal_dd("J1x", &StdI->J1[0][0], StdI->J1All, StdI->J[0][0]);
    StdFace_PrintVal_dd("J1y", &StdI->J1[1][1], StdI->J1All, StdI->J[1][1]);
    StdFace_PrintVal_dd("J1z", &StdI->J1[2][2], StdI->J1All, StdI->J[2][2]);
    StdFace_PrintVal_d("J1xy", &StdI->J1[0][1], StdI->J[0][1]);
    StdFace_PrintVal_d("J1xz", &StdI->J1[0][2], StdI->J[0][2]);
    StdFace_PrintVal_d("J1yz", &StdI->J1[1][2], StdI->J[1][2]);
    /**/
    StdFace_PrintVal_dd("J2yx", &StdI->J2[1][0], StdI->J2[0][1], StdI->J[1][0]);
    StdFace_PrintVal_dd("J2zx", &StdI->J2[2][0], StdI->J2[0][2], StdI->J[2][0]);
    StdFace_PrintVal_dd("J2zy", &StdI->J2[2][1], StdI->J2[1][2], StdI->J[2][1]);
    StdFace_PrintVal_dd("J2x", &StdI->J2[0][0], StdI->J2All, StdI->J[0][0]);
    StdFace_PrintVal_dd("J2y", &StdI->J2[1][1], StdI->J2All, StdI->J[1][1]);
    StdFace_PrintVal_dd("J2z", &StdI->J2[2][2], StdI->J2All, StdI->J[2][2]);
    StdFace_PrintVal_d("J2xy", &StdI->J2[0][1], StdI->J[0][1]);
    StdFace_PrintVal_d("J2xz", &StdI->J2[0][2], StdI->J[0][2]);
    StdFace_PrintVal_d("J2yz", &StdI->J2[1][2], StdI->J[1][2]);
    /**/
    StdFace_PrintVal_d("J'", &StdI->JpAll, 0.0);
    StdFace_PrintVal_d("J'x", &StdI->Jp[0][0], StdI->JpAll);
    StdFace_PrintVal_d("J'y", &StdI->Jp[1][1], StdI->JpAll);
    StdFace_PrintVal_d("J'z", &StdI->Jp[2][2], StdI->JpAll);
    StdFace_PrintVal_d("J'xy", &StdI->Jp[0][1], 0.0);
    StdFace_PrintVal_d("J'xz", &StdI->Jp[0][2], 0.0);
    StdFace_PrintVal_d("J'yz", &StdI->Jp[1][2], 0.0);
    StdFace_PrintVal_d("J'yx", &StdI->Jp[1][0], StdI->Jp[0][1]);
    StdFace_PrintVal_d("J'zx", &StdI->Jp[2][0], StdI->Jp[0][2]);
    StdFace_PrintVal_d("J'zy", &StdI->Jp[2][1], StdI->Jp[1][2]);
    /**/
    StdFace_NotUsed_d("K", StdI->K);
  }/*if (strcmp(StdI->model, "spin") == 0 )*/
  else {
    StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
    StdFace_PrintVal_d("U", &StdI->U, 0.0);
    StdFace_PrintVal_c("t", &StdI->t, 1.0);
    StdFace_PrintVal_c("t0", &StdI->t0, StdI->t);
    StdFace_PrintVal_c("t1", &StdI->t1, StdI->t);
    StdFace_PrintVal_c("t2", &StdI->t2, StdI->t);
    StdFace_PrintVal_d("V", &StdI->V, 0.0);
    StdFace_PrintVal_d("V0", &StdI->V0, StdI->V);
    StdFace_PrintVal_d("V1", &StdI->V1, StdI->V);
    StdFace_PrintVal_d("V2", &StdI->V2, StdI->V);
    StdFace_PrintVal_c("t'", &StdI->tp, 0.0);
    StdFace_PrintVal_d("V'", &StdI->Vp, 0.0);
    /**/

    if (strcmp(StdI->model, "hubbard") == 0 ) {
      StdFace_NotUsed_i("2S", StdI->S2);
      StdFace_NotUsed_d("J", StdI->JAll);
    }/*if (strcmp(StdI->model, "hubbard") == 0 )*/
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
    }/*if (model != "hubbard")*/

  }/*if (model != "spin")*/
   /*
   Local Spin
   */
  StdI->nsite = StdI->NsiteUC * StdI->NCell;
  if (strcmp(StdI->model, "kondo") == 0 ) StdI->nsite *= 2;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  /**/
  if (strcmp(StdI->model, "spin") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++) StdI->locspinflag[isite] = StdI->S2;
  else if (strcmp(StdI->model, "hubbard") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++) StdI->locspinflag[isite] = 0;
  else
    for (iL = 0; iL < StdI->nsite / 2; iL++) {
      StdI->locspinflag[iL] = StdI->S2;
      StdI->locspinflag[iL + StdI->nsite / 2] = 0;
    }
  /*
  The number of Transfer & Interaction
  */
  if (strcmp(StdI->model, "spin") == 0 ) {
    StdI->ntrans = 2 * StdI->NCell * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
    StdI->nintr = StdI->NCell * (2/*D*/ + 3/*J*/ + 6/*J'*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  else {
    StdI->ntrans = StdI->NCell * 2/*spin*/ * (2/*mu*/ + 6/*t*/ + 12/*t'*/);
    StdI->nintr = StdI->NCell * (2/*U*/ + 4 * (3/*V*/ + 6/*V'*/));

    if (strcmp(StdI->model, "kondo") == 0 )  StdI->nintr +=
      StdI->nsite / 2 * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  /**/
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double complex *)malloc(sizeof(double complex) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++) {
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
  for (kCell = 0; kCell < StdI->NCell; kCell++) {
    /**/
    iW = StdI->Cell[kCell][0];
    iL = StdI->Cell[kCell][1];
    /*
    Local term
    */
    isite = StdI->NsiteUC * kCell;
    if (strcmp(StdI->model, "kondo") == 0 ) isite += 2 * StdI->NCell;
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, isite);
      StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, isite + 1);
      StdFace_GeneralJ(StdI, StdI->D, StdI->S2, StdI->S2, isite, isite);
      StdFace_GeneralJ(StdI, StdI->D, StdI->S2, StdI->S2, isite + 1, isite + 1);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, StdI->mu, isite, isite);
      StdFace_Hopping(StdI, StdI->mu, isite + 1, isite + 1);
      StdFace_intr(StdI, StdI->U, isite, 0, isite, 0, isite, 1, isite, 1);
      StdFace_intr(StdI, StdI->U, isite + 1, 0, isite + 1, 0, isite + 1, 1, isite + 1, 1);
      /**/
      if (strcmp(StdI->model, "kondo") == 0 ) {
        jsite = StdI->NsiteUC * kCell;
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite, jsite);
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite + 1, jsite + 1);
      }/*if (strcmp(StdI->model, "kondo") == 0 )*/
    }
    /*
    Nearest neighbor intra cell 0 -> 1
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 0, 0, 0, 1, &isite, &jsite, 1);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J0, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, StdI->t0, isite, jsite);
      StdFace_Coulomb(StdI, StdI->V0, isite, jsite);
    }
    /*
    Nearest neighbor along W 1 -> 0
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 1, 0, 1, 0, &isite, &jsite, 1);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J1, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, StdI->t1, isite, jsite);
      StdFace_Coulomb(StdI, StdI->V1, isite, jsite);
    }
    /*
    Nearest neighbor along L 1 -> 0
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 0, 1, 1, 0, &isite, &jsite, 1);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J2, StdI->S2, StdI->S2, isite, jsite);
    }
    else {
      StdFace_Hopping(StdI, StdI->t2, isite, jsite);
      StdFace_Coulomb(StdI, StdI->V2, isite, jsite);
    }
    /*
    Second nearest neighbor along W 0 -> 0
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 1, 0, 0, 0, &isite, &jsite, 2);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, StdI->tp, isite, jsite);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
    /*
    Second nearest neighbor along W 1 -> 1
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 1, 0, 1, 1, &isite, &jsite, 2);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, StdI->tp, isite, jsite);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
    /*
    Second nearest neighbor along L 0 -> 0
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 0, 1, 0, 0, &isite, &jsite, 2);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, StdI->tp, isite, jsite);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
    /*
    Second nearest neighbor along L 1 -> 1
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 0, 1, 1, 1, &isite, &jsite, 2);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, StdI->tp, isite, jsite);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
    /*
    Second nearest neighbor along W-L 0 -> 0
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 1, - 1, 0, 0, &isite, &jsite, 2);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, StdI->tp, isite, jsite);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
    /*
    Second nearest neighbor along W - L 1 -> 1
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 1, -1, 1, 1, &isite, &jsite, 2);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, StdI->tp, isite, jsite);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
  }/*for (kCell = 0; kCell < StdI->NCell; kCell++)*/

  fprintf(fp, "plot x w l lw 0\n");
  fprintf(fp, "pause -1\n");
  fclose(fp);
}

/**
*
* Setup a Hamiltonian for the generalized Heisenberg model on a Heisenberg lattice
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_Honeycomb_Boost(struct StdIntList *StdI)
{
  int isite, ipivot;
  int ktrans, kintr;
  double LargeValue0, S;
  FILE *fp;

  StdI->NsiteUC = 2;
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
      0.25 * StdI->J0[0][0], 0.0, 0.25 * StdI->J0[0][1], 0.0, 0.25 * StdI->J0[0][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J0[0][1], 0.0, 0.25 * StdI->J0[1][1], 0.0, 0.25 * StdI->J0[1][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J0[0][2], 0.0, 0.25 * StdI->J0[1][2], 0.0, 0.25 * StdI->J0[2][2], 0.0);
  fprintf(fp, "# J 2\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1[0][0], 0.0, 0.25 * StdI->J1[0][1], 0.0, 0.25 * StdI->J1[0][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1[0][1], 0.0, 0.25 * StdI->J1[1][1], 0.0, 0.25 * StdI->J1[1][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1[0][2], 0.0, 0.25 * StdI->J1[1][2], 0.0, 0.25 * StdI->J1[2][2], 0.0);
  fprintf(fp, "# J 3\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2[0][0], 0.0, 0.25 * StdI->J2[0][1], 0.0, 0.25 * StdI->J2[0][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2[0][1], 0.0, 0.25 * StdI->J2[1][1], 0.0, 0.25 * StdI->J2[1][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2[0][2], 0.0, 0.25 * StdI->J2[1][2], 0.0, 0.25 * StdI->J2[2][2], 0.0);
  /*
   Topology
  */
  if (StdI->S2 != 1) {
    fprintf(stderr, "\n ERROR! S2 must be 1 in Boost. \n\n");
    exitMPI(-1);
  }
  StdI->ishift_nspin = 3;
  if (StdI->L < 2) {
    fprintf(stderr, "\n ERROR! L < 2 \n\n");
    exitMPI(-1);
  }
  if (StdI->W % StdI->ishift_nspin != 0) {
    fprintf(stderr, "\n ERROR! W %% %d != 0 \n\n", StdI->ishift_nspin);
    exitMPI(-1);
  }
  StdI->num_pivot = 2;
  if (StdI->W != 3) {
    fprintf(stderr, "DEBUG: W != 3\n");
    exitMPI(-1);
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

