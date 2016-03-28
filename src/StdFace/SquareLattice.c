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
 * Setup a Hamiltonian for the square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_Tetragonal(struct StdIntList *StdI, char *model)
{
  int isite, jsite;
  int iL, iW, jL, jW, kCell;
  int ktrans, kintr;
  FILE *fp;

  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  /*
   Initialize Cell
  */
  fp = fopen("lattice.gp", "w");
  /**/
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_RequiredVal_i("W", StdI->W);
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  StdFace_PrintVal_d("a0", &StdI->a0, StdI->a);
  StdFace_PrintVal_d("a1", &StdI->a1, StdI->a);
  StdFace_PrintVal_d("Wx", &StdI->Wx, StdI->a0);
  StdFace_PrintVal_d("Wy", &StdI->Wy, 0.0);
  StdFace_PrintVal_d("Lx", &StdI->Lx, 0.0);
  StdFace_PrintVal_d("Ly", &StdI->Ly, StdI->a1);
  StdFace_PrintVal_i("a0W", &StdI->a0W, StdI->W);
  StdFace_PrintVal_i("a0L", &StdI->a0L, 0);
  StdFace_PrintVal_i("a1W", &StdI->a1W, 0);
  StdFace_PrintVal_i("a1L", &StdI->a1L, StdI->L);
  /**/
  StdFace_InitSite2D(StdI, fp);
  /**/
  if (model == "spin") {
    StdFace_PrintVal_i("2S", &StdI->S2, 1);
    StdFace_PrintVal_d("h", &StdI->h, 0.0);
    StdFace_PrintVal_d("Gamma", &StdI->Gamma, 0.0);
    StdFace_PrintVal_d("D", &StdI->D[2][2], 0.0);
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
    if (StdI->J0All > 9999.0) {
      StdFace_PrintVal_d("J0", &StdI->J0All, StdI->JAll);
      StdFace_PrintVal_d("J0x", &StdI->J0[0][0], StdI->J[0][0]);
      StdFace_PrintVal_d("J0y", &StdI->J0[1][1], StdI->J[1][1]);
      StdFace_PrintVal_d("J0z", &StdI->J0[2][2], StdI->J[2][2]);
      StdFace_PrintVal_d("J0xy", &StdI->J0[0][1], StdI->J[0][1]);
      StdFace_PrintVal_d("J0xz", &StdI->J0[0][2], StdI->J[0][2]);
      StdFace_PrintVal_d("J0yz", &StdI->J0[1][2], StdI->J[1][2]);
      StdFace_PrintVal_d("J0yx", &StdI->J0[1][0], StdI->J[1][0]);
      StdFace_PrintVal_d("J0zx", &StdI->J0[2][0], StdI->J[2][0]);
      StdFace_PrintVal_d("J0zy", &StdI->J0[2][1], StdI->J[2][1]);
    }
    else {
      StdFace_PrintVal_d("J0", &StdI->J0All, StdI->JAll);
      StdFace_PrintVal_d("J0x", &StdI->J0[0][0], StdI->J0All);
      StdFace_PrintVal_d("J0y", &StdI->J0[1][1], StdI->J0All);
      StdFace_PrintVal_d("J0z", &StdI->J0[2][2], StdI->J0All);
      StdFace_PrintVal_d("J0xy", &StdI->J0[0][1], 0.0);
      StdFace_PrintVal_d("J0xz", &StdI->J0[0][2], 0.0);
      StdFace_PrintVal_d("J0yz", &StdI->J0[1][2], 0.0);
      StdFace_PrintVal_d("J0yx", &StdI->J0[1][0], StdI->J0[0][1]);
      StdFace_PrintVal_d("J0zx", &StdI->J0[2][0], StdI->J0[0][2]);
      StdFace_PrintVal_d("J0zy", &StdI->J0[2][1], StdI->J0[1][2]);
    }
    if (StdI->J1All > 9999.0) {
      StdFace_PrintVal_d("J1", &StdI->J1All, StdI->JAll);
      StdFace_PrintVal_d("J1x", &StdI->J1[0][0], StdI->J[0][0]);
      StdFace_PrintVal_d("J1y", &StdI->J1[1][1], StdI->J[1][1]);
      StdFace_PrintVal_d("J1z", &StdI->J1[2][2], StdI->J[2][2]);
      StdFace_PrintVal_d("J1xy", &StdI->J1[0][1], StdI->J[0][1]);
      StdFace_PrintVal_d("J1xz", &StdI->J1[0][2], StdI->J[0][2]);
      StdFace_PrintVal_d("J1yz", &StdI->J1[1][2], StdI->J[1][2]);
      StdFace_PrintVal_d("J1yx", &StdI->J1[1][0], StdI->J[1][0]);
      StdFace_PrintVal_d("J1zx", &StdI->J1[2][0], StdI->J[2][0]);
      StdFace_PrintVal_d("J1zy", &StdI->J1[2][1], StdI->J[2][1]);
    }
    else {
      StdFace_PrintVal_d("J1", &StdI->J1All, StdI->JAll);
      StdFace_PrintVal_d("J1x", &StdI->J1[0][0], StdI->J1All);
      StdFace_PrintVal_d("J1y", &StdI->J1[1][1], StdI->J1All);
      StdFace_PrintVal_d("J1z", &StdI->J1[2][2], StdI->J1All);
      StdFace_PrintVal_d("J1xy", &StdI->J1[0][1], 0.0);
      StdFace_PrintVal_d("J1xz", &StdI->J1[0][2], 0.0);
      StdFace_PrintVal_d("J1yz", &StdI->J1[1][2], 0.0);
      StdFace_PrintVal_d("J1yx", &StdI->J1[1][0], StdI->J1[0][1]);
      StdFace_PrintVal_d("J1zx", &StdI->J1[2][0], StdI->J1[0][2]);
      StdFace_PrintVal_d("J1zy", &StdI->J1[2][1], StdI->J1[1][2]);
    }
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
    StdFace_NotUsed_d("J2", StdI->J2All);
    StdFace_NotUsed_d("K", StdI->K);
  }/*if (model == "spin")*/
  else {
    StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
    StdFace_PrintVal_d("U", &StdI->U, 0.0);
    StdFace_PrintVal_c("t", &StdI->t, 1.0);
    StdFace_PrintVal_c("t0", &StdI->t0, StdI->t);
    StdFace_PrintVal_c("t1", &StdI->t1, StdI->t);
    StdFace_PrintVal_d("V", &StdI->V, 0.0);
    StdFace_PrintVal_d("V0", &StdI->V0, StdI->V);
    StdFace_PrintVal_d("V1", &StdI->V1, StdI->V);
    StdFace_PrintVal_c("t'", &StdI->tp, 0.0);
    StdFace_PrintVal_d("V'", &StdI->Vp, 0.0);
    /**/
    StdFace_NotUsed_c("t2", StdI->t2);
    StdFace_NotUsed_d("V2", StdI->V2);

    if (model == "hubbard") {
      StdFace_NotUsed_i("2S", StdI->S2);
      StdFace_NotUsed_d("J", StdI->JAll);
    }/*if (model == "hubbard")*/
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
  StdI->nsite = StdI->L * StdI->W;
  if (model == "kondo") StdI->nsite *= 2;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  /**/
  if(model == "spin")
    for (isite = 0; isite < StdI->nsite; isite++) StdI->locspinflag[isite] = StdI->S2;
  else if(model == "hubbard")
    for (isite = 0; isite < StdI->nsite; isite++) StdI->locspinflag[isite] = 0;
  else 
    for (iL = 0; iL < StdI->nsite / 2; iL++) {
      StdI->locspinflag[iL] = StdI->S2;
      StdI->locspinflag[iL + StdI->L * StdI->W] = 0;
    }
  /*
   The number of Transfer & Interaction
  */
  if (model == "spin") {
    StdI->ntrans = StdI->L * StdI->W * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
    StdI->nintr = StdI->L * StdI->W * (1/*D*/ + 2/*J*/ + 2/*J'*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  else {
    StdI->ntrans = StdI->L * StdI->W * 2/*spin*/ * (1/*mu*/ + 4/*t*/ + 4/*t'*/);
    StdI->nintr = StdI->L * StdI->W * (1/*U*/ + 4 * (2/*V*/ + 2/*V'*/));

    if (model == "kondo")  StdI->nintr += 
      StdI->nsite / 2 * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  /**/
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double *)malloc(sizeof(double) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    StdI->transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }
  /**/
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++) {
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  /*
   Set Transfer & Interaction
  */
  StdI->ntrans = 0;
  StdI->nintr = 0;
  for (kCell = 0; kCell < StdI->NCell; kCell++){
    /**/
    iW = StdI->Cell[kCell][0];
    iL = StdI->Cell[kCell][1];
    /*
     Local term
    */
    isite = kCell;
    if (model == "kondo") isite += StdI->NCell;
    /**/
    if (model == "spin") {
      StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, isite);
      StdFace_GeneralJ(StdI, StdI->D, StdI->S2, StdI->S2, isite, isite);
    }/*if (model == "spin")*/
    else {
      StdFace_Hopping(StdI, StdI->mu, isite, isite, "local");
      StdFace_intr(StdI, StdI->U, isite, 0, isite, 0, isite, 1, isite, 1);
      /**/
      if (model == "kondo") {
        jsite = kCell;
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite, jsite);
      }/*if (model == "kondo")*/
    }
    /*
     Nearest neighbor along W
    */
    StdFace_SetLabel(StdI, fp, iW, iL, iW + 1, iL, 0.0, 0.0, 0.0, 0.0, isite, &jsite, 1, 1, 0, model);
    /**/
    if (model == "spin") {
      StdFace_GeneralJ(StdI, StdI->J0, StdI->S2, StdI->S2, isite, jsite);
    }/*if (model == "spin")*/
    else {
      StdFace_Hopping(StdI, StdI->t0, isite, jsite, "hopp");
      StdFace_Coulomb(StdI, StdI->V0, isite, jsite);
    }
    /*
     Nearest neighbor along L
    */
    StdFace_SetLabel(StdI, fp, iW, iL, iW, iL + 1, 0.0, 0.0, 0.0, 0.0, isite, &jsite, 1, 1, 0, model);
    /**/
    if (model == "spin") {
      StdFace_GeneralJ(StdI, StdI->J1, StdI->S2, StdI->S2, isite, jsite);
    }
    else {
      StdFace_Hopping(StdI, StdI->t1, isite, jsite, "hopp");
      StdFace_Coulomb(StdI, StdI->V1, isite, jsite);
    }
    /*
     Second nearest neighbor 1
    */
    StdFace_SetLabel(StdI, fp, iW, iL, iW + 1, iL + 1, 0.0, 0.0, 0.0, 0.0, isite, &jsite, 2, 1, 0, model);
    /**/
    if (model == "spin") {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (model == "spin")*/
    else {
      StdFace_Hopping(StdI, StdI->tp, isite, jsite, "hopp");
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
    /*
     Second nearest neighbor 2
    */
    StdFace_SetLabel(StdI, fp, iW, iL, iW + 1, iL - 1, 0.0, 0.0, 0.0, 0.0, isite, &jsite, 2, 1, 0, model);
    /**/
    if (model == "spin") {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (model == "spin")*/
    else {
      StdFace_Hopping(StdI, StdI->tp, isite, jsite, "hopp");
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }/*if (model != "spin")*/
  }/*for (kCell = 0; kCell < StdI->NCell; kCell++)*/

  fprintf(fp, "plot x w l lw 0\n");
  fprintf(fp, "pause -1\n");
  fclose(fp);
}
