/*
HPhi-mVMC-StdFace - Common input generator
Copyright (C) 2015 The University of Tokyo

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
#include <string.h>

/**
*
* Setup a Hamiltonian for the Kagome lattice
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_Kagome(struct StdIntList *StdI, char *model)
{
  int isite, jsite, kCell;
  int iL, iW;
  int ktrans, kintr;
  FILE *fp;
  double complex phase;

  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  /*
  Initialize Cell
  */
  fp = fopen("lattice.gp", "w");
  /**/
  StdI->NsiteUC = 3;
  /**/
  fprintf(stdout, "  @ Lattice Size & Shape\n\n");
  
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  StdFace_PrintVal_d("a0", &StdI->a0, StdI->a);
  StdFace_PrintVal_d("a1", &StdI->a1, StdI->a);
  StdFace_PrintVal_d("Wx", &StdI->Wx, StdI->a0);
  StdFace_PrintVal_d("Wy", &StdI->Wy, 0.0);
  StdFace_PrintVal_d("Lx", &StdI->Lx, StdI->a1 * 0.5);
  StdFace_PrintVal_d("Ly", &StdI->Ly, StdI->a1 * 0.5 * sqrt(3.0));
  
  StdFace_InitSite2D(StdI, fp);
  StdI->tau[0][0] = 0.0; StdI->tau[0][1] = 0.0;
  StdI->tau[1][0] = 0.5; StdI->tau[1][1] = 0.0;
  StdI->tau[2][0] = 0.0; StdI->tau[2][1] = 0.5;
  /**/
  StdFace_PrintVal_d("phase0", &StdI->phase0, 0.0);
  StdFace_PrintVal_d("phase1", &StdI->phase1, 0.0);
  StdI->ExpPhase0 = cos(StdI->pi180 * StdI->phase0) + I*sin(StdI->pi180 * StdI->phase0);
  StdI->ExpPhase1 = cos(StdI->pi180 * StdI->phase1) + I*sin(StdI->pi180 * StdI->phase1);
  if (cabs(StdI->ExpPhase0 + 1.0) < 0.000001) StdI->AntiPeriod0 = 1;
  else StdI->AntiPeriod0 = 0;
  if (cabs(StdI->ExpPhase1 + 1.0) < 0.000001) StdI->AntiPeriod1 = 1;
  else StdI->AntiPeriod1 = 0;
  /**/
  fprintf(stdout, "\n  @ Hamiltonian \n\n");
  /**/
  StdFace_NotUsed_J("J1'", StdI->J1pAll, StdI->J1p);
  StdFace_NotUsed_J("J2'", StdI->J2pAll, StdI->J2p);
  StdFace_NotUsed_d("t1'", StdI->t1p);
  StdFace_NotUsed_d("t2'", StdI->t2p);
  StdFace_NotUsed_d("V1'", StdI->V1p);
  StdFace_NotUsed_d("V2'", StdI->V2p);
  StdFace_NotUsed_d("K", StdI->K);
  /**/
  if (strcmp(StdI->model, "spin") == 0 ) {
    StdFace_PrintVal_i("2S", &StdI->S2, 1);
    StdFace_PrintVal_d("h", &StdI->h, 0.0);
    StdFace_PrintVal_d("Gamma", &StdI->Gamma, 0.0);
    StdFace_PrintVal_d("D", &StdI->D[2][2], 0.0);
    StdFace_InputSpinNN(StdI, StdI->J0, StdI->J0All, "J0");
    StdFace_InputSpinNN(StdI, StdI->J1, StdI->J1All, "J1");
    StdFace_InputSpinNN(StdI, StdI->J2, StdI->J2All, "J2");
    StdFace_InputSpin(StdI, StdI->Jp, StdI->JpAll, "J'");
    /**/
    StdFace_NotUsed_d("mu", StdI->mu);
    StdFace_NotUsed_d("U", StdI->U);
    StdFace_NotUsed_c("t", StdI->t);
    StdFace_NotUsed_c("t0", StdI->t0);
    StdFace_NotUsed_c("t1", StdI->t1);
    StdFace_NotUsed_c("t2", StdI->t2);
    StdFace_NotUsed_c("t'", StdI->tp);
    StdFace_NotUsed_d("V", StdI->V);
    StdFace_NotUsed_d("V0", StdI->V0);
    StdFace_NotUsed_d("V1", StdI->V1);
    StdFace_NotUsed_d("V2", StdI->V2);
    StdFace_NotUsed_d("V'", StdI->Vp);
  }/*if (strcmp(StdI->model, "spin") == 0 )*/
  else {
    StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
    StdFace_PrintVal_d("U", &StdI->U, 0.0);
    StdFace_InputHopp(StdI, &StdI->t0, "t0");
    StdFace_InputHopp(StdI, &StdI->t1, "t1");
    StdFace_InputHopp(StdI, &StdI->t2, "t2");
    StdFace_PrintVal_c("t'", &StdI->tp, 0.0);
    StdFace_InputCoulombV(StdI, &StdI->V0, "V0");
    StdFace_InputCoulombV(StdI, &StdI->V1, "V1");
    StdFace_InputCoulombV(StdI, &StdI->V2, "V2");
    StdFace_PrintVal_d("V'", &StdI->Vp, 0.0);
    /**/
    StdFace_NotUsed_J("J0", StdI->J0All, StdI->J0);
    StdFace_NotUsed_J("J1", StdI->J1All, StdI->J1);
    StdFace_NotUsed_J("J2", StdI->J2All, StdI->J2);
    StdFace_NotUsed_J("J'", StdI->JpAll, StdI->Jp);
    StdFace_NotUsed_d("h", StdI->h);
    StdFace_NotUsed_d("Gamma", StdI->Gamma);
    StdFace_NotUsed_d("D", StdI->D[2][2]);

    if (strcmp(StdI->model, "hubbard") == 0 ) {
      StdFace_NotUsed_i("2S", StdI->S2);
      StdFace_NotUsed_J("J", StdI->JAll, StdI->J);
    }/*if (strcmp(StdI->model, "hubbard") == 0 )*/
    else {
      StdFace_PrintVal_i("2S", &StdI->S2, 1);
      StdFace_InputSpin(StdI, StdI->J, StdI->JAll, "J");
    }/*if (model != "hubbard")*/

  }/*if (model != "spin")*/
  fprintf(stdout, "\n  @ Numerical conditions\n\n");
  /*
   Local Spin
  */
  StdI->nsite = 3 * StdI->NCell;
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
    StdI->ntrans = StdI->nsite * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
    StdI->nintr = StdI->NCell * (StdI->NsiteUC/*D*/ + 6/*J*/ + 6/*J'*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  else {
    StdI->ntrans = StdI->NCell * 2/*spin*/ * (StdI->NsiteUC/*mu*/ + 12/*t*/ + 12/*t'*/);
    StdI->nintr = StdI->NCell * (StdI->NsiteUC/*U*/ + 4 * (6/*V*/ + 6/*V'*/));

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
  StdFace_MallocInteractions(StdI);/*debug*/
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
    if (strcmp(StdI->model, "kondo") == 0 ) isite += StdI->nsite / 2;
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, isite);
      StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, isite + 1);
      StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, isite + 2);
      StdFace_GeneralJ(StdI, StdI->D, StdI->S2, StdI->S2, isite, isite);
      StdFace_GeneralJ(StdI, StdI->D, StdI->S2, StdI->S2, isite + 1, isite + 1);
      StdFace_GeneralJ(StdI, StdI->D, StdI->S2, StdI->S2, isite + 2, isite + 2);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, StdI->mu, isite, isite, 0);
      StdFace_Hopping(StdI, StdI->mu, isite + 1, isite + 1, 0);
      StdFace_Hopping(StdI, StdI->mu, isite + 2, isite + 2, 0);
      StdI->Cintra[StdI->NCintra] = StdI->U; StdI->CintraIndx[StdI->NCintra][0] = isite; StdI->NCintra += 1;
      StdI->Cintra[StdI->NCintra] = StdI->U; StdI->CintraIndx[StdI->NCintra][0] = isite + 1; StdI->NCintra += 1;
      StdI->Cintra[StdI->NCintra] = StdI->U; StdI->CintraIndx[StdI->NCintra][0] = isite + 2; StdI->NCintra += 1;
      /**/
      if (strcmp(StdI->model, "kondo") == 0 ) {
        jsite = StdI->NsiteUC * kCell;
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite, jsite);
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite + 1, jsite + 1);
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite + 2, jsite + 2);
      }/*if (strcmp(StdI->model, "kondo") == 0 )*/
    }/*if (strcmp(StdI->model, "spin") != 0 )*/
    /*
    Nearest neighbor intra cell 0 -> 1
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 0, 0, 0, 1, &isite, &jsite, 1, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J2, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->t2, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V2, isite, jsite);
    }
    /*
    Nearest neighbor intra cell 0 -> 2
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 0, 0, 0, 2, &isite, &jsite, 1, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J1, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->t1, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V1, isite, jsite);
    }
    /*
    Nearest neighbor intra cell 1 -> 2
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 0, 0, 1, 2, &isite, &jsite, 1, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J0, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->t0, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V0, isite, jsite);
    }
    /*
    Nearest neighbor along W 1 -> 0
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 1, 0, 1, 0, &isite, &jsite, 1, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J2, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->t2, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V2, isite, jsite);
    }
    /*
    Nearest neighbor along L 2 -> 0
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 0, 1, 2, 0, &isite, &jsite, 1, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J1, StdI->S2, StdI->S2, isite, jsite);
    }
    else {
      StdFace_Hopping(StdI, phase * StdI->t1, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V1, isite, jsite);
    }
    /*
    Nearest neighbor along W-L 1 -> 2
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 1, - 1, 1, 2, &isite, &jsite, 1, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J0, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->t0, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V0, isite, jsite);
    }
    /*
    Second nearest neighbor along W 2 -> 0
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 1, 0, 2, 0, &isite, &jsite, 2, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->tp, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
    /*
    Second nearest neighbor along W 1 -> 2
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 1, 0, 1, 2, &isite, &jsite, 2, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->tp, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
    /*
    Second nearest neighbor along L 1 -> 0
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 0, 1, 1, 0, &isite, &jsite, 2, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->tp, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
    /*
    Second nearest neighbor along L 2 -> 1
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 0, 1, 2, 1, &isite, &jsite, 2, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->tp, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
    /*
    Second nearest neighbor along W-L 0 -> 2
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 1, - 1, 0, 2, &isite, &jsite, 2, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->tp, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
    /*
    Second nearest neighbor along L-W 0 -> 1
    */
    StdFace_SetLabel(StdI, fp, iW, iL, - 1, 1, 0, 1, &isite, &jsite, 2, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->tp, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
  }/*for (kCell = 0; kCell < StdI->NCell; kCell++)*/

  fprintf(fp, "plot \'-\' w d lc 7\n0.0 0.0\nend\npause -1\n");
  fclose(fp);
  StdFace_PrintGeometry(StdI);
}/*void StdFace_Kagome*/

#if defined(_HPhi)
/**
*
* Setup a Hamiltonian for the generalized Heisenberg model on a Heisenberg lattice
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_Kagome_Boost(struct StdIntList *StdI)
{
  int isite, ipivot, i1, i2;
  int kintr;
  FILE *fp;

  if (StdI->a0L != 0 || StdI->a1W != 0) {
    fprintf(stdout, "\nERROR ! (a0W, a0L, a1W, a1L) can not be used with SpinGCBoost.\n\n");
    StdFace_exit(-1);
  }
  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 3; i2++) {
      if (fabs(StdI->Jp[i1][i2]) > 1.0e-8) {
        fprintf(stdout, "\nERROR ! J' can not be used with SpinGCBoost.\n\n");
        StdFace_exit(-1);
      }
    }
  }
  /*
  Magnetic field
  */
  fp = fopen("boost.def", "w");
  fprintf(fp, "# Magnetic field\n");
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    -0.5 * StdI->Gamma, 0.0, -0.5 * StdI->h);
  /*
  Interaction
  */
  fprintf(fp, "%d  # Number of type of J\n", 3);
  fprintf(fp, "# J 0\n");
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J0[0][0], 0.25 * StdI->J0[0][1], 0.25 * StdI->J0[0][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J0[0][1], 0.25 * StdI->J0[1][1], 0.25 * StdI->J0[1][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J0[0][2], 0.25 * StdI->J0[1][2], 0.25 * StdI->J0[2][2]);
  fprintf(fp, "# J 1\n");
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1[0][0], 0.25 * StdI->J1[0][1], 0.25 * StdI->J1[0][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1[0][1], 0.25 * StdI->J1[1][1], 0.25 * StdI->J1[1][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1[0][2], 0.25 * StdI->J1[1][2], 0.25 * StdI->J1[2][2]);
  fprintf(fp, "# J 2\n");
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2[0][0], 0.25 * StdI->J2[0][1], 0.25 * StdI->J2[0][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2[0][1], 0.25 * StdI->J2[1][1], 0.25 * StdI->J2[1][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2[0][2], 0.25 * StdI->J2[1][2], 0.25 * StdI->J2[2][2]);  /*
  Topology
  */
  if (StdI->S2 != 1) {
    fprintf(stdout, "\n ERROR! S2 must be 1 in Boost. \n\n");
    StdFace_exit(-1);
  }
  StdI->ishift_nspin = 3;
  if (StdI->L < 2) {
    fprintf(stdout, "\n ERROR! L < 2 \n\n");
    StdFace_exit(-1);
  }
  if (StdI->W % StdI->ishift_nspin != 0) {
    fprintf(stdout, "\n ERROR! W %% %d != 0 \n\n", StdI->ishift_nspin);
    StdFace_exit(-1);
  }
  StdI->num_pivot = 4;
  if (StdI->W != 3) {
    fprintf(stdout, "DEBUG: W != 3\n");
    StdFace_exit(-1);
  }
  StdI->W = 9;
  fprintf(fp, "# W0  R0  StdI->num_pivot  StdI->ishift_nspin\n");
  fprintf(fp, "%d %d %d %d\n", StdI->W, StdI->L, StdI->num_pivot, StdI->ishift_nspin);

  StdI->list_6spin_star = (int **)malloc(sizeof(int*) * StdI->num_pivot);
  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    StdI->list_6spin_star[ipivot] = (int *)malloc(sizeof(int) * 7);
  }

  StdI->list_6spin_star[0][0] = 1; // num of J
  StdI->list_6spin_star[0][1] = 1;
  StdI->list_6spin_star[0][2] = 1;
  StdI->list_6spin_star[0][3] = 1;
  StdI->list_6spin_star[0][4] = 4;
  StdI->list_6spin_star[0][5] = 2;
  StdI->list_6spin_star[0][6] = -1; // flag

  StdI->list_6spin_star[1][0] = 6; // num of J
  StdI->list_6spin_star[1][1] = 1;
  StdI->list_6spin_star[1][2] = 1;
  StdI->list_6spin_star[1][3] = 1;
  StdI->list_6spin_star[1][4] = 6;
  StdI->list_6spin_star[1][5] = 7;
  StdI->list_6spin_star[1][6] = 1; // flag

  StdI->list_6spin_star[2][0] = 6; // num of J
  StdI->list_6spin_star[2][1] = 1;
  StdI->list_6spin_star[2][2] = 1;
  StdI->list_6spin_star[2][3] = 1;
  StdI->list_6spin_star[2][4] = 4;
  StdI->list_6spin_star[2][5] = 2;
  StdI->list_6spin_star[2][6] = 1; // flag

  StdI->list_6spin_star[3][0] = 5; // num of J
  StdI->list_6spin_star[3][1] = 1;
  StdI->list_6spin_star[3][2] = 1;
  StdI->list_6spin_star[3][3] = 1;
  StdI->list_6spin_star[3][4] = 4;
  StdI->list_6spin_star[3][5] = 2;
  StdI->list_6spin_star[3][6] = 1; // flag

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
  StdI->list_6spin_pair[0][1][0] = 4; //(2,1,1+2*j)=1
  StdI->list_6spin_pair[0][2][0] = 1; //(3,1,1+2*j)=2
  StdI->list_6spin_pair[0][3][0] = 2; //(4,1,1+2*j)=3
  StdI->list_6spin_pair[0][4][0] = 3; //(5,1,1+2*j)=4
  StdI->list_6spin_pair[0][5][0] = 5; //(6,1,1+2*j)=5
  StdI->list_6spin_pair[0][6][0] = 3; //(7,1,1+2*j)=3 ! type of J

  StdI->list_6spin_pair[1][0][0] = 0;
  StdI->list_6spin_pair[1][1][0] = 1;
  StdI->list_6spin_pair[1][2][0] = 2;
  StdI->list_6spin_pair[1][3][0] = 3;
  StdI->list_6spin_pair[1][4][0] = 4;
  StdI->list_6spin_pair[1][5][0] = 5;
  StdI->list_6spin_pair[1][6][0] = 3; // type of J
  StdI->list_6spin_pair[1][0][1] = 1;
  StdI->list_6spin_pair[1][1][1] = 2;
  StdI->list_6spin_pair[1][2][1] = 0;
  StdI->list_6spin_pair[1][3][1] = 3;
  StdI->list_6spin_pair[1][4][1] = 4;
  StdI->list_6spin_pair[1][5][1] = 5;
  StdI->list_6spin_pair[1][6][1] = 1; // type of J
  StdI->list_6spin_pair[1][0][2] = 0;
  StdI->list_6spin_pair[1][1][2] = 2;
  StdI->list_6spin_pair[1][2][2] = 1;
  StdI->list_6spin_pair[1][3][2] = 3;
  StdI->list_6spin_pair[1][4][2] = 4;
  StdI->list_6spin_pair[1][5][2] = 5;
  StdI->list_6spin_pair[1][6][2] = 2; // type of J
  StdI->list_6spin_pair[1][0][3] = 1;
  StdI->list_6spin_pair[1][1][3] = 3;
  StdI->list_6spin_pair[1][2][3] = 0;
  StdI->list_6spin_pair[1][3][3] = 2;
  StdI->list_6spin_pair[1][4][3] = 4;
  StdI->list_6spin_pair[1][5][3] = 5;
  StdI->list_6spin_pair[1][6][3] = 3; // type of J
  StdI->list_6spin_pair[1][0][4] = 2;
  StdI->list_6spin_pair[1][1][4] = 4;
  StdI->list_6spin_pair[1][2][4] = 0;
  StdI->list_6spin_pair[1][3][4] = 1;
  StdI->list_6spin_pair[1][4][4] = 3;
  StdI->list_6spin_pair[1][5][4] = 5;
  StdI->list_6spin_pair[1][6][4] = 2; // type of J
  StdI->list_6spin_pair[1][0][5] = 2;
  StdI->list_6spin_pair[1][1][5] = 5;
  StdI->list_6spin_pair[1][2][5] = 0;
  StdI->list_6spin_pair[1][3][5] = 1;
  StdI->list_6spin_pair[1][4][5] = 3;
  StdI->list_6spin_pair[1][5][5] = 4;
  StdI->list_6spin_pair[1][6][5] = 1; // type of J

  StdI->list_6spin_pair[2][0][0] = 0;
  StdI->list_6spin_pair[2][1][0] = 1;
  StdI->list_6spin_pair[2][2][0] = 2;
  StdI->list_6spin_pair[2][3][0] = 3;
  StdI->list_6spin_pair[2][4][0] = 4;
  StdI->list_6spin_pair[2][5][0] = 5;
  StdI->list_6spin_pair[2][6][0] = 3; // type of J
  StdI->list_6spin_pair[2][0][1] = 1;
  StdI->list_6spin_pair[2][1][1] = 2;
  StdI->list_6spin_pair[2][2][1] = 0;
  StdI->list_6spin_pair[2][3][1] = 3;
  StdI->list_6spin_pair[2][4][1] = 4;
  StdI->list_6spin_pair[2][5][1] = 5;
  StdI->list_6spin_pair[2][6][1] = 1; // type of J
  StdI->list_6spin_pair[2][0][2] = 0;
  StdI->list_6spin_pair[2][1][2] = 2;
  StdI->list_6spin_pair[2][2][2] = 1;
  StdI->list_6spin_pair[2][3][2] = 3;
  StdI->list_6spin_pair[2][4][2] = 4;
  StdI->list_6spin_pair[2][5][2] = 5;
  StdI->list_6spin_pair[2][6][2] = 2; // type of J
  StdI->list_6spin_pair[2][0][3] = 1;
  StdI->list_6spin_pair[2][1][3] = 3;
  StdI->list_6spin_pair[2][2][3] = 0;
  StdI->list_6spin_pair[2][3][3] = 2;
  StdI->list_6spin_pair[2][4][3] = 4;
  StdI->list_6spin_pair[2][5][3] = 5;
  StdI->list_6spin_pair[2][6][3] = 3; // type of J
  StdI->list_6spin_pair[2][0][4] = 2;
  StdI->list_6spin_pair[2][1][4] = 5;
  StdI->list_6spin_pair[2][2][4] = 0;
  StdI->list_6spin_pair[2][3][4] = 1;
  StdI->list_6spin_pair[2][4][4] = 3;
  StdI->list_6spin_pair[2][5][4] = 4;
  StdI->list_6spin_pair[2][6][4] = 2; // type of J
  StdI->list_6spin_pair[2][0][5] = 2;
  StdI->list_6spin_pair[2][1][5] = 4;
  StdI->list_6spin_pair[2][2][5] = 0;
  StdI->list_6spin_pair[2][3][5] = 1;
  StdI->list_6spin_pair[2][4][5] = 3;
  StdI->list_6spin_pair[2][5][5] = 5;
  StdI->list_6spin_pair[2][6][5] = 1; // type of J

  StdI->list_6spin_pair[3][0][0] = 0;
  StdI->list_6spin_pair[3][1][0] = 1;
  StdI->list_6spin_pair[3][2][0] = 2;
  StdI->list_6spin_pair[3][3][0] = 3;
  StdI->list_6spin_pair[3][4][0] = 4;
  StdI->list_6spin_pair[3][5][0] = 5;
  StdI->list_6spin_pair[3][6][0] = 3; // type of J
  StdI->list_6spin_pair[3][0][1] = 1;
  StdI->list_6spin_pair[3][1][1] = 2;
  StdI->list_6spin_pair[3][2][1] = 0;
  StdI->list_6spin_pair[3][3][1] = 3;
  StdI->list_6spin_pair[3][4][1] = 4;
  StdI->list_6spin_pair[3][5][1] = 5;
  StdI->list_6spin_pair[3][6][1] = 1; // type of J
  StdI->list_6spin_pair[3][0][2] = 0;
  StdI->list_6spin_pair[3][1][2] = 2;
  StdI->list_6spin_pair[3][2][2] = 1;
  StdI->list_6spin_pair[3][3][2] = 3;
  StdI->list_6spin_pair[3][4][2] = 4;
  StdI->list_6spin_pair[3][5][2] = 5;
  StdI->list_6spin_pair[3][6][2] = 2; // type of J
  StdI->list_6spin_pair[3][0][3] = 2;
  StdI->list_6spin_pair[3][1][3] = 5;
  StdI->list_6spin_pair[3][2][3] = 0;
  StdI->list_6spin_pair[3][3][3] = 1;
  StdI->list_6spin_pair[3][4][3] = 3;
  StdI->list_6spin_pair[3][5][3] = 4;
  StdI->list_6spin_pair[3][6][3] = 2; // type of J
  StdI->list_6spin_pair[3][0][4] = 2;
  StdI->list_6spin_pair[3][1][4] = 4;
  StdI->list_6spin_pair[3][2][4] = 0;
  StdI->list_6spin_pair[3][3][4] = 1;
  StdI->list_6spin_pair[3][4][4] = 3;
  StdI->list_6spin_pair[3][5][4] = 5;
  StdI->list_6spin_pair[3][6][4] = 1; // type of J

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
#endif
