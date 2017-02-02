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
 * Setup a Hamiltonian for the square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_Tetragonal(struct StdIntList *StdI, char *model)
{
  int isite, jsite;
  int iL, iW, kCell;
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
  StdI->NsiteUC = 1;
  /**/
  fprintf(stdout, "  @ Lattice Size & Shape\n\n");
  
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  StdFace_PrintVal_d("length[0]", &StdI->length[0], StdI->a);
  StdFace_PrintVal_d("length[1]", &StdI->length[1], StdI->a);
  StdFace_PrintVal_d("direct[0][0]", &StdI->direct[0][0], StdI->length[0]);
  StdFace_PrintVal_d("direct[0][1]", &StdI->direct[0][1], 0.0);
  StdFace_PrintVal_d("direct[1][0]", &StdI->direct[1][0], 0.0);
  StdFace_PrintVal_d("direct[1][1]", &StdI->direct[1][1], StdI->length[1]);
  
  StdFace_InitSite2D(StdI, fp);
  StdI->tau[0][0] = 0.0; StdI->tau[0][1] = 0.0;
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
  StdFace_NotUsed_J("J2", StdI->J2All, StdI->J2);
  StdFace_NotUsed_J("J1'", StdI->J1pAll, StdI->J1p);
  StdFace_NotUsed_J("J2'", StdI->J2pAll, StdI->J2p);
  StdFace_NotUsed_c("t2", StdI->t2);
  StdFace_NotUsed_d("t1'", StdI->t1p);
  StdFace_NotUsed_d("t2'", StdI->t2p);
  StdFace_NotUsed_d("V2", StdI->V2);
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
    StdFace_InputSpin(StdI, StdI->Jp, StdI->JpAll, "J'");
    /**/
    StdFace_NotUsed_d("mu", StdI->mu);
    StdFace_NotUsed_d("U", StdI->U);
    StdFace_NotUsed_c("t", StdI->t);
    StdFace_NotUsed_c("t0", StdI->t0);
    StdFace_NotUsed_c("t1", StdI->t1);
    StdFace_NotUsed_c("t'", StdI->tp);
    StdFace_NotUsed_d("V", StdI->V);
    StdFace_NotUsed_d("V0", StdI->V0);
    StdFace_NotUsed_d("V1", StdI->V1);
    StdFace_NotUsed_d("V'", StdI->Vp);
  }/*if (strcmp(StdI->model, "spin") == 0 )*/
  else {
    StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
    StdFace_PrintVal_d("U", &StdI->U, 0.0);
    StdFace_InputHopp(StdI, &StdI->t0, "t0");
    StdFace_InputHopp(StdI, &StdI->t1, "t1");
    StdFace_PrintVal_c("t'", &StdI->tp, 0.0);
    StdFace_InputCoulombV(StdI, &StdI->V0, "V0");
    StdFace_InputCoulombV(StdI, &StdI->V1, "V1");
    StdFace_PrintVal_d("V'", &StdI->Vp, 0.0);
    /**/
    StdFace_NotUsed_J("J0", StdI->J0All, StdI->J0);
    StdFace_NotUsed_J("J1", StdI->J1All, StdI->J1);
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
  StdI->nsite = StdI->NsiteUC * StdI->NCell;
  if (strcmp(StdI->model, "kondo") == 0 ) StdI->nsite *= 2;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  /**/
  if(strcmp(StdI->model, "spin") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++) StdI->locspinflag[isite] = StdI->S2;
  else if(strcmp(StdI->model, "hubbard") == 0 )
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
    StdI->nintr = StdI->NCell * (StdI->NsiteUC/*D*/ + 2/*J*/ + 2/*J'*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  else {
    StdI->ntrans = StdI->NCell * 2/*spin*/ * (StdI->NsiteUC/*mu*/ + 4/*t*/ + 4/*t'*/);
    StdI->nintr = StdI->NCell * (StdI->NsiteUC/*U*/ + 4 * (2/*V*/ + 2/*V'*/));

    if (strcmp(StdI->model, "kondo") == 0 )  StdI->nintr += 
      StdI->nsite / 2 * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
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
  StdFace_MallocInteractions(StdI);
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
    if (strcmp(StdI->model, "kondo") == 0 ) isite += StdI->NCell;
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, isite);
      StdFace_GeneralJ(StdI, StdI->D, StdI->S2, StdI->S2, isite, isite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, StdI->mu, isite, isite, 0);
      StdI->Cintra[StdI->NCintra] = StdI->U; StdI->CintraIndx[StdI->NCintra][0] = isite; StdI->NCintra += 1;
      /**/
      if (strcmp(StdI->model, "kondo") == 0 ) {
        jsite = kCell;
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite, jsite);
      }/*if (strcmp(StdI->model, "kondo") == 0 )*/
    }
    /*
     Nearest neighbor along W
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 1, 0, 0, 0, &isite, &jsite, 1, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J0, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->t0, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V0, isite, jsite);
    }
    /*
     Nearest neighbor along L
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 0, 1, 0, 0, &isite, &jsite, 1, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J1, StdI->S2, StdI->S2, isite, jsite);
    }
    else {
      StdFace_Hopping(StdI, phase * StdI->t1, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V1, isite, jsite);
    }
    /*
     Second nearest neighbor 1
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 1, 1, 0, 0, &isite, &jsite, 2, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->tp, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
    /*
     Second nearest neighbor 2
    */
    StdFace_SetLabel(StdI, fp, iW, iL, 1, -1, 0, 0, &isite, &jsite, 2, &phase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, phase * StdI->tp, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }/*if (model != "spin")*/
  }/*for (kCell = 0; kCell < StdI->NCell; kCell++)*/

  fprintf(fp, "plot \'-\' w d lc 7\n0.0 0.0\nend\npause -1\n");
  fclose(fp);
  StdFace_PrintGeometry(StdI);
}
