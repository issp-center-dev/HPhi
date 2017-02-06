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
 * Setup a Hamiltonian for the Pyrochlore structure
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_Pyrochlore(struct StdIntList *StdI, char *model)
{
  int isite, jsite;
  int iL, iW, iH, kCell;
  FILE *fp;
  double complex Cphase;

  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  /*
   Initialize Cell
  */
  StdI->NsiteUC = 4;
  /**/
  fprintf(stdout, "  @ Lattice Size & Shape\n\n");
  
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  StdFace_PrintVal_d("Wlength", &StdI->length[0], StdI->a);
  StdFace_PrintVal_d("Llength", &StdI->length[1], StdI->a);
  StdFace_PrintVal_d("Hlength", &StdI->length[2], StdI->a);
  StdFace_PrintVal_d("Wx", &StdI->direct[0][0], 0.0);
  StdFace_PrintVal_d("Wy", &StdI->direct[0][1], 0.5 * StdI->length[0]);
  StdFace_PrintVal_d("Wz", &StdI->direct[0][2], 0.5 * StdI->length[0]);
  StdFace_PrintVal_d("Lx", &StdI->direct[1][0], 0.5 * StdI->length[1]);
  StdFace_PrintVal_d("Ly", &StdI->direct[1][1], 0.0);
  StdFace_PrintVal_d("Lz", &StdI->direct[1][2], 0.5 * StdI->length[1]);
  StdFace_PrintVal_d("Hx", &StdI->direct[2][0], 0.5 * StdI->length[2]);
  StdFace_PrintVal_d("Hy", &StdI->direct[2][1], 0.5 * StdI->length[2]);
  StdFace_PrintVal_d("Hz", &StdI->direct[2][2], 0.0);

  StdFace_InitSite3D(StdI, fp);
  StdI->tau[0][0] = 0.0; StdI->tau[0][1] = 0.0; ; StdI->tau[0][2] = 0.0;
  StdI->tau[1][0] = 0.5; StdI->tau[1][1] = 0.0; ; StdI->tau[1][2] = 0.0;
  StdI->tau[2][0] = 0.0; StdI->tau[2][1] = 0.5; ; StdI->tau[2][2] = 0.0;
  StdI->tau[3][0] = 0.0; StdI->tau[3][1] = 0.0; ; StdI->tau[3][2] = 0.5;
  /**/
  StdFace_PrintVal_d("phase0", &StdI->phase[0], 0.0);
  StdFace_PrintVal_d("phase1", &StdI->phase[1], 0.0);
  StdFace_PrintVal_d("phase2", &StdI->phase[2], 0.0);
  /**/
  fprintf(stdout, "\n  @ Hamiltonian \n\n");
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
    StdFace_InputSpinNN(StdI, StdI->J0p, StdI->J0pAll, "J0'");
    StdFace_InputSpinNN(StdI, StdI->J1p, StdI->J1pAll, "J1'");
    StdFace_InputSpinNN(StdI, StdI->J2p, StdI->J2pAll, "J2'");
    StdFace_InputSpin(StdI, StdI->Jpp, StdI->JppAll, "J''");
    /**/
    StdFace_NotUsed_d("mu", StdI->mu);
    StdFace_NotUsed_d("U", StdI->U);
    StdFace_NotUsed_c("t", StdI->t);
    StdFace_NotUsed_c("t0", StdI->t0);
    StdFace_NotUsed_c("t1", StdI->t1);
    StdFace_NotUsed_c("t2", StdI->t2);
    StdFace_NotUsed_c("t'", StdI->tp);
    StdFace_NotUsed_c("t0'", StdI->t0p);
    StdFace_NotUsed_c("t1'", StdI->t1p);
    StdFace_NotUsed_c("t2'", StdI->t2p);
    StdFace_NotUsed_c("t''", StdI->tpp);
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
    StdFace_InputHopp(StdI, &StdI->t2, "t2");
    StdFace_InputHopp(StdI, &StdI->t0p, "t0'");
    StdFace_InputHopp(StdI, &StdI->t1p, "t1'");
    StdFace_InputHopp(StdI, &StdI->t2p, "t2'");
    StdFace_InputCoulombV(StdI, &StdI->V0, "V0");
    StdFace_InputCoulombV(StdI, &StdI->V1, "V1");
    StdFace_InputCoulombV(StdI, &StdI->V2, "V2");
    StdFace_InputCoulombV(StdI, &StdI->V0p, "V0'");
    StdFace_InputCoulombV(StdI, &StdI->V1p, "V1'");
    StdFace_InputCoulombV(StdI, &StdI->V2p, "V2'");
    /**/
    StdFace_NotUsed_J("J0", StdI->J0All, StdI->J0);
    StdFace_NotUsed_J("J1", StdI->J1All, StdI->J1);
    StdFace_NotUsed_J("J2", StdI->J2All, StdI->J2);
    StdFace_NotUsed_J("J0'", StdI->J0pAll, StdI->J0p);
    StdFace_NotUsed_J("J1'", StdI->J1pAll, StdI->J1p);
    StdFace_NotUsed_J("J2'", StdI->J2pAll, StdI->J2p);
    StdFace_NotUsed_J("J''", StdI->JppAll, StdI->Jpp);
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
    StdI->nintr = StdI->NCell * (StdI->NsiteUC/*D*/ + 12/*J*/ + 0/*J'*/ + 0/*J''*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  else {
    StdI->ntrans = StdI->NCell * 2/*spin*/ * (StdI->NsiteUC/*mu*/ + 24/*t*/ + 0/*t'*/ + 0/*t''*/);
    StdI->nintr = StdI->NCell * (StdI->NsiteUC/*U*/ + 4 * (12/*V*/ + 0/*V'*/ + 0/*V''*/));

    if (strcmp(StdI->model, "kondo") == 0 )  StdI->nintr += 
      StdI->nsite / 2 * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  /**/
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
    iH = StdI->Cell[kCell][2];
    /*
     (1) Local term
    */
    isite = StdI->NsiteUC * kCell;
    if (strcmp(StdI->model, "kondo") == 0) isite += StdI->nsite / 2;
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
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
      StdI->Cintra[StdI->NCintra] = StdI->U; StdI->CintraIndx[StdI->NCintra][0] = isite + 3; StdI->NCintra += 1;
      /**/
      if (strcmp(StdI->model, "kondo") == 0) {
        jsite = StdI->NsiteUC * kCell;
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite, jsite);
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite + 1, jsite + 1);
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite + 2, jsite + 2);
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite + 3, jsite + 3);
      }/*if (strcmp(StdI->model, "kondo") == 0 )*/
    }/*if (strcmp(StdI->model, "spin") != 0 )*/
    /*
     (2) Intra-Cell along W
    */
    StdFace_FindSite3d(StdI, iW, iL, iH, 0, 0, 0, 0, 1, &isite, &jsite, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J0, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t0, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V0, isite, jsite);
    }
    /*
     (3) Intra-Cell along L
    */
    StdFace_FindSite3d(StdI, iW, iL, iH, 0, 0, 0, 0, 2, &isite, &jsite, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J1, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t1, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V1, isite, jsite);
    }
    /*
     (4) Intra-Cell along H
    */
    StdFace_FindSite3d(StdI, iW, iL, iH, 0, 0, 0, 0, 3, &isite, &jsite, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J2, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t2, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V2, isite, jsite);
    }
    /*
     (5) Intra-Cell along L-H
    */
    StdFace_FindSite3d(StdI, iW, iL, iH, 0, 0, 0, 2, 3, &isite, &jsite, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J0p, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t0p, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V0p, isite, jsite);
    }
    /*
     (6) Intra-Cell along H-W
    */
    StdFace_FindSite3d(StdI, iW, iL, iH, 0, 0, 0, 3, 1, &isite, &jsite, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J1p, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t1p, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V1p, isite, jsite);
    }
    /*
     (7) Intra-Cell along W-L
    */
    StdFace_FindSite3d(StdI, iW, iL, iH, 0, 0, 0, 1, 2, &isite, &jsite, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J2p, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t2p, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V2p, isite, jsite);
    }
    /*
     (8) Inter-Cell along W
    */
    StdFace_FindSite3d(StdI, iW, iL, iH, 1, 0, 0, 1, 0, &isite, &jsite, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J0, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t0, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V0, isite, jsite);
    }
    /*
     (9) Inter-Cell along L
    */
    StdFace_FindSite3d(StdI, iW, iL, iH, 0, 1, 0, 2, 0, &isite, &jsite, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J1, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t1, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V1, isite, jsite);
    }
    /*
     (10) Inter-Cell along H
    */
    StdFace_FindSite3d(StdI, iW, iL, iH, 0, 0, 1, 3, 0, &isite, &jsite, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J2, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t2, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V2, isite, jsite);
    }
    /*
    (11) Inter-Cell along L-H
    */
    StdFace_FindSite3d(StdI, iW, iL, iH, 0, -1, 1, 3, 2, &isite, &jsite, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J0p, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t0p, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V0p, isite, jsite);
    }
    /*
    (12) Intra-Cell along H-W
    */
    StdFace_FindSite3d(StdI, iW, iL, iH, 1, 0, -1, 1, 3, &isite, &jsite, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J1p, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t1p, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V1p, isite, jsite);
    }
    /*
    (13) Intra-Cell along W-L
    */
    StdFace_FindSite3d(StdI, iW, iL, iH, -1, 1, 0, 2, 1, &isite, &jsite, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J2p, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t2p, isite, jsite, 1);
      StdFace_Coulomb(StdI, StdI->V2p, isite, jsite);
    }
  }/*for (kCell = 0; kCell < StdI->NCell; kCell++)*/

  StdFace_PrintXSF(StdI);
  StdFace_PrintGeometry(StdI);
}/*void StdFace_Pyrochlore*/
