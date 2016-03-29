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
#include <complex.h>
#include "../include/wrapperMPI.h"
#include <string.h>

/**
 *
 * Setup a Hamiltonian for the Hubbard model on a Chain lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_Chain(struct StdIntList *StdI, char *model)
{
  int isite, jsite;
  int ispin, iL;
  int ktrans, kintr;
  /**/
  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  /**/
  StdI->NsiteUC = 1;
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_NotUsed_i("W", StdI->W);
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  /**/
  StdFace_NotUsed_c("t0", StdI->t0);
  StdFace_NotUsed_c("t1", StdI->t1);
  StdFace_NotUsed_c("t2", StdI->t2);
  StdFace_NotUsed_d("V0", StdI->V0);
  StdFace_NotUsed_d("V1", StdI->V1);
  StdFace_NotUsed_d("V2", StdI->V2);
  StdFace_NotUsed_d("K", StdI->K);
  /**/
  if (strcmp(StdI->model, "spin") == 0 ) {
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
    StdFace_NotUsed_c("t", StdI->t);
    StdFace_NotUsed_c("t'", StdI->tp);
    StdFace_NotUsed_d("mu", StdI->mu);
  }
  else {

    StdFace_PrintVal_c("t", &StdI->t, 1.0);
    StdFace_PrintVal_c("t'", &StdI->tp, 0.0);
    StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
    StdFace_PrintVal_d("U", &StdI->U, 0.0);
    StdFace_PrintVal_d("V", &StdI->V, 0.0);
    StdFace_PrintVal_d("V'", &StdI->Vp, 0.0);

    if (strcmp(StdI->model, "hubbard") == 0 ) {
      StdFace_NotUsed_i("2S", StdI->S2);
      StdFace_NotUsed_d("J", StdI->JAll);
    }
    else if (strcmp(StdI->model, "kondo") == 0 ) {
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
  }
  /*
  Local Spin
  */
  StdI->nsite = StdI->L;
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
    StdI->ntrans = StdI->L * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
    StdI->nintr = StdI->L * (1/*D*/ + 1/*J*/ + 1/*J'*/) 
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  else {
    StdI->ntrans = StdI->L * 2/*spin*/ * (1/*mu*/ + 2/*t*/ + 2/*t'*/);
    StdI->nintr = StdI->L * (1/*U*/ + 4 * (1/*V*/ + 1/*V'*/));

    if(strcmp(StdI->model, "kondo") == 0 ) 
      StdI->nintr += StdI->nsite / 2 * (3 * 1 + 1) * (3 * StdI->S2 + 1);
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
  for (iL = 0; iL < StdI->L; iL++){

    isite = iL;
    if (strcmp(StdI->model, "kondo") == 0 ) isite += StdI->L;
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
        jsite = iL;
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite, jsite);
      }/*if (strcmp(StdI->model, "kondo") == 0 )*/
    }/*if (model != "spin")*/
    /*
    Nearest neighbor
   */
    jsite = (iL + 1) % StdI->L;
    if (strcmp(StdI->model, "kondo") == 0 ) jsite += StdI->L;
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J, StdI->S2, StdI->S2, isite, jsite);
    }
    else {
      StdFace_Hopping(StdI, StdI->t, isite, jsite);
      StdFace_Coulomb(StdI, StdI->V, isite, jsite);
    }
    /*
    Second nearest neighbor
    */
    jsite = (iL + 2) % StdI->L;
    if (strcmp(StdI->model, "kondo") == 0 ) jsite += StdI->L;
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }
    else {
      StdFace_Hopping(StdI, StdI->tp, isite, jsite);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
  }/*for (iL = 0; iL < StdI->L; iL++)*/
}

/**
*
* Setup a Hamiltonian for the generalized Heisenberg model on a Chain lattice
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_Chain_Boost(struct StdIntList *StdI)
{
  int isite, ipivot;
  int ktrans, kintr;
  double LargeValue0, S;
  FILE *fp;

  StdI->NsiteUC = 1;
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
    0.25 * StdI->J[0][0], 0.0, 0.25 * StdI->J[0][1], 0.0, 0.25 * StdI->J[0][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J[1][0], 0.0, 0.25 * StdI->J[1][1], 0.0, 0.25 * StdI->J[1][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->J[2][0], 0.0, 0.25 * StdI->J[2][1], 0.0, 0.25 * StdI->J[2][2], 0.0);
  fprintf(fp, "# J 2\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->Jp[0][0], 0.0, 0.25 * StdI->Jp[0][1], 0.0, 0.25 * StdI->Jp[0][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->Jp[1][0], 0.0, 0.25 * StdI->J[1][1], 0.0, 0.25 * StdI->Jp[1][2], 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->Jp[2][0], 0.0, 0.25 * StdI->Jp[2][1], 0.0, 0.25 * StdI->J[2][2], 0.0);
  /*
  Topology
  */
  if (StdI->S2 != 1) {
    fprintf(stderr, "\n ERROR! S2 must be 1 in Boost. \n\n");
    exitMPI(-1);
  }
  StdI->ishift_nspin = 4;
  if(StdI->L % 8 != 0){
    fprintf(stderr, "\n ERROR! L % 8 != 0 \n\n");
    exitMPI(-1);
  }
  StdI->W = StdI->L / 2;
  StdI->L = 2;
  StdI->num_pivot = StdI->W / 4;
/**/
  fprintf(fp, "# W0  R0  StdI->num_pivot  StdI->ishift_nspin\n");
  fprintf(fp, "%d %d %d %d\n", StdI->W, StdI->L, StdI->num_pivot, StdI->ishift_nspin);

  StdI->list_6spin_star = (int **)malloc(sizeof(int*) * StdI->num_pivot);
  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    StdI->list_6spin_star[ipivot] = (int *)malloc(sizeof(int) * 7);
  }

  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    StdI->list_6spin_star[ipivot][0] = 8; // num of J
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
    StdI->list_6spin_pair[ipivot][0][1] = 1;
    StdI->list_6spin_pair[ipivot][1][1] = 2;
    StdI->list_6spin_pair[ipivot][2][1] = 0;
    StdI->list_6spin_pair[ipivot][3][1] = 3;
    StdI->list_6spin_pair[ipivot][4][1] = 4;
    StdI->list_6spin_pair[ipivot][5][1] = 5;
    StdI->list_6spin_pair[ipivot][6][1] = 1; // type of J
    StdI->list_6spin_pair[ipivot][0][2] = 2;
    StdI->list_6spin_pair[ipivot][1][2] = 3;
    StdI->list_6spin_pair[ipivot][2][2] = 0;
    StdI->list_6spin_pair[ipivot][3][2] = 1;
    StdI->list_6spin_pair[ipivot][4][2] = 4;
    StdI->list_6spin_pair[ipivot][5][2] = 5;
    StdI->list_6spin_pair[ipivot][6][2] = 1; // type of J
    StdI->list_6spin_pair[ipivot][0][3] = 3;
    StdI->list_6spin_pair[ipivot][1][3] = 4;
    StdI->list_6spin_pair[ipivot][2][3] = 0;
    StdI->list_6spin_pair[ipivot][3][3] = 1;
    StdI->list_6spin_pair[ipivot][4][3] = 2;
    StdI->list_6spin_pair[ipivot][5][3] = 5;
    StdI->list_6spin_pair[ipivot][6][3] = 1; // type of J
    StdI->list_6spin_pair[ipivot][0][4] = 0;
    StdI->list_6spin_pair[ipivot][1][4] = 2;
    StdI->list_6spin_pair[ipivot][2][4] = 1;
    StdI->list_6spin_pair[ipivot][3][4] = 3;
    StdI->list_6spin_pair[ipivot][4][4] = 4;
    StdI->list_6spin_pair[ipivot][5][4] = 5;
    StdI->list_6spin_pair[ipivot][6][4] = 2; // type of J
    StdI->list_6spin_pair[ipivot][0][5] = 1;
    StdI->list_6spin_pair[ipivot][1][5] = 3;
    StdI->list_6spin_pair[ipivot][2][5] = 0;
    StdI->list_6spin_pair[ipivot][3][5] = 2;
    StdI->list_6spin_pair[ipivot][4][5] = 4;
    StdI->list_6spin_pair[ipivot][5][5] = 5;
    StdI->list_6spin_pair[ipivot][6][5] = 2; // type of J
    StdI->list_6spin_pair[ipivot][0][6] = 2;
    StdI->list_6spin_pair[ipivot][1][6] = 4;
    StdI->list_6spin_pair[ipivot][2][6] = 0;
    StdI->list_6spin_pair[ipivot][3][6] = 1;
    StdI->list_6spin_pair[ipivot][4][6] = 3;
    StdI->list_6spin_pair[ipivot][5][6] = 5;
    StdI->list_6spin_pair[ipivot][6][6] = 2; // type of J
    StdI->list_6spin_pair[ipivot][0][7] = 3;
    StdI->list_6spin_pair[ipivot][1][7] = 5;
    StdI->list_6spin_pair[ipivot][2][7] = 0;
    StdI->list_6spin_pair[ipivot][3][7] = 1;
    StdI->list_6spin_pair[ipivot][4][7] = 2;
    StdI->list_6spin_pair[ipivot][5][7] = 4;
    StdI->list_6spin_pair[ipivot][6][7] = 2; // type of J
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
