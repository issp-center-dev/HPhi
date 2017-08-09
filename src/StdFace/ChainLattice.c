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
/**@file
@brief Standard mode for the chain lattice
*/
#include "StdFace_vals.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "StdFace_ModelUtil.h"
#include <complex.h>
#include <string.h>

/**
@brief Setup a Hamiltonian for the Hubbard model on a Chain lattice
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_Chain(
  struct StdIntList *StdI//!<[inout]
)
{
  FILE *fp;
  int isite, jsite, ntransMax, nintrMax;
  int iL;
  double complex Cphase;
  
  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  /**@brief
  (1) Compute the shape of the super-cell and sites in the super-cell
  */
  fp = fopen("lattice.gp", "w");
  /**/
  StdI->NsiteUC = 1;
  /**/
  fprintf(stdout, "  @ Lattice Size & Shape\n\n");
  
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  StdFace_PrintVal_d("Wlength", &StdI->length[0], StdI->a);
  StdFace_PrintVal_d("Llength", &StdI->length[1], StdI->a);
  StdFace_PrintVal_d("Wx", &StdI->direct[0][0], StdI->length[0]);
  StdFace_PrintVal_d("Wy", &StdI->direct[0][1], 0.0);
  StdFace_PrintVal_d("Lx", &StdI->direct[1][0], 0.0);
  StdFace_PrintVal_d("Ly", &StdI->direct[1][1], StdI->length[1]);

  StdFace_PrintVal_d("phase0", &StdI->phase[0], 0.0);
  StdFace_NotUsed_d("phase1", StdI->phase[1]);
  StdI->phase[1] = StdI->phase[0];
  StdI->phase[0] = 0.0;
  /**/
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_NotUsed_i("W", StdI->W);
  StdI->W = 1;
  StdFace_InitSite(StdI, fp, 2);
  StdI->tau[0][0] = 0.0; StdI->tau[0][1] = 0.0; StdI->tau[0][2] = 0.0;
  /**@brief
  (2) check & store parameters of Hamiltonian
  */
  fprintf(stdout, "\n  @ Hamiltonian \n\n");
  StdFace_NotUsed_J("J1", StdI->J1All, StdI->J1);
  StdFace_NotUsed_J("J2", StdI->J2All, StdI->J2);
  StdFace_NotUsed_J("J1'", StdI->J1pAll, StdI->J1p);
  StdFace_NotUsed_J("J2'", StdI->J2pAll, StdI->J2p);
  StdFace_NotUsed_c("t1", StdI->t1);
  StdFace_NotUsed_c("t2", StdI->t2);
  StdFace_NotUsed_d("t1'", StdI->t1p);
  StdFace_NotUsed_d("t2'", StdI->t2p);
  StdFace_NotUsed_d("V1", StdI->V1);
  StdFace_NotUsed_d("V2", StdI->V2);
  StdFace_NotUsed_d("V1'", StdI->V1p);
  StdFace_NotUsed_d("V2'", StdI->V2p);
  StdFace_NotUsed_d("K", StdI->K);
  StdFace_PrintVal_d("h", &StdI->h, 0.0);
  StdFace_PrintVal_d("Gamma", &StdI->Gamma, 0.0);
  /**/
  if (strcmp(StdI->model, "spin") == 0 ) {
    StdFace_PrintVal_i("2S", &StdI->S2, 1);
    StdFace_PrintVal_d("D", &StdI->D[2][2], 0.0);
    StdFace_InputSpinNN(StdI, StdI->J0, StdI->J0All, "J0");
    StdFace_InputSpin(StdI, StdI->Jp, StdI->JpAll, "J'");
    /**/
    StdFace_NotUsed_d("mu", StdI->mu);
    StdFace_NotUsed_d("U", StdI->U);
    StdFace_NotUsed_c("t", StdI->t);
    StdFace_NotUsed_c("t0", StdI->t0);
    StdFace_NotUsed_c("t'", StdI->tp);
    StdFace_NotUsed_d("V", StdI->V);
    StdFace_NotUsed_d("V0", StdI->V0);
    StdFace_NotUsed_d("V'", StdI->Vp);
  }/*if (strcmp(StdI->model, "spin") == 0 )*/
  else {
    StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
    StdFace_PrintVal_d("U", &StdI->U, 0.0);
    StdFace_InputHopp(StdI, &StdI->t0, "t0");
    StdFace_PrintVal_c("t'", &StdI->tp, 0.0);
    StdFace_InputCoulombV(StdI, &StdI->V0, "V0");
    StdFace_PrintVal_d("V'", &StdI->Vp, 0.0);

    StdFace_NotUsed_J("J0", StdI->J0All, StdI->J0);
    StdFace_NotUsed_J("J'", StdI->JpAll, StdI->Jp);
    StdFace_NotUsed_d("D", StdI->D[2][2]);

    if (strcmp(StdI->model, "hubbard") == 0 ) {
      StdFace_NotUsed_i("2S", StdI->S2);
      StdFace_NotUsed_J("J", StdI->JAll, StdI->J);
    }
    else if (strcmp(StdI->model, "kondo") == 0 ) {
      StdFace_PrintVal_i("2S", &StdI->S2, 1);
      StdFace_InputSpin(StdI, StdI->J, StdI->JAll, "J");
    }
  }/*if (strcmp(StdI->model, "spin") != 0 )*/
  fprintf(stdout, "\n  @ Numerical conditions\n\n");
  /**@brief
  (3) Set local spin flag (StdIntList::locspinflag) and
  the number of sites (StdIntList::nsite)
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
  /**@brief
  (4) Compute the upper limit of the number of Transfer & Interaction and malloc them.
  */
  if (strcmp(StdI->model, "spin") == 0 ) {
    ntransMax = StdI->L * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
    nintrMax = StdI->L * (StdI->NsiteUC/*D*/ + 1/*J*/ + 1/*J'*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  else {
    ntransMax = StdI->L * 2/*spin*/ * (2 * StdI->NsiteUC/*mu+h+Gamma*/ + 2/*t*/ + 2/*t'*/);
    nintrMax = StdI->L * (StdI->NsiteUC/*U*/ + 4 * (1/*V*/ + 1/*V'*/));

    if (strcmp(StdI->model, "kondo") == 0) {
      ntransMax += StdI->L * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
      nintrMax += StdI->nsite / 2 * (3 * 1 + 1) * (3 * StdI->S2 + 1);
    }/*if (strcmp(StdI->model, "kondo") == 0)*/
  }
  /**/
  StdFace_MallocInteractions(StdI, ntransMax, nintrMax);
  /**@brief
  (5) Set Transfer & Interaction
  */
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
      StdFace_HubbardLocal(StdI, StdI->mu, -StdI->h, -StdI->Gamma, StdI->U, isite);
      if (strcmp(StdI->model, "kondo") == 0 ) {
        jsite = iL;
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite, jsite);
        StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, jsite);
      }/*if (strcmp(StdI->model, "kondo") == 0 )*/
    }/*if (model != "spin")*/
    /*
     Nearest neighbor
    */
    StdFace_SetLabel(StdI, fp, 0, iL, 0, 1, 0, 0, &isite, &jsite, 1, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J0, StdI->S2, StdI->S2, isite, jsite);
    }
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t0, isite, jsite);
      StdFace_Coulomb(StdI, StdI->V0, isite, jsite);
    }
    /*
     Second nearest neighbor
    */
    StdFace_SetLabel(StdI, fp, 0, iL, 0, 2, 0, 0, &isite, &jsite, 2, &Cphase);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }
    else {
      StdFace_Hopping(StdI, Cphase * StdI->tp, isite, jsite);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
  }/*for (iL = 0; iL < StdI->L; iL++)*/

  fprintf(fp, "plot \'-\' w d lc 7\n0.0 0.0\nend\npause -1\n");
  fclose(fp);
  StdFace_PrintGeometry(StdI);
}/*void StdFace_Chain*/

#if defined(_HPhi)
/**
@brief Setup a Hamiltonian for the generalized Heisenberg model on a Chain lattice
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_Chain_Boost(struct StdIntList *StdI)
{
  int isite, ipivot;
  int kintr;
  FILE *fp;

  StdI->NsiteUC = 1;
  /*
  Magnetic field
  */
  fp = fopen("boost.def", "w");
  fprintf(fp, "# Magnetic field\n");
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    -0.5 * StdI->Gamma, 0.0, -0.5 *StdI->h);
  /*
  Interaction
  */
  fprintf(fp, "%d  # Number of type of J\n", 2);
  fprintf(fp, "# J 1\n");
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J0[0][0], 0.25 * StdI->J0[0][1], 0.25 * StdI->J0[0][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J0[1][0], 0.25 * StdI->J0[1][1], 0.25 * StdI->J0[1][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J0[2][0], 0.25 * StdI->J0[2][1], 0.25 * StdI->J0[2][2]);
  fprintf(fp, "# J 2\n");
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->Jp[0][0], 0.25 * StdI->Jp[0][1], 0.25 * StdI->Jp[0][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->Jp[1][0], 0.25 * StdI->Jp[1][1], 0.25 * StdI->Jp[1][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->Jp[2][0], 0.25 * StdI->Jp[2][1], 0.25 * StdI->Jp[2][2]);
  /*
  Topology
  */
  if (StdI->S2 != 1) {
    fprintf(stdout, "\n ERROR! S2 must be 1 in Boost. \n\n");
    StdFace_exit(-1);
  }
  StdI->ishift_nspin = 4;
  if(StdI->L % 8 != 0){
    fprintf(stdout, "\n ERROR! L %% 8 != 0 \n\n");
    StdFace_exit(-1);
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
#endif
