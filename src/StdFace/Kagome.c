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
* Setup a Hamiltonian for the generalized Heisenberg model on a Heisenberg lattice
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void Spin_Kagome_Boost(
struct StdIntList *StdI,
  int Sz2 /**< [in] 2 * Total Sz */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int iL, iW, j;
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
  StdFace_PrintVal_d("Jxy", &StdI->Jxy, StdI->J);
  StdFace_PrintVal_d("Jx", &StdI->Jx, StdI->Jxy);
  StdFace_PrintVal_d("Jy", &StdI->Jy, StdI->Jxy);
  StdFace_PrintVal_d("Jz", &StdI->Jz, StdI->J);
  /**/
  StdFace_NotUsed_d("J'", StdI->Jp);
  StdFace_NotUsed_d("Jz'", StdI->Jzp);
  StdFace_NotUsed_d("Jxy'", StdI->Jxyp);
  StdFace_NotUsed_d("Jx'", StdI->Jxp);
  StdFace_NotUsed_d("Jy'", StdI->Jyp);
  StdFace_NotUsed_d("K", StdI->K);
  StdFace_NotUsed_d("J0", StdI->J0);
  StdFace_NotUsed_d("J1", StdI->J1);
  StdFace_NotUsed_d("J2", StdI->J2);
  StdFace_NotUsed_d("Jz0", StdI->Jz0);
  StdFace_NotUsed_d("Jz1", StdI->Jz1);
  StdFace_NotUsed_d("Jz2", StdI->Jz2);
  StdFace_NotUsed_d("Jxy0", StdI->Jxy0);
  StdFace_NotUsed_d("Jxy1", StdI->Jxy1);
  StdFace_NotUsed_d("Jxy2", StdI->Jxy2);
  StdFace_NotUsed_d("Jx0", StdI->Jx0);
  StdFace_NotUsed_d("Jx1", StdI->Jx1);
  StdFace_NotUsed_d("Jx2", StdI->Jx2);
  StdFace_NotUsed_d("Jy0", StdI->Jy0);
  StdFace_NotUsed_d("Jy1", StdI->Jy1);
  StdFace_NotUsed_d("Jy2", StdI->Jy2);
  /*
  Local Spin
  */
  StdI->nsite = StdI->L * StdI->W;
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
    LargeValue0 = (double)Sz2 / (double)(2 * StdI->nsite) * fabs(StdI->h) + S * fabs(StdI->D)
      + S * S * fabs(StdI->Gamma)
      + 4.0 / 2.0 * S * S * (fabs(StdI->Jx) + fabs(StdI->Jy) + fabs(StdI->Jz));
  }
  else {
    LargeValue0 = S * fabs(StdI->h) + S * fabs(StdI->D)
      + S * S * fabs(StdI->Gamma)
      + 4.0 / 2.0 * S * S * (fabs(StdI->Jx) + fabs(StdI->Jy) + fabs(StdI->Jz));
  }
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
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
  fprintf(fp, "%d  # Number of type of J\n", 1);
  fprintf(fp, "# J 1\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * StdI->Jx, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * StdI->Jy, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * StdI->Jz, 0.0);
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
  StdI->num_pivot = 4;
  if (StdI->W != 3) {
    fprintf(stderr, "DEBUG: W != 3\n");
    exit(-1);
  }
  StdI->W = 9;
  fprintf(fp, "# W0  R0  StdI->num_pivot  StdI->ishift_nspin\n");
  fprintf(fp, "%d %d %d %d\n", StdI->W, StdI->L, StdI->num_pivot, StdI->ishift_nspin);

  StdI->list_6spin_star = (int **)malloc(sizeof(int*) * StdI->num_pivot);
  for (j = 0; j < StdI->num_pivot; j++) {
    StdI->list_6spin_star[j] = (int *)malloc(sizeof(int) * 7);
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
  for (j = 0; j < StdI->num_pivot; j++) {
    fprintf(fp, "# pivot %d\n", j);
    for (iW = 0; iW < 7; iW++) {
      fprintf(fp, "%d ", StdI->list_6spin_star[j][iW]);
    }
    fprintf(fp, "\n");
  }

  StdI->list_6spin_pair = (int ***)malloc(sizeof(int**) * StdI->num_pivot);
  for (j = 0; j < StdI->num_pivot; j++) {
    StdI->list_6spin_pair[j] = (int **)malloc(sizeof(int*) * 7);
    for (iW = 0; iW < 7; iW++) {
      StdI->list_6spin_pair[j][iW] = (int *)malloc(sizeof(int) * StdI->list_6spin_star[j][0]);
    }
  }

  StdI->list_6spin_pair[0][0][0] = 0; //(1,1,1+2*j)=0 
  StdI->list_6spin_pair[0][1][0] = 4; //(2,1,1+2*j)=1
  StdI->list_6spin_pair[0][2][0] = 1; //(3,1,1+2*j)=2
  StdI->list_6spin_pair[0][3][0] = 2; //(4,1,1+2*j)=3
  StdI->list_6spin_pair[0][4][0] = 3; //(5,1,1+2*j)=4
  StdI->list_6spin_pair[0][5][0] = 5; //(6,1,1+2*j)=5
  StdI->list_6spin_pair[0][6][0] = 1; //(7,1,1+2*j)=3 ! type of J

  StdI->list_6spin_pair[1][0][0] = 0;
  StdI->list_6spin_pair[1][1][0] = 1;
  StdI->list_6spin_pair[1][2][0] = 2;
  StdI->list_6spin_pair[1][3][0] = 3;
  StdI->list_6spin_pair[1][4][0] = 4;
  StdI->list_6spin_pair[1][5][0] = 5;
  StdI->list_6spin_pair[1][6][0] = 1; // type of J
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
  StdI->list_6spin_pair[1][6][2] = 1; // type of J
  StdI->list_6spin_pair[1][0][3] = 1;
  StdI->list_6spin_pair[1][1][3] = 3;
  StdI->list_6spin_pair[1][2][3] = 0;
  StdI->list_6spin_pair[1][3][3] = 2;
  StdI->list_6spin_pair[1][4][3] = 4;
  StdI->list_6spin_pair[1][5][3] = 5;
  StdI->list_6spin_pair[1][6][3] = 1; // type of J
  StdI->list_6spin_pair[1][0][4] = 2;
  StdI->list_6spin_pair[1][1][4] = 4;
  StdI->list_6spin_pair[1][2][4] = 0;
  StdI->list_6spin_pair[1][3][4] = 1;
  StdI->list_6spin_pair[1][4][4] = 3;
  StdI->list_6spin_pair[1][5][4] = 5;
  StdI->list_6spin_pair[1][6][4] = 1; // type of J
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
  StdI->list_6spin_pair[2][6][0] = 1; // type of J
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
  StdI->list_6spin_pair[2][6][2] = 1; // type of J
  StdI->list_6spin_pair[2][0][3] = 1;
  StdI->list_6spin_pair[2][1][3] = 3;
  StdI->list_6spin_pair[2][2][3] = 0;
  StdI->list_6spin_pair[2][3][3] = 2;
  StdI->list_6spin_pair[2][4][3] = 4;
  StdI->list_6spin_pair[2][5][3] = 5;
  StdI->list_6spin_pair[2][6][3] = 1; // type of J
  StdI->list_6spin_pair[2][0][4] = 2;
  StdI->list_6spin_pair[2][1][4] = 5;
  StdI->list_6spin_pair[2][2][4] = 0;
  StdI->list_6spin_pair[2][3][4] = 1;
  StdI->list_6spin_pair[2][4][4] = 3;
  StdI->list_6spin_pair[2][5][4] = 4;
  StdI->list_6spin_pair[2][6][4] = 1; // type of J
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
  StdI->list_6spin_pair[3][6][0] = 1; // type of J
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
  StdI->list_6spin_pair[3][6][2] = 1; // type of J
  StdI->list_6spin_pair[3][0][3] = 2;
  StdI->list_6spin_pair[3][1][3] = 5;
  StdI->list_6spin_pair[3][2][3] = 0;
  StdI->list_6spin_pair[3][3][3] = 1;
  StdI->list_6spin_pair[3][4][3] = 3;
  StdI->list_6spin_pair[3][5][3] = 4;
  StdI->list_6spin_pair[3][6][3] = 1; // type of J
  StdI->list_6spin_pair[3][0][4] = 2;
  StdI->list_6spin_pair[3][1][4] = 4;
  StdI->list_6spin_pair[3][2][4] = 0;
  StdI->list_6spin_pair[3][3][4] = 1;
  StdI->list_6spin_pair[3][4][4] = 3;
  StdI->list_6spin_pair[3][5][4] = 5;
  StdI->list_6spin_pair[3][6][4] = 1; // type of J

  fprintf(fp, "# StdI->list_6spin_pair\n");
  for (j = 0; j < StdI->num_pivot; j++) {
    fprintf(fp, "# pivot %d\n", j);
    for (iL = 0; iL < StdI->list_6spin_star[j][0]; iL++) {
      for (iW = 0; iW < 7; iW++) {
        fprintf(fp, "%d ", StdI->list_6spin_pair[j][iW][iL]);
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);

  for (j = 0; j < StdI->num_pivot; j++) {
    free(StdI->list_6spin_star[j]);
  }
  free(StdI->list_6spin_star);

  for (j = 0; j < StdI->num_pivot; j++) {
    for (iW = 0; iW < 7; iW++) {
      free(StdI->list_6spin_pair[j][iW]);
    }
    free(StdI->list_6spin_pair[j]);
  }
  free(StdI->list_6spin_pair);

}
