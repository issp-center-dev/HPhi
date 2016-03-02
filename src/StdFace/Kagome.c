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
#include "../include/wrapperMPI.h"

/**
*
* Setup a Hamiltonian for the generalized Heisenberg model on a Heisenberg lattice
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void Spin_Kagome_Boost(
  int Sz2 /**< [in] 2 * Total Sz */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int iL, iW, j;
  int ktrans, kintr;
  double LargeValue0, S;
  FILE *fp;

  fprintf(stdoutMPI, "\n");
  fprintf(stdoutMPI, "#######  Parameter Summary  #######\n");
  fprintf(stdoutMPI, "\n");
  StdFace_RequiredVal_i("L", L);
  StdFace_RequiredVal_i("W", W);
  StdFace_PrintVal_d("a", &a, 1.0);
  /**/
  StdFace_PrintVal_i("2S", &S2, 1);
  StdFace_PrintVal_d("h", &h, 0.0);
  StdFace_PrintVal_d("Gamma", &Gamma, 0.0);
  StdFace_PrintVal_d("D", &D, 0.0);
  StdFace_PrintVal_d("J", &J, 0.0);
  StdFace_PrintVal_d("Jxy", &Jxy, J);
  StdFace_PrintVal_d("Jx", &Jx, Jxy);
  StdFace_PrintVal_d("Jy", &Jy, Jxy);
  StdFace_PrintVal_d("Jz", &Jz, J);
  /**/
  StdFace_NotUsed_d("J'", Jp);
  StdFace_NotUsed_d("Jz'", Jzp);
  StdFace_NotUsed_d("Jxy'", Jxyp);
  StdFace_NotUsed_d("Jx'", Jxp);
  StdFace_NotUsed_d("Jy'", Jyp);
  StdFace_NotUsed_d("K", K);
  StdFace_NotUsed_d("J0", J0);
  StdFace_NotUsed_d("J1", J1);
  StdFace_NotUsed_d("J2", J2);
  StdFace_NotUsed_d("Jz0", Jz0);
  StdFace_NotUsed_d("Jz1", Jz1);
  StdFace_NotUsed_d("Jz2", Jz2);
  StdFace_NotUsed_d("Jxy0", Jxy0);
  StdFace_NotUsed_d("Jxy1", Jxy1);
  StdFace_NotUsed_d("Jxy2", Jxy2);
  StdFace_NotUsed_d("Jx0", Jx0);
  StdFace_NotUsed_d("Jx1", Jx1);
  StdFace_NotUsed_d("Jx2", Jx2);
  StdFace_NotUsed_d("Jy0", Jy0);
  StdFace_NotUsed_d("Jy1", Jy1);
  StdFace_NotUsed_d("Jy2", Jy2);
  /*
  Local Spin
  */
  nsite = L * W;
  S2 = 1;
  locspinflag = (int *)malloc(sizeof(int) * nsite);
  for (isite = 0; isite < nsite; isite++)locspinflag[isite] = S2;
  /*
  Transfer
  */
  ntrans = 1;
  transindx = (int **)malloc(sizeof(int*) * ntrans);
  trans = (double *)malloc(sizeof(double) * ntrans);
  for (ktrans = 0; ktrans < ntrans; ktrans++) {
    transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }
  ntrans = 0;
  /*
  Interaction
  */
  nintr = 1;
  intrindx = (int **)malloc(sizeof(int*) * nintr);
  intr = (double *)malloc(sizeof(double) * nintr);
  for (kintr = 0; kintr < nintr; kintr++) {
    intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  nintr = 0;
  /*
  Set mTPQ parameter
  */
  S = (double)S2 * 0.5;
  if (lGC == 0) {
    LargeValue0 = (double)Sz2 / (double)(2 * nsite) * fabs(h) + S * fabs(D)
      + S * S * fabs(Gamma)
      + 4.0 / 2.0 * S * S * (fabs(Jx) + fabs(Jy) + fabs(Jz));
  }
  else {
    LargeValue0 = S * fabs(h) + S * fabs(D) + S * S * fabs(Gamma)
      + 4.0 / 2.0 * S * S * (fabs(Jx) + fabs(Jy) + fabs(Jz));
  }
  StdFace_PrintVal_d("LargeValue", &LargeValue, LargeValue0);
  /*
  Magnetic field
  */
  fp = fopenMPI("boost.def", "w");
  fprintf(fp, "# Magnetic field\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    -0.5 * Gamma, 0.0, 0.0, 0.0, -0.5 * h, 0.0);
  /*
  Interaction
  */
  fprintf(fp, "%d  # Number of type of J\n", 1);
  fprintf(fp, "# J 1\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * Jx, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * Jy, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * Jz, 0.0);
  /*
  Topology
  */
  if (S2 != 1) {
    fprintf(stderr, "\n ERROR! S2 must be 1 in Boost. \n\n");
    exitMPI(-1);
  }
  ishift_nspin = 3;
  if (L < 2) {
    fprintf(stderr, "\n ERROR! L < 2 \n\n");
    exitMPI(-1);
  }
  if (W % ishift_nspin != 0) {
    fprintf(stderr, "\n ERROR! W %% %d != 0 \n\n", ishift_nspin);
    exitMPI(-1);
  }
  num_pivot = 4;
  if (W != 9) {
    fprintf(stderr, "DEBUG: W != 9\n");
    exitMPI(-1);
  }
  fprintf(fp, "# W0  R0  num_pivot  ishift_nspin\n");
  fprintf(fp, "%d %d %d %d\n", W, L, num_pivot, ishift_nspin);

  list_6spin_star = (int **)malloc(sizeof(int*) * num_pivot);
  for (j = 0; j < num_pivot; j++) {
    list_6spin_star[j] = (int *)malloc(sizeof(int) * 7);
  }

  list_6spin_star[0][0] = 1; // num of J
  list_6spin_star[0][1] = 1;
  list_6spin_star[0][2] = 1;
  list_6spin_star[0][3] = 1;
  list_6spin_star[0][4] = 4;
  list_6spin_star[0][5] = 2;
  list_6spin_star[0][6] = -1; // flag

  list_6spin_star[1][0] = 6; // num of J
  list_6spin_star[1][1] = 1;
  list_6spin_star[1][2] = 1;
  list_6spin_star[1][3] = 1;
  list_6spin_star[1][4] = 6;
  list_6spin_star[1][5] = 7;
  list_6spin_star[1][6] = 1; // flag

  list_6spin_star[2][0] = 6; // num of J
  list_6spin_star[2][1] = 1;
  list_6spin_star[2][2] = 1;
  list_6spin_star[2][3] = 1;
  list_6spin_star[2][4] = 4;
  list_6spin_star[2][5] = 2;
  list_6spin_star[2][6] = 1; // flag

  list_6spin_star[3][0] = 5; // num of J
  list_6spin_star[3][1] = 1;
  list_6spin_star[3][2] = 1;
  list_6spin_star[3][3] = 1;
  list_6spin_star[3][4] = 4;
  list_6spin_star[3][5] = 2;
  list_6spin_star[3][6] = 1; // flag

  fprintf(fp, "# list_6spin_star\n");
  for (j = 0; j < num_pivot; j++) {
    fprintf(fp, "# pivot %d\n", j);
    for (iW = 0; iW < 7; iW++) {
      fprintf(fp, "%d ", list_6spin_star[j][iW]);
    }
    fprintf(fp, "\n");
  }

  list_6spin_pair = (int ***)malloc(sizeof(int**) * num_pivot);
  for (j = 0; j < num_pivot; j++) {
    list_6spin_pair[j] = (int **)malloc(sizeof(int*) * 7);
    for (iW = 0; iW < 7; iW++) {
      list_6spin_pair[j][iW] = (int *)malloc(sizeof(int) * list_6spin_star[j][0]);
    }
  }

  list_6spin_pair[0][0][0] = 0; //(1,1,1+2*j)=0 
  list_6spin_pair[0][1][0] = 4; //(2,1,1+2*j)=1
  list_6spin_pair[0][2][0] = 1; //(3,1,1+2*j)=2
  list_6spin_pair[0][3][0] = 2; //(4,1,1+2*j)=3
  list_6spin_pair[0][4][0] = 3; //(5,1,1+2*j)=4
  list_6spin_pair[0][5][0] = 5; //(6,1,1+2*j)=5
  list_6spin_pair[0][6][0] = 1; //(7,1,1+2*j)=3 ! type of J

  list_6spin_pair[1][0][0] = 0;
  list_6spin_pair[1][1][0] = 1;
  list_6spin_pair[1][2][0] = 2;
  list_6spin_pair[1][3][0] = 3;
  list_6spin_pair[1][4][0] = 4;
  list_6spin_pair[1][5][0] = 5;
  list_6spin_pair[1][6][0] = 1; // type of J
  list_6spin_pair[1][0][1] = 1;
  list_6spin_pair[1][1][1] = 2;
  list_6spin_pair[1][2][1] = 0;
  list_6spin_pair[1][3][1] = 3;
  list_6spin_pair[1][4][1] = 4;
  list_6spin_pair[1][5][1] = 5;
  list_6spin_pair[1][6][1] = 1; // type of J
  list_6spin_pair[1][0][2] = 0;
  list_6spin_pair[1][1][2] = 2;
  list_6spin_pair[1][2][2] = 1;
  list_6spin_pair[1][3][2] = 3;
  list_6spin_pair[1][4][2] = 4;
  list_6spin_pair[1][5][2] = 5;
  list_6spin_pair[1][6][2] = 1; // type of J
  list_6spin_pair[1][0][3] = 1;
  list_6spin_pair[1][1][3] = 3;
  list_6spin_pair[1][2][3] = 0;
  list_6spin_pair[1][3][3] = 2;
  list_6spin_pair[1][4][3] = 4;
  list_6spin_pair[1][5][3] = 5;
  list_6spin_pair[1][6][3] = 1; // type of J
  list_6spin_pair[1][0][4] = 2;
  list_6spin_pair[1][1][4] = 4;
  list_6spin_pair[1][2][4] = 0;
  list_6spin_pair[1][3][4] = 1;
  list_6spin_pair[1][4][4] = 3;
  list_6spin_pair[1][5][4] = 5;
  list_6spin_pair[1][6][4] = 1; // type of J
  list_6spin_pair[1][0][5] = 2;
  list_6spin_pair[1][1][5] = 5;
  list_6spin_pair[1][2][5] = 0;
  list_6spin_pair[1][3][5] = 1;
  list_6spin_pair[1][4][5] = 3;
  list_6spin_pair[1][5][5] = 4;
  list_6spin_pair[1][6][5] = 1; // type of J

  list_6spin_pair[2][0][0] = 0;
  list_6spin_pair[2][1][0] = 1;
  list_6spin_pair[2][2][0] = 2;
  list_6spin_pair[2][3][0] = 3;
  list_6spin_pair[2][4][0] = 4;
  list_6spin_pair[2][5][0] = 5;
  list_6spin_pair[2][6][0] = 1; // type of J
  list_6spin_pair[2][0][1] = 1;
  list_6spin_pair[2][1][1] = 2;
  list_6spin_pair[2][2][1] = 0;
  list_6spin_pair[2][3][1] = 3;
  list_6spin_pair[2][4][1] = 4;
  list_6spin_pair[2][5][1] = 5;
  list_6spin_pair[2][6][1] = 1; // type of J
  list_6spin_pair[2][0][2] = 0;
  list_6spin_pair[2][1][2] = 2;
  list_6spin_pair[2][2][2] = 1;
  list_6spin_pair[2][3][2] = 3;
  list_6spin_pair[2][4][2] = 4;
  list_6spin_pair[2][5][2] = 5;
  list_6spin_pair[2][6][2] = 1; // type of J
  list_6spin_pair[2][0][3] = 1;
  list_6spin_pair[2][1][3] = 3;
  list_6spin_pair[2][2][3] = 0;
  list_6spin_pair[2][3][3] = 2;
  list_6spin_pair[2][4][3] = 4;
  list_6spin_pair[2][5][3] = 5;
  list_6spin_pair[2][6][3] = 1; // type of J
  list_6spin_pair[2][0][4] = 2;
  list_6spin_pair[2][1][4] = 5;
  list_6spin_pair[2][2][4] = 0;
  list_6spin_pair[2][3][4] = 1;
  list_6spin_pair[2][4][4] = 3;
  list_6spin_pair[2][5][4] = 4;
  list_6spin_pair[2][6][4] = 1; // type of J
  list_6spin_pair[2][0][5] = 2;
  list_6spin_pair[2][1][5] = 4;
  list_6spin_pair[2][2][5] = 0;
  list_6spin_pair[2][3][5] = 1;
  list_6spin_pair[2][4][5] = 3;
  list_6spin_pair[2][5][5] = 5;
  list_6spin_pair[2][6][5] = 1; // type of J

  list_6spin_pair[3][0][0] = 0;
  list_6spin_pair[3][1][0] = 1;
  list_6spin_pair[3][2][0] = 2;
  list_6spin_pair[3][3][0] = 3;
  list_6spin_pair[3][4][0] = 4;
  list_6spin_pair[3][5][0] = 5;
  list_6spin_pair[3][6][0] = 1; // type of J
  list_6spin_pair[3][0][1] = 1;
  list_6spin_pair[3][1][1] = 2;
  list_6spin_pair[3][2][1] = 0;
  list_6spin_pair[3][3][1] = 3;
  list_6spin_pair[3][4][1] = 4;
  list_6spin_pair[3][5][1] = 5;
  list_6spin_pair[3][6][1] = 1; // type of J
  list_6spin_pair[3][0][2] = 0;
  list_6spin_pair[3][1][2] = 2;
  list_6spin_pair[3][2][2] = 1;
  list_6spin_pair[3][3][2] = 3;
  list_6spin_pair[3][4][2] = 4;
  list_6spin_pair[3][5][2] = 5;
  list_6spin_pair[3][6][2] = 1; // type of J
  list_6spin_pair[3][0][3] = 2;
  list_6spin_pair[3][1][3] = 5;
  list_6spin_pair[3][2][3] = 0;
  list_6spin_pair[3][3][3] = 1;
  list_6spin_pair[3][4][3] = 3;
  list_6spin_pair[3][5][3] = 4;
  list_6spin_pair[3][6][3] = 1; // type of J
  list_6spin_pair[3][0][4] = 2;
  list_6spin_pair[3][1][4] = 4;
  list_6spin_pair[3][2][4] = 0;
  list_6spin_pair[3][3][4] = 1;
  list_6spin_pair[3][4][4] = 3;
  list_6spin_pair[3][5][4] = 5;
  list_6spin_pair[3][6][4] = 1; // type of J

  fprintf(fp, "# list_6spin_pair\n");
  for (j = 0; j < num_pivot; j++) {
    fprintf(fp, "# pivot %d\n", j);
    for (iL = 0; iL < list_6spin_star[j][0]; iL++) {
      for (iW = 0; iW < 7; iW++) {
        fprintf(fp, "%d ", list_6spin_pair[j][iW][iL]);
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);

  for (j = 0; j < num_pivot; j++) {
    free(list_6spin_star[j]);
  }
  free(list_6spin_star);

  for (j = 0; j < num_pivot; j++) {
    for (iW = 0; iW < 7; iW++) {
      free(list_6spin_pair[j][iW]);
    }
    free(list_6spin_pair[j]);
  }
  free(list_6spin_pair);

}
