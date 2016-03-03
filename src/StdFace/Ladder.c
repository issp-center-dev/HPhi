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
 * Setup a Hamiltonian for the Hubbard model on a square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void FermionHubbard_Ladder(
  int nelec /**< [in] The number of electrons */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int ispin;
  int iL, iW;
  int ktrans, kintr;
  double LargeValue0;
  fprintf(stdoutMPI, "\n");
  fprintf(stdoutMPI, "#######  Parameter Summary  #######\n");
  fprintf(stdoutMPI, "\n");
  StdFace_RequiredVal_i("L", L);
  StdFace_RequiredVal_i("W", W);
  StdFace_PrintVal_d("a", &a, 1.0);
  /**/
  StdFace_PrintVal_d("mu", &mu, 0.0);
  StdFace_PrintVal_d("U", &U, 0.0);
  StdFace_PrintVal_d("t1", &t1, 1.0);
  StdFace_PrintVal_d("t1'", &t1p, 1.0);
  StdFace_PrintVal_d("t0", &t0, 1.0);
  StdFace_PrintVal_d("t2", &t2, 1.0);
  StdFace_PrintVal_d("t2", &t2p, 1.0);
  StdFace_PrintVal_d("V0", &V0, 0.0);
  StdFace_PrintVal_d("V1", &V1, 0.0);
  StdFace_PrintVal_d("V1'", &V1p, 0.0);
  StdFace_PrintVal_d("V2", &V2, 0.0);
  /**/
  StdFace_NotUsed_i("2S", S2);
  StdFace_NotUsed_d("t", t);
  StdFace_NotUsed_d("t'", tp);
  StdFace_NotUsed_d("t''", tpp);
  StdFace_NotUsed_d("V", V);
  StdFace_NotUsed_d("V'", Vp);
  StdFace_NotUsed_d("V''", Vpp);
  /*
  Local Spin
  */
  nsite = L * W;
  locspinflag = (int *)malloc(sizeof(int) * nsite);
  for (isite = 0; isite < nsite; isite++)locspinflag[isite] = 0;
  /*
  Transfer
  */
  ntrans = 2 * (5 * L*W + 6*L*(W - 1));
  transindx = (int **)malloc(sizeof(int*) * ntrans);
  trans = (double *)malloc(sizeof(double) * ntrans);
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  ktrans = 0;
  for (iL = 0; iL < L; iL++) {
    for (iW = 0; iW < W; iW++){
      isite = iW + iL * W;
      for (ispin = 0; ispin < 2; ispin++){

        StdFace_trans(&ktrans, mu, isite, ispin, isite, ispin);

        jsite = iW + ((iL + 1 + 2 * L) % L) * W;
        StdFace_trans(&ktrans, t1, isite, ispin, jsite, ispin);
        StdFace_trans(&ktrans, t1, jsite, ispin, isite, ispin);

        jsite = iW + ((iL + 2 + 2 * L) % L) * W;
        StdFace_trans(&ktrans, t1p, isite, ispin, jsite, ispin);
        StdFace_trans(&ktrans, t1p, jsite, ispin, isite, ispin);

        if (iW < W - 1){
          jsite = (iW + 1) + iL * W;
          StdFace_trans(&ktrans, t0, isite, ispin, jsite, ispin);
          StdFace_trans(&ktrans, t0, jsite, ispin, isite, ispin);

          jsite = (iW + 1) + ((iL + 1 + 2 * L) % L) * W;
          StdFace_trans(&ktrans, t2, isite, ispin, jsite, ispin);
          StdFace_trans(&ktrans, t2, jsite, ispin, isite, ispin);

          jsite = (iW + 1) + ((iL - 1 + 2 * L) % L) * W;
          StdFace_trans(&ktrans, t2p, isite, ispin, jsite, ispin);
          StdFace_trans(&ktrans, t2p, jsite, ispin, isite, ispin);
        }

      }
    }
  }
  /*
  Interaction
  */
  nintr = L*W + 4 * (L*(W - 1) + L*W + L*W + 2 * L*(W - 1));
  intrindx = (int **)malloc(sizeof(int*) * nintr);
  intr = (double *)malloc(sizeof(double) * nintr);
  for (kintr = 0; kintr < nintr; kintr++){
    intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  kintr = 0;
  for (iW = 0; iW < W; iW++){
    for (iL = 0; iL < L; iL++){

      isite = iW + iL * W;

      StdFace_intr(&kintr, U, isite, 0, isite, 0, isite, 1, isite, 1);

      jsite = iW + ((iL + 1 + 2 * L) % L) * W;
      StdFace_Coulomb(&kintr, V1, isite, jsite);

      jsite = iW + ((iL + 2 + 2 * L) % L) * W;
      StdFace_Coulomb(&kintr, V1p, isite, jsite);

      if (iW < W - 1){
        jsite = (iW + 1) + iL * W;
        StdFace_Coulomb(&kintr, V0, isite, jsite);

        jsite = (iW + 1) + ((iL + 1 + 2 * L) % L) * W;
        StdFace_Coulomb(&kintr, V2, isite, jsite);

        jsite = (iW + 1) + ((iL - 1 + 2 * L) % L) * W;
        StdFace_Coulomb(&kintr, V2p, isite, jsite);
      }
    }
  }
  /*
  Set mTPQ parameter
  */
  if (lGC == 0){
    LargeValue0 =  fabs(mu) * (double)nelec / (double)(L * W)
      + 2.0 * 2.0 * fabs(t1) + 2.0 * 2.0 * fabs(t1p) 
      + 2.0 * 2.0 * fabs(t0) + 2.0 * 2.0 * fabs(t2) + 2.0 * 2.0 * fabs(t2p)
      + fabs(U)
      + 4.0 * fabs(V1) + 4.0 * fabs(V1p) + 4.0 * fabs(V0)
      + 2.0 * 2.0 * fabs(V2) + 2.0 * 2.0 * fabs(V2p);
  }
  else{
    LargeValue0 = fabs(mu) * 2.0 
      + 2.0 * 2.0 * fabs(t1) + 2.0 * 2.0 * fabs(t1p)
      + 2.0 * 2.0 * fabs(t0) + 2.0 * 2.0 * fabs(t2) + 2.0 * 2.0 * fabs(t2p)
      + fabs(U) 
      + 4.0 * fabs(V1) + 4.0 * fabs(V1p) + 4.0 * fabs(V0)
      + 2.0 * 2.0 * fabs(V2) + 2.0 * 2.0 * fabs(V2p);
  }
  StdFace_PrintVal_d("LargeValue", &LargeValue, LargeValue0);
}

/**
 *
 * Setup a Hamiltonian for the generalized Heisenberg model on a square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void Spin_Ladder(
  int Sz2 /**< [in] 2 * Total Sz */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int iL, iW;
  int ktrans, kintr;
  double LargeValue0, S;

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
  StdFace_PrintVal_d("J1", &J1, 1.0);
  StdFace_PrintVal_d("J1'", &J1p, 1.0);
  StdFace_PrintVal_d("J0", &J0, 1.0);
  StdFace_PrintVal_d("J2", &J2, 1.0);
  StdFace_PrintVal_d("J2'", &J2p, 1.0);
  /**/
  StdFace_NotUsed_d("J", J);
  StdFace_NotUsed_d("J'", Jp);
  StdFace_NotUsed_d("J''", Jpp);
  StdFace_NotUsed_d("Jxy", Jxy);
  StdFace_NotUsed_d("Jx", Jx);
  StdFace_NotUsed_d("Jy", Jy);
  StdFace_NotUsed_d("Jz", Jz);
  StdFace_NotUsed_d("Jz0", Jz0);
  StdFace_NotUsed_d("Jz1", Jz1);
  StdFace_NotUsed_d("Jxy0", Jxy0);
  StdFace_NotUsed_d("Jxy1", Jxy1);
  StdFace_NotUsed_d("K", K);
  /*
  Local Spin
  */
  nsite = L * W;
  locspinflag = (int *)malloc(sizeof(int) * nsite);
  for (isite = 0; isite < nsite; isite++)locspinflag[isite] = S2;
  /*
  Transfer
  */
  ntrans = L * W * (S2 + 1 + 2 * S2);
  transindx = (int **)malloc(sizeof(int*) * ntrans);
  trans = (double *)malloc(sizeof(double) * ntrans);
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  ktrans = 0;
  for (isite = 0; isite < nsite; isite++){
    StdFace_MagLong(&ktrans, -h, isite);
    StdFace_MagTrans(&ktrans, -Gamma, isite);
  }
  /*
  Interaction
  */
  nintr = (3*L*W + 3*L*(W - 1)) * (S2 + 1) * (S2 + 1) + (2*L*W + 3*L*(W - 1))*S2*S2*2;
  intrindx = (int **)malloc(sizeof(int*) * nintr);
  intr = (double *)malloc(sizeof(double) * nintr);
  for (kintr = 0; kintr < nintr; kintr++){
    intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  kintr = 0;

  for (iL = 0; iL < L; iL++) {
    for (iW = 0; iW < W; iW++){

      isite = iW + iL * W;

      StdFace_SzSz(&kintr, S2, S2, D, isite, isite);
      
      jsite = iW + ((iL + 1 + 2 * L) % L) * W;
      StdFace_SzSz(&kintr, S2, S2, J1, isite, jsite);
      StdFace_exchange(&kintr, S2, S2, J1, isite, jsite);
      
      jsite = iW + ((iL + 2 + 2 * L) % L) * W;
      StdFace_SzSz(&kintr, S2, S2, J1p, isite, jsite);
      StdFace_exchange(&kintr, S2, S2, J1p, isite, jsite);
      
      
      if (iW < W - 1){
        jsite = (iW + 1) + iL * W;
        StdFace_SzSz(&kintr, S2, S2, J0, isite, jsite);
        StdFace_exchange(&kintr, S2, S2, J0, isite, jsite);

        jsite = (iW + 1) + ((iL + 1 + 2 * L) % L) * W;
        StdFace_SzSz(&kintr, S2, S2, J2, isite, jsite);
        StdFace_exchange(&kintr, S2, S2, J2, isite, jsite);

        jsite = (iW + 1) + ((iL - 1 + 2 * L) % L) * W;
        StdFace_SzSz(&kintr, S2, S2, J2p, isite, jsite);
        StdFace_exchange(&kintr, S2, S2, J2p, isite, jsite);
      }
    }
  }
  /*
  Set mTPQ parameter
  */
  S = (double)S2 * 0.5;
  if (lGC == 0){
    LargeValue0 = (double)Sz2 / (double)(2 * nsite) * fabs(h) + S * fabs(D) + 0.5 * S * fabs(Gamma)
    + S*S*fabs(J1) + S*S*fabs(J1p) + S*S*fabs(J0) + S*S*fabs(J2) + S*S*fabs(J2p);
  }
  else{
    LargeValue0 = S * fabs(h) + S * fabs(D) + 0.5 * S * fabs(Gamma)
      + S*S*fabs(J1) + S*S*fabs(J1p) + S*S*fabs(J0) + S*S*fabs(J2) + S*S*fabs(J2p);
   }
  StdFace_PrintVal_d("LargeValue", &LargeValue, LargeValue0);
}

/**
*
* Setup a Hamiltonian for the generalized Heisenberg model on a square lattice
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void Spin_Ladder_Boost(
  int Sz2 /**< [in] 2 * Total Sz */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int iL, iW;
  int ktrans, kintr;
  double LargeValue0, S;
  FILE *fp;
  int j;

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
  StdFace_PrintVal_d("J1", &J1, 1.0);
  StdFace_PrintVal_d("J1'", &J1p, 1.0);
  StdFace_PrintVal_d("J0", &J0, 1.0);
  StdFace_PrintVal_d("J2", &J2, 1.0);
  StdFace_PrintVal_d("J2'", &J2p, 1.0);
  /**/
  StdFace_NotUsed_d("J", J);
  StdFace_NotUsed_d("J'", Jp);
  StdFace_NotUsed_d("J''", Jpp);
  StdFace_NotUsed_d("Jxy", Jxy);
  StdFace_NotUsed_d("Jx", Jx);
  StdFace_NotUsed_d("Jy", Jy);
  StdFace_NotUsed_d("Jz", Jz);
  StdFace_NotUsed_d("Jz0", Jz0);
  StdFace_NotUsed_d("Jz1", Jz1);
  StdFace_NotUsed_d("Jxy0", Jxy0);
  StdFace_NotUsed_d("Jxy1", Jxy1);
  StdFace_NotUsed_d("K", K);
  /*
  Local Spin
  */
  nsite = L * W;
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
    LargeValue0 = (double)Sz2 / (double)(2 * nsite) * fabs(h) + S * fabs(D) + 0.5 * S * fabs(Gamma)
      + S*S*fabs(J1) + S*S*fabs(J1p) + S*S*fabs(J0) + S*S*fabs(J2) + S*S*fabs(J2p);
  }
  else {
    LargeValue0 = S * fabs(h) + S * fabs(D) + 0.5 * S * fabs(Gamma)
      + S*S*fabs(J1) + S*S*fabs(J1p) + S*S*fabs(J0) + S*S*fabs(J2) + S*S*fabs(J2p);
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
  fprintf(fp, "%d  # Number of type of J\n", 5);
  fprintf(fp, "# J 1 (inter chain, vertical)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * J0, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * J0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * J0, 0.0);
  fprintf(fp, "# J 2 (Nearest neighbor, along chain)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * J1, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * J1, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * J1, 0.0);
  fprintf(fp, "# J 3 (Second nearest neighbor, along chain)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * J1p, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * J1p, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * J1p, 0.0);
  fprintf(fp, "# J 4 (inter chain, diagonal1)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * J2, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * J2, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * J2, 0.0);
  fprintf(fp, "# J 5 (inter chain, diagonal2)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * J2p, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * J2p, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * J2p, 0.0);
  /*
  Topology
  */
  if (S2 != 1) {
    fprintf(stderr, "\n ERROR! S2 must be 1 in Boost. \n\n");
    exitMPI(-1);
  }
  ishift_nspin = 2;
  if (W != 2) {
    fprintf(stderr, "\n ERROR! W != 2 \n\n");
    exitMPI(-1);
  }
  if (L < 3) {
    fprintf(stderr, "\n ERROR! L < 3 \n\n");
    exitMPI(-1);
  }
  num_pivot = 1;
  /**/
  fprintf(fp, "# W0  R0  num_pivot  ishift_nspin\n");
  fprintf(fp, "%d %d %d %d\n", W, L, num_pivot, ishift_nspin);

  list_6spin_star = (int **)malloc(sizeof(int*) * num_pivot);
  for (j = 0; j < num_pivot; j++) {
    list_6spin_star[j] = (int *)malloc(sizeof(int) * 7);
  }

  list_6spin_star[0][0] = 7; // num of J
  list_6spin_star[0][1] = 1;
  list_6spin_star[0][2] = 1;
  list_6spin_star[0][3] = 1;
  list_6spin_star[0][4] = 1;
  list_6spin_star[0][5] = 1;
  list_6spin_star[0][6] = 1; // flag

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

  list_6spin_pair[0][0][0] = 0;
  list_6spin_pair[0][1][0] = 1;
  list_6spin_pair[0][2][0] = 2;
  list_6spin_pair[0][3][0] = 3;
  list_6spin_pair[0][4][0] = 4;
  list_6spin_pair[0][5][0] = 5;
  list_6spin_pair[0][6][0] = 1; // type of J
  list_6spin_pair[0][0][1] = 0;
  list_6spin_pair[0][1][1] = 2;
  list_6spin_pair[0][2][1] = 1;
  list_6spin_pair[0][3][1] = 3;
  list_6spin_pair[0][4][1] = 4;
  list_6spin_pair[0][5][1] = 5;
  list_6spin_pair[0][6][1] = 2; // type of J
  list_6spin_pair[0][0][2] = 1;
  list_6spin_pair[0][1][2] = 3;
  list_6spin_pair[0][2][2] = 0;
  list_6spin_pair[0][3][2] = 2;
  list_6spin_pair[0][4][2] = 4;
  list_6spin_pair[0][5][2] = 5;
  list_6spin_pair[0][6][2] = 2; // type of J
  list_6spin_pair[0][0][3] = 0;
  list_6spin_pair[0][1][3] = 4;
  list_6spin_pair[0][2][3] = 1;
  list_6spin_pair[0][3][3] = 2;
  list_6spin_pair[0][4][3] = 3;
  list_6spin_pair[0][5][3] = 5;
  list_6spin_pair[0][6][3] = 3; // type of J
  list_6spin_pair[0][0][4] = 1;
  list_6spin_pair[0][1][4] = 5;
  list_6spin_pair[0][2][4] = 0;
  list_6spin_pair[0][3][4] = 2;
  list_6spin_pair[0][4][4] = 3;
  list_6spin_pair[0][5][4] = 4;
  list_6spin_pair[0][6][4] = 3; // type of J
  list_6spin_pair[0][0][5] = 0;
  list_6spin_pair[0][1][5] = 3;
  list_6spin_pair[0][2][5] = 1;
  list_6spin_pair[0][3][5] = 2;
  list_6spin_pair[0][4][5] = 4;
  list_6spin_pair[0][5][5] = 5;
  list_6spin_pair[0][6][5] = 4; // type of J
  list_6spin_pair[0][0][6] = 1;
  list_6spin_pair[0][1][6] = 2;
  list_6spin_pair[0][2][6] = 0;
  list_6spin_pair[0][3][6] = 3;
  list_6spin_pair[0][4][6] = 4;
  list_6spin_pair[0][5][6] = 5;
  list_6spin_pair[0][6][6] = 5; // type of J

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

/**
 *
 * Setup a Hamiltonian for the Kondo lattice model on a square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void KondoLattice_Ladder(
  int nelec /**< [in] The number of valence electrons */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int ispin;
  int iL, iW;
  int ktrans, kintr;
  double LargeValue0, S;
  /**/
  fprintf(stdoutMPI, "\n");
  fprintf(stdoutMPI, "#######  Parameter Summary  #######\n");
  fprintf(stdoutMPI, "\n");
  StdFace_RequiredVal_i("L", L);
  StdFace_RequiredVal_i("W", W);
  StdFace_PrintVal_d("a", &a, 1.0);
  /**/
  StdFace_PrintVal_i("2S", &S2, 1);
  StdFace_PrintVal_d("mu", &mu, 0.0);
  StdFace_PrintVal_d("t0", &t0, 1.0);
  StdFace_PrintVal_d("t1", &t1, 1.0);
  StdFace_PrintVal_d("t1'", &t1p, 1.0);
  StdFace_PrintVal_d("t2", &t2, 1.0);
  StdFace_PrintVal_d("t2'", &t2p, 1.0);
  StdFace_PrintVal_d("J", &J, 0.0);
  /**/
  StdFace_NotUsed_d("U", U);
  StdFace_NotUsed_d("t", t);
  StdFace_NotUsed_d("t'", tp);
  StdFace_NotUsed_d("tpp", tpp);
  StdFace_NotUsed_d("V'", Vp);
  StdFace_NotUsed_d("Vpp", Vpp);
  StdFace_NotUsed_d("V0", V0);
  StdFace_NotUsed_d("V1", V1);
  StdFace_NotUsed_d("V2", V2);
  /*
  Local Spin
  */
  nsite = 2 * L * W;
  locspinflag = (int *)malloc(sizeof(int) * nsite);
  for (iL = 0; iL < L * W; iL++){
    locspinflag[iL] = S2;
    locspinflag[iL + L * W] = 0;
  }
  /*
  Transfer
  */
  ntrans = 2 * (5 * L*W + 6 * L*(W - 1));
  transindx = (int **)malloc(sizeof(int*) * ntrans);
  trans = (double *)malloc(sizeof(double) * ntrans);
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  ktrans = 0;
  for (iL = 0; iL < L; iL++) {
    for (iW = 0; iW < W; iW++){
      isite = 2 * L * W + iW + iL * W;
      for (ispin = 0; ispin < 2; ispin++){

        StdFace_trans(&ktrans, mu, isite, ispin, isite, ispin);

        jsite = L * W + iW + ((iL + 1 + 2 * L) % L) * W;
        StdFace_trans(&ktrans, t1, isite, ispin, jsite, ispin);
        StdFace_trans(&ktrans, t1, jsite, ispin, isite, ispin);

        jsite = L * W + iW + ((iL + 2 + 2 * L) % L) * W;
        StdFace_trans(&ktrans, t1p, isite, ispin, jsite, ispin);
        StdFace_trans(&ktrans, t1p, jsite, ispin, isite, ispin);

        if (iW < W - 1){
          jsite = L * W + (iW + 1) + iL * W;
          StdFace_trans(&ktrans, t0, isite, ispin, jsite, ispin);
          StdFace_trans(&ktrans, t0, jsite, ispin, isite, ispin);

          jsite = L * W + (iW + 1) + ((iL + 1 + 2 * L) % L) * W;
          StdFace_trans(&ktrans, t2, isite, ispin, jsite, ispin);
          StdFace_trans(&ktrans, t2, jsite, ispin, isite, ispin);

          jsite = L * W + (iW + 1) + ((iL - 1 + 2 * L) % L) * W;
          StdFace_trans(&ktrans, t2p, isite, ispin, jsite, ispin);
          StdFace_trans(&ktrans, t2p, jsite, ispin, isite, ispin);
        }

      }
    }
  }
  /*
  Interaction
  */
  nintr = L * W * ((S2 + 1) * (1 + 1) + 2 * S2 * 1);
  intrindx = (int **)malloc(sizeof(int*) * nintr);
  intr = (double *)malloc(sizeof(double) * nintr);
  for (kintr = 0; kintr < nintr; kintr++){
    intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  kintr = 0;
    for (iL = 0; iL < L; iL++){
      for (iW = 0; iW < W; iW++) {

      isite = iL + iW * L + L * W;
      jsite = iL + iW * L;

      StdFace_exchange(&kintr, 1, S2, J, isite, jsite);
      StdFace_SzSz(&kintr, 1, S2, J, isite, jsite);
    }
  }
  /*
  Set mTPQ parameter
  */
  S = (double)S2 * 0.5;
  if (lGC == 0){
    LargeValue0 = fabs(mu) * (double)nelec / (double)(L * W)
      + 2.0 * 2.0 * fabs(t1) + 2.0 * 2.0 * fabs(t1p)
      + 2.0 * 2.0 * fabs(t0) + 2.0 * 2.0 * fabs(t2) + 2.0 * 2.0 * fabs(t2p)
      + 0.5 * S * fabs(J);
  }
  else{
    LargeValue0 = fabs(mu) * 2.0 
      + 2.0 * 2.0 * fabs(t1) + 2.0 * 2.0 * fabs(t1p)
      + 2.0 * 2.0 * fabs(t0) + 2.0 * 2.0 * fabs(t2) + 2.0 * 2.0 * fabs(t2p)
      + 0.5 * S * fabs(J);
  }
  StdFace_PrintVal_d("LargeValue", &LargeValue, LargeValue0);
}
