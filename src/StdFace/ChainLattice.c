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
#include "../include/wrapperMPI.h"

/**
 *
 * Setup a Hamiltonian for the Hubbard model on a Chain lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void FermionHubbard_ChainLattice(
  int nelec /**< [in] The number of electrons */, 
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int ispin;
  int iL, iL2;
  int ktrans, kintr;
  double trans0, LargeValue0;
  fprintf(stdoutMPI, "\n");
  fprintf(stdoutMPI, "#######  Parameter Summary  #######\n");
  fprintf(stdoutMPI, "\n");
  StdFace_RequiredVal_i("L", L);
  StdFace_NotUsed_i("W", W);
  StdFace_PrintVal_d("a", &a, 1.0);
  /**/
  StdFace_PrintVal_d("mu", &mu, 0.0);
  StdFace_PrintVal_d("U", &U, 0.0);
  StdFace_PrintVal_d("t", &t, 1.0);
  StdFace_PrintVal_d("V", &V, 0.0);
  StdFace_PrintVal_d("t'", &tp, 0.0);
  StdFace_PrintVal_d("V'", &Vp, 0.0);
  /**/
  StdFace_NotUsed_i("2S",S2);
  StdFace_NotUsed_d("tpp", tpp);
  StdFace_NotUsed_d("t0", t0);
  StdFace_NotUsed_d("t1", t1);
  StdFace_NotUsed_d("t2", t2);
  StdFace_NotUsed_d("Vpp", Vpp);
  StdFace_NotUsed_d("V0", V0);
  StdFace_NotUsed_d("V1", V1);
  StdFace_NotUsed_d("V2", V2);
  /*
  Local Spin
  */
  nsite = L;
  locspinflag = (int *)malloc(sizeof(int) * nsite);
  for (isite = 0; isite < nsite; isite++)locspinflag[isite] = 0;
  /*
  Transfer
  */
  ntrans = nsite * 2 * 5;
  transindx = (int **)malloc(sizeof(int*) * ntrans);
  trans = (double *)malloc(sizeof(double) * ntrans);
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  ktrans = 0;
  for (iL = 0; iL < L; iL++){
    isite = iL;
    for (ispin = 0; ispin < 2; ispin++){
      for (iL2 = -2; iL2 <= 2; iL2++){

        if (iL2 == 0) trans0 = mu;
        else if(abs(iL2) == 1) trans0 = t;
        else trans0 = tp;

        jsite = (iL + iL2 + 2 * L) % L;
        StdFace_trans(&ktrans, trans0, isite, ispin, jsite, ispin);
      }
    }
  }
  /*
  Interaction
  */
  nintr = nsite * (1 + 4 * 2);
  intrindx = (int **)malloc(sizeof(int*) * nintr);
  intr = (double *)malloc(sizeof(double) * nintr);
  for (kintr = 0; kintr < nintr; kintr++){
    intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  kintr = 0;
  for (iL = 0; iL < L; iL++){

    isite = iL;

    StdFace_intr(&kintr, U, isite, 0, isite, 0, isite, 1, isite, 1);

    iL2 = (iL + 1) % L;
    jsite = iL2;
    StdFace_Coulomb(&kintr, V, isite, jsite);

    iL2 = (iL + 2) % L;
    jsite = iL2;
    StdFace_Coulomb(&kintr, Vp, isite, jsite);
  }
  /*
  Set mTPQ parameter
  */
  if (lGC == 0){
    LargeValue0 = fabs(mu) * (double)nelec / (double)(L * W) 
      + 2.0 * 2.0 * fabs(t) + 2.0 * 2.0 * fabs(tp)
      + fabs(U) + 2.0 * 2.0 * fabs(V) + 2.0 * 2.0 * fabs(Vp);
  }
  else{
    LargeValue0 = fabs(mu) * 2.0 + 2.0 * 2.0 * fabs(t) + 2.0 * 2.0 * fabs(tp)
        +fabs(U) + 2.0 * 2.0 * fabs(V) + 2.0 * 2.0 * fabs(Vp);
  }
  StdFace_PrintVal_d("LargeValue", &LargeValue, LargeValue0);
}

/**
  *
  * Setup a Hamiltonian for the generalized Heisenberg model on a Chain lattice
  *
  * @author Mitsuaki Kawamura (The University of Tokyo)
  */
void Spin_ChainLattice(
  int Sz2 /**< [in] 2 * Total Sz */, 
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int iL, iL2;
  int ktrans, kintr;
  double LargeValue0, S;

  fprintf(stdoutMPI, "\n");
  fprintf(stdoutMPI, "#######  Parameter Summary  #######\n");
  fprintf(stdoutMPI, "\n");
  StdFace_RequiredVal_i("L", L);
  StdFace_NotUsed_i("W", W);
  StdFace_PrintVal_d("a", &a, 1.0);
  /**/
  StdFace_PrintVal_i("2S", &S2, 1);
  StdFace_PrintVal_d("h", &h, 0.0);
  StdFace_PrintVal_d("Gamma", &Gamma, 0.0);
  StdFace_PrintVal_d("D", &D, 0.0);
  StdFace_PrintVal_d("J", &J, 1.0);
  StdFace_PrintVal_d("Jz", &Jz, J);
  StdFace_PrintVal_d("Jxy", &Jxy, J);
  StdFace_PrintVal_d("Jx", &Jx, Jxy);
  StdFace_PrintVal_d("Jy", &Jy, Jxy);
  StdFace_PrintVal_d("J'", &Jp, 0.0);
  StdFace_PrintVal_d("Jz'", &Jzp, Jp);
  StdFace_PrintVal_d("Jxy'", &Jxyp, Jp);
  StdFace_PrintVal_d("Jx'", &Jxp, Jxyp);
  StdFace_PrintVal_d("Jy'", &Jyp, Jxyp);
  Jxy = 0.5 * (Jx + Jy);
  Jxyp = 0.5 * (Jxp + Jyp);
  /**/
  StdFace_NotUsed_d("J0", J0);
  StdFace_NotUsed_d("J1", J1);
  StdFace_NotUsed_d("J2", J2);
  StdFace_NotUsed_d("Jz0", Jz0);
  StdFace_NotUsed_d("Jz1", Jz1);
  StdFace_NotUsed_d("Jxy0", Jxy0);
  StdFace_NotUsed_d("Jxy1", Jxy1);
  StdFace_NotUsed_d("K", K);
  /*
  Local Spin
  */
  nsite = L;
  locspinflag = (int *)malloc(sizeof(int) * nsite);
  for (isite = 0; isite < nsite; isite++)locspinflag[isite] = S2;
  /*
  Transfer
  */
  ntrans = L * (S2 + 1 + 2 * S2);
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
  nintr = L * ((S2 + 1) * (S2 + 1) * 3 + 2 * S2 * S2 * 4);
  intrindx = (int **)malloc(sizeof(int*) * nintr);
  intr = (double *)malloc(sizeof(double) * nintr);
  for (kintr = 0; kintr < nintr; kintr++){
    intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  kintr = 0;
  for (iL = 0; iL < L; iL++){

    isite = iL;
    StdFace_SzSz(&kintr, S2, S2, D, isite, isite);

    iL2 = (iL + 1) % L;
    jsite = iL2;
    StdFace_SzSz(&kintr, S2, S2, Jz, isite, jsite);
    StdFace_exchange(&kintr, S2, S2, Jxy, isite, jsite);
    StdFace_Kitaev(&kintr, S2, S2, 0.5 * (Jx - Jy), isite, jsite);

    iL2 = (iL + 2) % L;
    jsite = iL2;
    StdFace_SzSz(&kintr, S2, S2, Jzp, isite, jsite);
    StdFace_exchange(&kintr, S2, S2, Jxyp, isite, jsite);
    StdFace_Kitaev(&kintr, S2, S2, 0.5 * (Jxp - Jyp), isite, jsite);
  }
  /*
  Set mTPQ parameter
  */
  S = (double)S2 * 0.5;
  if (lGC == 0){
    LargeValue0 = (double)Sz2 / (double)(2 * nsite) * fabs(h) + S * fabs(D) + S * S * fabs(Gamma)
      + 2.0 / 2.0 * S * S * (fabs(Jx) + fabs(Jy) + fabs(Jz))
      + 2.0 / 2.0 * S * S * (fabs(Jxp) + fabs(Jyp) + fabs(Jzp));
  }
  else{
    LargeValue0 = S * fabs(h) + S * fabs(D) + S * S * fabs(Gamma)
      + 2.0 / 2.0 * S * S * (fabs(Jx) + fabs(Jy) + fabs(Jz))
      + 2.0 / 2.0 * S * S * (fabs(Jxp) + fabs(Jyp) + fabs(Jzp));
  }
  StdFace_PrintVal_d("LargeValue", &LargeValue, LargeValue0);
}

/**
*
* Setup a Hamiltonian for the generalized Heisenberg model on a Chain lattice
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void Spin_ChainLattice_Boost(
  int Sz2 /**< [in] 2 * Total Sz */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int iL, iL2;
  int ktrans, kintr;
  int j, iW;
  double LargeValue0, S;
  FILE *fp;

  fprintf(stdoutMPI, "\n");
  fprintf(stdoutMPI, "#######  Parameter Summary  #######\n");
  fprintf(stdoutMPI, "\n");
  StdFace_RequiredVal_i("L", L);
  StdFace_NotUsed_i("W", W);
  StdFace_PrintVal_d("a", &a, 1.0);
  /**/
  StdFace_PrintVal_i("2S", &S2, 1);
  StdFace_PrintVal_d("h", &h, 0.0);
  StdFace_PrintVal_d("Gamma", &Gamma, 0.0);
  StdFace_PrintVal_d("D", &D, 0.0);
  StdFace_PrintVal_d("J", &J, 1.0);
  StdFace_PrintVal_d("Jz", &Jz, J);
  StdFace_PrintVal_d("Jxy", &Jxy, J);
  StdFace_PrintVal_d("Jx", &Jx, Jxy);
  StdFace_PrintVal_d("Jy", &Jy, Jxy);
  StdFace_PrintVal_d("J'", &Jp, 0.0);
  StdFace_PrintVal_d("Jz'", &Jzp, Jp);
  StdFace_PrintVal_d("Jxy'", &Jxyp, Jp);
  StdFace_PrintVal_d("Jx'", &Jxp, Jxyp);
  StdFace_PrintVal_d("Jy'", &Jyp, Jxyp);
  /**/
  StdFace_NotUsed_d("J0", J0);
  StdFace_NotUsed_d("J1", J1);
  StdFace_NotUsed_d("J2", J2);
  StdFace_NotUsed_d("Jz0", Jz0);
  StdFace_NotUsed_d("Jz1", Jz1);
  StdFace_NotUsed_d("Jxy0", Jxy0);
  StdFace_NotUsed_d("Jxy1", Jxy1);
  StdFace_NotUsed_d("K", K);
  /*
  Local Spin
  */
  nsite = L;
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
    LargeValue0 = (double)Sz2 / (double)(2 * nsite) * fabs(h) + S * fabs(D) + S * S * fabs(Gamma)
      + 2.0 / 2.0 * S * S * (fabs(Jx) + fabs(Jy) + fabs(Jz))
      + 2.0 / 2.0 * S * S * (fabs(Jxp) + fabs(Jyp) + fabs(Jzp));
  }
  else {
    LargeValue0 = S * fabs(h) + S * fabs(D) + S * S * fabs(Gamma)
      + 2.0 / 2.0 * S * S * (fabs(Jx) + fabs(Jy) + fabs(Jz))
      + 2.0 / 2.0 * S * S * (fabs(Jxp) + fabs(Jyp) + fabs(Jzp));
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
  fprintf(fp, "%d  # Number of type of J\n", 2);
  fprintf(fp, "# J 1\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * Jx, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * Jy, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * Jz, 0.0);
  fprintf(fp, "# J 2\n");
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.25 * Jxp, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.25 * Jyp, 0.0, 0.0, 0.0);
  fprintf(fp, "%25.15e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
    0.0, 0.0, 0.0, 0.0, 0.25 * Jzp, 0.0);
  /*
  Topology
  */
  if (S2 != 1) {
    fprintf(stderr, "\n ERROR! S2 must be 1 in Boost. \n\n");
    exitMPI(-1);
  }
  ishift_nspin = 4;
  W = 4;
  if(L % 4 != 0){
    fprintf(stderr, "\n ERROR! L % 4 != 0 \n\n");
    exitMPI(-1);
  }
  L = L / 4;
  if (L < 2) {
    fprintf(stderr, "\n ERROR! L < 8 \n\n");
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

  list_6spin_star[0][0] = 8; // num of J
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

  list_6spin_pair[0][0][0] = 1; 
  list_6spin_pair[0][1][0] = 2;
  list_6spin_pair[0][2][0] = 0;
  list_6spin_pair[0][3][0] = 3;
  list_6spin_pair[0][4][0] = 4;
  list_6spin_pair[0][5][0] = 5; 
  list_6spin_pair[0][6][0] = 1; // type of J
  list_6spin_pair[0][0][1] = 2;
  list_6spin_pair[0][1][1] = 3;
  list_6spin_pair[0][2][1] = 0;
  list_6spin_pair[0][3][1] = 1;
  list_6spin_pair[0][4][1] = 4;
  list_6spin_pair[0][5][1] = 5;
  list_6spin_pair[0][6][1] = 1; // type of J
  list_6spin_pair[0][0][2] = 3;
  list_6spin_pair[0][1][2] = 4;
  list_6spin_pair[0][2][2] = 0;
  list_6spin_pair[0][3][2] = 1;
  list_6spin_pair[0][4][2] = 2;
  list_6spin_pair[0][5][2] = 5;
  list_6spin_pair[0][6][2] = 1; // type of J
  list_6spin_pair[0][0][3] = 4;
  list_6spin_pair[0][1][3] = 5;
  list_6spin_pair[0][2][3] = 0;
  list_6spin_pair[0][3][3] = 1;
  list_6spin_pair[0][4][3] = 2;
  list_6spin_pair[0][5][3] = 3;
  list_6spin_pair[0][6][3] = 1; // type of J
  list_6spin_pair[0][0][4] = 0;
  list_6spin_pair[0][1][4] = 2;
  list_6spin_pair[0][2][4] = 1;
  list_6spin_pair[0][3][4] = 3;
  list_6spin_pair[0][4][4] = 4;
  list_6spin_pair[0][5][4] = 5;
  list_6spin_pair[0][6][4] = 2; // type of J
  list_6spin_pair[0][0][5] = 1;
  list_6spin_pair[0][1][5] = 3;
  list_6spin_pair[0][2][5] = 0;
  list_6spin_pair[0][3][5] = 2;
  list_6spin_pair[0][4][5] = 4;
  list_6spin_pair[0][5][5] = 5;
  list_6spin_pair[0][6][5] = 2; // type of J
  list_6spin_pair[0][0][6] = 2;
  list_6spin_pair[0][1][6] = 4;
  list_6spin_pair[0][2][6] = 0;
  list_6spin_pair[0][3][6] = 1;
  list_6spin_pair[0][4][6] = 3;
  list_6spin_pair[0][5][6] = 5;
  list_6spin_pair[0][6][6] = 2; // type of J
  list_6spin_pair[0][0][7] = 3;
  list_6spin_pair[0][1][7] = 5;
  list_6spin_pair[0][2][7] = 0;
  list_6spin_pair[0][3][7] = 1;
  list_6spin_pair[0][4][7] = 2;
  list_6spin_pair[0][5][7] = 4;
  list_6spin_pair[0][6][7] = 2; // type of J

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
 * Setup a Hamiltonian for the Kondo lattice model on a Chain lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void KondoLattice_ChainLattice(
  int nelec /**< [in] The number of valence electrons */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int ispin;
  int iL;
  int ktrans, kintr;
  double LargeValue0, S;

  fprintf(stdoutMPI, "\n");
  fprintf(stdoutMPI, "#######  Parameter Summary  #######\n");
  fprintf(stdoutMPI, "\n");
  StdFace_RequiredVal_i("L", L);
  StdFace_NotUsed_i("W", W);
  StdFace_PrintVal_d("a", &a, 1.0);
  /**/
  StdFace_PrintVal_i("2S", &S2, 1);
  StdFace_PrintVal_d("mu", &mu, 0.0);
  StdFace_PrintVal_d("t", &t, 1.0);
  StdFace_PrintVal_d("J", &J, 0.0);
  /**/
  StdFace_NotUsed_d("U", U);
  StdFace_NotUsed_d("t'", tp);
  StdFace_NotUsed_d("tpp", tpp);
  StdFace_NotUsed_d("t0", t0);
  StdFace_NotUsed_d("t1", t1);
  StdFace_NotUsed_d("t2", t2);
  StdFace_NotUsed_d("V'", Vp);
  StdFace_NotUsed_d("Vpp", Vpp);
  StdFace_NotUsed_d("V0", V0);
  StdFace_NotUsed_d("V1", V1);
  StdFace_NotUsed_d("V2", V2);
  /*
  Local Spin
  */
  nsite = L * 2;
  locspinflag = (int *)malloc(sizeof(int) * nsite);
  for (iL = 0; iL < L; iL++){
    locspinflag[iL] = S2;
    locspinflag[iL + L] = 0;
  }
  /*
  Transfer
  */
  ntrans = L * 2 * 3;
  transindx = (int **)malloc(sizeof(int*) * ntrans);
  trans = (double *)malloc(sizeof(double) * ntrans);
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  ktrans = 0;
  for (iL = 0; iL < L; iL++){
    isite = L + iL;
    for (ispin = 0; ispin < 2; ispin++){

      StdFace_trans(&ktrans, mu, isite, ispin, isite, ispin);

      jsite = L + (iL + 1) % L;
      StdFace_trans(&ktrans, t, isite, ispin, jsite, ispin);
      StdFace_trans(&ktrans, t, jsite, ispin, isite, ispin);
    }
  }
  /*
  Interaction
  */
  nintr = L * ((S2 + 1) * (1 + 1) + 2 * S2 * 1);
  intrindx = (int **)malloc(sizeof(int*) * nintr);
  intr = (double *)malloc(sizeof(double) * nintr);
  for (kintr = 0; kintr < nintr; kintr++){
    intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  kintr = 0;
  for (iL = 0; iL < L; iL++){

    isite = iL + L;
    jsite = iL;

    StdFace_exchange(&kintr, 1, S2, J, isite, jsite);
    StdFace_SzSz(&kintr, 1, S2, J, isite, jsite);
  }
  /*
  Set mTPQ parameter
  */
  S = (double)S2 * 0.5;
  if (lGC == 0){
    LargeValue0 = fabs(mu) * (double)nelec / (double)(L * W) + 2.0 * 2.0 * fabs(t) + 0.5 * S * fabs(J);
  }
  else{
    LargeValue0 = fabs(mu) * 2.0 + 2.0 * 2.0 * fabs(t) + 0.5 * S * fabs(J);
  }
  StdFace_PrintVal_d("LargeValue", &LargeValue, LargeValue0);
}
