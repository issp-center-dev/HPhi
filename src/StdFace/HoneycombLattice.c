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
 * Setup a Hamiltonian for the Hubbard model on a Honeycomb lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void FermionHubbard_HoneycombLattice(
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
  StdFace_PrintVal_d("t", &t, 1.0);
  StdFace_PrintVal_d("t0", &t0, t);
  StdFace_PrintVal_d("t1", &t1, t);
  StdFace_PrintVal_d("t2", &t2, t);
  StdFace_PrintVal_d("V", &V, 0.0);
  StdFace_PrintVal_d("V0", &V0, V);
  StdFace_PrintVal_d("V1", &V1, V);
  StdFace_PrintVal_d("V2", &V2, V);
  /**/
  StdFace_NotUsed_i("2S", S2);
  StdFace_NotUsed_d("t'", tp);
  StdFace_NotUsed_d("tpp", tpp);
  StdFace_NotUsed_d("V'", Vp);
  StdFace_NotUsed_d("Vpp", Vpp);
  /*
  Local Spin
  */
  nsite = L * W * 2;
  locspinflag = (int *)malloc(sizeof(int) * nsite);
  for (isite = 0; isite < nsite; isite++)locspinflag[isite] = 0;
  /*
  Transfer
  */
  ntrans = L * W * 2 * 2 * 4;
  transindx = (int **)malloc(sizeof(int*) * ntrans);
  trans = (double *)malloc(sizeof(double) * ntrans);
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  ktrans = 0;
  for (iW = 0; iW < W; iW++){
    for (iL = 0; iL < L; iL++){
      for (ispin = 0; ispin < 2; ispin++){
        isite = (iL + iW * L) * 2;
        StdFace_trans(&ktrans, mu, isite, ispin, isite, ispin);
        jsite = (iL + iW * L) * 2 + 1;
        StdFace_trans(&ktrans, t0, isite, ispin, jsite, ispin);
        jsite = (((iL - 1 + 2 * L) % L) + ((iW + 0 + 2 * W) % W) * L) * 2 + 1;
        StdFace_trans(&ktrans, t1, isite, ispin, jsite, ispin);
        jsite = (((iL + 0 + 2 * L) % L) + ((iW - 1 + 2 * W) % W) * L) * 2 + 1;
        StdFace_trans(&ktrans, t2, isite, ispin, jsite, ispin);

        isite = (iL + iW * L) * 2 + 1;
        StdFace_trans(&ktrans, mu, isite, ispin, isite, ispin);
        jsite = (iL + iW * L) * 2;
        StdFace_trans(&ktrans, t0, isite, ispin, jsite, ispin);
        jsite = (((iL + 1 + 2 * L) % L) + ((iW + 0 + 2 * W) % W) * L) * 2;
        StdFace_trans(&ktrans, t1, isite, ispin, jsite, ispin);
        jsite = (((iL + 0 + 2 * L) % L) + ((iW + 1 + 2 * W) % W) * L) * 2;
        StdFace_trans(&ktrans, t2, isite, ispin, jsite, ispin);
      }
    }
  }
  /*
  Interaction
  */
  nintr = L * W * (2 + 4 * 3);
  intrindx = (int **)malloc(sizeof(int*) * nintr);
  intr = (double *)malloc(sizeof(double) * nintr);
  for (kintr = 0; kintr < nintr; kintr++){
    intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  kintr = 0;
  for (iW = 0; iW < W; iW++){
    for (iL = 0; iL < L; iL++){
      isite = (iL + iW * L) * 2;
      StdFace_intr(&kintr, U, isite, 0, isite, 0, isite, 1, isite, 1);
      isite = (iL + iW * L) * 2 + 1;
      StdFace_intr(&kintr, U, isite, 0, isite, 0, isite, 1, isite, 1);

      jsite = (iL + iW * L) * 2;
      StdFace_Coulomb(&kintr, V0, isite, jsite);
      jsite = (((iL + 1 + 2 * L) % L) + ((iW + 0 + 2 * W) % W) * L) * 2;
      StdFace_Coulomb(&kintr, V1, isite, jsite);
      jsite = (((iL + 0 + 2 * L) % L) + ((iW + 1 + 2 * W) % W) * L) * 2;
      StdFace_Coulomb(&kintr, V2, isite, jsite);
    }
  }
  /*
  Set mTPQ parameter
  */
  if (lGC == 0){
    LargeValue0 =  fabs(mu) * (double)nelec / (double)(2 * L * W)
      + 2.0 * 3.0 * fabs(t)
      + fabs(U) + 2.0 * 3.0 * fabs(V);
  }
  else{
    LargeValue0 = fabs(mu) * 2.0 + 2.0 * 3.0 * fabs(t)
      + fabs(U) + 2.0 * 3.0 * fabs(V);
  }
  StdFace_PrintVal_i("LargeValue", &LargeValue, (int)LargeValue0 + 1);
}

/**
 *
 * Setup a Hamiltonian for the generalized Heisenberg model on a Heisenberg lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void Spin_HoneycombLattice(
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
  StdFace_PrintVal_d("J", &J, 0.0);
  StdFace_PrintVal_d("J0", &J0, J);
  StdFace_PrintVal_d("J1", &J1, J);
  StdFace_PrintVal_d("J2", &J2, J);
  StdFace_PrintVal_d("Jz", &Jz, J);
  StdFace_PrintVal_d("Jz0", &Jz0, Jz);
  StdFace_PrintVal_d("Jz1", &Jz1, Jz);
  StdFace_PrintVal_d("Jz2", &Jz2, Jz);
  StdFace_PrintVal_d("Jxy", &Jxy, J);
  StdFace_PrintVal_d("Jxy0", &Jxy0, Jxy);
  StdFace_PrintVal_d("Jxy1", &Jxy1, Jxy);
  StdFace_PrintVal_d("Jxy2", &Jxy2, Jxy);
  StdFace_PrintVal_d("Jx", &Jx, Jxy);
  StdFace_PrintVal_d("Jx0", &Jx0, Jx);
  StdFace_PrintVal_d("Jx1", &Jx1, Jx);
  StdFace_PrintVal_d("Jx2", &Jx2, Jx);
  StdFace_PrintVal_d("Jy", &Jy, Jxy);
  StdFace_PrintVal_d("Jy0", &Jy0, Jy);
  StdFace_PrintVal_d("Jy1", &Jy1, Jy);
  StdFace_PrintVal_d("Jy2", &Jy2, Jy);
  Jxy = 0.5 * (Jx + Jy);
  Jxy0 = 0.5 * (Jx0 + Jy0);
  Jxy1 = 0.5 * (Jx1 + Jy1);
  Jxy2 = 0.5 * (Jx2 + Jy2);
  /**/
  StdFace_NotUsed_d("J'", Jp);
  StdFace_NotUsed_d("Jz'", Jzp);
  StdFace_NotUsed_d("Jxy'", Jxyp);
  StdFace_NotUsed_d("Jx'", Jxp);
  StdFace_NotUsed_d("Jy'", Jyp);
  StdFace_NotUsed_d("K", K);
  /*
  Local Spin
  */
  nsite = L * W * 2;
  locspinflag = (int *)malloc(sizeof(int) * nsite);
  for (isite = 0; isite < nsite; isite++)locspinflag[isite] = S2;
  /*
  Transfer
  */
  ntrans = L * W * 2 * (S2 + 1 + 2 * S2);
  transindx = (int **)malloc(sizeof(int*) * ntrans);
  trans = (double *)malloc(sizeof(double) * ntrans);
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  ktrans = 0;
  for (isite = 0; isite < nsite; isite++){
    StdFace_MagLong(&ktrans, h, isite);
    StdFace_MagTrans(&ktrans, Gamma, isite);
  }
  /*
  Interaction
  */
  nintr = L * W * ((S2 + 1) * (S2 + 1) * 5 + 2 * S2 * S2 * 6);
  intrindx = (int **)malloc(sizeof(int*) * nintr);
  intr = (double *)malloc(sizeof(double) * nintr);
  for (kintr = 0; kintr < nintr; kintr++){
    intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  kintr = 0;
  for (iW = 0; iW < W; iW++){
    for (iL = 0; iL < L; iL++){

      isite = (iL + iW * L) * 2;
      StdFace_SzSz(&kintr, S2, S2, D, isite, isite);
      isite = (iL + iW * L) * 2 + 1;
      StdFace_SzSz(&kintr, S2, S2, D, isite, isite);

      jsite = (iL + iW * L) * 2;
      StdFace_SzSz(&kintr, S2, S2, Jz0, isite, jsite);
      StdFace_exchange(&kintr, S2, S2, Jxy0, isite, jsite);
      StdFace_Kitaev(&kintr, S2, S2, 0.5 * (Jx0 - Jy0), isite, jsite);
      jsite = (((iL + 1 + 2 * L) % L) + ((iW + 0 + 2 * W) % W) * L) * 2;
      StdFace_SzSz(&kintr, S2, S2, Jz1, isite, jsite);
      StdFace_exchange(&kintr, S2, S2, Jxy1, isite, jsite);
      StdFace_Kitaev(&kintr, S2, S2, 0.5 * (Jx1 - Jy1), isite, jsite);
      jsite = (((iL + 0 + 2 * L) % L) + ((iW + 1 + 2 * W) % W) * L) * 2;
      StdFace_SzSz(&kintr, S2, S2, Jz2, isite, jsite);
      StdFace_exchange(&kintr, S2, S2, Jxy2, isite, jsite);
      StdFace_Kitaev(&kintr, S2, S2, 0.5 * (Jx2 - Jy2), isite, jsite);
    }
  }
  /*
  Set mTPQ parameter
  */
  S = (double)S2 * 0.5;
  if (lGC == 0){
    LargeValue0 = (double)Sz2 / (double)(2 * nsite) * fabs(h) + S * fabs(D) + S * S * fabs(Gamma)
      + 1.0 / 2.0 * S * S * (fabs(Jx0) + fabs(Jy0) + fabs(Jz0))
      + 1.0 / 2.0 * S * S * (fabs(Jx1) + fabs(Jy1) + fabs(Jz1))
      + 1.0 / 2.0 * S * S * (fabs(Jx2) + fabs(Jy2) + fabs(Jz2));
  }
  else{
    LargeValue0 = S * fabs(h) + S * fabs(D) + S * S * fabs(Gamma)
      + 1.0 / 2.0 * S * S * (fabs(Jx0) + fabs(Jy0) + fabs(Jz0))
      + 1.0 / 2.0 * S * S * (fabs(Jx1) + fabs(Jy1) + fabs(Jz1))
      + 1.0 / 2.0 * S * S * (fabs(Jx2) + fabs(Jy2) + fabs(Jz2));
  }
  StdFace_PrintVal_i("LargeValue", &LargeValue, (int)LargeValue0 + 1);
}

/**
 *
 * Setup a Hamiltonian for the Kondo lattice model on a Honeycomb lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void KondoLattice_HoneycombLattice(
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
  nsite = 4 * L * W;
  locspinflag = (int *)malloc(sizeof(int) * nsite);
  for (iL = 0; iL < L * W; iL++){
    locspinflag[4 * iL] = 0;
    locspinflag[4 * iL + 1] = S2;
    locspinflag[4 * iL + 2] = 0;
    locspinflag[4 * iL + 3] = S2;
  }
  /*
  Transfer
  */
  ntrans = L * W * 2 * 2 * 4;
  transindx = (int **)malloc(sizeof(int*) * ntrans);
  trans = (double *)malloc(sizeof(double) * ntrans);
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  ktrans = 0;
  for (iW = 0; iW < W; iW++){
    for (iL = 0; iL < L; iL++){
      for (ispin = 0; ispin < 2; ispin++){
        isite = (iL + iW * L) * 4;
        StdFace_trans(&ktrans, mu, isite, ispin, isite, ispin);
        jsite = (iL + iW * L) * 4 + 2;
        StdFace_trans(&ktrans, t0, isite, ispin, jsite, ispin);
        jsite = (((iL - 1 + 2 * L) % L) + ((iW + 0 + 2 * W) % W) * L) * 4 + 2;
        StdFace_trans(&ktrans, t1, isite, ispin, jsite, ispin);
        jsite = (((iL + 0 + 2 * L) % L) + ((iW - 1 + 2 * W) % W) * L) * 4 + 2;
        StdFace_trans(&ktrans, t2, isite, ispin, jsite, ispin);

        isite = (iL + iW * L) * 4 + 2;
        StdFace_trans(&ktrans, mu, isite, ispin, isite, ispin);
        jsite = (iL + iW * L) * 4;
        StdFace_trans(&ktrans, t0, isite, ispin, jsite, ispin);
        jsite = (((iL + 1 + 2 * L) % L) + ((iW + 0 + 2 * W) % W) * L) * 4;
        StdFace_trans(&ktrans, t1, isite, ispin, jsite, ispin);
        jsite = (((iL + 0 + 2 * L) % L) + ((iW + 1 + 2 * W) % W) * L) * 4;
        StdFace_trans(&ktrans, t2, isite, ispin, jsite, ispin);
      }
    }
  }
  /*
  Interaction
  */
  nintr = L * W * 2 * ((S2 + 1) * (1 + 1) + 2 * S2 * 1);
  intrindx = (int **)malloc(sizeof(int*) * nintr);
  intr = (double *)malloc(sizeof(double) * nintr);
  for (kintr = 0; kintr < nintr; kintr++){
    intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  kintr = 0;
  for (iW = 0; iW < W; iW++){
    for (iL = 0; iL < L; iL++){
      isite = 4 * (iL + iW * L);
      jsite = 4 * (iL + iW * L) + 1;
      StdFace_exchange(&kintr, 1, S2, J, isite, jsite);
      StdFace_SzSz(&kintr, 1, S2, J, isite, jsite);

      isite = 4 * (iL + iW * L) + 2;
      jsite = 4 * (iL + iW * L) + 3;
      StdFace_exchange(&kintr, 1, S2, J, isite, jsite);
      StdFace_SzSz(&kintr, 1, S2, J, isite, jsite);
    }
  }
  /*
  Set mTPQ parameter
  */
  S = (double)S2 * 0.5;
  if (lGC == 0){
    LargeValue0 = fabs(mu) * (double)nelec / (double)(2 * L * W) + 2.0 * 3.0 * fabs(t) + 0.5 * S * fabs(J);
  }
  else{
    LargeValue0 = fabs(mu) * 2.0 + 2.0 * 3.0 * fabs(t) + 0.5 * S * fabs(J);
  }
  StdFace_PrintVal_i("LargeValue", &LargeValue, (int)LargeValue0 + 1);
}
