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
 * Setup a Hamiltonian for the Hubbard model on a square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void FermionHubbard_SquareLattice(
  int nelec /**< [in] The number of electrons */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int ispin;
  int iL, iW, iL2, iW2;
  int ktrans, kintr;
  double trans0, LargeValue0;
  printf("\n");
  printf("#######  Parameter Summary  #######\n");
  printf("\n");
  StdFace_RequiredVal_i("L", L);
  StdFace_RequiredVal_i("W", W);
  StdFace_PrintVal_d("a", &a, 1.0);
  /**/
  StdFace_PrintVal_d("mu", &mu, 0.0);
  StdFace_PrintVal_d("U", &U, 0.0);
  StdFace_PrintVal_d("t", &t, 1.0);
  StdFace_PrintVal_d("V", &V, 0.0);
  StdFace_PrintVal_d("t'", &tp, 0.0);
  StdFace_PrintVal_d("V'", &Vp, 0.0);
  /**/
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
  nsite = L * W;
  locspinflag = (int *)malloc(sizeof(int) * nsite);
  for (isite = 0; isite < nsite; isite++)locspinflag[isite] = 1;
  /*
  Transfer
  */
  ntrans = L * W * 2 * 9;
  transindx = (int **)malloc(sizeof(int*) * ntrans);
  trans = (double *)malloc(sizeof(double) * ntrans);
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  ktrans = 0;
  for (iW = 0; iW < W; iW++){
    for (iL = 0; iL < L; iL++){
      isite = iL + iW * L;
      for (ispin = 0; ispin < 2; ispin++){
        for (iW2 = -1; iW2 <= 1; iW2++){
          for (iL2 = -1; iL2 <= 1; iL2++){

            if (iW2 == 0 && iL2 == 0) trans0 = mu;
            else if (iW2 == 0 || iL2 == 0) trans0 = t;
            else trans0 = tp;

            jsite = ((iL + iL2 + 2 * L) % L) + ((iW + iW2 + 2 * W) % W) * L;
            StdFace_trans(&ktrans, trans0, isite, ispin, jsite, ispin);
          }
        }
      }
    }
  }
  /*
  Interaction
  */
  nintr = L * W * (1 + 4 * 4);
  intrindx = (int **)malloc(sizeof(int*) * nintr);
  intr = (double *)malloc(sizeof(double) * nintr);
  for (kintr = 0; kintr < nintr; kintr++){
    intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  kintr = 0;
  for (iW = 0; iW < W; iW++){

    iW2 = (iW + 1) % W;

    for (iL = 0; iL < L; iL++){

      iL2 = (iL + 1) % L;

      isite = iL + iW * L;
      StdFace_intr(&kintr,U,isite,0,isite,0,isite,1,isite,1);

      isite = iL + iW * L;
      jsite = iL + iW2 * L;
      StdFace_Coulomb(&kintr, V, isite, jsite);

      jsite = iL2 + iW * L;
      StdFace_Coulomb(&kintr, V, isite, jsite);

      jsite = iL2 + iW2 * L;
      StdFace_Coulomb(&kintr, Vp, isite, jsite);

      isite = iL + iW2 * L;
      jsite = iL2 + iW * L;
      StdFace_Coulomb(&kintr, Vp, isite, jsite);
     }
  }
  /*
  Set mTPQ parameter
  */
  if (lGC == 0){
    LargeValue0 =  fabs(mu) * (double)nelec / (double)(L * W)
      + 2.0 * 4.0 * fabs(t) + 2.0 * 4.0 * fabs(tp)
      + fabs(U) + 2.0 * 4.0 * fabs(V) + 2.0 * 4.0 * fabs(Vp);
  }
  else{
    LargeValue0 = fabs(mu) * 2.0 + 2.0 * 4.0 * fabs(t) + 2.0 * 4.0 * fabs(tp)
      + fabs(U) + 2.0 * 4.0 * fabs(V) + 2.0 * 4.0 * fabs(Vp);
  }
  StdFace_PrintVal_i("LargeValue", &LargeValue, (int)LargeValue0 + 1);
}

/**
 *
 * Setup a Hamiltonian for the generalized Heisenberg model on a square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void Spin_SquareLattice(){
  int isite, jsite;
  int iL, iW, iW2, iL2;
  int ktrans, kintr;
  double LargeValue0;
  printf("\n");
  printf("#######  Parameter Summary  #######\n");
  printf("\n");
  StdFace_RequiredVal_i("L", L);
  StdFace_RequiredVal_i("W", W);
  StdFace_PrintVal_d("a", &a, 1.0);
  /**/
  StdFace_PrintVal_d("h", &h, 0.0);
  StdFace_PrintVal_d("Gamma", &Gamma, 0.0);
  StdFace_PrintVal_d("D", &D, 0.0);
  StdFace_PrintVal_d("J", &J, 0.0);
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
  nsite = L * W;
  locspinflag = (int *)malloc(sizeof(int) * nsite);
  for (isite = 0; isite < nsite; isite++)locspinflag[isite] = 0;
  /*
  Transfer
  */
  ntrans = L * W * 4;
  transindx = (int **)malloc(sizeof(int*) * ntrans);
  trans = (double *)malloc(sizeof(double) * ntrans);
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  ktrans = 0;
  for (isite = 0; isite < nsite; isite++){
    jsite = isite;
    StdFace_trans(&ktrans, - h * 0.5, isite, 0, jsite, 0);
    StdFace_trans(&ktrans,   h * 0.5, isite, 1, jsite, 1);
    StdFace_trans(&ktrans, 0.5 * Gamma, isite, 0, isite, 1);
    StdFace_trans(&ktrans, 0.5 * Gamma, isite, 1, isite, 0);
  }
  /*
  Interaction
  */
  nintr = L * W * (4 * 5 + 2 * 8);
  intrindx = (int **)malloc(sizeof(int*) * nintr);
  intr = (double *)malloc(sizeof(double) * nintr);
  for (kintr = 0; kintr < nintr; kintr++){
    intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  kintr = 0;
  for (iW = 0; iW < W; iW++){

    iW2 = (iW + 1) % W;

    for (iL = 0; iL < L; iL++){

      iL2 = (iL + 1) % L;

      isite = iL + iW * L;
      StdFace_SzSz(&kintr, D, isite, isite);

      jsite = iL + iW2 * L;
      StdFace_SzSz(&kintr, Jz, isite, jsite);
      StdFace_exchange(&kintr, Jxy, isite, jsite);
      StdFace_Kitaev(&kintr, 0.25 * (Jx - Jy), isite, jsite);

      jsite = iL2 + iW * L;
      StdFace_SzSz(&kintr, Jz, isite, jsite);
      StdFace_exchange(&kintr, Jxy, isite, jsite);
      StdFace_Kitaev(&kintr, 0.25 * (Jx - Jy), isite, jsite);

      jsite = iL2 + iW2 * L;
      StdFace_SzSz(&kintr, Jzp, isite, jsite);
      StdFace_exchange(&kintr, Jxyp, isite, jsite);
      StdFace_Kitaev(&kintr, 0.25 * (Jxp - Jyp), isite, jsite);

      isite = iL + iW2 * L;
      jsite = iL2 + iW * L;
      StdFace_SzSz(&kintr, Jzp, isite, jsite);
      StdFace_exchange(&kintr, Jxyp, isite, jsite);
      StdFace_Kitaev(&kintr, 0.25 * (Jxp - Jyp), isite, jsite);
    }
  }
  /*
  Set mTPQ parameter
  */
  LargeValue0 = 0.5 * fabs(h) + 0.5 * fabs(D) + 0.25 * fabs(Gamma)
    + 4.0 / 8.0 * (fabs(Jx) + fabs(Jy) + fabs(Jz))
    + 4.0 / 8.0 * (fabs(Jxp) + fabs(Jyp) + fabs(Jzp));
  StdFace_PrintVal_i("LargeValue", &LargeValue, (int)LargeValue0 + 1);
}

/**
 *
 * Setup a Hamiltonian for the Kondo lattice model on a square lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void KondoLattice_SquareLattice(
  int nelec /**< [in] The number of valence electrons */,
  int lGC /**< [in] 0 for Canonical ensemble, 1 for Grand Canonical */)
{
  int isite, jsite;
  int ispin;
  int iL, iW;
  int ktrans, kintr;
  double LargeValue0;
  /**/
  printf("\n");
  printf("#######  Parameter Summary  #######\n");
  printf("\n");
  StdFace_RequiredVal_i("L", L);
  StdFace_RequiredVal_i("W", W);
  StdFace_PrintVal_d("a", &a, 1.0);
  /**/
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
  nsite = 2 * L * W;
  locspinflag = (int *)malloc(sizeof(int) * nsite);
  for (iL = 0; iL < L * W; iL++){
    locspinflag[2 * iL] = 1;
    locspinflag[2 * iL + 1] = 0;
  }
  /*
  Transfer
  */
  ntrans = L * W * 2 * 5;
  transindx = (int **)malloc(sizeof(int*) * ntrans);
  trans = (double *)malloc(sizeof(double) * ntrans);
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }

  ktrans = 0;
  for (iW = 0; iW < W; iW++){
    for (iL = 0; iL < L; iL++){
      isite = 2 * (iL + iW * L);
      for (ispin = 0; ispin < 2; ispin++){

        StdFace_trans(&ktrans, mu, isite, ispin, isite, ispin);

        jsite = 2 * (((iL + 1) % L) + (iW * L));
        StdFace_trans(&ktrans, t, isite, ispin, jsite, ispin);
        StdFace_trans(&ktrans, t, jsite, ispin, isite, ispin);

        jsite = 2 * (iL + (((iW + 1) % W) * L));
        StdFace_trans(&ktrans, t, isite, ispin, jsite, ispin);
        StdFace_trans(&ktrans, t, jsite, ispin, isite, ispin);
      }
    }
  }
  /*
  Interaction
  */
  nintr = L * W * 6;
  intrindx = (int **)malloc(sizeof(int*) * nintr);
  intr = (double *)malloc(sizeof(double) * nintr);
  for (kintr = 0; kintr < nintr; kintr++){
    intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  kintr = 0;
  for (iW = 0; iW < W; iW++){
    for (iL = 0; iL < L; iL++){

      isite = 2 * (iL + iW * L);
      jsite = 2 * (iL + iW * L) + 1;

      StdFace_exchange(&kintr, J, isite, jsite);
      StdFace_SzSz(&kintr, J, isite, jsite);
    }
  }
  /*
  Set mTPQ parameter
  */
  if (lGC == 0){
    LargeValue0 = fabs(mu) * (double)nelec / (double)(L * W) + 2.0 * 4.0 * fabs(t) + 0.25 * fabs(J);
  }
  else{
    LargeValue0 = fabs(mu) * 2.0 + 2.0 * 4.0 * fabs(t) + 0.25 * fabs(J);
  }
  StdFace_PrintVal_i("LargeValue", &LargeValue, (int)LargeValue0 + 1);
}
