/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include "Common.h"
#include "wrapperMPI.h"

void CheckMPI(struct BindStruct *X)
{
  int isite, NDimInterPE, SmallDim, SpinNum, NsiteMPI;

  NsiteMPI = X->Def.Nsite;

  fprintf(stdout, "Debug %d %d \n", X->Def.iCalcModel, Hubbard);
  switch (X->Def.iCalcModel) {
  case HubbardGC: /****************************************************/
  case Hubbard:
  case HubbardNConserved:
  case Kondo:
  case KondoGC:

    /*
     Define local dimension
    */
    NDimInterPE = 1;
    for (isite = NsiteMPI; isite > 0; isite--) {
      if (NDimInterPE == nproc) {
        X->Def.Nsite = isite;
        break;
      }
      NDimInterPE *= 4;
    }
    
    if (isite == 0) {
      fprintf(stdout, "Error ! The number of PROCESS should be 4 exponent !\n");
      fprintf(stdout, "        The number of PROCESS : %d\n", nproc);
      exitMPI(-1);
    }

    X->Def.Tpow[2 * X->Def.Nsite] = 1;
    for (isite = 2 * X->Def.Nsite + 1; isite < 2 * NsiteMPI; isite++) 
      X->Def.Tpow[isite] = X->Def.Tpow[isite - 1] * 2;
 
    switch (X->Def.iCalcModel) {

    case Hubbard:

      /*X->Def.NupMPI = X->Def.Nup;
      X->Def.NdownMPI = X->Def.Ndown;*/
      /* Nup & Ndown should be differerent in each PE */
      SmallDim = myrank;
      for (isite = X->Def.Nsite; isite < NsiteMPI; isite++) {
        SpinNum = SmallDim % 4;
        SmallDim /= 4;
        if (SpinNum == 1 /*01*/) X->Def.Nup -= 1;
        else if (SpinNum == 2 /*10*/) X->Def.Ndown -= 1;
        else if (SpinNum == 3 /*11*/){
          X->Def.Nup -= 1;
          X->Def.Ndown -= 1;
        }
      } /*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/

      break;

    case HubbardNConserved:
      /*X->Def.NeMPI = X->Def.Ne;*/

      /* Ne should be different in each PE */
      SmallDim = myrank;
      for (isite = X->Def.Nsite; isite < NsiteMPI; isite++) {
        SpinNum = SmallDim % 4;
        SmallDim /= 4;
        if (SpinNum == 1 /*01*/ || SpinNum == 2 /*10*/) X->Def.Ne -= 1;
        else if (SpinNum == 3 /*11*/) X->Def.Ne -= 2;
      } /*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/

      break;

    case KondoGC:
    case Kondo:

      for (isite = X->Def.Nsite; isite < NsiteMPI; isite++)
        if (X->Def.LocSpn[isite] != ITINERANT) X->Def.NLocSpn -= 1;

      if (X->Def.iCalcModel == Kondo) {
        SmallDim = myrank;
        for (isite = X->Def.Nsite; isite < NsiteMPI; isite++) {
          SpinNum = SmallDim % 4;
          SmallDim /= 4;
          if (X->Def.LocSpn[isite] == ITINERANT) {
            if (SpinNum == 1 /*01*/) X->Def.Nup -= 1;
            else if (SpinNum == 2 /*10*/) X->Def.Ndown -= 1;
            else if (SpinNum == 3 /*11*/) {
              X->Def.Nup -= 1;
              X->Def.Ndown -= 1;
            }
          }
          else {
            if (SpinNum == 0) X->Def.Nup -= 1;
            else X->Def.Ndown -= 1;
          }
        }/*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/
      } /*if (X->Def.iCalcModel == Kondo)*/

      break;

    } /*switch (X->Def.iCalcModel)*/

    break;

  case SpinGC:/********************************************************/
  case Spin:

    if (X->Def.iFlgGeneralSpin == FALSE) {

      NDimInterPE = 1;
      for (isite = NsiteMPI; isite > 0; isite--) {
        if (NDimInterPE == nproc) {
          X->Def.Nsite = isite;
          break;
        }
        NDimInterPE *= 2;
      }/*for (isite = X->Def.NsiteMPI; isite > 0; isite--)*/

      if (isite == 0) {
        fprintf(stdout, "Error ! The number of PROCESS should be 2-exponent !\n");
        fprintf(stdout, "        The number of PROCESS : %d\n", nproc);
        exitMPI(-1);
      }

      X->Def.Tpow[X->Def.Nsite] = 1;
      for (isite = X->Def.Nsite + 1; isite < NsiteMPI; isite++)
        X->Def.Tpow[isite] = X->Def.Tpow[isite - 1] * 2;

      if (X->Def.iCalcModel == Spin) {
        /*X->Def.NeMPI = X->Def.Ne;*/

        /* Ne should be different in each PE */
        SmallDim = myrank;
        for (isite = X->Def.Nsite; isite < NsiteMPI; isite++) {
          SpinNum = SmallDim % 2;
          SmallDim /= 2;
          if (SpinNum == 0 /*01 up*/) X->Def.Ne -= 1;
          else /*10 down*/ X->Def.Ne += 1;
        }/*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/
      }/*if (X->Def.iCalcModel == Spin)*/

    }
    else{/* General Spin */
      NDimInterPE = 1;
      for (isite = NsiteMPI; isite > 0; isite--) {
        if (NDimInterPE == nproc) {
          X->Def.Nsite = isite;
          break;
        }
        NDimInterPE *= X->Def.SiteToBit[isite];
      }/*for (isite = X->Def.NsiteMPI; isite > 0; isite--)*/
      if (isite == 0) {
        fprintf(stdout, "Error ! The number of PROCESS is wrong !\n");
        fprintf(stdout, "        The number of PROCESS : %d\n", nproc);
        exitMPI(-1);
      }/**/

      X->Def.Tpow[X->Def.Nsite] = 1;
      for (isite = X->Def.Nsite + 1; isite < NsiteMPI; isite++)
        X->Def.Tpow[isite] = X->Def.Tpow[isite - 1] * X->Def.SiteToBit[isite - 1];
 
      if (X->Def.iCalcModel == Spin) {
        /*X->Def.Total2SzMPI = X->Def.Total2Sz;*/

        /* Ne should be different in each PE */
        SmallDim = myrank;
        for (isite = X->Def.Nsite; isite < NsiteMPI; isite++) {
          SpinNum = SmallDim % X->Def.SiteToBit[isite];
          SmallDim /= X->Def.SiteToBit[isite];

          X->Def.Total2Sz -= X->Def.SiteToBit[isite] - 1 - 2*SmallDim;
        }/*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/
      }/*if (X->Def.iCalcModel == Spin)*/
    }/*if (X->Def.iFlgGeneralSpin == TRUE)*/

    break;

  default:
    fprintf(stdout, "Error ! Wrong model !\n");
    exitMPI(-1);
  }

}
