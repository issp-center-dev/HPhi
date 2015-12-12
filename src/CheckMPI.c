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

int NsiteMPI;

/**
 *
 * Define the number of sites in each PE.
 * Reduce the number of electrons, total Sz by them in the inter process region 
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void CheckMPI(struct BindStruct *X/**< [inout] */)
{
  int isite, NDimInterPE, SmallDim, SpinNum;

  NsiteMPI = X->Def.Nsite;
  X->Def.NsiteMPI=NsiteMPI;
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
      } /*if (NDimInterPE == nproc)*/
      NDimInterPE *= 4;
    } /*for (isite = NsiteMPI; isite > 0; isite--)*/
    
    if (isite == 0) {
      fprintf(stderr, "Error ! The number of PROCESS should be 4 exponent !\n");
      fprintf(stderr, "        The number of PROCESS : %d\n", nproc);
      exitMPI(-1);
    } /*if (isite == 0)*/

    switch (X->Def.iCalcModel) /*2 (inner)*/ {

    case Hubbard:

      /*X->Def.NupMPI = X->Def.Nup;
      X->Def.NdownMPI = X->Def.Ndown;*/
      /* Nup & Ndown should be differerent in each PE */
      SmallDim = myrank;
      for (isite = X->Def.Nsite; isite < NsiteMPI; isite++) {
        SpinNum = SmallDim % 4;
        SmallDim /= 4;
        if (SpinNum == 1 /*01*/) {
          X->Def.Nup -= 1;
          X->Def.Ne -= 1;
        }
        else if (SpinNum == 2 /*10*/) {
          X->Def.Ndown -= 1;
          X->Def.Ne -= 1;
        }
        else if (SpinNum == 3 /*11*/){
          X->Def.Nup -= 1;
          X->Def.Ndown -= 1;
          X->Def.Ne -= 2;
        }
      } /*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/

      break;/*case Hubbard:*/

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

      break; /*case HubbardNConserved:*/

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
            if (SpinNum == 1 /*01*/) {
              X->Def.Nup -= 1;
              X->Def.Ne -= 1;
            }
            else if (SpinNum == 2 /*10*/) {
              X->Def.Ndown -= 1;
              X->Def.Ne -= 1;
            }
            else if (SpinNum == 3 /*11*/) {
              X->Def.Nup -= 1;
              X->Def.Ndown -= 1;
              X->Def.Ne -= 2;
            }
          }
          else {
            fprintf(stderr, "\n Stop because local spin in the inter process region\n");
            exitMPI(-1);
          }
        }/*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/
      } /*if (X->Def.iCalcModel == Kondo)*/

      break; /*case KondoGC, Kondo*/

    } /*switch (X->Def.iCalcModel) 2(inner)*/

    break; /*case HubbardGC, Hubbard, HubbardNConserved, Kondo, KondoGC:*/

  case SpinGC:/********************************************************/
  case Spin:

    if (X->Def.iFlgGeneralSpin == FALSE) {

      NDimInterPE = 1;
      for (isite = NsiteMPI; isite > 0; isite--) {
        if (NDimInterPE == nproc) {
          X->Def.Nsite = isite;
          break;
        }/*if (NDimInterPE == nproc)*/
        NDimInterPE *= 2;
      }/*for (isite = X->Def.NsiteMPI; isite > 0; isite--)*/

      if (isite == 0) {
        fprintf(stderr, "Error ! The number of PROCESS should be 2-exponent !\n");
        fprintf(stderr, "        The number of PROCESS : %d\n", nproc);
        exitMPI(-1);
      }/*if (isite == 0)*/

      if (X->Def.iCalcModel == Spin) {
        /*X->Def.NeMPI = X->Def.Ne;*/

        /* Ne should be different in each PE */
        SmallDim = myrank;
        for (isite = X->Def.Nsite; isite < NsiteMPI; isite++) {
          SpinNum = SmallDim % 2;
          SmallDim /= 2;
          if (SpinNum == 0) {
            X->Def.Ndown -= 1;
          }
          else {
            X->Def.Ne -= 1;
            X->Def.Nup -= 1;
          }
        }/*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/
      }/*if (X->Def.iCalcModel == Spin)*/

    } /*if (X->Def.iFlgGeneralSpin == FALSE)*/
    else{/* General Spin */
      NDimInterPE = 1;
      for (isite = NsiteMPI; isite > 0; isite--) {
        if (NDimInterPE == nproc) {
          X->Def.Nsite = isite;
          break;
        }/*if (NDimInterPE == nproc)*/
        NDimInterPE *= X->Def.SiteToBit[isite - 1];
      }/*for (isite = X->Def.NsiteMPI; isite > 0; isite--)*/
      if (isite == 0) {
        fprintf(stderr, "Error ! The number of PROCESS is wrong !\n");
        fprintf(stderr, "        The number of PROCESS : %d\n", nproc);
        exitMPI(-1);
      }/*if (isite == 0)*/

      if (X->Def.iCalcModel == Spin) {
        /*X->Def.Total2SzMPI = X->Def.Total2Sz;*/

        /* Ne should be different in each PE */
        SmallDim = myrank;
        for (isite = X->Def.Nsite; isite < NsiteMPI; isite++) {
          SpinNum = SmallDim % X->Def.SiteToBit[isite];
          SmallDim /= X->Def.SiteToBit[isite];

          X->Def.Total2Sz += X->Def.SiteToBit[isite] - 1 - 2*SpinNum;
        }/*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/
      }/*if (X->Def.iCalcModel == Spin)*/
    }/*if (X->Def.iFlgGeneralSpin == TRUE)*/

    break; /*case SpinGC, Spin*/

  default:
    fprintf(stderr, "Error ! Wrong model !\n");
    exitMPI(-1);
  }/*switch (X->Def.iCalcModel)*/

}/*void CheckMPI*/

/**
 *
 * Print infomation of MPI parallelization
 * Modify Tpow in the inter process region
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void CheckMPI_Summary(struct BindStruct *X/**< [inout] */) {

  int isite, iproc, SmallDim, SpinNum, Nelec;
  unsigned long int idimMPI;

  fprintf(stdoutMPI, "\n\n######  MPI site separation summary  ######\n\n");
  fprintf(stdoutMPI, "  INTRA process site\n");
  fprintf(stdoutMPI, "    Site    Bit\n");
  for (isite = 0; isite < X->Def.Nsite; isite++) {
    switch (X->Def.iCalcModel) {
    case HubbardGC: 
    case Hubbard:
    case HubbardNConserved:
    case Kondo:
    case KondoGC:

      fprintf(stdoutMPI, "    %4d    %4d\n", isite, 4);
      break;

    case Spin:
    case SpinGC:

      if (X->Def.iFlgGeneralSpin == FALSE) {
        fprintf(stdoutMPI, "    %4d    %4d\n", isite, 2);
      }/*if (X->Def.iFlgGeneralSpin == FALSE)*/
      else {
        fprintf(stdoutMPI, "    %4d    %4d\n", isite, X->Def.SiteToBit[isite]);
      }/*if (X->Def.iFlgGeneralSpin == TRUE)*/

      break;

    } /*switch (X->Def.iCalcModel)*/
  } /*for (isite = 0; isite < X->Def.Nsite; isite++)*/

  fprintf(stdoutMPI, "\n  INTER process site\n");
  fprintf(stdoutMPI, "    Site    Bit\n");
  for (isite = X->Def.Nsite; isite < NsiteMPI; isite++) {
    switch (X->Def.iCalcModel) {
    case HubbardGC:
    case Hubbard:
    case HubbardNConserved:
    case Kondo:
    case KondoGC:

      fprintf(stdoutMPI, "    %4d    %4d\n", isite, 4);
      break;

    case Spin:
    case SpinGC:

      if (X->Def.iFlgGeneralSpin == FALSE) {
        fprintf(stdoutMPI, "    %4d    %4d\n", isite, 2);
      }/*if (X->Def.iFlgGeneralSpin == FALSE) */
      else {
        fprintf(stdoutMPI, "    %4d    %4d\n", isite, X->Def.SiteToBit[isite]);
      }/*if (X->Def.iFlgGeneralSpin == TRUE) */

      break;

    }/*switch (X->Def.iCalcModel)*/
  }/*for (isite = X->Def.Nsite; isite < NsiteMPI; isite++)*/

  fprintf(stdoutMPI, "\n  Process element info\n");
  fprintf(stdoutMPI, "    Process       Dimension   Nup  Ndown  Nelec  Total2Sz   State\n");

  for (iproc = 0; iproc < nproc; iproc++) {

    fprintf(stdoutMPI, "    %7d", iproc);

    if (myrank == iproc) idimMPI = X->Check.idim_max;
    else idimMPI = 0;
    fprintf(stdoutMPI, " %15ld", SumMPI_li(idimMPI));

    if (myrank == iproc) Nelec = X->Def.Nup;
    else Nelec = 0;
    fprintf(stdoutMPI, "  %4d", SumMPI_i(Nelec));

    if (myrank == iproc) Nelec = X->Def.Ndown;
    else Nelec = 0;
    fprintf(stdoutMPI, "  %5d", SumMPI_i(Nelec));

    if (myrank == iproc) Nelec = X->Def.Ne;
    else Nelec = 0;
    fprintf(stdoutMPI, "  %5d", SumMPI_i(Nelec));

    if (myrank == iproc) Nelec = X->Def.Total2Sz;
    else Nelec = 0;
    fprintf(stdoutMPI, "  %8d   ", SumMPI_i(Nelec));
    /*
     Print the configuration in the inter process region of each PE
     as a binary (excepting general spin) format.
    */
    switch (X->Def.iCalcModel) {
    case HubbardGC: /****************************************************/
    case Hubbard:
    case HubbardNConserved:
    case Kondo:
    case KondoGC:

      SmallDim = iproc;
      for (isite = X->Def.Nsite; isite < NsiteMPI; isite++) {
        SpinNum = SmallDim % 4;
        SmallDim /= 4;
        if (SpinNum == 0) fprintf(stdoutMPI, "00");
        else if (SpinNum == 1) fprintf(stdoutMPI, "01");
        else if (SpinNum == 2) fprintf(stdoutMPI, "10");
        else if (SpinNum == 3) fprintf(stdoutMPI, "11");
      } /*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/

      break;

    case Spin:
    case SpinGC:

      SmallDim = iproc;
      if (X->Def.iFlgGeneralSpin == FALSE) {
        for (isite = X->Def.Nsite; isite < NsiteMPI; isite++) {
          SpinNum = SmallDim % 2;
          SmallDim /= 2;
          fprintf(stdoutMPI, "%1d", SpinNum);
        }/*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/
      }/*if (X->Def.iFlgGeneralSpin == FALSE)*/
      else {
        SmallDim = iproc;
        for (isite = X->Def.Nsite; isite < NsiteMPI; isite++) {
          SpinNum = SmallDim % (int)X->Def.SiteToBit[isite];
          SmallDim /= X->Def.SiteToBit[isite];
          fprintf(stdoutMPI, "%1d", SpinNum);
        }/*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/
      }/*if (X->Def.iFlgGeneralSpin == TRUE)*/

      break;
    
    }/*switch (X->Def.iCalcModel)*/
    fprintf(stdoutMPI, "\n");
  }/*for (iproc = 0; iproc < nproc; iproc++)*/

  fprintf(stdoutMPI, "\n   Total dimension : %ld\n\n", SumMPI_li(X->Check.idim_max));
  /*
    Reset Tpow[DefNsite], Tpow[DefNsite + 1] ... as inter process space
  */
  switch (X->Def.iCalcModel) {
  case HubbardGC: /****************************************************/
  case Hubbard:
  case HubbardNConserved:
  case Kondo:
  case KondoGC:

    X->Def.Tpow[2 * X->Def.Nsite] = 1;
    for (isite = 2 * X->Def.Nsite + 1; isite < 2 * NsiteMPI; isite++) 
      X->Def.Tpow[isite] = X->Def.Tpow[isite - 1] * 2;
 
    break;

  case SpinGC:/********************************************************/
  case Spin:

    if (X->Def.iFlgGeneralSpin == FALSE) {

      X->Def.Tpow[X->Def.Nsite] = 1;
      for (isite = X->Def.Nsite + 1; isite < NsiteMPI; isite++)
        X->Def.Tpow[isite] = X->Def.Tpow[isite - 1] * 2;

    }/*if (X->Def.iFlgGeneralSpin == FALSE)*/
    else{

      X->Def.Tpow[X->Def.Nsite] = 1;
      for (isite = X->Def.Nsite + 1; isite < NsiteMPI; isite++)
        X->Def.Tpow[isite] = X->Def.Tpow[isite - 1] * X->Def.SiteToBit[isite - 1];
 
    }/*if (X->Def.iFlgGeneralSpin == TRUE)*/
    break;
  } /*switch (X->Def.iCalcModel)*/
}/*void CheckMPI_Summary*/

