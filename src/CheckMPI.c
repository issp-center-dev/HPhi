/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

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
/**@file
@brief Compute total number of electrons, spins
*/
#include "Common.h"
#include "wrapperMPI.h"
/**
@brief Define the number of sites in each PE (DefineList.Nsite).
 Reduce the number of electrons (DefineList.Ne), 
 total Sz (DefineList.Total2Sz) by them in the inter process region 
@author Mitsuaki Kawamura (The University of Tokyo)
*/
int CheckMPI(struct BindStruct *X/**< [inout] */)
{
  int isite, NDimInterPE, SmallDim, SpinNum, ipivot, ishift, isiteMax, isiteMax0;

  /**@brief
  Branch for each model
  <ul>
  */
  X->Def.NsiteMPI = X->Def.Nsite;
  X->Def.Total2SzMPI = X->Def.Total2Sz;
  switch (X->Def.iCalcModel) {
  case HubbardGC: /****************************************************/
  case Hubbard:
  case HubbardNConserved:
  case Kondo:
  case KondoGC:

    /**@brief
     <li> For Hubbard & Kondo
     Define local dimension DefineList::Nsite</li>
     <ul>
    */
    NDimInterPE = 1;
    for (isite = X->Def.NsiteMPI; isite > 0; isite--) {
      if (NDimInterPE == nproc) {
        X->Def.Nsite = isite;
        break;
      } /*if (NDimInterPE == nproc)*/
      NDimInterPE *= 4;
    } /*for (isite = NsiteMPI; isite > 0; isite--)*/
    
    if (isite == 0) {
      fprintf(stdoutMPI, "%s", cErrNProcNumberHubbard);
      fprintf(stdoutMPI, cErrNProcNumber, nproc);
      NDimInterPE = 1;
      int ismallNproc=1;
      int ilargeNproc=1;
      for (isite = X->Def.NsiteMPI; isite > 0; isite--) {
        if (NDimInterPE > nproc) {
          ilargeNproc = NDimInterPE;
          if(isite >1)
            ismallNproc = NDimInterPE/4;
          break;
        }/*if (NDimInterPE > nproc)*/
        NDimInterPE *= 4;
      }/*for (isite = X->Def.NsiteMPI; isite > 0; isite--)*/
      fprintf(stdoutMPI, cErrNProcNumberSet,ismallNproc, ilargeNproc );
        return FALSE;
      //return FALSE;
    } /*if (isite == 0)*/

    switch (X->Def.iCalcModel) /*2 (inner)*/ {

    case Hubbard:
      /**@brief
      <li>For canonical Hubbard
      DefineList::Nup, DefineList::Ndown, and DefineList::Ne should be
      differerent in each PE.</li>
      */
      SmallDim = myrank;
      for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++) {
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
      /**@brief
      <li>For N-conserved canonical Hubbard
      DefineList::Ne should be differerent in each PE.</li>
      */
      SmallDim = myrank;
      for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++) {
        SpinNum = SmallDim % 4;
        SmallDim /= 4;
        if (SpinNum == 1 /*01*/ || SpinNum == 2 /*10*/) X->Def.Ne -= 1;
        else if (SpinNum == 3 /*11*/) X->Def.Ne -= 2;
      } /*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/

      break; /*case HubbardNConserved:*/

    case KondoGC:
    case Kondo:
      /**@brief
      <li>For canonical Kondo system
      DefineList::Nup, DefineList::Ndown, and DefineList::Ne should be
      differerent in each PE.</li>
      */
      for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)
        if (X->Def.LocSpn[isite] != ITINERANT) X->Def.NLocSpn -= 1;

      if (X->Def.iCalcModel == Kondo) {
        SmallDim = myrank;
        for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++) {
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
            fprintf(stdoutMPI, "\n Stop because local spin in the inter process region\n");
            return FALSE;
          }
        }/*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/
      } /*if (X->Def.iCalcModel == Kondo)*/

      break; /*case KondoGC, Kondo*/

    } /*switch (X->Def.iCalcModel) 2(inner)*/

    break; /*case HubbardGC, Hubbard, HubbardNConserved, Kondo, KondoGC:*/
    /**@brief</ul>*/
  case SpinGC:/********************************************************/
  case Spin:

    if (X->Def.iFlgGeneralSpin == FALSE) {
      /**@brief
      <li> For 1/2 Spin system,
      define local dimension DefineList::Nsite</li>
      */
      NDimInterPE = 1;
      for (isite = X->Def.NsiteMPI; isite > 0; isite--) {
        if (NDimInterPE == nproc) {
          X->Def.Nsite = isite;
          break;
        }/*if (NDimInterPE == nproc)*/
        NDimInterPE *= 2;
      }/*for (isite = X->Def.NsiteMPI; isite > 0; isite--)*/

      if (isite == 0) {
        fprintf(stdoutMPI, "%s", cErrNProcNumberSpin);
        fprintf(stdoutMPI, cErrNProcNumber, nproc);
      NDimInterPE = 1;
      int ismallNproc=1;
      int ilargeNproc=1;
      for (isite = X->Def.NsiteMPI; isite > 0; isite--) {
        if (NDimInterPE > nproc) {
          ilargeNproc = NDimInterPE;
          if(isite >1)
            ismallNproc = NDimInterPE/2;
          break;
        }/*if (NDimInterPE > nproc)*/
        NDimInterPE *= 2;
      }/*for (isite = X->Def.NsiteMPI; isite > 0; isite--)*/
      fprintf(stdoutMPI, cErrNProcNumberSet,ismallNproc, ilargeNproc );
        return FALSE;
      }/*if (isite == 0)*/

      if (X->Def.iCalcModel == Spin) {
        /*X->Def.NeMPI = X->Def.Ne;*/

        /* Ne should be different in each PE */
        SmallDim = myrank;
        for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++) {
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
      /**@brief
      <li> For general Spin system,
      define local dimension DefineList::Nsite</li>
      */
      NDimInterPE = 1;
      for (isite = X->Def.NsiteMPI; isite > 0; isite--) {
        if (NDimInterPE == nproc) {
          X->Def.Nsite = isite;
          break;
        }/*if (NDimInterPE == nproc)*/
        NDimInterPE *= X->Def.SiteToBit[isite - 1];
      }/*for (isite = X->Def.NsiteMPI; isite > 0; isite--)*/

      if (isite == 0) {
        fprintf(stdoutMPI, "%s", cErrNProcNumberGneralSpin);
        fprintf(stdoutMPI, cErrNProcNumber, nproc);
      NDimInterPE = 1;
      int ismallNproc=1;
      int ilargeNproc=1;
      for (isite = X->Def.NsiteMPI; isite > 0; isite--) {
        if (NDimInterPE > nproc) {
          ilargeNproc = NDimInterPE;
          if(isite >1)
            ismallNproc = NDimInterPE/X->Def.SiteToBit[isite - 2];
          break;
        }/*if (NDimInterPE > nproc)*/
        NDimInterPE *= X->Def.SiteToBit[isite - 1];
      }/*for (isite = X->Def.NsiteMPI; isite > 0; isite--)*/
      fprintf(stdoutMPI, cErrNProcNumberSet,ismallNproc, ilargeNproc );
        return FALSE;
      }/*if (isite == 0)*/

      if (X->Def.iCalcModel == Spin) {
        X->Def.Total2SzMPI = X->Def.Total2Sz;

        /* Ne should be different in each PE */
        SmallDim = myrank;
        for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++) {
          SpinNum = SmallDim % X->Def.SiteToBit[isite];
          SmallDim /= X->Def.SiteToBit[isite];

          X->Def.Total2Sz += X->Def.SiteToBit[isite] - 1 - 2*SpinNum;
        }/*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/
      }/*if (X->Def.iCalcModel == Spin)*/
    }/*if (X->Def.iFlgGeneralSpin == TRUE)*/

     /**@brief</ul>*/
    break; /*case SpinGC, Spin*/

  default:
    fprintf(stdoutMPI, "Error ! Wrong model !\n");
    return FALSE;
  }/*switch (X->Def.iCalcModel)*/

  /**@brief
   Check the number of processes for Boost
  */
  if (X->Boost.flgBoost == 1) {
    isiteMax = X->Boost.W0;
    ishift = 0;
    for (ipivot = 0; ipivot < X->Boost.num_pivot; ipivot++) {
      isiteMax0 = X->Boost.list_6spin_star[ipivot][1]
                + X->Boost.list_6spin_star[ipivot][2]
                + X->Boost.list_6spin_star[ipivot][3]
                + X->Boost.list_6spin_star[ipivot][4]
                + X->Boost.list_6spin_star[ipivot][5];
      if (ishift > 1) isiteMax0 = X->Def.NsiteMPI - isiteMax0 - 1 - ishift;
      else isiteMax0 = X->Def.NsiteMPI - isiteMax0 - 2;
      if (isiteMax0 < isiteMax) isiteMax = isiteMax0;
      if (X->Boost.list_6spin_star[ipivot][6] == 1) ishift += X->Boost.ishift_nspin;
    }/*for (ipivot = 0; ipivot < X->Boost.num_pivot; ipivot++)*/

    NDimInterPE = 1;
    for (isite = 0; isite < isiteMax; isite++) NDimInterPE *= 2;

    if (NDimInterPE < nproc) {
      fprintf(stderr, "\n Error ! in ReadDefFileIdxPara.\n");
      fprintf(stderr, "Too many MPI processes ! It should be <= %d. \n\n", NDimInterPE);
      exitMPI(-1);
    }/*if (NDimInterPE < nproc)*/
  }/*if (X->Boost.flgBoost == 1)*/

  return TRUE;
}/*void CheckMPI*/
/**
@brief Print infomation of MPI parallelization
Modify Definelist::Tpow in the inter process region
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void CheckMPI_Summary(struct BindStruct *X/**< [inout] */) {

  int isite, iproc, SmallDim, SpinNum, Nelec;
  unsigned long int idimMPI;

  if(X->Def.iFlgScaLAPACK == 0) {
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
            fprintf(stdoutMPI, "    %4d    %4ld\n", isite, X->Def.SiteToBit[isite]);
          }/*if (X->Def.iFlgGeneralSpin == TRUE)*/

              break;

      } /*switch (X->Def.iCalcModel)*/
    } /*for (isite = 0; isite < X->Def.Nsite; isite++)*/

    fprintf(stdoutMPI, "\n  INTER process site\n");
    fprintf(stdoutMPI, "    Site    Bit\n");
    for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++) {
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
            fprintf(stdoutMPI, "    %4d    %4ld\n", isite, X->Def.SiteToBit[isite]);
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

      if (myrank == iproc) {
        Nelec = X->Def.Ne; //X->Def.Nup
        if (X->Def.iCalcModel == Spin || X->Def.iCalcModel == SpinGC) Nelec += X->Def.Ndown;
      } else Nelec = 0;

      fprintf(stdoutMPI, "  %5d", SumMPI_i(Nelec));

      if (myrank == iproc) Nelec = X->Def.Total2Sz;
      else Nelec = 0;
      fprintf(stdoutMPI, "  %8d   ", SumMPI_i(Nelec));
      /**@brief
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
              for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++) {
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
                for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++) {
                  SpinNum = SmallDim % 2;
                  SmallDim /= 2;
                  fprintf(stdoutMPI, "%1d", SpinNum);
                }/*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/
              }/*if (X->Def.iFlgGeneralSpin == FALSE)*/
              else {
                SmallDim = iproc;
                for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++) {
                  SpinNum = SmallDim % (int) X->Def.SiteToBit[isite];
                  SmallDim /= X->Def.SiteToBit[isite];
                  fprintf(stdoutMPI, "%1d", SpinNum);
                }/*for (isite = X->Def.Nsite; isite < X->Def.NsiteMPI; isite++)*/
              }/*if (X->Def.iFlgGeneralSpin == TRUE)*/

              break;

      }/*switch (X->Def.iCalcModel)*/
      fprintf(stdoutMPI, "\n");
    }/*for (iproc = 0; iproc < nproc; iproc++)*/

    X->Check.idim_maxMPI = SumMPI_li(X->Check.idim_max);
    fprintf(stdoutMPI, "\n   Total dimension : %ld\n\n", X->Check.idim_maxMPI);
    if (X->Check.idim_maxMPI < 1) {
      fprintf(stdoutMPI, "ERROR! Total dimension < 1\n");
      exitMPI(-1);
    }
  }
  else{
    fprintf(stdoutMPI, "\n   Total dimension : %ld\n\n", X->Check.idim_max);
  }
  
  /**@brief
    Reset DefineList::Tpow[DefNsite], DefineList::Tpow[DefNsite + 1] ... 
    as inter process space
    For Hubbard & Kondo system, define DefineList::OrgTpow which is not
    affected by the number of processes.
  */
  switch (X->Def.iCalcModel) {
  case HubbardGC: /****************************************************/
  case Hubbard:
  case HubbardNConserved:
  case Kondo:
  case KondoGC:

    X->Def.Tpow[2 * X->Def.Nsite] = 1;
    for (isite = 2 * X->Def.Nsite + 1; isite < 2 * X->Def.NsiteMPI; isite++)
      X->Def.Tpow[isite] = X->Def.Tpow[isite - 1] * 2;

    X->Def.OrgTpow[0]=1;
    for (isite = 1; isite < 2 * X->Def.NsiteMPI; isite++)
      X->Def.OrgTpow[isite] = X->Def.OrgTpow[isite-1]*2;
    
    break;

  case SpinGC:/********************************************************/
  case Spin:

    if (X->Def.iFlgGeneralSpin == FALSE) {

      X->Def.Tpow[X->Def.Nsite] = 1;
      for (isite = X->Def.Nsite + 1; isite < X->Def.NsiteMPI; isite++)
        X->Def.Tpow[isite] = X->Def.Tpow[isite - 1] * 2;

    }/*if (X->Def.iFlgGeneralSpin == FALSE)*/
    else{

      X->Def.Tpow[X->Def.Nsite] = 1;
      for (isite = X->Def.Nsite + 1; isite < X->Def.NsiteMPI; isite++)
        X->Def.Tpow[isite] = X->Def.Tpow[isite - 1] * X->Def.SiteToBit[isite - 1];
 
    }/*if (X->Def.iFlgGeneralSpin == TRUE)*/
    break;
  } /*switch (X->Def.iCalcModel)*/
}/*void CheckMPI_Summary*/
