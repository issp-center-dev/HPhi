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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "StdFace_vals.h"
#include "StdFace_ModelUtil.h"
#include <complex.h>

#if defined(_HPhi)
/**
 * Set Largevalue for TPQ
 */
void StdFace_LargeValue(struct StdIntList *StdI) {
  int ktrans, kintr;
  double LargeValue0;

  LargeValue0 = 0.0;
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++) {
    LargeValue0 += cabs(StdI->trans[ktrans]);
  }
  for (kintr = 0; kintr < StdI->nintr; kintr++) {
    LargeValue0 += cabs(StdI->intr[kintr]);
  }
  for (kintr = 0; kintr < StdI->NCintra; kintr++) {
    LargeValue0 += fabs(StdI->Cintra[kintr]);
  }
  for (kintr = 0; kintr < StdI->NCinter; kintr++) {
    LargeValue0 += fabs(StdI->Cinter[kintr]);
  }
  for (kintr = 0; kintr < StdI->NEx; kintr++) {
    LargeValue0 += 2.0 * fabs(StdI->Ex[kintr]);
  }
  for (kintr = 0; kintr < StdI->NPairLift; kintr++) {
    LargeValue0 += 2.0 * fabs(StdI->PairLift[kintr]);
  }
  for (kintr = 0; kintr < StdI->NHund; kintr++) {
    LargeValue0 += 2.0 * fabs(StdI->Hund[kintr]);
  }
  LargeValue0 /= (double)StdI->nsite;
  StdFace_PrintVal_d("LargeValue", &StdI->LargeValue, LargeValue0);
}/*void StdFace_LargeValue*/
/**
 *
 * Print calcmod.def
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void PrintCalcMod(struct StdIntList *StdI)
{
  FILE *fp;
  int iCalcType, iCalcModel, iRestart, iCalcSpec, 
    iCalcEigenvec, iInitialVecTpye, InputEigenVec, OutputEigenVec;

  fprintf(stdout, "\n  @ CalcMod\n\n");
  /*
   Method
  */
  iCalcEigenvec = 0;
  if (strcmp(StdI->method, "****") == 0){
    fprintf(stdout, "ERROR ! Method is NOT specified !\n");
    StdFace_exit(-1);
  }
  else if (strcmp(StdI->method, "lanczos") == 0) iCalcType = 0;
  else if (strcmp(StdI->method, "lanczosenergy") == 0) { 
    iCalcType = 0; 
    iCalcEigenvec = 1;
  }
  else if (strcmp(StdI->method, "tpq") == 0) iCalcType = 1;
  else if (strcmp(StdI->method, "fulldiag") == 0 ||
    strcmp(StdI->method, "alldiag") == 0 ||
    strcmp(StdI->method, "direct") == 0 ) iCalcType = 2;
  else if (strcmp(StdI->method, "cg") == 0) iCalcType = 3;
  else{
    fprintf(stdout, "\n ERROR ! Unsupported Solver : %s\n", StdI->method);
    StdFace_exit(-1);
  }
  /*
   Model
  */
  if (strcmp(StdI->model, "hubbard") == 0) {
    if (StdI->lGC == 0)iCalcModel = 0;
    else iCalcModel = 3;
  }
  else if (strcmp(StdI->model, "spin") == 0) {
    if (StdI->lGC == 0)iCalcModel = 1;
    else iCalcModel = 4;
  }
  if (strcmp(StdI->model, "kondo") == 0) {
    if (StdI->lGC == 0)iCalcModel = 2;
    else iCalcModel = 5;
  }
  /*
    Restart
  */
  if (strcmp(StdI->Restart, "****") == 0) {
    strcpy(StdI->Restart, "none\0");
    fprintf(stdout, "          Restart = none        ######  DEFAULT VALUE IS USED  ######\n");
    iRestart = 0;
  }
  else {
    fprintf(stdout, "          Restart = %s\n", StdI->Restart);
    if (strcmp(StdI->Restart, "none") == 0) iRestart = 0;
    else if (strcmp(StdI->Restart, "save") == 0) iRestart = 1;
    else if (strcmp(StdI->Restart, "restartsave") == 0) iRestart = 2;
    else if (strcmp(StdI->Restart, "restart") == 0) iRestart = 3;
    else {
      fprintf(stdout, "\n ERROR ! Restart Mode : %s\n", StdI->Restart);
      StdFace_exit(-1);
    }
  }
  /*
  InitialVecType
  */
  if (strcmp(StdI->InitialVecType, "****") == 0) {
    strcpy(StdI->InitialVecType, "c\0");
    fprintf(stdout, "   InitialVecType = c           ######  DEFAULT VALUE IS USED  ######\n");
    iInitialVecTpye = 0;
  }
  else {
    fprintf(stdout, "   InitialVecType = %s\n", StdI->InitialVecType);
    if (strcmp(StdI->InitialVecType, "c") == 0) iInitialVecTpye = 0;
    else if (strcmp(StdI->InitialVecType, "r") == 0) iInitialVecTpye = 1;
    else {
      fprintf(stdout, "\n ERROR ! Restart Mode : %s\n", StdI->Restart);
      StdFace_exit(-1);
    }
  }
  /*
  EigenVecIO
  */
  InputEigenVec = 0;
  OutputEigenVec = 0;
  if (strcmp(StdI->EigenVecIO, "****") == 0) {
    strcpy(StdI->EigenVecIO, "none\0");
    fprintf(stdout, "       EigenVecIO = none        ######  DEFAULT VALUE IS USED  ######\n");
  }
  else {
    fprintf(stdout, "       EigenVecIO = %s\n", StdI->EigenVecIO);
    if (strcmp(StdI->EigenVecIO, "none") == 0) InputEigenVec = 0;
    else if (strcmp(StdI->EigenVecIO, "in") == 0) InputEigenVec = 1;
    else if (strcmp(StdI->EigenVecIO, "out") == 0) OutputEigenVec = 1;
    else if (strcmp(StdI->EigenVecIO, "inout") == 0) {
      InputEigenVec = 1;
      OutputEigenVec = 1;
    }
    else {
      fprintf(stdout, "\n ERROR ! EigenVecIO Mode : %s\n", StdI->Restart);
      StdFace_exit(-1);
    }
  }
  /*
  CalcSpec
  */
  if (strcmp(StdI->CalcSpec, "****") == 0) {
    strcpy(StdI->CalcSpec, "none\0");
    fprintf(stdout, "         CalcSpec = none        ######  DEFAULT VALUE IS USED  ######\n");
    iCalcSpec = 0;
  }
  else {
    fprintf(stdout, "         CalcSpec = %s\n", StdI->CalcSpec);
    if (strcmp(StdI->CalcSpec, "none") == 0) iCalcSpec = 0;
    else if (strcmp(StdI->CalcSpec, "normal") == 0) iCalcSpec = 1;
    else if (strcmp(StdI->CalcSpec, "noiteration") == 0) iCalcSpec = 2;
    else if (strcmp(StdI->CalcSpec, "save") == 0) iCalcSpec = 3;
    else if (strcmp(StdI->CalcSpec, "restart") == 0) iCalcSpec = 4;
    else if (strcmp(StdI->CalcSpec, "restartsave") == 0) iCalcSpec = 5;
    else {
      fprintf(stdout, "\n ERROR ! CalcSpec : %s\n", StdI->CalcSpec);
      StdFace_exit(-1);
    }
  }

  fp = fopen("calcmod.def", "w");
  fprintf(fp, "#CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 4:CG\n");
  fprintf(fp, "#CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC\n");
  fprintf(fp, "#Restart = 0:None, 1:Save, 2:Restart&Save, 3:Restart\n");
  fprintf(fp, "#CalcSpec = 0:None, 1:Normal, 2:No H*Phi, 3:Save, 4:Restart, 5:Restart&Save\n");
  fprintf(fp, "CalcType %3d\n", iCalcType);
  fprintf(fp, "CalcModel %3d\n", iCalcModel);
  fprintf(fp, "ReStart %3d\n", iRestart);
  fprintf(fp, "CalcSpec %3d\n", iCalcSpec);
  fprintf(fp, "CalcEigenVec %3d\n", iCalcEigenvec);
  fprintf(fp, "InitialVecType %3d\n", iInitialVecTpye);
  fprintf(fp, "InputEigenVec %3d\n", InputEigenVec);
  fprintf(fp, "OutputEigenVec %3d\n", OutputEigenVec);
  fflush(fp);
  fclose(fp);
  fprintf(stdout, "     calcmod.def is written.\n\n");
}/*static void PrintCalcMod*/
/**
 *
 * Print single.def or pair.def
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void PrintExcitation(struct StdIntList *StdI) {
  FILE *fp;
  int NumOp, spin[2][2], isite, ispin, icell, itau;
  double coef[2], pi, Cphase;
  double *fourier_r, *fourier_i;

  fourier_r = (double *)malloc(sizeof(double) * StdI->nsite);
  fourier_i = (double *)malloc(sizeof(double) * StdI->nsite);
  pi = acos(-1.0);

  fprintf(stdout, "\n  @ Spectrum\n\n");

  StdFace_PrintVal_d("SpectrumQW", &StdI->SpectrumQW, 0.0);
  StdFace_PrintVal_d("SpectrumQL", &StdI->SpectrumQL, 0.0);
  StdFace_PrintVal_d("SpectrumQH", &StdI->SpectrumQH, 0.0);

  if (strcmp(StdI->SpectrumType, "****") == 0) {
    strcpy(StdI->SpectrumType, "szsz\0");
    fprintf(stdout, "     SpectrumType = szsz        ######  DEFAULT VALUE IS USED  ######\n");
    NumOp = 2;
    coef[0] = 0.5;
    coef[1] = -0.5;
    spin[0][0] = 0;
    spin[0][1] = 0;
    spin[1][0] = 1;
    spin[1][1] = 1;
    StdI->SpectrumBody = 2;
  }
  else {
    fprintf(stdout, "     SpectrumType = %s\n", StdI->SpectrumType);
    if (strcmp(StdI->SpectrumType, "szsz") == 0) {
      NumOp = 2;
      coef[0] = 0.5;
      coef[1] = -0.5;
      spin[0][0] = 0;
      spin[0][1] = 0;
      spin[1][0] = 1;
      spin[1][1] = 1;
      StdI->SpectrumBody = 2;
    }
    else if (strcmp(StdI->SpectrumType, "s+s-") == 0) {
      NumOp = 1;
      coef[0] = 1.0;
      spin[0][0] = 0;
      spin[0][1] = 1;
      StdI->SpectrumBody = 2;
    }
    else if (strcmp(StdI->SpectrumType, "density") == 0) {
      NumOp = 2;
      coef[0] = 1,0;
      coef[1] = 1.0;
      spin[0][0] = 0;
      spin[0][1] = 0;
      spin[1][0] = 1;
      spin[1][1] = 1;
      StdI->SpectrumBody = 2;
    }
    else if (strcmp(StdI->SpectrumType, "up") == 0) {
      NumOp = 1;
      coef[0] = 1.0;
      spin[0][0] = 0;
      StdI->SpectrumBody = 1;
    }
    else if (strcmp(StdI->SpectrumType, "down") == 0) {
      NumOp = 1;
      coef[0] = 1.0;
      spin[0][0] = 1;
      StdI->SpectrumBody = 1;
    }
    else {
      fprintf(stdout, "\n ERROR ! SpectrumType : %s\n", StdI->SpectrumType);
      StdFace_exit(-1);
    }
  }

  if ((strcmp(StdI->model, "spin") == 0 && StdI->S2 > 1)
    || strcmp(StdI->model, "kondo") == 0) {
    printf("####################################\n");
    printf("###########  CAUTION  ##############\n");
    printf("####################################\n");
    printf("\n");
    printf(" For Kondo or S>1 system, excitation parameter file is NOT generated automatically.\n");
    printf(" Please write it by hand.\n");
  }/*if (StdI->S2 > 1 || strcmp(StdI->model, "kondo") == 0)*/

  isite = 0;
  for (icell = 0; icell < StdI->NCell; icell++) {
    for (itau = 0; itau < StdI->NsiteUC; itau++) {
      Cphase = (StdI->Cell[icell][0] + StdI->tau[itau][0])*StdI->SpectrumQW
             + (StdI->Cell[icell][1] + StdI->tau[itau][1])*StdI->SpectrumQL
             + (StdI->Cell[icell][2] + StdI->tau[itau][2])*StdI->SpectrumQH;
      fourier_r[isite] = cos(2.0*pi*Cphase);
      fourier_i[isite] = sin(2.0*pi*Cphase);
      isite += 1;
    }
  }

  if (StdI->SpectrumBody == 1) {
    fp = fopen("single.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NSingle %d\n", StdI->nsite * NumOp);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "============== Single Excitation ============\n");
    fprintf(fp, "=============================================\n");
    for (isite = 0; isite < StdI->nsite; isite++) {
      fprintf(fp, "%d %d 0 %25.15f %25.15f\n", isite, spin[0][0],
        fourier_r[isite]*coef[0], fourier_i[isite] * coef[0]);
    }
    fprintf(stdout, "      single.def is written.\n\n");
  }
  else {
    fp = fopen("pair.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NPair %d\n", StdI->nsite * NumOp);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "=============== Pair Excitation =============\n");
    fprintf(fp, "=============================================\n");
    for (isite = 0; isite < StdI->nsite; isite++) {
      for (ispin = 0; ispin < NumOp; ispin++) {
        fprintf(fp, "%d %d %d %d 0 %25.15f %25.15f\n", 
          isite, spin[ispin][0], isite, spin[ispin][1],
          fourier_r[isite] * coef[ispin], fourier_i[isite] * coef[ispin]);
      }
    }
    fprintf(stdout, "        pair.def is written.\n\n");
  }
  fflush(fp);
  fclose(fp);

  free(fourier_r);
  free(fourier_i);

}/*static void PrintExcitation()*/

#elif defined(_mVMC)
/**
 * Output Jastrow 
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void PrintJastrow(struct StdIntList *StdI) {
  FILE *fp;
  int isite, jsite, NJastrow, iJastrow, isite1, jsite1, iorb;
  int **Jastrow;

  Jastrow = (int **)malloc(sizeof(int*) * StdI->nsite);
  for (isite = 0; isite < StdI->nsite; isite++) {
    Jastrow[isite] = (int *)malloc(sizeof(int) * StdI->nsite);
    for (jsite = 0; jsite < StdI->nsite; jsite++) {
      Jastrow[isite][jsite] = StdI->Orb[isite][jsite];
    }/*for (jsite = 0; jsite < isite; jsite++)*/
  }/*for (isite = 0; isite < StdI->nsite; isite++)*/
  /*
   Symmetrize
  */
  for (iorb = 0; iorb < StdI->NOrb; iorb++) {
    for (isite = 0; isite < StdI->nsite; isite++) {
      for (jsite = 0; jsite < StdI->nsite; jsite++) {
        if (Jastrow[isite][jsite] == iorb) {
          Jastrow[jsite][isite] = Jastrow[isite][jsite];
        }
      }/*for (jsite = 0; jsite < isite; jsite++)*/
    }/*for (isite = 0; isite < StdI->nsite; isite++)*/
  }/*for (iorb = 0; iorb < StdI->NOrb; iorb++)*/
  /**/
  if (strcmp(StdI->model, "hubbard") == 0) NJastrow = 0;
  else NJastrow = -1;
  for (isite = 0; isite < StdI->nsite; isite++) {
    /*
     For Local spin
    */
    if (StdI->locspinflag[isite] != 0) {
      for (jsite = 0; jsite < StdI->nsite; jsite++) {
        Jastrow[isite][jsite] = -1;
        Jastrow[jsite][isite] = -1;
      }
      continue;
    }
    /**/
    for (jsite = 0; jsite < isite; jsite++) {
      if (Jastrow[isite][jsite] >= 0) {
        iJastrow = Jastrow[isite][jsite];
        NJastrow -= 1;
        for (isite1 = 0; isite1 < StdI->nsite; isite1++) {
          for (jsite1 = 0; jsite1 < StdI->nsite; jsite1++) {
            if (Jastrow[isite1][jsite1] == iJastrow)
              Jastrow[isite1][jsite1] = NJastrow;
          }/*for (jsite1 = 0; jsite1 < StdI->nsite; jsite1++)*/
        }/*for (isite1 = 0; isite1 < StdI->nsite; isite1++)*/
      }/*if (Jastrow[isite][jsite] >= 0)*/
    }/*for (jsite = 0; jsite < isite; jsite++)*/
  }/*for (isite = 0; isite < StdI->nsite; isite++)*/
  /**/
  NJastrow = -NJastrow;
  for (isite = 0; isite < StdI->nsite; isite++) {
    for (jsite = 0; jsite < StdI->nsite; jsite++) {
      Jastrow[isite][jsite] = -1 - Jastrow[isite][jsite];
    }/*for (jsite = 0; jsite < isite; jsite++)*/
  }/*for (isite = 0; isite < StdI->nsite; isite++)*/

  fp = fopen("jastrowidx.def", "w");
  fprintf(fp, "=============================================\n");
  fprintf(fp, "NJastrowIdx %10d\n", NJastrow);
  fprintf(fp, "ComplexType %10d\n", StdI->ComplexType);
  fprintf(fp, "=============================================\n");
  fprintf(fp, "=============================================\n");

  for (isite = 0; isite < StdI->nsite; isite++) {
    for (jsite = 0; jsite < StdI->nsite; jsite++) {
      if (isite == jsite) continue;
      fprintf(fp, "%5d  %5d  %5d\n", isite, jsite, Jastrow[isite][jsite]);
    }/*for (jsite = 0; jsite < isite; jsite++)*/
  }/*for (isite = 0; isite < StdI->nsite; isite++)*/

  if (strcmp(StdI->model, "hubbard") == 0) {
    for (iJastrow = 0; iJastrow < NJastrow; iJastrow++)
      fprintf(fp, "%5d  %5d\n", iJastrow, 1);
  }
  else if (strcmp(StdI->model, "spin") == 0) {
    fprintf(fp, "%5d  %5d\n", 0, 0);
  }
  else if (strcmp(StdI->model, "kondo") == 0) {
    fprintf(fp, "%5d  %5d\n", 0, 0);
    for (iJastrow = 1; iJastrow < NJastrow; iJastrow++)
      fprintf(fp, "%5d  %5d\n", iJastrow, 1);
  }
  else {
    printf("\nSomething wrong. \n\n");
    StdFace_exit(-1);
  }

  fflush(fp);
  fclose(fp);
  fprintf(stdout, "    jastrowidx.def is written.\n");

  for (isite = 0; isite < StdI->nsite; isite++) free(Jastrow[isite]);
  free(Jastrow);
}/*void PrintJastrow*/

 /*
 * Output Jastrow
 */
void PrintOrb(struct StdIntList *StdI) {
  FILE *fp;
  int isite, jsite, iOrb;

  fp = fopen("orbitalidx.def", "w");
  fprintf(fp, "=============================================\n");
  fprintf(fp, "NOrbitalIdx %10d\n", StdI->NOrb);
  fprintf(fp, "ComplexType %10d\n", StdI->ComplexType);
  fprintf(fp, "=============================================\n");
  fprintf(fp, "=============================================\n");

  for (isite = 0; isite < StdI->nsite; isite++) {
    for (jsite = 0; jsite < StdI->nsite; jsite++) {
      if (StdI->AntiPeriod[0] == 1 || StdI->AntiPeriod[1] == 1) {
        fprintf(fp, "%5d  %5d  %5d  %5d\n", isite, jsite, StdI->Orb[isite][jsite], StdI->AntiOrb[isite][jsite]);
      }
      else {
        fprintf(fp, "%5d  %5d  %5d\n", isite, jsite, StdI->Orb[isite][jsite]);
      }
    }/*for (jsite = 0; jsite < isite; jsite++)*/
  }/*for (isite = 0; isite < StdI->nsite; isite++)*/

  for (iOrb = 0; iOrb < StdI->NOrb; iOrb++)
    fprintf(fp, "%5d  %5d\n", iOrb, 1);

  fflush(fp);
  fclose(fp);
  fprintf(stdout, "    orbitalidx.def is written.\n");

  for (isite = 0; isite < StdI->nsite; isite++) free(StdI->Orb[isite]);
  free(StdI->Orb);
}/*void PrintJastrow*/
/**
* Output .def file for Gutzwiller
*
*/
static void PrintGutzwiller(struct StdIntList *StdI)
{
  FILE *fp;
  int isite, jsite, NGutzwiller, iGutz;
  int *Gutz;

  Gutz = (int *)malloc(sizeof(int) * StdI->nsite);
  for (isite = 0; isite < StdI->nsite; isite++) Gutz[isite] = StdI->Orb[isite][isite];

  if (strcmp(StdI->model, "hubbard") == 0) NGutzwiller = 0;
  else NGutzwiller = -1;
  for (isite = 0; isite < StdI->nsite; isite++) {
    /*
    For Local spin
    */
    if (StdI->locspinflag[isite] != 0) {
      Gutz[isite] = -1;
      continue;
    }
    /**/
    if (Gutz[isite] >= 0) {
      iGutz = Gutz[isite];
      NGutzwiller -= 1;
      for (jsite = 0; jsite < StdI->nsite; jsite++) {
        if (Gutz[jsite] == iGutz)
          Gutz[jsite] = NGutzwiller;
      }/*for (jsite = 0; jsite < StdI->nsite; jsite++)*/
    }/*if (Gutz[isite] >= 0)*/
  }/*for (isite = 0; isite < StdI->nsite; isite++)*/
  /**/
  NGutzwiller = -NGutzwiller;
  for (isite = 0; isite < StdI->nsite; isite++) {
    Gutz[isite] = -1 - Gutz[isite];
  }/*for (isite = 0; isite < StdI->nsite; isite++)*/

  fp = fopen("gutzwilleridx.def", "w");
  fprintf(fp, "=============================================\n");
  fprintf(fp, "NGutzwillerIdx %10d\n", NGutzwiller);
  fprintf(fp, "ComplexType %10d\n", StdI->ComplexType);
  fprintf(fp, "=============================================\n");
  fprintf(fp, "=============================================\n");

  for (isite = 0; isite < StdI->nsite; isite++)
    fprintf(fp, "%5d  %5d\n", isite, Gutz[isite]);

  if (strcmp(StdI->model, "hubbard") == 0) {
    for (iGutz = 0; iGutz < NGutzwiller; iGutz++)
      fprintf(fp, "%5d  %5d\n", iGutz, 1);
  }
  else if (strcmp(StdI->model, "spin") == 0) {
    fprintf(fp, "%5d  %5d\n", 0, 0);
  }
  else if (strcmp(StdI->model, "kondo") == 0) {
    fprintf(fp, "%5d  %5d\n", 0, 0);
    for (iGutz = 1; iGutz < NGutzwiller; iGutz++)
      fprintf(fp, "%5d  %5d\n", iGutz, 1);
  }
  else {
    printf("\nSomething wrong. \n\n");
    StdFace_exit(-1);
  }
  fflush(fp);
  fclose(fp);
  fprintf(stdout, "    gutzwilleridx.def is written.\n");

  free(Gutz);
}/*static void PrintGutzwiller*/
#endif

/**
*
* Clear grobal variables in the standard mode
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
static void StdFace_ResetVals(struct StdIntList *StdI) {
  int i, j;
  double NaN_d;

  NaN_d = 0.0 / 0.0;
  StdI->NaN_i = 2147483647;
  /**/
  StdI->a = NaN_d;
  for (i = 0; i < 3; i++) StdI->length[i] = NaN_d;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      StdI->box[i][j] = StdI->NaN_i;
  StdI->Gamma = NaN_d;
  StdI->h = NaN_d;
  StdI->Height = StdI->NaN_i;
  StdI->JAll = NaN_d;
  StdI->JpAll = NaN_d;
  StdI->JppAll = NaN_d;
  StdI->J0All = NaN_d;
  StdI->J0pAll = NaN_d;
  StdI->J1All = NaN_d;
  StdI->J1pAll = NaN_d;
  StdI->J2All = NaN_d;
  StdI->J2pAll = NaN_d;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      StdI->J[i][j] = NaN_d;
      StdI->Jp[i][j] = NaN_d;
      StdI->Jpp[i][j] = NaN_d;
      StdI->J0[i][j] = NaN_d;
      StdI->J0p[i][j] = NaN_d;
      StdI->J1[i][j] = NaN_d;
      StdI->J1p[i][j] = NaN_d;
      StdI->J2[i][j] = NaN_d;
      StdI->J2p[i][j] = NaN_d;
      StdI->D[i][j] = 0.0;
    }
  }
  StdI->D[2][2] = NaN_d;
  StdI->K = NaN_d;
  StdI->L = StdI->NaN_i;
  for (i = 0; i < 3; i++) 
    for (j = 0; j < 3; j++)
      StdI->direct[i][j] = NaN_d;
  StdI->mu = NaN_d;
  StdI->S2 = StdI->NaN_i;
  StdI->t = NaN_d;
  StdI->tp = NaN_d;
  StdI->tpp = NaN_d;
  StdI->t0 = NaN_d;
  StdI->t0p = NaN_d;
  StdI->t1 = NaN_d;
  StdI->t1p = NaN_d;
  StdI->t2 = NaN_d;
  StdI->t2p = NaN_d;
  StdI->U = NaN_d;
  StdI->V = NaN_d;
  StdI->Vp = NaN_d;
  StdI->Vpp = NaN_d;
  StdI->V0 = NaN_d;
  StdI->V0p = NaN_d;
  StdI->V1 = NaN_d;
  StdI->V1p = NaN_d;
  StdI->V2 = NaN_d;
  StdI->V2p = NaN_d;
  StdI->W = StdI->NaN_i;
  StdI->phase[0] = NaN_d;
  StdI->phase[1] = NaN_d;
  StdI->phase[2] = NaN_d;
  StdI->pi180 = 0.01745329251994329576;/*Pi/180*/

  StdI->nelec = StdI->NaN_i;
  StdI->Sz2 = StdI->NaN_i;
  strcpy(StdI->model, "****\0");
  strcpy(StdI->lattice, "****\0");
  strcpy(StdI->outputmode, "****\0");
  strcpy(StdI->CDataFileHead, "****\0");
  strcpy(StdI->W90_geom, "****\0");
  strcpy(StdI->W90_hr, "****\0");
  StdI->W90_cutoff = NaN_d;
#if defined(_HPhi)
  StdI->LargeValue = NaN_d;
  StdI->OmegaMax = NaN_d;
  StdI->OmegaMin = NaN_d;
  StdI->OmegaIm = NaN_d;
  StdI->Nomega = StdI->NaN_i;
  StdI->SpectrumQW = NaN_d;
  StdI->SpectrumQL = NaN_d;
  StdI->SpectrumQH = NaN_d;
  strcpy(StdI->method, "****\0");
  strcpy(StdI->Restart, "****\0");
  strcpy(StdI->EigenVecIO, "****\0");
  strcpy(StdI->InitialVecType, "****\0");
  strcpy(StdI->CalcSpec, "****\0");
  strcpy(StdI->SpectrumType, "****\0");
  StdI->FlgTemp = 1;
  StdI->Lanczos_max = StdI->NaN_i;
  StdI->initial_iv = StdI->NaN_i;
  StdI->nvec = StdI->NaN_i;
  StdI->exct = StdI->NaN_i;
  StdI->LanczosEps = StdI->NaN_i;
  StdI->LanczosTarget = StdI->NaN_i;
  StdI->NumAve = StdI->NaN_i;
  StdI->ExpecInterval = StdI->NaN_i;
#elif defined(_mVMC)
  strcpy(StdI->CParaFileHead, "****\0");
  StdI->NVMCCalMode = StdI->NaN_i;
  StdI->NLanczosMode = StdI->NaN_i;
  StdI->NDataIdxStart = StdI->NaN_i;
  StdI->NDataQtySmp = StdI->NaN_i;
  StdI->NSPGaussLeg = StdI->NaN_i;
  StdI->NMPTrans = StdI->NaN_i;
  StdI->NSROptItrStep = StdI->NaN_i;
  StdI->NSROptItrSmp = StdI->NaN_i;
  StdI->DSROptRedCut = NaN_d;
  StdI->DSROptStaDel = NaN_d;
  StdI->DSROptStepDt = NaN_d;
  StdI->NVMCWarmUp = StdI->NaN_i;
  StdI->NVMCInterval = StdI->NaN_i;
  StdI->NVMCSample = StdI->NaN_i;
  StdI->NExUpdatePath = StdI->NaN_i;
  StdI->RndSeed = StdI->NaN_i;
  StdI->NSplitSize = StdI->NaN_i;
  StdI->NStore = StdI->NaN_i;
  StdI->ComplexType = StdI->NaN_i;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      StdI->boxsub[i][j] = StdI->NaN_i;
  StdI->Hsub = StdI->NaN_i;
  StdI->Lsub = StdI->NaN_i;
  StdI->Wsub = StdI->NaN_i;
#endif
}/*static void StdFace_ResetVals*/
/*
 * Make all characters lower
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void Text2Lower(char *value /**< [inout] Keyword or value*/){
  char value2;
  int valuelen, ii;

  valuelen = strlen(value);
  for (ii = 0; ii < valuelen; ii++) {
    value2 = tolower(value[ii]);
    value[ii] = value2;
  }
}/*static void Text2Lower*/
/**
 *
 * Remove : space etc. from keyword and value in an iput file
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void TrimSpaceQuote(char *value /**< [inout] Keyword or value*/){
  char value2[256];
  int valuelen, valuelen2, ii;

  valuelen = strlen(value);
  valuelen2 = 0;
  for (ii = 0; ii < valuelen; ii++){
    if (value[ii] != ' ' &&
      value[ii] != ':' &&
      value[ii] != ';' &&
      value[ii] != '\"' &&
      value[ii] != '\b' &&
      value[ii] != '\\' &&
      value[ii] != '\v' &&
      value[ii] != '\n' &&
      value[ii] != '\0'){
      value2[valuelen2] = value[ii];
      valuelen2++;
    }
  }

  strncpy(value, value2, valuelen2);
  value[valuelen2] = '\0';

}/*static void TrimSpaceQuote*/
/**
 *
 * Store an input value into the valiable (string)
 * If duplicated, HPhi will stop.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void StoreWithCheckDup_s(
  char *keyword /**< [in] keyword read from the input file*/, 
  char *valuestring /**< [in] value read from the input file*/, 
  char *value /**< [out] */)
{
  if (strcmp(value, "****") != 0){
    fprintf(stdout, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    StdFace_exit(-1);
  }
  else{
    strcpy(value, valuestring);
  }
}/*static void StoreWithCheckDup_s*/
 /**
 *
 * Store an input value into the valiable (string) 
 * Force string lower.
 * If duplicated, HPhi will stop.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void StoreWithCheckDup_sl(
  char *keyword /**< [in] keyword read from the input file*/,
  char *valuestring /**< [in] value read from the input file*/,
  char *value /**< [out] */)
{
  if (strcmp(value, "****") != 0) {
    fprintf(stdout, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    StdFace_exit(-1);
  }
  else {
    strcpy(value, valuestring);
    Text2Lower(value);
  }
}/*static void StoreWithCheckDup_sl*/
/**
 *
 * Store an input value into the valiable (integer)
 * If duplicated, HPhi will stop.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void StoreWithCheckDup_i(
  char *keyword /**< [in] keyword read from the input file*/,
  char *valuestring /**< [in] value read from the input file*/,
  int *value /**< [out] */)
{
  int NaN_i = 2147483647;

  if (*value != NaN_i){
    fprintf(stdout, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    StdFace_exit(-1);
  }
  else{
    sscanf(valuestring, "%d", value);
  }
}/*static void StoreWithCheckDup_i*/

/**
 *
 * Store an input value into the valiable (double)
 * If duplicated, HPhi will stop.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void StoreWithCheckDup_d(
  char *keyword /**< [in] keyword read from the input file*/,
  char *valuestring /**< [in] value read from the input file*/,
  double *value /**< [out] */)
{

  if (isnan(*value) == 0){
    fprintf(stdout, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    StdFace_exit(-1);
  }
  else{
    sscanf(valuestring, "%lf", value);
  }

}/*static void StoreWithCheckDup_d*/

/**
*
* Store an input value into the valiable (Double complex)
* If duplicated, HPhi will stop.
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
static void StoreWithCheckDup_c(
  char *keyword /**< [in] keyword read from the input file*/,
  char *valuestring /**< [in] value read from the input file*/,
  double complex *value /**< [out] */)
{
  int num;
  char *valuestring_r, *valuestring_i;
  double value_r, value_i;

  if (isnan(creal(*value)) == 0) {
    fprintf(stdout, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    StdFace_exit(-1);
  }
  else {

    if (valuestring[0] == ',') {
      valuestring_r = NULL;
      valuestring_i = strtok(valuestring, ",");
    }
    else {
      valuestring_r = strtok(valuestring, ",");
      valuestring_i = strtok(NULL, ",");
    }
    
    if (valuestring_r == NULL) {
      *value = 0.0;
    }
    else {
      num = sscanf(valuestring_r, "%lf", &value_r);
      if (num == 1) *value = value_r;
      else *value = 0.0;
    }

    if (valuestring_i == NULL) {
      *value += I * 0.0;
    }
    else {
        num = sscanf(valuestring_i, "%lf", &value_i);
      if (num == 1) *value += I * value_i;
      else *value += I * 0.0;
    }
  }
}/*static void StoreWithCheckDup_c*/

/**
*
* Print the locspin file
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void PrintLocSpin(struct StdIntList *StdI) {
  FILE *fp;
  int isite, nlocspin;

  nlocspin = 0;
  for (isite = 0; isite < StdI->nsite; isite++)
    if (StdI->locspinflag[isite] != 0) nlocspin = nlocspin + 1;

  fp = fopen("locspn.def", "w");
  fprintf(fp, "================================ \n");
  fprintf(fp, "NlocalSpin %5d  \n", nlocspin);
  fprintf(fp, "================================ \n");
  fprintf(fp, "========i_0LocSpn_1IteElc ====== \n");
  fprintf(fp, "================================ \n");

  for (isite = 0; isite < StdI->nsite; isite++)
    fprintf(fp, "%5d  %5d\n", isite, StdI->locspinflag[isite]);

  fflush(fp);
  fclose(fp);
  fprintf(stdout, "    locspn.def is written.\n");
}/*void PrintLocSpin*/
/**
*
* Print the transfer file
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintTrans(struct StdIntList *StdI){
  FILE *fp;
  int jtrans, ktrans, ntrans0;

  for (jtrans = 0; jtrans < StdI->ntrans; jtrans++){
    for (ktrans = jtrans + 1; ktrans < StdI->ntrans; ktrans++){
      if (StdI->transindx[jtrans][0] == StdI->transindx[ktrans][0]
        && StdI->transindx[jtrans][1] == StdI->transindx[ktrans][1]
        && StdI->transindx[jtrans][2] == StdI->transindx[ktrans][2]
        && StdI->transindx[jtrans][3] == StdI->transindx[ktrans][3]){
        StdI->trans[jtrans] = StdI->trans[jtrans] + StdI->trans[ktrans];
        StdI->trans[ktrans] = 0.0;
      }
    }/*for (ktrans = jtrans + 1; ktrans < StdI->ntrans; ktrans++)*/
  }/*for (jtrans = 0; jtrans < StdI->ntrans; jtrans++)*/

  ntrans0 = 0;
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    if (cabs(StdI->trans[ktrans]) > 0.000001) ntrans0 = ntrans0 + 1;
  }
  if (ntrans0 != 0 || StdI->lBoost == 1) StdI->Ltrans = 1;
  else StdI->Ltrans = 0;

  if(StdI->Ltrans == 1){
    fp = fopen("trans.def", "w");
    fprintf(fp, "======================== \n");
    fprintf(fp, "NTransfer %7d  \n", ntrans0);
    fprintf(fp, "======================== \n");
    fprintf(fp, "========i_j_s_tijs====== \n");
    fprintf(fp, "======================== \n");

    ntrans0 = 0;
    for (ktrans = 0; ktrans < StdI->ntrans; ktrans++) {
      if (cabs(StdI->trans[ktrans]) > 0.000001)
        fprintf(fp, "%5d %5d %5d %5d %25.15f %25.15f\n",
          StdI->transindx[ktrans][0], StdI->transindx[ktrans][1],
          StdI->transindx[ktrans][2], StdI->transindx[ktrans][3],
          creal(StdI->trans[ktrans]), cimag(StdI->trans[ktrans]));
    }

    fflush(fp);
    fclose(fp);
    fprintf(stdout, "      trans.def is written.\n");
  }/*if (StdI->Ltrans == 1)*/
}/*static void PrintTrans*/
/**
 *
 * Print namelist.def  
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void PrintNamelist(struct StdIntList *StdI){
  FILE *fp;

  fp = fopen("namelist.def", "w");
  fprintf(                         fp, "         ModPara  modpara.def\n");
  fprintf(                         fp, "         LocSpin  locspn.def\n");
  if (StdI->Ltrans == 1) fprintf(  fp, "           Trans  trans.def\n");
  if (StdI->LCintra == 1) fprintf( fp, "    CoulombIntra  coulombintra.def\n");
  if (StdI->LCinter == 1) fprintf( fp, "    CoulombInter  coulombinter.def\n");
  if (StdI->LHund == 1)fprintf(    fp, "            Hund  hund.def\n");
  if (StdI->LEx == 1)fprintf(      fp, "        Exchange  exchange.def\n");
  if (StdI->LPairLift == 1)fprintf(fp, "        PairLift  pairlift.def\n");
  if (StdI->Lintr == 1)fprintf(    fp, "        InterAll  interall.def\n");
  if (StdI->ioutputmode != 0) {
    fprintf(                       fp, "        OneBodyG  greenone.def\n");
    fprintf(                       fp, "        TwoBodyG  greentwo.def\n");
  }
#if defined(_HPhi)
  fprintf(                         fp, "         CalcMod  calcmod.def\n");
  if(StdI->SpectrumBody == 1) 
    fprintf(                       fp, "SingleExcitation  single.def\n");
  else fprintf(                    fp, "  PairExcitation  pair.def\n");
  fprintf(                         fp, "     SpectrumVec  %s_eigenvec_0\n",
                                   StdI->CDataFileHead);
  if (StdI->lBoost == 1) fprintf(  fp, "           Boost  boost.def\n");
#elif defined(_mVMC)
  fprintf(                         fp, "      Gutzwiller  gutzwilleridx.def\n");
  fprintf(                         fp, "         Jastrow  jastrowidx.def\n");
  fprintf(                         fp, "         Orbital  orbitalidx.def\n");
  fprintf(                         fp, "        TransSym  qptransidx.def\n");
#endif
  
  fflush(fp);
  fclose(fp);
  fprintf(stdout, "    namelist.def is written.\n");
}/*static void PrintNamelist*/
/**
 *
 * Print modpara.def
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void PrintModPara(struct StdIntList *StdI)
{
  FILE *fp;

  fp = fopen("modpara.def", "w");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "Model_Parameters   0\n");
  fprintf(fp, "--------------------\n");
#if defined(_HPhi)
  fprintf(fp, "HPhi_Cal_Parameters\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "CDataFileHead  %s\n", StdI->CDataFileHead);
  fprintf(fp, "CParaFileHead  zqp\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "Nsite          %-5d\n", StdI->nsite);
  if (StdI->Sz2 != StdI->NaN_i) fprintf(fp, "2Sz            %-5d\n", StdI->Sz2);
  if (StdI->nelec != StdI->NaN_i) fprintf(fp, "Ncond          %-5d\n", StdI->nelec);
  fprintf(fp, "Lanczos_max    %-5d\n", StdI->Lanczos_max);
  fprintf(fp, "initial_iv     %-5d\n", StdI->initial_iv);
  if(StdI->nvec != StdI->NaN_i) fprintf(fp, "nvec           %-5d\n", StdI->nvec);
  fprintf(fp, "exct           %-5d\n", StdI->exct);
  fprintf(fp, "LanczosEps     %-5d\n", StdI->LanczosEps);
  fprintf(fp, "LanczosTarget  %-5d\n", StdI->LanczosTarget);
  fprintf(fp, "LargeValue     %-25.15e\n", StdI->LargeValue);
  fprintf(fp, "NumAve         %-5d\n", StdI->NumAve);
  fprintf(fp, "ExpecInterval  %-5d\n", StdI->ExpecInterval);
  fprintf(fp, "NOmega         %-5d\n", StdI->Nomega);
  fprintf(fp, "OmegaMax       %-25.15e\n", StdI->OmegaMax);
  fprintf(fp, "OmegaMin       %-25.15e\n", StdI->OmegaMin);
  fprintf(fp, "OmegaIm        %-25.15e\n", StdI->OmegaIm);
#elif defined(_mVMC)
  fprintf(fp, "VMC_Cal_Parameters\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "CDataFileHead  %s\n", StdI->CDataFileHead);
  fprintf(fp, "CParaFileHead  %s\n", StdI->CParaFileHead);
  fprintf(fp, "--------------------\n");
  fprintf(fp, "NVMCCalMode    %d\n", StdI->NVMCCalMode);
  /*fprintf(fp, "NLanczosMode   %d\n", StdI->NLanczosMode);*/
  fprintf(fp, "--------------------\n");
  fprintf(fp, "NDataIdxStart  %d\n", StdI->NDataIdxStart);
  fprintf(fp, "NDataQtySmp    %d\n", StdI->NDataQtySmp);
  fprintf(fp, "--------------------\n");
  fprintf(fp, "Nsite          %d\n", StdI->nsite);
  fprintf(fp, "Nelectron      %d\n", StdI->nelec);
  fprintf(fp, "NSPGaussLeg    %d\n", StdI->NSPGaussLeg);
  fprintf(fp, "NSPStot        %d\n", StdI->Sz2);
  fprintf(fp, "NMPTrans       %d\n", StdI->NMPTrans);
  fprintf(fp, "NSROptItrStep  %d\n", StdI->NSROptItrStep);
  fprintf(fp, "NSROptItrSmp   %d\n", StdI->NSROptItrSmp);
  fprintf(fp, "DSROptRedCut   %.10f\n", StdI->DSROptRedCut);
  fprintf(fp, "DSROptStaDel   %.10f\n", StdI->DSROptStaDel);
  fprintf(fp, "DSROptStepDt   %.10f\n", StdI->DSROptStepDt);
  fprintf(fp, "NVMCWarmUp     %d\n", StdI->NVMCWarmUp);
  fprintf(fp, "NVMCInterval  %d\n", StdI->NVMCInterval);
  fprintf(fp, "NVMCSample     %d\n", StdI->NVMCSample);
  fprintf(fp, "NExUpdatePath  %d\n", StdI->NExUpdatePath);
  fprintf(fp, "RndSeed        %d\n", StdI->RndSeed);
  fprintf(fp, "NSplitSize     %d\n", StdI->NSplitSize);
  fprintf(fp, "NStore         %d\n", StdI->NStore);
#endif

  fflush(fp);
  fclose(fp);
  fprintf(stdout, "     modpara.def is written.\n");
}/*static void PrintModPara*/

/**
 *
 * Print greenone.def
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void Print1Green(struct StdIntList *StdI)
{
  FILE *fp;
  int ngreen, igreen, store;
  int isite, jsite, ispin, jspin, SiMax, SjMax;
  int **greenindx;
  /*
   Set Indices of correlation functions
  */
  ngreen = 0;
  if (StdI->ioutputmode != 0) {
    for (store = 0; store < 2; store++) {

      if (store == 1) {
        greenindx = (int **)malloc(sizeof(int*) * (ngreen + 1));
        for (igreen = 0; igreen < ngreen; igreen++) {
          greenindx[igreen] = (int *)malloc(sizeof(int) * 4);
        }
        ngreen = 0;
      }/*if (store == 1)*/

      for (isite = 0; isite < StdI->nsite; isite++) {

        if (StdI->locspinflag[isite] == 0) SiMax = 1;
        else SiMax = StdI->locspinflag[isite];

        for (ispin = 0; ispin <= SiMax; ispin++) {
          for (jsite = 0; jsite < StdI->nsite; jsite++) {

            if (StdI->locspinflag[jsite] == 0) SjMax = 1;
            else SjMax = StdI->locspinflag[jsite];

            for (jspin = 0; jspin <= SjMax; jspin++) {

              if (isite != jsite &&
                (StdI->locspinflag[isite] != 0 && StdI->locspinflag[jsite] != 0)) continue;

              if (StdI->ioutputmode == 2 || ispin == jspin)
              {
                if (store == 1) {
                  greenindx[ngreen][0] = isite;
                  greenindx[ngreen][1] = ispin;
                  greenindx[ngreen][2] = jsite;
                  greenindx[ngreen][3] = jspin;
                }
                ngreen++;
              }

            }/*for (jspin = 0; jspin <= SjMax; jspin++)*/
          }/*for (jsite = 0; jsite < StdI->nsite; jsite++)*/
        }/*for (ispin = 0; ispin <= SiMax; ispin++)*/
      }/*for (isite = 0; isite < StdI->nsite; isite++)*/
    }

    fp = fopen("greenone.def", "w");
    fprintf(fp, "===============================\n");
    fprintf(fp, "NCisAjs %10d\n", ngreen);
    fprintf(fp, "===============================\n");
    fprintf(fp, "======== Green functions ======\n");
    fprintf(fp, "===============================\n");
    for (igreen = 0; igreen < ngreen; igreen++) {
    fprintf(fp,   "%5d %5d %5d %5d\n",
      greenindx[igreen][0], greenindx[igreen][1], greenindx[igreen][2], greenindx[igreen][3]);
    }
    fflush(fp);
    fclose(fp);

    fprintf(stdout, "    greenone.def is written.\n");

    for (igreen = 0; igreen < ngreen; igreen++) {
      free(greenindx[igreen]);
    }
    free(greenindx);

  }/*if (StdI->ioutputmode != 0) */
}/*static void Print1Green*/

/**
 *
 * Print greentwo.def
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void Print2Green(struct StdIntList *StdI) {
  FILE *fp;
  int ngreen, store, igreen;
  int site1, site2, site3, site4;
  int spin1, spin2, spin3, spin4;
  int S1Max, S2Max, S3Max, S4Max;
  int **greenindx;
  /*
   Set Indices of correlation functions
  */
  ngreen = 0;
  if (StdI->ioutputmode != 0) {
    for (store = 0; store < 2; store++) {

      if (store == 1) {
        greenindx = (int **)malloc(sizeof(int*) * (ngreen + 1));
        for (igreen = 0; igreen < ngreen; igreen++) {
          greenindx[igreen] = (int *)malloc(sizeof(int) * 8);
        }
        ngreen = 0;
      }/*if (store == 1)*/

      for (site1 = 0; site1 < StdI->nsite; site1++) {

        if (StdI->locspinflag[site1] == 0) S1Max = 1;
        else S1Max = StdI->locspinflag[site1];

        for (spin1 = 0; spin1 <= S1Max; spin1++) {
          for (site2 = 0; site2 < StdI->nsite; site2++) {

            if (StdI->locspinflag[site1] != 0 && StdI->locspinflag[site2] != 0
              && site1 != site2) continue;

            if (StdI->locspinflag[site2] == 0) S2Max = 1;
            else S2Max = StdI->locspinflag[site2];

            for (spin2 = 0; spin2 <= S2Max; spin2++) {
              for (site3 = 0; site3 < StdI->nsite; site3++) {

                if (StdI->locspinflag[site3] == 0) S3Max = 1;
                else S3Max = StdI->locspinflag[site3];

                for (spin3 = 0; spin3 <= S3Max; spin3++) {
                  for (site4 = 0; site4 < StdI->nsite; site4++) {

                    if (StdI->locspinflag[site3] != 0 && StdI->locspinflag[site4] != 0
                      && site3 != site4) continue;

                    if (StdI->locspinflag[site4] == 0) S4Max = 1;
                    else S4Max = StdI->locspinflag[site4];

                    for (spin4 = 0; spin4 <= S4Max; spin4++) {

                      if (StdI->ioutputmode == 2 ||
                        (site1 == site2 && site3 == site4 &&
                          spin1 - spin2 + spin3 - spin4 == 0)) {
                        if (store == 1) {
#if defined(_mVMC)
                          if (spin1 != spin2 || spin3 != spin4) 
                          {
                            greenindx[ngreen][0] = site1;
                            greenindx[ngreen][1] = spin1;
                            greenindx[ngreen][2] = site4;
                            greenindx[ngreen][3] = spin4;
                            greenindx[ngreen][4] = site3;
                            greenindx[ngreen][5] = spin3;
                            greenindx[ngreen][6] = site2;
                            greenindx[ngreen][7] = spin2;
                          }
                          else 
#endif
                          {
                            greenindx[ngreen][0] = site1;
                            greenindx[ngreen][1] = spin1;
                            greenindx[ngreen][2] = site2;
                            greenindx[ngreen][3] = spin2;
                            greenindx[ngreen][4] = site3;
                            greenindx[ngreen][5] = spin3;
                            greenindx[ngreen][6] = site4;
                            greenindx[ngreen][7] = spin4;
                          }
                        }
                        ngreen++;
                      }

                    }/*for (spin4 = 0; spin4 <= S4Max; spin4++)*/
                  }/*for (site4 = 0; site4 < StdI->nsite; site4++)*/
                }/*for (spin3 = 0; spin3 <= S3Max; spin3++*/
              }/*for (site3 = 0; site3 < StdI->nsite; site3++)*/
            }/*for (spin2 = 0; spin2 <= S2Max; spin2++)*/
          }/*for (site2 = 0; site2 < StdI->nsite; site2++)*/
        }/*for (spin1 = 0; spin1 <= S1Max; spin1++)*/
      }/*for (site1 = 0; site1 < StdI->nsite; site1++)*/

    }/*for (store = 0; store < 2; store++)*/

    fp = fopen("greentwo.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NCisAjsCktAltDC %10d\n", ngreen);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "======== Green functions for Sq AND Nq ======\n");
    fprintf(fp, "=============================================\n");
    for (igreen = 0; igreen < ngreen; igreen++) {
      fprintf(fp, "%5d %5d %5d %5d %5d %5d %5d %5d\n",
        greenindx[igreen][0], greenindx[igreen][1], greenindx[igreen][2], greenindx[igreen][3],
        greenindx[igreen][4], greenindx[igreen][5], greenindx[igreen][6], greenindx[igreen][7]);
    }
    fflush(fp);
    fclose(fp);

    fprintf(stdout, "    greentwo.def is written.\n");

    for (igreen = 0; igreen < ngreen; igreen++) {
      free(greenindx[igreen]);
    }
    free(greenindx);

  }/*if (StdI->ioutputmode != 0)*/
}/*void Print2Green(struct StdIntList *StdI)*/
/**
 *
 * Stop HPhi if unsupported model is read 
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void UnsupportedSystem(
  char *model /**< [in]*/, 
  char *lattice /**< [in]*/)
{
  fprintf(stdout, "\nSorry, specified combination, \n");
  fprintf(stdout, "    MODEL : %s  \n", model);
  fprintf(stdout, "  LATTICE : %s, \n", lattice);
  fprintf(stdout, "is unsupported in the STANDARD MODE...\n");
  fprintf(stdout, "Please use the EXPART MODE, or write a NEW FUNCTION and post us.\n");
  StdFace_exit(-1);
}/*static void UnsupportedSystem*/
/**
 *
 * Verify outputmode
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void CheckOutputMode(struct StdIntList *StdI)
{
  /*
  Form for Correlation function
  */
  if (strcmp(StdI->outputmode, "non") == 0
    || strcmp(StdI->outputmode, "none") == 0
    || strcmp(StdI->outputmode, "off") == 0) {
    StdI->ioutputmode = 0;
    fprintf(stdout, "      ioutputmode = %-10d\n", StdI->ioutputmode);
  }
  else if (strcmp(StdI->outputmode, "cor") == 0
    || strcmp(StdI->outputmode, "corr") == 0
    || strcmp(StdI->outputmode, "correlation") == 0) {
    StdI->ioutputmode = 1;
    fprintf(stdout, "      ioutputmode = %-10d\n", StdI->ioutputmode);
  }
  else if (strcmp(StdI->outputmode, "****") == 0) {
    StdI->ioutputmode = 1;
    fprintf(stdout, "      ioutputmode = %-10d  ######  DEFAULT VALUE IS USED  ######\n", StdI->ioutputmode);
  }
  else if (strcmp(StdI->outputmode, "raw") == 0
    || strcmp(StdI->outputmode, "all") == 0
    || strcmp(StdI->outputmode, "full") == 0) {
    StdI->ioutputmode = 2;
    fprintf(stdout, "      ioutputmode = %-10d\n", StdI->ioutputmode);
  }
  else{
    fprintf(stdout, "\n ERROR ! Unsupported OutPutMode : %s\n", StdI->outputmode);
    StdFace_exit(-1);
  }
}/*static void CheckOutputMode*/
/**
 *
 * Summary numerical parameter check the combination of
 * the number of sites, total spin, the number of electrons
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void CheckModPara(struct StdIntList *StdI)
{
#if defined(_mVMC)
  int NSym2;
#endif

  if (strcmp(StdI->CDataFileHead, "****") == 0) {
    strcpy(StdI->CDataFileHead, "zvo\0");
    fprintf(stdout, "    CDataFileHead = %-12s######  DEFAULT VALUE IS USED  ######\n", StdI->CDataFileHead);
  }
  else fprintf(stdout, "    CDataFileHead = %-s\n", StdI->CDataFileHead);

  /**/
#if defined(_HPhi)
  StdFace_PrintVal_i("Lanczos_max", &StdI->Lanczos_max, 2000);
  StdFace_PrintVal_i("initial_iv", &StdI->initial_iv, -1);
  /*StdFace_PrintVal_i("nvec", &StdI->nvec, 1);*/
  StdFace_PrintVal_i("exct", &StdI->exct, 1);
  StdFace_PrintVal_i("LanczosEps", &StdI->LanczosEps, 14);
  StdFace_PrintVal_i("LanczosTarget", &StdI->LanczosTarget, 2);
  if(StdI->LanczosTarget < StdI->exct) StdI->LanczosTarget = StdI->exct;
  StdFace_PrintVal_i("NumAve", &StdI->NumAve, 5);
  StdFace_PrintVal_i("ExpecInterval", &StdI->ExpecInterval, 20);
  StdFace_PrintVal_i("NOmega", &StdI->Nomega, 200);
  StdFace_PrintVal_d("OmegaMax", &StdI->OmegaMax, StdI->LargeValue*StdI->nsite);
  StdFace_PrintVal_d("OmegaMin", &StdI->OmegaMin, -StdI->LargeValue*StdI->nsite);
  StdFace_PrintVal_d("OmegaIm", &StdI->OmegaIm, 0.01* (int)StdI->LargeValue);
#elif defined(_mVMC)
  if (strcmp(StdI->CParaFileHead, "****") == 0) {
    strcpy(StdI->CParaFileHead, "zqp\0");
    fprintf(stdout, "    CParaFileHead = %-12s######  DEFAULT VALUE IS USED  ######\n", StdI->CParaFileHead);
  }
  else fprintf(stdout, "    CParaFileHead = %-s\n", StdI->CParaFileHead);
  
  StdFace_PrintVal_i("NVMCCalMode", &StdI->NVMCCalMode, 0);
  StdFace_PrintVal_i("NLanczosMode", &StdI->NLanczosMode, 0);
  StdFace_PrintVal_i("NDataIdxStart", &StdI->NDataIdxStart, 1);

  if (StdI->NVMCCalMode == 0) StdFace_NotUsed_i("NDataQtySmp", StdI->NDataQtySmp);
  /*else*/StdFace_PrintVal_i("NDataQtySmp", &StdI->NDataQtySmp, 1);

  StdFace_PrintVal_i("NSPGaussLeg", &StdI->NSPGaussLeg, 8);
  if (StdI->AntiPeriod[0] == 1 || StdI->AntiPeriod[1] == 1) NSym2 = -StdI->NSym;
  else NSym2 = StdI->NSym;
  StdFace_PrintVal_i("NMPTrans", &StdI->NMPTrans, NSym2);

  StdFace_PrintVal_i("NSROptItrStep", &StdI->NSROptItrStep, 1000);
  
  if (StdI->NVMCCalMode == 1) StdFace_NotUsed_i("NSROptItrSmp", StdI->NSROptItrSmp);
  /*else*/ StdFace_PrintVal_i("NSROptItrSmp", &StdI->NSROptItrSmp, StdI->NSROptItrStep/10);

  StdFace_PrintVal_i("NVMCWarmUp", &StdI->NVMCWarmUp, 10);
  StdFace_PrintVal_i("NVMCInterval", &StdI->NVMCInterval, 1);
  StdFace_PrintVal_i("NVMCSample", &StdI->NVMCSample, 1000);

  if (strcmp(StdI->model, "hubbard") == 0) StdI->NExUpdatePath = 0;
  else if (strcmp(StdI->model, "spin") == 0) StdI->NExUpdatePath = 2;
  else if (strcmp(StdI->model, "kondo") == 0) StdI->NExUpdatePath = 1;
  fprintf(stdout, "  %15s = %-10d\n", "NExUpdatePath", StdI->NExUpdatePath);

  StdFace_PrintVal_i("RndSeed", &StdI->RndSeed, 123456789);
  StdFace_PrintVal_i("NSplitSize", &StdI->NSplitSize, 1);
  StdFace_PrintVal_i("NStore", &StdI->NStore, 0);

  StdFace_PrintVal_d("DSROptRedCut", &StdI->DSROptRedCut, 0.001);
  StdFace_PrintVal_d("DSROptStaDel", &StdI->DSROptStaDel, 0.02);
  StdFace_PrintVal_d("DSROptStepDt", &StdI->DSROptStepDt, 0.02);
#endif
  /**/
  if (strcmp(StdI->model, "hubbard") == 0){
    if (StdI->lGC == 0) StdFace_RequiredVal_i("nelec", StdI->nelec);
    else {
      StdFace_NotUsed_i("nelec", StdI->nelec);
      StdFace_NotUsed_i("2Sz", StdI->Sz2);
    }
#if defined(_mVMC)
    if (StdI->lGC == 0) {
      StdFace_PrintVal_i("2Sz", &StdI->Sz2, 0);
      if (StdI->nelec % 2 != 0) {
        printf("\nERROR ! nelec should be an even number !\n\n");
        StdFace_exit(-1);
      }
      else {
        StdI->nelec = StdI->nelec / 2;
      }
    }
#endif
  }
  else if (strcmp(StdI->model, "spin") == 0) {
    StdFace_NotUsed_i("nelec", StdI->nelec);
    if (StdI->lGC == 0) StdFace_RequiredVal_i("2Sz", StdI->Sz2);
    else StdFace_NotUsed_i("2Sz", StdI->Sz2);
#if defined(_mVMC)
    StdI->nelec = StdI->nsite / 2;
#endif
  }
  else if (strcmp(StdI->model, "kondo") == 0) {
    if (StdI->lGC == 0) StdFace_RequiredVal_i("nelec", StdI->nelec);
    else {
      StdFace_NotUsed_i("nelec", StdI->nelec);
      StdFace_NotUsed_i("2Sz", StdI->Sz2);
    }
#if defined(_mVMC)
    if (StdI->lGC == 0) {
      StdFace_PrintVal_i("2Sz", &StdI->Sz2, 0);
      if ((StdI->nelec + StdI->nsite / 2) % 2 != 0) {
        printf("\nERROR ! nelec should be an even number !\n\n");
        StdFace_exit(-1);
      }
      else {
        StdI->nelec = (StdI->nelec + StdI->nsite / 2) / 2;
      }
    }
#endif
  }
}/*static void CheckModPara*/
/**
 *
 * Output .def file for Specific interaction
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void PrintInteractions(struct StdIntList *StdI)
{
  FILE *fp;
  int nintr0, kintr, jintr;
  /*
   Coulomb INTRA
  */
  nintr0 = 0;
  for (kintr = 0; kintr < StdI->NCintra; kintr++) {
    if (fabs(StdI->Cintra[kintr]) > 0.000001) nintr0 = nintr0 + 1;
  }
  if (nintr0 == 0 || StdI->lBoost == 1) StdI->LCintra = 0;
  else StdI->LCintra = 1;

  if (StdI->LCintra == 1) {
    fp = fopen("coulombintra.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NCoulombIntra %10d\n", nintr0);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "================== CoulombIntra ================\n");
    fprintf(fp, "=============================================\n");
    for (kintr = 0; kintr < StdI->NCintra; kintr++) {
      if (fabs(StdI->Cintra[kintr]) > 0.000001)
        fprintf(fp, "%5d %25.15f\n",
          StdI->CintraIndx[kintr][0], StdI->Cintra[kintr]);
    }
    fflush(fp);
    fclose(fp);
    fprintf(stdout, "    coulombintra.def is written.\n");
  }/*if (StdI->LCintra == 1)*/
  /*
  Coulomb INTER
  */
  nintr0 = 0;
  for (kintr = 0; kintr < StdI->NCinter; kintr++) {
    if (fabs(StdI->Cinter[kintr]) > 0.000001) nintr0 = nintr0 + 1;
  }
  if (nintr0 == 0 || StdI->lBoost == 1) StdI->LCinter = 0;
  else StdI->LCinter = 1;

  if (StdI->LCinter == 1) {
    fp = fopen("coulombinter.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NCoulombInter %10d\n", nintr0);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "================== CoulombInter ================\n");
    fprintf(fp, "=============================================\n");
    for (kintr = 0; kintr < StdI->NCinter; kintr++) {
      if (fabs(StdI->Cinter[kintr]) > 0.000001)
        fprintf(fp, "%5d %5d %25.15f\n",
          StdI->CinterIndx[kintr][0], StdI->CinterIndx[kintr][1], StdI->Cinter[kintr]);
    }
    fflush(fp);
    fclose(fp);
    fprintf(stdout, "    coulombinter.def is written.\n");
  }/*if (StdI->LCinter == 1)*/
  /*
  Hund
  */
  nintr0 = 0;
  for (kintr = 0; kintr < StdI->NHund; kintr++) {
    if (fabs(StdI->Hund[kintr]) > 0.000001) nintr0 = nintr0 + 1;
  }
  if (nintr0 == 0 || StdI->lBoost == 1) StdI->LHund = 0;
  else StdI->LHund = 1;

  if (StdI->LHund == 1) {
    fp = fopen("hund.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NHund %10d\n", nintr0);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "=============== Hund coupling ===============\n");
    fprintf(fp, "=============================================\n");
    for (kintr = 0; kintr < StdI->NHund; kintr++) {
      if (fabs(StdI->Hund[kintr]) > 0.000001)
        fprintf(fp, "%5d %5d %25.15f\n",
          StdI->HundIndx[kintr][0], StdI->HundIndx[kintr][1], StdI->Hund[kintr]);
    }
    fflush(fp);
    fclose(fp);
    fprintf(stdout, "    hund.def is written.\n");
  }/*if (StdI->LHund == 1)*/
  /*
  Exchange
  */
  nintr0 = 0;
  for (kintr = 0; kintr < StdI->NEx; kintr++) {
    if (fabs(StdI->Ex[kintr]) > 0.000001) nintr0 = nintr0 + 1;
  }
  if (nintr0 == 0 || StdI->lBoost == 1) StdI->LEx = 0;
  else StdI->LEx = 1;

  if (StdI->LEx == 1) {
    fp = fopen("exchange.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NExchange %10d\n", nintr0);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "====== ExchangeCoupling coupling ============\n");
    fprintf(fp, "=============================================\n");
    for (kintr = 0; kintr < StdI->NEx; kintr++) {
      if (fabs(StdI->Ex[kintr]) > 0.000001)
        fprintf(fp, "%5d %5d %25.15f\n",
          StdI->ExIndx[kintr][0], StdI->ExIndx[kintr][1], StdI->Ex[kintr]);
    }
    fflush(fp);
    fclose(fp);
    fprintf(stdout, "    exchange.def is written.\n");
  }
  /*
    PairLift
  */
  nintr0 = 0;
  for (kintr = 0; kintr < StdI->NPairLift; kintr++) {
    if (fabs(StdI->PairLift[kintr]) > 0.000001) nintr0 = nintr0 + 1;
  }
  if (nintr0 == 0 || StdI->lBoost == 1) StdI->LPairLift = 0;
  else StdI->LPairLift = 1;

  if (StdI->LPairLift == 1) {
    fp = fopen("pairlift.def", "w");
    fprintf(fp, "=============================================\n");
    fprintf(fp, "NExchange %10d\n", nintr0);
    fprintf(fp, "=============================================\n");
    fprintf(fp, "====== ExchangeCoupling coupling ============\n");
    fprintf(fp, "=============================================\n");
    for (kintr = 0; kintr < StdI->NPairLift; kintr++) {
      if (fabs(StdI->PairLift[kintr]) > 0.000001)
        fprintf(fp, "%5d %5d %25.15f\n",
          StdI->PLIndx[kintr][0], StdI->PLIndx[kintr][1], StdI->PairLift[kintr]);
    }
    fflush(fp);
    fclose(fp);
    fprintf(stdout, "    exchange.def is written.\n");
  }
  /*
   InterAll
  */
  for (jintr = 0; jintr < StdI->nintr; jintr++) {
    for (kintr = jintr + 1; kintr < StdI->nintr; kintr++) {
      if (
        (StdI->intrindx[jintr][0] == StdI->intrindx[kintr][0]
          && StdI->intrindx[jintr][1] == StdI->intrindx[kintr][1]
          && StdI->intrindx[jintr][2] == StdI->intrindx[kintr][2]
          && StdI->intrindx[jintr][3] == StdI->intrindx[kintr][3]
          && StdI->intrindx[jintr][4] == StdI->intrindx[kintr][4]
          && StdI->intrindx[jintr][5] == StdI->intrindx[kintr][5]
          && StdI->intrindx[jintr][6] == StdI->intrindx[kintr][6]
          && StdI->intrindx[jintr][7] == StdI->intrindx[kintr][7])
        ||
        (StdI->intrindx[jintr][0] == StdI->intrindx[kintr][4]
          && StdI->intrindx[jintr][1] == StdI->intrindx[kintr][5]
          && StdI->intrindx[jintr][2] == StdI->intrindx[kintr][6]
          && StdI->intrindx[jintr][3] == StdI->intrindx[kintr][7]
          && StdI->intrindx[jintr][4] == StdI->intrindx[kintr][0]
          && StdI->intrindx[jintr][5] == StdI->intrindx[kintr][1]
          && StdI->intrindx[jintr][6] == StdI->intrindx[kintr][2]
          && StdI->intrindx[jintr][7] == StdI->intrindx[kintr][3])
        ) {
        StdI->intr[jintr] = StdI->intr[jintr] + StdI->intr[kintr];
        StdI->intr[kintr] = 0.0;
      }
      else if (
        (StdI->intrindx[jintr][0] == StdI->intrindx[kintr][4]
          && StdI->intrindx[jintr][1] == StdI->intrindx[kintr][5]
          && StdI->intrindx[jintr][2] == StdI->intrindx[kintr][2]
          && StdI->intrindx[jintr][3] == StdI->intrindx[kintr][3]
          && StdI->intrindx[jintr][4] == StdI->intrindx[kintr][0]
          && StdI->intrindx[jintr][5] == StdI->intrindx[kintr][1]
          && StdI->intrindx[jintr][6] == StdI->intrindx[kintr][6]
          && StdI->intrindx[jintr][7] == StdI->intrindx[kintr][7])
        ||
        (StdI->intrindx[jintr][0] == StdI->intrindx[kintr][0]
          && StdI->intrindx[jintr][1] == StdI->intrindx[kintr][1]
          && StdI->intrindx[jintr][2] == StdI->intrindx[kintr][6]
          && StdI->intrindx[jintr][3] == StdI->intrindx[kintr][7]
          && StdI->intrindx[jintr][4] == StdI->intrindx[kintr][4]
          && StdI->intrindx[jintr][5] == StdI->intrindx[kintr][5]
          && StdI->intrindx[jintr][6] == StdI->intrindx[kintr][2]
          && StdI->intrindx[jintr][7] == StdI->intrindx[kintr][3])
        ) {
        StdI->intr[jintr] = StdI->intr[jintr] - StdI->intr[kintr];
        StdI->intr[kintr] = 0.0;
      }
    }/*for (kintr = jintr + 1; kintr < StdI->nintr; kintr++)*/
  }/*for (jintr = 0; jintr < StdI->nintr; jintr++)*/

  for (jintr = 0; jintr < StdI->nintr; jintr++) {
    for (kintr = jintr + 1; kintr < StdI->nintr; kintr++) {
      if (StdI->intrindx[jintr][6] == StdI->intrindx[kintr][4]
        && StdI->intrindx[jintr][7] == StdI->intrindx[kintr][5]
        && StdI->intrindx[jintr][4] == StdI->intrindx[kintr][6]
        && StdI->intrindx[jintr][5] == StdI->intrindx[kintr][7]
        && StdI->intrindx[jintr][2] == StdI->intrindx[kintr][0]
        && StdI->intrindx[jintr][3] == StdI->intrindx[kintr][1]
        && StdI->intrindx[jintr][0] == StdI->intrindx[kintr][2]
        && StdI->intrindx[jintr][1] == StdI->intrindx[kintr][3]
        ) {
        StdI->intrindx[kintr][0] = StdI->intrindx[jintr][6];
        StdI->intrindx[kintr][1] = StdI->intrindx[jintr][7];
        StdI->intrindx[kintr][2] = StdI->intrindx[jintr][4];
        StdI->intrindx[kintr][3] = StdI->intrindx[jintr][5];
        StdI->intrindx[kintr][4] = StdI->intrindx[jintr][2];
        StdI->intrindx[kintr][5] = StdI->intrindx[jintr][3];
        StdI->intrindx[kintr][6] = StdI->intrindx[jintr][0];
        StdI->intrindx[kintr][7] = StdI->intrindx[jintr][1];
      }
      else if (
        (StdI->intrindx[jintr][6] == StdI->intrindx[kintr][4]
          && StdI->intrindx[jintr][7] == StdI->intrindx[kintr][5]
          && StdI->intrindx[jintr][4] == StdI->intrindx[kintr][2]
          && StdI->intrindx[jintr][5] == StdI->intrindx[kintr][3]
          && StdI->intrindx[jintr][2] == StdI->intrindx[kintr][0]
          && StdI->intrindx[jintr][3] == StdI->intrindx[kintr][1]
          && StdI->intrindx[jintr][0] == StdI->intrindx[kintr][6]
          && StdI->intrindx[jintr][1] == StdI->intrindx[kintr][7])
        ||
        (StdI->intrindx[jintr][6] == StdI->intrindx[kintr][0]
          && StdI->intrindx[jintr][7] == StdI->intrindx[kintr][1]
          && StdI->intrindx[jintr][4] == StdI->intrindx[kintr][6]
          && StdI->intrindx[jintr][5] == StdI->intrindx[kintr][7]
          && StdI->intrindx[jintr][2] == StdI->intrindx[kintr][4]
          && StdI->intrindx[jintr][3] == StdI->intrindx[kintr][5]
          && StdI->intrindx[jintr][0] == StdI->intrindx[kintr][2]
          && StdI->intrindx[jintr][1] == StdI->intrindx[kintr][3])
        ) {
        StdI->intrindx[kintr][0] = StdI->intrindx[jintr][6];
        StdI->intrindx[kintr][1] = StdI->intrindx[jintr][7];
        StdI->intrindx[kintr][2] = StdI->intrindx[jintr][4];
        StdI->intrindx[kintr][3] = StdI->intrindx[jintr][5];
        StdI->intrindx[kintr][4] = StdI->intrindx[jintr][2];
        StdI->intrindx[kintr][5] = StdI->intrindx[jintr][3];
        StdI->intrindx[kintr][6] = StdI->intrindx[jintr][0];
        StdI->intrindx[kintr][7] = StdI->intrindx[jintr][1];

        StdI->intr[kintr] = -StdI->intr[kintr];
      }
    }/*for (kintr = jintr + 1; kintr < StdI->nintr; kintr++)*/
  }/*for (jintr = 0; jintr < StdI->nintr; jintr++)*/

  for (jintr = 0; jintr < StdI->nintr; jintr++) {

    if (
      (StdI->intrindx[jintr][0] == StdI->intrindx[jintr][4]
        && StdI->intrindx[jintr][1] == StdI->intrindx[jintr][5]) ||
        (StdI->intrindx[jintr][2] == StdI->intrindx[jintr][6]
          && StdI->intrindx[jintr][3] == StdI->intrindx[jintr][7])
      ) {

      if (!(
        (StdI->intrindx[jintr][0] == StdI->intrindx[jintr][2]
          && StdI->intrindx[jintr][1] == StdI->intrindx[jintr][3])
        ||
        (StdI->intrindx[jintr][0] == StdI->intrindx[jintr][6]
          && StdI->intrindx[jintr][1] == StdI->intrindx[jintr][7])
        ||
        (StdI->intrindx[jintr][4] == StdI->intrindx[jintr][2]
          && StdI->intrindx[jintr][5] == StdI->intrindx[jintr][3])
        ||
        (StdI->intrindx[jintr][4] == StdI->intrindx[jintr][6]
          && StdI->intrindx[jintr][5] == StdI->intrindx[jintr][7])
        ))
        StdI->intr[jintr] = 0.0;
    }
  }/*for (jintr = 0; jintr < StdI->nintr; jintr++)*/

  nintr0 = 0;
  for (kintr = 0; kintr < StdI->nintr; kintr++) {
    if (cabs(StdI->intr[kintr]) > 0.000001) nintr0 = nintr0 + 1;
  }
  if (nintr0 == 0 || StdI->lBoost == 1) StdI->Lintr = 0;
  else StdI->Lintr = 1;

  if (StdI->Lintr == 1) {
    fp = fopen("interall.def", "w");
    fprintf(fp, "====================== \n");
    fprintf(fp, "NInterAll %7d  \n", nintr0);
    fprintf(fp, "====================== \n");
    fprintf(fp, "========zInterAll===== \n");
    fprintf(fp, "====================== \n");

    if (StdI->lBoost == 0) {
      nintr0 = 0;
      for (kintr = 0; kintr < StdI->nintr; kintr++) {
        if (cabs(StdI->intr[kintr]) > 0.000001)
          fprintf(fp, "%5d %5d %5d %5d %5d %5d %5d %5d %25.15f  %25.15f\n",
            StdI->intrindx[kintr][0], StdI->intrindx[kintr][1],
            StdI->intrindx[kintr][2], StdI->intrindx[kintr][3],
            StdI->intrindx[kintr][4], StdI->intrindx[kintr][5],
            StdI->intrindx[kintr][6], StdI->intrindx[kintr][7],
            creal(StdI->intr[kintr]), cimag(StdI->intr[kintr]));
      }/*for (kintr = 0; kintr < StdI->nintr; kintr++)*/
    }/* if (StdI->lBoost == 0)*/

    fflush(fp);
    fclose(fp);
    fprintf(stdout, "    interall.def is written.\n");
  }
}/*static void PrintInteractions*/
/**
*
* Main routine for the standard mode
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_main(char *fname  /**< [in] Input file name for the standard mode */) {

  struct StdIntList StdI;

  FILE *fp;
  int ktrans, kintr;
  char ctmpline[256];
  char *keyword, *value;

  fprintf(stdout, "\n######  Standard Intarface Mode STARTS  ######\n");
  if ((fp = fopen(fname, "r")) == NULL) {
    fprintf(stdout, "\n  ERROR !  Cannot open input file %s !\n\n", fname);
    StdFace_exit(-1);
  }
  else {
    fprintf(stdout, "\n  Open Standard-Mode Inputfile %s \n\n", fname);
  }

  StdFace_ResetVals(&StdI);

  while (fgets(ctmpline, 256, fp) != NULL) {

    TrimSpaceQuote(ctmpline);
    if (strncmp(ctmpline, "//", 2) == 0) {
      fprintf(stdout, "  Skipping a line.\n");
      continue;
    }
    else if (ctmpline[0] == '\0') {
      fprintf(stdout, "  Skipping a line.\n");
      continue;
    }
    keyword = strtok(ctmpline, "=");
    value = strtok(NULL, "=");
    if (value == NULL) {
      fprintf(stdout, "\n  ERROR !  \"=\" is NOT found !\n\n");
      StdFace_exit(-1);
    }
    Text2Lower(keyword);
    fprintf(stdout, "  KEYWORD : %-20s | VALUE : %s \n", keyword, value);

    if (strcmp(keyword, "a") == 0) StoreWithCheckDup_d(keyword, value, &StdI.a);
    else if (strcmp(keyword, "a0h") == 0) StoreWithCheckDup_i(keyword, value, &StdI.box[0][2]);
    else if (strcmp(keyword, "a0l") == 0) StoreWithCheckDup_i(keyword, value, &StdI.box[0][1]);
    else if (strcmp(keyword, "a0w") == 0) StoreWithCheckDup_i(keyword, value, &StdI.box[0][0]);
    else if (strcmp(keyword, "a1h") == 0) StoreWithCheckDup_i(keyword, value, &StdI.box[1][2]);
    else if (strcmp(keyword, "a1l") == 0) StoreWithCheckDup_i(keyword, value, &StdI.box[1][1]);
    else if (strcmp(keyword, "a1w") == 0) StoreWithCheckDup_i(keyword, value, &StdI.box[1][0]);
    else if (strcmp(keyword, "a2h") == 0) StoreWithCheckDup_i(keyword, value, &StdI.box[2][2]);
    else if (strcmp(keyword, "a2l") == 0) StoreWithCheckDup_i(keyword, value, &StdI.box[2][1]);
    else if (strcmp(keyword, "a2w") == 0) StoreWithCheckDup_i(keyword, value, &StdI.box[2][0]);
    else if (strcmp(keyword, "d") == 0) StoreWithCheckDup_d(keyword, value, &StdI.D[2][2]);
    else if (strcmp(keyword, "gamma") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Gamma);
    else if (strcmp(keyword, "h") == 0) StoreWithCheckDup_d(keyword, value, &StdI.h);
    else if (strcmp(keyword, "height") == 0) StoreWithCheckDup_i(keyword, value, &StdI.Height);
    else if (strcmp(keyword, "hlength") == 0) StoreWithCheckDup_d(keyword, value, &StdI.length[2]);
    else if (strcmp(keyword, "hx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.direct[2][0]);
    else if (strcmp(keyword, "hy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.direct[2][1]);
    else if (strcmp(keyword, "hz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.direct[2][2]);
    else if (strcmp(keyword, "j") == 0) StoreWithCheckDup_d(keyword, value, &StdI.JAll);
    else if (strcmp(keyword, "jx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J[0][0]);
    else if (strcmp(keyword, "jxy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J[0][1]);
    else if (strcmp(keyword, "jxz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J[0][2]);
    else if (strcmp(keyword, "jy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J[1][1]);
    else if (strcmp(keyword, "jyx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J[1][0]);
    else if (strcmp(keyword, "jyz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J[1][2]);
    else if (strcmp(keyword, "jz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J[2][2]);
    else if (strcmp(keyword, "jzx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J[2][0]);
    else if (strcmp(keyword, "jzy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J[2][1]);
    else if (strcmp(keyword, "j0") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0All);
    else if (strcmp(keyword, "j0x") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0[0][0]);
    else if (strcmp(keyword, "j0xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0[0][1]);
    else if (strcmp(keyword, "j0xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0[0][2]);
    else if (strcmp(keyword, "j0y") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0[1][1]);
    else if (strcmp(keyword, "j0yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0[1][0]);
    else if (strcmp(keyword, "j0yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0[1][2]);
    else if (strcmp(keyword, "j0z") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0[2][2]);
    else if (strcmp(keyword, "j0zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0[2][0]);
    else if (strcmp(keyword, "j0zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0[2][1]);
    else if (strcmp(keyword, "j0'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0pAll);
    else if (strcmp(keyword, "j0'x") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0p[0][0]);
    else if (strcmp(keyword, "j0'xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0p[0][1]);
    else if (strcmp(keyword, "j0'xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0p[0][2]);
    else if (strcmp(keyword, "j0'y") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0p[1][1]);
    else if (strcmp(keyword, "j0'yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0p[1][0]);
    else if (strcmp(keyword, "j0'yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0p[1][2]);
    else if (strcmp(keyword, "j0'z") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0p[2][2]);
    else if (strcmp(keyword, "j0'zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0p[2][0]);
    else if (strcmp(keyword, "j0'zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0p[2][1]);
    else if (strcmp(keyword, "j1") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1All);
    else if (strcmp(keyword, "j1x") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1[0][0]);
    else if (strcmp(keyword, "j1xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1[0][1]);
    else if (strcmp(keyword, "j1xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1[0][2]);
    else if (strcmp(keyword, "j1y") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1[1][1]);
    else if (strcmp(keyword, "j1yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1[1][0]);
    else if (strcmp(keyword, "j1yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1[1][2]);
    else if (strcmp(keyword, "j1z") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1[2][2]);
    else if (strcmp(keyword, "j1zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1[2][0]);
    else if (strcmp(keyword, "j1zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1[2][1]);
    else if (strcmp(keyword, "j1'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1pAll);
    else if (strcmp(keyword, "j1'x") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1p[0][0]);
    else if (strcmp(keyword, "j1'xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1p[0][1]);
    else if (strcmp(keyword, "j1'xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1p[0][2]);
    else if (strcmp(keyword, "j1'y") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1p[1][1]);
    else if (strcmp(keyword, "j1'yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1p[1][0]);
    else if (strcmp(keyword, "j1'yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1p[1][2]);
    else if (strcmp(keyword, "j1'z") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1p[2][2]);
    else if (strcmp(keyword, "j1'zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1p[2][0]);
    else if (strcmp(keyword, "j1'zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1p[2][1]);
    else if (strcmp(keyword, "j2") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2All);
    else if (strcmp(keyword, "j2x") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2[0][0]);
    else if (strcmp(keyword, "j2xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2[0][1]);
    else if (strcmp(keyword, "j2xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2[0][2]);
    else if (strcmp(keyword, "j2y") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2[1][1]);
    else if (strcmp(keyword, "j2yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2[1][0]);
    else if (strcmp(keyword, "j2yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2[1][2]);
    else if (strcmp(keyword, "j2z") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2[2][2]);
    else if (strcmp(keyword, "j2zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2[2][0]);
    else if (strcmp(keyword, "j2zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2[2][1]);
    else if (strcmp(keyword, "j2'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2pAll);
    else if (strcmp(keyword, "j2'x") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2p[0][0]);
    else if (strcmp(keyword, "j2'xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2p[0][1]);
    else if (strcmp(keyword, "j2'xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2p[0][2]);
    else if (strcmp(keyword, "j2'y") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2p[1][1]);
    else if (strcmp(keyword, "j2'yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2p[1][0]);
    else if (strcmp(keyword, "j2'yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2p[1][2]);
    else if (strcmp(keyword, "j2'z") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2p[2][2]);
    else if (strcmp(keyword, "j2'zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2p[2][0]);
    else if (strcmp(keyword, "j2'zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2p[2][1]);
    else if (strcmp(keyword, "j'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.JpAll);
    else if (strcmp(keyword, "j'x") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jp[0][0]);
    else if (strcmp(keyword, "j'xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jp[0][1]);
    else if (strcmp(keyword, "j'xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jp[0][2]);
    else if (strcmp(keyword, "j'y") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jp[1][1]);
    else if (strcmp(keyword, "j'yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jp[1][0]);
    else if (strcmp(keyword, "j'yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jp[1][2]);
    else if (strcmp(keyword, "j'z") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jp[2][2]);
    else if (strcmp(keyword, "j'zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jp[2][0]);
    else if (strcmp(keyword, "j'zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jp[2][1]);
    else if (strcmp(keyword, "j''") == 0) StoreWithCheckDup_d(keyword, value, &StdI.JppAll);
    else if (strcmp(keyword, "j''x") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jpp[0][0]);
    else if (strcmp(keyword, "j''xy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jpp[0][1]);
    else if (strcmp(keyword, "j''xz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jpp[0][2]);
    else if (strcmp(keyword, "j''y") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jpp[1][1]);
    else if (strcmp(keyword, "j''yx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jpp[1][0]);
    else if (strcmp(keyword, "j''yz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jpp[1][2]);
    else if (strcmp(keyword, "j''z") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jpp[2][2]);
    else if (strcmp(keyword, "j''zx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jpp[2][0]);
    else if (strcmp(keyword, "j''zy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jpp[2][1]);
    else if (strcmp(keyword, "k") == 0) StoreWithCheckDup_d(keyword, value, &StdI.K);
    else if (strcmp(keyword, "l") == 0) StoreWithCheckDup_i(keyword, value, &StdI.L);
    else if (strcmp(keyword, "lattice") == 0) StoreWithCheckDup_sl(keyword, value, StdI.lattice);
    else if (strcmp(keyword, "llength") == 0) StoreWithCheckDup_d(keyword, value, &StdI.length[1]);
    else if (strcmp(keyword, "lx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.direct[1][0]);
    else if (strcmp(keyword, "ly") == 0) StoreWithCheckDup_d(keyword, value, &StdI.direct[1][1]);
    else if (strcmp(keyword, "lz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.direct[1][2]);
    else if (strcmp(keyword, "model") == 0) StoreWithCheckDup_sl(keyword, value, StdI.model);
    else if (strcmp(keyword, "mu") == 0) StoreWithCheckDup_d(keyword, value, &StdI.mu);
    else if (strcmp(keyword, "nelec") == 0) StoreWithCheckDup_i(keyword, value, &StdI.nelec);
    else if (strcmp(keyword, "outputmode") == 0) StoreWithCheckDup_sl(keyword, value, StdI.outputmode);
    else if (strcmp(keyword, "phase0") == 0) StoreWithCheckDup_d(keyword, value, &StdI.phase[0]);
    else if (strcmp(keyword, "phase1") == 0) StoreWithCheckDup_d(keyword, value, &StdI.phase[1]);
    else if (strcmp(keyword, "phase2") == 0) StoreWithCheckDup_d(keyword, value, &StdI.phase[2]);
    else if (strcmp(keyword, "t") == 0) StoreWithCheckDup_c(keyword, value, &StdI.t);
    else if (strcmp(keyword, "t0") == 0) StoreWithCheckDup_c(keyword, value, &StdI.t0);
    else if (strcmp(keyword, "t0'") == 0) StoreWithCheckDup_c(keyword, value, &StdI.t0p);
    else if (strcmp(keyword, "t1") == 0) StoreWithCheckDup_c(keyword, value, &StdI.t1);
    else if (strcmp(keyword, "t1'") == 0) StoreWithCheckDup_c(keyword, value, &StdI.t1p);
    else if (strcmp(keyword, "t2") == 0) StoreWithCheckDup_c(keyword, value, &StdI.t2);
    else if (strcmp(keyword, "t2'") == 0) StoreWithCheckDup_c(keyword, value, &StdI.t2p);
    else if (strcmp(keyword, "t'") == 0) StoreWithCheckDup_c(keyword, value, &StdI.tp);
    else if (strcmp(keyword, "t''") == 0) StoreWithCheckDup_c(keyword, value, &StdI.tpp);
    else if (strcmp(keyword, "u") == 0) StoreWithCheckDup_d(keyword, value, &StdI.U);
    else if (strcmp(keyword, "v") == 0) StoreWithCheckDup_d(keyword, value, &StdI.V);
    else if (strcmp(keyword, "v0") == 0) StoreWithCheckDup_d(keyword, value, &StdI.V0);
    else if (strcmp(keyword, "v0'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.V0p);
    else if (strcmp(keyword, "v1") == 0) StoreWithCheckDup_d(keyword, value, &StdI.V1);
    else if (strcmp(keyword, "v1'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.V1p);
    else if (strcmp(keyword, "v2") == 0) StoreWithCheckDup_d(keyword, value, &StdI.V2);
    else if (strcmp(keyword, "v2p") == 0) StoreWithCheckDup_d(keyword, value, &StdI.V2);
    else if (strcmp(keyword, "v'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Vp);
    else if (strcmp(keyword, "v''") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Vpp);
    else if (strcmp(keyword, "w") == 0) StoreWithCheckDup_i(keyword, value, &StdI.W);
    else if (strcmp(keyword, "wlength") == 0) StoreWithCheckDup_d(keyword, value, &StdI.length[0]);
    else if (strcmp(keyword, "wx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.direct[0][0]);
    else if (strcmp(keyword, "wy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.direct[0][1]);
    else if (strcmp(keyword, "wz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.direct[0][2]);
    else if (strcmp(keyword, "w90_cutoff") == 0) StoreWithCheckDup_d(keyword, value, &StdI.W90_cutoff);
    else if (strcmp(keyword, "w90_geom") == 0) StoreWithCheckDup_s(keyword, value, StdI.W90_geom);
    else if (strcmp(keyword, "w90_hr") == 0) StoreWithCheckDup_s(keyword, value, StdI.W90_hr);
    else if (strcmp(keyword, "2sz") == 0) StoreWithCheckDup_i(keyword, value, &StdI.Sz2);

#if defined(_HPhi)
    else if (strcmp(keyword, "calcspec") == 0) StoreWithCheckDup_sl(keyword, value, StdI.CalcSpec);
    else if (strcmp(keyword, "exct") == 0) StoreWithCheckDup_i(keyword, value, &StdI.exct);
    else if (strcmp(keyword, "eigenvecio") == 0) StoreWithCheckDup_sl(keyword, value, StdI.EigenVecIO);
    else if (strcmp(keyword, "expecinterval") == 0) StoreWithCheckDup_i(keyword, value, &StdI.ExpecInterval);
    else if (strcmp(keyword, "cdatafilehead") == 0) StoreWithCheckDup_s(keyword, value, StdI.CDataFileHead);
    else if (strcmp(keyword, "flgtemp") == 0) StoreWithCheckDup_i(keyword, value, &StdI.FlgTemp);
    else if (strcmp(keyword, "initialvectype") == 0) StoreWithCheckDup_sl(keyword, value, StdI.InitialVecType);
    else if (strcmp(keyword, "initial_iv") == 0) StoreWithCheckDup_i(keyword, value, &StdI.initial_iv);
    else if (strcmp(keyword, "lanczoseps") == 0) StoreWithCheckDup_i(keyword, value, &StdI.LanczosEps);
    else if (strcmp(keyword, "lanczostarget") == 0) StoreWithCheckDup_i(keyword, value, &StdI.LanczosTarget);
    else if (strcmp(keyword, "lanczos_max") == 0) StoreWithCheckDup_i(keyword, value, &StdI.Lanczos_max);
    else if (strcmp(keyword, "largevalue") == 0) StoreWithCheckDup_d(keyword, value, &StdI.LargeValue);
    else if (strcmp(keyword, "method") == 0) StoreWithCheckDup_sl(keyword, value, StdI.method);
    else if (strcmp(keyword, "nomega") == 0) StoreWithCheckDup_i(keyword, value, &StdI.Nomega);
    else if (strcmp(keyword, "numave") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NumAve);
    else if (strcmp(keyword, "nvec") == 0) StoreWithCheckDup_i(keyword, value, &StdI.nvec);
    else if (strcmp(keyword, "omegamax") == 0) StoreWithCheckDup_d(keyword, value, &StdI.OmegaMax);
    else if (strcmp(keyword, "omegamin") == 0) StoreWithCheckDup_d(keyword, value, &StdI.OmegaMin);
    else if (strcmp(keyword, "omegaim") == 0) StoreWithCheckDup_d(keyword, value, &StdI.OmegaIm);
    else if (strcmp(keyword, "restart") == 0) StoreWithCheckDup_sl(keyword, value, StdI.Restart);
    else if (strcmp(keyword, "spectrumqh") == 0) StoreWithCheckDup_d(keyword, value, &StdI.SpectrumQH);
    else if (strcmp(keyword, "spectrumql") == 0) StoreWithCheckDup_d(keyword, value, &StdI.SpectrumQL);
    else if (strcmp(keyword, "spectrumqw") == 0) StoreWithCheckDup_d(keyword, value, &StdI.SpectrumQW);
    else if (strcmp(keyword, "spectrumtype") == 0) StoreWithCheckDup_sl(keyword, value, StdI.SpectrumType);
    else if (strcmp(keyword, "2s") == 0) StoreWithCheckDup_i(keyword, value, &StdI.S2);
#elif defined(_mVMC)
    else if (strcmp(keyword, "a0hsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI.boxsub[0][2]);
    else if (strcmp(keyword, "a0lsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI.boxsub[0][1]);
    else if (strcmp(keyword, "a0wsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI.boxsub[0][0]);
    else if (strcmp(keyword, "a1hsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI.boxsub[1][2]);
    else if (strcmp(keyword, "a1lsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI.boxsub[1][1]);
    else if (strcmp(keyword, "a1wsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI.boxsub[1][0]);
    else if (strcmp(keyword, "a2hsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI.boxsub[2][2]);
    else if (strcmp(keyword, "a2lsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI.boxsub[2][1]);
    else if (strcmp(keyword, "a2wsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI.boxsub[2][0]);
    else if (strcmp(keyword, "complextype") == 0) StoreWithCheckDup_i(keyword, value, &StdI.ComplexType);
    else if (strcmp(keyword, "cparafilehead") == 0) StoreWithCheckDup_s(keyword, value, StdI.CParaFileHead);
    else if (strcmp(keyword, "dsroptredcut") == 0) StoreWithCheckDup_d(keyword, value, &StdI.DSROptRedCut);
    else if (strcmp(keyword, "dsroptstadel") == 0) StoreWithCheckDup_d(keyword, value, &StdI.DSROptStaDel);
    else if (strcmp(keyword, "dsroptstepdt") == 0) StoreWithCheckDup_d(keyword, value, &StdI.DSROptStepDt);
    else if (strcmp(keyword, "hsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI.Hsub);
    else if (strcmp(keyword, "lsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI.Lsub);
    else if (strcmp(keyword, "nvmccalmode") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NVMCCalMode);
    else if (strcmp(keyword, "ndataidxstart") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NDataIdxStart);
    else if (strcmp(keyword, "ndataqtysmp") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NDataQtySmp);
    else if (strcmp(keyword, "nlanczosmode") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NLanczosMode);
    else if (strcmp(keyword, "nmptrans") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NMPTrans);
    else if (strcmp(keyword, "nspgaussleg") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NSPGaussLeg);
    else if (strcmp(keyword, "nsplitsize") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NSplitSize);
    else if (strcmp(keyword, "nsroptitrsmp") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NSROptItrSmp);
    else if (strcmp(keyword, "nsroptitrstep") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NSROptItrStep);
    else if (strcmp(keyword, "nstore") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NStore);
    else if (strcmp(keyword, "nvmcinterval") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NVMCInterval);
    else if (strcmp(keyword, "nvmcsample") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NVMCSample);
    else if (strcmp(keyword, "nvmcwarmup") == 0) StoreWithCheckDup_i(keyword, value, &StdI.NVMCWarmUp);
    else if (strcmp(keyword, "rndseed") == 0) StoreWithCheckDup_i(keyword, value, &StdI.RndSeed);
    else if (strcmp(keyword, "wsub") == 0) StoreWithCheckDup_i(keyword, value, &StdI.Wsub);
#endif
    else {
      fprintf(stdout, "ERROR ! Unsupported Keyword !\n");
      StdFace_exit(-1);
    }
  }
  fflush(fp);
  fclose(fp);
  /*
  Check the model
  */
  StdI.lGC = 0;
  StdI.lBoost = 0;
  if (strcmp(StdI.model, "fermionhubbard") == 0
    || strcmp(StdI.model, "hubbard") == 0)
    strcpy(StdI.model, "hubbard\0");
  else if(strcmp(StdI.model, "fermionhubbardgc") == 0
    || strcmp(StdI.model, "hubbardgc") == 0) {
    strcpy(StdI.model, "hubbard\0");
    StdI.lGC = 1;
  }
  else if (strcmp(StdI.model, "spin") == 0)
    strcpy(StdI.model, "spin\0");
  else if (strcmp(StdI.model, "spingc") == 0) {
    strcpy(StdI.model, "spin\0");
    StdI.lGC = 1;
  }
#if defined(_HPhi)
  else if(strcmp(StdI.model, "spingcboost") == 0) {
    strcpy(StdI.model, "spin\0");
    StdI.lGC = 1;
    StdI.lBoost = 1;
  }
#endif
  else if (strcmp(StdI.model, "kondolattice") == 0
    || strcmp(StdI.model, "kondo") == 0) {
    strcpy(StdI.model, "kondo\0");
  }
  else if(strcmp(StdI.model, "kondolatticegc") == 0
    || strcmp(StdI.model, "kondogc") == 0) {
    strcpy(StdI.model, "kondo\0");
    StdI.lGC = 1;
  }
  else UnsupportedSystem(StdI.model, StdI.lattice);

#if defined(_mVMC)
  if (StdI.lGC == 1) {
    printf("\nSorry, Grandcanonical has not been supported !\n\n");
    StdFace_exit(-1);
  }
#endif
  /*
  Generate Hamiltonian definition files
  */
  if (strcmp(StdI.lattice, "chain") == 0
    || strcmp(StdI.lattice, "chainlattice") == 0) StdFace_Chain(&StdI, StdI.model);
  else if (strcmp(StdI.lattice, "facecenteredorthorhombic") == 0
    || strcmp(StdI.lattice, "fcorthorhombic") == 0
    || strcmp(StdI.lattice, "fco") == 0) StdFace_FCOrtho(&StdI, StdI.model);
  else if (strcmp(StdI.lattice, "honeycomb") == 0
    || strcmp(StdI.lattice, "honeycomblattice") == 0) StdFace_Honeycomb(&StdI, StdI.model);
  else if (strcmp(StdI.lattice, "kagome") == 0
    || strcmp(StdI.lattice, "kagomelattice") == 0) StdFace_Kagome(&StdI, StdI.model);
  else if (strcmp(StdI.lattice, "ladder") == 0
    || strcmp(StdI.lattice, "ladderlattice") == 0) StdFace_Ladder(&StdI, StdI.model);
  else if (strcmp(StdI.lattice, "orthorhombic") == 0
    || strcmp(StdI.lattice, "simpleorthorhombic") == 0) StdFace_Orthorhombic(&StdI, StdI.model);
  else if (strcmp(StdI.lattice, "pyrochlore") == 0) StdFace_Pyrochlore(&StdI, StdI.model);
  else if (strcmp(StdI.lattice, "tetragonal") == 0
    || strcmp(StdI.lattice, "tetragonallattice") == 0
    || strcmp(StdI.lattice, "square") == 0
    || strcmp(StdI.lattice, "squarelattice") == 0) StdFace_Tetragonal(&StdI, StdI.model);
  else if (strcmp(StdI.lattice, "triangular") == 0
    || strcmp(StdI.lattice, "triangularlattice") == 0) StdFace_Triangular(&StdI, StdI.model);
  else if (strcmp(StdI.lattice, "wannier90") == 0) StdFace_Wannier90(&StdI, StdI.model);
  else UnsupportedSystem(StdI.model, StdI.lattice);
  /**/
#if defined(_HPhi)
  StdFace_LargeValue(&StdI);
  /*
  Generate Hamiltonian for Boost
  */
  if (StdI.lBoost == 1) {
    if (strcmp(StdI.lattice, "chain") == 0
      || strcmp(StdI.lattice, "chainlattice") == 0) StdFace_Chain_Boost(&StdI);
    else if (strcmp(StdI.lattice, "honeycomb") == 0
      || strcmp(StdI.lattice, "honeycomblattice") == 0) StdFace_Honeycomb_Boost(&StdI);
    else if (strcmp(StdI.lattice, "kagome") == 0
      || strcmp(StdI.lattice, "kagomelattice") == 0) StdFace_Kagome_Boost(&StdI);
    else if (strcmp(StdI.lattice, "ladder") == 0
      || strcmp(StdI.lattice, "ladderlattice") == 0) StdFace_Ladder_Boost(&StdI);
    else UnsupportedSystem(StdI.model, StdI.lattice);
  }
#endif
  /**/
  fprintf(stdout, "\n");
  fprintf(stdout, "######  Print Expert input files  ######\n");
  fprintf(stdout, "\n");
  PrintLocSpin(&StdI);
  PrintTrans(&StdI);
  PrintInteractions(&StdI);
#if defined(_HPhi)
  PrintExcitation(&StdI);
  PrintCalcMod(&StdI);
#elif defined(_mVMC)
  StdFace_PrintVal_i("ComplexType", &StdI.ComplexType, 0);
  StdFace_generate_orb(&StdI);
  StdFace_Proj(&StdI);
  PrintJastrow(&StdI);
  PrintGutzwiller(&StdI);
  PrintOrb(&StdI);
#endif
  CheckModPara(&StdI);
  PrintModPara(&StdI);
  CheckOutputMode(&StdI);
  Print1Green(&StdI);
  Print2Green(&StdI);
  PrintNamelist(&StdI);
  /*
  Finalize All
  */
  free(StdI.locspinflag);
  for (ktrans = 0; ktrans < StdI.ntrans; ktrans++) {
    free(StdI.transindx[ktrans]);
  }
  free(StdI.transindx);
  free(StdI.trans);
  for (kintr = 0; kintr < StdI.nintr; kintr++) {
    free(StdI.intrindx[kintr]);
  }
  free(StdI.intrindx);
  free(StdI.intr);

  fprintf(stdout, "\n######  Input files are generated.  ######\n\n");

}/*void StdFace_main*/
