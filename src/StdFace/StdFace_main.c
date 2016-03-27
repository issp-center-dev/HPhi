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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "StdFace_vals.h"
#include "StdFace_ModelUtil.h"
#include <complex.h>

static void TrimSpaceQuote(char *value);
static void StoreWithCheckDup_s(char *keyword, char *valuestring, char *value);
static void StoreWithCheckDup_i(char *keyword, char *valuestring, int *value);
static void StoreWithCheckDup_d(char *keyword, char *valuestring, double *value);
static void StdFace_ResetVals(struct StdIntList *StdI);

static void PrintLocSpin(struct StdIntList *StdI);
static void PrintTrans(struct StdIntList *StdI);
static void PrintInter(struct StdIntList *StdI);
static void PrintNamelist(char *model);
static void PrintCalcMod(char *method, int FlgTemp, char *model,int ioutputmode);
static void PrintModPara(struct StdIntList *StdI, int Sz2, int nelec, 
  int Lanczos_max, int initial_iv, int nvec,
  int exct,  int LanczosEps, int LanczosTarget,
  int NumAve, int ExpecInterval, char* filehead);
static void Print1Green(struct StdIntList *StdI, int ioutputmode);
static void Print2Green(struct StdIntList *StdI,int ioutputmode);

static int CheckOutputMode(char* outputmode);

static void UnsupportedSystem(char *model, char *lattice);

static void CheckModPara(char* model,
  int *nelec, int *Sz2, int *Lanczos_max, int *initial_iv, int *nvec,
  int *exct, int *LanczosEps, int *LanczosTarget,
  int *NumAve, int *ExpecInterval, char* filehead);

/**
 *
 * Main routine for the standard mode
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_main(char *fname  /**< [in] Input file name for the standard mode */){

  struct StdIntList StdI;

  FILE *fp;
  int ktrans, kintr;
  int FlgTemp, Lanczos_max, initial_iv, nvec, exct, 
    LanczosEps, LanczosTarget, 
    NumAve, ExpecInterval, Sz2, nelec, ioutputmode;
  char ctmpline[256];
  char *keyword, *value;
  char model[256], lattice[256], method[256], outputmode[256], filehead[256];

  fprintf(stdout, "\n######  Standard Intarface Mode STARTS  ######\n");
  if ((fp = fopen(fname, "r")) == NULL){
    fprintf(stderr, "\n  ERROR !  Cannot open input file %s !\n\n", fname);
    exit(-1);
  }
  else{
    fprintf(stdout, "\n  Open Standard-Mode Inputfile %s \n\n", fname);
  }

  strcpy(model, "****");
  strcpy(lattice, "****");
  strcpy(method, "****");
  strcpy(outputmode, "****");
  strcpy(filehead, "****");
  FlgTemp = 1;
  nelec       = 9999;
  Sz2         = 9999;
  Lanczos_max = 9999;
  initial_iv    = 9999;
  nvec          = 9999;
  exct          = 9999;
  LanczosEps    = 9999;
  LanczosTarget = 9999;
  NumAve        = 9999;
  ExpecInterval = 9999;
  StdFace_ResetVals(&StdI);

  while (fgets(ctmpline, 256, fp) != NULL){

    TrimSpaceQuote(ctmpline);
    if (strncmp(ctmpline, "//", 2) == 0){
      fprintf(stdout, "  Skipping a line.\n");
      continue;
    }
    else if (ctmpline[0] == '\0'){
      fprintf(stdout, "  Skipping a line.\n");
      continue;
    }
    keyword = strtok(ctmpline, "=");
    value = strtok(NULL, "=");
    if (value == NULL){
      fprintf(stderr, "\n  ERROR !  \"=\" is NOT found !\n\n");
      exit(-1);
    }
    TrimSpaceQuote(keyword);
    TrimSpaceQuote(value);
    fprintf(stdout, "  KEYWORD : %-20s | VALUE : %s \n", keyword, value);

    if (strcmp(keyword, "h") == 0) StoreWithCheckDup_d(keyword, value, &StdI.h);
    else if (strcmp(keyword, "d") == 0) StoreWithCheckDup_d(keyword, value, &StdI.D);
    else if (strcmp(keyword, "exct") == 0) StoreWithCheckDup_i(keyword, value, &exct);
    else if (strcmp(keyword, "expecinterval") == 0) StoreWithCheckDup_i(keyword, value, &ExpecInterval);
    else if (strcmp(keyword, "filehead") == 0) StoreWithCheckDup_s(keyword, value, filehead);
    else if (strcmp(keyword, "flgtemp") == 0) StoreWithCheckDup_i(keyword, value, &FlgTemp);
    else if (strcmp(keyword, "gamma") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Gamma);
    else if (strcmp(keyword, "initial_iv") == 0) StoreWithCheckDup_i(keyword, value, &initial_iv);
    else if (strcmp(keyword, "j") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J);
    else if (strcmp(keyword, "j0") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J0);
    else if (strcmp(keyword, "j1") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1);
    else if (strcmp(keyword, "j1'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J1p);
    else if (strcmp(keyword, "j2") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2);
    else if (strcmp(keyword, "j2'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.J2p);
    else if (strcmp(keyword, "jx") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jx);
    else if (strcmp(keyword, "jx0") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jx0);
    else if (strcmp(keyword, "jx1") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jx1);
    else if (strcmp(keyword, "jx2") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jx2);
    else if (strcmp(keyword, "jxy") == 0) StoreWithCheckDup_c(keyword, value, &StdI.J[0][1]);
    else if (strcmp(keyword, "jxy0") == 0) StoreWithCheckDup_c(keyword, value, &StdI.J[0][1]0);
    else if (strcmp(keyword, "jxy1") == 0) StoreWithCheckDup_c(keyword, value, &StdI.J[0][1]1);
    else if (strcmp(keyword, "jxy2") == 0) StoreWithCheckDup_c(keyword, value, &StdI.J[0][1]2);
    else if (strcmp(keyword, "jxy'") == 0) StoreWithCheckDup_c(keyword, value, &StdI.J[0][1]p);
    else if (strcmp(keyword, "jx'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jxp);
    else if (strcmp(keyword, "jy") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jy);
    else if (strcmp(keyword, "jy0") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jy0);
    else if (strcmp(keyword, "jy1") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jy1);
    else if (strcmp(keyword, "jy2") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jy2);
    else if (strcmp(keyword, "jy'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jyp);
    else if (strcmp(keyword, "jz") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jz);
    else if (strcmp(keyword, "jz0") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jz0);
    else if (strcmp(keyword, "jz1") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jz1);
    else if (strcmp(keyword, "jz2") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jz2);
    else if (strcmp(keyword, "jz'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jzp);
    else if (strcmp(keyword, "j'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jp);
    else if (strcmp(keyword, "j''") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Jp);
    else if (strcmp(keyword, "k") == 0) StoreWithCheckDup_d(keyword, value, &StdI.K);
    else if (strcmp(keyword, "l") == 0) StoreWithCheckDup_i(keyword, value, &StdI.L);
    else if (strcmp(keyword, "lanczoseps") == 0) StoreWithCheckDup_i(keyword, value, &LanczosEps);
    else if (strcmp(keyword, "lanczostarget") == 0) StoreWithCheckDup_i(keyword, value, &LanczosTarget);
    else if (strcmp(keyword, "lanczos_max") == 0) StoreWithCheckDup_i(keyword, value, &Lanczos_max);
    else if (strcmp(keyword, "largevalue") == 0) StoreWithCheckDup_d(keyword, value, &StdI.LargeValue);
    else if (strcmp(keyword, "lattice") == 0) StoreWithCheckDup_s(keyword, value, lattice);
    else if (strcmp(keyword, "method") == 0) StoreWithCheckDup_s(keyword, value, method);
    else if (strcmp(keyword, "model") == 0) StoreWithCheckDup_s(keyword, value, model);
    else if (strcmp(keyword, "outputmode") == 0) StoreWithCheckDup_s(keyword, value, outputmode);
    else if (strcmp(keyword, "mu") == 0) StoreWithCheckDup_d(keyword, value, &StdI.mu);
    else if (strcmp(keyword, "nelec") == 0) StoreWithCheckDup_i(keyword, value, &nelec);
    else if (strcmp(keyword, "numave") == 0) StoreWithCheckDup_i(keyword, value, &NumAve);
    else if (strcmp(keyword, "nvec") == 0) StoreWithCheckDup_i(keyword, value, &nvec);
    else if (strcmp(keyword, "2sz") == 0) StoreWithCheckDup_i(keyword, value, &Sz2);
    else if (strcmp(keyword, "2s") == 0) StoreWithCheckDup_i(keyword, value, &StdI.S2);
    else if (strcmp(keyword, "t") == 0) StoreWithCheckDup_c(keyword, value, &StdI.t);
    else if (strcmp(keyword, "t0") == 0) StoreWithCheckDup_c(keyword, value, &StdI.t0);
    else if (strcmp(keyword, "t1") == 0) StoreWithCheckDup_c(keyword, value, &StdI.t1);
    else if (strcmp(keyword, "t1'") == 0) StoreWithCheckDup_c(keyword, value, &StdI.t1p);
    else if (strcmp(keyword, "t2") == 0) StoreWithCheckDup_c(keyword, value, &StdI.t2);
    else if (strcmp(keyword, "t'") == 0) StoreWithCheckDup_c(keyword, value, &StdI.tp);
    else if (strcmp(keyword, "t''") == 0) StoreWithCheckDup_c(keyword, value, &StdI.tp);
    else if (strcmp(keyword, "u") == 0) StoreWithCheckDup_d(keyword, value, &StdI.U);
    else if (strcmp(keyword, "v") == 0) StoreWithCheckDup_d(keyword, value, &StdI.V);
    else if (strcmp(keyword, "v0") == 0) StoreWithCheckDup_d(keyword, value, &StdI.V0);
    else if (strcmp(keyword, "v1") == 0) StoreWithCheckDup_d(keyword, value, &StdI.V1);
    else if (strcmp(keyword, "v1'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.V1p);
    else if (strcmp(keyword, "v2") == 0) StoreWithCheckDup_d(keyword, value, &StdI.V2);
    else if (strcmp(keyword, "v'") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Vp);
    else if (strcmp(keyword, "v''") == 0) StoreWithCheckDup_d(keyword, value, &StdI.Vp);
    else if (strcmp(keyword, "w") == 0) StoreWithCheckDup_i(keyword, value, &StdI.W);
    else {
      fprintf(stderr, "ERROR ! Unsupported Keyword !\n");
      exit(-1);
    }
  }
  fclose(fp);
  /*
   Generate Hamiltonian definition files
  */
  if (strcmp(model, "fermionhubbard") == 0){
    if (strcmp(lattice, "squarelattice") == 0) FermionHubbard_SquareLattice(&StdI, nelec, 0);
    else if (strcmp(lattice, "chainlattice") == 0) FermionHubbard_ChainLattice(&StdI, nelec, 0);
    else if (strcmp(lattice, "triangularlattice") == 0) FermionHubbard_TriangularLattice(&StdI, nelec, 0);
    else if (strcmp(lattice, "honeycomblattice") == 0) FermionHubbard_HoneycombLattice(&StdI, nelec, 0);
    else if (strcmp(lattice, "ladder") == 0) FermionHubbard_Ladder(&StdI, nelec, 0);
    else UnsupportedSystem(model, lattice);
  }
  else if (strcmp(model, "fermionhubbardgc") == 0
    || strcmp(model, "hubbardgc") == 0){
    if (strcmp(lattice, "squarelattice") == 0) FermionHubbard_SquareLattice(&StdI, nelec, 1);
    else if (strcmp(lattice, "chainlattice") == 0) FermionHubbard_ChainLattice(&StdI, nelec, 1);
    else if (strcmp(lattice, "triangularlattice") == 0) FermionHubbard_TriangularLattice(&StdI, nelec, 1);
    else if (strcmp(lattice, "honeycomblattice") == 0) FermionHubbard_HoneycombLattice(&StdI, nelec, 1);
    else if (strcmp(lattice, "ladder") == 0) FermionHubbard_Ladder(&StdI, nelec, 1);
    else UnsupportedSystem(model, lattice);
  }
  else if (strcmp(model, "spin") == 0){
    if (strcmp(lattice, "squarelattice") == 0) Spin_SquareLattice(&StdI, Sz2, 0);
    else if (strcmp(lattice, "chainlattice") == 0) Spin_ChainLattice(&StdI, Sz2, 0);
    else if (strcmp(lattice, "triangularlattice") == 0) Spin_TriangularLattice(&StdI, Sz2, 0);
    else if (strcmp(lattice, "honeycomblattice") == 0) Spin_HoneycombLattice(&StdI, Sz2, 0);
    else if (strcmp(lattice, "ladder") == 0) Spin_Ladder(&StdI, Sz2, 0);
    else UnsupportedSystem(model, lattice);
  }
  else if (strcmp(model, "spingc") == 0){
    if (strcmp(lattice, "squarelattice") == 0) Spin_SquareLattice(&StdI, Sz2, 1);
    else if (strcmp(lattice, "chainlattice") == 0) Spin_ChainLattice(&StdI, Sz2, 1);
    else if (strcmp(lattice, "triangularlattice") == 0) Spin_TriangularLattice(&StdI, Sz2, 1);
    else if (strcmp(lattice, "honeycomblattice") == 0) Spin_HoneycombLattice(&StdI, Sz2, 1);
    else if (strcmp(lattice, "ladder") == 0) Spin_Ladder(&StdI, Sz2, 1);
    else UnsupportedSystem(model, lattice);
  }
  else if (strcmp(model, "spingcboost") == 0) {
    if (strcmp(lattice, "squarelattice") == 0) UnsupportedSystem(model, lattice);
    else if (strcmp(lattice, "chainlattice") == 0) Spin_ChainLattice_Boost(&StdI, Sz2, 1);
    else if (strcmp(lattice, "triangularlattice") == 0) UnsupportedSystem(model, lattice);
    else if (strcmp(lattice, "honeycomblattice") == 0) Spin_HoneycombLattice_Boost(&StdI, Sz2, 1);
    else if (strcmp(lattice, "kagome") == 0) Spin_Kagome_Boost(&StdI, Sz2, 1);
    else if (strcmp(lattice, "ladder") == 0) Spin_Ladder_Boost(&StdI, Sz2, 1);
    else UnsupportedSystem(model, lattice);
  }
  else if (strcmp(model, "kondolattice") == 0){
    if (strcmp(lattice, "squarelattice") == 0) KondoLattice_SquareLattice(&StdI, nelec, 0);
    else if (strcmp(lattice, "chainlattice") == 0) KondoLattice_ChainLattice(&StdI, nelec, 0);
    else if (strcmp(lattice, "triangularlattice") == 0) KondoLattice_TriangularLattice(&StdI, nelec, 0);
    else if (strcmp(lattice, "honeycomblattice") == 0) KondoLattice_HoneycombLattice(&StdI, nelec, 0);
    else if (strcmp(lattice, "ladder") == 0) KondoLattice_Ladder(&StdI, nelec, 0);
    else UnsupportedSystem(model, lattice);
  }
  else if (strcmp(model, "kondolatticegc") == 0
    || strcmp(model, "kondogc") == 0){
    if (strcmp(lattice, "squarelattice") == 0) KondoLattice_SquareLattice(&StdI, nelec, 1);
    else if (strcmp(lattice, "chainlattice") == 0) KondoLattice_ChainLattice(&StdI, nelec, 1);
    else if (strcmp(lattice, "triangularlattice") == 0) KondoLattice_TriangularLattice(&StdI, nelec, 1);
    else if (strcmp(lattice, "honeycomblattice") == 0) KondoLattice_HoneycombLattice(&StdI, nelec, 1);
    else if (strcmp(lattice, "ladder") == 0) KondoLattice_Ladder(&StdI, nelec, 1);
    else UnsupportedSystem(model, lattice);
  }
  else UnsupportedSystem(model,lattice);
  /**/
  fprintf(stdout, "\n");
  /**/
  CheckModPara(model,
    &nelec, &Sz2, &Lanczos_max, &initial_iv, &nvec,
    &exct, &LanczosEps, &LanczosTarget,
    &NumAve, &ExpecInterval,filehead);
  ioutputmode = CheckOutputMode(outputmode);
  /**/
  fprintf(stdout, "\n");
  fprintf(stdout, "######  Print Expert input files  ######\n");
  fprintf(stdout, "\n");
  PrintLocSpin(&StdI);
  PrintTrans(&StdI);
  PrintInter(&StdI);
  PrintNamelist(model);
  PrintCalcMod(method, FlgTemp, model, ioutputmode);
  PrintModPara(&StdI,Sz2, nelec, Lanczos_max, initial_iv, nvec, exct,
    LanczosEps, LanczosTarget,NumAve, ExpecInterval,filehead);  
  Print1Green(&StdI, ioutputmode);
  Print2Green(&StdI, ioutputmode);
  /*
   Finalize All
  */
  free(StdI.locspinflag);
  for (ktrans = 0; ktrans < StdI.ntrans; ktrans++){
    free(StdI.transindx[ktrans]);
  }
  free(StdI.transindx);
  free(StdI.trans);
  for (kintr = 0; kintr < StdI.nintr; kintr++){
    free(StdI.intrindx[kintr]);
  }
  free(StdI.intrindx);
  free(StdI.intr);
 
  fprintf(stdout, "\n######  Input files are generated.  ######\n\n");

}

/**
 *
 * Clear grobal variables in the standard mode
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void StdFace_ResetVals(struct StdIntList *StdI){
  int i, j;
  /*
  Parameters for LATTICE
  */
  StdI->a = 9999.9;
  StdI->L = 9999;
  StdI->W = 9999;
  /*
  Parameters for MODEL
  */
  StdI->mu = 9999.9;
  StdI->t = 9999.9;
  StdI->tp = 9999.9;
  StdI->t0 = 9999.9;
  StdI->t1 = 9999.9;
  StdI->t1p = 9999.9;
  StdI->t2 = 9999.9;
  StdI->U = 9999.9;
  StdI->V = 9999.9;
  StdI->Vp = 9999.9;
  StdI->V0 = 9999.9;
  StdI->V1 = 9999.9;
  StdI->V1p = 9999.9;
  StdI->V2 = 9999.9;
  StdI->JAll = 9999.9;
  StdI->JpAll = 9999.9;
  StdI->J0All = 9999.9;
  StdI->J1All = 9999.9;
  StdI->J1pAll = 9999.9;
  StdI->J2All = 9999.9;
  StdI->J2pAll = 9999.9;
  /**/
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      StdI->J[i][j] = 9999.9;
      StdI->Jp[i][j] = 9999.9;
      StdI->J0[i][j] = 9999.9;
      StdI->J1[i][j] = 9999.9;
      StdI->J1p[i][j] = 9999.9;
      StdI->J2[i][j] = 9999.9;
      StdI->J2p[i][j] = 9999.9;
      StdI->D[i][j] = 0.0;
    }
  }
  StdI->h = 9999.9;
  StdI->Gamma = 9999.9;
  StdI->K = 9999.9;
  StdI->D[2][2] = 9999.9;
  /**/
  StdI->LargeValue = 9999.9;
  StdI->S2 = 9999;
}

/**
 *
 * Remove , : space etc. from keyword and value in an iput file
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
      value[ii] != ',' &&
      value[ii] != ':' &&
      value[ii] != ';' &&
      value[ii] != '\"' &&
      value[ii] != '\b' &&
      value[ii] != '\\' &&
      value[ii] != '\v' &&
      value[ii] != '\n' &&
      value[ii] != '\0'){
      value2[valuelen2] = tolower(value[ii]);
      valuelen2++;
    }
  }

  strncpy(value, value2, valuelen2);
  value[valuelen2] = '\0';

}

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
    fprintf(stderr, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    exit(-1);
  }
  else{
    strcpy(value, valuestring);
  }
}

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
  if (*value != 9999){
    fprintf(stderr, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    exit(-1);
  }
  else{
    sscanf(valuestring, "%d", value);
  }
}

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

  if (*value != 9999.9){
    fprintf(stderr, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    exit(-1);
  }
  else{
    sscanf(valuestring, "%lf", value);
  }

}

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
  _Dcomplex *value /**< [out] */)
{
  int num;
  double value_r, value_i;

  if (creal(*value) != 9999.9) {
    fprintf(stderr, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    exit(-1);
  }
  else {
    num = sscanf(valuestring, "%lf %lf", &value_r, &value_i);
    if (num == 1) {

    }
    else if (num == 2) {
      *value = _Cbuild(value_r, value_i);
    }
    else {
      fprintf(stderr, "\nERROR! in reading complex number.\n\n");
      exit(-1);
    }
  }

}

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

  fp = fopen("zlocspn.def", "w");
  fprintf(fp, "================================ \n");
  fprintf(fp, "NlocalSpin %5d  \n", nlocspin);
  fprintf(fp, "================================ \n");
  fprintf(fp, "========i_0LocSpn_1IteElc ====== \n");
  fprintf(fp, "================================ \n");

  for (isite = 0; isite < StdI->nsite; isite++)
    fprintf(fp, "%5d  %5d\n", isite, StdI->locspinflag[isite]);

  fclose(fp);
  fprintf(stdout, "    zlocspin.def is written.\n");
}

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
    }
  }

  ntrans0 = 0;
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    if (fabs(StdI->trans[ktrans]) > 0.000001) ntrans0 = ntrans0 + 1;
  }

  fp = fopen("zTrans.def", "w");
  fprintf(fp, "======================== \n");
  fprintf(fp, "NTransfer %7d  \n", ntrans0);
  fprintf(fp, "======================== \n");
  fprintf(fp, "========i_j_s_tijs====== \n");
  fprintf(fp, "======================== \n");

  ntrans0 = 0;
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    if (fabs(StdI->trans[ktrans]) > 0.000001)
      fprintf(fp, "%5d %5d %5d %5d %25.15f  0.000000\n",
        StdI->transindx[ktrans][0], StdI->transindx[ktrans][1], 
        StdI->transindx[ktrans][2], StdI->transindx[ktrans][3],
        StdI->trans[ktrans]);
  }

  fclose(fp);
  fprintf(stdout, "      zTrans.def is written.\n");
}

/**
*
* Print zInterAll.def
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintInter(struct StdIntList *StdI){
  FILE *fp;
  int jintr, kintr, nintr0;

  for (jintr = 0; jintr < StdI->nintr; jintr++){
    for (kintr = jintr + 1; kintr < StdI->nintr; kintr++){
      if ( 
        (StdI->intrindx[jintr][0] == StdI->intrindx[kintr][0]
        && StdI->intrindx[jintr][1] == StdI->intrindx[kintr][1]
        && StdI->intrindx[jintr][2] == StdI->intrindx[kintr][2]
        && StdI->intrindx[jintr][3] == StdI->intrindx[kintr][3]
        && StdI->intrindx[jintr][4] == StdI->intrindx[kintr][4]
        && StdI->intrindx[jintr][5] == StdI->intrindx[kintr][5]
        && StdI->intrindx[jintr][6] == StdI->intrindx[kintr][6]
        && StdI->intrindx[jintr][7] == StdI->intrindx[kintr][7] )
        ||
        (StdI->intrindx[jintr][0] == StdI->intrindx[kintr][4]
        && StdI->intrindx[jintr][1] == StdI->intrindx[kintr][5]
        && StdI->intrindx[jintr][2] == StdI->intrindx[kintr][6]
        && StdI->intrindx[jintr][3] == StdI->intrindx[kintr][7]
        && StdI->intrindx[jintr][4] == StdI->intrindx[kintr][0]
        && StdI->intrindx[jintr][5] == StdI->intrindx[kintr][1]
        && StdI->intrindx[jintr][6] == StdI->intrindx[kintr][2]
        && StdI->intrindx[jintr][7] == StdI->intrindx[kintr][3])
        ){
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
        ){
        StdI->intr[jintr] = StdI->intr[jintr] - StdI->intr[kintr];
        StdI->intr[kintr] = 0.0;
      }
    }
  }

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
    }
  }

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
  }
 
  nintr0 = 0;
  for (kintr = 0; kintr < StdI->nintr; kintr++){
    if (fabs(StdI->intr[kintr]) > 0.000001) nintr0 = nintr0 + 1;
  }

  fp = fopen("zInterAll.def", "w");
  fprintf(fp, "====================== \n");
  fprintf(fp, "NInterAll %7d  \n", nintr0);
  fprintf(fp, "====================== \n");
  fprintf(fp, "========zInterAll===== \n");
  fprintf(fp, "====================== \n");

  nintr0 = 0;
  for (kintr = 0; kintr < StdI->nintr; kintr++){
    if (fabs(StdI->intr[kintr]) > 0.000001)
      fprintf(fp, "%5d %5d %5d %5d %5d %5d %5d %5d %25.15f  0.000000\n",
        StdI->intrindx[kintr][0], StdI->intrindx[kintr][1], 
        StdI->intrindx[kintr][2], StdI->intrindx[kintr][3],
        StdI->intrindx[kintr][4], StdI->intrindx[kintr][5], 
        StdI->intrindx[kintr][6], StdI->intrindx[kintr][7],
        StdI->intr[kintr]);
  }

  fclose(fp);
  fprintf(stdout, "   zInterAll.def is written.\n");
}

/**
 *
 * Print namelist.def  
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void PrintNamelist(char *model){
  FILE *fp;

  fp = fopen("namelist.def", "w");
  fprintf(fp, "CalcMod calcmod.def\n");
  fprintf(fp, "ModPara modpara.def\n");
  fprintf(fp, "LocSpin zlocspn.def\n");
  fprintf(fp, "Trans zTrans.def\n");
  fprintf(fp, "InterAll zInterAll.def\n");
  fprintf(fp, "OneBodyG greenone.def\n");
  fprintf(fp, "TwoBodyG greentwo.def\n");

  if (strcmp(model, "spingcboost") == 0) 
    fprintf(fp, "Boost boost.def\n");

  fclose(fp);
  fprintf(stdout, "    namelist.def is written.\n");
}

/**
 *
 * Print calcmod.def
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void PrintCalcMod(
  char *method /**< [in]*/, 
  int FlgTemp /**< [in]*/,
  char *model /**< [in]*/,
  int ioutputmode /**< [in]*/)
{
  FILE *fp;
  int iCalcType, iCalcModel, ioutputmode2;

  if (strcmp(method, "****") == 0){
    fprintf(stderr, "ERROR ! Method is NOT specified !\n");
    exit(-1);
  }
  else if (strcmp(method, "lanczos") == 0) iCalcType = 0;
  else if (strcmp(method, "tpq") == 0) iCalcType = 1;
  else if (strcmp(method, "fulldiag") == 0 || 
    strcmp(method, "alldiag") == 0 ||
    strcmp(method, "direct") == 0 ) iCalcType = 2;
  else{
    fprintf(stderr, "\n ERROR ! Unsupported Solver : %s\n", method);
    exit(-1);
  }

  if (strcmp(model, "fermionhubbard") == 0) iCalcModel = 0;
  else if (strcmp(model, "spin") == 0) iCalcModel = 1;
  else if (strcmp(model, "kondolattice") == 0) iCalcModel = 2;
  else if (strcmp(model, "fermionhubbardgc") == 0
    || strcmp(model, "hubbardgc") == 0) iCalcModel = 3;
  else if (strcmp(model, "spingc") == 0
    || strcmp(model, "spingcboost") == 0) iCalcModel = 4;
  else if (strcmp(model, "kondolatticegc") == 0
    || strcmp(model, "kondogc") == 0) iCalcModel = 5;
  else{
    fprintf(stderr, "\n ERROR ! Unsupported Model : %s\n", model);
    exit(-1);
  }

  if (ioutputmode == 2) ioutputmode2 = 0;
  else ioutputmode2 = ioutputmode;

  fp = fopen("calcmod.def", "w");
  fprintf(fp, "#CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag\n");
  fprintf(fp, "#FlgFiniteTemperature= 0:Zero temperature, 1:Finite temperature. This parameter is active only for CalcType=2.\n");
  fprintf(fp, "#CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC \n");
  fprintf(fp, "CalcType %3d\n", iCalcType);
  fprintf(fp, "FlgFiniteTemperature %3d\n", FlgTemp);
  fprintf(fp, "CalcModel %3d\n", iCalcModel);
  fprintf(fp, "OutputMode %3d\n", ioutputmode2);
  fclose(fp);
  fprintf(stdout, "     calcmod.def is written.\n");
}

/**
 *
 * Print modpara.def
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void PrintModPara(
  struct StdIntList *StdI,
  int Sz2 /**< [in]*/,
  int nelec /**< [in]*/,
  int Lanczos_max /**< [in]*/,
  int initial_iv /**< [in]*/,
  int nvec /**< [in]*/,
  int exct /**< [in]*/,
  int LanczosEps /**< [in]*/,
  int LanczosTarget /**< [in]*/,
  int NumAve /**< [in]*/,
  int ExpecInterval /**< [in]*/,
  char* filehead /**< [in]*/)
{
  FILE *fp;

  fp = fopen("modpara.def", "w");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "Model_Parameters   0\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "HPhi_Cal_Parameters\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "CDataFileHead  %s\n", filehead);
  fprintf(fp, "CParaFileHead  zqp\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "Nsite          %-5d\n", StdI->nsite);
  if (Sz2 != 9999) fprintf(fp, "2Sz            %-5d\n", Sz2);
  if (nelec != 9999) fprintf(fp, "Ncond          %-5d\n", nelec);
  fprintf(fp, "Lanczos_max    %-5d\n", Lanczos_max);
  fprintf(fp, "initial_iv     %-5d\n", initial_iv);
  fprintf(fp, "nvec           %-5d\n", nvec);
  fprintf(fp, "exct           %-5d\n", exct);
  fprintf(fp, "LanczosEps     %-5d\n", LanczosEps);
  fprintf(fp, "LanczosTarget  %-5d\n", LanczosTarget);
  fprintf(fp, "LargeValue     %-25.15e\n", StdI->LargeValue);
  fprintf(fp, "NumAve         %-5d\n", NumAve);
  fprintf(fp, "ExpecInterval  %-5d\n", ExpecInterval);

  fclose(fp);
  fprintf(stdout, "     modpara.def is written.\n");
}

/**
 *
 * Print greenone.def
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void Print1Green(struct StdIntList *StdI,
  int ioutputmode /**< [in]*/){
  FILE *fp;
  int ngreen, igreen, isite, jsite, ispin,jspin, SiMax, SjMax;
  int **greenindx;

  if (ioutputmode == 0){
    ngreen = 0;
  }
  else if(ioutputmode == 1){
    ngreen = 0;
    for (isite = 0; isite < StdI->nsite; isite++){
      if (StdI->locspinflag[isite] == 0) ngreen = ngreen + 2;
      else ngreen = ngreen + (StdI->locspinflag[isite] + 1);
    }
  }
  else{
    ngreen = 0;
    for (isite = 0; isite < StdI->nsite; isite++){
      if (StdI->locspinflag[isite] == 0) ngreen = ngreen + 2;
      else ngreen = ngreen + (StdI->locspinflag[isite] + 1);
    }
    ngreen = ngreen * ngreen;
  }
  greenindx = (int **)malloc(sizeof(int*) * (ngreen + 1));
  for (igreen = 0; igreen < ngreen; igreen++){
    greenindx[igreen] = (int *)malloc(sizeof(int) * 4);
  }
  if (ioutputmode == 1){
    igreen = 0;
    for (isite = 0; isite < StdI->nsite; isite++){

      if (StdI->locspinflag[isite] == 0) SiMax = 1;
      else SiMax = StdI->locspinflag[isite];
      
      for (ispin = 0; ispin <= SiMax; ispin++){
        greenindx[igreen][0] = isite;
        greenindx[igreen][1] = ispin;
        greenindx[igreen][2] = isite;
        greenindx[igreen][3] = ispin;
        igreen++;
      }
    }
  }
  else if (ioutputmode == 2){
    igreen = 0;
    for (isite = 0; isite < StdI->nsite; isite++){

      if (StdI->locspinflag[isite] == 0) SiMax = 1;
      else SiMax = StdI->locspinflag[isite];

      for (ispin = 0; ispin <= SiMax; ispin++){
        for (jsite = 0; jsite < StdI->nsite; jsite++){

          if (StdI->locspinflag[jsite] == 0) SjMax = 1;
          else SjMax = StdI->locspinflag[jsite];

          for (jspin = 0; jspin <= SjMax; jspin++){
            greenindx[igreen][0] = isite;
            greenindx[igreen][1] = ispin;
            greenindx[igreen][2] = jsite;
            greenindx[igreen][3] = jspin;
            igreen++;
          }
        }
      }
    }
  }

  fp = fopen("greenone.def", "w");
  fprintf(fp, "===============================\n");
  fprintf(fp, "NCisAjs %10d\n", ngreen);
  fprintf(fp, "===============================\n");
  fprintf(fp, "======== Green functions ======\n");
  fprintf(fp, "===============================\n");
  for (igreen = 0; igreen < ngreen; igreen++){
    fprintf(fp,"%5d %5d %5d %5d\n",
      greenindx[igreen][0], greenindx[igreen][1], greenindx[igreen][2], greenindx[igreen][3]);
  }
  fclose(fp);

  fprintf(stdout, "    greenone.def is written.\n");
  //[s] free
  for (igreen = 0; igreen < ngreen; igreen++){
    free(greenindx[igreen]);
  }
  free(greenindx);
  //[e] free
}

/**
 *
 * Print greentwo.def
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void Print2Green(struct StdIntList *StdI,
  int ioutputmode /**< [in]*/){
  FILE *fp;
  int ngreen, igreen;
  int site1, site2, site3, site4;
  int spin1, spin2, spin3, spin4;
  int S1Max, S2Max, S3Max, S4Max;
  int **greenindx;

  if (ioutputmode == 0){
    ngreen = 0;
  }
  else if (ioutputmode == 1){
    ngreen = 0;
    for (site1 = 0; site1 < StdI->nsite; site1++){
      if (StdI->locspinflag[site1] == 0) ngreen = ngreen + 2;
      else ngreen = ngreen + (StdI->locspinflag[site1] + 1);
    }
    ngreen = ngreen * ngreen;
  }
  else{
    ngreen = 0;
    for (site1 = 0; site1 < StdI->nsite; site1++){
      if (StdI->locspinflag[site1] == 0) ngreen = ngreen + 2;
      else ngreen = ngreen + (StdI->locspinflag[site1] + 1);
    }
    ngreen = ngreen * ngreen * ngreen * ngreen;
  }
  greenindx = (int **)malloc(sizeof(int*) * (ngreen + 1));
  for (igreen = 0; igreen < ngreen; igreen++){
    greenindx[igreen] = (int *)malloc(sizeof(int) * 8);
  }
  if (ioutputmode == 1){
    igreen = 0;
    for (site1 = 0; site1 < StdI->nsite; site1++){

      if (StdI->locspinflag[site1] == 0) S1Max = 1;
      else S1Max = StdI->locspinflag[site1];

      for (spin1 = 0; spin1 <= S1Max; spin1++){
        for (site2 = 0; site2 < StdI->nsite; site2++){

          if (StdI->locspinflag[site2] == 0) S2Max = 1;
          else S2Max = StdI->locspinflag[site2];

          for (spin2 = 0; spin2 <= S2Max; spin2++){
            greenindx[igreen][0] = site1;
            greenindx[igreen][1] = spin1;
            greenindx[igreen][2] = site1;
            greenindx[igreen][3] = spin1;
            greenindx[igreen][4] = site2;
            greenindx[igreen][5] = spin2;
            greenindx[igreen][6] = site2;
            greenindx[igreen][7] = spin2;
           igreen++;
          }
        }
      }
    }
   }
  else if (ioutputmode == 2){
    igreen = 0;
    for (site1 = 0; site1 < StdI->nsite; site1++){

      if (StdI->locspinflag[site1] == 0) S1Max = 1;
      else S1Max = StdI->locspinflag[site1];

      for (spin1 = 0; spin1 <= S1Max; spin1++){
        for (site2 = 0; site2 < StdI->nsite; site2++){

          if (StdI->locspinflag[site1] != 0 && StdI->locspinflag[site2] != 0 
            && site1 != site2) continue;

          if (StdI->locspinflag[site2] == 0) S2Max = 1;
          else S2Max = StdI->locspinflag[site2];

          for (spin2 = 0; spin2 <= S2Max; spin2++){
            for (site3 = 0; site3 < StdI->nsite; site3++){

              if (StdI->locspinflag[site3] == 0) S3Max = 1;
              else S3Max = StdI->locspinflag[site3];

              for (spin3 = 0; spin3 <= S3Max; spin3++){
                for (site4 = 0; site4 < StdI->nsite; site4++){

                  if (StdI->locspinflag[site3] != 0 && StdI->locspinflag[site4] != 0
                    && site3 != site4) continue;

                  if (StdI->locspinflag[site4] == 0) S4Max = 1;
                  else S4Max = StdI->locspinflag[site4];

                  for (spin4 = 0; spin4 <= S4Max; spin4++){
                    greenindx[igreen][0] = site1;
                    greenindx[igreen][1] = spin1;
                    greenindx[igreen][2] = site2;
                    greenindx[igreen][3] = spin2;
                    greenindx[igreen][4] = site3;
                    greenindx[igreen][5] = spin3;
                    greenindx[igreen][6] = site4;
                    greenindx[igreen][7] = spin4;
                    igreen++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  ngreen = igreen;

  fp = fopen("greentwo.def", "w");
  fprintf(fp, "=============================================\n");
  fprintf(fp, "NCisAjsCktAltDC %10d\n", ngreen);
  fprintf(fp, "=============================================\n");
  fprintf(fp, "======== Green functions for Sq AND Nq ======\n");
  fprintf(fp, "=============================================\n");
  for (igreen = 0; igreen < ngreen; igreen++){
    fprintf(fp,"%5d %5d %5d %5d %5d %5d %5d %5d\n",
      greenindx[igreen][0], greenindx[igreen][1], greenindx[igreen][2], greenindx[igreen][3],
      greenindx[igreen][4], greenindx[igreen][5], greenindx[igreen][6], greenindx[igreen][7]);
  }
  fclose(fp);

  fprintf(stdout, "    greentwo.def is written.\n");
  //[s] free
  for (igreen = 0; igreen < ngreen; igreen++){
    free(greenindx[igreen]);
  }
  free(greenindx);
  //[e] free
}

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
  fprintf(stderr, "\nSorry, specified combination, \n");
  fprintf(stderr, "    MODEL : %s  \n", model);
  fprintf(stderr, "  LATTICE : %s, \n", lattice);
  fprintf(stderr, "is unsupported in the STANDARD MODE...\n");
  fprintf(stderr, "Please use the EXPART MODE, or write a NEW FUNCTION and post us.\n");
  exit(-1);
}

/**
 *
 * Verify outputmode
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static int CheckOutputMode(char* outputmode /**< [in]*/){
  int ioutputmode;
  /*
  Form for Correlation function
  */
  if (strcmp(outputmode, "non") == 0
    || strcmp(outputmode, "none") == 0
    || strcmp(outputmode, "off") == 0) {
    ioutputmode = 0;
    fprintf(stdout, "      ioutputmode = %-10d\n", ioutputmode);
  }
  else if (strcmp(outputmode, "cor") == 0
    || strcmp(outputmode, "corr") == 0
    || strcmp(outputmode, "correlation") == 0) {
    ioutputmode = 1;
    fprintf(stdout, "      ioutputmode = %-10d\n", ioutputmode);
  }
  else if (strcmp(outputmode, "****") == 0) {
    ioutputmode = 1;
    fprintf(stdout, "      ioutputmode = %-10d  ######  DEFAULT VALUE IS USED  ######\n", ioutputmode);
  }
  else if (strcmp(outputmode, "raw") == 0
    || strcmp(outputmode, "all") == 0
    || strcmp(outputmode, "full") == 0) {
    ioutputmode = 2;
    fprintf(stdout, "      ioutputmode = %-10d\n", ioutputmode);
  }
  else{
    fprintf(stderr, "\n ERROR ! Unsupported OutPutMode : %s\n", outputmode);
    exit(-1);
  }
  return(ioutputmode);
}

/**
 *
 * Summary numerical parameter check the combination of
 * the number of sites, total spin, the number of electrons
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void CheckModPara(
  char* model /**< [inout]*/,
  int *nelec /**< [inout]*/,
  int *Sz2 /**< [inout]*/,
  int *Lanczos_max /**< [inout]*/,
  int *initial_iv /**< [inout]*/,
  int *nvec /**< [inout]*/,
  int *exct /**< [inout]*/,
  int *LanczosEps /**< [inout]*/,
  int *LanczosTarget /**< [inout]*/,
  int *NumAve /**< [inout]*/,
  int *ExpecInterval /**< [inout]*/,
  char* filehead /**< [inout]*/)
{
  if (strcmp(filehead, "****") == 0) {
    strcpy(filehead, "zvo");
    fprintf(stdout, "         filehead = %-12s######  DEFAULT VALUE IS USED  ######\n", filehead);
  }
  else fprintf(stdout, "         filehead = %-s\n", filehead);
  /**/
  StdFace_PrintVal_i("Lanczos_max", Lanczos_max, 2000);
  StdFace_PrintVal_i("initial_iv", initial_iv, 1);
  StdFace_PrintVal_i("nvec", nvec, 1);
  StdFace_PrintVal_i("exct", exct, 1);
  StdFace_PrintVal_i("LanczosEps", LanczosEps, 14);
  StdFace_PrintVal_i("LanczosTarget", LanczosTarget, 2);
  StdFace_PrintVal_i("NumAve", NumAve, 5);
  StdFace_PrintVal_i("ExpecInterval", ExpecInterval, 20);
  /**/
  if (strcmp(model, "fermionhubbard") == 0) {

    StdFace_RequiredVal_i("nelec", *nelec);

  }
  else if (strcmp(model, "spin") == 0) {

    StdFace_RequiredVal_i("2Sz", *Sz2);
    StdFace_NotUsed_i("nelec", *nelec);
  }
  else if (strcmp(model, "kondolattice") == 0) {

    StdFace_RequiredVal_i("nelec", *nelec);

  }
  else if (strcmp(model, "fermionhubbardgc") == 0
    || strcmp(model, "hubbardgc") == 0
    || strcmp(model, "spingc") == 0
    || strcmp(model, "spingcboost") == 0
    || strcmp(model, "kondolatticegc") == 0
    || strcmp(model, "kondogc") == 0){
 
    StdFace_NotUsed_i("nelec", *nelec);
    StdFace_NotUsed_i("2Sz", *Sz2);
  }
  else{
    fprintf(stderr, "\n ERROR ! Unsupported Model !\n");
    exit(-1);
  }


}
