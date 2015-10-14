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
#include "../include/wrapperMPI.h"

static void TrimSpaceQuote(char *value);
static void StoreWithCheckDup_s(char *keyword, char *valuestring, char *value);
static void StoreWithCheckDup_i(char *keyword, char *valuestring, int *value);
static void StoreWithCheckDup_d(char *keyword, char *valuestring, double *value);
static void StdFace_ResetVals();

static void PrintLocSpin();
static void PrintTrans();
static void PrintInter();
static void PrintNamelist();
static void PrintCalcMod(char *method, int FlgTemp, char *model,int ioutputmode);
static void PrintModPara(int nup, int ndown, int Lanczos_max, int initial_iv, int nvec, 
  int exct,  int LanczosEps, int LanczosTarget, int WRITE, int READ,
  int NumAve, int ExpecInterval, char* filehead);
static void Print1Green(int ioutputmode);
static void Print2Green(int ioutputmode);

static int CheckOutputMode(char* outputmode);

static void UnsupportedSystem(char *model, char *lattice);

static void CheckModPara(int *nup, int *ndown, char* model,
  int *nelec, int *Sz2, int *Lanczos_max, int *initial_iv, int *nvec,
  int *exct, int *LanczosEps, int *LanczosTarget, int *WRITE, int *READ,
  int *NumAve, int *ExpecInterval, char* filehead);

/**
 *
 * Main routine for the standard mode
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_main(char *fname  /**< [in] Input file name for the standard mode */){
  FILE *fp;
  int ktrans, kintr;
  int FlgTemp, Lanczos_max, initial_iv, nvec, exct, 
    LanczosEps, LanczosTarget, WRITE, READ, nup, ndown, 
    NumAve, ExpecInterval, Sz2, nelec, ioutputmode;
  char ctmpline[256];
  char *keyword, *value;
  char model[256], lattice[256], method[256], outputmode[256], filehead[256];

  fprintf(stdoutMPI, "\n######  Standard Intarface Mode STARTS  ######\n");
  if ((fp = fopenMPI(fname, "r")) == NULL){
    fprintf(stdoutMPI, "\n  ERROR !  Cannot open input file %s !\n\n", fname);
    exitMPI(-1);
  }
  else{
    fprintf(stdoutMPI, "\n  Open Standard-Mode Inputfile %s \n\n", fname);
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
  //WRITE         = 9999;
  //READ          = 9999;
  NumAve        = 9999;
  ExpecInterval = 9999;
  StdFace_ResetVals();

  while (fgetsMPI(ctmpline, 256, fp) != NULL){

    TrimSpaceQuote(ctmpline);
    if (strncmp(ctmpline, "//", 2) == 0){
      fprintf(stdoutMPI, "  Skipping a line.\n");
      continue;
    }
    else if (ctmpline[0] == '\0'){
      fprintf(stdoutMPI, "  Skipping a line.\n");
      continue;
    }
    keyword = strtok(ctmpline, "=");
    value = strtok(NULL, "=");
    if (value == NULL){
      fprintf(stdoutMPI, "\n  ERROR !  \"=\" is NOT found !\n\n");
      exitMPI(-1);
    }
    TrimSpaceQuote(keyword);
    TrimSpaceQuote(value);
    fprintf(stdoutMPI, "  KEYWORD : %-20s | VALUE : %s \n", keyword, value);

    if (strcmp(keyword, "h") == 0) StoreWithCheckDup_d(keyword, value, &h);
    else if (strcmp(keyword, "exct") == 0) StoreWithCheckDup_i(keyword, value, &exct);
    else if (strcmp(keyword, "expecinterval") == 0) StoreWithCheckDup_i(keyword, value, &ExpecInterval);
    else if (strcmp(keyword, "filehead") == 0) StoreWithCheckDup_s(keyword, value, filehead);
    else if (strcmp(keyword, "flgtemp") == 0) StoreWithCheckDup_i(keyword, value, &FlgTemp);
    else if (strcmp(keyword, "g") == 0) StoreWithCheckDup_d(keyword, value, &D);
    else if (strcmp(keyword, "gamma") == 0) StoreWithCheckDup_d(keyword, value, &Gamma);
    else if (strcmp(keyword, "initial_iv") == 0) StoreWithCheckDup_i(keyword, value, &initial_iv);
    else if (strcmp(keyword, "j") == 0) StoreWithCheckDup_d(keyword, value, &J);
    else if (strcmp(keyword, "j0") == 0) StoreWithCheckDup_d(keyword, value, &J0);
    else if (strcmp(keyword, "j1") == 0) StoreWithCheckDup_d(keyword, value, &J1);
    else if (strcmp(keyword, "j2") == 0) StoreWithCheckDup_d(keyword, value, &J2);
    else if (strcmp(keyword, "jx") == 0) StoreWithCheckDup_d(keyword, value, &Jx);
    else if (strcmp(keyword, "jx0") == 0) StoreWithCheckDup_d(keyword, value, &Jx0);
    else if (strcmp(keyword, "jx1") == 0) StoreWithCheckDup_d(keyword, value, &Jx1);
    else if (strcmp(keyword, "jx2") == 0) StoreWithCheckDup_d(keyword, value, &Jx2);
    else if (strcmp(keyword, "jxy") == 0) StoreWithCheckDup_d(keyword, value, &Jxy);
    else if (strcmp(keyword, "jxy0") == 0) StoreWithCheckDup_d(keyword, value, &Jxy0);
    else if (strcmp(keyword, "jxy1") == 0) StoreWithCheckDup_d(keyword, value, &Jxy1);
    else if (strcmp(keyword, "jxy2") == 0) StoreWithCheckDup_d(keyword, value, &Jxy2);
    else if (strcmp(keyword, "jxy'") == 0) StoreWithCheckDup_d(keyword, value, &Jxyp);
    else if (strcmp(keyword, "jx'") == 0) StoreWithCheckDup_d(keyword, value, &Jxp);
    else if (strcmp(keyword, "jy") == 0) StoreWithCheckDup_d(keyword, value, &Jy);
    else if (strcmp(keyword, "jy0") == 0) StoreWithCheckDup_d(keyword, value, &Jy0);
    else if (strcmp(keyword, "jy1") == 0) StoreWithCheckDup_d(keyword, value, &Jy1);
    else if (strcmp(keyword, "jy2") == 0) StoreWithCheckDup_d(keyword, value, &Jy2);
    else if (strcmp(keyword, "jy'") == 0) StoreWithCheckDup_d(keyword, value, &Jyp);
    else if (strcmp(keyword, "jz") == 0) StoreWithCheckDup_d(keyword, value, &Jz);
    else if (strcmp(keyword, "jz0") == 0) StoreWithCheckDup_d(keyword, value, &Jz0);
    else if (strcmp(keyword, "jz1") == 0) StoreWithCheckDup_d(keyword, value, &Jz1);
    else if (strcmp(keyword, "jz2") == 0) StoreWithCheckDup_d(keyword, value, &Jz2);
    else if (strcmp(keyword, "jz'") == 0) StoreWithCheckDup_d(keyword, value, &Jzp);
    else if (strcmp(keyword, "j'") == 0) StoreWithCheckDup_d(keyword, value, &Jp);
    else if (strcmp(keyword, "j''") == 0) StoreWithCheckDup_d(keyword, value, &Jp);
    else if (strcmp(keyword, "k") == 0) StoreWithCheckDup_d(keyword, value, &K);
    else if (strcmp(keyword, "l") == 0) StoreWithCheckDup_i(keyword, value, &L);
    else if (strcmp(keyword, "lanczoseps") == 0) StoreWithCheckDup_i(keyword, value, &LanczosEps);
    else if (strcmp(keyword, "lanczostarget") == 0) StoreWithCheckDup_i(keyword, value, &LanczosTarget);
    else if (strcmp(keyword, "lanczos_max") == 0) StoreWithCheckDup_i(keyword, value, &Lanczos_max);
    else if (strcmp(keyword, "largevalue") == 0) StoreWithCheckDup_i(keyword, value, &LargeValue);
    else if (strcmp(keyword, "lattice") == 0) StoreWithCheckDup_s(keyword, value, lattice);
    else if (strcmp(keyword, "method") == 0) StoreWithCheckDup_s(keyword, value, method);
    else if (strcmp(keyword, "model") == 0) StoreWithCheckDup_s(keyword, value, model);
    else if (strcmp(keyword, "outputmode") == 0) StoreWithCheckDup_s(keyword, value, outputmode);
    else if (strcmp(keyword, "mu") == 0) StoreWithCheckDup_d(keyword, value, &mu);
    else if (strcmp(keyword, "nelec") == 0) StoreWithCheckDup_i(keyword, value, &nelec);
    else if (strcmp(keyword, "numave") == 0) StoreWithCheckDup_i(keyword, value, &NumAve);
    else if (strcmp(keyword, "nvec") == 0) StoreWithCheckDup_i(keyword, value, &nvec);
    else if (strcmp(keyword, "2sz") == 0) StoreWithCheckDup_i(keyword, value, &Sz2);
    else if (strcmp(keyword, "t") == 0) StoreWithCheckDup_d(keyword, value, &t);
    else if (strcmp(keyword, "t0") == 0) StoreWithCheckDup_d(keyword, value, &t0);
    else if (strcmp(keyword, "t1") == 0) StoreWithCheckDup_d(keyword, value, &t1);
    else if (strcmp(keyword, "t2") == 0) StoreWithCheckDup_d(keyword, value, &t2);
    else if (strcmp(keyword, "t'") == 0) StoreWithCheckDup_d(keyword, value, &tp);
    else if (strcmp(keyword, "t''") == 0) StoreWithCheckDup_d(keyword, value, &tp);
    else if (strcmp(keyword, "u") == 0) StoreWithCheckDup_d(keyword, value, &U);
    else if (strcmp(keyword, "v") == 0) StoreWithCheckDup_d(keyword, value, &V);
    else if (strcmp(keyword, "v0") == 0) StoreWithCheckDup_d(keyword, value, &V0);
    else if (strcmp(keyword, "v1") == 0) StoreWithCheckDup_d(keyword, value, &V1);
    else if (strcmp(keyword, "v2") == 0) StoreWithCheckDup_d(keyword, value, &V2);
    else if (strcmp(keyword, "v'") == 0) StoreWithCheckDup_d(keyword, value, &Vp);
    else if (strcmp(keyword, "v''") == 0) StoreWithCheckDup_d(keyword, value, &Vp);
    else if (strcmp(keyword, "w") == 0) StoreWithCheckDup_i(keyword, value, &W);
    else {
      fprintf(stdoutMPI, "ERROR ! Unsupported Keyword !\n");
      exitMPI(-1);
    }
  }
  fclose(fp);
  /*
   Generate Hamiltonian definition files
  */
  if (strcmp(model, "fermionhubbard") == 0){
    if (strcmp(lattice, "squarelattice") == 0) FermionHubbard_SquareLattice(nelec, 0);
    else if (strcmp(lattice, "chainlattice") == 0) FermionHubbard_ChainLattice(nelec, 0);
    else if (strcmp(lattice, "triangularlattice") == 0) FermionHubbard_TriangularLattice(nelec, 0);
    else if (strcmp(lattice, "honeycomblattice") == 0) FermionHubbard_HoneycombLattice(nelec, 0);
    else UnsupportedSystem(model, lattice);
  }
  else if (strcmp(model, "fermionhubbardgc") == 0
    || strcmp(model, "hubbardgc") == 0){
    if (strcmp(lattice, "squarelattice") == 0) FermionHubbard_SquareLattice(nelec, 1);
    else if (strcmp(lattice, "chainlattice") == 0) FermionHubbard_ChainLattice(nelec, 1);
    else if (strcmp(lattice, "triangularlattice") == 0) FermionHubbard_TriangularLattice(nelec, 1);
    else if (strcmp(lattice, "honeycomblattice") == 0) FermionHubbard_HoneycombLattice(nelec, 1);
    else UnsupportedSystem(model, lattice);
  }
  else if (strcmp(model, "spin") == 0
    || strcmp(model, "spingc") == 0){
    if (strcmp(lattice, "squarelattice") == 0) Spin_SquareLattice();
    else if (strcmp(lattice, "chainlattice") == 0) Spin_ChainLattice();
    else if (strcmp(lattice, "triangularlattice") == 0) Spin_TriangularLattice();
    else if (strcmp(lattice, "honeycomblattice") == 0) Spin_HoneycombLattice();
    else UnsupportedSystem(model, lattice);
  }
  else if (strcmp(model, "kondolattice") == 0){
    if (strcmp(lattice, "squarelattice") == 0) KondoLattice_SquareLattice(nelec, 0);
    else if (strcmp(lattice, "chainlattice") == 0) KondoLattice_ChainLattice(nelec, 0);
    else if (strcmp(lattice, "triangularlattice") == 0) KondoLattice_TriangularLattice(nelec, 0);
    else if (strcmp(lattice, "honeycomblattice") == 0) KondoLattice_HoneycombLattice(nelec, 0);
    else UnsupportedSystem(model, lattice);
  }
  else if (strcmp(model, "kondolatticegc") == 0
    || strcmp(model, "kondogc") == 0){
    if (strcmp(lattice, "squarelattice") == 0) KondoLattice_SquareLattice(nelec, 1);
    else if (strcmp(lattice, "chainlattice") == 0) KondoLattice_ChainLattice(nelec, 1);
    else if (strcmp(lattice, "triangularlattice") == 0) KondoLattice_TriangularLattice(nelec, 1);
    else if (strcmp(lattice, "honeycomblattice") == 0) KondoLattice_HoneycombLattice(nelec, 1);
    else UnsupportedSystem(model, lattice);
  }
  else UnsupportedSystem(model,lattice);
  /**/
  fprintf(stdoutMPI, "\n");
  /**/
  CheckModPara(&nup, &ndown, model,
    &nelec, &Sz2, &Lanczos_max, &initial_iv, &nvec,
    &exct, &LanczosEps, &LanczosTarget, &WRITE, &READ,
    &NumAve, &ExpecInterval,filehead);
  ioutputmode = CheckOutputMode(outputmode);
  /**/
  fprintf(stdoutMPI, "\n");
  fprintf(stdoutMPI, "######  Print Expart input files  ######\n");
  fprintf(stdoutMPI, "\n");
  PrintLocSpin();
  PrintTrans();
  PrintInter();
  PrintNamelist();
  PrintCalcMod(method, FlgTemp, model, ioutputmode);
  PrintModPara(nup, ndown,  Lanczos_max, initial_iv, nvec,  exct,
    LanczosEps, LanczosTarget, WRITE, READ,
    NumAve, ExpecInterval,filehead);  
  Print1Green(ioutputmode);
  Print2Green(ioutputmode);
  /*
   Finalize All
  */
  free(locspinflag);
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    free(transindx[ktrans]);
  }
  free(transindx);
  free(trans);
  for (kintr = 0; kintr < nintr; kintr++){
    free(intrindx[kintr]);
  }
  free(intrindx);
  free(intr);

  fprintf(stdoutMPI, "\n######  Input files are generated.  ######\n\n");

}

/**
 *
 * Clear grobal variables in the standard mode
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void StdFace_ResetVals(){
  /*
  Parameters for LATTICE
  */
  a = 9999.9;
  L = 9999;
  W = 9999;
  /*
  Parameters for MODEL
  */
  mu = 9999.9;
  t = 9999.9;
  tp = 9999.9;
  tpp = 9999.9;
  t0 = 9999.9;
  t1 = 9999.9;
  t2 = 9999.9;
  U = 9999.9;
  V = 9999.9;
  Vp = 9999.9;
  Vpp = 9999.9;
  V0 = 9999.9;
  V1 = 9999.9;
  V2 = 9999.9;
  J = 9999.9;
  Jp = 9999.9;
  Jpp = 9999.9;
  J0 = 9999.9;
  J1 = 9999.9;
  J2 = 9999.9;
  /**/
  Jx = 9999.9;
  Jy = 9999.9;
  Jz = 9999.9;
  Jxy = 9999.9;
  Jx0 = 9999.9;
  Jy0 = 9999.9;
  Jz0 = 9999.9;
  Jxy0 = 9999.9;
  Jx1 = 9999.9;
  Jy1 = 9999.9;
  Jz1 = 9999.9;
  Jxy1 = 9999.9;
  Jx2 = 9999.9;
  Jy2 = 9999.9;
  Jz2 = 9999.9;
  Jxy2 = 9999.9;
  Jxp = 9999.9;
  Jyp = 9999.9;
  Jzp = 9999.9;
  Jxyp = 9999.9;
  h = 9999.9;
  Gamma = 9999.9;
  D = 9999.9;
  K = 9999.9;
  /**/
  LargeValue = 9999;
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
    fprintf(stdoutMPI, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    exitMPI(-1);
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
    fprintf(stdoutMPI, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    exitMPI(-1);
  }
  else{
    sscanf(valuestring, "%d", value);
  }
}

/**
 *
 * Store an input value into the valiable (real)
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
    fprintf(stdoutMPI, "ERROR !  Keyword %s is duplicated ! \n", keyword);
    exitMPI(-1);
  }
  else{
    sscanf(valuestring, "%lf", value);
  }

}

/**
*
* Print the locspin file
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void PrintLocSpin(){
  FILE *fp;
  int isite, nlocspin;

  nlocspin = 0;
  for (isite = 0; isite < nsite; isite++)
    if (locspinflag[isite] == 0) nlocspin = nlocspin + 1;

  fp = fopenMPI("zlocspn.def", "w");
  fprintf(fp, "================================ \n");
  fprintf(fp, "NlocalSpin %5d  \n", nlocspin);
  fprintf(fp, "================================ \n");
  fprintf(fp, "========i_0LocSpn_1IteElc ====== \n");
  fprintf(fp, "================================ \n");

  for (isite = 0; isite < nsite; isite++)
    fprintf(fp, "%5d  %5d\n", isite, locspinflag[isite]);

  fclose(fp);
  fprintf(stdoutMPI, "    zlocspin.def is written.\n");
}

/**
*
* Print the transfer file
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintTrans(){
  FILE *fp;
  int jtrans, ktrans, ntrans0;

  for (jtrans = 0; jtrans < ntrans; jtrans++){
    for (ktrans = jtrans + 1; ktrans < ntrans; ktrans++){
      if ( transindx[jtrans][0] == transindx[ktrans][0]
        && transindx[jtrans][1] == transindx[ktrans][1]
        && transindx[jtrans][2] == transindx[ktrans][2]
        && transindx[jtrans][3] == transindx[ktrans][3]){
        trans[jtrans] = trans[jtrans] + trans[ktrans];
        trans[ktrans] = 0.0;
      }
    }
  }

  ntrans0 = 0;
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    if (fabs(trans[ktrans]) > 0.000001) ntrans0 = ntrans0 + 1;
  }

  fp = fopenMPI("zTrans.def", "w");
  fprintf(fp, "======================== \n");
  fprintf(fp, "NTransfer %7d  \n", ntrans0);
  fprintf(fp, "======================== \n");
  fprintf(fp, "========i_j_s_tijs====== \n");
  fprintf(fp, "======================== \n");

  ntrans0 = 0;
  for (ktrans = 0; ktrans < ntrans; ktrans++){
    if (fabs(trans[ktrans]) > 0.000001)
      fprintf(fp, "%5d %5d %5d %5d %10.6f  0.000000\n",
      transindx[ktrans][0], transindx[ktrans][1], transindx[ktrans][2], transindx[ktrans][3],
      trans[ktrans]);
  }

  fclose(fp);
  fprintf(stdoutMPI, "      zTrans.def is written.\n");
}

/**
*
* Print zInterAll.def
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
static void PrintInter(){
  FILE *fp;
  int jintr, kintr, nintr0;

  for (jintr = 0; jintr < nintr; jintr++){
    for (kintr = jintr + 1; kintr < nintr; kintr++){
      if ( intrindx[jintr][0] == intrindx[kintr][0]
        && intrindx[jintr][1] == intrindx[kintr][1]
        && intrindx[jintr][2] == intrindx[kintr][2]
        && intrindx[jintr][3] == intrindx[kintr][3]
        && intrindx[jintr][4] == intrindx[kintr][4]
        && intrindx[jintr][5] == intrindx[kintr][5]
        && intrindx[jintr][6] == intrindx[kintr][6]
        && intrindx[jintr][7] == intrindx[kintr][7]){
        intr[jintr] = intr[jintr] + intr[kintr];
        intr[kintr] = 0.0;
      }
    }
  }

  nintr0 = 0;
  for (kintr = 0; kintr < nintr; kintr++){
    if (fabs(intr[kintr]) > 0.000001) nintr0 = nintr0 + 1;
  }

  fp = fopenMPI("zInterAll.def", "w");
  fprintf(fp, "====================== \n");
  fprintf(fp, "NInterAll %7d  \n", nintr0);
  fprintf(fp, "====================== \n");
  fprintf(fp, "========zInterAll===== \n");
  fprintf(fp, "====================== \n");

  nintr0 = 0;
  for (kintr = 0; kintr < nintr; kintr++){
    if (fabs(intr[kintr]) > 0.000001)
      fprintf(fp, "%5d %5d %5d %5d %5d %5d %5d %5d %10.6f  0.000000\n",
      intrindx[kintr][0], intrindx[kintr][1], intrindx[kintr][2], intrindx[kintr][3],
      intrindx[kintr][4], intrindx[kintr][5], intrindx[kintr][6], intrindx[kintr][7],
      intr[kintr]);
  }

  fclose(fp);
  fprintf(stdoutMPI, "   zInterAll.def is written.\n");
}

/**
 *
 * Print namelist.def  
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void PrintNamelist(){
  FILE *fp;

  fp = fopenMPI("namelist.def", "w");
  fprintf(fp, "CalcMod calcmod.def\n");
  fprintf(fp, "ModPara modpara.def\n");
  fprintf(fp, "LocSpin zlocspn.def\n");
  fprintf(fp, "Trans zTrans.def\n");
  fprintf(fp, "InterAll zInterAll.def\n");
  fprintf(fp, "OneBodyG greenone.def\n");
  fprintf(fp, "TwoBodyG greentwo.def\n");

  fclose(fp);
  fprintf(stdoutMPI, "    namelist.def is written.\n");
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
    fprintf(stdoutMPI, "ERROR ! Method is NOT specified !\n");
    exitMPI(-1);
  }
  else if (strcmp(method, "lanczos") == 0) iCalcType = 0;
  else if (strcmp(method, "tpq") == 0) iCalcType = 1;
  else if (strcmp(method, "fulldiag") == 0 || 
    strcmp(method, "alldiag") == 0 ||
    strcmp(method, "direct") == 0 ) iCalcType = 2;
  else{
    fprintf(stdoutMPI, "\n ERROR ! Unsupported Solver : %s\n", method);
    exitMPI(-1);
  }

  if (strcmp(model, "fermionhubbard") == 0) iCalcModel = 0;
  else if (strcmp(model, "spin") == 0) iCalcModel = 1;
  else if (strcmp(model, "kondolattice") == 0) iCalcModel = 2;
  else if (strcmp(model, "fermionhubbardgc") == 0
    || strcmp(model, "hubbardgc") == 0) iCalcModel = 3;
  else if (strcmp(model, "spingc") == 0) iCalcModel = 4;
  else if (strcmp(model, "kondolatticegc") == 0
    || strcmp(model, "kondogc") == 0) iCalcModel = 5;
  else{
    fprintf(stdoutMPI, "\n ERROR ! Unsupported Model : %s\n", model);
    exitMPI(-1);
  }

  if (ioutputmode == 2) ioutputmode2 = 0;
  else ioutputmode2 = ioutputmode;

  fp = fopenMPI("calcmod.def", "w");
  fprintf(fp, "#CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag\n");
  fprintf(fp, "#FlgFiniteTemperature= 0:Zero temperature, 1:Finite temperature. This parameter is active only for CalcType=2.\n");
  fprintf(fp, "#CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC \n");
  fprintf(fp, "CalcType %3d\n", iCalcType);
  fprintf(fp, "FlgFiniteTemperature %3d\n", FlgTemp);
  fprintf(fp, "CalcModel %3d\n", iCalcModel);
  fprintf(fp, "OutputMode %3d\n", ioutputmode2);
  fclose(fp);
  fprintf(stdoutMPI, "     calcmod.def is written.\n");
}

/**
 *
 * Print modpara.def
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void PrintModPara(
  int nup /**< [in]*/,
  int ndown /**< [in]*/,
  int Lanczos_max /**< [in]*/,
  int initial_iv /**< [in]*/,
  int nvec /**< [in]*/,
  int exct /**< [in]*/,
  int LanczosEps /**< [in]*/,
  int LanczosTarget /**< [in]*/,
  int WRITE /**< [in]*/, 
  int READ /**< [in]*/,
  int NumAve /**< [in]*/,
  int ExpecInterval /**< [in]*/,
  char* filehead /**< [in]*/)
{
  FILE *fp;

  fp = fopenMPI("modpara.def", "w");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "Model_Parameters   0\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "HPhi_Cal_Parameters\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "CDataFileHead  %s\n", filehead);
  fprintf(fp, "CParaFileHead  zqp\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "Nsite          %-5d\n", nsite);
  fprintf(fp, "Nup            %-5d\n", nup);
  fprintf(fp, "Ndown          %-5d\n", ndown);
  fprintf(fp, "Lanczos_max    %-5d\n", Lanczos_max);
  fprintf(fp, "initial_iv     %-5d\n", initial_iv);
  fprintf(fp, "nvec           %-5d\n", nvec);
  fprintf(fp, "exct           %-5d\n", exct);
  fprintf(fp, "LanczosEps     %-5d\n", LanczosEps);
  fprintf(fp, "LanczosTarget  %-5d\n", LanczosTarget);
  //fprintf(fp, "WRITE          %-5d\n", WRITE);
  //fprintf(fp, "READ           %-5d\n", READ);
  fprintf(fp, "LargeValue     %-5d\n", LargeValue);
  fprintf(fp, "NumAve         %-5d\n", NumAve);
  fprintf(fp, "ExpecInterval  %-5d\n", ExpecInterval);

  fclose(fp);
  fprintf(stdoutMPI, "     modpara.def is written.\n");
}

/**
 *
 * Print greenone.def
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void Print1Green(int ioutputmode /**< [in]*/){
  FILE *fp;
  int ngreen, igreen, isite, jsite, ispin,jspin;
  int **greenindx;

  if (ioutputmode == 0){
    ngreen = 0;
  }
  else if(ioutputmode == 1){
    ngreen = nsite * 2;
  }
  else{
    ngreen = nsite * 2 * nsite * 2;
  }
  greenindx = (int **)malloc(sizeof(int*) * (ngreen + 1));
  for (igreen = 0; igreen < ngreen; igreen++){
    greenindx[igreen] = (int *)malloc(sizeof(int) * 4);
  }
  if (ioutputmode == 1){
    igreen = 0;
    for (isite = 0; isite < nsite; isite++){
      for (ispin = 0; ispin < 2; ispin++){
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
    for (isite = 0; isite < nsite; isite++){
      for (ispin = 0; ispin < 2; ispin++){
        for (jsite = 0; jsite < nsite; jsite++){
          for (jspin = 0; jspin < 2; jspin++){
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

  fp = fopenMPI("greenone.def", "w");
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

  fprintf(stdoutMPI, "    greenone.def is written.\n");
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
static void Print2Green(int ioutputmode /**< [in]*/){
  FILE *fp;
  int ngreen, igreen;
  int site1, site2, site3, site4;
  int spin1, spin2, spin3, spin4;
  int **greenindx;

  if (ioutputmode == 0){
    ngreen = 0;
  }
  else if (ioutputmode == 1){
    ngreen = nsite * 2 * nsite * 2;
  }
  else{
    ngreen = nsite * 2 * nsite * 2 * nsite * 2 * nsite * 2;
  }
  greenindx = (int **)malloc(sizeof(int*) * (ngreen + 1));
  for (igreen = 0; igreen < ngreen; igreen++){
    greenindx[igreen] = (int *)malloc(sizeof(int) * 8);
  }
  if (ioutputmode == 1){
    igreen = 0;
    for (site1 = 0; site1 < nsite; site1++){
      for (spin1 = 0; spin1 < 2; spin1++){
        for (site2 = 0; site2 < nsite; site2++){
          for (spin2 = 0; spin2 < 2; spin2++){
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
    for (site1 = 0; site1 < nsite; site1++){
      for (spin1 = 0; spin1 < 2; spin1++){
        for (site2 = 0; site2 < nsite; site2++){
          for (spin2 = 0; spin2 < 2; spin2++){
            for (site3 = 0; site3 < nsite; site3++){
              for (spin3 = 0; spin3 < 2; spin3++){
                for (site4 = 0; site4 < nsite; site4++){
                  for (spin4 = 0; spin4 < 2; spin4++){
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

  fp = fopenMPI("greentwo.def", "w");
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

  fprintf(stdoutMPI, "    greentwo.def is written.\n");
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
  fprintf(stdoutMPI, "\nSorry, specified combination, \n");
  fprintf(stdoutMPI, "    MODEL : %s  \n", model);
  fprintf(stdoutMPI, "  LATTICE : %s, \n", lattice);
  fprintf(stdoutMPI, "is unsupported in the STANDARD MODE...\n");
  fprintf(stdoutMPI, "Please use the EXPART MODE, or write a NEW FUNCTION and post us.\n");
  exitMPI(-1);
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
    fprintf(stdoutMPI, "      ioutputmode = %-10d\n", ioutputmode);
  }
  else if (strcmp(outputmode, "cor") == 0
    || strcmp(outputmode, "corr") == 0
    || strcmp(outputmode, "correlation") == 0) {
    ioutputmode = 1;
    fprintf(stdoutMPI, "      ioutputmode = %-10d\n", ioutputmode);
  }
  else if (strcmp(outputmode, "****") == 0) {
    ioutputmode = 1;
    fprintf(stdoutMPI, "      ioutputmode = %-10d  ######  DEFAULT VALUE IS USED  ######\n", ioutputmode);
  }
  else if (strcmp(outputmode, "raw") == 0
    || strcmp(outputmode, "all") == 0
    || strcmp(outputmode, "full") == 0) {
    ioutputmode = 2;
    fprintf(stdoutMPI, "      ioutputmode = %-10d\n", ioutputmode);
  }
  else{
    fprintf(stdoutMPI, "\n ERROR ! Unsupported OutPutMode : %s\n", outputmode);
    exitMPI(-1);
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
  int *nup /**< [out] The number of electrons of down spin*/, 
  int *ndown /**< [out] The number of electrons of up spin*/,
  char* model /**< [inout]*/,
  int *nelec /**< [inout]*/,
  int *Sz2 /**< [inout]*/,
  int *Lanczos_max /**< [inout]*/,
  int *initial_iv /**< [inout]*/,
  int *nvec /**< [inout]*/,
  int *exct /**< [inout]*/,
  int *LanczosEps /**< [inout]*/,
  int *LanczosTarget /**< [inout]*/,
  int *WRITE /**< [inout]*/,
  int *READ /**< [inout]*/,
  int *NumAve /**< [inout]*/,
  int *ExpecInterval /**< [inout]*/,
  char* filehead /**< [inout]*/)
{
  if (strcmp(filehead, "****") == 0) {
    strcpy(filehead, "zvo");
    fprintf(stdoutMPI, "         filehead = %-12s######  DEFAULT VALUE IS USED  ######\n", filehead);
  }
  else fprintf(stdoutMPI, "         filehead = %-s\n", filehead);
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

    StdFace_RequiredVal_i("2Sz", *Sz2);
    StdFace_RequiredVal_i("nelec", *nelec);

    if (abs(*Sz2) > nsite){
      fprintf(stdoutMPI, "\n ERROR ! abs(2 * Sz) > nsite in Hubbard model ! \n");
      exitMPI(-1);
    }
    else if (*nelec > 2 * nsite){
      fprintf(stdoutMPI, "\n ERROR ! Nelec > 2 * nsite in Hubbard model ! \n");
      exitMPI(-1);
    }
    else if ((*nelec + *Sz2) % 2 != 0){
      fprintf(stdoutMPI, "\n ERROR ! (nelec + 2 * Sz) %% 2 != 0 in Hubbard model ! \n");
      exitMPI(-1);
    }
    else if (*nelec <= nsite && abs(*Sz2) > *nelec){
      fprintf(stdoutMPI, "\n ERROR ! nelec <= nsite && 2 * |Sz| > nelec in Hubbard model ! \n");
      exitMPI(-1);
    }
    else if (*nelec > nsite && abs(*Sz2) > 2 * nsite - *nelec){
      fprintf(stdoutMPI, "\n ERROR ! nelec > nsite && 2 * |Sz| > 2 * nsite - nelec in Hubbard model ! \n");
      exitMPI(-1);
    }
    else {
      *nup = (*nelec + *Sz2) / 2;
      *ndown = (*nelec - *Sz2) / 2;
    }
  }
  else if (strcmp(model, "spin") == 0) {

    StdFace_RequiredVal_i("2Sz", *Sz2);
    StdFace_NotUsed_i("nelec", *nelec);

    if (abs(*Sz2) > nsite){
      fprintf(stdoutMPI, "\n ERROR ! abs(2 * Sz) > nsite in Spin model ! \n");
      exitMPI(-1);
    }
    else if ((nsite + *Sz2) % 2 != 0){
      fprintf(stdoutMPI, "\n ERROR ! (nsite + 2 * Sz) %% 2 != 0 in Spin model ! \n");
      exitMPI(-1);
    }
    else{
      *nup = (nsite + *Sz2) / 2;
      *ndown = (nsite - *Sz2) / 2;
    }
  }
  else if (strcmp(model, "kondolattice") == 0) {

    StdFace_RequiredVal_i("2Sz", *Sz2);
    StdFace_RequiredVal_i("nelec", *nelec);

    if (abs(*Sz2) > nsite){
      fprintf(stdoutMPI, "\n ERROR ! abs(2 * Sz) > nsite in Hubbard model ! \n");
      exitMPI(-1);
    }
    else if (*nelec > nsite){
      fprintf(stdoutMPI, "\n ERROR ! Nelec_cond / 2 + Nelec_loc > nsite in Kondo model ! \n");
      exitMPI(-1);
    }
    else if ((*nelec + nsite / 2 + *Sz2) % 2 != 0){
      fprintf(stdoutMPI, "\n ERROR ! (nelec_cond + nelec_loc + 2 * Sz) %% 2 != 0 in Kondo model ! \n");
      exitMPI(-1);
    }
    else if (*nelec <= nsite / 2 && abs(*Sz2) > *nelec + nsite / 2){
      fprintf(stdoutMPI, "\n ERROR ! nelec_cond <= nsite / 2 && 2 * |Sz| > nelec_cond + nelec_loc in Kondo model ! \n");
      exitMPI(-1);
    }
    else if (*nelec > nsite / 2 && abs(*Sz2) > nsite / 2 * 3 - *nelec){
      fprintf(stdoutMPI, "\n ERROR ! nelec_cond > nsite / 2 && abs(Sz2) > nsite / 2 * 3 - nelec in Kondo model ! \n");
      exitMPI(-1);
    }
    else {
      *nup = (*nelec + nsite / 2 + *Sz2) / 2;
      *ndown = (*nelec + nsite / 2 - *Sz2) / 2;
    }

  }
  else if (strcmp(model, "fermionhubbardgc") == 0
    || strcmp(model, "hubbardgc") == 0
    || strcmp(model, "spingc") == 0
    || strcmp(model, "kondolatticegc") == 0
    || strcmp(model, "kondogc") == 0){
 
    StdFace_NotUsed_i("nelec", *nelec);
    StdFace_NotUsed_i("2Sz", *Sz2);

    *nup = nsite / 2;
    *ndown = nsite - *nup;
  }
  else{
    fprintf(stdoutMPI, "\n ERROR ! Unsupported Model !\n");
    exitMPI(-1);
  }


}
