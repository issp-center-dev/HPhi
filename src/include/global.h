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

#ifndef HPHI_GLOBAL_H
#define HPHI_GLOBAL_H

#include <complex.h>
#include <math.h>
#define D_FileNameMax 256

double complex *v0; 
double complex *v1;
//test
double complex *vg;

double *alpha,*beta;
double complex **vec;
double *list_Diagonal; 
double *list_D_num;
long unsigned int *list_1;
long unsigned  int *list_2_1;
long unsigned  int *list_2_2;
long unsigned  int *list_jb;
int *list_3;

/*[s] For Lanczos */
//double *eigen_vec;
int     initial_mode;
/*[e] For Lanczos */

/*[s] For TPQ*/
double LargeValue;
int    NumAve, ExpecInterval;
long int global_iv;
int step_i, step_spin;
double **All_S, **All_C, **All_Sz;
double global_spn2;
double global_norm, global_1st_norm;
/*[e] For TPQ*/

/*[s] For All Diagonalization*/
double *list_num_up,*list_num_down;
double complex**Ham;
double complex **L_vec; 
//double *eigen_vec_2;
/*[e] For All Diagonalization*/

/*OutputPath Name*/
const char* cParentOutputFolder;

//For TimeKeep
const char* cFileNameTimeKeep;
const char* cFileNameSzTimeKeep;

//For Check
const char* cFileNameCheckCoulombIntra;
const char* cFileNameCheckChemi;
const char* cFileNameCheckInterU;
const char* cFileNameCheckHund;
const char* cFileNameCheckInterAll;
const char* cFileNameCheckMemory;
const char* cFileNameCheckSdim;

//For EDTrans
const char* cFileNameWarningOnTransfer;

//For Lanczos
const char* cFileNameLanczosStep;
const char* cFileNameEnergy_Lanczos;
const char* cFileNameEigenvalue_Lanczos;
const char* cFileNameEnergy_CG;
const char* cFileName1BGreen_Lanczos;
const char* cFileName1BGreen_CG;
const char* cFileName2BGreen_Lanczos;
const char* cFileName2BGreen_CG;
const char* cFileNameTimeEV_CG;
const char* cFileNameListModel;
const char* cFileNameListKondo;

//For TPQ
const char* cFileNameSSRand;
const char* cFileNameTPQStep;
const char* cFileNameNormRand;
const char* cFileName1BGreen_TPQ;
const char* cFileName2BGreen_TPQ;

//For FullDiag
const char* cFileNamePhys_FullDiag;
const char* cFileNamePhys_FullDiag_GC;
const char* cFileName1BGreen_FullDiag;
const char* cFileName2BGreen_FullDiag;

//For Error
const char* cFileNameErrorSz;

/********************************************************************/
/********************************************************************/
double eps;
double eps_CG;
double eps_Lanczos; 
double eps_Bisec;
double eps_Energy;
double eps_CheckImag0;
double dShiftBeta;
double eps_vec12;

#endif /* HPHI_GLOBAL_H */
