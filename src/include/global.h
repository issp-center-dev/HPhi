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

#ifndef HPHI_GLOBAL_H
#define HPHI_GLOBAL_H

#include <complex.h>
#include <math.h>
#define D_FileNameMax 256
#define MPIFALSE -1
#define FALSE 0
#define TRUE 1

/**
 * Kind of electron. For ver.0.1, ITINERANT=1. From ver.0.2, ITINERANT=0.
 **/
#define ITINERANT 0
#define LOCSPIN 1

double complex *v0; 
double complex *v1;
double complex *v2;
double complex *v1buf;

//For calcSpectrum
//[TODO]: This array can be omitted
double complex *v1Org;

//test
double complex *vg;

double *alpha,*beta;
double complex **vec;
double *list_Diagonal; 
double *list_D_num;
long unsigned int *list_1;
long unsigned int *list_1buf;
long unsigned int *list_2_1;
long unsigned int *list_2_2;
long unsigned int *HilbertNumToSz;

/*[s] For Spectrum */
long unsigned int *list_1_org;
long unsigned int *list_1buf_org;
long unsigned int *list_2_1_org;
long unsigned int *list_2_2_org;
long unsigned int *list_jb_org;
long unsigned int *list_2_1_Sz_org;
long unsigned int *list_2_2_Sz_org;
/*[e] For Spectrum */

/*[s] For Lanczos */
//double *eigen_vec;
int     initial_mode;
/*[e] For Lanczos */

/*[s] For TPQ*/
double LargeValue;
int    NumAve;
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
const char* cFileNameOutputEigen;
const char* cFileNameInputEigen;
const char* cFileNameCalcDynamicalGreen;
const char* cFileNameTridiagonalMatrixComponents;
const char* cFileNameLanczosOutputVector;
const char* cFileNameLanczosInputVector;

//For TPQ
const char* cFileNameSSRand;
const char* cFileNameTPQStep;
const char* cFileNameNormRand;
const char* cFileNameFlctRand;
const char* cFileName1BGreen_TPQ;
const char* cFileName2BGreen_TPQ;
const char* cFileNameOutputVector;
const char* cFileNameInputVector;

//For FullDiag
const char* cFileNamePhys_FullDiag;
const char* cFileNamePhys_FullDiag_GC;
const char* cFileName1BGreen_FullDiag;
const char* cFileName2BGreen_FullDiag;
const char* cFileNamePhys_FullDiag_Ham;

//For Spectrum
const char* cFileNameOutputRestartVec;

//For Time Evolution
const char* cFileNameInputVec;


//For Error
const char* cFileNameErrorSz;

//For Timer
//int NTimer=100;
double *Timer, *TimerStart;


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
