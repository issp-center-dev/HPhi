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
#pragma once
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
const char* cParentOutputFolder = "./output/";

//For TimeKeep
const char* cFileNameTimeKeep="%s_TimeKeeper.dat";
const char* cFileNameSzTimeKeep="%s_sz_TimeKeeper.dat";

//For Check
const char* cFileNameCheckCoulombIntra="CHECK_CoulombIntra.dat";
const char* cFileNameCheckChemi="CHECK_Chemi.dat";
const char* cFileNameCheckInterU="CHECK_INTER_U.dat";
const char* cFileNameCheckHund="CHECK_Hund.dat";
const char* cFileNameCheckInterAll="CHECK_InterAll.dat";
const char* cFileNameCheckMemory="CHECK_Memory.dat";
const char* cFileNameCheckSdim="CHECK_Sdim.dat";

//For EDTrans
const char* cFileNameWarningOnTransfer="WarningOnTransfer.dat";

//For Lanczos
const char* cFileNameLanczosStep="%s_Lanczos_Step.dat";
const char* cFileNameEnergy_Lanczos= "%s_energy.dat";
const char* cFileNameEigenvalue_Lanczos= "Eigenvalue.dat";
const char* cFileNameEnergy_CG="%s_energy.dat";
const char* cFileName1BGreen_Lanczos="%s_cisajs.dat";
const char* cFileName1BGreen_CG="%s_cisajs.dat";
const char* cFileName2BGreen_Lanczos="%s_cisajscktalt.dat";
const char* cFileName2BGreen_CG="%s_cisajscktalt.dat";
const char* cFileNameTimeEV_CG="Time_EigenVector.dat";
const char* cFileNameListModel="ListForModel_Ns%d_Nup%dNdown%d.dat";
const char* cFileNameListKondo="ListForKondo_Ns%d_Ncond%d.dat";

//For TPQ
const char* cFileNameSSRand="SS_rand%d.dat";
const char* cFileNameTPQStep="%s_TPQ_Step.dat";
const char* cFileNameNormRand="Norm_rand%d.dat";
const char* cFileName1BGreen_TPQ="%s_cisajs_set%dstep%d.dat";
const char* cFileName2BGreen_TPQ="%s_cisajscktalt_set%dstep%d.dat";

//For FullDiag
const char* cFileNamePhys_FullDiag="%s_phys_Nup%d_Ndown%d.dat";
const char* cFileNamePhys_FullDiag_GC="%s_phys.dat";
const char* cFileName1BGreen_FullDiag="%s_cisajs_eigen%d.dat";
const char* cFileName2BGreen_FullDiag="%s_cisajscktalt_eigen%d.dat";

//For Error
const char* cFileNameErrorSz="Err_sz.dat";


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
