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

#include "global.h"

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
const char* cFileName3BGreen_Lanczos="%s_ThreeBody.dat";
const char* cFileName4BGreen_Lanczos="%s_FourBody.dat";
const char* cFileName6BGreen_Lanczos="%s_SixBody.dat";
const char* cFileNameTimeEV_CG="Time_EigenVector.dat";
const char* cFileNameListModel="ListForModel_Ns%d_Nup%dNdown%d.dat";
const char* cFileNameOutputEigen="%s_eigenvec_%d_rank_%d.dat";
const char* cFileNameInputEigen="%s_eigenvec_%d_rank_%d.dat";
const char* cFileNameCalcDynamicalGreen="%s_DynamicalGreen.dat";
const char* cFileNameTridiagonalMatrixComponents="%s_TMComponents.dat";


//For TPQ
const char* cFileNameSSRand="SS_rand%d.dat";
const char* cFileNameTPQStep="%s_Time_TPQ_Step.dat";
const char* cFileNameNormRand="Norm_rand%d.dat";
const char* cFileNameFlctRand="Flct_rand%d.dat";
const char* cFileName1BGreen_TPQ="%s_cisajs_set%dstep%d.dat";
const char* cFileName2BGreen_TPQ="%s_cisajscktalt_set%dstep%d.dat";
const char* cFileName3BGreen_TPQ="%s_ThreeBody_set%dstep%d.dat";
const char* cFileName4BGreen_TPQ="%s_FourBody_set%dstep%d.dat";
const char* cFileName6BGreen_TPQ="%s_SixBody_set%dstep%d.dat";
const char* cFileName1BGreen_TE="%s_cisajs_step%d.dat";
const char* cFileName2BGreen_TE="%s_cisajscktalt_step%d.dat";
const char* cFileName3BGreen_TE="%s_ThreeBody_step%d.dat";
const char* cFileName4BGreen_TE="%s_FourBody_step%d.dat";
const char* cFileName6BGreen_TE="%s_SixBody_step%d.dat";
const char* cFileNameOutputVector="tmpvec_set%d_rank_%d.dat";
const char* cFileNameInputVector="tmpvec_set%d_rank_%d.dat";

//Fot Time evolution
const char* cFileNameTEStep="%s_Time_TE_Step.dat";
const char* cFileNameSS="SS.dat";
const char* cFileNameNorm="Norm.dat";
const char* cFileNameFlct="Flct.dat";

//For FullDiag
const char* cFileNamePhys_FullDiag="%s_phys_Nup%d_Ndown%d.dat";
const char* cFileNamePhys_FullDiag_GC="%s_phys.dat";
const char* cFileName1BGreen_FullDiag="%s_cisajs_eigen%d.dat";
const char* cFileName2BGreen_FullDiag="%s_cisajscktalt_eigen%d.dat";
const char* cFileName3BGreen_FullDiag="%s_ThreeBody_eigen%d.dat";
const char* cFileName4BGreen_FullDiag="%s_FourBody_eigen%d.dat";
const char* cFileName6BGreen_FullDiag="%s_SixBody_eigen%d.dat";
const char* cFileNamePhys_FullDiag_Ham="%s_Ham.dat";

//For Spectrum
const char* cFileNameOutputRestartVec="%s_recalcvec_rank_%d.dat";
const char* cFileNameOutputExcitedVec="%s_excitedvec_rank_%d.dat";
//For Error
const char* cFileNameErrorSz="Err_sz.dat";

double complex *v0 = 0;
double complex *v1 = 0;
double complex *v2 = 0;
double complex *v1buf = 0;
double complex *v1Org = 0;
double complex *vg=0;
double *alpha = 0;
double *beta = 0;
double complex **vec = 0;
double *list_Diagonal = 0;
long unsigned int *list_1 = 0;
long unsigned int *list_1buf = 0;
long unsigned int *list_2 = 0;
long unsigned int *list_2_1 = 0;
long unsigned int *list_2_2 = 0;

/* Spectrum */
long unsigned int *list_1_org = 0;
long unsigned int *list_1buf_org = 0;
long unsigned int *list_2_1_org = 0;
long unsigned int *list_2_2_org = 0;

/* Lanczos */
int initial_mode = 0;

/* TPQ */
double LargeValue = 1.0e30;
int NumAve = 1;
int step_i = 1;
double global_norm = 0.0;
double global_1st_norm = 0.0;
int step_spin = 1;

/* All Diagonalization */
double complex **Ham = 0;
double complex **L_vec = 0;
#ifdef _SCALAPACK
double complex *Z_vec=0;
int descZ_vec[9] = {0};
#endif

/* Timer */
double *Timer = 0;
double *TimerStart = 0;

/* epsilons */
double eps = 1e-10;
double eps_CG = 1e-10;
double eps_Lanczos = 1e-10;
double eps_Energy = 1e-10;
double eps_CheckImag0 = 1e-10;

/* MPI */
int nproc = 1;
int myrank = 0;
int nthreads = 1;
FILE *stdoutMPI = 0;
