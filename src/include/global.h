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
#include <stdio.h>
#define D_FileNameMax 256
#define MPIFALSE -1
#define FALSE 0
#define TRUE 1

/**
 * Kind of electron. For ver.0.1, ITINERANT=1. From ver.0.2, ITINERANT=0.
 **/
#define ITINERANT 0
#define LOCSPIN 1

double complex *v0;  /**< A vector after multiplying Hamiltonian, @f$ v_0 = H v_1@f$.*/
double complex *v1;  /**< A vector before multiplying Hamiltonian, @f$ v_0 = H v_1@f$.*/
double complex *v2;  /**< A temporary vector for time evolution calculation, @f$ v2 = H*v1 = H^coef |psi(t)>@f$.*/
double complex *v1buf; /**< A temporary vector for MPI. */

//[s] For calcSpectrum
double complex *v1Org; /**< An input vector to calculate spectrum function.*/
double complex *vg; /**< A vector used in the CG mode.*/
//[e] For calcSpectrum

double *alpha,*beta; /**< Tridiagonal components used in Lanczos mode.*/
double complex **vec; /**< Eigen vectors.*/
double *list_Diagonal; /**< list for diagonal components.*/
long unsigned int *list_1; /**< list of getting real-space configuration for canonical state*/
long unsigned int *list_1buf;/**< list of getting real-space configuration for canonical state across processes*/
long unsigned int *list_2_1;/**< list to get index of list_1*/
long unsigned int *list_2_2;/**< list to get index of list_1*/

/*[s] For Spectrum */
long unsigned int *list_1_org; /**< list of getting real-space configuration for canonical state before excitation*/
long unsigned int *list_1buf_org;/**< list of getting real-space configuration for canonical state before excitation across processes*/
long unsigned int *list_2_1_org;/**< list to get index of list_1_org*/
long unsigned int *list_2_2_org;/**< list to get index of list_1_org*/
/*[e] For Spectrum */

/*[s] For Lanczos */
int     initial_mode;/**< mode to get initial state (0: use same random generator for MPI, 1: use each random generator for MPI)*/
/*[e] For Lanczos */

/*[s] For TPQ*/
double LargeValue;/**< constant value l for TPQ calculation.*/
int    NumAve;/**< Average number for TPQ calculation*/
int step_i;/**< step for TPQ calculation*/
double global_norm;/**< norm before normalization for TPQ calculation*/
double global_1st_norm;/**< 1-st norm for TPQ calculation*/
int step_spin;/**< output step for TE calculation.*/
/*[e] For TPQ*/

/*[s] For All Diagonalization*/
double complex**Ham; /**> Hamiltonian for full diagonalization. */
double complex **L_vec;/**> eigen vectors*/
/*[e] For All Diagonalization*/

const char* cParentOutputFolder; /**> Path to output results*/

//For TimeKeep
const char* cFileNameTimeKeep; /**> Name of the file to output calculation processing time*/
const char* cFileNameSzTimeKeep;/**> Name of the file to output sz calculation time*/

//For Check
const char* cFileNameCheckCoulombIntra;/**> Name of the file to check Coulomb Intra interactions.*/
const char* cFileNameCheckChemi;/**> Name of the file to check Chemical potentials.*/
const char* cFileNameCheckInterU;/**> Name of the file to check Inter U.*/
const char* cFileNameCheckHund;/**>Name of the file to check Hund interactions. */
const char* cFileNameCheckInterAll;/**> Name of the file to check InterAll interactions.*/
const char* cFileNameCheckMemory;/**> Name of the file to check memory.*/
const char* cFileNameCheckSdim;/**> Name of the file to check dimension.*/

//For EDTrans
const char* cFileNameWarningOnTransfer; /**> Name of the file to output warning of transfer integrals.*/

//For Lanczos
const char* cFileNameLanczosStep;/**> Name of the file to output Lanczos step.*/
const char* cFileNameEnergy_Lanczos;/**> Name of the file to output energies.*/
const char* cFileNameEigenvalue_Lanczos;/**> Name of the file to output eigen values.*/
const char* cFileNameEnergy_CG;/**> Name of the file to output energies obtained by CG method.*/
const char* cFileName1BGreen_Lanczos;/**> Name of the file to output One-Body Green's functions obtained by Lanczos method.*/
const char* cFileName1BGreen_CG;/**> Name of the file to output One-Body Green's functions obtained by CG method.*/
const char* cFileName2BGreen_Lanczos;/**> Name of the file to output Two-Body Green's functions obtained by Lanczos method.*/
const char* cFileName2BGreen_CG;/**> Name of the file to output Two-Body Green's functions obtained by CG method.*/
const char* cFileNameTimeEV_CG;/**> Name of the file to output time for getting eigen vector by CG method.*/
const char* cFileNameListModel;/**> Name of the file to output list.*/
const char* cFileNameOutputEigen;/**> Name of the file to output eigen vector.*/
const char* cFileNameInputEigen;/**> Name of the file to input eigen vector.*/
const char* cFileNameCalcDynamicalGreen;/**> Name of the file to output dynamical Green's function.*/
const char* cFileNameTridiagonalMatrixComponents;/**> Name of the file to output tridiagonal matrix components.*/

//For TPQ
const char* cFileNameSSRand;/**> Name of the SS_rand file.*/
const char* cFileNameTPQStep;/**> Name of the Time_TPQ_Step file.*/
const char* cFileNameNormRand;/**> Name of the NormRand file.*/
const char* cFileNameFlctRand;/**> Name of the Flct file.*/
const char* cFileName1BGreen_TPQ;/**> Name of the file to output one-body Green's functions for TPQ calculation.*/
const char* cFileName2BGreen_TPQ;/**> Name of the file to output two-body Green's functions for TPQ calculation.*/
const char* cFileNameOutputVector;/**> Name of the file to output TPQ vector.*/
const char* cFileNameInputVector;/**> Name of the file to input TPQ vector.*/

//For Time evolution
const char* cFileNameTEStep; /**> Name of the Time_TE_Step file.*/

//For FullDiag
const char* cFileNamePhys_FullDiag;/**> Name of the file to output physical values for canonical ensemble.*/
const char* cFileNamePhys_FullDiag_GC;/**> Name of the file to output physical values for grand canonical ensemble.*/
const char* cFileName1BGreen_FullDiag;/**> Name of the file to output one-body Green's functions for Full diagonalization.*/
const char* cFileName2BGreen_FullDiag;/**> Name of the file to output two-body Green's functions for Full diagonalization.*/
const char* cFileNamePhys_FullDiag_Ham;/**> Name of the file to output Hamiltonian for Full diagonalization.*/

//For Spectrum
const char* cFileNameOutputRestartVec;/**> Name of the file to output restart vector for spectrum mode.*/

//For Error
const char* cFileNameErrorSz;/**> Name of the file to output Error message for Sz calculation.*/

//For Timer
double *Timer; /**> The procedure execution time.*/
double *TimerStart;/**> Timer when the procedure starts.*/

/********************************************************************/
/********************************************************************/
double eps; /**> epsilon used in getting spectrum by Lanczos method and Lanczos eigenvector by CG method.*/
double eps_CG;/**> epsilon used in getting Lanczos eigenvector by CG method.*/
double eps_Lanczos;/**> epsilon used in LOBPCG, BiCG and Lanczos eigen value.*/
double eps_Energy;/**> epsilon for energy*/
double eps_CheckImag0;/**> epsilon for checking values of one-body and two-body interactions.*/

/*
 Variables for the MPI parallelism
*/
int nproc;//!< Number of processors, defined in InitializeMPI()
int myrank;//!< Process ID, defined in InitializeMPI()
int nthreads;//!< Number of Threads, defined in InitializeMPI()
FILE *stdoutMPI;/**<@brief File pointer to the standard output
                defined in InitializeMPI()*/

#endif /* HPHI_GLOBAL_H */

/**
@page page_variable Global variables and Data structure

HPhi uses global variables. List of them can be found in global.h

Sometimes, we pass variables function as
@code
func(&(X.Bind.Def))
@endcode
This C-structure can be found in struct.h
*/