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

extern double complex *v0;  /**< A vector after multiplying Hamiltonian, @f$ v_0 = H v_1@f$.*/
extern double complex *v1;  /**< A vector before multiplying Hamiltonian, @f$ v_0 = H v_1@f$.*/
extern double complex *v2;  /**< A temporary vector for time evolution calculation, @f$ v2 = H*v1 = H^coef |psi(t)>@f$.*/
extern double complex *v1buf; /**< A temporary vector for MPI. */

//[s] For calcSpectrum
extern double complex *v1Org; /**< An input vector to calculate spectrum function.*/
extern double complex *vg; /**< A vector used in the CG mode.*/
//[e] For calcSpectrum

extern double *alpha,*beta; /**< Tridiagonal components used in Lanczos mode.*/
extern double complex **vec; /**< Eigen vectors.*/
extern double *list_Diagonal; /**< list for diagonal components.*/
extern long unsigned int *list_1; /**< list of getting real-space configuration for canonical state*/
extern long unsigned int *list_1buf;/**< list of getting real-space configuration for canonical state across processes*/
extern long unsigned int *list_2_1;/**< list to get index of list_1*/
extern long unsigned int *list_2_2;/**< list to get index of list_1*/

/*[s] For Spectrum */
extern long unsigned int *list_1_org; /**< list of getting real-space configuration for canonical state before excitation*/
extern long unsigned int *list_1buf_org;/**< list of getting real-space configuration for canonical state before excitation across processes*/
extern long unsigned int *list_2_1_org;/**< list to get index of list_1_org*/
extern long unsigned int *list_2_2_org;/**< list to get index of list_1_org*/
/*[e] For Spectrum */

/*[s] For Lanczos */
extern int     initial_mode;/**< mode to get initial state (0: use same random generator for MPI, 1: use each random generator for MPI)*/
/*[e] For Lanczos */

/*[s] For TPQ*/
extern double LargeValue;/**< constant value l for TPQ calculation.*/
extern int    NumAve;/**< Average number for TPQ calculation*/
extern int step_i;/**< step for TPQ calculation*/
extern double global_norm;/**< norm before normalization for TPQ calculation*/
extern double global_1st_norm;/**< 1-st norm for TPQ calculation*/
extern int step_spin;/**< output step for TE calculation.*/
/*[e] For TPQ*/

/*[s] For All Diagonalization*/
extern double complex**Ham; /**> Hamiltonian for full diagonalization. */
extern double complex **L_vec;/**> eigen vectors*/
#ifdef _SCALAPACK
extern double complex *Z_vec; /**> distributed matrix of eigen vector*/
extern int descZ_vec[9]; /*descriptor for Z_vec*/
#endif
/*[e] For All Diagonalization*/

extern const char* cParentOutputFolder; /**> Path to output results*/

//For TimeKeep
extern const char* cFileNameTimeKeep; /**> Name of the file to output calculation processing time*/
extern const char* cFileNameSzTimeKeep;/**> Name of the file to output sz calculation time*/

//For Check
extern const char* cFileNameCheckCoulombIntra;/**> Name of the file to check Coulomb Intra interactions.*/
extern const char* cFileNameCheckChemi;/**> Name of the file to check Chemical potentials.*/
extern const char* cFileNameCheckInterU;/**> Name of the file to check Inter U.*/
extern const char* cFileNameCheckHund;/**>Name of the file to check Hund interactions. */
extern const char* cFileNameCheckInterAll;/**> Name of the file to check InterAll interactions.*/
extern const char* cFileNameCheckMemory;/**> Name of the file to check memory.*/
extern const char* cFileNameCheckSdim;/**> Name of the file to check dimension.*/

//For EDTrans
extern const char* cFileNameWarningOnTransfer; /**> Name of the file to output warning of transfer integrals.*/

//For Lanczos
extern const char* cFileNameLanczosStep;/**> Name of the file to output Lanczos step.*/
extern const char* cFileNameEnergy_Lanczos;/**> Name of the file to output energies.*/
extern const char* cFileNameEigenvalue_Lanczos;/**> Name of the file to output eigen values.*/
extern const char* cFileNameEnergy_CG;/**> Name of the file to output energies obtained by CG method.*/
extern const char* cFileName1BGreen_Lanczos;/**> Name of the file to output One-Body Green's functions obtained by Lanczos method.*/
extern const char* cFileName1BGreen_CG;/**> Name of the file to output One-Body Green's functions obtained by CG method.*/
extern const char* cFileName2BGreen_Lanczos;/**> Name of the file to output Two-Body Green's functions obtained by Lanczos method.*/
extern const char* cFileName3BGreen_Lanczos;/**> Name of the file to output Three-Body Green's functions obtained by Lanczos method.*/
extern const char* cFileName4BGreen_Lanczos;/**> Name of the file to output Four-Body Green's functions obtained by Lanczos method.*/
extern const char* cFileName6BGreen_Lanczos;/**> Name of the file to output Six-Body Green's functions obtained by Lanczos method.*/
extern const char* cFileName2BGreen_CG;/**> Name of the file to output Two-Body Green's functions obtained by CG method.*/
extern const char* cFileNameTimeEV_CG;/**> Name of the file to output time for getting eigen vector by CG method.*/
extern const char* cFileNameListModel;/**> Name of the file to output list.*/
extern const char* cFileNameOutputEigen;/**> Name of the file to output eigen vector.*/
extern const char* cFileNameInputEigen;/**> Name of the file to input eigen vector.*/
extern const char* cFileNameCalcDynamicalGreen;/**> Name of the file to output dynamical Green's function.*/
extern const char* cFileNameTridiagonalMatrixComponents;/**> Name of the file to output tridiagonal matrix components.*/

//For TPQ
extern const char* cFileNameSSRand;/**> Name of the SS_rand file.*/
extern const char* cFileNameTPQStep;/**> Name of the Time_TPQ_Step file.*/
extern const char* cFileNameNormRand;/**> Name of the NormRand file.*/
extern const char* cFileNameFlctRand;/**> Name of the Flct file.*/
extern const char* cFileName1BGreen_TPQ;/**> Name of the file to output one-body Green's functions for TPQ calculation.*/
extern const char* cFileName2BGreen_TPQ;/**> Name of the file to output two-body Green's functions for TPQ calculation.*/
extern const char* cFileName3BGreen_TPQ;/**> Name of the file to output three-body Green's functions for TPQ calculation.*/
extern const char* cFileName4BGreen_TPQ;/**> Name of the file to output four-body Green's functions for TPQ calculation.*/
extern const char* cFileName6BGreen_TPQ;/**> Name of the file to output six-body Green's functions for TPQ calculation.*/
extern const char* cFileName1BGreen_TE;/**> Name of the file to output one-body Green's functions for Time Evolution calculation.*/
extern const char* cFileName2BGreen_TE;/**> Name of the file to output two-body Green's functions for Time Evolution calculation.*/
extern const char* cFileName3BGreen_TE;/**> Name of the file to output three-body Green's functions for Time Evolution calculation.*/
extern const char* cFileName4BGreen_TE;/**> Name of the file to output four-body Green's functions for Time Evolution calculation.*/
extern const char* cFileName6BGreen_TE;/**> Name of the file to output six-body Green's functions for Time Evolution calculation.*/
extern const char* cFileNameOutputVector;/**> Name of the file to output TPQ vector.*/
extern const char* cFileNameInputVector;/**> Name of the file to input TPQ vector.*/

//For Time evolution
extern const char* cFileNameTEStep; /**> Name of the Time_TE_Step file.*/
extern const char* cFileNameSS;/**> Name of the SS file.*/
extern const char* cFileNameNorm;/**> Name of the Norm file.*/
extern const char* cFileNameFlct;/**> Name of the Flct file.*/

//For FullDiag
extern const char* cFileNamePhys_FullDiag;/**> Name of the file to output physical values for canonical ensemble.*/
extern const char* cFileNamePhys_FullDiag_GC;/**> Name of the file to output physical values for grand canonical ensemble.*/
extern const char* cFileName1BGreen_FullDiag;/**> Name of the file to output one-body Green's functions for Full diagonalization.*/
extern const char* cFileName2BGreen_FullDiag;/**> Name of the file to output two-body Green's functions for Full diagonalization.*/
extern const char* cFileName3BGreen_FullDiag;/**> Name of the file to output three-body Green's functions for Full diagonalization.*/
extern const char* cFileName4BGreen_FullDiag;/**> Name of the file to output four-body Green's functions for Full diagonalization.*/
extern const char* cFileName6BGreen_FullDiag;/**> Name of the file to output six-body Green's functions for Full diagonalization.*/
extern const char* cFileNamePhys_FullDiag_Ham;/**> Name of the file to output Hamiltonian for Full diagonalization.*/

//For Spectrum
extern const char* cFileNameOutputRestartVec;/**> Name of the file to output restart vector for spectrum mode.*/
extern const char* cFileNameOutputExcitedVec;/**> Name of the file to output the excited vector for spectrum mode.*/

//For Error
extern const char* cFileNameErrorSz;/**> Name of the file to output Error message for Sz calculation.*/

//For Timer
extern double *Timer; /**> The procedure execution time.*/
extern double *TimerStart;/**> Timer when the procedure starts.*/

/********************************************************************/
/********************************************************************/
extern double eps; /**> epsilon used in getting spectrum by Lanczos method and Lanczos eigenvector by CG method.*/
extern double eps_CG;/**> epsilon used in getting Lanczos eigenvector by CG method.*/
extern double eps_Lanczos;/**> epsilon used in LOBPCG, BiCG and Lanczos eigen value.*/
extern double eps_Energy;/**> epsilon for energy*/
extern double eps_CheckImag0;/**> epsilon for checking values of one-body and two-body interactions.*/

/*
 Variables for the MPI parallelism
*/
extern int nproc;//!< Number of processors, defined in InitializeMPI()
extern int myrank;//!< Process ID, defined in InitializeMPI()
extern int nthreads;//!< Number of Threads, defined in InitializeMPI()
extern FILE *stdoutMPI;/**<@brief File pointer to the standard output
                defined in InitializeMPI()*/

#endif /* HPHI_GLOBAL_H */

/**
@page page_variable Global variables and Data structure

In HPhi, global variables are used. List of them can be found in global.h

Sometimes, we pass variables to the function as
@code
func(&(X.Bind.Def))
@endcode
This C-structure is defined in struct.h.
*/

