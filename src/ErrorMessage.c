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

#include "ErrorMessage.h"

//! Error Message in HPhiMain.c
int iErrCodeMem=1;

char *cErrNameList="Error: *.out(*.exe) NameListFile [OptParaFile]\n";
char *cErrOutput="Error: Fail to make output folder in current directory. \n";
char *cErrDefFile="Error: Definition files(*.def) are incomplete.\n";
char *cErrNvec="Error: nvec should be smaller than exct are incorrect.\n";
char *cErrnvec="Error: nvec should be smaller than exct are incorrect.\n";
char *cErrnvecShow="Error: nvec = %d, exct=%d.\n";
char *cErrIndices="Error: Indices and Parameters of Definition files(*.def) are incomplete.\n";
char *cErrArgv="Error: %s is incorrect. set '-e' or '-s'.\n";

//! Error Message in HPhiTrans.c
char *cErrTransfer="WarningOnTransfer.dat";
char *cErrDoubleCounting="double conuntings in transfers: i=%d j=%d spni %d spnj %d  \n";
char *cErrChemicalPotential="possible error in chemical potentaial: usually, # of chemical potential = 2*(# of sites)  \n";
char *cErrLargeMem="Error: Fail for memory allocation. see error code %d in user manual. \n";


//! Error Message in readdef.c
char *cErrReadDefFile="Error: %s (Broken file or Not exist)\n";
char *cErrDefFileFormat="Error: incorrect format= %s. \n";
char *cErrNLoc ="Error: Ne=Nup+Ndown must be (Ne >= NLocalSpin).\n";
char *cErrDefFileParam="Error: In %s, wrong parameter name:%s \n";
char *cErrCalcType="Error in %s\n CalcType: 0: Lanczos Method, 1: Thermal Pure Quantum State Method, 2: Full Diagonalization Method, 3: Calculation Spectrum mode.\n";
char *cErrOutputMode="Error in %s\n OutputMode: \n 0: calc one body green function and two body green functions,\n 1: calc one body green function and two body green functions and correlatinos for charge and spin.\n";
char *cErrCalcEigenVec="Error in %s\n CalcEigenVec: \n 0: Lanczos+CG method,\n 1: Lanczos method.\n";
char *cErrOutputHam="Error in %s\n OutputHam: \n 0: not output Hamiltonian,\n 1: output Hamiltonian.\n";
char *cErrOutputHamForFullDiag="Error in %s\n OutputHam is only defined for FullDiag mode, CalcType=2.\n";
char *cErrInputHam="Error in %s\n InputHam: \n 0: not input Hamiltonian,\n 1: input Hamiltonian.\n";
char *cErrInputOutputHam="Error in %s\n OutputHam=1 and InputHam=1.\n";
char *cErrCalcModel="Error in %s\n CalcModel: \n 0: Hubbard, 1: Spin, 2: Kondo, 3: HubbardGC, 4: SpinGC, 5:KondoGC.\n";
char *cErrFiniteTemp="Error in %s\n FlgFiniteTemperature: Finite Temperature, 1: Zero Temperature.\n";
char *cErrSetIniVec="Error in %s\n InitialVecType: \n 0: complex type,\n 1: real type.\n";
char *cErrRestart="Error in %s\n Restart: \n 0: not restart (default).\n 1: output a restart vector.\n 2: input a restart vector and output a new restart vector.\n 3: input a restart vector.\n";
char *cErrCUDA="Error in %s\n NGPU: NGPU must be greater than 0.\n";
char *cErrScaLAPACK="Error in %s\n ScaLAPACK: \n 0: Use LAPACK for FullDiag mode,\n 1: Use ScaLAPACK for FullDiag mode.\n";

char *cErrNcond= "Error in %s\n Ncond must be greater than 0.\n ";
char *cErrNsite= "Error in %s\n Nsite must be positive value.\n ";
char *cErrNumAve= "Error in %s\n NumAve must be positive value.\n ";
char *cErrExpecInterval= "Error in %s\n ExpecInterval must be positive value.\n ";
char *cErrLanczos_max="Error in %s\n Lanczos_max must be positive value.\n ";
char *cErrLanczos_eps="Error in %s\n Lanczos_eps must be positive value.\n ";
char *cErrLanczosExct="Error in %s\n exct=%d must be greater than 1 and smaller than nvec=%d.\n ";
char *cErrLanczosTarget="Error in %s\n LanczosTarget=%d must be greater than exct=%d.\n ";


char *cErrKW="Error: Wrong keywords '%s' in %s.\n";
char *cErrKW_ShowList="Choose Keywords as follows: \n";
char *cErrKW_Same="Error: Same keywords exist in %s.\n";
char *cErrKW_InCorPair="Error: keyword and filename must be set as a pair in %s.\n";

char *cErrMakeDef="Error: Need to make a def file for %s.\n";
char *cErrIncorrectDef="Error: incorrect file.\n";
char *cErrNonHermiteTrans="Error: NonHermite (i, spni, j, spnj) = (%d,  %d, %d, %d), trans_re= %lf, trans_im= %lf.\n";
char *cErrNonHermiteTransForAll= "Error: NonHermite Pair exists in transfer integral. \n";
char *cErrNonHermiteInterAll="Error: NonHermite (i, spni, j, spnj, k, spnk, l, spnl) = (%d, %d, %d, %d, %d, %d, %d, %d), InterAll_re= %lf, InterAll_im= %lf . \n";
char *cErrNonConservedInterAll="Error: This operator breaks Sz Symmetry (i, spni, j, spnj, k, spnk, l, spnl) = (%d, %d, %d, %d, %d, %d, %d, %d), InterAll_re= %lf, InterAll_im= %lf . \n";
char *cErrNonHermiteInterAllForAll= "Error: NonHermite Pair exists in InterAll. \n";
char *cErrIncorrectFormatForKondoInt= "Error: Site component of (i, j, k, l) =(%d, %d, %d, %d) is incorrect.\n";
char *cErrIncorrectFormatForKondoTrans= "Error: Site component of (i, j) =(%d, %d) is incorrect.\n";


char *cErrIncorrectFormatForSpinTrans= "Error: Site component of (i, j) =(%d, %d) is incorrect.\n";
char *cWarningIncorrectFormatForSpin= "Warning: Site component of (i, j, k, l) =(%d, %d, %d, %d) is not correct; i=j and k=l must be satisfied. \n";
char *cWarningIncorrectFormatForSpin2= "Warning: Site component of (i, j) =(%d, %d) is ignored.\n";
char *cErrIncorrectFormatInter= "Error: Use only InterAll for setteing interactions for general spin.\n";
char *cErrIncorrectSpinIndexForInter="Error: Spin index is incorrect for interactions defined in InterAll file.\n";
char *cErrIncorrectSpinIndexForTrans="Error: Spin index is incorrect for transfers defined in Trans file.\n";

//! Error Message in CheckMPI.c
char *cErrNProcNumberHubbard = "Error ! The number of PROCESS should be 4-exponent !\n";
char *cErrNProcNumberSpin = "Error ! The number of PROCESS should be 2-exponent !\n";
char *cErrNProcNumberGneralSpin = "Error ! The number of PROCESS is wrong !\n";
char *cErrNProcNumber = "        The number of PROCESS : %d\n";
char *cErrNProcNumberSet = "        Set the number of PROCESS as %d or %d.\n";


//! Error Message in diagonal calc.c
char *cErrNoModel ="Error: CalcModel %d is incorrect.\n";
char *cErrNoHilbertSpace = "Error: There is no target of Hilbert space.\n";

//! Error Message in bitcalc.c
char *cErrSiteNumber = "Error: Total Site Number is incorrect.\n";


//! Error Message in mltiply.c
char *cErrMltiply="ERROR IN mltiply.c \n";


//! Error Message in FileIO.c
char *cErrFIOpen ="FileOpenError: %s.\n";

//! Error Message in sz.c
char* cErrSz="Error: in sz. \n";
char* cErrSz_NoFile="No file. Please set READ=0.\n";
char* cErrSz_NoFile_Show=" %s does not exist. \n";
char* cErrSz_ShowDim="imax = %ld, Check.idim_max=%ld \n";
char* cErrSz_OutFile="Caution!!  Error in sz !!!! idim_max is not correct \n";
