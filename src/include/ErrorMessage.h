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

#ifndef HPHI_ERRORMESSAGE_H
#define HPHI_ERRORMESSAGE_H

//! Error Message in HPhiMain.c
int iErrCodeMem;

char *cErrNameList;
char *cErrOutput;
char *cErrDefFile;
char *cErrNvec;
char *cErrnvec;
char *cErrnvecShow;
char *cErrIndices;
char *cErrArgv;

//! Error Message in HPhiTrans.c
char *cErrTransfer;
char *cErrDoubleCounting;
char *cErrChemicalPotential;
char *cErrLargeMem;


//! Error Message in readdef.c
char *cErrReadDefFile;
char *cErrDefFileFormat;
char *cErrNLoc;
char *cErrDefFileParam;
char *cErrCalcType;
char *cErrOutputMode;
char *cErrCalcModel;
char *cErrCalcEigenVec;
char *cErrSetIniVec;
char *cErrOutputHam;
char *cErrInputHam;
char *cErrInputOutputHam;
char *cErrOutputHamForFullDiag;
char *cErrRestart;
char *cErrFiniteTemp;
char *cErrCUDA;
char *cErrScaLAPACK;

char *cErrKW;
char *cErrKW_ShowList;
char *cErrKW_Same;
char *cErrKW_InCorPair;

char *cErrNsite;
char *cErrNcond;
char *cErrNumAve;
char *cErrExpecInterval;
char *cErrLanczos_max;
char *cErrLanczos_eps;
char *cErrLanczosTarget;
char *cErrLanczosExct;

char *cErrMakeDef;
char *cErrIncorrectDef;
char *cErrNonHermiteTrans;
char *cErrNonHermiteTransForAll;
char *cErrNonHermiteInterAll;
char *cErrNonConservedInterAll;
char *cErrNonHermiteInterAllForAll;
char *cErrIncorrectFormatForKondoInt;
char *cErrIncorrectFormatForKondoTrans;
char *cErrIncorrectFormatInter;
char *cErrIncorrectSpinIndexForInter;
char *cErrIncorrectSpinIndexForTrans;

char *cErrIncorrectFormatForSpinTrans;
char *cWarningIncorrectFormatForSpin;
char *cWarningIncorrectFormatForSpin2;


//! Error Message in CheckMPI.c
char *cErrNProcNumberHubbard;
char *cErrNProcNumberSpin;
char *cErrNProcNumberGneralSpin;
char *cErrNProcNumber;
char *cErrNProcNumberSet;

//! Error Message in diagonal calc.c
char *cErrNoModel;
char *cErrNoHilbertSpace;

//! Error Message in bitcalc.c
char *cErrSiteNumber;

//! Error Message in mltiply.c
char *cErrMltiply;


//! Error Message in FileIO.c
char *cErrFIOpen;

//! Error Message in sz.c
char* cErrSz;
char* cErrSz_NoFile;
char* cErrSz_NoFile_Show;
char* cErrSz_ShowDim;
char* cErrSz_OutFile;

#endif /* HPHI_ERRORMESSAGE_H */
