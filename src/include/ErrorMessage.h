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
extern int iErrCodeMem;

extern char *cErrNameList;
extern char *cErrOutput;
extern char *cErrDefFile;
extern char *cErrNvec;
extern char *cErrnvec;
extern char *cErrnvecShow;
extern char *cErrIndices;
extern char *cErrArgv;

//! Error Message in HPhiTrans.c
extern char *cErrTransfer;
extern char *cErrDoubleCounting;
extern char *cErrChemicalPotential;
extern char *cErrLargeMem;


//! Error Message in readdef.c
extern char *cErrReadDefFile;
extern char *cErrDefFileFormat;
extern char *cErrNLoc;
extern char *cErrDefFileParam;
extern char *cErrCalcType;
extern char *cErrOutputMode;
extern char *cErrCalcModel;
extern char *cErrCalcEigenVec;
extern char *cErrSetIniVec;
extern char *cErrOutputHam;
extern char *cErrInputHam;
extern char *cErrInputOutputHam;
extern char *cErrOutputHamForFullDiag;
extern char *cErrRestart;
extern char *cErrFiniteTemp;
extern char *cErrCUDA;
extern char *cErrScaLAPACK;

extern char *cErrKW;
extern char *cErrKW_ShowList;
extern char *cErrKW_Same;
extern char *cErrKW_InCorPair;

extern char *cErrNsite;
extern char *cErrNcond;
extern char *cErrNumAve;
extern char *cErrExpecInterval;
extern char *cErrLanczos_max;
extern char *cErrLanczos_eps;
extern char *cErrLanczosTarget;
extern char *cErrLanczosExct;

extern char *cErrMakeDef;
extern char *cErrIncorrectDef;
extern char *cErrNonHermiteTrans;
extern char *cErrNonHermiteTransForAll;
extern char *cErrNonHermiteInterAll;
extern char *cErrNonConservedInterAll;
extern char *cErrNonHermiteInterAllForAll;
extern char *cErrIncorrectFormatForKondoInt;
extern char *cErrIncorrectFormatForKondoTrans;
extern char *cErrIncorrectFormatInter;
extern char *cErrIncorrectSpinIndexForInter;
extern char *cErrIncorrectSpinIndexForTrans;

extern char *cErrIncorrectFormatForSpinTrans;
extern char *cWarningIncorrectFormatForSpin;
extern char *cWarningIncorrectFormatForSpin2;


//! Error Message in CheckMPI.c
extern char *cErrNProcNumberHubbard;
extern char *cErrNProcNumberSpin;
extern char *cErrNProcNumberGneralSpin;
extern char *cErrNProcNumber;
extern char *cErrNProcNumberSet;

//! Error Message in diagonal calc.c
extern char *cErrNoModel;
extern char *cErrNoHilbertSpace;

//! Error Message in bitcalc.c
extern char *cErrSiteNumber;

//! Error Message in mltiply.c
extern char *cErrMltiply;


//! Error Message in FileIO.c
extern char *cErrFIOpen;

//! Error Message in sz.c
extern char* cErrSz;
extern char* cErrSz_NoFile;
extern char* cErrSz_NoFile_Show;
extern char* cErrSz_ShowDim;
extern char* cErrSz_OutFile;

#endif /* HPHI_ERRORMESSAGE_H */
