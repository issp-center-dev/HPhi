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

/**
 * @file   ProgressMessage.c
 *
 * @brief  File for defining progress messages.
 *
 *
 */

#include "ProgressMessage.h"
const char* cProFinishDefFiles ="\n######  Definition files are correct.  ######\n\n";
const char* cProFinishDefCheck ="\n######  Indices and Parameters of Definition files(*.def) are complete.  ######\n\n";
const char* cProFinishAlloc ="\n######  LARGE ALLOCATE FINISH !  ######\n\n";
const char* cProStartCalcSz = "  Start: Calculate HilbertNum for fixed Sz. \n";
const char* cProEndCalcSz = "  End  : Calculate HilbertNum for fixed Sz. \n\n";
const char* cProStartOutputList = "  Start: output list. \n";
const char* cProEndOutputList = "  End  : output list. \n";
const char* cProStartCalcDiag = "  Start: Calculate diagaonal components of Hamiltonian. \n";
const char* cProEndCalcDiag = "  End  : Calculate diagaonal components of Hamiltonian. \n\n";

const char* cLogLanczos_EigenValueStart="######  Start: Calculate Lanczos Eigenvalue.  ######\n\n";
const char* cLogLanczos_EigenValueEnd="\n######  End  : Calculate Lanczos EigenValue.  ######\n\n";
const char* cLogLanczos_EigenVecStart="######  Start: Calculate Lanczos Eigenvec.  ######\n\n";
const char* cLogLanczos_EigenVecEnd="\n######  End  : Calculate Lanczos EigenVec.  ######\n\n";


const char* cLogTPQ_Start="######  Start: TPQCalculation.  ######\n\n";
const char* cLogTPQ_End="######  End  : TPQCalculation.  ######\n\n";


const char* cLogFullDiag_SetHam_Start="######  Start: Setting Hamiltonian.  ######\n\n";
const char* cLogFullDiag_SetHam_End="######  End  : Setting Hamiltonian.  ######\n\n";

const char* cLogFullDiag_InputHam_Start="######  Start: Input Hamiltonian.  ######\n\n";
const char* cLogFullDiag_InputHam_End="######  End  : Input Hamiltonian.  ######\n\n";


const char* cLogFullDiag_OutputHam_Start="######  Start: Output Hamiltonian.  ######\n\n";
const char* cLogFullDiag_OutputHam_End="######  End  : Output Hamiltonian.  ######\n\n";


const char* cLogFullDiag_Start="######  Start: Diagonalization.  ######\n\n";
const char* cLogFullDiag_End="######  End  : Diagonalization.  ######\n\n";

const char* cLogFullDiag_ExpecValue_Start="######  Start: Calc Expected value.  ######\n\n";
const char* cLogFullDiag_ExpecValue_End="######  End  : Calc Expected value.  ######\n\n";
const char* cLogFinish ="######  Finish Calculation.  ######\n";

const char* cLogTEM_Start="######  Start: TimeEvolution by Taylor Expansion.  ######\n\n";
const char* cLogTEM_End="######  End  : TimeEvolution by Taylor Expansion.  ######\n\n";

