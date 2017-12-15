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
const char* cFileName1BGreen_TE="%s_cisajs_step%d.dat";
const char* cFileName2BGreen_TE="%s_cisajscktalt_step%d.dat";
const char* cFileNameOutputVector="tmpvec_set%d_rank_%d.dat";
const char* cFileNameInputVector="tmpvec_set%d_rank_%d.dat";

//Fot Time evolution
const char* cFileNameTEStep="%s_Time_TE_Step.dat";

//For FullDiag
const char* cFileNamePhys_FullDiag="%s_phys_Nup%d_Ndown%d.dat";
const char* cFileNamePhys_FullDiag_GC="%s_phys.dat";
const char* cFileName1BGreen_FullDiag="%s_cisajs_eigen%d.dat";
const char* cFileName2BGreen_FullDiag="%s_cisajscktalt_eigen%d.dat";
const char* cFileNamePhys_FullDiag_Ham="%s_Ham.dat";

//For Spectrum
const char* cFileNameOutputRestartVec="%s_recalcvec_rank_%d.dat";

//For Error
const char* cFileNameErrorSz="Err_sz.dat";
