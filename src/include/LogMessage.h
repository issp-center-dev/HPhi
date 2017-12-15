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

#ifndef HPHI_LOGMESSAGE_H
#define HPHI_LOGMESSAGE_H

//readdef.c
const char* cReadFileNamelist;
const char* cReadFile;
const char* cReadDefStart;
const char* cReadDefFinish;

//sz.c
const char* cStateLocSpin;
const char* cStateNupNdown;
const char* cInitalSz;
const char* cOMPSzStart;
const char* cOMPSzMid;
const char* cOMPSzFinish;
const char* cReadSzStart;
const char* cReadSzEnd;

//CalcByLanczos.c
const char* cReadEigenVecStart;
const char* cReadEigenVecFinish;
const char* cOutputEigenVecStart;
const char* cOutputEigenVecFinish;


//CG_EigenVector.c
const char* cLogCG_EigenVecStart;
const char* cLogCG_EigenVecEnd;
const char* cCG_EigenVecStart;
const char* cCG_EigenVecFinish;

//diagonalcalc.c
const char* cDiagonalCalcStart;
const char* cDiagonalCalcFinish;

//calcspectrum.c
const char* c_InputEigenVectorStart;
const char* c_InputEigenVectorEnd;
const char* c_CalcExcitedStateStart;
const char* c_CalcExcitedStateEnd;
const char* c_CalcSpectrumStart;
const char* c_CalcSpectrumEnd;
const char* c_GetTridiagonalStart;
const char* c_GetTridiagonalEnd;
const char* c_CalcSpectrumFromTridiagonalStart;
const char* c_CalcSpectrumFromTridiagonalEnd;
const char* c_OutputSpectrumRecalcvecStart;
const char* c_OutputSpectrumRecalcvecEnd;
const char* c_InputSpectrumRecalcvecStart;
const char* c_InputSpectrumRecalcvecEnd;

//calcspectrum in Lanczos_Eigenvalue.c
const char* c_Lanczos_SpectrumStep;


//FirstMultiply.c, Multiply.c
const char* cTPQStep;
const char* cTPQStepEnd;

//CalcByTEM.c
const char* cTEStep;
const char* cTEStepEnd;


//Lanczos_EigenValue.c
const char* cLogLanczos_EigenValueNotConverged;
const char* cLanczos_EigenValueStart;
const char* cLanczos_EigenValueStep;
const char* cLanczos_EigenValueFinish;

//Lanczos_EigenVector.c
const char* cLogLanczos_EigenVectorStart;
const char* cLogLanczos_EigenVectorEnd;
const char* cLanczos_EigenVectorStart;
const char* cLanczos_EigenVectorFinish;

//expec.c
const char* cExpecStart;
const char* cExpecEnd;
const char* cTPQExpecStart;
const char* cTPQExpecEnd;

//expec_cisajs.c
const char* cLogLanczosExpecOneBodyGStart;
const char* cLogCGExpecOneBodyGStart;
const char* cLogLanczosExpecOneBodyGEnd;
const char* cLogCGExpecOneBodyGEnd;

const char* cLanczosExpecOneBodyGFinish;
const char* cLanczosExpecOneBodyGStart;
const char* cTPQExpecOneBodyGStart;
const char* cTPQExpecOneBodyGFinish;
const char* cCGExpecOneBodyGStart;
const char* cCGExpecOneBodyGFinish;
const char* cTEExpecOneBodyGStart;
const char* cTEExpecOneBodyGFinish;

//expec_cisajucktaltdc.c
const char*  cLogLanczosExpecTwoBodyGStart;
const char*  cLogLanczosExpecTwoBodyGFinish;
const char*  cLanczosExpecTwoBodyGStart;
const char*  cLogCGExpecTwoBodyGFinish;
const char*  cLanczosExpecTwoBodyGFinish;
const char*  cCGExpecTwoBodyGStart;
const char*  cCGExpecTwoBodyGFinish;
const char*  cTPQExpecTwoBodyGStart;
const char*  cTPQExpecTwoBodyGFinish;
const char*  cTEExpecTwoBodyGStart;
const char*  cTEExpecTwoBodyGFinish;

//expec_energy.c
const char* cLogExpecEnergyStart;
const char* cLogExpecEnergyEnd;

//CalcByTPQ.c
const char* cLogTPQRand;
const char* cLogSSRand;
const char* cLogNormRand;
const char* cLogFlctRand;
const char* cLogTPQStep;
const char* cLogTPQEnd;

//CalcByTEM.c
const char* cLogTEStep;
const char* cLogSS;
const char* cLogNorm;
const char* cLogFlct;

const char* cLogInputVecStart;
const char* cLogInputVecFinish;
const char* cLogOutputVecStart;
const char* cLogOutputVecFinish;
const char* cOutputVecStart;
const char* cOutputVecFinish;


//FirstMultiply.c
const char* cLogCheckInitComplex;
const char* cLogCheckInitReal;
#endif /* HPHI_LOGMESSAGE_H */
