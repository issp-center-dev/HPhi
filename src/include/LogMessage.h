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
extern const char* cReadFileNamelist;
extern const char* cReadFile;
extern const char* cReadDefStart;
extern const char* cReadDefFinish;

//sz.c
extern const char* cStateLocSpin;
extern const char* cStateNupNdown;
extern const char* cInitalSz;
extern const char* cOMPSzStart;
extern const char* cOMPSzMid;
extern const char* cOMPSzFinish;
extern const char* cReadSzStart;
extern const char* cReadSzEnd;

//CalcByLanczos.c
extern const char* cReadEigenVecStart;
extern const char* cReadEigenVecFinish;
extern const char* cOutputEigenVecStart;
extern const char* cOutputEigenVecFinish;


//CG_EigenVector.c
extern const char* cLogCG_EigenVecStart;
extern const char* cLogCG_EigenVecEnd;
extern const char* cCG_EigenVecStart;
extern const char* cCG_EigenVecFinish;

//diagonalcalc.c
extern const char* cDiagonalCalcStart;
extern const char* cDiagonalCalcFinish;

//calcspectrum.c
extern const char* c_InputEigenVectorStart;
extern const char* c_InputEigenVectorEnd;
extern const char* c_CalcExcitedStateStart;
extern const char* c_CalcExcitedStateEnd;
extern const char* c_CalcSpectrumStart;
extern const char* c_CalcSpectrumEnd;
extern const char* c_GetTridiagonalStart;
extern const char* c_GetTridiagonalEnd;
extern const char* c_CalcSpectrumFromTridiagonalStart;
extern const char* c_CalcSpectrumFromTridiagonalEnd;
extern const char* c_OutputSpectrumRecalcvecStart;
extern const char* c_OutputSpectrumRecalcvecEnd;
extern const char* c_InputSpectrumRecalcvecStart;
extern const char* c_InputSpectrumRecalcvecEnd;

//calcspectrum in Lanczos_Eigenvalue.c
extern const char* c_Lanczos_SpectrumStep;

//FirstMultiply.c, Multiply.c
extern const char* cTPQStep;
extern const char* cTPQStepEnd;

//CalcByTEM.c
extern const char* cTEStep;
extern const char* cTEStepEnd;

//Lanczos_EigenValue.c
extern const char* cLogLanczos_EigenValueNotConverged;
extern const char* cLanczos_EigenValueStart;
extern const char* cLanczos_EigenValueStep;
extern const char* cLanczos_EigenValueFinish;

//Lanczos_EigenVector.c
extern const char* cLogLanczos_EigenVectorStart;
extern const char* cLogLanczos_EigenVectorEnd;
extern const char* cLanczos_EigenVectorStart;
extern const char* cLanczos_EigenVectorFinish;

//expec.c
extern const char* cExpecStart;
extern const char* cExpecEnd;
extern const char* cTPQExpecStart;
extern const char* cTPQExpecEnd;

//expec_cisajs.c
extern const char* cLogLanczosExpecOneBodyGStart;
extern const char* cLogCGExpecOneBodyGStart;
extern const char* cLogLanczosExpecOneBodyGEnd;
extern const char* cLogCGExpecOneBodyGEnd;

extern const char* cLanczosExpecOneBodyGFinish;
extern const char* cLanczosExpecOneBodyGStart;
extern const char* cTPQExpecOneBodyGStart;
extern const char* cTPQExpecOneBodyGFinish;
extern const char* cCGExpecOneBodyGStart;
extern const char* cCGExpecOneBodyGFinish;
extern const char* cTEExpecOneBodyGStart;
extern const char* cTEExpecOneBodyGFinish;

//expec_cisajucktaltdc.c
extern const char*  cLogLanczosExpecTwoBodyGStart;
extern const char*  cLogLanczosExpecTwoBodyGFinish;
extern const char*  cLanczosExpecTwoBodyGStart;
extern const char*  cLogCGExpecTwoBodyGFinish;
extern const char*  cLanczosExpecTwoBodyGFinish;
extern const char*  cCGExpecTwoBodyGStart;
extern const char*  cCGExpecTwoBodyGFinish;
extern const char*  cTPQExpecTwoBodyGStart;
extern const char*  cTPQExpecTwoBodyGFinish;
extern const char*  cTEExpecTwoBodyGStart;
extern const char*  cTEExpecTwoBodyGFinish;

//expec_energy.c
extern const char* cLogExpecEnergyStart;
extern const char* cLogExpecEnergyEnd;

//CalcByTPQ.c
extern const char* cLogTPQRand;
extern const char* cLogSSRand;
extern const char* cLogNormRand;
extern const char* cLogFlctRand;
extern const char* cLogTPQStep;
extern const char* cLogTPQEnd;

//CalcByTEM.c
extern const char* cLogTEStep;
extern const char* cLogSS;
extern const char* cLogNorm;
extern const char* cLogFlct;

extern const char* cLogInputVecStart;
extern const char* cLogInputVecFinish;
extern const char* cLogOutputVecStart;
extern const char* cLogOutputVecFinish;
extern const char* cOutputVecStart;
extern const char* cOutputVecFinish;

//FirstMultiply.c
extern const char* cLogCheckInitComplex;
extern const char* cLogCheckInitReal;
#endif /* HPHI_LOGMESSAGE_H */
