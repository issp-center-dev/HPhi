/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima */

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

//sz.c
const char* cStateLocSpin;
const char* cInitalSz;
const char* cOMPSzStart;
const char* cOMPSzMid;
const char* cOMPSzFinish;
const char* cSingleSzStart;
const char* cSingleSzFinish;
const char* cReadSzStart;
const char* cReadSzEnd;

//CG_EigenVector.c
const char* cLogCG_EigenVecStart;
const char* cLogCG_EigenVecEnd;
const char* cCG_EigenVecFinish;

//diagonalcalc.c
const char* cDiagonalCalcFinish;


//FirstMultiply.c, Multiply.c
const char* cTPQStep;
const char* cTPQStepEnd;

//Lanczos_EigenValue.c
const char* cLogLanczos_EigenValueStart;
const char* cLogLanczos_EigenValueNotConverged;
const char* cLogLanczos_EigenValueEnd;
const char* cLanczos_EigenValueStart;
const char* cLanczos_EigenValueStep;
const char* cLanczos_EigenValueFinish;

//Lanczos_EigenVector.c
const char* cLogLanczos_EigenVectorStart;
const char* cLogLanczos_EigenVectorEnd;
const char* cLanczos_EigenVectorFinish;

//expec.c
const char* cExpecStart;
const char* cExpecEnd;

//expec_cisajs.c
const char* cLogLanczosExpecOneBodyGStart;
const char* cLogCGExpecOneBodyGStart;
const char* cLogLanczosExpecOneBodyGEnd;
const char* cLogCGExpecOneBodyGEnd;

const char* cLanczosExpecOneBodyGFinish;
const char* cTPQExpecOneBodyGStart;
const char* cTPQExpecOneBodyGFinish;
const char* cExpecOneBodyGFinish;

//expec_cisajucktaltdc.c
const char*  cLogLanczosExpecTwoBodyGFinish;
const char*  cLogCGExpecTwoBodyGFinish;
const char*  cLanczosExpecTwoBodyGFinish;
const char*  cCGExpecTwoBodyGFinish;
const char*  cTPQExpecTwoBodyGStart;
const char*  cTPQExpecTwoBodyGFinish;

//expec_energy.c
const char* cLogExpecEnergyStart;
const char* cLogLanczosExpecEnergyEnd;
const char* cLogCGExpecEnergyEnd;
const char* cLogEnergy;
const char* cLanczosExpecEnergyEnd;
const char* cCGExpecEnergyEnd;

#endif /* HPHI_LOGMESSAGE_H */
