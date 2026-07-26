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
#pragma once
#include "Common.h"

int expec_energy_flct(struct BindStruct *X);

int expec_energy_flct_Hubbard(struct BindStruct *X);

int expec_energy_flct_HubbardGC(struct BindStruct *X);

int expec_energy_flct_HalfSpinGC(struct BindStruct *X);

int expec_energy_flct_GeneralSpinGC(struct BindStruct *X);

int expec_energy_flct_HalfSpin(struct BindStruct *X);

int expec_energy_flct_GeneralSpin(struct BindStruct *X);

int expec_energy_flct_SpinlessFermionGC(struct BindStruct *X);

/* Per-basis-state ("per-k") fluctuation coefficient helpers (phase 3c Task 3).
 * Each returns the RAW quantity the corresponding evaluator loop computes for
 * basis index k (the 1-based loop index j == CSR row + 1): D(k) the doublon
 * count, N(k) the particle number, S(k) the raw 2*Sz bit sum. The caller
 * applies the same scalings it applies today. Shared verbatim by the Mode-1
 * evaluators and the CSR trace collector so no algebra is reimplemented. */
void EnergyFlctCoeff_Hubbard(struct BindStruct *X, long int k,
                             double *D, double *N, double *S);
void EnergyFlctCoeff_HubbardGC(struct BindStruct *X, long int k,
                               double *D, double *N, double *S);
void EnergyFlctCoeff_HalfSpinGC(struct BindStruct *X, long int k, double *S);
void EnergyFlctCoeff_GeneralSpinGC(struct BindStruct *X, long int k, double *S);
