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

#ifndef HPHI_DEFCOMMON_H
#define HPHI_DEFCOMMON_H

/*!< CalcType */
#define NUM_CALCTYPE 5
#define Lanczos 0 /*!< CalcType is Exact Diagonalization method.*/
#define TPQCalc 1 /*!< CalcType is TPQ calculation.*/
#define FullDiag 2 /*!< CalcType is Full Diagonalization method.*/
#define CG 3 /*!< CalcType is CG method */
#define TimeEvolution 4 /*!< CalcType is Time Evolution method*/

/*!< CalcModel */
#define NUM_CALCMODEL 6 /*!< Number of model types defined by CalcModel in calcmodfile. Note: HubbardNConserved is not explicitly defined in calcmod file and thus not counted. SpinlessFermion and SpinlessFermionGC are not yet supported*/
#define Hubbard 0 /*!< CalcModel is Hubbard model.*/
#define Spin 1 /*!< CalcModel is Spin system.*/
#define Kondo 2 /*!< CalcModel is Kondo model.*/
#define HubbardGC 3 /*!< CalcModel is GrandCanonical Hubbard model.*/
#define SpinGC 4 /*!< CalcModel is GrandCanonical Spin system.*/
#define KondoGC 5 /*!< CalcModel is GrandCanonical Kondo model.*/
#define HubbardNConserved 6 /*!< CalcModel is Hubbard model under particle number conserved. This symmetry is automatically introduced by not defining 2Sz in a modpara file.*/
#define SpinlessFermion 7 /*!< CalcModel is GrandCanonical Spinless fermion model.*/
#define SpinlessFermionGC 8 /*!< CalcModel is GrandCanonical Spinless fermionGC model.*/

/*!< OutputMode */
#define NUM_OUTPUTMODE 2 /*!< Number of output mode.*/
#define RAWMODE  0 /*!< calc one body green function and two body green functions.*/
#define CORRMODE 1 /*!< calc one body green function and two body green functions and correlatinos for charge and spin.*/
#define NUM_OUTPUTHAM 2 /*!< Number of output Hamiltonian mode */

/*!< InputMode */
#define NUM_INPUTHAM 2 /*!< Number of input Hamiltonian mode */

/*!< CalcEigenVector */
#define NUM_SETINITAILVEC 2 /*!< Number of setting type of initial vectors.*/
#define NUM_CALCEIGENVEC 2 /*!< Number of calculating eigenvector mode.*/
#define CALCVEC_LANCZOSCG  0 /*!< Lanczos + CG method*/
#define CALCVEC_LANCZOS 1 /*!< Lanczos method*/
#define CALCVEC_NOT -1 /*!< eigenvector is not calculated*/

/*!< CalcSpectrum */
#define CALCSPEC_NOT 0
#define RECALC_NOT 1
#define RECALC_FROM_TMComponents 2
#define RECALC_OUTPUT_TMComponents_VEC 3
#define RECALC_FROM_TMComponents_VEC 4
#define RECALC_INOUT_TMComponents_VEC 5
#define CALCSPEC_SCRATCH 6

/*!< ReStartVector */
#define NUM_RESTART 4
#define RESTART_NOT 0
#define RESTART_OUT 1
#define RESTART_INOUT 2
#define RESTART_IN 3

#endif /* HPHI_DEFCOMMON_H */
