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
#pragma once
#define PI 3.14159265358979 

/*!< CalcType */
#define NUM_CALCTYPE 3
#define Lanczos 0 /*!< CalcType is Exact Diagonalization method.*/
#define TPQCalc 1 /*!< CalcType is TPQ calculation.*/
#define FullDiag 2 /*!< CalcType is Full Diagonalization method.*/

/*!< CalcModel */
#define NUM_CALCMODEL 6 /*!< Number of model types defined by CalcModel.*/
#define Hubbard 0 /*!< CalcModel is Hubbard model.*/
#define Spin 1 /*!< CalcModel is Spin system.*/
#define Kondo 2 /*!< CalcModel is Kondo model.*/
#define HubbardGC 3 /*!< CalcModel is GrandCanonical Hubbard model.*/
#define SpinGC 4 /*!< CalcModel is GrandCanonical Spin system.*/
#define KondoGC 5 /*!< CalcModel is GrandCanonical Kondo model.*/

/*!< OutputMode */
#define NUM_OUTPUTMODE 2 /*!< Number of output mode.*/
#define RAWMODE  0 /*!< calc one body green function and two body green functions.*/
#define CORRMODE 1 /*!< calc one body green function and two body green functions and correlatinos for charge and spin.*/
