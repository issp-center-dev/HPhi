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
#include "Common.h"

int CalcSpectrum(
                 struct EDMainCalStruct *X
);

/** @brief A set of excitation operators (ket A or bra B) to apply to |phi>.
    Lets GetExcitedState build either A|phi> (ket / X->Def fields) or B|phi>
    (bra / X->Def *Bra fields) without mutating DefineList. */
typedef struct {
  unsigned int NSingle;       /**< number of single excitation operators */
  int **Single;               /**< [NSingle][3] = {site, spin, type} */
  double complex *ParaSingle; /**< [NSingle] coefficients */
  unsigned int NPair;         /**< number of pair excitation operators */
  int **Pair;                 /**< [NPair][5] = {site1, spin1, site2, spin2, type} */
  double complex *ParaPair;   /**< [NPair] coefficients */
} ExcitationOperatorSet;

int GetExcitedState(
                struct BindStruct *X,
                const ExcitationOperatorSet *op,
                double complex *tmp_v0,
                double complex *tmp_v1
);


int MakeExcitedList(
                struct BindStruct *X,
                  int *iFlgListModifed
                );

int ReSetList(struct BindStruct *X);

int SetOmega
(
 struct DefineList *X
 );
