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
int Lanczos_EigenValue(struct BindStruct *X);
int Lanczos_GetTridiagonalMatrixComponents(struct BindStruct *X, double *alpha, double *beta, double complex *_v1, unsigned long int *Lanczos_step);

int ReadInitialVector(struct BindStruct *X, double complex* tmp_v0, double complex *tmp_v1, unsigned long int *liLanczosStp_vec);

int OutputLanczosVector(struct BindStruct *X, double complex* tmp_v0, double complex *tmp_v1, unsigned long int liLanczosStp_vec);

void SetInitialVector(struct BindStruct *X, double complex* tmp_v0, double complex *tmp_v1);

int ReadTMComponents(
        struct BindStruct *X,
        double *_dnorm,
        unsigned long int *i_max,
        const int iFlg
);

int OutputTMComponents(
        struct BindStruct *X,
        double *_alpha,
        double *_beta,
        double _dnorm,
        int liLanczosStp
);

