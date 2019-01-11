/*
HPhi-mVMC-StdFace - Common input generator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//
// Created by Kazuyoshi Yoshimi on 2019-01-09.
//

#ifndef MVMC_SETMEMORY_H
#define MVMC_SETMEMORY_H

#include <stdlib.h>
#include <string.h>
#include <complex.h>
///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \param A [in,out] Array to allocate
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
int *i_1d_allocate(const int N);

///
/// \brief Function to free 1d array (int)
/// \param A Pointer of 1d array A
void free_i_1d_allocate(int *A);


///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
int **i_2d_allocate(int N, int M);
///
/// \brief Function to free 2d array (int)
/// \param A Pointer of 2d array A
void free_i_2d_allocate(int **A);


///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \param A [in,out] Array to allocate
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
double *d_1d_allocate(const int N);

///
/// \brief Function to free 1d array (double)
/// \param A Pointer of 1d array A
void free_d_1d_allocate(double *A);


///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
double **d_2d_allocate(int N, int M);
///
/// \brief Function to free 2d array (double)
/// \param A Pointer of 2d array A
void free_d_2d_allocate(double **A);


///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
complex double *cd_1d_allocate(int N);

///
/// \brief Function to free 1d array (complex double)
/// \param A Pointer of 2d array A
void free_cd_1d_allocate(double complex*A);

///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
complex double **cd_2d_allocate(int N, int M);

///
/// \brief Function to free 2d array (complex double)
/// \param A Pointer of 2d array A
void free_cd_2d_allocate(double complex**A);

//
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
double complex***cd_3d_allocate(int N, int M, int L);
///
/// \brief Function to free 3d array (complex double)
/// \param A A pointer of 3d array A
void free_cd_3d_allocate(double complex***A);

#endif //MVMC_SETMEMORY_H
