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

#include "setmemory.h"

///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \param A [in,out] Array to allocate
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
unsigned int *ui_1d_allocate(const long unsigned int N){
    unsigned int *A;
    A     = (unsigned int*)calloc((N),sizeof(unsigned int));
    return A;
}

///
/// \brief Function to free 1d array (int)
/// \param A Pointer of 1d array A
void free_ui_1d_allocate(unsigned int *A){
    free(A);
}

///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \param A [in,out] Array to allocate
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
long unsigned int *lui_1d_allocate(const long unsigned int N){
    long unsigned int *A;
    A     = (long unsigned int*)calloc((N),sizeof(long unsigned int));
    return A;
}

///
/// \brief Function to free 1d array (int)
/// \param A Pointer of 1d array A
void free_lui_1d_allocate(long unsigned int *A){
    free(A);
}

///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \param A [in,out] Array to allocate
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
long int *li_1d_allocate(const long unsigned int N){
    long int *A;
    A     = (long  int*)calloc((N),sizeof(long int));
    return A;
}

///
/// \brief Function to free 1d array (int)
/// \param A Pointer of 1d array A
void free_li_1d_allocate(long int *A){
    free(A);
}

///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
long int **li_2d_allocate(const long unsigned int N, const long unsigned int M) {
    long int **A;
    long unsigned int int_i;
    A = (long int **) calloc((N) , sizeof(long int *));
    A[0] = (long int *) calloc((M * N) ,sizeof(long int));
    for (int_i = 0; int_i < N; int_i++) {
        A[int_i] = A[0] + int_i * M;
    }
    //memset(A[0], 0, sizeof(long int)*M*N);
    return A;
}
///
/// \brief Function to free 2d array (int)
/// \param A Pointer of 2d array A
void free_li_2d_allocate(long int **A){
    free(A[0]);
    free(A);
}


///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \param A [in,out] Array to allocate
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
int *i_1d_allocate(const long unsigned int N){
    int *A;
    A     = (int*)calloc((N),sizeof(int));
    //memset(A, 0, sizeof(int)*N);
    return A;
}

///
/// \brief Function to free 1d array (int)
/// \param A Pointer of 1d array A
void free_i_1d_allocate(int *A){
    free(A);
}


///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
int **i_2d_allocate(const long unsigned int N, const long unsigned int M) {
    int **A;
    long unsigned int int_i;
    A = (int **) calloc((N) , sizeof(int *));
    A[0] = (int *) calloc((M * N) , sizeof(int));
    for (int_i = 0; int_i < N; int_i++) {
        A[int_i] = A[0] + int_i * M;
    }
    //memset(A[0], 0, sizeof(int)*M*N);
    return A;
}
///
/// \brief Function to free 2d array (int)
/// \param A Pointer of 2d array A
void free_i_2d_allocate(int **A){
    free(A[0]);
    free(A);
}

/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
int***i_3d_allocate(const long unsigned int N, const long unsigned int M, const long unsigned int L){
    long unsigned int int_i, int_j;
    int*** A;
    A     = (int***)calloc((N),sizeof(int**));
    A[0]  = (int**)calloc((M*N),sizeof(int*));
    A[0][0] = (int*)calloc((L*M*N),sizeof(int));
    for(int_i=0;int_i<N; int_i++) {
        A[int_i] = A[0] + int_i*M;
        for(int_j = 0; int_j<M; int_j++){
            A[int_i][int_j]= A[0][0] + int_i*M*L + int_j*L;
        }
    }
    return A;
}

///
/// \brief Function to free 3d array (int)
/// \param A A pointer of 3d array A
void free_i_3d_allocate(int ***A){
    free(A[0][0]);
    free(A[0]);
    free(A);
}

///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \param A [in,out] Array to allocate
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
double *d_1d_allocate(const long unsigned int N){
    double *A;
    A     = (double*)calloc((N),sizeof(double));
    return A;
}
///
/// \brief Function to free 1d array (double)
/// \param A Pointer of 1d array A
void free_d_1d_allocate(double *A){
    free(A);
}


///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
double **d_2d_allocate(const long unsigned int N, const long unsigned int M){
    long unsigned int int_i;
    double **A;
    A     = (double**)calloc((N),sizeof(double*));
    A[0]  = (double*)calloc((M*N),sizeof(double));
    for(int_i=0;int_i<N;int_i++){
        A[int_i] = A[0] + int_i*M;
    }
    return A;
}

///
/// \brief Function to free 2d array (double)
/// \param A Pointer of 2d array A
void free_d_2d_allocate(double **A){
    free(A[0]);
    free(A);
}

///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \param A [in,out] Array to allocate
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
double complex *cd_1d_allocate(const long unsigned int N){
    double complex*A;
    A     = (double complex*)calloc((N),sizeof(double complex));
    return A;
}
///
/// \brief Function to free 1d array (double complex)
/// \param A Pointer of 1d array A
void free_cd_1d_allocate(double complex *A){
    free(A);
}


///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
complex double **cd_2d_allocate(const long unsigned int N, const long unsigned int M){
    long unsigned int int_i;
    complex double **A;
    A     = (complex double**)calloc((N),sizeof(complex double));
    A[0]  = (complex double*)calloc((M*N),sizeof(complex double));
    for(int_i=0;int_i<N;int_i++){
        A[int_i] = A[0]+int_i*M;
    }
    return A;
}

///
/// \brief Function to free 2d array (complex double)
/// \param A Pointer of 2d array A
void free_cd_2d_allocate(double complex**A){
    free(A[0]);
    free(A);
}

/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \return A Pointer to array A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
double complex***cd_3d_allocate(const long unsigned int N, const long unsigned int M, const long unsigned int L){
    long unsigned int int_i, int_j;
    double complex***A;
    A     = (double complex***)calloc((N),sizeof(double complex**));
    A[0]  = (double complex**)calloc((M*N),sizeof(double complex*));
    A[0][0] = (double complex*)calloc((L*M*N),sizeof(double complex));
    for(int_i=0;int_i<N; int_i++) {
        A[int_i] = A[0] + int_i*M;
        for(int_j = 0; int_j<M; int_j++){
            A[int_i][int_j]= A[0][0] + int_i*M*L + int_j*L;
        }
    }
    return A;
}

///
/// \brief Function to free 3d array (complex double)
/// \param A A pointer of 3d array A
void free_cd_3d_allocate(double complex***A){
    free(A[0][0]);
    free(A[0]);
    free(A);
}
