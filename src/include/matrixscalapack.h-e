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
/*=================================================================================================*/

#pragma once
#ifdef _SCALAPACK
#include "mpi.h"
#include "global.h"
#include <complex.h>
#include <stdlib.h>

void blacs_pinfo_(int *mypnum, int *nprocs); 
void blacs_get_(const int *context, const int *request, int *value ); 
void blacs_gridinfo_(const int *ictxt, int *nprow, int *npcol, int *myrow, int *mycol);
void blacs_gridinit_(const int *ictxt, char *order, int *nprow, int *npcol );
void blacs_gridexit_(const int *ictxt);

void descinit_(int *desc, const long int *m, const long int *n, const long int *mb, 
               const long int *nb, const int *irsrc, const int *icsrc, 
               const int *ixtxt, const int *lld, int *info);
void pzelset_(double complex *A, const long int *ia, const long int *ja, 
              const int *descA, const double complex *alpha);
void pzelget_(char *scope, char *top, double complex *alpha, 
              const double complex *A, const long int *ia, const long int *ja, 
              const int *descA);
long int numroc_(const long int *n, const long int *nb, const int *iproc, 
            const int *isrcproc, const int *nprocs);
void pzheev_(char *jobz, char *uplo, const long int *n, double complex *a, 
             const int *ia, const int *ja, int *desca, double *w, 
             double complex *z, const int *iz, const int *jz, int *descz, 
             double complex *work, const long int *lwork, double complex *rwork, 
             const long int *lrwork, int *info );

long int GetBlockSize(long int Msize, long int nproc);
long int GetPArrayIndex(long int i, long int np, long int nb);
long int GetLocalIndex(long int i, long int np, long int nb);
long int GetGlobalIndex(long int il, long int p, long int np, long int nb);
long int MatToRank(long int i, long int j, long int nprow, long int npcol, long int nb);
long int GetMatRawInRank(long int lj, long int rank, long int npcol, long int nb);
long int *GetMatElementInRank(long int i, long int j, long int nprow, long int npcol, long int nb);
void DivMat(long int m, long int n, double complex Aorgmn, double complex *A, int *desca);
void GetEigenVector(long int i, long int m, double complex *Z, int *descZ, double complex *vec);
int diag_scalapack_cmp(long int xNsize, double complex **A,
                       double complex *r, double complex *Z, int *descZ);

extern int use_scalapack;
#endif

