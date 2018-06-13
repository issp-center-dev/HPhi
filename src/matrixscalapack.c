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
#include "matrixscalapack.h"

#ifdef _SCALAPACK
int use_scalapack = 0;

/**
 * @file matrixscalapack.c
 * @version 3.1
 * 
 * @brief File for diagonalization using scalapack with functions for getting indices of global and local array
 */

long int GetBlockSize(long int Msize, long int nproc) {
  long int block_size = 16;
  if(Msize*Msize/nproc > block_size*block_size)
    return block_size;
  return 1;
}

/**
 * @brief A function to get processor array index from global array index
 * @param i global array index
 * @param np processor array dimension
 * @param nb block size
 * @return processor array index
 * @version 3.1 
 */
long int GetPArrayIndex(long int i, long int np, long int nb) {
  return (i/nb)%np;
}

/**
 * @brief A function to get local array index from global array index
 * @param i global array index
 * @param np processor array dimension
 * @param nb block size
 * @return local array index
 * @version 3.1 
 */
long int GetLocalIndex(long int i, long int np, long int nb) {
  return (i/(np*nb))*nb + i%nb;
}

/**
 * @brief A function to get global array index from local array index and processor array index
 * @param il local array index
 * @param p processor array index
 * @param np processor array dimension
 * @param nb block size
 * @return global array index
 * @version 3.1 
 */
long int GetGlobalIndex(long int il, long int p, long int np, long int nb){
  return ((il/nb)*np+p)*nb + il%nb;
}

/**
 * @brief A function to get rank of processor from indices of global matrix
 * @param i index of global matrix
 * @param j index of global matrix
 * @param nprow processor array dimension for row
 * @param npcol processor array dimension for column
 * @param nb block size
 * @return rank of processor
 * @version 3.1 
 */
long int MatToRank(long int i, long int j, long int nprow, long int npcol, long int nb){
  long int iproc, jproc;
  iproc = GetPArrayIndex(i, nprow, nb);
  jproc = GetPArrayIndex(j, npcol, nb);
  return iproc+jproc*nprow;
}

/**
 * @brief A function to get column index of global matrix from given rank of processor and column index of local matrix
 * @param lj column index of local matrix
 * @param rank rank of processor
 * @param npcol processor array dimension for column
 * @param nb block size
 * @return column index of global matrix
 * @version 3.1 
 */
long int GetMatRawInRank(long int lj, long int rank, long int npcol, long int nb){
  long int pcol = rank/npcol;
  return GetGlobalIndex(lj, pcol, npcol, nb);
}

/**
 * @brief A function to get indices of local matrix from indices (need to free memory after this function used)
 * @param i index of local matrix
 * @param j index of local matrix
 * @param nprow processor array dimension for row
 * @param npcol processor array dimension for column
 * @param nb block size
 * @return indices array of local matrix
 * @version 3.1 
 */
long int *GetMatElementInRank(long int i, long int j, long int nprow, long int npcol, long int nb){
  long int *ij;
  ij = malloc(2*sizeof(int));
  ij[0] = GetLocalIndex(i, nprow, nb);
  ij[1] = GetLocalIndex(j, npcol, nb);
  return ij;
}

/**
 * @brief A function to divide matrix
 * @param m index of column of matrix
 * @param n index of row of matrix
 * @param Aorgmn value of matrix in (i, j)
 * @param[in, out] A distribution matrix
 * @param desca descriptor of distribution matrix A
 * @version 3.1 
 */
void DivMat(long int m, long int n, double complex Aorgmn, double complex *A, int *desca){
  long int mp = m+1, np = n+1;
  pzelset_(A, &mp, &np, desca, &Aorgmn);
}

/**
 * @brief A function to get eigenvector from distributed matrix
 * @param i index of eigenvector
 * @param m size of eigenvector
 * @param Z distribution matrix of eigenvector
 * @param descZ descriptor for Z
 * @param[in, out] vec eigenvector
 * @version 3.1 
 */
void GetEigenVector(long int i, long int m, double complex *Z, int *descZ, double complex *vec) {
  double complex alpha;
  long int j, ip, jp;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ip = i+1;
  for(j=0; j<m; j++){
    jp = j+1;
    pzelget_("A", " ", &alpha, Z, &jp, &ip, descZ);
    if(rank==0) {
      vec[j] = alpha;
    }
  }
}

/**
 * @brief A function for diagonalization using scalapack
 * @param xNsize size of matrix
 * @param A input matrix
 * @param[in, out] r eigenvalue
 * @param[in, out] Z distribution matrix of eigenvector
 * @param[in, out] descZ descriptor for Z
 * @return this returns 0 when it finished normally
 * @version 3.1 
 */
int diag_scalapack_cmp(long int xNsize, double complex **A, 
                       double complex *r, double complex *Z, int *descZ) {
  const int i_one=1, i_zero=0;
  const long int i_negone=-1;
  const double zero=0.0, one=1.0;
  long int m, n, mb, nb;
  int nprow, npcol;
  int myrow, mycol, info, lld;
  long int mp, nq;
  int ictxt;
  complex double *A_distr, *work, *rwork;
  double *W;
  int descA_distr[9];
  int rank, size, iam, nprocs;
  long int lwork, lrwork;
  int dims[2]={0,0};
  long int i, j, ip, jp;
  m=n=xNsize;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size,2,dims);
  nprow=dims[0]; npcol=dims[1];
 
  blacs_pinfo_(&iam, &nprocs); 
  blacs_get_((int *)&i_negone, &i_zero, &ictxt);
  blacs_gridinit_(&ictxt, "R", &nprow, &npcol);
  blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);
 
  mb = GetBlockSize(m, size);
  nb = GetBlockSize(n, size);

  mp = numroc_(&m, &mb, &myrow, &i_zero, &nprow);
  nq = numroc_(&n, &nb, &mycol, &i_zero, &npcol);
  W = malloc(n*sizeof(double));
  A_distr = malloc(mp*nq*sizeof(complex double));

  lld = (mp>0) ? mp : 1;
  descinit_(descA_distr, &m, &n, &mb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info);
  descinit_(descZ, &m, &n, &mb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      DivMat(i, j, A[i][j], A_distr, descA_distr);
    }
  }

  double complex wkopt, rwkopt;
  pzheev_("V", "U", &n, A_distr, &i_one, &i_one, descA_distr, W, Z, &i_one, &i_one, descZ, &wkopt, &i_negone, &rwkopt, &i_negone, &info);

  lwork = (long int)wkopt;
  lrwork = (long int)rwkopt;
  work = malloc(lwork*sizeof(complex double));
  rwork = malloc(lrwork*sizeof(complex double));


  pzheev_("V", "U", &n, A_distr, &i_one, &i_one, descA_distr, W, Z, &i_one, &i_one, descZ, work, &lwork, rwork, &lrwork, &info);

  if(rank == 0){
    for(i=0; i<n; i++){
      r[i] = W[i];
    }
  }

  free(A_distr);
  free(work);
  free(rwork);
  free(W);

  use_scalapack = 1;

  return 0;
}

#endif
