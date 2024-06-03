/*
HPhi  -  Quantum Lattice Model Simulator
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
/**@file
@brief MPI wrapper for init, finalize, bcast, etc.
*/
#ifdef MPI
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wrapperMPI.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <complex.h>
#include "splash.h"
#include "global.h"

/**
@brief MPI initialization wrapper
Process ID (::myrank), Number of processes (::nproc), 
Number of threads (::nthreads), and pointer to the standard output
(::stdoutMPI) are specified here.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void InitializeMPI(int argc, char *argv[]){
#ifdef MPI
  int ierr;
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(ierr != 0) exitMPI(ierr);
#else
  nproc = 1;
  myrank = 0;
#endif
  if (myrank == 0) stdoutMPI = stdout;
  else stdoutMPI = fopen("/dev/null", "w");
  splash();

#pragma omp parallel default(none) shared(nthreads)
#pragma omp master
#ifdef _OPENMP
  nthreads = omp_get_num_threads();
#else
  nthreads=1;
#endif
  fprintf(stdoutMPI, "\n\n#####  Parallelization Info.  #####\n\n");
  fprintf(stdoutMPI, "  OpenMP threads : %d\n", nthreads);
  fprintf(stdoutMPI, "  MPI PEs : %d \n\n", nproc);
}/*void InitializeMPI(int argc, char *argv[])*/
/**
@brief MPI Finitialization wrapper
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void FinalizeMPI(){
#ifdef MPI
  int ierr;
  ierr = MPI_Finalize();
  if (ierr != 0) fprintf(stderr, "\n  MPI_Finalize() = %d\n\n", ierr);
#endif
  if (myrank != 0) fclose(stdoutMPI);
}
/**
@brief MPI Abortation wrapper
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void exitMPI(
  int errorcode//!<[in] Error-code to be returned as that of this program
)
{
  fflush(stdout);
#ifdef MPI
  fprintf(stdout,"\n\n #######  [HPhi] You DO NOT have to WORRY about the following MPI-ERROR MESSAGE.  #######\n\n");
  int ierr;
  ierr = MPI_Abort(MPI_COMM_WORLD, errorcode);
  ierr = MPI_Finalize();
  if (ierr != 0) fprintf(stderr, "\n  MPI_Finalize() = %d\n\n", ierr);
#endif
  exit(errorcode);
}/*void exitMPI*/
/**
@brief MPI file I/O (open) wrapper.
Only the root node (::myrank = 0) should be open/read/write (small) parameter files.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
FILE* fopenMPI(
  const char* FileName,//!<[in] Input/output file
  const char* mode//!<[in] "w", "r", etc.
){
  FILE* fp;

  if (myrank == 0) fp = fopen(FileName, mode);
  else fp = fopen("/dev/null", "w");

  return fp;
}/*FILE* fopenMPI*/
/**
@brief MPI file I/O (get a line, fgets) wrapper.
Only the root node (::myrank = 0) reads and broadcast string.
@return The same as that of fgets
@author Mitsuaki Kawamura (The University of Tokyo)
*/
char* fgetsMPI(
  char* InputString,//!<[out] read line.
  int maxcount,//!<[in] Length of string
  FILE* fp//!<[in] file pointer
){
  int inull;
  char *ctmp;

  ctmp = InputString;
  inull = 0;
  if (myrank == 0) {
    ctmp = fgets(InputString, maxcount, fp);
    if (ctmp == NULL){
      inull = 1;
    }
    
    while(*InputString == '\n' || strncmp(InputString, "#", 1)==0){
      ctmp = fgets(InputString, maxcount, fp);
      if (ctmp == NULL){
        inull=1;
        break;
      }
    }
  }
#ifdef MPI
  MPI_Bcast(InputString, maxcount, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&inull, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  if (myrank != 0 && inull == 1) {
    ctmp = NULL;
  }

  return ctmp;
}/*char* fgetsMPI*/
/**
@brief MPI barrier wrapper.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void BarrierMPI(){
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}/*void BarrierMPI()*/
/**
@brief MPI wrapper function to obtain maximum unsigned
long integer across processes.
@return Maximum value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
unsigned long int MaxMPI_li(
  unsigned long int idim//!<[in] Value to be maximized
){
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &idim, 1,
    MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  if(ierr != 0) exitMPI(-1);
#endif
  return(idim);
}/*unsigned long int MaxMPI_li*/
/**
@brief MPI wrapper function to obtain maximum Double
across processes.
@return Maximum value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double MaxMPI_d(
  double dvalue//!<[in] Value to be maximized
){
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &dvalue, 1,
    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if(ierr != 0) exitMPI(-1);
#endif
  return(dvalue);
}/*double MaxMPI_d*/
/**
@brief MPI wrapper function to obtain sum of Double
complex across processes.
@return Sumed value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double complex SumMPI_dc(
  double complex norm//!<[in] Value to be summed
){
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &norm, 1,
    MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  if(ierr != 0) exitMPI(-1);
#endif
  return(norm);
}/*double complex SumMPI_dc*/
/**
@brief MPI wrapper function to obtain sum of Double
across processes.
@return Sumed value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double SumMPI_d(
  double norm//!<[in] Value to be summed
){
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &norm, 1,
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
  if(ierr != 0) exitMPI(-1);
#endif
  return(norm);
}/*double SumMPI_d*/
/**
@brief MPI wrapper function to obtain sum of Double array
across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void SumMPI_dv(
  int nnorm,
  double *norm//!<[in] Value to be summed
) {
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, norm, nnorm,
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
  if (ierr != 0) exitMPI(-1);
#endif
}/*void SumMPI_dv*/
/**
@brief MPI wrapper function to obtain sum of Double array
across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void SumMPI_cv(
  int nnorm,
  double complex *norm//!<[in] Value to be summed
) {
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, norm, nnorm,
    MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  if (ierr != 0) exitMPI(-1);
#endif
}/*void SumMPI_cv*/
/**
@brief MPI wrapper function to obtain sum of unsigned
long integer across processes.
@return Sumed value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
unsigned long int SumMPI_li(
  unsigned long int idim//!<[in] Value to be summed
){
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &idim, 1,
    MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  if(ierr != 0) exitMPI(-1);
#endif
  return(idim);
}/*unsigned long int SumMPI_li*/
/**
@brief MPI wrapper function to obtain sum of
integer across processes.
@return Sumed value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
int SumMPI_i(
  int idim//!<[in] Value to be summed
) {
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &idim, 1,
                       MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(ierr != 0) exitMPI(-1);
#endif
  return(idim);
}/*int SumMPI_i*/
/**
@brief MPI wrapper function to broadcast unsigned long
integer across processes.
@return Broadcasted value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
unsigned long int BcastMPI_li(
  int root,//!<[in] The source process of the broadcast
  unsigned long int idim//!<[in] Value to be broadcasted
) {
  unsigned long int idim0;
  idim0 = idim;
#ifdef MPI
    MPI_Bcast(&idim0, 1, MPI_UNSIGNED_LONG, root, MPI_COMM_WORLD);
#endif
  return(idim0);
}/*unsigned long int BcastMPI_li*/
/**
@brief MPI wrapper function to broadcast integer across processes.
@return Broadcasted value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
int BcastMPI_i(
  int root,//!<[in] The source process of the broadcast
  int nsub//!<[in] Value to be broadcasted
) {
  int nsub0;
  nsub0 = nsub;
#ifdef MPI
  MPI_Bcast(&nsub0, 1, MPI_INT, root, MPI_COMM_WORLD);
#endif
  return(nsub0);
}/*int BcastMPI_i*/
/**
@brief MPI wrapper function to broadcast double precision vector.
       And store it in place.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void BcastMPI_dv(
  int root,//!<[in] The source process of the broadcast
  int nlen,//!<[in] Length of broadcasted vector
  double *vector//!<[inout] Broadcasted vector
) {
#ifdef MPI
  MPI_Bcast(vector, nlen, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD);
#endif
}/*void BcastMPI_dv*/
/**
@brief MPI wrapper function to broadcast double complex vector.
       And store it in place.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void BcastMPI_cv(
  int root,//!<[in] The source process of the broadcast
  int nlen,//!<[in] Length of broadcasted vector
  double complex* vector//!<[inout] Broadcasted vector
) {
#ifdef MPI
  MPI_Bcast(vector, nlen, MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD);
#endif
}/*void BcastMPI_cv*/
/**
@brief Compute norm of process-distributed vector
@f$|{\bf v}_1|^2@f$
@return Norm @f$|{\bf v}_1|^2@f$
*/
double NormMPI_dc(
  unsigned long int idim,//!<[in] Local dimension of vector
  double complex *_v1//!<[in] [idim] vector to be producted
){
  double complex cdnorm=0;
  double dnorm =0;
  unsigned long int i;
  //DEBUG
#pragma omp parallel for default(none) private(i) firstprivate(myrank) shared(_v1, idim) reduction(+: cdnorm)
  for(i=1;i<=idim;i++){
    cdnorm += conj(_v1[i])*_v1[i];
  }
#ifdef MPI
  cdnorm = SumMPI_dc(cdnorm);
#endif
  dnorm=creal(cdnorm);
  dnorm=sqrt(dnorm);

  return dnorm;
}/*double NormMPI_dc*/
/**
@brief Compute norm of process-distributed vector
@f$|{\bf v}_1|^2@f$
@return Norm @f$|{\bf v}_1|^2@f$
*/
void NormMPI_dv(
  unsigned long int ndim,//!<[in] Local dimension of vector
  int nstate,
  double complex **_v1,//!<[in] [idim] vector to be producted
  double *dnorm
) {
  unsigned long int idim;
  int istate;

  for (istate = 0; istate < nstate; istate++) dnorm[istate] = 0.0;
  for (idim = 1; idim <= ndim; idim++) {
    for (istate = 0; istate < nstate; istate++) {
      dnorm[istate] += conj(_v1[idim][istate])*_v1[idim][istate];
    }
  }
  SumMPI_dv(nstate, dnorm);
  for (istate = 0; istate < nstate; istate++) dnorm[istate] = sqrt(dnorm[istate]);
}/*double NormMPI_cv*/
/**
@brief Compute conjugate scaler product of process-distributed vector
@f${\bf v}_1^* \cdot {\bf v}_2@f$
@return Conjugate scaler product @f${\bf v}_1^* \cdot {\bf v}_2@f$
*/
double complex VecProdMPI(
  long unsigned int ndim,//!<[in] Local dimension of vector
  double complex *v1,//!<[in] [ndim] vector to be producted
  double complex *v2//!<[in] [ndim] vector to be producted
){
  long unsigned int idim;
  double complex prod;

  prod = 0.0;
#pragma omp parallel for default(none) shared(v1,v2,ndim) private(idim) reduction(+: prod)
  for (idim = 1; idim <= ndim; idim++) prod += conj(v1[idim]) * v2[idim];
  prod = SumMPI_dc(prod);

  return(prod);
}/*double complex VecProdMPI*/
/**
@brief Compute conjugate scaler product of process-distributed vector
@f${\bf v}_1^* \cdot {\bf v}_2@f$
*/
void MultiVecProdMPI(
  long unsigned int ndim,//!<[in] Local dimension of vector
  int nstate,
  double complex **v1,//!<[in] [ndim] vector to be producted
  double complex **v2,//!<[in] [ndim] vector to be producted
  double complex *prod
) {
  long unsigned int idim;
  int istate;

  for (istate = 0; istate < nstate; istate++) prod[istate] = 0.0;
  for (idim = 1; idim <= ndim; idim++) {
    for (istate = 0; istate < nstate; istate++) {
      prod[istate] += conj(v1[idim][istate])*v2[idim][istate];
    }
  }
  SumMPI_cv(nstate, prod);
}/*void MultiVecProdMPI*/
/**
@brief Wrapper of MPI_Sendrecv for double complex number.
When we pass a message longer than 2^31-1 
(max of int: 2147483647), we need to divide it.
*/
void SendRecv_cv(
  int origin,
  unsigned long int nMsgS,
  unsigned long int nMsgR,
  double complex *vecs,
  double complex *vecr
) {
#ifdef MPI
  int ierr, two31m1 = 2147483647, modMsg, nMsgS2, nMsgR2;
  unsigned long int nMsg, nnMsg, iMsg, sMsgR, sMsgS;
  MPI_Status statusMPI;

  if (nMsgS > nMsgR) nMsg = nMsgS;
  else nMsg = nMsgR;
  nnMsg = nMsg / two31m1;
  modMsg = nMsg % two31m1;
  if (modMsg != 0) nnMsg += 1;

  sMsgS = 0;
  sMsgR = 0;
  for (iMsg = 0; iMsg < nnMsg; iMsg++) {
    nMsgS2 = nMsgS / nnMsg;
    nMsgR2 = nMsgR / nnMsg;
    if (iMsg < nMsgS % nnMsg) nMsgS2 += 1;
    if (iMsg < nMsgR % nnMsg) nMsgR2 += 1;

    ierr = MPI_Sendrecv(&vecs[sMsgS], nMsgS2, MPI_DOUBLE_COMPLEX, origin, 0,
                        &vecr[sMsgR], nMsgR2, MPI_DOUBLE_COMPLEX, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    sMsgS += nMsgS2;
    sMsgR += nMsgR2;
  }
#endif
}/*void SendRecv_cv*/
/**
@brief Wrapper of MPI_Sendrecv for long unsigned integer number.
When we pass a message longer than 2^31-1
(max of int: 2147483647), we need to divide it.
*/
void SendRecv_iv(
  int origin,
  unsigned long int nMsgS,
  unsigned long int nMsgR,
  unsigned long int *vecs,
  unsigned long int *vecr
) {
#ifdef MPI
  int ierr, two31m1 = 2147483647, modMsg, nMsgS2, nMsgR2;
  unsigned long int nMsg, nnMsg, iMsg, sMsgR, sMsgS;
  MPI_Status statusMPI;

  if (nMsgS > nMsgR) nMsg = nMsgS;
  else nMsg = nMsgR;
  nnMsg = nMsg / two31m1;
  modMsg = nMsg % two31m1;
  if (modMsg != 0) nnMsg += 1;

  sMsgS = 0;
  sMsgR = 0;
  for (iMsg = 0; iMsg < nnMsg; iMsg++) {
    nMsgS2 = nMsgS / nnMsg;
    nMsgR2 = nMsgR / nnMsg;
    if (iMsg < nMsgS % nnMsg) nMsgS2 += 1;
    if (iMsg < nMsgR % nnMsg) nMsgR2 += 1;

    ierr = MPI_Sendrecv(&vecs[sMsgS], nMsgS2, MPI_UNSIGNED_LONG, origin, 0,
                        &vecr[sMsgR], nMsgR2, MPI_UNSIGNED_LONG, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    sMsgS += nMsgS2;
    sMsgR += nMsgR2;
  }
#endif
}/*void SendRecv_iv*/
/**
@brief Wrapper of MPI_Sendrecv for long unsigned integer number.
*/
unsigned long int SendRecv_i(
  int origin,
  unsigned long int isend
) {
#ifdef MPI
  int ierr;
  MPI_Status statusMPI;
  unsigned long int ircv;
  ierr = MPI_Sendrecv(&isend, 1, MPI_UNSIGNED_LONG, origin, 0,
                      &ircv,  1, MPI_UNSIGNED_LONG, origin, 0,
                      MPI_COMM_WORLD, &statusMPI);
  if (ierr != 0) exitMPI(ierr);
  return ircv;
#else
  return isend;
#endif
}/*void SendRecv_i*/
