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
  int ierr;

#ifdef MPI
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
  int ierr;
#ifdef MPI
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
  int ierr;
  fflush(stdout);
#ifdef MPI
  fprintf(stdout,"\n\n #######  [HPhi] You DO NOT have to WORRY about the following MPI-ERROR MESSAGE.  #######\n\n");
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
