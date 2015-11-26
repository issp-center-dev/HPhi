/*
HPhi  -  Quantum Lattice Model Simulator
Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima

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
#ifdef MPI
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include "wrapperMPI.h"
#include <omp.h>
#include <math.h>
#include <complex.h>

/**
 *
 * MPI initialization wrapper
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void InitializeMPI(int argc, char *argv[]){
  int ierr, nthreads;

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

#pragma omp parallel default(none) shared(nthreads)
#pragma omp master
  nthreads = omp_get_num_threads();

  fprintf(stdoutMPI, "\n\n#####  Parallelization Info.  #####\n\n");
  fprintf(stdoutMPI, "  OpenMP threads : %d\n", nthreads);
  fprintf(stdoutMPI, "  MPI PEs : %d \n\n", nproc);
}

/**
 *
 * MPI Finitialization wrapper
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
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
*
* MPI Abortation wrapper
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void exitMPI(int errorcode /**< [in]*/)
{
  int ierr;
#ifdef MPI
  fprintf(stderr,"\n\n #######  [HPhi] You DO NOT have to WORRY about the following MPI-ERROR MESSAGE.  #######\n\n");
  fprintf(stdout,"\n\n #######  [HPhi] You DO NOT have to WORRY about the following MPI-ERROR MESSAGE.  #######\n\n");
  ierr = MPI_Abort(MPI_COMM_WORLD, errorcode);
  ierr = MPI_Finalize();
  if (ierr != 0) fprintf(stderr, "\n  MPI_Finalize() = %d\n\n", ierr);
#endif
  exit(errorcode);
}

/**
 *
 * MPI file I/O (open) wrapper
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
FILE* fopenMPI(
  const char* FileName /**< [in] Input/output file*/, 
  const char* mode /**< [in] "w", "r", etc. */){
  FILE* fp;

  if (myrank == 0) fp = fopen(FileName, mode);
  else fp = fopen("/dev/null", "w");

  return fp;
}

/**
 *
 * MPI file I/O (get a line) wrapper
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
char* fgetsMPI(
  char* InputString, /**< [out] read line.*/
  int maxcount /**< [in] Length of string*/,
  FILE* fp /**< [in] file pointer*/)
{
  int inull;
  char *ctmp;

  ctmp = InputString;
  inull = 0;
  if (myrank == 0) {
    ctmp = fgets(InputString, maxcount, fp);
    if (ctmp == NULL){
      inull = 1;
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
}

void BarrierMPI(){
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

unsigned long int MaxMPI_li(unsigned long int idim)
{
  int ierr;
#ifdef MPI
  ierr = MPI_Allreduce(MPI_IN_PLACE, &idim, 1,
    MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
#endif
  return(idim);
}

double complex SumMPI_dc(double complex norm)
{
  int ierr;
#ifdef MPI
  ierr = MPI_Allreduce(MPI_IN_PLACE, &norm, 1,
    MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
  return(norm);
}

double SumMPI_d(double norm)
{
  int ierr;
#ifdef MPI
  ierr = MPI_Allreduce(MPI_IN_PLACE, &norm, 1,
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
#endif
  return(norm);
}

unsigned long int SumMPI_li(unsigned long int idim)
{
  int ierr;
#ifdef MPI
  ierr = MPI_Allreduce(MPI_IN_PLACE, &idim, 1,
    MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  return(idim);
}

int SumMPI_i(int idim)
{
  int ierr;
#ifdef MPI
  ierr = MPI_Allreduce(MPI_IN_PLACE, &idim, 1,
    MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  return(idim);
}

