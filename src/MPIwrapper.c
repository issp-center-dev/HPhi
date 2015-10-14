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
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "MPIwrapper.h"

/**
 *
 * MPI initialization wrapper
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void InitializeMPI(int argc, char *argv[]){
  int ierr;

#ifdef MPI
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(ierr != 0) AbortMPI(ierr);
#else
  nproc = 1;
  myrank = 0;
#endif
  if (myrank == 0) stdoutMPI = stdout;
  else stdoutMPI = fopen("/dev/null", "w");
}

/**
 *
 * MPI Abortation wrapper
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void AbortMPI(int errorcode /**< [in]*/)
{
  int ierr;
#ifdef MPI
  ierr = MPI_Abort(MPI_COMM_WORLD,errorcode);
  ierr = MPI_Finalize();
  if (ierr != 0) fprintf(stdoutMPI, "\n  MPI_Finalize() = %d\n\n", ierr);
#endif
  exit(errorcode);
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
  if (ierr != 0) fprintf(stdoutMPI, "\n  MPI_Finalize() = %d\n\n", ierr);
#endif
}

/**
 *
 * MPI file I/O (open) wrapper
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
FILE* fopenMPI(
  char* FileName /**< [in] Input/output file*/, 
  char* mode /**< [in] "w", "r", etc. */){
  FILE* fp;

  if (myrank == 0) fp = fopen(FileName, mode);
  else fp = fopen("/dev/null", "w");

  return fp;
}

/**
 *
 * MPI file I/O (close) wrapper
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void  fcloseMPI(
  FILE* fp /**< [in]*/){
  if (myrank == 0) fclose(fp);
}

/**
 *
 * MPI file I/O (get a line) wrapper
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void fgetsMPI(
  FILE* fp /**< [in] file pointer*/, 
  char* InputString /**< [out] read line. The length must be 256*/)
{
  fgets(InputString, 256, fp);
#ifdef MPI
  MPI_Bcast(InputString, 256, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
}

