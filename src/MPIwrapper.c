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
 * Main routine for the standard mode
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void InitializeMPI(int argc, char *argv[]){
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0){
    stdoutMPI = stdout;
  }
  else{
    stdoutMPI = fopen("/dev/null", "w");
  }
#else
  nproc = 1;
  myrank = 0;
  stdoutMPI = stdout;
#endif
}

void AbortMPI(int errorcode /**< [in]*/)
{
  int ierr;
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD,errorcode,ierr);
#else
  exit(errorcode);
#endif
}

void FinalizeMPI(){
#ifdef MPI
  MPI_Finalize();
#endif
}

void fopenMPI(char* FileName, FILE* fp, char* mode){
  if (myrank == 0){
    if ((fp = fopen(FileName, mode)) == NULL){
      fprintf(stdout,"\n  ERROR !  Cannot open %s !\n\n", FileName);
      AbortMPI(-1);
    }
  }
  else{
    fp = 0;
  }
}

void 

void fgetsMPI(FILE* fp, char* InputString)
{
  fgets(InputString, 256, fp);
#ifdef MPI
  MPI_Bcast(InputString, 256, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
}