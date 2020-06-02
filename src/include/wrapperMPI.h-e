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
@brief Global variables related to MPI and OpenMP
*/
#ifndef HPHI_WRAPPER_H
#define HPHI_WRAPPER_H
#include <complex.h>

void InitializeMPI(int argc, char *argv[]);
void FinalizeMPI();
void exitMPI(int errorcode);
FILE* fopenMPI(const char* FileName, const char* mode);
char* fgetsMPI(char* InputString, int maxcount,FILE* fp);
void BarrierMPI();
unsigned long int MaxMPI_li(unsigned long int idim);
double MaxMPI_d(double dvalue);
double complex SumMPI_dc(double complex norm);
double SumMPI_d(double norm);
unsigned long int SumMPI_li(unsigned long int idim);
int SumMPI_i(int idim);
unsigned long int BcastMPI_li(int root, unsigned long int idim);
double NormMPI_dc(unsigned long int idim, double complex *_v1);
double complex VecProdMPI(long unsigned int ndim, double complex *v1, double complex *v2);

#endif
