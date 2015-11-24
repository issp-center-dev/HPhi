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
#ifndef HPHI_WRAPPER_H
#define HPHI_WRAPPER_H

int nproc, myrank;
FILE *stdoutMPI;

void InitializeMPI(int argc, char *argv[]);
void FinalizeMPI();
void exitMPI(int errorcode);
FILE* fopenMPI(const char* FileName, const char* mode);
char* fgetsMPI(char* InputString, int maxcount,FILE* fp);
void BarrierMPI();
unsigned long int RedduceMaxMPI_li(unsigned long int idim);

#endif
