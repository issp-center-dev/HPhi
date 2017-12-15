#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int OutputRealEigenValue(int xNsize, double *ene, char *filename);
int OutputCmpEigenValue(int xNsize, complex double *ene, char *filename);
int OutputRealEigenVec(int xNsize, const int nene, double **vec, const int nproc, char *filename);
int OutputCmpEigenVec(int xNsize, const int nene, complex double **vec, const int nproc, char *filename);
int InputRealEigenValue(int xNsize, double *ene, char *filename);
int InputCmpEigenValue(int xNsize, complex double *ene, char *filename);
int InputRealEigenVec(int xNsize, const int nene, double **vec, const int nproc, char *filename);
int InputCmpEigenVec(int xNsize, const int nene, complex double **vec, const int nproc, char *filename);
