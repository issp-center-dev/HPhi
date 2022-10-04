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
/*-------------------------------------------------------------*/
#include "eigenIO.h"

int comp(const void *p, const void *q) {
  if(*(double *)p > *(double *)q) return 1;
  if(*(double *)p < *(double *)q) return -1;
  return 0;
}

int OutputRealEigenValue(int xNsize, double *ene, char *filename) {
  FILE *fp = NULL;
  int i;
  double *buffene = (double *) malloc(xNsize * sizeof(double));

  for(i = 0; i < xNsize; i++)
    buffene[i] = ene[i];

  qsort(buffene, xNsize, sizeof(double), comp);

  fp = fopen(filename, "wb+");
  if(fp == NULL) {
    free(buffene);
    return -1;
  }

  fwrite(buffene, sizeof(double), xNsize, fp);

  fclose(fp);
  free(buffene);

  return 0;
}

int OutputCmpEigenValue(int xNsize, complex double *ene, char *filename) {
  FILE *fp = NULL;
  int i;
  double *buffene = (double *) malloc(xNsize * sizeof(double));

  for(i = 0; i < xNsize; i++)
    buffene[i] = creal(ene[i]);

  qsort(buffene, xNsize, sizeof(double), comp);

  fp = fopen(filename, "wb+");
  if(fp == NULL) {
    free(buffene);
    return -1;
  }

  fwrite(buffene, sizeof(double), xNsize, fp);

  fclose(fp);
  free(buffene);

  return 0;
}

int OutputRealEigenVec(int xNsize, const int nene, double **vec, const int nproc, char *filename) {
  FILE *fp = NULL;

  fp = fopen(filename, "wb+");
  if(fp == NULL) {
    return -1;
  }

  fwrite(vec[nene], sizeof(double), xNsize, fp);

  fclose(fp);

  return 0;
}

int OutputCmpEigenVec(int xNsize, const int nene, complex double **vec, const int nproc, char *filename) {
  FILE *fp = NULL;

  fp = fopen(filename, "wb+");
  if(fp == NULL) {
    return -1;
  }

  fwrite(vec[nene], sizeof(complex double), xNsize, fp);

  fclose(fp);

  return 0;
}

int InputRealEigenValue(int xNsize, double *ene, char *filename) {
  FILE *fp = NULL;

  fp = fopen(filename, "rb+");
  if(fp == NULL) {
    return -1;
  }

  fread(ene, sizeof(double), xNsize, fp);

  fclose(fp);

  return 0;
}

int InputCmpEigenValue(int xNsize, complex double *ene, char *filename) {
  FILE *fp = NULL;
  fp = fopen(filename, "rb+");
  if(fp == NULL) {
    return -1;
  }

  fread(ene, sizeof(complex double), xNsize, fp);

  fclose(fp);

  return 0;
}

int InputRealEigenVec(int xNsize, const int nene, double **vec, const int nproc, char *filename) {
  FILE *fp = NULL;

  fp = fopen(filename, "rb+");
  if(fp == NULL) {
    return -1;
  }

  fread(vec[nene], sizeof(double), xNsize, fp);

  fclose(fp);

  return 0;
}

int InputCmpEigenVec(int xNsize, const int nene, complex double **vec, const int nproc, char *filename) {
  FILE *fp = NULL;

  fp = fopen(filename, "rb+");
  if(fp == NULL) {
    return -1;
  }

  fread(vec[nene], sizeof(complex double), xNsize, fp);

  fclose(fp);

  return 0;
}

