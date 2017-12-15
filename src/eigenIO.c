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
  int i = 0;

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

