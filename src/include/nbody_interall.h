/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

#pragma once
#include "Common.h"

int ParseNBodyInterAllLine(
  const char *line,
  unsigned int *N,
  int **factors,
  double *re,
  double *im
);

int ValidateNBodyInterAllScope(const struct DefineList *D);
int NormalizeNBodyInterAllTerms(struct DefineList *D);
int CheckNBodyInterAllSpinConservation(const struct DefineList *D);
int ClassifyNBodyInterAllTerms(struct DefineList *D);
int CheckNBodyInterAllHermitePairs(const struct DefineList *D);

int ApplyNBodyInterAllSpinGC(
  const struct BindStruct *X,
  unsigned int term_index,
  unsigned long int local_in,
  int rank_in,
  unsigned long int *local_out,
  int *rank_out,
  double complex *matrix_element
);

int SetDiagonalNBodyInterAllSpinGC(struct BindStruct *X);
int SetDiagonalNBodyInterAllHubbardGC(struct BindStruct *X);
int MultiplyNBodyInterAllSpinGC(
  struct BindStruct *X,
  double complex *tmp_v0,
  double complex *tmp_v1
);
int MultiplyNBodyInterAllHubbardGC(
  struct BindStruct *X,
  double complex *tmp_v0,
  double complex *tmp_v1
);
int AddNBodyInterAllToHamSpinGC(struct BindStruct *X);
int AddNBodyInterAllToHamHubbardGC(struct BindStruct *X);
