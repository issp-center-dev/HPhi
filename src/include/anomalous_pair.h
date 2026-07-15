/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

#pragma once
#include "Common.h"

int ParseAnomalousTermLine(
  const char *line,
  int term[5],
  double *re,
  double *im
);

int ParseAnomalousGLine(
  const char *line,
  int term[5]
);

int ValidateAnomalousTermScope(const struct DefineList *D);
int ValidateAnomalousGScope(const struct DefineList *D);
int CheckAnomalousTermHermitePairs(const struct DefineList *D);

int ApplyAnomalousPairHubbardGC(
  const struct DefineList *D,
  const int term[5],
  unsigned long int local_in,
  int rank_in,
  unsigned long int *local_out,
  int *rank_out,
  int *sign
);

int MultiplyAnomalousTermHubbardGC(
  struct BindStruct *X,
  double complex *tmp_v0,
  double complex *tmp_v1
);

int AddAnomalousTermToHamHubbardGC(struct BindStruct *X);
int expec_anomalousg(struct BindStruct *X, double complex *vec);
