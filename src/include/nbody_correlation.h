/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

#pragma once
#include "Common.h"

int ParseNBodyGLine(
  const char *line,
  unsigned int *N,
  int **factors
);

int ValidateNBodyGScope(const struct DefineList *D);
int NormalizeNBodyGTerms(struct DefineList *D);
int CheckNBodyGSpinConservation(const struct DefineList *D);
int CheckNBodyGHubbardConservation(const struct DefineList *D);
int expec_nbodyg(struct BindStruct *X, double complex *vec);
