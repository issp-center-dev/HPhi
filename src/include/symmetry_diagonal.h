#ifndef HPHI_SYMMETRY_DIAGONAL_H
#define HPHI_SYMMETRY_DIAGONAL_H

#include "Common.h"

struct DefineList;

int EvaluateSymmetryStateDiagonal(
    const struct DefineList *def,
    unsigned long int state,
    double *diagonal);

#endif /* HPHI_SYMMETRY_DIAGONAL_H */
