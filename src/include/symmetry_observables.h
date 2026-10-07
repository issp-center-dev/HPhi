#ifndef HPHI_SYMMETRY_OBSERVABLES_H
#define HPHI_SYMMETRY_OBSERVABLES_H

#include "Common.h"

struct BindStruct;

int EvaluateSymmetrySpinGCMoments(struct BindStruct *X,
                                  const double complex *vec);

#endif
