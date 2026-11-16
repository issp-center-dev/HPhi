#ifndef HPHI_SYMMETRY_OBSERVABLES_H
#define HPHI_SYMMETRY_OBSERVABLES_H

#include "Common.h"

struct BindStruct;

int EvaluateSymmetrySpinGCMoments(struct BindStruct *X,
                                  const double complex *vec);
/* N and N^2 of a SpinlessFermionGC sector state; spin moments are zero. */
int EvaluateSymmetrySpinlessGCMoments(struct BindStruct *X,
                                      const double complex *vec);

#endif
