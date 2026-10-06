#ifndef HPHI_TEST_SYMMETRY_SPINGC_PROBE_H
#define HPHI_TEST_SYMMETRY_SPINGC_PROBE_H
#include <complex.h>
struct BindStruct;
int SpinGCProbeBeforeSolver(struct BindStruct *X);
int SpinGCProbeFinalVector(const struct BindStruct *X,
                           const double complex *vec);
#endif
