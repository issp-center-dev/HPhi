#ifndef HPHI_TEST_SYMMETRY_PROBE_H
#define HPHI_TEST_SYMMETRY_PROBE_H
#include <complex.h>
struct BindStruct;
int SymmetryProbeBeforeSolver(struct BindStruct *X);
int SymmetryProbeFinalVector(const struct BindStruct *X, const double complex *vector);
int SymmetryProbeWriteVector(const struct BindStruct *X, const double complex *vector,
                             const char *stage, int sample, int step, double prenorm);
#endif
