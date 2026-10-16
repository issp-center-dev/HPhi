#ifndef HPHI_SYMMETRY_CORRELATION_H
#define HPHI_SYMMETRY_CORRELATION_H

#include <stddef.h>
#include <stdint.h>
#include "Common.h"

struct BindStruct;
struct DefineList;

#ifndef HPHI_SYMMETRY_CORRELATION_BLOCK_TRANSITIONS
#define HPHI_SYMMETRY_CORRELATION_BLOCK_TRANSITIONS UINT64_C(8388608)
#endif

struct SymmetryCorrelationOperator {
  unsigned int factors;
  const int *index;
};

struct SymmetryCorrelationStats {
  unsigned long long waves;
  unsigned long long transitions;
  unsigned long long unique_keys;
  unsigned long long peak_bytes;
  unsigned int threads;
};

int SymmetryCorrelationOrbitCount(const struct DefineList *def,
                                  const struct SymmetryCorrelationOperator *ops,
                                  size_t count,
                                  size_t *orbit_count,
                                  size_t *member_total);
int SymmetryCorrelationExpectationWithStats(
    const struct BindStruct *X, const double complex *vec,
    const struct SymmetryCorrelationOperator *ops, size_t count,
    double complex *values, struct SymmetryCorrelationStats *stats);
int SymmetryCorrelationExpectation(
    const struct BindStruct *X, const double complex *vec,
    const struct SymmetryCorrelationOperator *ops, size_t count,
    double complex *values);

#endif
