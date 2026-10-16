#ifndef HPHI_SYMMETRY_TE_H
#define HPHI_SYMMETRY_TE_H
#include "Common.h"
#include <stdint.h>

/* Owns only replacement Transfer/InterAll arrays; the base definition is borrowed. */
struct SymmetryTEHamiltonian {
  struct DefineList base;
  int **transfer, **interall;
  int *transfer_data, *interall_data;
  double complex *transfer_values, *interall_values;
  uint64_t *digests;
};
int SymmetryTEIsDynamic(const struct DefineList *def);
double SymmetryTETime(const struct DefineList *def, unsigned int step);
int ValidateSymmetryTESchedule(struct BindStruct *X, struct SymmetryTEHamiltonian *view);
int SelectSymmetryTEHamiltonian(struct BindStruct *X, struct SymmetryTEHamiltonian *view,
                                unsigned int step, double time);
int RebuildSymmetryTEPlan(struct BindStruct *X);
void FreeSymmetryTEHamiltonian(struct SymmetryTEHamiltonian *view);
#endif
