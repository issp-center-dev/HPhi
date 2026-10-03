#ifndef HPHI_SYMMETRY_TERMS_H
#define HPHI_SYMMETRY_TERMS_H

#include <complex.h>
struct DefineList;

/* Ordered product of one or two c^dagger c factors (local matrix units for
 * Spin). Coefficients already include HPhi's Transfer/Hund sign convention. */
struct SymmetryTerm {
  unsigned int factors;
  int index[8];
  double complex value;
};
typedef int (*SymmetryTermCallback)(const struct SymmetryTerm *, void *);

/* kind: -1 all, 0 diagonal, 1 off-diagonal. Uses the reader's split InterAll
 * arrays and EDChemi, never the original InterAll as well. */
int EnumerateSymmetryTerms(const struct DefineList *def, int kind,
                          SymmetryTermCallback callback, void *context);
/* Return 1 for a nonzero matrix element, 0 for zero, -1 for invalid input. */
int ApplySymmetryTerm(const struct DefineList *def,
                      const struct SymmetryTerm *term, unsigned long state,
                      unsigned long *out, double complex *value);
int ValidateSymmetryTerms(const struct DefineList *def);
int SymmetryUsesExtendedTerms(const struct DefineList *def);

#endif
