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
/* Apply an ordered product of factors from right to left. The index array has
 * 4*factors entries (site_out, spin_out, site_in, spin_in per factor). No
 * fixed factor limit is imposed; only overflow of 4*factors is rejected.
 * Return 1 for a nonzero element, 0 for zero, and -1 for invalid input. */
int ApplySymmetryFactors(const struct DefineList *def, unsigned int factors,
                         const int *index, unsigned long state,
                         unsigned long *out, double *sign);
/* The fixed-size SymmetryTerm wrapper retains its one/two-factor contract and
 * multiplies the returned sign by term->value. */
int ApplySymmetryTerm(const struct DefineList *def,
                      const struct SymmetryTerm *term, unsigned long state,
                      unsigned long *out, double complex *value);
int ValidateSymmetryTerms(const struct DefineList *def);
int SymmetryUsesExtendedTerms(const struct DefineList *def);

#endif
