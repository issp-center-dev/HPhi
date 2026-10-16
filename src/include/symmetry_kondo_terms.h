#ifndef HPHI_SYMMETRY_KONDO_TERMS_H
#define HPHI_SYMMETRY_KONDO_TERMS_H

#include <complex.h>
struct DefineList;
struct SymmetryTerm;

/* Independent local basis {I,E01,E10,E11}; identity sites are omitted.
 * The creation and annihilation strings are separately ascending. */
struct SymmetryKondoMonomial {
  unsigned int nlocal, ncreate, nannihilate;
  int local_site[2], local_out[2], local_in[2];
  int create[2], annihilate[2]; /* global 2*site+spin */
  double complex value;
};
typedef int (*SymmetryKondoMonomialCallback)(
    const struct SymmetryKondoMonomial *, void *);

/* NULL permutation means identity. Return 0 on success, -1 on invalid input,
 * allocation failure, or any nonzero callback result. */
int CanonicalizeSymmetryKondoTerm(const struct DefineList *def,
    const struct SymmetryTerm *term, const int *permutation,
    SymmetryKondoMonomialCallback callback, void *context);
int ValidateSymmetryKondoTerms(const struct DefineList *def);

#endif
