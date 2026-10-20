#ifndef HPHI_SYMMETRY_TERMS_H
#define HPHI_SYMMETRY_TERMS_H

#include <complex.h>
struct DefineList;

/* BILINEARS: ordered product of one or two c^dagger c factors (local matrix
 * units for Spin). The pair kinds are c_i c_j and c^dagger_i c^dagger_j of
 * AnomalousTerm (types 0 and 1), with factors=1 and index (i, si, j, sj).
 * Coefficients already include HPhi's Transfer/Hund sign convention. */
enum SymmetryTermKind {
  SYMMETRY_TERM_BILINEARS = 0,
  SYMMETRY_TERM_ANNIHILATE_PAIR,
  SYMMETRY_TERM_CREATE_PAIR
};
struct SymmetryTerm {
  unsigned int factors;
  int index[8];
  double complex value;
  enum SymmetryTermKind kind; /* last, so positional initializers stay bilinear */
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
/* Apply c_i c_j (create=0) or c^dagger_i c^dagger_j (create=1) of a GC
 * fermion model from the right, so the j operator acts first. index is
 * (site_i, spin_i, site_j, spin_j) of two different orbitals. Return values
 * are those of ApplySymmetryFactors. */
int ApplySymmetryPair(const struct DefineList *def, int create, const int *index,
                      unsigned long state, unsigned long *out, double *sign);
/* The fixed-size SymmetryTerm wrapper retains its one/two-factor contract and
 * multiplies the returned sign by term->value. */
int ApplySymmetryTerm(const struct DefineList *def,
                      const struct SymmetryTerm *term, unsigned long state,
                      unsigned long *out, double complex *value);
int ValidateSymmetryTerms(const struct DefineList *def);
int SymmetryUsesExtendedTerms(const struct DefineList *def);
/* Itinerant spinful fermions without local spins, with fixed Sz (Hubbard,
 * tJ), with only N fixed (HubbardNConserved, tJNConserved) or with N free
 * (HubbardGC). */
int IsSymmetrySpinfulFermionModel(int model);
/* Spinful models whose physical states exclude doubly occupied sites. */
int IsSymmetryTJModel(int model);
/* One orbital per site, with fixed N (SpinlessFermion) or N free. */
int IsSymmetrySpinlessFermionModel(int model);
/* Fermion sectors without a fixed particle number; ParticleParity may
 * select even or odd N. */
int IsSymmetryFermionGCModel(int model);

#endif
