#ifndef HPHI_SYMMETRY_STATE_ENUMERATOR_H
#define HPHI_SYMMETRY_STATE_ENUMERATOR_H

#include <limits.h>

#include "Common.h"

struct DefineList;

#define HPHI_SYMMETRY_STATE_WORD_BITS (CHAR_BIT * sizeof(unsigned long int))

struct SymmetryStateEnumerator {
  int model;
  unsigned int nsite;
  unsigned int nup;
  unsigned int ndown;
  unsigned int bit_count;
  unsigned long local_site_mask;
  unsigned int ncond;
  unsigned int prefix_local[HPHI_SYMMETRY_STATE_WORD_BITS / 2U + 1U];
  unsigned int prefix_conduction[HPHI_SYMMETRY_STATE_WORD_BITS / 2U + 1U];
  unsigned long int raw_dim;
  unsigned long int
      binomial[HPHI_SYMMETRY_STATE_WORD_BITS + 1U]
              [HPHI_SYMMETRY_STATE_WORD_BITS + 1U];
};

/* Validated, nonzero physical dimension bounded by LONG_MAX. */
int ComputeSymmetryKondoDimension(const struct DefineList *def,
                                  unsigned long *dimension);
/* Counts need no LocSpn array. Impossible residuals return count=0;
 * invalid models/site counts or checked-integer overflow return -1. */
int CountSymmetryKondoCompletions(int model, unsigned int local_sites,
    unsigned int conduction_sites, int nup, int ndown, int ncond,
    unsigned long *count);

int InitSymmetryStateEnumerator(
    const struct DefineList *def,
    unsigned long int expected_raw_dim,
    struct SymmetryStateEnumerator *enumerator);
int SymmetryStateEnumeratorStateAt(
    const struct SymmetryStateEnumerator *enumerator,
    unsigned long int raw_index,
    unsigned long int *state);

#endif /* HPHI_SYMMETRY_STATE_ENUMERATOR_H */
