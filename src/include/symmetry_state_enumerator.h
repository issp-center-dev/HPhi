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
  unsigned long int raw_dim;
  unsigned long int
      binomial[HPHI_SYMMETRY_STATE_WORD_BITS + 1U]
              [HPHI_SYMMETRY_STATE_WORD_BITS + 1U];
};

int InitSymmetryStateEnumerator(
    const struct DefineList *def,
    unsigned long int expected_raw_dim,
    struct SymmetryStateEnumerator *enumerator);
int SymmetryStateEnumeratorStateAt(
    const struct SymmetryStateEnumerator *enumerator,
    unsigned long int raw_index,
    unsigned long int *state);

#endif /* HPHI_SYMMETRY_STATE_ENUMERATOR_H */
