#include <limits.h>
#include <string.h>

#include "DefCommon.h"
#include "symmetry_state_enumerator.h"
#include "struct.h"

static int initialize_binomial_table(
    struct SymmetryStateEnumerator *enumerator)
{
  unsigned int n, k;
  if (enumerator == NULL ||
      enumerator->nsite > HPHI_SYMMETRY_STATE_WORD_BITS) {
    return -1;
  }
  memset(enumerator->binomial, 0, sizeof(enumerator->binomial));
  enumerator->binomial[0][0] = 1UL;
  for (n = 1U; n <= enumerator->nsite; n++) {
    enumerator->binomial[n][0] = 1UL;
    enumerator->binomial[n][n] = 1UL;
    for (k = 1U; k < n; k++) {
      if (enumerator->binomial[n - 1U][k] >
          ULONG_MAX - enumerator->binomial[n - 1U][k - 1U]) {
        return -1;
      }
      enumerator->binomial[n][k] =
          enumerator->binomial[n - 1U][k] +
          enumerator->binomial[n - 1U][k - 1U];
    }
  }
  return 0;
}

static unsigned long int enumerator_binomial(
    const struct SymmetryStateEnumerator *enumerator,
    unsigned int n,
    unsigned int k)
{
  if (k > n) return 0UL;
  return enumerator->binomial[n][k];
}

static int checked_product(unsigned long int lhs,
                           unsigned long int rhs,
                           unsigned long int *value)
{
  if (value == NULL || (rhs != 0UL && lhs > ULONG_MAX / rhs)) return -1;
  *value = lhs * rhs;
  return 0;
}

int InitSymmetryStateEnumerator(
    const struct DefineList *def,
    unsigned long int expected_raw_dim,
    struct SymmetryStateEnumerator *enumerator)
{
  const unsigned int word_bits =
      (unsigned int)(CHAR_BIT * sizeof(unsigned long int));
  struct SymmetryStateEnumerator initialized;
  unsigned long int up_dim, down_dim, raw_dim;
  if (def == NULL || enumerator == NULL || def->Nsite == 0U) return -1;
  memset(&initialized, 0, sizeof(initialized));
  initialized.model = def->iCalcModel;
  initialized.nsite = def->Nsite;
  initialized.nup = def->Nup;
  initialized.ndown = def->Ndown;

  switch (def->iCalcModel) {
  case Spin:
    if (def->iFlgGeneralSpin != FALSE ||
        def->Nsite > word_bits ||
        def->Nup > def->Nsite ||
        (def->iFlgSzConserved != TRUE &&
         (def->Ndown > def->Nsite ||
          def->Nup != def->Nsite - def->Ndown))) {
      return -1;
    }
    initialized.bit_count = def->Nsite;
    break;
  case SpinlessFermion:
    if (def->Ne > def->Nsite || def->Nsite > word_bits) return -1;
    initialized.nup = def->Ne;
    initialized.ndown = 0U;
    initialized.bit_count = def->Nsite;
    break;
  case Hubbard:
    if (def->Nsite > word_bits / 2U ||
        def->Nup > def->Nsite ||
        def->Ndown > def->Nsite ||
        def->Ne != def->Nup + def->Ndown) {
      return -1;
    }
    initialized.bit_count = 2U * def->Nsite;
    break;
  default:
    return -1;
  }
  if (initialize_binomial_table(&initialized) != 0) return -1;
  up_dim = enumerator_binomial(
      &initialized, initialized.nsite, initialized.nup);
  if (initialized.model == Hubbard) {
    down_dim = enumerator_binomial(
        &initialized, initialized.nsite, initialized.ndown);
    if (checked_product(up_dim, down_dim, &raw_dim) != 0) return -1;
  } else {
    raw_dim = up_dim;
  }
  if (raw_dim == 0UL || raw_dim != expected_raw_dim) return -1;
  initialized.raw_dim = raw_dim;
  *enumerator = initialized;
  return 0;
}

static int fixed_popcount_state_at(
    const struct SymmetryStateEnumerator *enumerator,
    unsigned long int rank,
    unsigned long int *state)
{
  unsigned int position = enumerator->nsite;
  unsigned int remaining = enumerator->nup;
  unsigned long int result = 0UL;
  while (position-- > 0U) {
    unsigned long int zero_count =
        enumerator_binomial(enumerator, position, remaining);
    if (rank >= zero_count) {
      result |= 1UL << position;
      rank -= zero_count;
      if (remaining == 0U) return -1;
      remaining--;
    }
  }
  if (rank != 0UL || remaining != 0U) return -1;
  *state = result;
  return 0;
}

static int hubbard_state_at(
    const struct SymmetryStateEnumerator *enumerator,
    unsigned long int rank,
    unsigned long int *state)
{
  unsigned int position = enumerator->bit_count;
  unsigned int remaining_up = enumerator->nup;
  unsigned int remaining_down = enumerator->ndown;
  unsigned long int result = 0UL;
  while (position-- > 0U) {
    unsigned int lower_up = (position + 1U) / 2U;
    unsigned int lower_down = position / 2U;
    unsigned int *remaining =
        (position % 2U == 0U) ? &remaining_up : &remaining_down;
    unsigned long int up_count, down_count, zero_count;
    up_count = enumerator_binomial(enumerator, lower_up, remaining_up);
    down_count = enumerator_binomial(enumerator, lower_down, remaining_down);
    if (checked_product(up_count, down_count, &zero_count) != 0) {
      return -1;
    }
    if (rank >= zero_count) {
      result |= 1UL << position;
      rank -= zero_count;
      if (*remaining == 0U) return -1;
      (*remaining)--;
    }
  }
  if (rank != 0UL || remaining_up != 0U || remaining_down != 0U) return -1;
  *state = result;
  return 0;
}

int SymmetryStateEnumeratorStateAt(
    const struct SymmetryStateEnumerator *enumerator,
    unsigned long int raw_index,
    unsigned long int *state)
{
  const unsigned int word_bits =
      (unsigned int)(CHAR_BIT * sizeof(unsigned long int));
  unsigned long int result;
  int status;
  if (enumerator == NULL || state == NULL ||
      enumerator->raw_dim == 0UL ||
      raw_index == 0UL || raw_index > enumerator->raw_dim) {
    return -1;
  }
  if (enumerator->model == Hubbard) {
    if (enumerator->nsite == 0U ||
        enumerator->nsite > word_bits / 2U ||
        enumerator->bit_count != 2U * enumerator->nsite ||
        enumerator->nup > enumerator->nsite ||
        enumerator->ndown > enumerator->nsite) {
      return -1;
    }
    status = hubbard_state_at(enumerator, raw_index - 1UL, &result);
  } else if (enumerator->model == Spin ||
             enumerator->model == SpinlessFermion) {
    if (enumerator->nsite == 0U ||
        enumerator->nsite > word_bits ||
        enumerator->bit_count != enumerator->nsite ||
        enumerator->nup > enumerator->nsite) {
      return -1;
    }
    status = fixed_popcount_state_at(enumerator, raw_index - 1UL, &result);
  } else {
    return -1;
  }
  if (status != 0) return -1;
  *state = result;
  return 0;
}
