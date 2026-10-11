#include <limits.h>
#include <string.h>

#include "DefCommon.h"
#include "symmetry_state_enumerator.h"
#include "struct.h"
#include "symmetry_kondo.h"

static int initialize_binomial_table(
    struct SymmetryStateEnumerator *enumerator, unsigned int maximum)
{
  unsigned int n, k;
  if (enumerator == NULL ||
      maximum > HPHI_SYMMETRY_STATE_WORD_BITS) {
    return -1;
  }
  memset(enumerator->binomial, 0, sizeof(enumerator->binomial));
  enumerator->binomial[0][0] = 1UL;
  for (n = 1U; n <= maximum; n++) {
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


/* Dimension arithmetic has a signed-index bound; Pascal's triangle above
 * separately checks the unsigned intermediate coefficients. */
static int dimension_product(unsigned long a, unsigned long b, unsigned long *result)
{
  if (a > LONG_MAX || b > LONG_MAX || (b && a > (unsigned long)LONG_MAX / b))
    return -1;
  *result = a * b;
  return 0;
}

static unsigned long kondo_binomial(const struct SymmetryStateEnumerator *e,
                                    unsigned int n, int k)
{
  return k < 0 || (unsigned int)k > n ? 0UL : e->binomial[n][k];
}

static int kondo_completions(const struct SymmetryStateEnumerator *e, unsigned int l,
                             unsigned int c, int up, int down, int nc,
                             unsigned long *count)
{
  unsigned long total = 0, a, b, term;
  unsigned int u, exponent;
  if (e->model == Kondo) {
    if (up < 0 || down < 0 || (unsigned int)up > l + c || (unsigned int)down > l + c) {
      *count = 0;
      return 0;
    }
    for (u = 0; u <= l; ++u) {
      a = kondo_binomial(e, c, up - (int)u);
      b = kondo_binomial(e, c, down - (int)l + (int)u);
      if (!a || !b)
        continue;
      if (dimension_product(a, b, &term) ||
          dimension_product(term, e->binomial[l][u], &term) ||
          term > (unsigned long)LONG_MAX || total > (unsigned long)LONG_MAX - term)
        return -1;
      total += term;
    }
  } else {
    total = e->model == KondoGC ? 1UL : kondo_binomial(e, 2 * c, nc);
    exponent = l + (e->model == KondoGC ? 2 * c : 0U);
    for (u = 0; u < exponent; ++u)
      if (dimension_product(total, 2UL, &total))
        return -1;
  }
  *count = total;
  return 0;
}

int CountSymmetryKondoCompletions(int model, unsigned int local_sites,
                                  unsigned int conduction_sites, int nup, int ndown,
                                  int ncond, unsigned long *count)
{
  struct SymmetryStateEnumerator e;
  unsigned int maximum;
  if (!count || !IsSymmetryKondoModel(model) ||
      local_sites > HPHI_SYMMETRY_STATE_WORD_BITS ||
      conduction_sites > HPHI_SYMMETRY_STATE_WORD_BITS / 2U)
    return -1;
  memset(&e, 0, sizeof(e));
  e.model = model;
  maximum = local_sites > 2 * conduction_sites ? local_sites : 2 * conduction_sites;
  if (initialize_binomial_table(&e, maximum))
    return -1;
  return kondo_completions(&e, local_sites, conduction_sites, nup, ndown, ncond, count);
}

int ComputeSymmetryKondoDimension(const struct DefineList *def, unsigned long *dimension)
{
  unsigned long result;
  if (!dimension || ValidateSymmetryKondoSpace(def))
    return -1;
  if (CountSymmetryKondoCompletions(def->iCalcModel, def->NLocSpn,
                                    def->Nsite - def->NLocSpn, def->Nup, def->Ndown,
                                    def->NCond, &result) ||
      !result)
    return -1;
  *dimension = result;
  return 0;
}

static int initialize_kondo(const struct DefineList *def, unsigned long expected,
                            struct SymmetryStateEnumerator *e)
{
  unsigned int site, maximum, l = def->NLocSpn, c = def->Nsite - l;
  unsigned long dimension;
  if (ComputeSymmetryKondoDimension(def, &dimension) || dimension != expected ||
      SymmetryKondoLocalMask(def, &e->local_site_mask))
    return -1;
  e->bit_count = 2 * def->Nsite;
  e->ncond = def->NCond;
  e->raw_dim = dimension;
  maximum = l > 2 * c ? l : 2 * c;
  if (initialize_binomial_table(e, maximum))
    return -1;
  for (site = 0; site < e->nsite; ++site) {
    unsigned int local = (unsigned int)((e->local_site_mask >> site) & 1UL);
    e->prefix_local[site + 1] = e->prefix_local[site] + local;
    e->prefix_conduction[site + 1] = e->prefix_conduction[site] + 1U - local;
  }
  return 0;
}

static int kondo_state_at(const struct SymmetryStateEnumerator *e, unsigned long rank,
                          unsigned long *state)
{
  unsigned int site = e->nsite;
  int up = (int)e->nup, down = (int)e->ndown, nc = (int)e->ncond;
  unsigned long word = 0;
  while (site-- > 0U) {
    unsigned int digit, local = (unsigned int)((e->local_site_mask >> site) & 1UL);
    for (digit = 0; digit < 4; ++digit) {
      unsigned long completions;
      int next_up = up - (int)(digit & 1U), next_down = down - (int)(digit >> 1);
      int next_nc = nc - (local ? 0 : (int)((digit & 1U) + (digit >> 1)));
      if (local && (digit == 0 || digit == 3))
        continue;
      if (kondo_completions(e, e->prefix_local[site], e->prefix_conduction[site], next_up,
                            next_down, next_nc, &completions))
        return -1;
      if (rank > completions) {
        rank -= completions;
        continue;
      }
      word |= (unsigned long)digit << (2 * site);
      up = next_up;
      down = next_down;
      nc = next_nc;
      break;
    }
    if (digit == 4)
      return -1;
  }
  if (rank != 1 || (e->model == Kondo && (up || down)) ||
      (e->model == KondoNConserved && nc))
    return -1;
  *state = word;
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

  if (IsSymmetryKondoModel(def->iCalcModel)) {
    if (initialize_kondo(def, expected_raw_dim, &initialized)) return -1;
    *enumerator = initialized;
    return 0;
  }

  if (def->iCalcModel == SpinGC) {
    if (def->iFlgGeneralSpin != FALSE || def->Nsite >= word_bits) return -1;
    initialized.bit_count = def->Nsite;
    initialized.nup = initialized.ndown = 0U;
    initialized.raw_dim = 1UL << def->Nsite;
    if (initialized.raw_dim != expected_raw_dim) return -1;
    *enumerator = initialized;
    return 0;
  }

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
  case tJ:
    if (def->Nsite > word_bits / 2U ||
        def->Nup > def->Nsite ||
        def->Ndown > def->Nsite ||
        def->Ne != def->Nup + def->Ndown ||
        (def->iCalcModel == tJ && def->Ne > def->Nsite)) {
      return -1;
    }
    initialized.bit_count = 2U * def->Nsite;
    break;
  default:
    return -1;
  }
  if (initialize_binomial_table(&initialized, initialized.nsite) != 0) return -1;
  up_dim = enumerator_binomial(
      &initialized, initialized.nsite, initialized.nup);
  if (initialized.model == Hubbard || initialized.model == tJ) {
    down_dim = enumerator_binomial(
        &initialized, initialized.nsite - (initialized.model == tJ ? initialized.nup : 0U),
        initialized.ndown);
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

static int tj_completion_count(const struct SymmetryStateEnumerator *enumerator,
                                unsigned int sites, unsigned int up,
                                unsigned int down, unsigned long *count)
{
  if (up > sites || down > sites - up) { *count = 0UL; return 0; }
  return checked_product(enumerator_binomial(enumerator, sites, up),
                         enumerator_binomial(enumerator, sites - up, down), count);
}

/* Lexicographic base-4 state order, with local digits 0 (empty), 1 (up),
 * 2 (down). Count completions instead of constructing a Hubbard list and
 * removing doubly occupied states. This also covers empty/full sectors. */
static int tj_state_at(const struct SymmetryStateEnumerator *enumerator,
                       unsigned long rank, unsigned long *state)
{
  unsigned int site = enumerator->nsite;
  unsigned int up = enumerator->nup, down = enumerator->ndown;
  unsigned long result = 0UL;
  while (site-- > 0U) {
    unsigned long empty_count, up_count = 0UL;
    if (tj_completion_count(enumerator, site, up, down, &empty_count)) return -1;
    if (rank < empty_count) continue;
    rank -= empty_count;
    if (up && tj_completion_count(enumerator, site, up - 1U, down, &up_count)) return -1;
    if (rank < up_count) {
      result |= 1UL << (2U * site);
      --up;
    } else {
      rank -= up_count;
      if (!down) return -1;
      result |= 1UL << (2U * site + 1U);
      --down;
    }
  }
  if (rank || up || down) return -1;
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
  if (IsSymmetryKondoModel(enumerator->model)) {
    if (!enumerator->nsite || enumerator->nsite > word_bits / 2U ||
        (enumerator->model == KondoGC && enumerator->nsite == word_bits / 2U) ||
        enumerator->bit_count != 2U * enumerator->nsite) return -1;
    return kondo_state_at(enumerator, raw_index, state);
  }
  if (enumerator->model == SpinGC) {
    if (enumerator->nsite == 0U || enumerator->nsite >= word_bits ||
        enumerator->bit_count != enumerator->nsite ||
        enumerator->raw_dim != (1UL << enumerator->nsite)) {
      return -1;
    }
    *state = raw_index - 1UL;
    return 0;
  }
  if (enumerator->model == Hubbard || enumerator->model == tJ) {
    if (enumerator->nsite == 0U ||
        enumerator->nsite > word_bits / 2U ||
        enumerator->bit_count != 2U * enumerator->nsite ||
        enumerator->nup > enumerator->nsite ||
        enumerator->ndown > enumerator->nsite ||
        (enumerator->model == tJ && enumerator->nup > enumerator->nsite - enumerator->ndown)) {
      return -1;
    }
    status = enumerator->model == tJ
        ? tj_state_at(enumerator, raw_index - 1UL, &result)
        : hubbard_state_at(enumerator, raw_index - 1UL, &result);
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
