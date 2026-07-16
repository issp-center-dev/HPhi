#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symmetry_basis.h"
#include "struct.h"

FILE *stdoutMPI = NULL;
long unsigned int *list_1 = NULL;
long unsigned int *list_2_1 = NULL;
long unsigned int *list_2_2 = NULL;
double *list_Diagonal = NULL;
int g_tj_odd_split_guard_enabled = 0;
long unsigned int g_tj_odd_split_up_mask = 0;
long unsigned int g_tj_odd_split_down_mask = 0;
static unsigned long int test_raw_dim = 0;

static int perm_storage[6][6];
static int anti_storage[6][6];
static int *perm_rows[6];
static int *anti_rows[6];
static double complex chars_storage[6];
static int exchange_storage[6][2];
static int *exchange_rows[6];
static double exchange_params[6];

static int popcount_ulong(unsigned long int x)
{
  int count = 0;
  while (x != 0UL) {
    count += (int)(x & 1UL);
    x >>= 1U;
  }
  return count;
}

static unsigned long int setup_fixed_sz_basis(unsigned int nsite,
                                              unsigned int nup)
{
  unsigned long int state;
  unsigned long int dim = 0;
  unsigned long int limit = 1UL << nsite;
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
  for (state = 0; state < limit; state++) {
    if ((unsigned int)popcount_ulong(state) == nup) dim++;
  }
  list_1 = (unsigned long int *)calloc(dim + 1, sizeof(unsigned long int));
  list_Diagonal = (double *)calloc(dim + 1, sizeof(double));
  if (list_1 == NULL || list_Diagonal == NULL) {
    fprintf(stderr, "failed to allocate fixed-Sz basis\n");
    exit(1);
  }
  dim = 0;
  for (state = 0; state < limit; state++) {
    if ((unsigned int)popcount_ulong(state) == nup) {
      dim++;
      list_1[dim] = state;
    }
  }
  test_raw_dim = dim;
  return dim;
}

static unsigned long int raw_index_for_state(unsigned long int state)
{
  unsigned long int raw;
  for (raw = 1; raw <= test_raw_dim; raw++) {
    if (list_1[raw] == state) return raw;
  }
  return 0;
}

static double spin_ising_pair_energy(unsigned long int state,
                                     unsigned int site0,
                                     unsigned int site1,
                                     double coupling)
{
  unsigned long int bit0 = (state >> site0) & 1UL;
  unsigned long int bit1 = (state >> site1) & 1UL;
  if (bit0 == bit1) return 0.25 * coupling;
  return -0.25 * coupling;
}

static void set_ising_ring_diagonal(unsigned int nsite, double coupling)
{
  unsigned long int raw;
  for (raw = 1; raw <= test_raw_dim; raw++) {
    unsigned int site;
    double diagonal = 0.0;
    for (site = 0; site < nsite; site++) {
      diagonal += spin_ising_pair_energy(list_1[raw], site, (site + 1U) % nsite, coupling);
    }
    list_Diagonal[raw] = diagonal;
  }
}

static void setup_cyclic_def(struct DefineList *def,
                             unsigned int nsite,
                             unsigned int momentum_index)
{
  unsigned int g, s;
  const double pi = 3.141592653589793238462643383279502884;
  memset(def, 0, sizeof(*def));
  for (g = 0; g < nsite; g++) {
    perm_rows[g] = perm_storage[g];
    anti_rows[g] = anti_storage[g];
  }
  def->iFlgSymmetryBasis = TRUE;
  def->Nsite = nsite;
  def->NSymTrans = nsite;
  def->SymTrans = perm_rows;
  def->SymTransAnti = anti_rows;
  def->SymTransChar = chars_storage;
  for (g = 0; g < nsite; g++) {
    double angle = -2.0 * pi * (double)momentum_index * (double)g / (double)nsite;
    chars_storage[g] = cos(angle) + I * sin(angle);
    for (s = 0; s < nsite; s++) {
      perm_storage[g][s] = (int)((s + g) % nsite);
      anti_storage[g][s] = 1;
    }
  }
}

static void setup_c4_def(struct DefineList *def)
{
  setup_cyclic_def(def, 4, 0);
}

static void setup_exchange_ring(struct DefineList *def, unsigned int nsite)
{
  unsigned int site;
  def->NExchangeCoupling = nsite;
  def->ExchangeCoupling = exchange_rows;
  def->ParaExchangeCoupling = exchange_params;
  for (site = 0; site < nsite; site++) {
    exchange_rows[site] = exchange_storage[site];
    exchange_storage[site][0] = (int)site;
    exchange_storage[site][1] = (int)((site + 1U) % nsite);
    exchange_params[site] = 1.0;
  }
}

static void setup_bind(struct BindStruct *X,
                       unsigned int nsite,
                       unsigned int nup,
                       unsigned int momentum_index)
{
  memset(X, 0, sizeof(*X));
  setup_cyclic_def(&X->Def, nsite, momentum_index);
  setup_exchange_ring(&X->Def, nsite);
  X->Check.idim_max = setup_fixed_sz_basis(nsite, nup);
}

static double complex orbit_coeff_sum_for_state(unsigned long int raw_state,
                                                const struct DefineList *def,
                                                unsigned long int target_state)
{
  unsigned int g;
  double complex sum = 0.0;
  for (g = 0; g < def->NSymTrans; g++) {
    unsigned long int moved = SymmetryApplyToSpinBits(raw_state, def->SymTrans[g], def->Nsite);
    if (moved == target_state) sum += conj(def->SymTransChar[g]);
  }
  return sum;
}

static int apply_exchange_halfspin_test(unsigned long int state,
                                        int site0,
                                        int site1,
                                        unsigned long int *out_state)
{
  unsigned long int b0 = (state >> (unsigned int)site0) & 1UL;
  unsigned long int b1 = (state >> (unsigned int)site1) & 1UL;
  if (b0 == b1) return 0;
  *out_state = state ^ (1UL << (unsigned int)site0) ^ (1UL << (unsigned int)site1);
  return 1;
}

static double complex projected_coeff_for_state(const struct BindStruct *X,
                                                unsigned long int basis_index,
                                                unsigned long int state)
{
  unsigned int g;
  const struct SymmetryBasisVector *basis = &X->Sym->basis[basis_index];
  double complex sum = 0.0;
  for (g = 0; g < X->Def.NSymTrans; g++) {
    unsigned long int moved = SymmetryApplyToSpinBits(basis->rep_state,
                                                      X->Def.SymTrans[g],
                                                      X->Def.Nsite);
    if (moved == state) sum += conj(X->Def.SymTransChar[g]);
  }
  return sum / basis->norm;
}

static double complex raw_reference_matrix_element(const struct BindStruct *X,
                                                   unsigned long int alpha,
                                                   unsigned long int beta)
{
  unsigned long int raw;
  double complex value = 0.0;
  for (raw = 1; raw <= X->Check.idim_max; raw++) {
    unsigned int term;
    unsigned long int state = list_1[raw];
    double complex beta_coeff = projected_coeff_for_state(X, beta, state);
    if (cabs(beta_coeff) < 1.0e-12) continue;
    if (list_Diagonal != NULL) {
      double complex alpha_coeff = projected_coeff_for_state(X, alpha, state);
      value += conj(alpha_coeff) * list_Diagonal[raw] * beta_coeff;
    }
    for (term = 0; term < X->Def.NExchangeCoupling; term++) {
      unsigned long int out_state;
      if (apply_exchange_halfspin_test(state,
                                       X->Def.ExchangeCoupling[term][0],
                                       X->Def.ExchangeCoupling[term][1],
                                       &out_state) != 0) {
        double complex alpha_coeff = projected_coeff_for_state(X, alpha, out_state);
        value += conj(alpha_coeff) * X->Def.ParaExchangeCoupling[term] * beta_coeff;
      }
    }
  }
  return value;
}

static double complex canonicalized_matrix_element(const struct BindStruct *X,
                                                   unsigned long int alpha,
                                                   unsigned long int beta)
{
  unsigned int term;
  double complex value = 0.0;
  if (alpha == beta && X->Sym->sym_diagonal != NULL) {
    value += X->Sym->sym_diagonal[beta];
  }
  for (term = 0; term < X->Def.NExchangeCoupling; term++) {
    unsigned long int out_state;
    if (apply_exchange_halfspin_test(X->Sym->basis[beta].rep_state,
                                     X->Def.ExchangeCoupling[term][0],
                                     X->Def.ExchangeCoupling[term][1],
                                     &out_state) != 0) {
      struct SymmetryCanonicalResult result;
      if (SymmetryCanonicalizeSpinState(X, out_state, &result) != 0) {
        fprintf(stderr, "canonicalize returned an internal error\n");
        exit(1);
      }
      if (result.found != 0 && result.basis_index == alpha) {
        double norm_factor = X->Sym->basis[alpha].norm / X->Sym->basis[beta].norm;
        value += X->Def.ParaExchangeCoupling[term] * result.phase * norm_factor;
      }
    }
  }
  return value;
}

static void assert_ulong_eq(unsigned long int got,
                            unsigned long int expected,
                            const char *label)
{
  if (got != expected) {
    fprintf(stderr, "%s: got %lu expected %lu\n", label, got, expected);
    exit(1);
  }
}

static void assert_int_eq(int got, int expected, const char *label)
{
  if (got != expected) {
    fprintf(stderr, "%s: got %d expected %d\n", label, got, expected);
    exit(1);
  }
}

static void assert_complex_close(double complex got,
                                 double complex expected,
                                 double tol,
                                 const char *label)
{
  if (cabs(got - expected) > tol) {
    fprintf(stderr, "%s: got %.16e%+.16ei expected %.16e%+.16ei\n",
            label, creal(got), cimag(got), creal(expected), cimag(expected));
    exit(1);
  }
}

static void assert_symmetry_dim(unsigned int nsite,
                                unsigned int nup,
                                unsigned int momentum_index,
                                unsigned long int expected_dim,
                                const char *label)
{
  struct BindStruct X;
  setup_bind(&X, nsite, nup, momentum_index);
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: BuildSymmetryBasis failed\n", label);
    exit(1);
  }
  assert_ulong_eq(X.Sym->dim, expected_dim, label);
  FreeSymmetryBasis(X.Sym);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

static void assert_canonicalized_matrix_matches_raw(unsigned int nsite,
                                                    unsigned int nup,
                                                    unsigned int momentum_index,
                                                    double diagonal_coupling,
                                                    const char *label)
{
  struct BindStruct X;
  unsigned long int alpha, beta;
  setup_bind(&X, nsite, nup, momentum_index);
  if (diagonal_coupling != 0.0) set_ising_ring_diagonal(nsite, diagonal_coupling);
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: BuildSymmetryBasis failed\n", label);
    exit(1);
  }
  for (beta = 1; beta <= X.Sym->dim; beta++) {
    for (alpha = 1; alpha <= X.Sym->dim; alpha++) {
      double complex raw_value = raw_reference_matrix_element(&X, alpha, beta);
      double complex canonical_value = canonicalized_matrix_element(&X, alpha, beta);
      assert_complex_close(canonical_value, raw_value, 1.0e-10, label);
    }
  }
  FreeSymmetryBasis(X.Sym);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

static void assert_orbit_diagonal_is_representative(unsigned int nsite,
                                                    unsigned int nup,
                                                    unsigned int momentum_index,
                                                    double diagonal_coupling,
                                                    const char *label)
{
  struct BindStruct X;
  unsigned long int beta;
  setup_bind(&X, nsite, nup, momentum_index);
  set_ising_ring_diagonal(nsite, diagonal_coupling);
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: BuildSymmetryBasis failed\n", label);
    exit(1);
  }
  for (beta = 1; beta <= X.Sym->dim; beta++) {
    unsigned int g;
    double representative_diagonal = X.Sym->sym_diagonal[beta];
    for (g = 0; g < X.Def.NSymTrans; g++) {
      unsigned long int moved = SymmetryApplyToSpinBits(X.Sym->basis[beta].rep_state,
                                                        X.Def.SymTrans[g],
                                                        X.Def.Nsite);
      unsigned long int raw = raw_index_for_state(moved);
      assert_int_eq(raw != 0, 1, label);
      assert_complex_close(list_Diagonal[raw], representative_diagonal, 1.0e-12, label);
    }
  }
  FreeSymmetryBasis(X.Sym);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

int main(void)
{
  int shift4[4] = {1, 2, 3, 0};
  stdoutMPI = stderr;
  assert_ulong_eq(SymmetryApplyToSpinBits(0x1UL, shift4, 4), 0x2UL, "single bit shift");
  assert_ulong_eq(SymmetryApplyToSpinBits(0x9UL, shift4, 4), 0x3UL, "wrap shift");
  assert_ulong_eq(SymmetryApplyToSpinBits(0x6UL, shift4, 4), 0xcUL, "two bit shift");
  assert_int_eq(1, 1, "unit harness still running");
  {
    struct DefineList def;
    setup_c4_def(&def);
    assert_int_eq(ValidateSymmetryGroupInput(&def), 0, "C4 group validates");
  }
  {
    struct DefineList def;
    setup_c4_def(&def);
    def.NSymTrans = 2;
    assert_int_eq(ValidateSymmetryGroupInput(&def), -1, "generator-only C4 rejects");
  }
  {
    struct DefineList def;
    setup_c4_def(&def);
    def.SymTransAnti[1][0] = -1;
    assert_int_eq(ValidateSymmetryGroupInput(&def), -1, "anti boundary rejects");
  }
  {
    struct DefineList def;
    setup_c4_def(&def);
    assert_int_eq(cabs(orbit_coeff_sum_for_state(0x3UL, &def, 0x3UL) - 1.0) < 1.0e-12, 1,
                  "orbit coeff identity");
    assert_int_eq(cabs(orbit_coeff_sum_for_state(0x3UL, &def, 0x6UL) - 1.0) < 1.0e-12, 1,
                  "orbit coeff shifted state");
  }
  assert_symmetry_dim(4, 2, 0, 2, "C4 k=0 sector dimension");
  assert_symmetry_dim(4, 2, 1, 1, "C4 k=pi/2 sector dimension");
  assert_symmetry_dim(6, 3, 3, 4, "C6 k=pi sector dimension");
  assert_symmetry_dim(6, 3, 1, 3, "C6 k=pi/3 sector dimension");
  {
    struct BindStruct X;
    unsigned long int raw;
    setup_bind(&X, 4, 2, 0);
    if (BuildSymmetryBasis(&X) != 0) {
      fprintf(stderr, "C4 k=0 BuildSymmetryBasis failed\n");
      exit(1);
    }
    for (raw = 1; raw <= X.Check.idim_max; raw++) {
      struct SymmetryCanonicalResult result;
      unsigned long int expected_rep = list_1[raw];
      unsigned int g;
      for (g = 0; g < X.Def.NSymTrans; g++) {
        unsigned long int moved = SymmetryApplyToSpinBits(list_1[raw],
                                                          X.Def.SymTrans[g],
                                                          X.Def.Nsite);
        if (moved < expected_rep) expected_rep = moved;
      }
      assert_int_eq(SymmetryCanonicalizeSpinState(&X, list_1[raw], &result), 0,
                    "C4 k=0 canonicalize return");
      assert_int_eq(result.found, 1, "C4 k=0 canonicalize found");
      assert_ulong_eq(X.Sym->basis[result.basis_index].rep_state, expected_rep,
                      "C4 k=0 canonicalize representative");
    }
    FreeSymmetryBasis(X.Sym);
    free(list_1);
    free(list_Diagonal);
    list_1 = NULL;
    list_Diagonal = NULL;
  }
  {
    struct BindStruct X;
    struct SymmetryCanonicalResult result;
    setup_bind(&X, 4, 2, 1);
    if (BuildSymmetryBasis(&X) != 0) {
      fprintf(stderr, "C4 k=pi/2 BuildSymmetryBasis failed\n");
      exit(1);
    }
    assert_int_eq(SymmetryCanonicalizeSpinState(&X, 0x5UL, &result), 0,
                  "C4 k=pi/2 period-2 canonicalize return");
    assert_int_eq(result.found, 0, "C4 k=pi/2 period-2 orbit absent");
    assert_int_eq(SymmetryCanonicalizeSpinState(&X, 0xaUL, &result), 0,
                  "C4 k=pi/2 shifted period-2 canonicalize return");
    assert_int_eq(result.found, 0, "C4 k=pi/2 shifted period-2 orbit absent");
    FreeSymmetryBasis(X.Sym);
    free(list_1);
    free(list_Diagonal);
    list_1 = NULL;
    list_Diagonal = NULL;
  }
  {
    struct BindStruct X;
    struct SymmetryCanonicalResult result;
    unsigned long int rep = 0x7UL;
    unsigned long int shifted;
    setup_bind(&X, 6, 3, 1);
    if (BuildSymmetryBasis(&X) != 0) {
      fprintf(stderr, "C6 k=pi/3 BuildSymmetryBasis failed\n");
      exit(1);
    }
    shifted = SymmetryApplyToSpinBits(rep, X.Def.SymTrans[1], X.Def.Nsite);
    assert_int_eq(SymmetryCanonicalizeSpinState(&X, shifted, &result), 0,
                  "C6 k=pi/3 canonicalize return");
    assert_int_eq(result.found, 1, "C6 k=pi/3 canonicalize found");
    assert_ulong_eq(X.Sym->basis[result.basis_index].rep_state, rep,
                    "C6 k=pi/3 canonicalize representative");
    assert_complex_close(result.phase, X.Def.SymTransChar[1], 1.0e-12,
                         "C6 k=pi/3 canonicalize phase");
    FreeSymmetryBasis(X.Sym);
    free(list_1);
    free(list_Diagonal);
    list_1 = NULL;
    list_Diagonal = NULL;
  }
  assert_canonicalized_matrix_matches_raw(4, 2, 0, 0.0,
                                          "C4 k=0 canonicalized matrix matches raw reference");
  assert_canonicalized_matrix_matches_raw(6, 3, 1, 0.0,
                                          "C6 k=pi/3 canonicalized matrix matches raw reference");
  assert_orbit_diagonal_is_representative(6, 3, 1, 1.0,
                                          "C6 k=pi/3 Ising diagonal is orbit-invariant");
  assert_canonicalized_matrix_matches_raw(6, 3, 1, 1.0,
                                          "C6 k=pi/3 Ising canonicalized matrix matches raw reference");
  return 0;
}
