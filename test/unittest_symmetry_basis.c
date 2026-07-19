#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "DefCommon.h"
#include "symmetry_basis.h"
#include "symmetry_matvec_plan.h"
#include "struct.h"

FILE *stdoutMPI = NULL;
int nproc = 1;
int myrank = 0;
long unsigned int *list_1 = NULL;
long unsigned int *list_2_1 = NULL;
long unsigned int *list_2_2 = NULL;
double *list_Diagonal = NULL;
int g_tj_odd_split_guard_enabled = 0;
long unsigned int g_tj_odd_split_up_mask = 0;
long unsigned int g_tj_odd_split_down_mask = 0;
static unsigned long int test_raw_dim = 0;

void StartTimer(int timer_id)
{
  (void)timer_id;
}

void StopTimer(int timer_id)
{
  (void)timer_id;
}

int BcastMPI_i(int root, int value)
{
  (void)root;
  return value;
}

static int perm_storage[6][6];
static int anti_storage[6][6];
static int *perm_rows[6];
static int *anti_rows[6];
static double complex chars_storage[6];
static int exchange_storage[6][2];
static int *exchange_rows[6];
static double exchange_params[6];
static int transfer_storage[24][4];
static int *transfer_rows[24];
static double complex transfer_params[24];
static int coulomb_intra_storage[8][1];
static int *coulomb_intra_rows[8];
static double coulomb_intra_params[8];
static int coulomb_storage[12][2];
static int *coulomb_rows[12];
static double coulomb_params[12];

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

static int count_hubbard_spin(unsigned long int state,
                              unsigned int nsite,
                              unsigned int spin)
{
  unsigned int site;
  int count = 0;
  for (site = 0; site < nsite; site++) {
    count += (int)((state >> (2U * site + spin)) & 1UL);
  }
  return count;
}

static unsigned long int setup_hubbard_basis(unsigned int nsite,
                                             unsigned int nup,
                                             unsigned int ndown)
{
  unsigned long int state;
  unsigned long int dim = 0;
  unsigned long int limit = 1UL << (2U * nsite);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
  for (state = 0; state < limit; state++) {
    if ((unsigned int)count_hubbard_spin(state, nsite, 0) == nup &&
        (unsigned int)count_hubbard_spin(state, nsite, 1) == ndown) {
      dim++;
    }
  }
  list_1 = (unsigned long int *)calloc(dim + 1, sizeof(unsigned long int));
  list_Diagonal = (double *)calloc(dim + 1, sizeof(double));
  if (list_1 == NULL || list_Diagonal == NULL) {
    fprintf(stderr, "failed to allocate Hubbard basis\n");
    exit(1);
  }
  dim = 0;
  for (state = 0; state < limit; state++) {
    if ((unsigned int)count_hubbard_spin(state, nsite, 0) == nup &&
        (unsigned int)count_hubbard_spin(state, nsite, 1) == ndown) {
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

static double spinless_coulomb_pair_energy(unsigned long int state,
                                           unsigned int site0,
                                           unsigned int site1,
                                           double coupling)
{
  unsigned long int bit0 = (state >> site0) & 1UL;
  unsigned long int bit1 = (state >> site1) & 1UL;
  return (bit0 != 0UL && bit1 != 0UL) ? coupling : 0.0;
}

static void set_spinless_coulomb_ring_diagonal(unsigned int nsite, double coupling)
{
  unsigned long int raw;
  for (raw = 1; raw <= test_raw_dim; raw++) {
    unsigned int site;
    double diagonal = 0.0;
    for (site = 0; site < nsite; site++) {
      diagonal += spinless_coulomb_pair_energy(list_1[raw], site, (site + 1U) % nsite, coupling);
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
  def->iCalcModel = Spin;
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

static void setup_spinless_bind(struct BindStruct *X,
                                unsigned int nsite,
                                unsigned int ne,
                                unsigned int momentum_index)
{
  memset(X, 0, sizeof(*X));
  setup_cyclic_def(&X->Def, nsite, momentum_index);
  X->Def.iCalcModel = SpinlessFermion;
  X->Def.Ne = ne;
  X->Def.Nup = ne;
  X->Def.Ndown = 0;
  X->Check.idim_max = setup_fixed_sz_basis(nsite, ne);
}

static void setup_spinless_transfer_ring(struct DefineList *def, unsigned int nsite)
{
  unsigned int site;
  def->EDNTransfer = 2U * nsite;
  def->EDGeneralTransfer = transfer_rows;
  def->EDParaGeneralTransfer = transfer_params;
  for (site = 0; site < nsite; site++) {
    unsigned int next = (site + 1U) % nsite;
    unsigned int forward = 2U * site;
    unsigned int reverse = forward + 1U;
    transfer_rows[forward] = transfer_storage[forward];
    transfer_storage[forward][0] = (int)site;
    transfer_storage[forward][1] = 0;
    transfer_storage[forward][2] = (int)next;
    transfer_storage[forward][3] = 0;
    transfer_params[forward] = 1.0;
    transfer_rows[reverse] = transfer_storage[reverse];
    transfer_storage[reverse][0] = (int)next;
    transfer_storage[reverse][1] = 0;
    transfer_storage[reverse][2] = (int)site;
    transfer_storage[reverse][3] = 0;
    transfer_params[reverse] = 1.0;
  }
}

static void setup_hubbard_bind(struct BindStruct *X,
                               unsigned int nsite,
                               unsigned int nup,
                               unsigned int ndown,
                               unsigned int momentum_index)
{
  memset(X, 0, sizeof(*X));
  setup_cyclic_def(&X->Def, nsite, momentum_index);
  X->Def.iCalcModel = Hubbard;
  X->Def.Nup = nup;
  X->Def.Ndown = ndown;
  X->Def.Ne = nup + ndown;
  X->Check.idim_max = setup_hubbard_basis(nsite, nup, ndown);
}

static void setup_hubbard_transfer_ring(struct DefineList *def, unsigned int nsite)
{
  unsigned int site, spin;
  unsigned int row = 0;
  def->NTransfer = 4U * nsite;
  def->EDNTransfer = 4U * nsite;
  def->GeneralTransfer = transfer_rows;
  def->EDGeneralTransfer = transfer_rows;
  def->ParaGeneralTransfer = transfer_params;
  def->EDParaGeneralTransfer = transfer_params;
  for (site = 0; site < nsite; site++) {
    unsigned int next = (site + 1U) % nsite;
    for (spin = 0; spin < 2U; spin++) {
      transfer_rows[row] = transfer_storage[row];
      transfer_storage[row][0] = (int)site;
      transfer_storage[row][1] = (int)spin;
      transfer_storage[row][2] = (int)next;
      transfer_storage[row][3] = (int)spin;
      transfer_params[row] = 1.0;
      row++;

      transfer_rows[row] = transfer_storage[row];
      transfer_storage[row][0] = (int)next;
      transfer_storage[row][1] = (int)spin;
      transfer_storage[row][2] = (int)site;
      transfer_storage[row][3] = (int)spin;
      transfer_params[row] = 1.0;
      row++;
    }
  }
}

static void set_hubbard_coulomb_intra_diagonal(unsigned int nsite, double coupling)
{
  unsigned long int raw;
  for (raw = 1; raw <= test_raw_dim; raw++) {
    unsigned int site;
    double diagonal = 0.0;
    for (site = 0; site < nsite; site++) {
      unsigned long int up = (list_1[raw] >> (2U * site)) & 1UL;
      unsigned long int down = (list_1[raw] >> (2U * site + 1U)) & 1UL;
      diagonal += coupling * (double)(up * down);
    }
    list_Diagonal[raw] = diagonal;
  }
}

static void setup_hubbard_coulomb_intra(struct DefineList *def,
                                        unsigned int nsite,
                                        double coupling)
{
  unsigned int site;
  def->NCoulombIntra = nsite;
  def->CoulombIntra = coulomb_intra_rows;
  def->ParaCoulombIntra = coulomb_intra_params;
  for (site = 0; site < nsite; site++) {
    coulomb_intra_rows[site] = coulomb_intra_storage[site];
    coulomb_intra_storage[site][0] = (int)site;
    coulomb_intra_params[site] = coupling;
  }
  set_hubbard_coulomb_intra_diagonal(nsite, coupling);
}

static void setup_spinless_coulomb_ring(struct DefineList *def, unsigned int nsite, double coupling)
{
  unsigned int site;
  def->NCoulombInter = nsite;
  def->CoulombInter = coulomb_rows;
  def->ParaCoulombInter = coulomb_params;
  for (site = 0; site < nsite; site++) {
    coulomb_rows[site] = coulomb_storage[site];
    coulomb_storage[site][0] = (int)site;
    coulomb_storage[site][1] = (int)((site + 1U) % nsite);
    coulomb_params[site] = coupling;
  }
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

static unsigned long int mask_between_sites_test(unsigned int site0, unsigned int site1)
{
  unsigned int lo = site0 < site1 ? site0 : site1;
  unsigned int hi = site0 < site1 ? site1 : site0;
  if (hi <= lo + 1U) return 0UL;
  return (1UL << hi) - (1UL << (lo + 1U));
}

static int apply_spinless_hopping_hermite_test(unsigned long int state,
                                               unsigned int site1,
                                               unsigned int site2,
                                               double complex trans,
                                               unsigned long int *out_state,
                                               double complex *hval)
{
  unsigned long int mask1 = 1UL << site1;
  unsigned long int mask2 = 1UL << site2;
  unsigned long int occupied1 = state & mask1;
  unsigned long int occupied2 = state & mask2;
  int sgn = (popcount_ulong(state & mask_between_sites_test(site1, site2)) % 2 == 0) ? 1 : -1;
  if (site1 == site2) return 0;
  if ((occupied1 == 0UL && occupied2 == 0UL) ||
      (occupied1 != 0UL && occupied2 != 0UL)) {
    return 0;
  }
  *out_state = state ^ mask1 ^ mask2;
  if (occupied1 != 0UL && occupied2 == 0UL) {
    *hval = (double)sgn * conj(trans);
  } else {
    *hval = (double)sgn * trans;
  }
  return 1;
}

static int apply_hubbard_hopping_hermite_test(unsigned long int state,
                                              unsigned int site1,
                                              unsigned int spin1,
                                              unsigned int site2,
                                              unsigned int spin2,
                                              double complex trans,
                                              unsigned long int *out_state,
                                              double complex *hval)
{
  if (spin1 > 1U || spin2 > 1U) return 0;
  return apply_spinless_hopping_hermite_test(state,
                                             2U * site1 + spin1,
                                             2U * site2 + spin2,
                                             trans, out_state, hval);
}

static double complex projected_coeff_for_state(const struct BindStruct *X,
                                                unsigned long int basis_index,
                                                unsigned long int state)
{
  unsigned int g;
  const struct SymmetryBasisVector *basis = &X->Sym->basis[basis_index];
  double complex sum = 0.0;
  for (g = 0; g < X->Def.NSymTrans; g++) {
    struct SymmetryTransformResult moved;
    if (SymmetryApplyToState(&X->Def, basis->rep_state, g, &moved) != 0) {
      fprintf(stderr, "state transform returned an internal error\n");
      exit(1);
    }
    if (moved.state == state) sum += conj(X->Def.SymTransChar[g]) * moved.amplitude;
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
    if (X->Def.iCalcModel == Spin) {
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
    } else if (X->Def.iCalcModel == SpinlessFermion) {
      for (term = 0; term < X->Def.EDNTransfer; term += 2U) {
        unsigned long int out_state;
        double complex hval;
        double complex trans = -X->Def.EDParaGeneralTransfer[term];
        if (apply_spinless_hopping_hermite_test(state,
                                                (unsigned int)X->Def.EDGeneralTransfer[term][0],
                                                (unsigned int)X->Def.EDGeneralTransfer[term][2],
                                                trans, &out_state, &hval) != 0) {
          double complex alpha_coeff = projected_coeff_for_state(X, alpha, out_state);
          value += conj(alpha_coeff) * hval * beta_coeff;
        }
      }
    } else if (X->Def.iCalcModel == Hubbard) {
      for (term = 0; term < X->Def.EDNTransfer; term += 2U) {
        unsigned long int out_state;
        double complex hval;
        double complex trans = -X->Def.EDParaGeneralTransfer[term];
        if (apply_hubbard_hopping_hermite_test(state,
                                               (unsigned int)X->Def.EDGeneralTransfer[term][0],
                                               (unsigned int)X->Def.EDGeneralTransfer[term][1],
                                               (unsigned int)X->Def.EDGeneralTransfer[term][2],
                                               (unsigned int)X->Def.EDGeneralTransfer[term][3],
                                               trans, &out_state, &hval) != 0) {
          double complex alpha_coeff = projected_coeff_for_state(X, alpha, out_state);
          value += conj(alpha_coeff) * hval * beta_coeff;
        }
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
  if (X->Def.iCalcModel == Spin) {
    for (term = 0; term < X->Def.NExchangeCoupling; term++) {
      unsigned long int out_state;
      if (apply_exchange_halfspin_test(X->Sym->basis[beta].rep_state,
                                       X->Def.ExchangeCoupling[term][0],
                                       X->Def.ExchangeCoupling[term][1],
                                       &out_state) != 0) {
        struct SymmetryCanonicalResult result;
        if (SymmetryCanonicalizeState(X, out_state, &result) != 0) {
          fprintf(stderr, "canonicalize returned an internal error\n");
          exit(1);
        }
        if (result.found != 0 && result.basis_index == alpha) {
          double norm_factor = X->Sym->basis[alpha].norm / X->Sym->basis[beta].norm;
          value += X->Def.ParaExchangeCoupling[term] * result.phase * norm_factor;
        }
      }
    }
  } else if (X->Def.iCalcModel == SpinlessFermion) {
    for (term = 0; term < X->Def.EDNTransfer; term += 2U) {
      unsigned long int out_state;
      double complex hval;
      double complex trans = -X->Def.EDParaGeneralTransfer[term];
      if (apply_spinless_hopping_hermite_test(X->Sym->basis[beta].rep_state,
                                              (unsigned int)X->Def.EDGeneralTransfer[term][0],
                                              (unsigned int)X->Def.EDGeneralTransfer[term][2],
                                              trans, &out_state, &hval) != 0) {
        struct SymmetryCanonicalResult result;
        if (SymmetryCanonicalizeState(X, out_state, &result) != 0) {
          fprintf(stderr, "canonicalize returned an internal error\n");
          exit(1);
        }
        if (result.found != 0 && result.basis_index == alpha) {
          double norm_factor = X->Sym->basis[alpha].norm / X->Sym->basis[beta].norm;
          value += hval * result.phase * norm_factor;
        }
      }
    }
  } else if (X->Def.iCalcModel == Hubbard) {
    for (term = 0; term < X->Def.EDNTransfer; term += 2U) {
      unsigned long int out_state;
      double complex hval;
      double complex trans = -X->Def.EDParaGeneralTransfer[term];
      if (apply_hubbard_hopping_hermite_test(X->Sym->basis[beta].rep_state,
                                             (unsigned int)X->Def.EDGeneralTransfer[term][0],
                                             (unsigned int)X->Def.EDGeneralTransfer[term][1],
                                             (unsigned int)X->Def.EDGeneralTransfer[term][2],
                                             (unsigned int)X->Def.EDGeneralTransfer[term][3],
                                             trans, &out_state, &hval) != 0) {
        struct SymmetryCanonicalResult result;
        if (SymmetryCanonicalizeState(X, out_state, &result) != 0) {
          fprintf(stderr, "canonicalize returned an internal error\n");
          exit(1);
        }
        if (result.found != 0 && result.basis_index == alpha) {
          double norm_factor = X->Sym->basis[alpha].norm / X->Sym->basis[beta].norm;
          value += hval * result.phase * norm_factor;
        }
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

struct LegacyVectorContext {
  double complex *output;
  double complex input_amp;
  unsigned long int dim;
};

static int accumulate_legacy_vector_entry(unsigned long int out_index,
                                          double complex coefficient,
                                          void *context)
{
  struct LegacyVectorContext *legacy = (struct LegacyVectorContext *)context;
  if (out_index == 0UL || out_index > legacy->dim) return -1;
  legacy->output[out_index] += coefficient * legacy->input_amp;
  return 0;
}

static void assert_plan_matches_canonicalized_matrix(struct BindStruct *X,
                                                     int require_duplicate,
                                                     const char *label)
{
  unsigned long int alpha, beta;
  size_t p;
  size_t matrix_size;
  int duplicate_found = 0;
  double difference_norm2 = 0.0;
  double legacy_norm2 = 0.0;
  double complex expected_prdct = 0.0;
  double complex plan_prdct = 0.0;
  double complex *dense;
  double complex *input;
  double complex *legacy_output;
  double complex *output;
  unsigned int *multiplicity;
  struct SymmetryMatvecPlan *plan;

  if (ActivateSymmetryBasisDimension(X) != 0 || BuildSymmetryMatvecPlan(X) != 0) {
    fprintf(stderr, "%s: plan setup failed\n", label);
    exit(1);
  }
  plan = X->Sym->matvec_plan;
  assert_int_eq(plan != NULL && plan->ready == TRUE, 1, label);
  assert_ulong_eq(plan->dim, X->Sym->dim, label);
  assert_ulong_eq(plan->local_offset, 0UL, label);
  assert_ulong_eq(plan->local_dim, X->Sym->dim, label);

  if (plan->dim > SIZE_MAX / plan->dim) {
    fprintf(stderr, "%s: dense matrix size overflow\n", label);
    exit(1);
  }
  matrix_size = (size_t)plan->dim * (size_t)plan->dim;
  dense = (double complex *)calloc(matrix_size, sizeof(*dense));
  multiplicity = (unsigned int *)calloc(matrix_size, sizeof(*multiplicity));
  input = (double complex *)calloc((size_t)plan->dim + 1U, sizeof(*input));
  legacy_output = (double complex *)calloc((size_t)plan->dim + 1U,
                                           sizeof(*legacy_output));
  output = (double complex *)calloc((size_t)plan->dim + 1U, sizeof(*output));
  if (dense == NULL || multiplicity == NULL || input == NULL ||
      legacy_output == NULL || output == NULL) {
    fprintf(stderr, "%s: dense test allocation failed\n", label);
    exit(1);
  }

  for (alpha = 1UL; alpha <= plan->dim; alpha++) {
    unsigned long int local_row = alpha - 1UL;
    for (p = plan->row_ptr[local_row]; p < plan->row_ptr[local_row + 1UL]; p++) {
      size_t index = (size_t)(alpha - 1UL) * (size_t)plan->dim +
                     (size_t)(plan->col_index[p] - 1UL);
      dense[index] += plan->values[p];
      multiplicity[index]++;
      if (multiplicity[index] > 1U) duplicate_found = 1;
    }
  }
  if (require_duplicate != 0) assert_int_eq(duplicate_found, 1, label);

  for (alpha = 1UL; alpha <= plan->dim; alpha++) {
    for (beta = 1UL; beta <= plan->dim; beta++) {
      double complex plan_value = dense[(size_t)(alpha - 1UL) * (size_t)plan->dim +
                                        (size_t)(beta - 1UL)];
      double complex expected = canonicalized_matrix_element(X, alpha, beta);
      double complex transpose = dense[(size_t)(beta - 1UL) * (size_t)plan->dim +
                                       (size_t)(alpha - 1UL)];
      assert_complex_close(plan_value, expected, 1.0e-12, label);
      assert_complex_close(plan_value, conj(transpose), 1.0e-12, label);
    }
  }

  for (beta = 1UL; beta <= plan->dim; beta++) {
    input[beta] = 0.125 * (double)beta + I * 0.0625 * (double)(beta + 1UL);
  }
  {
    struct LegacyVectorContext legacy;
    legacy.output = legacy_output;
    legacy.dim = plan->dim;
    for (beta = 1UL; beta <= plan->dim; beta++) {
      legacy.input_amp = input[beta];
      if (SymmetryEnumerateColumn(X, beta, accumulate_legacy_vector_entry,
                                  &legacy) != 0) {
        fprintf(stderr, "%s: legacy vector scan failed\n", label);
        exit(1);
      }
    }
  }
  if (ApplySymmetryMatvecPlan(X, output, input, &plan_prdct) != 0) {
    fprintf(stderr, "%s: plan apply failed\n", label);
    exit(1);
  }
  for (alpha = 1UL; alpha <= plan->dim; alpha++) {
    double complex expected = 0.0;
    for (beta = 1UL; beta <= plan->dim; beta++) {
      expected += dense[(size_t)(alpha - 1UL) * (size_t)plan->dim +
                        (size_t)(beta - 1UL)] * input[beta];
    }
    assert_complex_close(output[alpha], expected, 1.0e-12, label);
    assert_complex_close(output[alpha], legacy_output[alpha], 1.0e-12, label);
    difference_norm2 += pow(cabs(output[alpha] - legacy_output[alpha]), 2.0);
    legacy_norm2 += pow(cabs(legacy_output[alpha]), 2.0);
    expected_prdct += conj(input[alpha]) * expected;
  }
  if (sqrt(difference_norm2) / fmax(sqrt(legacy_norm2), 1.0e-300) > 1.0e-12) {
    fprintf(stderr, "%s: plan/legacy relative L2 error exceeds 1e-12\n", label);
    exit(1);
  }
  assert_complex_close(plan_prdct, expected_prdct, 1.0e-12, label);

  X->Sym->local_offset++;
  assert_int_eq(ApplySymmetryMatvecPlan(X, output, input, &plan_prdct), -1,
                "plan snapshot guard");
  X->Sym->local_offset--;

  free(dense);
  free(multiplicity);
  free(input);
  free(legacy_output);
  free(output);
}

static void assert_spin_plan(unsigned int nsite,
                             unsigned int nup,
                             unsigned int momentum_index,
                             double diagonal_coupling,
                             const char *label)
{
  struct BindStruct X;
  setup_bind(&X, nsite, nup, momentum_index);
  if (diagonal_coupling != 0.0) set_ising_ring_diagonal(nsite, diagonal_coupling);
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: BuildSymmetryBasis failed\n", label);
    exit(1);
  }
  assert_plan_matches_canonicalized_matrix(&X, 1, label);
  FreeSymmetryBasis(X.Sym);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

static void assert_spin_diagonal_only_plan(const char *label)
{
  struct BindStruct X;
  setup_bind(&X, 6, 3, 1);
  X.Def.NExchangeCoupling = 0U;
  set_ising_ring_diagonal(6, 0.37);
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: BuildSymmetryBasis failed\n", label);
    exit(1);
  }
  assert_plan_matches_canonicalized_matrix(&X, 0, label);
  FreeSymmetryBasis(X.Sym);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

static void assert_spinless_plan(unsigned int nsite,
                                 unsigned int ne,
                                 unsigned int momentum_index,
                                 double density_coupling,
                                 const char *label)
{
  struct BindStruct X;
  setup_spinless_bind(&X, nsite, ne, momentum_index);
  setup_spinless_transfer_ring(&X.Def, nsite);
  if (density_coupling != 0.0) {
    setup_spinless_coulomb_ring(&X.Def, nsite, density_coupling);
    set_spinless_coulomb_ring_diagonal(nsite, density_coupling);
  }
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: BuildSymmetryBasis failed\n", label);
    exit(1);
  }
  assert_plan_matches_canonicalized_matrix(&X, 0, label);
  FreeSymmetryBasis(X.Sym);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

static void assert_hubbard_plan(unsigned int nsite,
                                unsigned int nup,
                                unsigned int ndown,
                                unsigned int momentum_index,
                                double coulomb_intra,
                                const char *label)
{
  struct BindStruct X;
  setup_hubbard_bind(&X, nsite, nup, ndown, momentum_index);
  setup_hubbard_transfer_ring(&X.Def, nsite);
  if (coulomb_intra != 0.0) {
    setup_hubbard_coulomb_intra(&X.Def, nsite, coulomb_intra);
  }
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: BuildSymmetryBasis failed\n", label);
    exit(1);
  }
  assert_plan_matches_canonicalized_matrix(&X, 0, label);
  FreeSymmetryBasis(X.Sym);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

static void assert_zero_row_plan(const char *label)
{
  struct BindStruct X;
  double complex input[2] = {0.0, 1.0};
  double complex output[1] = {0.0};
  double complex prdct = 1.0;
  setup_bind(&X, 4, 2, 1);
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: BuildSymmetryBasis failed\n", label);
    exit(1);
  }
  nproc = 2;
  myrank = 1;
  if (ActivateSymmetryBasisDimension(&X) != 0 || BuildSymmetryMatvecPlan(&X) != 0) {
    fprintf(stderr, "%s: zero-row plan setup failed\n", label);
    exit(1);
  }
  assert_ulong_eq(X.Sym->local_dim, 0UL, label);
  assert_int_eq(X.Sym->matvec_plan != NULL, 1, label);
  assert_ulong_eq((unsigned long int)X.Sym->matvec_plan->nnz, 0UL, label);
  assert_int_eq(ApplySymmetryMatvecPlan(&X, output, input, &prdct), 0, label);
  assert_complex_close(prdct, 0.0, 1.0e-12, label);
  nproc = 1;
  myrank = 0;
  FreeSymmetryBasis(X.Sym);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
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

static void assert_spinless_symmetry_dim(unsigned int nsite,
                                         unsigned int ne,
                                         unsigned int momentum_index,
                                         unsigned long int expected_dim,
                                         const char *label)
{
  struct BindStruct X;
  setup_spinless_bind(&X, nsite, ne, momentum_index);
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

static void assert_hubbard_symmetry_dim(unsigned int nsite,
                                        unsigned int nup,
                                        unsigned int ndown,
                                        unsigned int momentum_index,
                                        unsigned long int expected_dim,
                                        const char *label)
{
  struct BindStruct X;
  setup_hubbard_bind(&X, nsite, nup, ndown, momentum_index);
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

static void assert_spinless_canonicalized_matrix_matches_raw(unsigned int nsite,
                                                             unsigned int ne,
                                                             unsigned int momentum_index,
                                                             double density_coupling,
                                                             const char *label)
{
  struct BindStruct X;
  unsigned long int alpha, beta;
  setup_spinless_bind(&X, nsite, ne, momentum_index);
  setup_spinless_transfer_ring(&X.Def, nsite);
  if (density_coupling != 0.0) {
    setup_spinless_coulomb_ring(&X.Def, nsite, density_coupling);
    set_spinless_coulomb_ring_diagonal(nsite, density_coupling);
  }
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

static void assert_hubbard_canonicalized_matrix_matches_raw(unsigned int nsite,
                                                            unsigned int nup,
                                                            unsigned int ndown,
                                                            unsigned int momentum_index,
                                                            double coulomb_intra,
                                                            const char *label)
{
  struct BindStruct X;
  unsigned long int alpha, beta;
  setup_hubbard_bind(&X, nsite, nup, ndown, momentum_index);
  setup_hubbard_transfer_ring(&X.Def, nsite);
  if (coulomb_intra != 0.0) {
    setup_hubbard_coulomb_intra(&X.Def, nsite, coulomb_intra);
  }
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

static void assert_representative_hash_matches_basis(unsigned int nsite,
                                                     unsigned int nup,
                                                     unsigned int momentum_index,
                                                     const char *label)
{
  struct BindStruct X;
  unsigned long int beta;
  unsigned long int slot;
  unsigned long int occupied = 0;
  struct SymmetryCanonicalResult result;
  setup_bind(&X, nsite, nup, momentum_index);
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: BuildSymmetryBasis failed\n", label);
    exit(1);
  }
  assert_int_eq(X.Sym->capacity >= X.Sym->dim, 1, label);
  assert_int_eq(X.Sym->rep_hash_size > X.Sym->dim, 1, label);
  assert_int_eq(X.Sym->rep_hash_keys != NULL, 1, label);
  assert_int_eq(X.Sym->rep_hash_values != NULL, 1, label);

  for (slot = 0; slot < X.Sym->rep_hash_size; slot++) {
    unsigned long int basis_index = X.Sym->rep_hash_values[slot];
    if (basis_index == 0UL) continue;
    occupied++;
    assert_int_eq(basis_index <= X.Sym->dim, 1, label);
    assert_ulong_eq(X.Sym->basis[basis_index].rep_state,
                    X.Sym->rep_hash_keys[slot],
                    label);
  }
  assert_ulong_eq(occupied, X.Sym->dim, label);

  for (beta = 1; beta <= X.Sym->dim; beta++) {
    assert_int_eq(SymmetryCanonicalizeSpinState(&X,
                                                X.Sym->basis[beta].rep_state,
                                                &result),
                  0,
                  label);
    assert_int_eq(result.found, 1, label);
    assert_ulong_eq(result.basis_index, beta, label);
  }

  assert_int_eq(SymmetryCanonicalizeSpinState(&X, 0UL, &result), 0, label);
  assert_int_eq(result.found, 0, label);

  FreeSymmetryBasis(X.Sym);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

static void assert_hash_probe_lookup_handles_collision(const char *label)
{
  struct BindStruct X;
  struct SymmetryBasisRuntime sym;
  struct SymmetryBasisVector basis[2];
  unsigned long int keys[2];
  unsigned long int values[2];
  int identity_perm[4] = {0, 1, 2, 3};
  int identity_anti[4] = {1, 1, 1, 1};
  int *perm_rows_one[1];
  int *anti_rows_one[1];
  double complex chars_one[1];
  unsigned long int target_state = 0x5UL;
  unsigned long int dummy_state = 0x6UL;
  int scenario;

  memset(&X, 0, sizeof(X));
  memset(&sym, 0, sizeof(sym));
  memset(basis, 0, sizeof(basis));
  perm_rows_one[0] = identity_perm;
  anti_rows_one[0] = identity_anti;
  chars_one[0] = 1.0;
  X.Def.iFlgSymmetryBasis = TRUE;
  X.Def.iCalcModel = Spin;
  X.Def.Nsite = 4;
  X.Def.NSymTrans = 1;
  X.Def.SymTrans = perm_rows_one;
  X.Def.SymTransAnti = anti_rows_one;
  X.Def.SymTransChar = chars_one;
  X.Sym = &sym;
  sym.enabled = TRUE;
  sym.dim = 1;
  sym.basis = basis;
  sym.rep_hash_size = 2;
  sym.rep_hash_keys = keys;
  sym.rep_hash_values = values;
  basis[1].rep_state = target_state;

  for (scenario = 0; scenario < 2; scenario++) {
    struct SymmetryCanonicalResult result;
    if (scenario == 0) {
      keys[0] = dummy_state;
      values[0] = 2UL;
      keys[1] = target_state;
      values[1] = 1UL;
    } else {
      keys[0] = target_state;
      values[0] = 1UL;
      keys[1] = dummy_state;
      values[1] = 2UL;
    }
    assert_int_eq(SymmetryCanonicalizeSpinState(&X, target_state, &result), 0, label);
    assert_int_eq(result.found, 1, label);
    assert_ulong_eq(result.basis_index, 1UL, label);
  }
}

int main(void)
{
  int shift4[4] = {1, 2, 3, 0};
  stdoutMPI = stderr;
  assert_ulong_eq(SymmetryApplyToSpinBits(0x1UL, shift4, 4), 0x2UL, "single bit shift");
  assert_ulong_eq(SymmetryApplyToSpinBits(0x9UL, shift4, 4), 0x3UL, "wrap shift");
  assert_ulong_eq(SymmetryApplyToSpinBits(0x6UL, shift4, 4), 0xcUL, "two bit shift");
  {
    struct DefineList def;
    struct SymmetryTransformResult moved;
    setup_cyclic_def(&def, 4, 0);
    def.iCalcModel = SpinlessFermion;
    assert_int_eq(SymmetryApplyToState(&def, 0x3UL, 1, &moved), 0,
                  "spinless adjacent translate");
    assert_ulong_eq(moved.state, 0x6UL, "spinless adjacent translated state");
    assert_complex_close(moved.amplitude, 1.0, 1.0e-12,
                         "spinless adjacent translation sign");
    assert_int_eq(SymmetryApplyToState(&def, 0x9UL, 1, &moved), 0,
                  "spinless wrap translate");
    assert_ulong_eq(moved.state, 0x3UL, "spinless wrap translated state");
    assert_complex_close(moved.amplitude, -1.0, 1.0e-12,
                         "spinless wrap translation sign");
  }
  {
    struct DefineList def;
    struct SymmetryTransformResult moved;
    setup_cyclic_def(&def, 4, 0);
    def.iCalcModel = Hubbard;
    assert_int_eq(SymmetryApplyToState(&def, (1UL << 0) | (1UL << 6), 1, &moved), 0,
                  "hubbard same-spin wrap translate");
    assert_ulong_eq(moved.state, (1UL << 0) | (1UL << 2),
                    "hubbard same-spin wrap translated state");
    assert_complex_close(moved.amplitude, -1.0, 1.0e-12,
                         "hubbard same-spin wrap translation sign");
    assert_int_eq(SymmetryApplyToState(&def, (1UL << 0) | (1UL << 7), 1, &moved), 0,
                  "hubbard mixed-spin wrap translate");
    assert_ulong_eq(moved.state, (1UL << 1) | (1UL << 2),
                    "hubbard mixed-spin wrap translated state");
    assert_complex_close(moved.amplitude, -1.0, 1.0e-12,
                         "hubbard mixed-spin wrap translation sign");
  }
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
  assert_spinless_symmetry_dim(4, 2, 0, 1,
                               "SpinlessFermion C4 k=0 sector dimension");
  assert_spinless_symmetry_dim(4, 2, 1, 2,
                               "SpinlessFermion C4 k=pi/2 sector dimension");
  assert_hubbard_symmetry_dim(4, 1, 1, 0, 4,
                              "Hubbard C4 k=0 sector dimension");
  assert_hubbard_symmetry_dim(4, 1, 1, 1, 4,
                              "Hubbard C4 k=pi/2 sector dimension");
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
  assert_spinless_canonicalized_matrix_matches_raw(4, 2, 1, 0.0,
                                                   "SpinlessFermion C4 k=pi/2 matrix matches raw reference");
  assert_spinless_canonicalized_matrix_matches_raw(4, 2, 1, 0.25,
                                                   "SpinlessFermion C4 k=pi/2 CoulombInter matrix matches raw reference");
  assert_hubbard_canonicalized_matrix_matches_raw(4, 1, 1, 0, 0.0,
                                                  "Hubbard C4 k=0 matrix matches raw reference");
  assert_hubbard_canonicalized_matrix_matches_raw(4, 1, 1, 1, 0.5,
                                                  "Hubbard C4 k=pi/2 CoulombIntra matrix matches raw reference");
  assert_hubbard_canonicalized_matrix_matches_raw(4, 2, 1, 1, 0.0,
                                                  "Hubbard C4 k=pi/2 same-spin matrix matches raw reference");
  assert_orbit_diagonal_is_representative(6, 3, 1, 1.0,
                                          "C6 k=pi/3 Ising diagonal is orbit-invariant");
  assert_canonicalized_matrix_matches_raw(6, 3, 1, 1.0,
                                          "C6 k=pi/3 Ising canonicalized matrix matches raw reference");
  assert_spin_plan(6, 3, 0, 1.0,
                   "C6 k=0 Spin local-row plan matches canonicalized matrix");
  assert_spin_plan(6, 3, 1, 1.0,
                   "C6 k=pi/3 Spin local-row plan matches canonicalized matrix");
  assert_spin_diagonal_only_plan(
      "C6 k=pi/3 diagonal-only local-row plan matches canonicalized matrix");
  assert_spinless_plan(6, 3, 0, 0.25,
                       "SpinlessFermion C6 k=0 local-row plan matches canonicalized matrix");
  assert_spinless_plan(6, 3, 1, 0.25,
                       "SpinlessFermion C6 k=pi/3 local-row plan matches canonicalized matrix");
  assert_hubbard_plan(4, 2, 2, 0, 0.5,
                      "Hubbard C4 k=0 local-row plan matches canonicalized matrix");
  assert_hubbard_plan(4, 2, 2, 1, 0.5,
                      "Hubbard C4 k=pi/2 local-row plan matches canonicalized matrix");
  assert_zero_row_plan("local-row plan supports zero-row rank");
  assert_representative_hash_matches_basis(6, 3, 1,
                                           "C6 k=pi/3 representative hash matches basis");
  assert_hash_probe_lookup_handles_collision("representative hash probing handles collisions");
  return 0;
}
