#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "DefCommon.h"
#include "symmetry_basis.h"
#include "symmetry_diagonal.h"
#include "symmetry_distribution.h"
#include "symmetry_matvec_plan.h"
#include "symmetry_state_enumerator.h"
#include "symmetry_vector_halo.h"
#include "struct.h"

#ifdef MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

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

int SumMPI_i(int value)
{
  return value;
}

unsigned long int SumMPI_li(unsigned long int value)
{
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
static int hund_storage[12][2];
static int *hund_rows[12];
static double hund_params[12];

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

static void set_ising_ring_diagonal(struct DefineList *def,
                                    unsigned int nsite,
                                    double coupling)
{
  unsigned long int raw;
  unsigned int site;
  def->NIsingCoupling = nsite;
  def->NCoulombInter = nsite;
  def->CoulombInter = coulomb_rows;
  def->ParaCoulombInter = coulomb_params;
  def->NHundCoupling = nsite;
  def->HundCoupling = hund_rows;
  def->ParaHundCoupling = hund_params;
  for (site = 0U; site < nsite; site++) {
    unsigned int next = (site + 1U) % nsite;
    coulomb_rows[site] = coulomb_storage[site];
    coulomb_storage[site][0] = (int)site;
    coulomb_storage[site][1] = (int)next;
    coulomb_params[site] = -coupling / 4.0;
    hund_rows[site] = hund_storage[site];
    hund_storage[site][0] = (int)site;
    hund_storage[site][1] = (int)next;
    hund_params[site] = -coupling / 2.0;
  }
  for (raw = 1; raw <= test_raw_dim; raw++) {
    double diagonal = 0.0;
    for (site = 0U; site < nsite; site++) {
      diagonal += coulomb_params[site];
    }
    for (site = 0U; site < nsite; site++) {
      unsigned int next = (site + 1U) % nsite;
      if (((list_1[raw] >> site) & 1UL) ==
          ((list_1[raw] >> next) & 1UL)) {
        diagonal += -hund_params[site];
      }
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
  X->Def.Nup = nup;
  X->Def.Ndown = nsite - nup;
  X->Def.Ne = nup;
  X->Def.iFlgSzConserved = TRUE;
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
  if (alpha == beta) {
    value += X->Sym->basis[beta].diagonal;
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

static int state_matches_sector(int model,
                                unsigned long int state,
                                unsigned int nsite,
                                unsigned int nup,
                                unsigned int ndown)
{
  if (model == Hubbard) {
    return (unsigned int)count_hubbard_spin(state, nsite, 0U) == nup &&
           (unsigned int)count_hubbard_spin(state, nsite, 1U) == ndown;
  }
  return (unsigned int)popcount_ulong(state) == nup;
}

static void assert_state_enumerator_sequence(int model,
                                             unsigned int nsite,
                                             unsigned int nup,
                                             unsigned int ndown,
                                             int thread_count,
                                             const char *label)
{
  struct DefineList def;
  struct SymmetryStateEnumerator enumerator;
  unsigned long int limit = 1UL << (model == Hubbard ? 2U * nsite : nsite);
  unsigned long int state, raw, dim = 0UL;
  unsigned long int *states;
  int failed = 0;
  memset(&def, 0, sizeof(def));
  def.iCalcModel = model;
  def.Nsite = nsite;
  def.Nup = nup;
  def.Ndown = ndown;
  def.Ne = model == SpinlessFermion ? nup : nup + ndown;
  if (model == Spin) def.iFlgSzConserved = TRUE;
  for (state = 0UL; state < limit; state++) {
    if (state_matches_sector(model, state, nsite, nup, ndown)) dim++;
  }
  assert_int_eq(InitSymmetryStateEnumerator(&def, dim, &enumerator), 0, label);
  assert_ulong_eq(enumerator.raw_dim, dim, label);
  states = (unsigned long int *)calloc(dim + 1UL, sizeof(*states));
  if (states == NULL) {
    fprintf(stderr, "%s: allocation failed\n", label);
    exit(1);
  }
#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(thread_count);
#pragma omp parallel for schedule(static) reduction(|:failed)
#else
  (void)thread_count;
#endif
  for (raw = 1UL; raw <= dim; raw++) {
    if (SymmetryStateEnumeratorStateAt(&enumerator, raw, &states[raw]) != 0) {
      failed = 1;
    }
  }
  assert_int_eq(failed, 0, label);
  raw = 0UL;
  for (state = 0UL; state < limit; state++) {
    if (!state_matches_sector(model, state, nsite, nup, ndown)) continue;
    raw++;
    assert_ulong_eq(states[raw], state, label);
    if (raw > 1UL) assert_int_eq(states[raw] > states[raw - 1UL], 1, label);
  }
  assert_ulong_eq(raw, dim, label);
  free(states);
}

static void assert_state_enumerator_exact(void)
{
  const unsigned int word_bits =
      (unsigned int)(CHAR_BIT * sizeof(unsigned long int));
  struct DefineList def;
  struct SymmetryStateEnumerator enumerator;
  unsigned int nsite, nup, ndown;
  unsigned long int state = 123UL;
  int thread_count;
#ifdef _OPENMP
  int saved_dynamic = omp_get_dynamic();
  int saved_threads = omp_get_max_threads();
#endif
  for (thread_count = 1; thread_count <= 4; thread_count += 3) {
    for (nsite = 1U; nsite <= 8U; nsite++) {
      for (nup = 0U; nup <= nsite; nup++) {
        assert_state_enumerator_sequence(
            Spin, nsite, nup, nsite - nup, thread_count,
            "Spin raw-state enumerator matches numeric list_1 order");
        assert_state_enumerator_sequence(
            SpinlessFermion, nsite, nup, 0U, thread_count,
            "Spinless raw-state enumerator matches numeric list_1 order");
      }
    }
    for (nsite = 1U; nsite <= 5U; nsite++) {
      for (nup = 0U; nup <= nsite; nup++) {
        for (ndown = 0U; ndown <= nsite; ndown++) {
          assert_state_enumerator_sequence(
              Hubbard, nsite, nup, ndown, thread_count,
              "Hubbard raw-state enumerator matches numeric list_1 order");
        }
      }
    }
  }

  memset(&def, 0, sizeof(def));
  def.iCalcModel = Spin;
  def.Nsite = 4U;
  def.Nup = 2U;
  def.Ndown = 2U;
  def.iFlgSzConserved = TRUE;
  assert_int_eq(InitSymmetryStateEnumerator(&def, 6UL, &enumerator), 0,
                "valid enumerator initializes");
  assert_int_eq(SymmetryStateEnumeratorStateAt(&enumerator, 0UL, &state), -1,
                "raw index zero rejects");
  assert_int_eq(SymmetryStateEnumeratorStateAt(&enumerator, 7UL, &state), -1,
                "raw index above dimension rejects");
  assert_int_eq(SymmetryStateEnumeratorStateAt(NULL, 1UL, &state), -1,
                "null enumerator rejects");
  assert_int_eq(SymmetryStateEnumeratorStateAt(&enumerator, 1UL, NULL), -1,
                "null state output rejects");
  assert_ulong_eq(state, 123UL, "failed enumeration preserves output");
  enumerator.bit_count++;
  assert_int_eq(SymmetryStateEnumeratorStateAt(&enumerator, 1UL, &state), -1,
                "inconsistent enumerator metadata rejects");
  assert_ulong_eq(state, 123UL, "invalid enumerator preserves output");
  assert_int_eq(InitSymmetryStateEnumerator(&def, 5UL, &enumerator), -1,
                "dimension mismatch rejects");
  def.iFlgGeneralSpin = TRUE;
  assert_int_eq(InitSymmetryStateEnumerator(&def, 6UL, &enumerator), -1,
                "general Spin rejects");
  def.iFlgGeneralSpin = FALSE;
  def.iCalcModel = SpinGC;
  assert_int_eq(InitSymmetryStateEnumerator(&def, 16UL, &enumerator), -1,
                "unsupported model rejects");
  def.iCalcModel = Spin;
  def.Nsite = word_bits + 1U;
  def.Nup = 0U;
  def.Ndown = def.Nsite;
  assert_int_eq(InitSymmetryStateEnumerator(&def, 1UL, &enumerator), -1,
                "Spin state wider than unsigned long rejects");
  def.iCalcModel = Hubbard;
  def.Nsite = word_bits / 2U + 1U;
  def.Nup = 0U;
  def.Ndown = 0U;
  def.Ne = 0U;
  assert_int_eq(InitSymmetryStateEnumerator(&def, 1UL, &enumerator), -1,
                "Hubbard state wider than unsigned long rejects");

  memset(&def, 0, sizeof(def));
  def.iCalcModel = Spin;
  def.Nsite = word_bits;
  def.Nup = word_bits;
  def.Ndown = 0U;
  def.iFlgSzConserved = TRUE;
  assert_int_eq(InitSymmetryStateEnumerator(&def, 1UL, &enumerator), 0,
                "full Spin sector at word width initializes");
  assert_int_eq(SymmetryStateEnumeratorStateAt(&enumerator, 1UL, &state), 0,
                "full Spin sector at word width enumerates");
  assert_ulong_eq(state, ULONG_MAX, "full Spin word-width state is exact");
  def.iCalcModel = Hubbard;
  def.Nsite = word_bits / 2U;
  def.Nup = def.Nsite;
  def.Ndown = def.Nsite;
  def.Ne = def.Nup + def.Ndown;
  assert_int_eq(InitSymmetryStateEnumerator(&def, 1UL, &enumerator), 0,
                "full Hubbard sector at word width initializes");
  assert_int_eq(SymmetryStateEnumeratorStateAt(&enumerator, 1UL, &state), 0,
                "full Hubbard sector at word width enumerates");
  assert_ulong_eq(state, ULONG_MAX, "full Hubbard word-width state is exact");
  assert_int_eq(InitSymmetryStateEnumerator(NULL, 1UL, &enumerator), -1,
                "null definition rejects");
  assert_int_eq(InitSymmetryStateEnumerator(&def, 1UL, NULL), -1,
                "null enumerator output rejects");
#ifdef _OPENMP
  omp_set_num_threads(saved_threads);
  omp_set_dynamic(saved_dynamic);
#endif
}

static void assert_double_bitwise(double got,
                                  double expected,
                                  const char *label)
{
  uint64_t got_bits;
  uint64_t expected_bits;
  memcpy(&got_bits, &got, sizeof(got_bits));
  memcpy(&expected_bits, &expected, sizeof(expected_bits));
  if (got_bits != expected_bits) {
    fprintf(stderr, "%s: got %016llx expected %016llx\n",
            label,
            (unsigned long long)got_bits,
            (unsigned long long)expected_bits);
    exit(1);
  }
}

static double reference_symmetry_diagonal(const struct DefineList *def,
                                          unsigned long int state)
{
  unsigned int index;
  double value = 0.0;
  if (def->iCalcModel == Spin) {
    for (index = 0U; index < def->NCoulombInter; index++) {
      value += def->ParaCoulombInter[index];
    }
    for (index = 0U; index < def->NHundCoupling; index++) {
      unsigned int site0 = (unsigned int)def->HundCoupling[index][0];
      unsigned int site1 = (unsigned int)def->HundCoupling[index][1];
      if (((state >> site0) & 1UL) == ((state >> site1) & 1UL)) {
        value += -def->ParaHundCoupling[index];
      }
    }
  } else if (def->iCalcModel == SpinlessFermion) {
    for (index = 0U; index < def->NCoulombInter; index++) {
      unsigned int site0 = (unsigned int)def->CoulombInter[index][0];
      unsigned int site1 = (unsigned int)def->CoulombInter[index][1];
      if (((state >> site0) & 1UL) != 0UL &&
          ((state >> site1) & 1UL) != 0UL) {
        value += def->ParaCoulombInter[index];
      }
    }
  } else if (def->iCalcModel == Hubbard) {
    for (index = 0U; index < def->NCoulombIntra; index++) {
      unsigned int site = (unsigned int)def->CoulombIntra[index][0];
      if (((state >> (2U * site)) & 1UL) != 0UL &&
          ((state >> (2U * site + 1U)) & 1UL) != 0UL) {
        value += def->ParaCoulombIntra[index];
      }
    }
  }
  return value;
}

static void assert_state_diagonal_sequence(const struct DefineList *def,
                                           int thread_count,
                                           const char *label)
{
  unsigned int bit_count =
      def->iCalcModel == Hubbard ? 2U * def->Nsite : def->Nsite;
  unsigned long int limit = 1UL << bit_count;
  unsigned long int state;
  double *values = (double *)calloc((size_t)limit, sizeof(*values));
  int failed = 0;
  if (values == NULL) {
    fprintf(stderr, "%s: allocation failed\n", label);
    exit(1);
  }
#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(thread_count);
#pragma omp parallel for schedule(static) reduction(|:failed)
#else
  (void)thread_count;
#endif
  for (state = 0UL; state < limit; state++) {
    if (EvaluateSymmetryStateDiagonal(def, state, &values[state]) != 0) {
      failed = 1;
    }
  }
  assert_int_eq(failed, 0, label);
  for (state = 0UL; state < limit; state++) {
    assert_double_bitwise(
        values[state], reference_symmetry_diagonal(def, state), label);
  }
  free(values);
}

static void assert_state_diagonal_exact(void)
{
  const unsigned int word_bits =
      (unsigned int)(CHAR_BIT * sizeof(unsigned long int));
  struct DefineList spin_def;
  struct DefineList spinless_def;
  struct DefineList hubbard_def;
  struct DefineList invalid;
  int spin_coulomb_storage[6][2];
  int spin_hund_storage[6][2];
  int *spin_coulomb_rows[6];
  int *spin_hund_rows[6];
  double spin_coulomb_parameters[6];
  double spin_hund_parameters[6];
  int spinless_storage[6][2];
  int *spinless_rows[6];
  double spinless_parameters[6] = {
      0.125, -0.375, 0.2, 0.0, 1.125, -0.0625};
  int hubbard_storage[6][1];
  int *hubbard_rows[6];
  double hubbard_parameters[6] = {
      0.25, -0.5, 1.125, 0.0625, -0.125, 0.75};
  const int spin_pairs[6][2] = {
      {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}, {0, 2}};
  const int spinless_pairs[6][2] = {
      {0, 1}, {2, 1}, {2, 3}, {4, 3}, {4, 0}, {1, 4}};
  const double ising_parameters[6] = {
      0.3, -0.7, 0.125, 1.1, -0.2, 0.05};
  const int hubbard_sites[6] = {0, 1, 2, 3, 1, 3};
  unsigned int index;
  double diagonal = 19.25;
  int thread_count;
#ifdef _OPENMP
  int saved_dynamic = omp_get_dynamic();
  int saved_threads = omp_get_max_threads();
#endif

  memset(&spin_def, 0, sizeof(spin_def));
  spin_def.iCalcModel = Spin;
  spin_def.Nsite = 5U;
  spin_def.iFlgSzConserved = TRUE;
  spin_def.NIsingCoupling = 6U;
  spin_def.NCoulombInter = 6U;
  spin_def.NHundCoupling = 6U;
  spin_def.CoulombInter = spin_coulomb_rows;
  spin_def.ParaCoulombInter = spin_coulomb_parameters;
  spin_def.HundCoupling = spin_hund_rows;
  spin_def.ParaHundCoupling = spin_hund_parameters;
  for (index = 0U; index < 6U; index++) {
    spin_coulomb_rows[index] = spin_coulomb_storage[index];
    spin_hund_rows[index] = spin_hund_storage[index];
    spin_coulomb_storage[index][0] = spin_pairs[index][0];
    spin_coulomb_storage[index][1] = spin_pairs[index][1];
    spin_hund_storage[index][0] = spin_pairs[index][0];
    spin_hund_storage[index][1] = spin_pairs[index][1];
    spin_coulomb_parameters[index] = -ising_parameters[index] / 4.0;
    spin_hund_parameters[index] = -ising_parameters[index] / 2.0;
  }

  memset(&spinless_def, 0, sizeof(spinless_def));
  spinless_def.iCalcModel = SpinlessFermion;
  spinless_def.Nsite = 5U;
  spinless_def.NCoulombInter = 6U;
  spinless_def.CoulombInter = spinless_rows;
  spinless_def.ParaCoulombInter = spinless_parameters;
  for (index = 0U; index < 6U; index++) {
    spinless_rows[index] = spinless_storage[index];
    spinless_storage[index][0] = spinless_pairs[index][0];
    spinless_storage[index][1] = spinless_pairs[index][1];
  }

  memset(&hubbard_def, 0, sizeof(hubbard_def));
  hubbard_def.iCalcModel = Hubbard;
  hubbard_def.Nsite = 4U;
  hubbard_def.NCoulombIntra = 6U;
  hubbard_def.CoulombIntra = hubbard_rows;
  hubbard_def.ParaCoulombIntra = hubbard_parameters;
  for (index = 0U; index < 6U; index++) {
    hubbard_rows[index] = hubbard_storage[index];
    hubbard_storage[index][0] = hubbard_sites[index];
  }

  for (thread_count = 1; thread_count <= 4; thread_count += 3) {
    assert_state_diagonal_sequence(
        &spin_def, thread_count,
        "Spin single-state diagonal matches diagonalcalc semantics");
    assert_state_diagonal_sequence(
        &spinless_def, thread_count,
        "Spinless single-state diagonal matches diagonalcalc semantics");
    assert_state_diagonal_sequence(
        &hubbard_def, thread_count,
        "Hubbard single-state diagonal matches diagonalcalc semantics");
  }

  invalid = spinless_def;
  assert_int_eq(EvaluateSymmetryStateDiagonal(NULL, 0UL, &diagonal), -1,
                "null diagonal definition rejects");
  assert_int_eq(EvaluateSymmetryStateDiagonal(&invalid, 0UL, NULL), -1,
                "null diagonal output rejects");
  assert_int_eq(
      EvaluateSymmetryStateDiagonal(
          &invalid, 1UL << invalid.Nsite, &diagonal),
      -1, "state bits outside model width reject");
  invalid.iCalcModel = SpinGC;
  assert_int_eq(EvaluateSymmetryStateDiagonal(&invalid, 0UL, &diagonal), -1,
                "unsupported diagonal model rejects");
  invalid = spinless_def;
  invalid.CoulombInter = NULL;
  assert_int_eq(EvaluateSymmetryStateDiagonal(&invalid, 0UL, &diagonal), -1,
                "missing diagonal term storage rejects");
  invalid = spinless_def;
  invalid.EDNChemi = 1U;
  assert_int_eq(EvaluateSymmetryStateDiagonal(&invalid, 0UL, &diagonal), -1,
                "unsupported diagonal term family rejects");
  invalid = spin_def;
  invalid.iFlgGeneralSpin = TRUE;
  assert_int_eq(EvaluateSymmetryStateDiagonal(&invalid, 0UL, &diagonal), -1,
                "general Spin diagonal rejects");
  invalid = spin_def;
  invalid.NHundCoupling--;
  assert_int_eq(EvaluateSymmetryStateDiagonal(&invalid, 0UL, &diagonal), -1,
                "incomplete Ising expansion rejects");
  invalid = hubbard_def;
  hubbard_storage[0][0] = (int)hubbard_def.Nsite;
  assert_int_eq(EvaluateSymmetryStateDiagonal(&invalid, 0UL, &diagonal), -1,
                "diagonal term site outside model rejects");
  hubbard_storage[0][0] = hubbard_sites[0];
  assert_double_bitwise(diagonal, 19.25,
                        "failed diagonal evaluation preserves output");

  memset(&invalid, 0, sizeof(invalid));
  invalid.iCalcModel = Spin;
  invalid.Nsite = word_bits;
  assert_int_eq(
      EvaluateSymmetryStateDiagonal(&invalid, ULONG_MAX, &diagonal), 0,
      "word-width Spin state evaluates");
  assert_double_bitwise(diagonal, 0.0,
                        "empty word-width Spin diagonal is zero");
  invalid.iCalcModel = Hubbard;
  invalid.Nsite = word_bits / 2U;
  assert_int_eq(
      EvaluateSymmetryStateDiagonal(&invalid, ULONG_MAX, &diagonal), 0,
      "word-width Hubbard state evaluates");
  assert_double_bitwise(diagonal, 0.0,
                        "empty word-width Hubbard diagonal is zero");
#ifdef _OPENMP
  omp_set_num_threads(saved_threads);
  omp_set_dynamic(saved_dynamic);
#endif
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

static void assert_complex_bitwise(double complex got,
                                   double complex expected,
                                   const char *label)
{
  uint64_t got_bits[2];
  uint64_t expected_bits[2];
  memcpy(got_bits, &got, sizeof(got_bits));
  memcpy(expected_bits, &expected, sizeof(expected_bits));
  if (memcmp(got_bits, expected_bits, sizeof(got_bits)) != 0) {
    fprintf(stderr,
            "%s: got %016llx %016llx expected %016llx %016llx\n",
            label,
            (unsigned long long)got_bits[0],
            (unsigned long long)got_bits[1],
            (unsigned long long)expected_bits[0],
            (unsigned long long)expected_bits[1]);
    exit(1);
  }
}

static uint64_t symmetry_plan_column_slot(
    const struct SymmetryMatvecPlan *plan, size_t column)
{
  if (plan->column_slot_width == SYMMETRY_COLUMN_U32 &&
      plan->column_slot32 != NULL) {
    return (uint64_t)plan->column_slot32[column];
  }
  if (plan->column_slot_width == SYMMETRY_COLUMN_U64 &&
      plan->column_slot64 != NULL) {
    return plan->column_slot64[column];
  }
  fprintf(stderr, "column slot storage is not initialized\n");
  exit(1);
}

static int apply_direct_flat_global(
    const struct SymmetryMatvecPlan *plan,
    double complex *output,
    const double complex *input,
    double complex *prdct_out)
{
  unsigned long int local_row;
  double complex prdct = 0.0;
  if (plan == NULL || output == NULL || input == NULL || prdct_out == NULL ||
      plan->columns_remapped == TRUE ||
      (plan->nnz > 0U &&
       (plan->col_index == NULL || plan->values == NULL))) {
    return -1;
  }
#pragma omp parallel for default(none) schedule(static) reduction(+:prdct) \
  shared(plan, output, input)
  for (local_row = 0UL; local_row < plan->local_dim; local_row++) {
    size_t p;
    double complex sum = 0.0;
    unsigned long int global_alpha =
        plan->local_offset + local_row + 1UL;
    for (p = plan->row_ptr[local_row];
         p < plan->row_ptr[local_row + 1UL]; p++) {
      sum += plan->values[p] * input[plan->col_index[p]];
    }
    output[local_row + 1UL] += sum;
    prdct += conj(input[global_alpha]) * sum;
  }
  *prdct_out = prdct;
  return 0;
}

static int apply_direct_flat_halo(
    const struct SymmetryMatvecPlan *plan,
    double complex *output,
    const double complex *input,
    double complex *prdct_out)
{
  unsigned long int local_row;
  double complex prdct = 0.0;
  int apply_error = 0;
  size_t slot_count;
  if (plan == NULL || output == NULL || input == NULL || prdct_out == NULL ||
      plan->columns_remapped != TRUE ||
      (plan->nnz > 0U && plan->values == NULL) ||
      (size_t)plan->local_dim > SIZE_MAX - plan->halo.ghost_count) {
    return -1;
  }
  if ((plan->column_slot_width == SYMMETRY_COLUMN_U32 &&
       (plan->column_slot64 != NULL ||
        (plan->nnz > 0U && plan->column_slot32 == NULL))) ||
      (plan->column_slot_width == SYMMETRY_COLUMN_U64 &&
       (plan->column_slot32 != NULL ||
        (plan->nnz > 0U && plan->column_slot64 == NULL))) ||
      (plan->column_slot_width != SYMMETRY_COLUMN_U32 &&
       plan->column_slot_width != SYMMETRY_COLUMN_U64)) {
    return -1;
  }
  slot_count = (size_t)plan->local_dim + plan->halo.ghost_count;
  if (plan->column_slot_width == SYMMETRY_COLUMN_U32) {
#pragma omp parallel for default(none) schedule(static) \
  reduction(+:prdct) reduction(|:apply_error) \
  shared(plan, output, input, slot_count)
    for (local_row = 0UL; local_row < plan->local_dim; local_row++) {
      size_t p;
      double complex sum = 0.0;
      for (p = plan->row_ptr[local_row];
           p < plan->row_ptr[local_row + 1UL]; p++) {
        size_t slot = (size_t)plan->column_slot32[p];
        double complex input_value;
        if (slot >= slot_count) {
          apply_error = 1;
          continue;
        }
        input_value =
            slot < (size_t)plan->local_dim
                ? input[slot + 1U]
                : plan->halo.ghost_values[slot - (size_t)plan->local_dim];
        sum += plan->values[p] * input_value;
      }
      output[local_row + 1UL] += sum;
      prdct += conj(input[local_row + 1UL]) * sum;
    }
  } else {
#pragma omp parallel for default(none) schedule(static) \
  reduction(+:prdct) reduction(|:apply_error) \
  shared(plan, output, input, slot_count)
    for (local_row = 0UL; local_row < plan->local_dim; local_row++) {
      size_t p;
      double complex sum = 0.0;
      for (p = plan->row_ptr[local_row];
           p < plan->row_ptr[local_row + 1UL]; p++) {
        uint64_t raw_slot = plan->column_slot64[p];
        size_t slot;
        double complex input_value;
        if (raw_slot >= (uint64_t)slot_count) {
          apply_error = 1;
          continue;
        }
        slot = (size_t)raw_slot;
        input_value =
            slot < (size_t)plan->local_dim
                ? input[slot + 1U]
                : plan->halo.ghost_values[slot - (size_t)plan->local_dim];
        sum += plan->values[p] * input_value;
      }
      output[local_row + 1UL] += sum;
      prdct += conj(input[local_row + 1UL]) * sum;
    }
  }
  if (apply_error != 0) return -1;
  *prdct_out = prdct;
  return 0;
}

static void reference_fermion_permutation(unsigned long int state,
                                          const int *perm,
                                          unsigned int nsite,
                                          unsigned int orbitals_per_site,
                                          unsigned long int *out_state,
                                          int *sign)
{
  unsigned int mapped[sizeof(unsigned long int) * CHAR_BIT];
  unsigned int count = 0U;
  unsigned int inversions = 0U;
  unsigned int orb;
  unsigned int i, j;
  unsigned int norb = nsite * orbitals_per_site;
  *out_state = 0UL;
  for (orb = 0U; orb < norb; orb++) {
    if ((state & (1UL << orb)) != 0UL) {
      unsigned int site = orb / orbitals_per_site;
      unsigned int flavor = orb % orbitals_per_site;
      unsigned int target = orbitals_per_site * (unsigned int)perm[site] +
                            flavor;
      mapped[count++] = target;
      *out_state |= 1UL << target;
    }
  }
  for (i = 0U; i < count; i++) {
    for (j = i + 1U; j < count; j++) {
      if (mapped[i] > mapped[j]) inversions++;
    }
  }
  *sign = (inversions & 1U) == 0U ? 1 : -1;
}

static void assert_spin_permutation_states(const int *perm,
                                           unsigned int nsite,
                                           const char *label)
{
  unsigned long int state;
  unsigned long int limit = 1UL << nsite;
  for (state = 0UL; state < limit; state++) {
    unsigned long int expected = 0UL;
    unsigned int site;
    for (site = 0U; site < nsite; site++) {
      if ((state & (1UL << site)) != 0UL) {
        expected |= 1UL << (unsigned int)perm[site];
      }
    }
    if (SymmetryApplyToSpinBits(state, perm, nsite) != expected) {
      fprintf(stderr, "%s: state=%#lx expected=%#lx\n",
              label, state, expected);
      exit(1);
    }
  }
}

static void assert_fermion_permutation_states(const int *perm,
                                              unsigned int nsite,
                                              unsigned int orbitals_per_site,
                                              int model,
                                              const char *label)
{
  struct DefineList def;
  struct SymmetryTransformResult result;
  int *perm_rows_local[1];
  unsigned int norb = nsite * orbitals_per_site;
  unsigned long int state;
  unsigned long int limit = 1UL << norb;
  memset(&def, 0, sizeof(def));
  perm_rows_local[0] = (int *)perm;
  def.Nsite = nsite;
  def.NSymTrans = 1U;
  def.SymTrans = perm_rows_local;
  def.iCalcModel = model;
  for (state = 0UL; state < limit; state++) {
    unsigned long int expected_state;
    int expected_sign;
    reference_fermion_permutation(state, perm, nsite, orbitals_per_site,
                                  &expected_state, &expected_sign);
    if (SymmetryApplyToState(&def, state, 0U, &result) != 0 ||
        result.state != expected_state ||
        result.amplitude != (double)expected_sign) {
      fprintf(stderr, "%s: state=%#lx expected_state=%#lx expected_sign=%d\n",
              label, state, expected_state, expected_sign);
      exit(1);
    }
  }
}

static void enumerate_fermion_permutations(int *perm,
                                           int *used,
                                           unsigned int nsite,
                                           unsigned int depth)
{
  unsigned int target;
  if (depth == nsite) {
    assert_spin_permutation_states(perm, nsite,
                                   "spin exhaustive set-bit permutation");
    assert_fermion_permutation_states(perm, nsite, 1U, SpinlessFermion,
                                      "spinless exhaustive permutation parity");
    assert_fermion_permutation_states(perm, nsite, 2U, Hubbard,
                                      "Hubbard exhaustive permutation parity");
    return;
  }
  for (target = 0U; target < nsite; target++) {
    if (used[target] != 0) continue;
    used[target] = 1;
    perm[depth] = (int)target;
    enumerate_fermion_permutations(perm, used, nsite, depth + 1U);
    used[target] = 0;
  }
}

static void assert_exhaustive_fermion_permutation_parity(void)
{
  int perm[4] = {0, 0, 0, 0};
  int used[4] = {0, 0, 0, 0};
  enumerate_fermion_permutations(perm, used, 4U, 0U);
}

static void assert_fermion_parity_word_boundary(void)
{
  const unsigned int word_bits =
      (unsigned int)(sizeof(unsigned long int) * CHAR_BIT);
  struct DefineList def;
  struct SymmetryTransformResult result;
  int permutation[32];
  int *perm_rows_local[1];
  unsigned long int state;
  unsigned long int expected_state;
  int expected_sign;
  unsigned int site;
  if (word_bits < 64U) return;
  for (site = 0U; site < 32U; site++) permutation[site] = 31 - (int)site;
  memset(&def, 0, sizeof(def));
  perm_rows_local[0] = permutation;
  def.Nsite = 32U;
  def.NSymTrans = 1U;
  def.SymTrans = perm_rows_local;
  def.iCalcModel = Hubbard;
  state = (1UL << 1U) | (1UL << 62U);
  reference_fermion_permutation(state, permutation, 32U, 2U,
                                &expected_state, &expected_sign);
  assert_int_eq(SymmetryApplyToState(&def, state, 0U, &result), 0,
                "Hubbard parity handles mapped orbital 63");
  assert_ulong_eq(result.state, expected_state,
                  "Hubbard parity boundary transformed state");
  assert_int_eq((int)result.amplitude, expected_sign,
                "Hubbard parity boundary sign");
}

static void assert_spin_permutation_word_boundary(void)
{
  const unsigned int word_bits =
      (unsigned int)(sizeof(unsigned long int) * CHAR_BIT);
  int permutation[64];
  unsigned int site;
  if (word_bits < 64U) return;
  for (site = 0U; site < 64U; site++) permutation[site] = 63 - (int)site;
  assert_ulong_eq(SymmetryApplyToSpinBits(1UL, permutation, 64U),
                  1UL << 63U, "Spin set-bit maps to bit 63");
  assert_ulong_eq(SymmetryApplyToSpinBits(1UL << 63U, permutation, 64U),
                  1UL, "Spin set-bit reads bit 63");
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
  size_t expected_row_nnz_max = 0U;
  size_t matrix_size;
  int duplicate_found = 0;
  double difference_norm2 = 0.0;
  double legacy_norm2 = 0.0;
  double complex expected_prdct = 0.0;
  double complex flat_prdct = 0.0;
  double complex plan_prdct = 0.0;
  double complex *dense;
  double complex *input;
  double complex *flat_output;
  double complex *legacy_output;
  double complex *output;
  unsigned int *multiplicity;
  struct SymmetryMatvecBlockView view;
  struct SymmetryMatvecPlan *plan;
#ifdef _OPENMP
  int saved_dynamic = omp_get_dynamic();
  int saved_threads = omp_get_max_threads();
  /*
   * Compare the two kernels with one deterministic reduction partition.
   * The dedicated OpenMP regression below exercises block-view build/apply
   * with both one and four threads.
   */
  omp_set_dynamic(0);
  omp_set_num_threads(1);
#endif

  if (setenv("HPHI_SYMMETRY_VECTOR_EXCHANGE", "allgather", 1) != 0 ||
      ActivateSymmetryBasisDimension(X) != 0 ||
      BuildSymmetryMatvecPlan(X) != 0) {
    fprintf(stderr, "%s: plan setup failed\n", label);
    exit(1);
  }
  unsetenv("HPHI_SYMMETRY_VECTOR_EXCHANGE");
  plan = X->Sym->matvec_plan;
  assert_int_eq(plan != NULL && plan->ready == TRUE, 1, label);
  assert_ulong_eq(
      (unsigned long int)SymmetryMatvecPlanBlockCount(plan), 1UL, label);
  assert_ulong_eq(
      (unsigned long int)SymmetryMatvecPlanBlockCount(NULL), 0UL, label);
  assert_int_eq(SymmetryMatvecPlanGetBlockView(plan, 0U, &view), 0, label);
  assert_int_eq(SymmetryMatvecPlanGetBlockView(plan, 1U, &view), -1, label);
  assert_int_eq(SymmetryMatvecPlanGetBlockView(NULL, 0U, &view), -1, label);
  assert_int_eq(SymmetryMatvecPlanGetBlockView(plan, 0U, NULL), -1, label);
  assert_int_eq(SymmetryMatvecPlanGetBlockView(plan, 0U, &view), 0, label);
  assert_ulong_eq(plan->dim, X->Sym->dim, label);
  assert_ulong_eq(plan->local_offset, 0UL, label);
  assert_ulong_eq(plan->local_dim, X->Sym->dim, label);
  assert_ulong_eq((unsigned long int)plan->local_column_nnz,
                  (unsigned long int)plan->nnz, label);
  assert_ulong_eq((unsigned long int)plan->remote_column_nnz, 0UL, label);
  assert_int_eq(plan->halo.request_layout_ready, TRUE, label);
  assert_int_eq(plan->halo.ready, TRUE, label);
  assert_ulong_eq((unsigned long int)plan->halo.ghost_count, 0UL, label);
  assert_ulong_eq((unsigned long int)plan->halo.send_value_count, 0UL, label);
  assert_ulong_eq((unsigned long int)plan->halo.incoming_peer_count, 0UL,
                  label);
  assert_ulong_eq((unsigned long int)plan->halo.outgoing_peer_count, 0UL,
                  label);
  assert_int_eq(plan->halo.schedule_checksum != 0ULL, 1, label);
  assert_ulong_eq((unsigned long int)plan->column_slot_width, 32UL, label);
  assert_int_eq(plan->nnz == 0U || plan->col_index != NULL, 1, label);
  assert_int_eq(plan->column_slot32 == NULL, 1, label);
  assert_int_eq(plan->column_slot64 == NULL, 1, label);
  assert_ulong_eq(
      (unsigned long int)plan->column_storage_bytes,
      (unsigned long int)(plan->nnz * sizeof(*plan->col_index)), label);
  assert_ulong_eq(
      (unsigned long int)plan->matrix_storage_bytes,
      (unsigned long int)(
          ((size_t)plan->local_dim + 1U) * sizeof(*plan->row_ptr) +
          plan->nnz *
              (sizeof(*plan->col_index) + sizeof(*plan->values))),
      label);
  assert_ulong_eq(
      (unsigned long int)plan->allgather_nonlocal_values_per_call, 0UL, label);
  assert_ulong_eq(
      (unsigned long int)plan->allgather_payload_bytes_per_call, 0UL, label);
  assert_ulong_eq(view.local_row_begin, 0UL, label);
  assert_ulong_eq(view.local_row_count, plan->local_dim, label);
  assert_ulong_eq((unsigned long int)view.nnz,
                  (unsigned long int)plan->nnz, label);
  assert_int_eq(view.row_ptr == plan->row_ptr, 1, label);
  assert_int_eq(view.global_columns == plan->col_index, 1, label);
  assert_int_eq(view.column_slot32 == NULL, 1, label);
  assert_int_eq(view.column_slot64 == NULL, 1, label);
  assert_int_eq(view.values == plan->values, 1, label);
  assert_ulong_eq((unsigned long int)view.row_ptr[0], 0UL, label);
  assert_ulong_eq(
      (unsigned long int)view.row_ptr[view.local_row_count],
      (unsigned long int)view.nnz, label);

  if (plan->dim > SIZE_MAX / plan->dim) {
    fprintf(stderr, "%s: dense matrix size overflow\n", label);
    exit(1);
  }
  matrix_size = (size_t)plan->dim * (size_t)plan->dim;
  dense = (double complex *)calloc(matrix_size, sizeof(*dense));
  multiplicity = (unsigned int *)calloc(matrix_size, sizeof(*multiplicity));
  input = (double complex *)calloc((size_t)plan->dim + 1U, sizeof(*input));
  flat_output = (double complex *)calloc(
      (size_t)plan->dim + 1U, sizeof(*flat_output));
  legacy_output = (double complex *)calloc((size_t)plan->dim + 1U,
                                           sizeof(*legacy_output));
  output = (double complex *)calloc((size_t)plan->dim + 1U, sizeof(*output));
  if (dense == NULL || multiplicity == NULL || input == NULL ||
      flat_output == NULL ||
      legacy_output == NULL || output == NULL) {
    fprintf(stderr, "%s: dense test allocation failed\n", label);
    exit(1);
  }

  for (alpha = 1UL; alpha <= plan->dim; alpha++) {
    unsigned long int local_row = alpha - 1UL;
    size_t row_nnz = plan->row_ptr[local_row + 1UL] - plan->row_ptr[local_row];
    if (row_nnz > expected_row_nnz_max) expected_row_nnz_max = row_nnz;
    for (p = plan->row_ptr[local_row]; p < plan->row_ptr[local_row + 1UL]; p++) {
      size_t index = (size_t)(alpha - 1UL) * (size_t)plan->dim +
                     (size_t)(plan->col_index[p] - 1UL);
      dense[index] += plan->values[p];
      multiplicity[index]++;
      if (multiplicity[index] > 1U) duplicate_found = 1;
    }
  }
  assert_ulong_eq((unsigned long int)plan->row_nnz_max,
                  (unsigned long int)expected_row_nnz_max, label);
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
  if (apply_direct_flat_global(plan, flat_output, input, &flat_prdct) != 0) {
    fprintf(stderr, "%s: direct-flat apply failed\n", label);
    exit(1);
  }
  if (ApplySymmetryMatvecPlan(X, output, input, &plan_prdct) != 0) {
    fprintf(stderr, "%s: plan apply failed\n", label);
    exit(1);
  }
  assert_int_eq(
      memcmp(output, flat_output,
             ((size_t)plan->dim + 1U) * sizeof(*output)) == 0,
      1, "direct-flat/global block-view output is bitwise identical");
  assert_complex_bitwise(
      plan_prdct, flat_prdct,
      "direct-flat/global block-view prdct is bitwise identical");
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

  memset(output, 0, ((size_t)X->Sym->dim + 1U) * sizeof(*output));
  memset(flat_output, 0,
         ((size_t)X->Sym->dim + 1U) * sizeof(*flat_output));
  flat_prdct = 0.0;
  plan_prdct = 0.0;
  if (BuildSymmetryMatvecPlan(X) != 0) {
    fprintf(stderr, "%s: default serial halo plan setup failed\n", label);
    exit(1);
  }
  plan = X->Sym->matvec_plan;
  assert_int_eq(plan != NULL && plan->ready == TRUE, 1, label);
  assert_int_eq(plan->columns_remapped, TRUE, label);
  assert_int_eq(X->Sym->vector_exchange_mode,
                SYMMETRY_VECTOR_EXCHANGE_HALO, label);
  assert_int_eq(X->Sym->mpi_full_v1 == NULL, 1, label);
  assert_int_eq(plan->col_index == NULL, 1, label);
  assert_int_eq(plan->nnz == 0U || plan->column_slot32 != NULL, 1, label);
  assert_int_eq(plan->column_slot64 == NULL, 1, label);
  assert_int_eq(plan->halo.ghost_global_index == NULL, 1, label);
  assert_ulong_eq(
      (unsigned long int)SymmetryMatvecPlanBlockCount(plan), 1UL, label);
  assert_int_eq(SymmetryMatvecPlanGetBlockView(plan, 0U, &view), 0, label);
  assert_ulong_eq(view.local_row_begin, 0UL, label);
  assert_ulong_eq(view.local_row_count, plan->local_dim, label);
  assert_ulong_eq((unsigned long int)view.nnz,
                  (unsigned long int)plan->nnz, label);
  assert_int_eq(view.row_ptr == plan->row_ptr, 1, label);
  assert_int_eq(view.global_columns == NULL, 1, label);
  assert_int_eq(view.column_slot32 == plan->column_slot32, 1, label);
  assert_int_eq(view.column_slot64 == plan->column_slot64, 1, label);
  assert_int_eq(view.values == plan->values, 1, label);
  assert_ulong_eq((unsigned long int)view.row_ptr[0], 0UL, label);
  assert_ulong_eq(
      (unsigned long int)view.row_ptr[view.local_row_count],
      (unsigned long int)view.nnz, label);
  assert_ulong_eq(
      (unsigned long int)plan->column_storage_bytes,
      (unsigned long int)(plan->nnz * sizeof(*plan->column_slot32)), label);
  assert_ulong_eq(
      (unsigned long int)plan->matrix_storage_bytes,
      (unsigned long int)(
          ((size_t)plan->local_dim + 1U) * sizeof(*plan->row_ptr) +
          plan->nnz *
              (sizeof(*plan->column_slot32) + sizeof(*plan->values))),
      label);
  for (p = 0U; p < plan->nnz; p++) {
    assert_int_eq(symmetry_plan_column_slot(plan, p) <
                      (uint64_t)plan->local_dim +
                          (uint64_t)plan->halo.ghost_count,
                  1, label);
  }
  assert_int_eq(ExchangeSymmetryVectorHalo(&plan->halo, input), 0, label);
  if (apply_direct_flat_halo(
          plan, flat_output, input, &flat_prdct) != 0) {
    fprintf(stderr, "%s: direct-flat halo apply failed\n", label);
    exit(1);
  }
  assert_int_eq(ApplySymmetryMatvecPlan(X, output, input, &plan_prdct), -1,
                "global-column apply rejects remapped plan");
  assert_int_eq(
      ApplySymmetryMatvecPlanHalo(X, output, input, &plan_prdct), 0, label);
  assert_int_eq(
      memcmp(output, flat_output,
             ((size_t)plan->dim + 1U) * sizeof(*output)) == 0,
      1, "direct-flat/halo block-view output is bitwise identical");
  assert_complex_bitwise(
      plan_prdct, flat_prdct,
      "direct-flat/halo block-view prdct is bitwise identical");
  for (alpha = 1UL; alpha <= plan->dim; alpha++) {
    assert_complex_close(output[alpha], legacy_output[alpha], 1.0e-12, label);
  }
  assert_complex_close(plan_prdct, expected_prdct, 1.0e-12, label);
  assert_ulong_eq((unsigned long int)plan->halo.exchange_calls, 1UL, label);

  free(dense);
  free(multiplicity);
  free(input);
  free(flat_output);
  free(legacy_output);
  free(output);
#ifdef _OPENMP
  omp_set_num_threads(saved_threads);
  omp_set_dynamic(saved_dynamic);
#endif
}

static void assert_spin_plan(unsigned int nsite,
                             unsigned int nup,
                             unsigned int momentum_index,
                             double diagonal_coupling,
                             const char *label)
{
  struct BindStruct X;
  setup_bind(&X, nsite, nup, momentum_index);
  if (diagonal_coupling != 0.0)
    set_ising_ring_diagonal(&X.Def, nsite, diagonal_coupling);
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
  set_ising_ring_diagonal(&X.Def, 6, 0.37);
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

#ifdef _OPENMP
static int basis_vector_fields_equal(
    const struct SymmetryBasisVector *lhs,
    const struct SymmetryBasisVector *rhs)
{
  return lhs->rep_state == rhs->rep_state &&
         lhs->orbit_size == rhs->orbit_size &&
         lhs->stabilizer_size == rhs->stabilizer_size &&
         memcmp(&lhs->norm, &rhs->norm, sizeof(lhs->norm)) == 0 &&
         memcmp(&lhs->stabilizer_character_sum,
                &rhs->stabilizer_character_sum,
                sizeof(lhs->stabilizer_character_sum)) == 0 &&
         memcmp(&lhs->diagonal, &rhs->diagonal,
                sizeof(lhs->diagonal)) == 0;
}

static void assert_parallel_basis_matches_serial(const char *label)
{
  struct BindStruct X;
  struct SymmetryBasisVector *serial_basis;
  unsigned long int serial_dim;
  unsigned long int index;
  unsigned long long serial_raw_states;
  unsigned long long serial_candidates;
  unsigned long long serial_survivors;
  unsigned long long serial_transform_calls;
  int saved_dynamic = omp_get_dynamic();
  int saved_threads = omp_get_max_threads();

  omp_set_dynamic(0);
  omp_set_num_threads(1);
  setup_hubbard_bind(&X, 6, 3, 3, 1);
  setup_hubbard_coulomb_intra(&X.Def, 6, 0.5);
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: serial basis setup failed\n", label);
    exit(1);
  }
  serial_dim = X.Sym->dim;
  serial_raw_states = X.Sym->basis_raw_states;
  serial_candidates = X.Sym->basis_representative_candidates;
  serial_survivors = X.Sym->basis_compatible_survivors;
  serial_transform_calls = X.Sym->basis_transform_calls;
  serial_basis = (struct SymmetryBasisVector *)calloc(
      (size_t)serial_dim + 1U, sizeof(*serial_basis));
  if (serial_basis == NULL) {
    fprintf(stderr, "%s: serial basis snapshot allocation failed\n", label);
    exit(1);
  }
  for (index = 1UL; index <= serial_dim; index++) {
    serial_basis[index] = X.Sym->basis[index];
  }
  FreeSymmetryBasis(X.Sym);
  X.Sym = NULL;

  omp_set_num_threads(4);
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: parallel basis setup failed\n", label);
    exit(1);
  }
  assert_ulong_eq(X.Sym->dim, serial_dim, label);
  assert_int_eq(X.Sym->basis_raw_states == serial_raw_states, 1, label);
  assert_int_eq(X.Sym->basis_representative_candidates == serial_candidates,
                1, label);
  assert_int_eq(X.Sym->basis_compatible_survivors == serial_survivors,
                1, label);
  assert_int_eq(X.Sym->basis_transform_calls == serial_transform_calls,
                1, label);
  for (index = 1UL; index <= serial_dim; index++) {
    assert_int_eq(basis_vector_fields_equal(&X.Sym->basis[index],
                                            &serial_basis[index]),
                  1, label);
  }

  free(serial_basis);
  FreeSymmetryBasis(X.Sym);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
  omp_set_num_threads(saved_threads);
  omp_set_dynamic(saved_dynamic);
}

static void assert_parallel_plan_matches_serial(const char *label)
{
  struct BindStruct X;
  struct SymmetryMatvecPlan *parallel_plan;
  unsigned long int index;
  size_t row_ptr_count;
  size_t serial_nnz;
  size_t serial_row_nnz_max;
  size_t vector_count;
  size_t *serial_row_ptr;
  unsigned long int *serial_col_index;
  double complex parallel_prdct = 0.0;
  double complex serial_prdct = 0.0;
  double complex *input;
  double complex *parallel_output;
  double complex *serial_output;
  double complex *serial_values;
  int saved_dynamic = omp_get_dynamic();
  int saved_threads = omp_get_max_threads();

  omp_set_dynamic(0);
  omp_set_num_threads(1);
  setup_hubbard_bind(&X, 6, 3, 3, 1);
  setup_hubbard_transfer_ring(&X.Def, 6);
  setup_hubbard_coulomb_intra(&X.Def, 6, 0.5);
  if (setenv("HPHI_SYMMETRY_VECTOR_EXCHANGE", "allgather", 1) != 0 ||
      BuildSymmetryBasis(&X) != 0 ||
      ActivateSymmetryBasisDimension(&X) != 0 ||
      BuildSymmetryMatvecPlan(&X) != 0) {
    fprintf(stderr, "%s: serial plan setup failed\n", label);
    exit(1);
  }

  row_ptr_count = (size_t)X.Sym->matvec_plan->local_dim + 1U;
  serial_nnz = X.Sym->matvec_plan->nnz;
  serial_row_nnz_max = X.Sym->matvec_plan->row_nnz_max;
  serial_row_ptr = (size_t *)malloc(row_ptr_count * sizeof(*serial_row_ptr));
  serial_col_index = (unsigned long int *)malloc(
      serial_nnz * sizeof(*serial_col_index));
  serial_values = (double complex *)malloc(serial_nnz * sizeof(*serial_values));
  vector_count = (size_t)X.Sym->matvec_plan->local_dim + 1U;
  input = (double complex *)calloc(vector_count, sizeof(*input));
  serial_output =
      (double complex *)calloc(vector_count, sizeof(*serial_output));
  parallel_output =
      (double complex *)calloc(vector_count, sizeof(*parallel_output));
  if (serial_row_ptr == NULL ||
      input == NULL || serial_output == NULL || parallel_output == NULL ||
      (serial_nnz > 0U &&
       (serial_col_index == NULL || serial_values == NULL))) {
    fprintf(stderr, "%s: serial plan snapshot allocation failed\n", label);
    exit(1);
  }
  memcpy(serial_row_ptr, X.Sym->matvec_plan->row_ptr,
         row_ptr_count * sizeof(*serial_row_ptr));
  if (serial_nnz > 0U) {
    memcpy(serial_col_index, X.Sym->matvec_plan->col_index,
           serial_nnz * sizeof(*serial_col_index));
    memcpy(serial_values, X.Sym->matvec_plan->values,
           serial_nnz * sizeof(*serial_values));
  }
  for (index = 1UL; index <= X.Sym->matvec_plan->local_dim; index++) {
    input[index] =
        0.125 * (double)index + I * 0.0625 * (double)(index + 1UL);
  }
  if (ApplySymmetryMatvecPlan(
          &X, serial_output, input, &serial_prdct) != 0) {
    fprintf(stderr, "%s: serial block-view apply failed\n", label);
    exit(1);
  }

  omp_set_num_threads(4);
  if (BuildSymmetryMatvecPlan(&X) != 0) {
    fprintf(stderr, "%s: parallel plan setup failed\n", label);
    exit(1);
  }
  parallel_plan = X.Sym->matvec_plan;
  unsetenv("HPHI_SYMMETRY_VECTOR_EXCHANGE");
  assert_ulong_eq((unsigned long int)parallel_plan->nnz,
                  (unsigned long int)serial_nnz, label);
  assert_ulong_eq((unsigned long int)parallel_plan->row_nnz_max,
                  (unsigned long int)serial_row_nnz_max, label);
  assert_int_eq(memcmp(parallel_plan->row_ptr, serial_row_ptr,
                       row_ptr_count * sizeof(*serial_row_ptr)) == 0,
                1, label);
  if (serial_nnz > 0U) {
    assert_int_eq(memcmp(parallel_plan->col_index, serial_col_index,
                         serial_nnz * sizeof(*serial_col_index)) == 0,
                  1, label);
    assert_int_eq(memcmp(parallel_plan->values, serial_values,
                         serial_nnz * sizeof(*serial_values)) == 0,
                  1, label);
  }
  if (ApplySymmetryMatvecPlan(
          &X, parallel_output, input, &parallel_prdct) != 0) {
    fprintf(stderr, "%s: parallel block-view apply failed\n", label);
    exit(1);
  }
  assert_int_eq(
      memcmp(serial_output, parallel_output,
             vector_count * sizeof(*serial_output)) == 0,
      1, "OpenMP 1/4 block-view output is bitwise identical");
  assert_complex_close(
      parallel_prdct, serial_prdct, 1.0e-12,
      "OpenMP 1/4 block-view prdct is numerically identical");

  free(serial_row_ptr);
  free(serial_col_index);
  free(serial_values);
  free(input);
  free(serial_output);
  free(parallel_output);
  FreeSymmetryBasis(X.Sym);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
  omp_set_num_threads(saved_threads);
  omp_set_dynamic(saved_dynamic);
}
#endif

static void assert_halo_plan_exact(
    const struct SymmetryVectorHaloPlan *got,
    const struct SymmetryVectorHaloPlan *expected,
    const char *label)
{
  size_t rank_bytes;
  assert_int_eq(got->request_layout_ready, expected->request_layout_ready,
                label);
  assert_int_eq(got->ready, expected->ready, label);
  assert_int_eq(got->nrank, expected->nrank, label);
  assert_int_eq(got->rank, expected->rank, label);
  assert_ulong_eq(got->dim, expected->dim, label);
  assert_ulong_eq(got->local_offset, expected->local_offset, label);
  assert_ulong_eq(got->local_dim, expected->local_dim, label);
  assert_ulong_eq((unsigned long int)got->ghost_count,
                  (unsigned long int)expected->ghost_count, label);
  assert_ulong_eq((unsigned long int)got->send_value_count,
                  (unsigned long int)expected->send_value_count, label);
  assert_ulong_eq((unsigned long int)got->incoming_peer_count,
                  (unsigned long int)expected->incoming_peer_count, label);
  assert_ulong_eq((unsigned long int)got->outgoing_peer_count,
                  (unsigned long int)expected->outgoing_peer_count, label);
  assert_ulong_eq((unsigned long int)got->max_recv_from_peer,
                  (unsigned long int)expected->max_recv_from_peer, label);
  assert_ulong_eq((unsigned long int)got->max_send_to_peer,
                  (unsigned long int)expected->max_send_to_peer, label);
  assert_ulong_eq((unsigned long int)got->topology_scratch_bytes,
                  (unsigned long int)expected->topology_scratch_bytes, label);
  assert_ulong_eq((unsigned long int)got->schedule_bytes,
                  (unsigned long int)expected->schedule_bytes, label);
  assert_ulong_eq((unsigned long int)got->runtime_buffer_bytes,
                  (unsigned long int)expected->runtime_buffer_bytes, label);
  assert_int_eq(got->schedule_checksum == expected->schedule_checksum,
                1, label);
  rank_bytes = (size_t)got->nrank * sizeof(*got->recv_counts);
  assert_int_eq(memcmp(got->recv_counts, expected->recv_counts,
                       rank_bytes) == 0, 1, label);
  assert_int_eq(memcmp(got->recv_displs, expected->recv_displs,
                       rank_bytes) == 0, 1, label);
  assert_int_eq(memcmp(got->send_counts, expected->send_counts,
                       rank_bytes) == 0, 1, label);
  assert_int_eq(memcmp(got->send_displs, expected->send_displs,
                       rank_bytes) == 0, 1, label);
  if (got->ghost_count > 0U) {
    assert_int_eq(
        memcmp(got->ghost_global_index, expected->ghost_global_index,
               got->ghost_count * sizeof(*got->ghost_global_index)) == 0,
        1, label);
  }
  if (got->send_value_count > 0U) {
    assert_int_eq(
        memcmp(got->send_local_index, expected->send_local_index,
               got->send_value_count * sizeof(*got->send_local_index)) == 0,
        1, label);
  }
}

#ifdef MPI
static void assert_mpi_column_spans_exact(const char *label)
{
  const unsigned long int dim = 10UL;
  struct SymmetryVectorHaloPlan flat_halo;
  struct SymmetryVectorHaloPlan split_halo;
  struct SymmetryMatvecPlan flat_plan;
  struct SymmetryMatvecPlan split_plan;
  struct SymmetryBasisRuntime sym;
  struct BindStruct X;
  struct SymmetryGlobalColumnSpan flat_span;
  struct SymmetryGlobalColumnSpan split_spans[7];
  unsigned long int base;
  unsigned long int remainder;
  unsigned long int local_offset;
  unsigned long int local_dim;
  unsigned long int local_row;
  unsigned long int *columns = NULL;
  unsigned long int *flat_columns = NULL;
  unsigned long int *split_columns = NULL;
  size_t *row_ptr = NULL;
  double complex *values = NULL;
  double complex *local_vector = NULL;
  double complex *flat_output = NULL;
  double complex *split_output = NULL;
  double complex flat_prdct = 0.0;
  double complex split_prdct = 0.0;
  size_t flat_local_columns = 0U;
  size_t flat_remote_columns = 0U;
  size_t split_local_columns = 0U;
  size_t split_remote_columns = 0U;
  size_t nnz;
  size_t first_count;
  size_t middle_count;
  size_t last_count;
  size_t index;
#ifdef _OPENMP
  int saved_dynamic = omp_get_dynamic();
  int saved_threads = omp_get_max_threads();
  omp_set_dynamic(0);
  omp_set_num_threads(1);
#endif

  memset(&flat_halo, 0, sizeof(flat_halo));
  memset(&split_halo, 0, sizeof(split_halo));
  memset(&flat_plan, 0, sizeof(flat_plan));
  memset(&split_plan, 0, sizeof(split_plan));
  memset(&sym, 0, sizeof(sym));
  memset(&X, 0, sizeof(X));
  base = dim / (unsigned long int)nproc;
  remainder = dim % (unsigned long int)nproc;
  local_dim =
      base + ((unsigned long int)myrank < remainder ? 1UL : 0UL);
  local_offset =
      base * (unsigned long int)myrank +
      ((unsigned long int)myrank < remainder
           ? (unsigned long int)myrank
           : remainder);
  nnz = (size_t)local_dim * 3U;
  row_ptr = (size_t *)calloc((size_t)local_dim + 1U, sizeof(*row_ptr));
  local_vector = (double complex *)calloc(
      (size_t)local_dim + 1U, sizeof(*local_vector));
  flat_output = (double complex *)calloc(
      (size_t)local_dim + 1U, sizeof(*flat_output));
  split_output = (double complex *)calloc(
      (size_t)local_dim + 1U, sizeof(*split_output));
  if (nnz > 0U) {
    columns = (unsigned long int *)malloc(nnz * sizeof(*columns));
    flat_columns =
        (unsigned long int *)malloc(nnz * sizeof(*flat_columns));
    split_columns =
        (unsigned long int *)malloc(nnz * sizeof(*split_columns));
    values = (double complex *)malloc(nnz * sizeof(*values));
  }
  if (row_ptr == NULL || local_vector == NULL || flat_output == NULL ||
      split_output == NULL ||
      (nnz > 0U && (columns == NULL || flat_columns == NULL ||
                    split_columns == NULL || values == NULL))) {
    fprintf(stderr, "%s: MPI column-span allocation failed\n", label);
    exit(1);
  }
  for (local_row = 0UL; local_row < local_dim; local_row++) {
    unsigned long int global_row = local_offset + local_row + 1UL;
    size_t row_begin = (size_t)local_row * 3U;
    row_ptr[local_row] = row_begin;
    columns[row_begin] = global_row;
    columns[row_begin + 1U] = global_row == dim ? 1UL : global_row + 1UL;
    columns[row_begin + 2U] = global_row;
    local_vector[local_row + 1UL] =
        0.25 * (double)global_row +
        I * 0.125 * (double)(global_row + 1UL);
  }
  row_ptr[local_dim] = nnz;
  for (index = 0U; index < nnz; index++) {
    values[index] =
        1.0 + 0.03125 * (double)(index + 1U) +
        I * 0.015625 * (double)(index + 2U);
  }
  if (nnz > 0U) {
    memcpy(flat_columns, columns, nnz * sizeof(*columns));
    memcpy(split_columns, columns, nnz * sizeof(*columns));
  }

  first_count = nnz < 2U ? nnz : 2U;
  last_count = nnz > first_count ? 1U : 0U;
  middle_count = nnz - first_count - last_count;
  flat_span.columns = columns;
  flat_span.count = nnz;
  split_spans[0].columns = NULL;
  split_spans[0].count = 0U;
  split_spans[1].columns = first_count > 0U ? columns : NULL;
  split_spans[1].count = first_count;
  split_spans[2].columns = NULL;
  split_spans[2].count = 0U;
  split_spans[3].columns =
      middle_count > 0U ? columns + first_count : NULL;
  split_spans[3].count = middle_count;
  split_spans[4].columns = NULL;
  split_spans[4].count = 0U;
  split_spans[5].columns =
      last_count > 0U ? columns + nnz - last_count : NULL;
  split_spans[5].count = last_count;
  split_spans[6].columns = NULL;
  split_spans[6].count = 0U;

  assert_int_eq(
      BuildSymmetryVectorHaloPlan(
          &flat_halo, dim, local_offset, local_dim,
          &flat_span, 1U, nproc, myrank,
          &flat_local_columns, &flat_remote_columns),
      0, label);
  assert_int_eq(
      BuildSymmetryVectorHaloPlan(
          &split_halo, dim, local_offset, local_dim,
          split_spans, sizeof(split_spans) / sizeof(split_spans[0]),
          nproc, myrank, &split_local_columns, &split_remote_columns),
      0, label);
  assert_ulong_eq((unsigned long int)split_local_columns,
                  (unsigned long int)flat_local_columns, label);
  assert_ulong_eq((unsigned long int)split_remote_columns,
                  (unsigned long int)flat_remote_columns, label);
  assert_halo_plan_exact(&split_halo, &flat_halo, label);
  assert_int_eq(
      ExchangeSymmetryVectorHalo(&flat_halo, local_vector), 0, label);
  assert_int_eq(
      ExchangeSymmetryVectorHalo(&split_halo, local_vector), 0, label);
  if (flat_halo.ghost_count > 0U) {
    assert_int_eq(
        memcmp(split_halo.ghost_values, flat_halo.ghost_values,
               flat_halo.ghost_count *
                   sizeof(*flat_halo.ghost_values)) == 0,
        1, label);
  }

  flat_plan.ready = TRUE;
  flat_plan.columns_remapped = FALSE;
  flat_plan.dim = dim;
  flat_plan.local_offset = local_offset;
  flat_plan.local_dim = local_dim;
  flat_plan.block_count = 1U;
  flat_plan.nnz = nnz;
  flat_plan.row_ptr = row_ptr;
  flat_plan.col_index = flat_columns;
  flat_plan.values = values;
  flat_plan.halo = flat_halo;
  split_plan = flat_plan;
  split_plan.col_index = split_columns;
  split_plan.halo = split_halo;
  assert_int_eq(RemapSymmetryMatvecPlanColumns(&flat_plan), 0, label);
  assert_int_eq(RemapSymmetryMatvecPlanColumns(&split_plan), 0, label);
  assert_ulong_eq((unsigned long int)split_plan.column_slot_width,
                  (unsigned long int)flat_plan.column_slot_width, label);
  if (flat_plan.column_slot_width == SYMMETRY_COLUMN_U32 && nnz > 0U) {
    assert_int_eq(
        memcmp(split_plan.column_slot32, flat_plan.column_slot32,
               nnz * sizeof(*flat_plan.column_slot32)) == 0,
        1, label);
  } else if (nnz > 0U) {
    assert_int_eq(
        memcmp(split_plan.column_slot64, flat_plan.column_slot64,
               nnz * sizeof(*flat_plan.column_slot64)) == 0,
        1, label);
  }

  sym.enabled = TRUE;
  sym.dim = dim;
  sym.local_offset = local_offset;
  sym.local_dim = local_dim;
  X.Sym = &sym;
  sym.matvec_plan = &flat_plan;
  assert_int_eq(
      ApplySymmetryMatvecPlanHalo(
          &X, flat_output, local_vector, &flat_prdct),
      0, label);
  sym.matvec_plan = &split_plan;
  assert_int_eq(
      ApplySymmetryMatvecPlanHalo(
          &X, split_output, local_vector, &split_prdct),
      0, label);
  assert_int_eq(
      memcmp(split_output, flat_output,
             ((size_t)local_dim + 1U) * sizeof(*flat_output)) == 0,
      1, label);
  assert_complex_bitwise(split_prdct, flat_prdct, label);

  free(flat_plan.column_slot32);
  free(flat_plan.column_slot64);
  free(split_plan.column_slot32);
  free(split_plan.column_slot64);
  free(columns);
  free(row_ptr);
  free(values);
  free(local_vector);
  free(flat_output);
  free(split_output);
  FreeSymmetryVectorHaloPlan(&flat_halo);
  FreeSymmetryVectorHaloPlan(&split_halo);
#ifdef _OPENMP
  omp_set_num_threads(saved_threads);
  omp_set_dynamic(saved_dynamic);
#endif
}
#endif

static void assert_vector_owner_and_request_layout(const char *label)
{
  const unsigned long int columns[] = {5UL, 1UL, 1UL, 4UL,
                                       8UL, 10UL, 8UL};
  const unsigned long int local_columns_only[] = {5UL, 6UL, 7UL};
  const unsigned long int remote_columns_only[] = {1UL, 4UL, 8UL, 10UL};
  const unsigned long int expected_ghosts[] = {1UL, 4UL, 8UL, 10UL};
  const struct SymmetryGlobalColumnSpan flat_span = {
      columns, sizeof(columns) / sizeof(columns[0])};
  const struct SymmetryGlobalColumnSpan split_spans[] = {
      {NULL, 0U},
      {columns, 2U},
      {NULL, 0U},
      {columns + 2U, 3U},
      {columns + 5U, 2U},
      {NULL, 0U}};
  const struct SymmetryGlobalColumnSpan local_span = {
      local_columns_only,
      sizeof(local_columns_only) / sizeof(local_columns_only[0])};
  const struct SymmetryGlobalColumnSpan remote_span = {
      remote_columns_only,
      sizeof(remote_columns_only) / sizeof(remote_columns_only[0])};
  const struct SymmetryGlobalColumnSpan invalid_span = {NULL, 1U};
  struct SymmetryVectorHaloPlan first;
  struct SymmetryVectorHaloPlan second;
  struct SymmetryVectorHaloPlan empty;
  size_t local_columns = 0U;
  size_t remote_columns = 0U;
  size_t index;

  memset(&first, 0, sizeof(first));
  memset(&second, 0, sizeof(second));
  memset(&empty, 0, sizeof(empty));
  assert_int_eq(SymmetryVectorOwnerOfGlobalIndex(10UL, 3, 1UL), 0, label);
  assert_int_eq(SymmetryVectorOwnerOfGlobalIndex(10UL, 3, 4UL), 0, label);
  assert_int_eq(SymmetryVectorOwnerOfGlobalIndex(10UL, 3, 5UL), 1, label);
  assert_int_eq(SymmetryVectorOwnerOfGlobalIndex(10UL, 3, 7UL), 1, label);
  assert_int_eq(SymmetryVectorOwnerOfGlobalIndex(10UL, 3, 8UL), 2, label);
  assert_int_eq(SymmetryVectorOwnerOfGlobalIndex(10UL, 3, 10UL), 2, label);
  assert_int_eq(SymmetryVectorOwnerOfGlobalIndex(2UL, 4, 1UL), 0, label);
  assert_int_eq(SymmetryVectorOwnerOfGlobalIndex(2UL, 4, 2UL), 1, label);
  assert_int_eq(SymmetryVectorOwnerOfGlobalIndex(2UL, 4, 0UL), -1, label);
  assert_int_eq(SymmetryVectorOwnerOfGlobalIndex(2UL, 4, 3UL), -1, label);

  assert_int_eq(
      BuildSymmetryVectorHaloPlan(
          &first, 10UL, 4UL, 3UL, &flat_span, 1U, 3, 1,
          &local_columns, &remote_columns),
      0, label);
  assert_ulong_eq((unsigned long int)local_columns, 1UL, label);
  assert_ulong_eq((unsigned long int)remote_columns, 6UL, label);
  assert_int_eq(first.request_layout_ready, TRUE, label);
  assert_int_eq(first.ready, FALSE, label);
  assert_ulong_eq((unsigned long int)first.ghost_count, 4UL, label);
  assert_ulong_eq((unsigned long int)first.incoming_peer_count, 2UL, label);
  assert_ulong_eq((unsigned long int)first.max_recv_from_peer, 2UL, label);
  assert_int_eq(first.recv_counts[0], 2, label);
  assert_int_eq(first.recv_counts[1], 0, label);
  assert_int_eq(first.recv_counts[2], 2, label);
  assert_int_eq(first.recv_displs[0], 0, label);
  assert_int_eq(first.recv_displs[1], 2, label);
  assert_int_eq(first.recv_displs[2], 2, label);
  for (index = 0U; index < first.ghost_count; index++) {
    assert_ulong_eq(first.ghost_global_index[index],
                    expected_ghosts[index], label);
  }

  local_columns = 0U;
  remote_columns = 0U;
  assert_int_eq(
      BuildSymmetryVectorHaloPlan(
          &second, 10UL, 4UL, 3UL, split_spans,
          sizeof(split_spans) / sizeof(split_spans[0]), 3, 1,
          &local_columns, &remote_columns),
      0, label);
  assert_ulong_eq((unsigned long int)local_columns, 1UL, label);
  assert_ulong_eq((unsigned long int)remote_columns, 6UL, label);
  assert_halo_plan_exact(&second, &first, label);
  FreeSymmetryVectorHaloPlan(&first);
  FreeSymmetryVectorHaloPlan(&second);

  assert_int_eq(
      BuildSymmetryVectorHaloPlan(
          &empty, 10UL, 4UL, 3UL, NULL, 0U, 3, 1,
          &local_columns, &remote_columns),
      0, label);
  assert_ulong_eq((unsigned long int)local_columns, 0UL, label);
  assert_ulong_eq((unsigned long int)remote_columns, 0UL, label);
  assert_ulong_eq((unsigned long int)empty.ghost_count, 0UL, label);
  FreeSymmetryVectorHaloPlan(&empty);

  assert_int_eq(
      BuildSymmetryVectorHaloPlan(
          &empty, 10UL, 4UL, 3UL, &local_span, 1U, 3, 1,
          &local_columns, &remote_columns),
      0, label);
  assert_ulong_eq((unsigned long int)local_columns, 3UL, label);
  assert_ulong_eq((unsigned long int)remote_columns, 0UL, label);
  assert_ulong_eq((unsigned long int)empty.ghost_count, 0UL, label);
  FreeSymmetryVectorHaloPlan(&empty);

  assert_int_eq(
      BuildSymmetryVectorHaloPlan(
          &empty, 10UL, 4UL, 3UL, &remote_span, 1U, 3, 1,
          &local_columns, &remote_columns),
      0, label);
  assert_ulong_eq((unsigned long int)local_columns, 0UL, label);
  assert_ulong_eq((unsigned long int)remote_columns, 4UL, label);
  assert_ulong_eq((unsigned long int)empty.ghost_count, 4UL, label);
  FreeSymmetryVectorHaloPlan(&empty);

  assert_int_eq(
      BuildSymmetryVectorHaloPlan(
          &empty, 10UL, 4UL, 3UL, &invalid_span, 1U, 3, 1,
          &local_columns, &remote_columns),
      -1, label);

  memset(&first, 0, sizeof(first));
  assert_int_eq(
      BuildSymmetryVectorHaloPlan(
          &first, 10UL, 3UL, 3UL, &flat_span, 1U, 3, 1,
          &local_columns, &remote_columns),
      -1, label);
}

static void assert_remote_topology_plan(const char *label)
{
  struct BindStruct X;
  struct SymmetryMatvecPlan *plan;
  setup_bind(&X, 6, 3, 1);
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: BuildSymmetryBasis failed\n", label);
    exit(1);
  }
  nproc = 4;
  myrank = 0;
  if (setenv("HPHI_SYMMETRY_VECTOR_EXCHANGE", "allgather", 1) != 0 ||
      ActivateSymmetryBasisDimension(&X) != 0 ||
      BuildSymmetryMatvecPlan(&X) != 0) {
    fprintf(stderr, "%s: remote topology plan setup failed\n", label);
    exit(1);
  }
  unsetenv("HPHI_SYMMETRY_VECTOR_EXCHANGE");
  plan = X.Sym->matvec_plan;
  assert_int_eq(plan != NULL && plan->ready == TRUE, 1, label);
  assert_ulong_eq(
      (unsigned long int)(plan->local_column_nnz + plan->remote_column_nnz),
      (unsigned long int)plan->nnz, label);
  assert_int_eq(plan->remote_column_nnz > 0U, 1, label);
  assert_int_eq(plan->halo.request_layout_ready, TRUE, label);
  assert_int_eq(plan->halo.ready, FALSE, label);
  assert_int_eq(plan->halo.ghost_count > 0U, 1, label);
  assert_int_eq(plan->halo.ghost_count <= plan->remote_column_nnz, 1, label);
  assert_int_eq(plan->halo.incoming_peer_count > 0U, 1, label);
  assert_int_eq(plan->halo.max_recv_from_peer > 0U, 1, label);
  assert_ulong_eq((unsigned long int)plan->halo.send_value_count, 0UL, label);
  assert_ulong_eq((unsigned long int)plan->halo.outgoing_peer_count, 0UL,
                  label);
  assert_int_eq(plan->halo.topology_scratch_bytes > 0U, 1, label);
  assert_int_eq(plan->halo.schedule_checksum != 0ULL, 1, label);
  assert_ulong_eq((unsigned long int)plan->column_slot_width, 32UL, label);
  assert_ulong_eq(
      (unsigned long int)plan->allgather_nonlocal_values_per_call,
      X.Sym->dim - X.Sym->local_dim, label);
  assert_ulong_eq(
      (unsigned long int)plan->allgather_payload_bytes_per_call,
      (X.Sym->dim - X.Sym->local_dim) * sizeof(double complex), label);
  nproc = 1;
  myrank = 0;
  FreeSymmetryBasis(X.Sym);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

static void assert_mixed_column_remap(const char *label)
{
  const unsigned long int initial_columns[] =
      {5UL, 1UL, 4UL, 8UL, 10UL, 7UL};
  const unsigned long int expected_slots[] = {0UL, 3UL, 4UL, 5UL, 6UL, 2UL};
  unsigned long int missing_ghost_column[] = {9UL};
  unsigned long int ghosts[] = {1UL, 4UL, 8UL, 10UL};
  unsigned long int *columns;
  struct SymmetryMatvecPlan plan;
  struct SymmetryMatvecPlan invalid_plan;
  size_t index;
  columns = (unsigned long int *)malloc(sizeof(initial_columns));
  if (columns == NULL) {
    fprintf(stderr, "%s: column allocation failed\n", label);
    exit(1);
  }
  memcpy(columns, initial_columns, sizeof(initial_columns));
  memset(&plan, 0, sizeof(plan));
  plan.dim = 10UL;
  plan.local_offset = 4UL;
  plan.local_dim = 3UL;
  plan.nnz = sizeof(initial_columns) / sizeof(initial_columns[0]);
  plan.col_index = columns;
  plan.halo.ghost_count = sizeof(ghosts) / sizeof(ghosts[0]);
  plan.halo.ghost_global_index = ghosts;
  assert_int_eq(RemapSymmetryMatvecPlanColumns(&plan), 0, label);
  assert_int_eq(plan.columns_remapped, TRUE, label);
  assert_int_eq(plan.col_index == NULL, 1, label);
  assert_ulong_eq(
      (unsigned long int)plan.column_slot_width,
      (unsigned long int)SYMMETRY_COLUMN_U32, label);
  assert_int_eq(plan.column_slot32 != NULL, 1, label);
  assert_int_eq(plan.column_slot64 == NULL, 1, label);
  for (index = 0U; index < plan.nnz; index++) {
    assert_ulong_eq((unsigned long int)symmetry_plan_column_slot(&plan, index),
                    expected_slots[index], label);
  }
  assert_int_eq(RemapSymmetryMatvecPlanColumns(&plan), -1,
                "column remap rejects a second in-place remap");

  memset(&invalid_plan, 0, sizeof(invalid_plan));
  invalid_plan.dim = plan.dim;
  invalid_plan.local_offset = plan.local_offset;
  invalid_plan.local_dim = plan.local_dim;
  invalid_plan.nnz = 1U;
  invalid_plan.col_index = missing_ghost_column;
  invalid_plan.halo.ghost_count = plan.halo.ghost_count;
  invalid_plan.halo.ghost_global_index = ghosts;
  assert_int_eq(RemapSymmetryMatvecPlanColumns(&invalid_plan), -1,
                "column remap rejects a missing ghost index");
  free(plan.column_slot32);
  free(plan.column_slot64);
}

static void assert_column_slot_width_boundaries(const char *label)
{
  struct SymmetryMatvecPlan plan32;
  unsigned long int *column32 =
      (unsigned long int *)malloc(sizeof(*column32));
  if (column32 == NULL) {
    fprintf(stderr, "%s: 32-bit boundary allocation failed\n", label);
    exit(1);
  }
  *column32 = (unsigned long int)UINT32_MAX;
  memset(&plan32, 0, sizeof(plan32));
  plan32.dim = (unsigned long int)UINT32_MAX;
  plan32.local_dim = (unsigned long int)UINT32_MAX;
  plan32.nnz = 1U;
  plan32.col_index = column32;
  assert_int_eq(RemapSymmetryMatvecPlanColumns(&plan32), 0, label);
  assert_ulong_eq(
      (unsigned long int)plan32.column_slot_width,
      (unsigned long int)SYMMETRY_COLUMN_U32, label);
  assert_ulong_eq(
      (unsigned long int)symmetry_plan_column_slot(&plan32, 0U),
      (unsigned long int)UINT32_MAX - 1UL, label);
  free(plan32.column_slot32);
  free(plan32.column_slot64);

#if ULONG_MAX > UINT32_MAX
  {
    struct SymmetryMatvecPlan plan64;
    unsigned long int *column64 =
        (unsigned long int *)malloc(sizeof(*column64));
    unsigned long int count64 = (unsigned long int)UINT32_MAX + 2UL;
    if (column64 == NULL) {
      fprintf(stderr, "%s: 64-bit boundary allocation failed\n", label);
      exit(1);
    }
    *column64 = count64;
    memset(&plan64, 0, sizeof(plan64));
    plan64.dim = count64;
    plan64.local_dim = count64;
    plan64.nnz = 1U;
    plan64.col_index = column64;
    assert_int_eq(RemapSymmetryMatvecPlanColumns(&plan64), 0, label);
    assert_ulong_eq(
        (unsigned long int)plan64.column_slot_width,
        (unsigned long int)SYMMETRY_COLUMN_U64, label);
    assert_ulong_eq(
        (unsigned long int)symmetry_plan_column_slot(&plan64, 0U),
        (unsigned long int)UINT32_MAX + 1UL, label);
    free(plan64.column_slot32);
    free(plan64.column_slot64);
  }
#endif
}

static void assert_u64_column_slot_apply(const char *label)
{
#if ULONG_MAX > UINT32_MAX
  struct BindStruct X;
  struct SymmetryBasisRuntime sym;
  struct SymmetryMatvecBlockView view;
  struct SymmetryMatvecPlan plan;
  size_t row_ptr[] = {0U, 1U};
  uint64_t column_slot64[] = {0U};
  double complex values[] = {2.0 - I};
  double complex ghost_value = 3.0 + 4.0 * I;
  double complex input[] = {0.0, 1.0 + 2.0 * I};
  double complex output[] = {0.0, 0.0};
  double complex prdct = 0.0;
  double complex expected = values[0] * input[1];
  memset(&X, 0, sizeof(X));
  memset(&sym, 0, sizeof(sym));
  memset(&plan, 0, sizeof(plan));
  X.Sym = &sym;
  sym.dim = (unsigned long int)UINT32_MAX + 1UL;
  sym.local_dim = 1UL;
  sym.matvec_plan = &plan;
  plan.ready = TRUE;
  plan.block_count = 1U;
  plan.columns_remapped = TRUE;
  plan.dim = sym.dim;
  plan.local_dim = sym.local_dim;
  plan.nnz = 1U;
  plan.row_ptr = row_ptr;
  plan.column_slot_width = SYMMETRY_COLUMN_U64;
  plan.column_slot64 = column_slot64;
  plan.values = values;
  plan.halo.ready = TRUE;
  plan.halo.ghost_count = (size_t)UINT32_MAX;
  plan.halo.ghost_values = &ghost_value;
  assert_int_eq(
      SymmetryMatvecPlanGetBlockView(&plan, 0U, &view), 0, label);
  assert_int_eq(view.global_columns == NULL, 1, label);
  assert_int_eq(view.column_slot32 == NULL, 1, label);
  assert_int_eq(view.column_slot64 == column_slot64, 1, label);
  assert_int_eq(view.values == values, 1, label);
  assert_int_eq(
      ApplySymmetryMatvecPlanHalo(&X, output, input, &prdct), 0, label);
  assert_complex_close(output[1], expected, 1.0e-12, label);
  assert_complex_close(prdct, conj(input[1]) * expected, 1.0e-12, label);

  column_slot64[0] = 1U;
  output[1] = 0.0;
  prdct = 0.0;
  expected = values[0] * ghost_value;
  assert_int_eq(
      ApplySymmetryMatvecPlanHalo(&X, output, input, &prdct), 0, label);
  assert_complex_close(output[1], expected, 1.0e-12, label);
  assert_complex_close(prdct, conj(input[1]) * expected, 1.0e-12, label);

  column_slot64[0] = (uint64_t)UINT32_MAX + 1ULL;
  output[1] = 0.0;
  prdct = 0.0;
  assert_int_eq(
      ApplySymmetryMatvecPlanHalo(&X, output, input, &prdct), -1,
      "64-bit column apply rejects an out-of-range slot");
#else
  (void)label;
#endif
}

static void assert_zero_row_plan(const char *label)
{
  struct BindStruct X;
  struct SymmetryMatvecBlockView view;
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
  if (setenv("HPHI_SYMMETRY_VECTOR_EXCHANGE", "allgather", 1) != 0 ||
      ActivateSymmetryBasisDimension(&X) != 0 ||
      BuildSymmetryMatvecPlan(&X) != 0) {
    fprintf(stderr, "%s: zero-row plan setup failed\n", label);
    exit(1);
  }
  unsetenv("HPHI_SYMMETRY_VECTOR_EXCHANGE");
  assert_ulong_eq(X.Sym->local_dim, 0UL, label);
  assert_int_eq(X.Sym->matvec_plan != NULL, 1, label);
  assert_ulong_eq(
      (unsigned long int)SymmetryMatvecPlanBlockCount(
          X.Sym->matvec_plan),
      1UL, label);
  assert_int_eq(
      SymmetryMatvecPlanGetBlockView(X.Sym->matvec_plan, 0U, &view),
      0, label);
  assert_ulong_eq(view.local_row_begin, 0UL, label);
  assert_ulong_eq(view.local_row_count, 0UL, label);
  assert_ulong_eq((unsigned long int)view.nnz, 0UL, label);
  assert_int_eq(view.row_ptr == X.Sym->matvec_plan->row_ptr, 1, label);
  assert_int_eq(view.row_ptr != NULL, 1, label);
  assert_ulong_eq((unsigned long int)view.row_ptr[0], 0UL, label);
  assert_ulong_eq((unsigned long int)X.Sym->matvec_plan->nnz, 0UL, label);
  assert_ulong_eq((unsigned long int)X.Sym->matvec_plan->row_nnz_max, 0UL, label);
  assert_ulong_eq(
      (unsigned long int)X.Sym->matvec_plan->remote_column_nnz, 0UL, label);
  assert_int_eq(X.Sym->matvec_plan->halo.request_layout_ready, TRUE, label);
  assert_int_eq(X.Sym->matvec_plan->halo.ready, FALSE, label);
  assert_ulong_eq(
      (unsigned long int)X.Sym->matvec_plan->halo.ghost_count, 0UL, label);
  assert_ulong_eq(
      (unsigned long int)X.Sym->matvec_plan->column_slot_width, 32UL, label);
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

static void assert_basis_ownership_accessors(const char *label)
{
  struct BindStruct X;
  struct SymmetryBasisRuntime sym;
  struct SymmetryBasisVector basis[5];
  const struct SymmetryBasisVector *entry;
  double raw_diagonal[3] = {0.0, -41.25, 73.5};
  double *saved_list_diagonal = list_Diagonal;
  double diagonal = 0.0;

  memset(&X, 0, sizeof(X));
  memset(&sym, 0, sizeof(sym));
  memset(basis, 0, sizeof(basis));
  sym.enabled = TRUE;
  sym.dim = 4UL;
  sym.capacity = 4UL;
  sym.local_offset = 1UL;
  sym.local_dim = 2UL;
  sym.basis = basis;
  basis[1].rep_state = 0x1UL;
  basis[2].rep_state = 0x2UL;
  basis[2].diagonal = 17.25;
  basis[3].rep_state = 0x4UL;
  basis[3].diagonal = -9.5;
  basis[4].rep_state = 0x8UL;

  entry = SymmetryBasisLocalEntry(&sym, 1UL);
  assert_int_eq(entry == &basis[2], 1, label);
  entry = SymmetryBasisLocalEntry(&sym, 2UL);
  assert_int_eq(entry == &basis[3], 1, label);
  assert_int_eq(SymmetryBasisLocalEntry(&sym, 0UL) == NULL, 1, label);
  assert_int_eq(SymmetryBasisLocalEntry(&sym, 3UL) == NULL, 1, label);
  assert_int_eq(SymmetryBasisLocalEntry(NULL, 1UL) == NULL, 1, label);

  entry = SymmetryBasisReplicatedGlobalEntry(&sym, 1UL);
  assert_int_eq(entry == &basis[1], 1, label);
  entry = SymmetryBasisReplicatedGlobalEntry(&sym, 4UL);
  assert_int_eq(entry == &basis[4], 1, label);
  assert_int_eq(
      SymmetryBasisReplicatedGlobalEntry(&sym, 0UL) == NULL, 1, label);
  assert_int_eq(
      SymmetryBasisReplicatedGlobalEntry(&sym, 5UL) == NULL, 1, label);
  sym.capacity = 1UL;
  assert_int_eq(SymmetryBasisLocalEntry(&sym, 1UL) == NULL, 1, label);
  assert_int_eq(
      SymmetryBasisReplicatedGlobalEntry(&sym, 2UL) == NULL, 1, label);
  sym.capacity = 4UL;

  sym.local_dim = 0UL;
  assert_int_eq(SymmetryBasisLocalEntry(&sym, 1UL) == NULL, 1, label);
  sym.local_dim = 2UL;
  sym.local_offset = 4UL;
  assert_int_eq(SymmetryBasisLocalEntry(&sym, 1UL) == NULL, 1, label);
  sym.local_offset = 1UL;
  sym.basis = NULL;
  assert_int_eq(SymmetryBasisLocalEntry(&sym, 1UL) == NULL, 1, label);
  assert_int_eq(
      SymmetryBasisReplicatedGlobalEntry(&sym, 1UL) == NULL, 1, label);
  sym.basis = basis;
  sym.enabled = FALSE;
  assert_int_eq(SymmetryBasisLocalEntry(&sym, 1UL) == NULL, 1, label);
  assert_int_eq(
      SymmetryBasisReplicatedGlobalEntry(&sym, 1UL) == NULL, 1, label);
  sym.enabled = TRUE;

  list_Diagonal = raw_diagonal;
  X.Def.iFlgSymmetryBasis = TRUE;
  X.Check.idim_max = sym.local_dim;
  X.Sym = &sym;
  assert_int_eq(GetOwnedHamiltonianDiagonal(&X, 1UL, &diagonal), 0,
                label);
  assert_complex_close(diagonal, basis[2].diagonal, 0.0, label);
  assert_int_eq(GetOwnedHamiltonianDiagonal(&X, 2UL, &diagonal), 0,
                label);
  assert_complex_close(diagonal, basis[3].diagonal, 0.0, label);
  assert_int_eq(GetOwnedHamiltonianDiagonal(&X, 0UL, &diagonal), -1,
                label);
  assert_int_eq(GetOwnedHamiltonianDiagonal(&X, 3UL, &diagonal), -1,
                label);
  assert_int_eq(GetOwnedHamiltonianDiagonal(&X, 1UL, NULL), -1, label);

  X.Def.iFlgSymmetryBasis = FALSE;
  X.Sym = NULL;
  X.Check.idim_max = 2UL;
  assert_int_eq(GetOwnedHamiltonianDiagonal(&X, 1UL, &diagonal), 0,
                label);
  assert_complex_close(diagonal, raw_diagonal[1], 0.0, label);
  assert_int_eq(GetOwnedHamiltonianDiagonal(&X, 2UL, &diagonal), 0,
                label);
  assert_complex_close(diagonal, raw_diagonal[2], 0.0, label);
  assert_int_eq(GetOwnedHamiltonianDiagonal(&X, 3UL, &diagonal), -1,
                label);
  assert_int_eq(GetOwnedHamiltonianDiagonal(NULL, 1UL, &diagonal), -1,
                label);
  list_Diagonal = NULL;
  assert_int_eq(GetOwnedHamiltonianDiagonal(&X, 1UL, &diagonal), -1,
                label);
  list_Diagonal = saved_list_diagonal;
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
  assert_int_eq(X.Sym->basis_transform_calls <
                    X.Sym->basis_raw_states * X.Def.NSymTrans +
                    X.Sym->basis_representative_candidates * X.Def.NSymTrans,
                1, label);
  assert_ulong_eq((unsigned long int)X.Sym->basis_raw_states,
                  test_raw_dim, label);
  assert_int_eq(X.Sym->basis_thread_count >= 1U, 1, label);
  assert_int_eq(X.Sym->basis_thread_transform_calls_max <=
                    X.Sym->basis_transform_calls,
                1, label);
  assert_ulong_eq(X.Sym->dim, expected_dim, label);
  assert_ulong_eq((unsigned long int)X.Sym->basis_compatible_survivors,
                  expected_dim, label);
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

static void discard_raw_basis_storage(void)
{
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

static void assert_streaming_basis_without_raw_lists(void)
{
  struct BindStruct X;

  setup_bind(&X, 6U, 3U, 1U);
  set_ising_ring_diagonal(&X.Def, 6U, 0.37);
  discard_raw_basis_storage();
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "Spin streaming basis build failed\n");
    exit(1);
  }
  assert_ulong_eq(X.Sym->dim, 3UL,
                  "Spin basis builds without list_1/list_Diagonal");
  FreeSymmetryBasis(X.Sym);

  setup_spinless_bind(&X, 4U, 2U, 1U);
  setup_spinless_coulomb_ring(&X.Def, 4U, 0.25);
  discard_raw_basis_storage();
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "Spinless streaming basis build failed\n");
    exit(1);
  }
  assert_ulong_eq(X.Sym->dim, 2UL,
                  "Spinless basis builds without list_1/list_Diagonal");
  FreeSymmetryBasis(X.Sym);

  setup_hubbard_bind(&X, 4U, 1U, 1U, 0U);
  setup_hubbard_coulomb_intra(&X.Def, 4U, 0.5);
  discard_raw_basis_storage();
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "Hubbard streaming basis build failed\n");
    exit(1);
  }
  assert_ulong_eq(X.Sym->dim, 4UL,
                  "Hubbard basis builds without list_1/list_Diagonal");
  FreeSymmetryBasis(X.Sym);
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
  if (diagonal_coupling != 0.0)
    set_ising_ring_diagonal(&X.Def, nsite, diagonal_coupling);
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
  set_ising_ring_diagonal(&X.Def, nsite, diagonal_coupling);
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: BuildSymmetryBasis failed\n", label);
    exit(1);
  }
  for (beta = 1; beta <= X.Sym->dim; beta++) {
    unsigned int g;
    double representative_diagonal = X.Sym->basis[beta].diagonal;
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

static int c1_basis_vector_fields_equal(
    const struct SymmetryBasisVector *lhs,
    const struct SymmetryBasisVector *rhs)
{
  return lhs->rep_state == rhs->rep_state &&
      lhs->orbit_size == rhs->orbit_size &&
      lhs->stabilizer_size == rhs->stabilizer_size &&
      memcmp(&lhs->norm, &rhs->norm, sizeof(lhs->norm)) == 0 &&
      memcmp(&lhs->stabilizer_character_sum,
             &rhs->stabilizer_character_sum,
             sizeof(lhs->stabilizer_character_sum)) == 0 &&
      memcmp(&lhs->diagonal, &rhs->diagonal,
             sizeof(lhs->diagonal)) == 0;
}

static void assert_rank_local_basis_run_contract(const char *label)
{
  struct BindStruct X;
  struct SymmetryBasisRuntime local_sym;
  struct SymmetryBasisRun run = {NULL, 0UL, 0UL};
  struct SymmetryBasisRun invalid_run;
  struct SymmetryBasisVector *invalid_entries;
  unsigned long int *presence;
  unsigned long int global_count;
  unsigned long long global_raw_states;
  unsigned long int index;
  int empty_rank_count;

  setup_spinless_bind(&X, 4U, 2U, 1U);
  memset(&local_sym, 0, sizeof(local_sym));
  invalid_entries = (struct SymmetryBasisVector *)calloc(
      1U, sizeof(*invalid_entries));
  if (invalid_entries == NULL) {
    fprintf(stderr, "%s: invalid-run fixture allocation failed\n", label);
    exit(1);
  }
  invalid_run.entries = invalid_entries;
  invalid_run.count = 0UL;
  invalid_run.capacity = 0UL;
  assert_int_eq(
      BuildRankLocalSymmetryBasisRun(&X, &local_sym, &invalid_run),
      -1, "invalid rank-local run is rejected");
  assert_int_eq(invalid_run.entries == invalid_entries &&
                    invalid_run.count == 0UL &&
                    invalid_run.capacity == 0UL,
                1, "failed rank-local build preserves caller ownership");
  free(invalid_entries);

  if (BuildRankLocalSymmetryBasisRun(&X, &local_sym, &run) != 0) {
    fprintf(stderr, "%s: rank-local run build failed\n", label);
    exit(1);
  }
  assert_int_eq(run.entries != NULL, 1, label);
  assert_ulong_eq(run.capacity, run.count + 1UL, label);
  assert_ulong_eq(run.entries[0].rep_state, 0UL, label);
  assert_int_eq(run.entries[0].orbit_size == 0U &&
                    run.entries[0].stabilizer_size == 0U,
                1, label);
  assert_complex_close(run.entries[0].norm, 0.0, 0.0, label);
  assert_complex_close(run.entries[0].stabilizer_character_sum,
                       0.0, 0.0, label);
  assert_complex_close(run.entries[0].diagonal, 0.0, 0.0, label);

  global_count = run.count;
  global_raw_states = local_sym.basis_raw_states;
  empty_rank_count = run.count == 0UL ? 1 : 0;
#ifdef MPI
  if (nproc > 1) {
    unsigned long int reduced_count = 0UL;
    unsigned long long reduced_raw_states = 0ULL;
    int reduced_empty_ranks = 0;
    if (MPI_Allreduce(&global_count, &reduced_count, 1, MPI_UNSIGNED_LONG,
                      MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&global_raw_states, &reduced_raw_states, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&empty_rank_count, &reduced_empty_ranks, 1, MPI_INT,
                      MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
      fprintf(stderr, "%s: rank-local count reduction failed\n", label);
      exit(1);
    }
    global_count = reduced_count;
    global_raw_states = reduced_raw_states;
    empty_rank_count = reduced_empty_ranks;
  }
#endif
  assert_ulong_eq((unsigned long int)global_raw_states,
                  X.Check.idim_max, label);

  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: replicated compatibility build failed\n", label);
    exit(1);
  }
  assert_ulong_eq(global_count, X.Sym->dim, label);
  assert_int_eq(local_sym.basis_raw_states == X.Sym->basis_raw_states &&
                    local_sym.basis_representative_candidates ==
                        X.Sym->basis_representative_candidates &&
                    local_sym.basis_compatible_survivors ==
                        X.Sym->basis_compatible_survivors &&
                    local_sym.basis_transform_calls ==
                        X.Sym->basis_transform_calls &&
                    local_sym.basis_state_enumerator_calls ==
                        X.Sym->basis_state_enumerator_calls &&
                    local_sym.basis_diagonal_evaluator_calls ==
                        X.Sym->basis_diagonal_evaluator_calls,
                1, label);
  if (nproc > 1 && (unsigned long int)nproc > X.Sym->dim) {
    assert_int_eq(
        empty_rank_count >= nproc - (int)X.Sym->dim, 1, label);
  }

  presence = (unsigned long int *)calloc(
      (size_t)X.Sym->dim + 1U, sizeof(*presence));
  if (presence == NULL) {
    fprintf(stderr, "%s: presence allocation failed\n", label);
    exit(1);
  }
  for (index = 1UL; index <= run.count; index++) {
    unsigned long int beta;
    int found = FALSE;
    for (beta = 1UL; beta <= X.Sym->dim; beta++) {
      if (run.entries[index].rep_state !=
          X.Sym->basis[beta].rep_state) {
        continue;
      }
      assert_int_eq(c1_basis_vector_fields_equal(
                        &run.entries[index], &X.Sym->basis[beta]),
                    1, label);
      presence[beta]++;
      found = TRUE;
      break;
    }
    assert_int_eq(found, TRUE, label);
  }
#ifdef MPI
  if (nproc > 1) {
    if (X.Sym->dim + 1UL > (unsigned long int)INT_MAX ||
        MPI_Allreduce(MPI_IN_PLACE, presence, (int)(X.Sym->dim + 1UL),
                      MPI_UNSIGNED_LONG, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      fprintf(stderr, "%s: rank-local presence reduction failed\n", label);
      exit(1);
    }
  }
#endif
  for (index = 1UL; index <= X.Sym->dim; index++) {
    assert_ulong_eq(presence[index], 1UL, label);
  }

  free(presence);
  FreeSymmetryBasisRun(&run);
  assert_int_eq(run.entries == NULL && run.count == 0UL &&
                    run.capacity == 0UL,
                1, label);
  FreeSymmetryBasis(X.Sym);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

static void assert_empty_rank_local_basis_run(const char *label)
{
  struct BindStruct X;
  struct SymmetryBasisRuntime local_sym;
  struct SymmetryBasisRun run = {NULL, 0UL, 0UL};
  unsigned long int global_count;

  setup_bind(&X, 4U, 0U, 1U);
  memset(&local_sym, 0, sizeof(local_sym));
  if (BuildRankLocalSymmetryBasisRun(&X, &local_sym, &run) != 0) {
    fprintf(stderr, "%s: empty rank-local run build failed\n", label);
    exit(1);
  }
  assert_ulong_eq(run.count, 0UL, label);
  assert_ulong_eq(run.capacity, 1UL, label);
  assert_int_eq(run.entries != NULL, 1, label);
  global_count = run.count;
#ifdef MPI
  if (nproc > 1) {
    unsigned long int reduced_count = 1UL;
    if (MPI_Allreduce(&global_count, &reduced_count, 1, MPI_UNSIGNED_LONG,
                      MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
      fprintf(stderr, "%s: empty-run count reduction failed\n", label);
      exit(1);
    }
    global_count = reduced_count;
  }
#endif
  assert_ulong_eq(global_count, 0UL, label);
  assert_int_eq(BuildSymmetryBasis(&X), -1,
                "replicated wrapper rejects zero-dimensional sector");
  assert_int_eq(X.Sym == NULL, 1, label);
  FreeSymmetryBasisRun(&run);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

int main(int argc, char **argv)
{
#ifdef MPI
  if (argc == 2 && strcmp(argv[1], "--mpi-rank-local-run") == 0) {
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS ||
        MPI_Comm_size(MPI_COMM_WORLD, &nproc) != MPI_SUCCESS ||
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank) != MPI_SUCCESS) {
      fprintf(stderr, "MPI rank-local run test initialization failed\n");
      return 1;
    }
    stdoutMPI = stderr;
    assert_rank_local_basis_run_contract(
        "rank-local run union matches replicated compatibility basis");
    assert_empty_rank_local_basis_run(
        "rank-local distribution accepts an empty global sector");
    if (myrank == 0) {
      fprintf(stdout,
              "rank-local symmetry basis run gate: PASS (%d MPI ranks)\n",
              nproc);
    }
    if (MPI_Finalize() != MPI_SUCCESS) return 1;
    return 0;
  } else if (argc == 2 && strcmp(argv[1], "--mpi-state-enumerator") == 0) {
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS ||
        MPI_Comm_size(MPI_COMM_WORLD, &nproc) != MPI_SUCCESS ||
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank) != MPI_SUCCESS) {
      fprintf(stderr, "MPI state-enumerator test initialization failed\n");
      return 1;
    }
    stdoutMPI = stderr;
    assert_state_enumerator_exact();
    assert_state_diagonal_exact();
    if (myrank == 0) {
      fprintf(stdout,
              "allocation-free state primitive gate: PASS (%d MPI ranks)\n",
              nproc);
    }
    if (MPI_Finalize() != MPI_SUCCESS) return 1;
    return 0;
  } else if (argc == 2 && strcmp(argv[1], "--mpi-column-spans") == 0) {
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS ||
        MPI_Comm_size(MPI_COMM_WORLD, &nproc) != MPI_SUCCESS ||
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank) != MPI_SUCCESS) {
      fprintf(stderr, "MPI column-span test initialization failed\n");
      return 1;
    }
    assert_mpi_column_spans_exact(
        "single-span and multi-span MPI halo plans are exact");
    if (myrank == 0) {
      fprintf(stdout,
              "single/multi-span halo exact gate: PASS (%d MPI ranks)\n",
              nproc);
    }
    if (MPI_Finalize() != MPI_SUCCESS) return 1;
    return 0;
  }
#else
  (void)argc;
  (void)argv;
#endif
  int shift4[4] = {1, 2, 3, 0};
  stdoutMPI = stderr;
  assert_state_enumerator_exact();
  assert_state_diagonal_exact();
  assert_exhaustive_fermion_permutation_parity();
  assert_fermion_parity_word_boundary();
  assert_spin_permutation_word_boundary();
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
  assert_rank_local_basis_run_contract(
      "serial rank-local run matches replicated compatibility basis");
  assert_empty_rank_local_basis_run(
      "serial rank-local distribution accepts an empty global sector");
  assert_basis_ownership_accessors(
      "basis ownership accessors enforce local/global ranges");
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
  assert_streaming_basis_without_raw_lists();
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
#ifdef _OPENMP
  assert_parallel_basis_matches_serial(
      "Hubbard basis fields are identical for one and four OpenMP threads");
  assert_parallel_plan_matches_serial(
      "Hubbard plan CSR is identical for one and four OpenMP threads");
#endif
  assert_vector_owner_and_request_layout(
      "halo owner mapping and request layout are deterministic");
  assert_remote_topology_plan(
      "local-row plan reports remote-column topology without changing CSR");
  assert_mixed_column_remap(
      "mixed local/remote columns remap to bounded local/ghost slots");
  assert_column_slot_width_boundaries(
      "column slot storage switches at the UINT32_MAX boundary");
  assert_u64_column_slot_apply(
      "64-bit column slot apply preserves values and bounds checks");
  assert_zero_row_plan("local-row plan supports zero-row rank");
  assert_representative_hash_matches_basis(6, 3, 1,
                                           "C6 k=pi/3 representative hash matches basis");
  assert_hash_probe_lookup_handles_collision("representative hash probing handles collisions");
  return 0;
}
