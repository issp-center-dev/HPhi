#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "DefCommon.h"
#include "mltplySpinSym.h"
#include "symmetry_basis.h"
#include "symmetry_diagonal.h"
#include "symmetry_directory.h"
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

static int perm_storage[8][8];
static int anti_storage[8][8];
static int *perm_rows[8];
static int *anti_rows[8];
static double complex chars_storage[8];
static int exchange_storage[8][2];
static int *exchange_rows[8];
static double exchange_params[8];
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

static int representative_batch_stats_equal(
    const struct SymmetryRepresentativeBatchStats *left,
    const struct SymmetryRepresentativeBatchStats *right)
{
  /* Field-wise comparison avoids reading implementation-defined padding. */
  return left->directory_batch_calls ==
             right->directory_batch_calls &&
      left->directory_request_entries_sent ==
          right->directory_request_entries_sent &&
      left->directory_request_entries_received ==
          right->directory_request_entries_received &&
      left->directory_found_entries ==
          right->directory_found_entries &&
      left->directory_not_found_entries ==
          right->directory_not_found_entries &&
      left->directory_lookup_probe_count ==
          right->directory_lookup_probe_count &&
      left->directory_lookup_max_probe ==
          right->directory_lookup_max_probe &&
      left->directory_owner_peer_count_max ==
          right->directory_owner_peer_count_max &&
      left->directory_requester_peer_count_max ==
          right->directory_requester_peer_count_max &&
      left->directory_exchange_used_chunked ==
          right->directory_exchange_used_chunked &&
      left->directory_exchange_message_byte_limit ==
          right->directory_exchange_message_byte_limit &&
      left->directory_exchange_max_message_bytes ==
          right->directory_exchange_max_message_bytes &&
      left->directory_exchange_send_messages ==
          right->directory_exchange_send_messages &&
      left->directory_exchange_recv_messages ==
          right->directory_exchange_recv_messages &&
      left->directory_batch_temporary_peak_bytes ==
          right->directory_batch_temporary_peak_bytes &&
      left->directory_batch_memory_warning_byte_threshold ==
          right->directory_batch_memory_warning_byte_threshold &&
      left->directory_batch_memory_byte_limit ==
          right->directory_batch_memory_byte_limit;
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
  struct SymmetryMatvecBlock flat_block;
  struct SymmetryMatvecBlock split_block;
  struct SymmetryMatvecPlan flat_plan;
  struct SymmetryMatvecPlan split_plan;
  struct SymmetryBasisRuntime sym;
  struct BindStruct X;
  struct SymmetryGlobalColumnSpan flat_span;
  struct SymmetryGlobalColumnSpan split_spans[7];
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
  memset(&flat_block, 0, sizeof(flat_block));
  memset(&split_block, 0, sizeof(split_block));
  memset(&flat_plan, 0, sizeof(flat_plan));
  memset(&split_plan, 0, sizeof(split_plan));
  memset(&sym, 0, sizeof(sym));
  memset(&X, 0, sizeof(X));
  assert_int_eq(
      SymmetryBlockRange(
          dim, myrank, nproc, &local_offset, &local_dim),
      0, label);
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
  flat_plan.blocks = &flat_block;
  flat_plan.nnz = nnz;
  flat_block.local_row_count = local_dim;
  flat_block.nnz = nnz;
  flat_block.row_ptr = row_ptr;
  flat_block.global_columns = flat_columns;
  flat_block.values = values;
  flat_plan.halo = flat_halo;
  split_plan = flat_plan;
  split_plan.blocks = &split_block;
  split_block = flat_block;
  split_block.global_columns = split_columns;
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

  free(flat_block.column_slot32);
  free(flat_block.column_slot64);
  free(split_block.column_slot32);
  free(split_block.column_slot64);
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

static void assert_vector_owner_inverse_contract(const char *label)
{
  unsigned long int dim;
  int nrank;
  for (nrank = 1; nrank <= 24; nrank++) {
    for (dim = 0UL; dim <= 200UL; dim++) {
      int rank;
      unsigned long int covered = 0UL;
      for (rank = 0; rank < nrank; rank++) {
        unsigned long int offset;
        unsigned long int count;
        unsigned long int local_index;
        assert_int_eq(
            SymmetryBlockRange(dim, rank, nrank, &offset, &count),
            0, label);
        assert_ulong_eq(offset, covered, label);
        for (local_index = 0UL; local_index < count; local_index++) {
          assert_int_eq(
              SymmetryVectorOwnerOfGlobalIndex(
                  dim, nrank, offset + local_index + 1UL),
              rank, label);
        }
        covered = offset + count;
      }
      assert_ulong_eq(covered, dim, label);
      assert_int_eq(
          SymmetryVectorOwnerOfGlobalIndex(dim, nrank, 0UL),
          -1, label);
      if (dim < ULONG_MAX) {
        assert_int_eq(
            SymmetryVectorOwnerOfGlobalIndex(dim, nrank, dim + 1UL),
            -1, label);
      }
    }
  }
}

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
  assert_vector_owner_inverse_contract(label);
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
  struct SymmetryMatvecBlock block;
  struct SymmetryMatvecBlock invalid_block;
  struct SymmetryMatvecPlan plan;
  struct SymmetryMatvecPlan invalid_plan;
  size_t index;
  columns = (unsigned long int *)malloc(sizeof(initial_columns));
  if (columns == NULL) {
    fprintf(stderr, "%s: column allocation failed\n", label);
    exit(1);
  }
  memcpy(columns, initial_columns, sizeof(initial_columns));
  memset(&block, 0, sizeof(block));
  memset(&plan, 0, sizeof(plan));
  plan.dim = 10UL;
  plan.local_offset = 4UL;
  plan.local_dim = 3UL;
  plan.block_count = 1U;
  plan.blocks = &block;
  plan.nnz = sizeof(initial_columns) / sizeof(initial_columns[0]);
  block.local_row_count = plan.local_dim;
  block.nnz = plan.nnz;
  block.global_columns = columns;
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

  memset(&invalid_block, 0, sizeof(invalid_block));
  memset(&invalid_plan, 0, sizeof(invalid_plan));
  invalid_plan.dim = plan.dim;
  invalid_plan.local_offset = plan.local_offset;
  invalid_plan.local_dim = plan.local_dim;
  invalid_plan.block_count = 1U;
  invalid_plan.blocks = &invalid_block;
  invalid_plan.nnz = 1U;
  invalid_block.local_row_count = invalid_plan.local_dim;
  invalid_block.nnz = invalid_plan.nnz;
  invalid_block.global_columns = missing_ghost_column;
  invalid_plan.halo.ghost_count = plan.halo.ghost_count;
  invalid_plan.halo.ghost_global_index = ghosts;
  assert_int_eq(RemapSymmetryMatvecPlanColumns(&invalid_plan), -1,
                "column remap rejects a missing ghost index");
  free(block.column_slot32);
  free(block.column_slot64);
}

static void assert_column_slot_width_boundaries(const char *label)
{
  struct SymmetryMatvecBlock block32;
  struct SymmetryMatvecPlan plan32;
  unsigned long int *column32 =
      (unsigned long int *)malloc(sizeof(*column32));
  if (column32 == NULL) {
    fprintf(stderr, "%s: 32-bit boundary allocation failed\n", label);
    exit(1);
  }
  *column32 = (unsigned long int)UINT32_MAX;
  memset(&block32, 0, sizeof(block32));
  memset(&plan32, 0, sizeof(plan32));
  plan32.dim = (unsigned long int)UINT32_MAX;
  plan32.local_dim = (unsigned long int)UINT32_MAX;
  plan32.block_count = 1U;
  plan32.blocks = &block32;
  plan32.nnz = 1U;
  block32.local_row_count = plan32.local_dim;
  block32.nnz = plan32.nnz;
  block32.global_columns = column32;
  assert_int_eq(RemapSymmetryMatvecPlanColumns(&plan32), 0, label);
  assert_ulong_eq(
      (unsigned long int)plan32.column_slot_width,
      (unsigned long int)SYMMETRY_COLUMN_U32, label);
  assert_ulong_eq(
      (unsigned long int)symmetry_plan_column_slot(&plan32, 0U),
      (unsigned long int)UINT32_MAX - 1UL, label);
  free(block32.column_slot32);
  free(block32.column_slot64);

#if ULONG_MAX > UINT32_MAX
  {
    struct SymmetryMatvecBlock block64;
    struct SymmetryMatvecPlan plan64;
    unsigned long int *column64 =
        (unsigned long int *)malloc(sizeof(*column64));
    unsigned long int count64 = (unsigned long int)UINT32_MAX + 2UL;
    if (column64 == NULL) {
      fprintf(stderr, "%s: 64-bit boundary allocation failed\n", label);
      exit(1);
    }
    *column64 = count64;
    memset(&block64, 0, sizeof(block64));
    memset(&plan64, 0, sizeof(plan64));
    plan64.dim = count64;
    plan64.local_dim = count64;
    plan64.block_count = 1U;
    plan64.blocks = &block64;
    plan64.nnz = 1U;
    block64.local_row_count = plan64.local_dim;
    block64.nnz = plan64.nnz;
    block64.global_columns = column64;
    assert_int_eq(RemapSymmetryMatvecPlanColumns(&plan64), 0, label);
    assert_ulong_eq(
        (unsigned long int)plan64.column_slot_width,
        (unsigned long int)SYMMETRY_COLUMN_U64, label);
    assert_ulong_eq(
        (unsigned long int)symmetry_plan_column_slot(&plan64, 0U),
        (unsigned long int)UINT32_MAX + 1UL, label);
    free(block64.column_slot32);
    free(block64.column_slot64);
  }
#endif
}

static void assert_u64_column_slot_apply(const char *label)
{
#if ULONG_MAX > UINT32_MAX
  struct BindStruct X;
  struct SymmetryBasisRuntime sym;
  struct SymmetryMatvecBlock block;
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
  memset(&block, 0, sizeof(block));
  memset(&plan, 0, sizeof(plan));
  X.Sym = &sym;
  sym.dim = (unsigned long int)UINT32_MAX + 1UL;
  sym.local_dim = 1UL;
  sym.matvec_plan = &plan;
  plan.ready = TRUE;
  plan.block_count = 1U;
  plan.blocks = &block;
  plan.columns_remapped = TRUE;
  plan.dim = sym.dim;
  plan.local_dim = sym.local_dim;
  plan.nnz = 1U;
  block.local_row_count = plan.local_dim;
  block.nnz = plan.nnz;
  block.row_ptr = row_ptr;
  plan.column_slot_width = SYMMETRY_COLUMN_U64;
  block.column_slot64 = column_slot64;
  block.values = values;
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

static void assert_owned_multi_block_plan(const char *label)
{
  const unsigned long int columns0_data[] = {1UL, 3UL, 2UL};
  const unsigned long int columns1_data[] = {4UL, 1UL, 3UL};
  size_t row_ptr0[] = {0U, 2U, 3U};
  size_t row_ptr1[] = {0U, 1U, 3U};
  double complex values0[] = {2.0, -0.5 + 0.25 * I, 1.25};
  double complex values1[] = {0.75 - 0.5 * I, -1.0, 0.5 + 0.5 * I};
  double complex input[] = {
      0.0, 1.0 + 0.5 * I, -0.25 + 0.75 * I,
      0.5 - I, 1.5 + 0.25 * I};
  double complex global_output[5] = {0.0};
  double complex halo_output[5] = {0.0};
  double complex expected[5] = {0.0};
  double complex global_prdct = 0.0;
  double complex halo_prdct = 0.0;
  double complex expected_prdct = 0.0;
  struct SymmetryMatvecBlock blocks[2];
  struct SymmetryMatvecBlockView view;
  struct SymmetryMatvecPlan plan;
  struct SymmetryBasisRuntime sym;
  struct BindStruct X;
  unsigned long int *columns0;
  unsigned long int *columns1;
  unsigned long int row;
  size_t p;
#ifdef _OPENMP
  int saved_dynamic = omp_get_dynamic();
  int saved_threads = omp_get_max_threads();
  omp_set_dynamic(0);
  omp_set_num_threads(1);
#endif

  columns0 =
      (unsigned long int *)malloc(sizeof(columns0_data));
  columns1 =
      (unsigned long int *)malloc(sizeof(columns1_data));
  if (columns0 == NULL || columns1 == NULL) {
    fprintf(stderr, "%s: multi-block column allocation failed\n", label);
    exit(1);
  }
  memcpy(columns0, columns0_data, sizeof(columns0_data));
  memcpy(columns1, columns1_data, sizeof(columns1_data));
  memset(blocks, 0, sizeof(blocks));
  memset(&plan, 0, sizeof(plan));
  memset(&sym, 0, sizeof(sym));
  memset(&X, 0, sizeof(X));
  blocks[0].local_row_count = 2UL;
  blocks[0].nnz = 3U;
  blocks[0].row_ptr = row_ptr0;
  blocks[0].global_columns = columns0;
  blocks[0].values = values0;
  blocks[1].local_row_begin = 2UL;
  blocks[1].local_row_count = 2UL;
  blocks[1].nnz = 3U;
  blocks[1].row_ptr = row_ptr1;
  blocks[1].global_columns = columns1;
  blocks[1].values = values1;
  plan.ready = TRUE;
  plan.dim = 4UL;
  plan.local_dim = 4UL;
  plan.block_count = 2U;
  plan.blocks = blocks;
  plan.nnz = 6U;
  plan.halo.ready = TRUE;
  sym.enabled = TRUE;
  sym.dim = plan.dim;
  sym.local_dim = plan.local_dim;
  sym.matvec_plan = &plan;
  X.Sym = &sym;

  assert_ulong_eq(
      (unsigned long int)SymmetryMatvecPlanBlockCount(&plan), 2UL, label);
  assert_int_eq(
      SymmetryMatvecPlanGetBlockView(&plan, 0U, &view), 0, label);
  assert_ulong_eq(view.local_row_begin, 0UL, label);
  assert_ulong_eq(view.local_row_count, 2UL, label);
  assert_int_eq(view.row_ptr == row_ptr0, 1, label);
  assert_int_eq(
      SymmetryMatvecPlanGetBlockView(&plan, 1U, &view), 0, label);
  assert_ulong_eq(view.local_row_begin, 2UL, label);
  assert_ulong_eq(view.local_row_count, 2UL, label);
  assert_int_eq(view.row_ptr == row_ptr1, 1, label);
  assert_int_eq(
      SymmetryMatvecPlanGetBlockView(&plan, 2U, &view), -1, label);

  for (row = 0UL; row < plan.local_dim; row++) {
    const struct SymmetryMatvecBlock *block =
        row < 2UL ? &blocks[0] : &blocks[1];
    unsigned long int block_row = row - block->local_row_begin;
    double complex sum = 0.0;
    for (p = block->row_ptr[block_row];
         p < block->row_ptr[block_row + 1UL]; p++) {
      sum += block->values[p] * input[block->global_columns[p]];
    }
    expected[row + 1UL] = sum;
    expected_prdct += conj(input[row + 1UL]) * sum;
  }
  assert_int_eq(
      ApplySymmetryMatvecPlan(
          &X, global_output, input, &global_prdct),
      0, label);
  for (row = 1UL; row <= plan.local_dim; row++) {
    assert_complex_bitwise(global_output[row], expected[row], label);
  }
  assert_complex_bitwise(global_prdct, expected_prdct, label);

  blocks[1].local_row_begin = 3UL;
  assert_int_eq(
      RemapSymmetryMatvecPlanColumns(&plan), -1,
      "multi-block remap rejects a row gap");
  blocks[1].local_row_begin = 2UL;
  columns1[0] = 5UL;
  assert_int_eq(
      RemapSymmetryMatvecPlanColumns(&plan), -1,
      "multi-block remap failure is atomic across blocks");
  assert_int_eq(
      blocks[0].column_slot32 == NULL &&
          blocks[0].column_slot64 == NULL &&
          blocks[1].column_slot32 == NULL &&
          blocks[1].column_slot64 == NULL &&
          blocks[0].global_columns == columns0 &&
          blocks[1].global_columns == columns1 &&
          plan.columns_remapped == FALSE,
      1, label);
  columns1[0] = 4UL;
  assert_int_eq(RemapSymmetryMatvecPlanColumns(&plan), 0, label);
  assert_int_eq(plan.row_ptr == NULL && plan.col_index == NULL &&
                    plan.column_slot32 == NULL &&
                    plan.column_slot64 == NULL && plan.values == NULL,
                1, label);
  assert_int_eq(
      ApplySymmetryMatvecPlanHalo(
          &X, halo_output, input, &halo_prdct),
      0, label);
  for (row = 1UL; row <= plan.local_dim; row++) {
    assert_complex_bitwise(halo_output[row], expected[row], label);
  }
  assert_complex_bitwise(halo_prdct, expected_prdct, label);

  free(blocks[0].column_slot32);
  free(blocks[0].column_slot64);
  free(blocks[1].column_slot32);
  free(blocks[1].column_slot64);
#ifdef _OPENMP
  omp_set_num_threads(saved_threads);
  omp_set_dynamic(saved_dynamic);
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

enum C5ReferenceModel {
  C5_REFERENCE_SPIN = 0,
  C5_REFERENCE_SPINLESS = 1,
  C5_REFERENCE_HUBBARD = 2
};

static void setup_c5_reference_bind(
    struct BindStruct *X,
    enum C5ReferenceModel model)
{
  if (model == C5_REFERENCE_SPIN) {
    setup_bind(X, 6U, 3U, 1U);
    set_ising_ring_diagonal(&X->Def, 6U, 0.37);
  } else if (model == C5_REFERENCE_SPINLESS) {
    setup_spinless_bind(X, 4U, 2U, 1U);
    setup_spinless_transfer_ring(&X->Def, 4U);
    setup_spinless_coulomb_ring(&X->Def, 4U, 0.25);
  } else {
    setup_hubbard_bind(X, 4U, 1U, 1U, 0U);
    setup_hubbard_transfer_ring(&X->Def, 4U);
    setup_hubbard_coulomb_intra(&X->Def, 4U, 0.5);
  }
}

static void assert_representative_discovery_model(
    enum C5ReferenceModel model,
    const char *label)
{
  struct BindStruct X;
  struct SymmetryBasisRuntime *saved_sym;
  struct SymmetryRepresentativeResult representative;
  struct SymmetryRepresentativeResult fallback;
  struct SymmetryRepresentativeResult zero_result;
  unsigned int saved_enabled;
  unsigned long int raw;
  unsigned long int member_count = 0UL;
  unsigned long int nonmember_count = 0UL;
  unsigned long int identity_count = 0UL;
  unsigned long int nontrivial_count = 0UL;
  unsigned long int stabilizer_count = 0UL;

  setup_c5_reference_bind(&X, model);
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: basis build failed\n", label);
    exit(1);
  }

  memset(&zero_result, 0, sizeof(zero_result));
  memset(&representative, 0xa5, sizeof(representative));
  assert_int_eq(
      SymmetryFindRepresentative(NULL, 0UL, &representative), -1, label);
  assert_int_eq(
      memcmp(&representative, &zero_result, sizeof(representative)), 0,
      label);
  assert_int_eq(
      SymmetryFindRepresentative(&X, 0UL, NULL), -1, label);

  saved_sym = X.Sym;
  X.Sym = NULL;
  memset(&representative, 0xa5, sizeof(representative));
  assert_int_eq(
      SymmetryFindRepresentative(&X, 0UL, &representative), -1, label);
  assert_int_eq(
      memcmp(&representative, &zero_result, sizeof(representative)), 0,
      label);
  X.Sym = saved_sym;

  saved_enabled = (unsigned int)X.Sym->enabled;
  X.Sym->enabled = FALSE;
  memset(&representative, 0xa5, sizeof(representative));
  assert_int_eq(
      SymmetryFindRepresentative(&X, 0UL, &representative), -1, label);
  assert_int_eq(
      memcmp(&representative, &zero_result, sizeof(representative)), 0,
      label);
  X.Sym->enabled = (int)saved_enabled;

  assert_int_eq(
      SymmetryFindRepresentative(&X, 0UL, &representative), 0, label);
  assert_ulong_eq(representative.rep_state, 0UL, label);
  assert_int_eq(representative.op_rep_to_state, 0U, label);
  assert_complex_bitwise(representative.phase, 1.0, label);
#ifdef HPHI_SYMMETRY_CANONICAL_VERIFY
  {
    unsigned int saved_inverse = X.Sym->group_inverse[0];
    X.Sym->group_inverse[0] = 1U;
    memset(&representative, 0xa5, sizeof(representative));
    assert_int_eq(
        SymmetryFindRepresentative(&X, 0UL, &representative), -1, label);
    assert_int_eq(
        memcmp(&representative, &zero_result, sizeof(representative)), 0,
        label);
    X.Sym->group_inverse[0] = saved_inverse;
  }
#endif

  for (raw = 1UL; raw <= X.Check.idim_max; raw++) {
    struct SymmetryCanonicalResult canonical;
    unsigned long int state = list_1[raw];
    unsigned long int expected_rep = state;
    unsigned long int expected_basis_index = 0UL;
    unsigned int expected_op = UINT_MAX;
    double complex expected_phase = 0.0;
    unsigned int g;
    unsigned int stabilizer_size = 0U;
    unsigned long int beta;

    for (g = 0U; g < X.Def.NSymTrans; g++) {
      struct SymmetryTransformResult moved;
      assert_int_eq(
          SymmetryApplyToState(&X.Def, state, g, &moved), 0, label);
      if (moved.state < expected_rep) expected_rep = moved.state;
      if (moved.state == state) stabilizer_size++;
    }
    for (g = 0U; g < X.Def.NSymTrans; g++) {
      struct SymmetryTransformResult moved;
      assert_int_eq(
          SymmetryApplyToState(&X.Def, expected_rep, g, &moved), 0,
          label);
      if (moved.state == state) {
        expected_op = g;
        expected_phase = X.Def.SymTransChar[g] * moved.amplitude;
        break;
      }
    }
    assert_int_eq(expected_op != UINT_MAX, 1, label);
    if (expected_op == 0U) {
      identity_count++;
    } else {
      nontrivial_count++;
    }
    if (stabilizer_size > 1U) stabilizer_count++;

    assert_int_eq(
        SymmetryFindRepresentative(&X, state, &representative), 0,
        label);
    assert_ulong_eq(representative.rep_state, expected_rep, label);
    assert_int_eq(representative.op_rep_to_state, expected_op, label);
    assert_complex_bitwise(representative.phase, expected_phase, label);

    for (beta = 1UL; beta <= X.Sym->dim; beta++) {
      if (X.Sym->basis[beta].rep_state == expected_rep) {
        expected_basis_index = beta;
        break;
      }
    }
    assert_int_eq(
        SymmetryCanonicalizeState(&X, state, &canonical), 0, label);
    if (expected_basis_index != 0UL) {
      member_count++;
      assert_int_eq(canonical.found, TRUE, label);
      assert_ulong_eq(canonical.basis_index, expected_basis_index, label);
      assert_int_eq(
          canonical.op_rep_to_state, representative.op_rep_to_state,
          label);
      assert_complex_bitwise(
          canonical.phase, representative.phase, label);
    } else {
      nonmember_count++;
      assert_int_eq(canonical.found, FALSE, label);
      assert_ulong_eq(canonical.basis_index, 0UL, label);
      assert_int_eq(canonical.op_rep_to_state, 0U, label);
      assert_complex_bitwise(canonical.phase, 0.0, label);
    }
  }
  assert_int_eq(member_count > 0UL, 1, label);
  assert_int_eq(identity_count > 0UL, 1, label);
  assert_int_eq(nontrivial_count > 0UL, 1, label);
  if (model == C5_REFERENCE_SPIN) {
    assert_int_eq(nonmember_count > 0UL, 1, label);
    assert_int_eq(stabilizer_count > 0UL, 1, label);
  }

  assert_int_eq(
      SymmetryFindRepresentative(
          &X, list_1[X.Check.idim_max], &representative),
      0, label);
  {
    unsigned int *saved_group_inverse = X.Sym->group_inverse;
    X.Sym->group_inverse = NULL;
    assert_int_eq(
        SymmetryFindRepresentative(
            &X, list_1[X.Check.idim_max], &fallback),
        0, label);
    X.Sym->group_inverse = saved_group_inverse;
  }
  assert_ulong_eq(fallback.rep_state, representative.rep_state, label);
  assert_int_eq(
      fallback.op_rep_to_state, representative.op_rep_to_state, label);
  assert_complex_bitwise(fallback.phase, representative.phase, label);

  FreeSymmetryBasis(X.Sym);
  X.Sym = NULL;
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

static int c5_noop_entry(
    unsigned long int out_index,
    double complex coefficient,
    void *context)
{
  (void)out_index;
  (void)coefficient;
  (void)context;
  return 0;
}

struct C5DirectoryTarget {
  unsigned long int rep_state;
  unsigned long int global_beta;
  double norm;
};

static int compare_c5_directory_target(const void *lhs, const void *rhs)
{
  const struct C5DirectoryTarget *left =
      (const struct C5DirectoryTarget *)lhs;
  const struct C5DirectoryTarget *right =
      (const struct C5DirectoryTarget *)rhs;
  if (left->rep_state < right->rep_state) return -1;
  if (left->rep_state > right->rep_state) return 1;
  return 0;
}

static int c5_transition(
    const struct BindStruct *X,
    enum C5ReferenceModel model,
    unsigned long int state,
    unsigned int term,
    unsigned long int *out_state,
    double complex *hval)
{
  if (X == NULL || out_state == NULL || hval == NULL) return -1;
  if (model == C5_REFERENCE_SPIN) {
    int applied = apply_exchange_halfspin_test(
        state, X->Def.ExchangeCoupling[term][0],
        X->Def.ExchangeCoupling[term][1], out_state);
    if (applied == TRUE) *hval = X->Def.ParaExchangeCoupling[term];
    return applied;
  }
  if (model == C5_REFERENCE_SPINLESS) {
    double complex trans = -X->Def.EDParaGeneralTransfer[term];
    return apply_spinless_hopping_hermite_test(
        state, (unsigned int)X->Def.EDGeneralTransfer[term][0],
        (unsigned int)X->Def.EDGeneralTransfer[term][2],
        trans, out_state, hval);
  }
  {
    double complex trans = -X->Def.EDParaGeneralTransfer[term];
    return apply_hubbard_hopping_hermite_test(
        state, (unsigned int)X->Def.EDGeneralTransfer[term][0],
        (unsigned int)X->Def.EDGeneralTransfer[term][1],
        (unsigned int)X->Def.EDGeneralTransfer[term][2],
        (unsigned int)X->Def.EDGeneralTransfer[term][3],
        trans, out_state, hval);
  }
}

static int c5_transition_state(
    const struct BindStruct *X,
    enum C5ReferenceModel model,
    unsigned long int state,
    unsigned int term,
    unsigned long int *out_state)
{
  double complex hval;
  return c5_transition(
      X, model, state, term, out_state, &hval);
}

static struct C5DirectoryTarget *collect_c5_directory_targets(
    const struct BindStruct *X,
    enum C5ReferenceModel model,
    unsigned long int *target_count,
    int *found_count,
    int *interior_miss_count,
    int *range_miss_count,
    const char *label)
{
  struct C5DirectoryTarget *targets;
  unsigned long int count = 0UL;
  unsigned long int maximum;
  unsigned long int beta;
  unsigned int term_count;
  unsigned int term_step;
  if (X == NULL || X->Sym == NULL || X->Sym->dim == 0UL ||
      target_count == NULL || found_count == NULL ||
      interior_miss_count == NULL || range_miss_count == NULL) {
    fprintf(stderr, "%s: invalid transition target collector input\n", label);
    exit(1);
  }
  if (model == C5_REFERENCE_SPIN) {
    term_count = X->Def.NExchangeCoupling;
    term_step = 1U;
  } else {
    term_count = X->Def.EDNTransfer;
    term_step = 2U;
  }
  if (term_count == 0U ||
      X->Sym->dim > ULONG_MAX / (unsigned long int)term_count) {
    fprintf(stderr, "%s: transition target capacity overflow\n", label);
    exit(1);
  }
  maximum = X->Sym->dim * (unsigned long int)term_count;
  if (maximum > (unsigned long int)(SIZE_MAX / sizeof(*targets))) {
    fprintf(stderr, "%s: transition target byte overflow\n", label);
    exit(1);
  }
  targets = (struct C5DirectoryTarget *)malloc(
      (size_t)maximum * sizeof(*targets));
  if (targets == NULL) {
    fprintf(stderr, "%s: transition target allocation failed\n", label);
    exit(1);
  }
  for (beta = 1UL; beta <= X->Sym->dim; beta++) {
    unsigned int term;
    for (term = 0U; term < term_count; term += term_step) {
      struct SymmetryRepresentativeResult representative;
      struct SymmetryCanonicalResult canonical;
      unsigned long int out_state;
      int applied = c5_transition_state(
          X, model, X->Sym->basis[beta].rep_state, term, &out_state);
      if (applied < 0) {
        fprintf(stderr, "%s: Hamiltonian transition failed\n", label);
        exit(1);
      }
      if (applied == 0) continue;
      if (SymmetryFindRepresentative(
              X, out_state, &representative) != 0 ||
          SymmetryCanonicalizeState(X, out_state, &canonical) != 0) {
        fprintf(stderr, "%s: transition canonicalization failed\n", label);
        exit(1);
      }
      targets[count].rep_state = representative.rep_state;
      if (canonical.found != 0) {
        if (canonical.basis_index == 0UL ||
            canonical.basis_index > X->Sym->dim ||
            X->Sym->basis[canonical.basis_index].rep_state !=
                representative.rep_state) {
          fprintf(stderr, "%s: transition canonical result mismatch\n",
                  label);
          exit(1);
        }
        targets[count].global_beta = canonical.basis_index;
        targets[count].norm = X->Sym->basis[canonical.basis_index].norm;
      } else {
        targets[count].global_beta = 0UL;
        targets[count].norm = 0.0;
      }
      count++;
    }
  }
  qsort(targets, (size_t)count, sizeof(*targets),
        compare_c5_directory_target);
  {
    unsigned long int input;
    unsigned long int output = 0UL;
    for (input = 0UL; input < count; input++) {
      if (output > 0UL &&
          targets[output - 1UL].rep_state == targets[input].rep_state) {
        if (targets[output - 1UL].global_beta !=
                targets[input].global_beta ||
            memcmp(&targets[output - 1UL].norm, &targets[input].norm,
                   sizeof(targets[input].norm)) != 0) {
          fprintf(stderr, "%s: duplicate transition oracle mismatch\n",
                  label);
          exit(1);
        }
        continue;
      }
      targets[output++] = targets[input];
    }
    count = output;
  }
  *found_count = 0;
  *interior_miss_count = 0;
  *range_miss_count = 0;
  for (beta = 0UL; beta < count; beta++) {
    if (targets[beta].global_beta != 0UL) {
      (*found_count)++;
    } else if (targets[beta].rep_state >
                   X->Sym->basis[1].rep_state &&
               targets[beta].rep_state <
                   X->Sym->basis[X->Sym->dim].rep_state) {
      (*interior_miss_count)++;
    } else {
      (*range_miss_count)++;
    }
  }
  if (count == 0UL || *found_count == 0) {
    fprintf(stderr, "%s: transition oracle lacks found entries\n", label);
    exit(1);
  }
  *target_count = count;
  return targets;
}

static void assert_c5_directory_transition_batch(
    const struct BindStruct *X,
    enum C5ReferenceModel model,
    const struct C5DirectoryTarget *targets,
    unsigned long int target_count,
    const char *label)
{
  unsigned long int *keys;
  unsigned long int *global_beta;
  double *norm;
  unsigned long int index;
  if (target_count > (unsigned long int)(SIZE_MAX / sizeof(*keys)) ||
      target_count > (unsigned long int)(SIZE_MAX / sizeof(*norm))) {
    fprintf(stderr, "%s: directory output capacity overflow\n", label);
    exit(1);
  }
  keys = (unsigned long int *)malloc((size_t)target_count * sizeof(*keys));
  global_beta =
      (unsigned long int *)malloc((size_t)target_count *
                                  sizeof(*global_beta));
  norm = (double *)malloc((size_t)target_count * sizeof(*norm));
  if (keys == NULL || global_beta == NULL || norm == NULL) {
    fprintf(stderr, "%s: directory output allocation failed\n", label);
    exit(1);
  }
  for (index = 0UL; index < target_count; index++) {
    keys[index] = targets[index].rep_state;
    global_beta[index] = ULONG_MAX;
    norm[index] = -1.0;
  }
  if (SymmetryBasisRepresentativeDirectoryReady(X->Sym) != TRUE ||
      SymmetryResolveRepresentativeBatch(
          X->Sym->representative_directory,
          keys, (uint64_t)target_count,
          global_beta, norm) != 0) {
    fprintf(stderr, "%s: distributed transition batch failed\n", label);
    exit(1);
  }
  for (index = 0UL; index < target_count; index++) {
    if (global_beta[index] != targets[index].global_beta ||
        memcmp(&norm[index], &targets[index].norm,
               sizeof(norm[index])) != 0) {
      fprintf(stderr,
              "%s: transition batch mismatch for model %d key %lu\n",
              label, (int)model, keys[index]);
      exit(1);
    }
  }
  free(keys);
  free(global_beta);
  free(norm);
}

static int unresolved_request_contains(
    const struct SymmetryUnresolvedMatvecBlock *block,
    unsigned long int key)
{
  size_t left = 0U;
  size_t right = block->request_count;
  while (left < right) {
    size_t middle = left + (right - left) / 2U;
    if (block->request_keys[middle] < key) {
      left = middle + 1U;
    } else {
      right = middle;
    }
  }
  return left < block->request_count &&
      block->request_keys[left] == key;
}

static void assert_c2_unresolved_blocks(
    const struct BindStruct *X,
    enum C5ReferenceModel model,
    const char *label)
{
  unsigned long int local_row_begin = 0UL;
  unsigned int term_count =
      model == C5_REFERENCE_SPIN
          ? X->Def.NExchangeCoupling : X->Def.EDNTransfer;
  unsigned int term_step =
      model == C5_REFERENCE_SPIN ? 1U : 2U;
  if (X == NULL || X->Sym == NULL ||
      X->Sym->basis_layout != SYMMETRY_BASIS_DISTRIBUTED) {
    fprintf(stderr, "%s: unresolved block fixture is not distributed\n",
            label);
    exit(1);
  }
  do {
    struct SymmetryUnresolvedMatvecBlock block;
    struct SymmetryUnresolvedMatvecBlock failed;
    unsigned long int block_row;
    unsigned long int local_row_count =
        X->Sym->local_dim - local_row_begin;
    size_t saved_peak;
    if (local_row_count > 2UL) local_row_count = 2UL;
    memset(&block, 0, sizeof(block));
    memset(&failed, 0, sizeof(failed));
    if (BuildSymmetryUnresolvedMatvecBlock(
            X, local_row_begin, local_row_count,
            HPHI_SYMMETRY_PLAN_BLOCK_MEMORY_BYTES, &block) != 0) {
      fprintf(stderr, "%s: unresolved block build failed\n", label);
      exit(1);
    }
    assert_ulong_eq(block.local_row_begin, local_row_begin, label);
    assert_ulong_eq(block.local_row_count, local_row_count, label);
    assert_int_eq(block.row_ptr != NULL, 1, label);
    assert_ulong_eq(
        (unsigned long int)block.transition_count,
        (unsigned long int)block.row_ptr[local_row_count], label);
    assert_ulong_eq(
        (unsigned long int)block.offdiagonal_count,
        (unsigned long int)(block.transition_count -
                            (size_t)local_row_count),
        label);
    assert_int_eq(
        block.storage_bytes <= block.temporary_peak_bytes &&
            block.temporary_peak_bytes <=
                (size_t)block.memory_byte_limit,
        1, label);
    if (block.offdiagonal_count == 0U) {
      assert_int_eq(
          block.request_count == 0U && block.request_keys == NULL,
          1, label);
    } else {
      size_t request;
      assert_int_eq(block.request_keys != NULL, 1, label);
      assert_int_eq(
          block.request_count > 0U &&
              block.request_count <= block.offdiagonal_count,
          1, label);
      for (request = 1U; request < block.request_count; request++) {
        assert_int_eq(
            block.request_keys[request - 1U] <
                block.request_keys[request],
            1, label);
      }
    }
    for (block_row = 0UL;
         block_row < local_row_count;
         block_row++) {
      unsigned long int local_index =
          local_row_begin + block_row + 1UL;
      const struct SymmetryBasisVector *source =
          SymmetryBasisLocalEntry(X->Sym, local_index);
      size_t transition = block.row_ptr[block_row];
      size_t end = block.row_ptr[block_row + 1UL];
      unsigned int term;
      double complex expected_diagonal;
      if (source == NULL || transition >= end) {
        fprintf(stderr, "%s: unresolved row source is invalid\n", label);
        exit(1);
      }
      expected_diagonal = source->diagonal;
      assert_int_eq(
          block.transitions[transition].is_diagonal, TRUE, label);
      assert_ulong_eq(
          block.transitions[transition].rep_state, 0UL, label);
      assert_complex_bitwise(
          block.transitions[transition].phased_hval,
          expected_diagonal, label);
      transition++;
      for (term = 0U; term < term_count; term += term_step) {
        struct SymmetryRepresentativeResult representative;
        unsigned long int out_state;
        double complex hval;
        double complex expected_phased_hval;
        int applied = c5_transition(
            X, model, source->rep_state, term, &out_state, &hval);
        if (applied < 0) {
          fprintf(stderr, "%s: unresolved transition oracle failed\n",
                  label);
          exit(1);
        }
        if (applied == 0) continue;
        if (transition >= end ||
            SymmetryFindRepresentative(
                X, out_state, &representative) != 0) {
          fprintf(stderr, "%s: unresolved transition order mismatch\n",
                  label);
          exit(1);
        }
        expected_phased_hval = hval * representative.phase;
        assert_int_eq(
            block.transitions[transition].is_diagonal, FALSE, label);
        assert_ulong_eq(
            block.transitions[transition].rep_state,
            representative.rep_state, label);
        assert_complex_bitwise(
            block.transitions[transition].phased_hval,
            expected_phased_hval, label);
        assert_int_eq(
            unresolved_request_contains(
                &block, representative.rep_state),
            TRUE, label);
        transition++;
      }
      assert_ulong_eq(
          (unsigned long int)transition,
          (unsigned long int)end, label);
    }
    saved_peak = block.temporary_peak_bytes;
    FreeSymmetryUnresolvedMatvecBlock(&block);
    assert_int_eq(
        block.row_ptr == NULL && block.transitions == NULL &&
            block.request_keys == NULL &&
            block.transition_count == 0U &&
            block.offdiagonal_count == 0U &&
            block.request_count == 0U,
        1, label);
    assert_int_eq(saved_peak > 0U, 1, label);
    assert_int_eq(
        BuildSymmetryUnresolvedMatvecBlock(
            X, local_row_begin, local_row_count,
            (uint64_t)(saved_peak - 1U), &failed),
        -1, "unresolved block enforces its memory cap");
    assert_int_eq(
        failed.row_ptr == NULL && failed.transitions == NULL &&
            failed.request_keys == NULL &&
            failed.transition_count == 0U &&
            failed.offdiagonal_count == 0U &&
            failed.request_count == 0U,
        1, "failed unresolved build preserves empty output");
    local_row_begin += local_row_count;
  } while (local_row_begin < X->Sym->local_dim);
}

static void assert_c3_plan_matches_replicated(
    struct BindStruct *X,
    const size_t *reference_row_ptr,
    const unsigned long int *reference_columns,
    const double complex *reference_values,
    unsigned long int reference_local_offset,
    unsigned long int reference_local_dim,
    size_t reference_nnz,
    const char *label)
{
  struct SymmetryRepresentativeBatchStats before;
  struct SymmetryRepresentativeBatchStats after;
  struct SymmetryDistributedMatvecPlanOptions options;
  struct SymmetryRepresentativeBatchOptions directory_options;
  int pass;
  memset(&options, 0, sizeof(options));
  memset(&directory_options, 0, sizeof(directory_options));
  directory_options.corrupt_response_rank = -1;
  options.block_memory_byte_limit =
      HPHI_SYMMETRY_PLAN_BLOCK_MEMORY_BYTES;

  assert_int_eq(
      GetSymmetryRepresentativeDirectoryBatchStats(
          X->Sym->representative_directory, &before),
      0, label);
  options.local_rows_per_block = 1UL;
  options.block_memory_byte_limit = 1U;
  assert_int_eq(
      BuildSymmetryDistributedMatvecPlanWithOptions(
          X, &options),
      -1, "distributed plan enforces collective block memory cap");
  assert_int_eq(X->Sym->matvec_plan == NULL, 1, label);
  assert_int_eq(
      GetSymmetryRepresentativeDirectoryBatchStats(
          X->Sym->representative_directory, &after),
      0, label);
  assert_int_eq(
      representative_batch_stats_equal(&after, &before),
      1, "pre-batch block cap failure preserves directory stats");

  options.block_memory_byte_limit =
      HPHI_SYMMETRY_PLAN_BLOCK_MEMORY_BYTES;
  options.directory_options = &directory_options;
  directory_options.corrupt_response_rank = 0;
  assert_int_eq(
      BuildSymmetryDistributedMatvecPlanWithOptions(
          X, &options),
      -1, "injected directory failure aborts collective plan build");
  assert_int_eq(X->Sym->matvec_plan == NULL, 1, label);
  assert_int_eq(
      GetSymmetryRepresentativeDirectoryBatchStats(
          X->Sym->representative_directory, &after),
      0, label);
  assert_int_eq(
      representative_batch_stats_equal(&after, &before),
      1, "failed directory wave preserves published directory stats");
  directory_options.corrupt_response_rank = -1;
  options.directory_options = NULL;
  for (pass = 0; pass < 2; pass++) {
    const struct SymmetryMatvecPlan *plan;
    unsigned long int covered_rows = 0UL;
    size_t flattened_nnz = 0U;
    size_t block_index;
    uint64_t expected_rounds;
    uint64_t expected_local_rounds;
    int rank;
    if (pass == 0) {
      options.local_rows_per_block = 1UL;
      options.directory_options = NULL;
    } else {
      options.local_rows_per_block = 2UL;
      directory_options.chunk_limit = 1U;
      directory_options.force_chunked = 1;
      directory_options.debug_echo = 1;
      options.directory_options = &directory_options;
    }
    expected_local_rounds =
        (uint64_t)(X->Sym->local_dim /
                   options.local_rows_per_block);
    if (X->Sym->local_dim %
            options.local_rows_per_block != 0UL) {
      expected_local_rounds++;
    }
    expected_rounds = 0U;
    for (rank = 0; rank < nproc; rank++) {
      unsigned long int rank_local_dim =
          X->Sym->rank_offsets[rank + 1] -
          X->Sym->rank_offsets[rank];
      uint64_t rank_rounds =
          (uint64_t)(rank_local_dim /
                     options.local_rows_per_block);
      if (rank_local_dim %
              options.local_rows_per_block != 0UL) {
        rank_rounds++;
      }
      if (rank_rounds > expected_rounds) {
        expected_rounds = rank_rounds;
      }
    }
    assert_int_eq(
        GetSymmetryRepresentativeDirectoryBatchStats(
            X->Sym->representative_directory, &before),
        0, label);
    if (BuildSymmetryDistributedMatvecPlanWithOptions(
            X, &options) != 0) {
      fprintf(stderr, "%s: distributed plan build failed pass %d\n",
              label, pass);
      exit(1);
    }
    assert_int_eq(
        GetSymmetryRepresentativeDirectoryBatchStats(
            X->Sym->representative_directory, &after),
        0, label);
    assert_int_eq(
        after.directory_batch_calls -
                before.directory_batch_calls ==
            expected_rounds,
        1, "one directory batch is used per local block wave");
    if (expected_local_rounds < expected_rounds) {
      assert_int_eq(
          after.directory_batch_calls -
                  before.directory_batch_calls ==
              expected_rounds,
          1, "rank without a local block participates in zero-request waves");
    }
    plan = X->Sym->matvec_plan;
    assert_int_eq(plan != NULL && plan->ready == TRUE, 1, label);
    assert_int_eq(plan->columns_remapped, FALSE, label);
    assert_int_eq(plan->halo.ready, FALSE, label);
    assert_ulong_eq(plan->dim, X->Sym->dim, label);
    assert_ulong_eq(plan->local_offset, reference_local_offset, label);
    assert_ulong_eq(plan->local_dim, reference_local_dim, label);
    assert_ulong_eq((unsigned long int)plan->nnz,
                    (unsigned long int)reference_nnz, label);
    assert_ulong_eq(
        (unsigned long int)plan->build_local_wave_count,
        (unsigned long int)expected_local_rounds, label);
    assert_ulong_eq(
        (unsigned long int)plan->build_max_wave_count,
        (unsigned long int)expected_rounds, label);
    assert_int_eq(
        plan->build_temporary_peak_bytes > 0U &&
            plan->build_memory_warning_byte_threshold ==
                (uint64_t)HPHI_SYMMETRY_MEMORY_WARN_BYTES &&
            plan->build_memory_byte_limit ==
                (uint64_t)HPHI_SYMMETRY_PLAN_BLOCK_MEMORY_BYTES &&
            plan->build_temporary_peak_bytes <=
                (size_t)plan->build_memory_byte_limit,
        1, "distributed plan memory policy stats mismatch");
    for (block_index = 0U;
         block_index < SymmetryMatvecPlanBlockCount(plan);
         block_index++) {
      struct SymmetryMatvecBlockView view;
      unsigned long int block_row;
      assert_int_eq(
          SymmetryMatvecPlanGetBlockView(
              plan, block_index, &view),
          0, label);
      assert_ulong_eq(view.local_row_begin, covered_rows, label);
      for (block_row = 0UL;
           block_row < view.local_row_count;
           block_row++) {
        unsigned long int local_row =
            covered_rows + block_row;
        size_t expected_begin = reference_row_ptr[local_row];
        size_t expected_end = reference_row_ptr[local_row + 1UL];
        size_t actual_begin = view.row_ptr[block_row];
        size_t actual_end = view.row_ptr[block_row + 1UL];
        size_t row_nnz = expected_end - expected_begin;
        assert_ulong_eq(
            (unsigned long int)(actual_end - actual_begin),
            (unsigned long int)row_nnz, label);
        if (row_nnz > 0U) {
          assert_int_eq(
              memcmp(
                  &view.global_columns[actual_begin],
                  &reference_columns[expected_begin],
                  row_nnz * sizeof(*reference_columns)) == 0,
              1, "distributed/global replicated columns are bitwise exact");
          assert_int_eq(
              memcmp(
                  &view.values[actual_begin],
                  &reference_values[expected_begin],
                  row_nnz * sizeof(*reference_values)) == 0,
              1, "distributed/global replicated values are bitwise exact");
        }
        flattened_nnz += row_nnz;
      }
      covered_rows += view.local_row_count;
    }
    assert_ulong_eq(covered_rows, reference_local_dim, label);
    assert_ulong_eq(
        (unsigned long int)flattened_nnz,
        (unsigned long int)reference_nnz, label);
    FreeSymmetryMatvecPlan(X->Sym->matvec_plan);
    X->Sym->matvec_plan = NULL;
  }
}

static double complex c4_fixed_input_value(
    unsigned long int global_beta)
{
  return 0.125 * (double)global_beta +
      I * 0.0625 * (double)(global_beta + 1UL);
}

static void assert_c4_distributed_solver_plan(
    struct BindStruct *X,
    const size_t *reference_row_ptr,
    const unsigned long int *reference_columns,
    const double complex *reference_values,
    unsigned long int reference_local_offset,
    unsigned long int reference_local_dim,
    const char *label)
{
  struct SymmetryDistributedMatvecPlanOptions options;
  double complex *input;
  double complex *output;
  double complex *expected;
  double complex expected_prdct = 0.0;
  unsigned long int local_row;
  size_t vector_count =
      reference_local_dim == 0UL
          ? 1U : (size_t)reference_local_dim + 1U;
  memset(&options, 0, sizeof(options));
  options.local_rows_per_block = 1UL;
  options.block_memory_byte_limit =
      HPHI_SYMMETRY_PLAN_BLOCK_MEMORY_BYTES;
  input = (double complex *)calloc(vector_count, sizeof(*input));
  output = (double complex *)calloc(vector_count, sizeof(*output));
  expected = (double complex *)calloc(vector_count, sizeof(*expected));
  if (input == NULL || output == NULL || expected == NULL) {
    fprintf(stderr, "%s: C4 fixed-vector allocation failed\n", label);
    exit(1);
  }
  for (local_row = 0UL;
       local_row < reference_local_dim; local_row++) {
    size_t p;
    unsigned long int global_alpha =
        reference_local_offset + local_row + 1UL;
    input[local_row + 1UL] =
        c4_fixed_input_value(global_alpha);
    for (p = reference_row_ptr[local_row];
         p < reference_row_ptr[local_row + 1UL]; p++) {
      expected[local_row + 1UL] +=
          reference_values[p] *
          c4_fixed_input_value(reference_columns[p]);
    }
    expected_prdct +=
        conj(input[local_row + 1UL]) *
        expected[local_row + 1UL];
  }
  if (BuildSymmetryDistributedMatvecPlanForSolverWithOptions(
          X, &options) != 0) {
    fprintf(stderr, "%s: C4 solver plan build failed\n", label);
    exit(1);
  }
  assert_int_eq(
      X->Sym->matvec_plan != NULL &&
          X->Sym->matvec_plan->ready == TRUE &&
          X->Sym->matvec_plan->columns_remapped == TRUE &&
          X->Sym->matvec_plan->halo.ready == TRUE,
      1, label);
  assert_int_eq(
      X->Sym->matvec_mode == SYMMETRY_MATVEC_MODE_PLAN &&
          X->Sym->vector_exchange_mode ==
              SYMMETRY_VECTOR_EXCHANGE_HALO,
      1, label);
  assert_int_eq(
      X->Sym->mpi_full_v1 == NULL &&
          X->Sym->mpi_recvcounts == NULL &&
          X->Sym->mpi_displs == NULL,
      1, "distributed solver plan allocates no full vector state");
  assert_int_eq(
      X->Sym->representative_directory == NULL &&
          X->Sym->representative_directory_stats_ready == TRUE &&
          X->Sym->representative_directory_heavy_storage_released == TRUE,
      1, "distributed solver releases representative directory storage");
  assert_int_eq(
      X->Sym->representative_directory_info.splitter_bytes > 0U &&
          (X->Sym->local_dim == 0UL
               ? X->Sym->representative_directory_index_stats
                         .table_bytes == 0U
               : X->Sym->representative_directory_index_stats
                         .table_bytes > 0U) &&
          X->Sym->representative_directory_batch_stats
                  .directory_batch_calls > 0U,
      1, "distributed solver retains lightweight directory statistics");
  assert_int_eq(
      X->Sym->matvec_plan->halo.ghost_global_index == NULL,
      1, "distributed solver releases remap-only ghost global indices");
  assert_ulong_eq(
      (unsigned long int)SymmetryMatvecPlanBlockCount(
          X->Sym->matvec_plan),
      reference_local_dim == 0UL ? 1UL : reference_local_dim,
      label);
  X->Large.prdct = 0.0;
  assert_int_eq(mltplySpinSym(X, output, input), 0, label);
  for (local_row = 0UL;
       local_row < reference_local_dim; local_row++) {
    assert_complex_bitwise(
        output[local_row + 1UL],
        expected[local_row + 1UL],
        "distributed halo Hv is bitwise exact");
  }
  assert_complex_close(
      X->Large.prdct, expected_prdct, 1.0e-13,
      "distributed halo local prdct");
  memset(output, 0, vector_count * sizeof(*output));
  X->Large.prdct = 0.0;
  assert_int_eq(mltplySpinSym(X, output, input), 0,
                "released directory supports repeated distributed matvec");
  for (local_row = 0UL;
       local_row < reference_local_dim; local_row++) {
    assert_complex_bitwise(
        output[local_row + 1UL],
        expected[local_row + 1UL],
        "repeated distributed halo Hv is bitwise exact");
  }
  assert_complex_close(
      X->Large.prdct, expected_prdct, 1.0e-13,
      "repeated distributed halo local prdct");
  assert_int_eq(
      X->Sym->matvec_plan->matvec_calls == 2ULL &&
          X->Sym->matvec_plan->halo.exchange_calls == 2ULL &&
          X->Sym->matvec_plan->input_allgather_calls == 0ULL,
      1, label);
  FreeSymmetryMatvecPlan(X->Sym->matvec_plan);
  X->Sym->matvec_plan = NULL;
  free(input);
  free(output);
  free(expected);
}

static void assert_c6_full_basis_directory_batch(
    struct SymmetryBasisRuntime *sym,
    const struct SymmetryBasisVector *reference_basis,
    const char *label)
{
  struct SymmetryRepresentativeBatchOptions options;
  unsigned long int *keys;
  unsigned long int *global_beta;
  double *norm;
  unsigned long int missing_key = 0UL;
  unsigned long int index;
  int pass;
  if (sym == NULL || reference_basis == NULL || sym->dim == 0UL ||
      sym->dim > (unsigned long int)(SIZE_MAX / sizeof(*keys)) ||
      sym->dim > (unsigned long int)(SIZE_MAX / sizeof(*norm))) {
    fprintf(stderr, "%s: invalid full-basis directory fixture\n", label);
    exit(1);
  }
  keys = (unsigned long int *)malloc((size_t)sym->dim * sizeof(*keys));
  global_beta =
      (unsigned long int *)malloc((size_t)sym->dim *
                                  sizeof(*global_beta));
  norm = (double *)malloc((size_t)sym->dim * sizeof(*norm));
  if (keys == NULL || global_beta == NULL || norm == NULL) {
    fprintf(stderr, "%s: full-basis directory allocation failed\n", label);
    exit(1);
  }
  for (index = 1UL; index <= sym->dim; index++) {
    keys[index - 1UL] = reference_basis[index].rep_state;
  }
  memset(&options, 0, sizeof(options));
  options.corrupt_response_rank = -1;
  for (pass = 0; pass < 2; pass++) {
    int status;
    for (index = 0UL; index < sym->dim; index++) {
      global_beta[index] = ULONG_MAX;
      norm[index] = -1.0;
    }
    if (pass == 0) {
      status = SymmetryResolveRepresentativeBatch(
          sym->representative_directory, keys, (uint64_t)sym->dim,
          global_beta, norm);
    } else {
      options.force_chunked = 1;
      options.chunk_limit = 1U;
      options.debug_echo = 1;
      status = SymmetryResolveRepresentativeBatchWithOptions(
          sym->representative_directory, keys, (uint64_t)sym->dim,
          global_beta, norm, &options);
    }
    if (status != 0) {
      fprintf(stderr, "%s: full-basis directory batch failed\n", label);
      exit(1);
    }
    for (index = 1UL; index <= sym->dim; index++) {
      if (global_beta[index - 1UL] != index ||
          memcmp(&norm[index - 1UL], &reference_basis[index].norm,
                 sizeof(norm[index - 1UL])) != 0) {
        fprintf(stderr,
                "%s: full-basis directory mismatch at beta %lu pass %d\n",
                label, index, pass);
        exit(1);
      }
    }
  }
  for (index = 1UL; index <= sym->dim; index++) {
    if (missing_key < reference_basis[index].rep_state) break;
    if (missing_key == reference_basis[index].rep_state) {
      if (missing_key == ULONG_MAX) {
        fprintf(stderr, "%s: missing-key search overflow\n", label);
        exit(1);
      }
      missing_key++;
    }
  }
  {
    unsigned long int missing_beta = ULONG_MAX;
    double missing_norm = -1.0;
    const double positive_zero = 0.0;
    if (SymmetryResolveRepresentativeBatch(
            sym->representative_directory, &missing_key, 1U,
            &missing_beta, &missing_norm) != 0 ||
        missing_beta != 0UL ||
        memcmp(&missing_norm, &positive_zero,
               sizeof(missing_norm)) != 0) {
      fprintf(stderr, "%s: nonexistent representative mismatch\n", label);
      exit(1);
    }
  }
  free(keys);
  free(global_beta);
  free(norm);
}

static struct SymmetryRepresentativeDirectory *
build_c6_zero_dimension_directory(const char *label)
{
  struct SymmetryRepresentativeDirectory *directory = NULL;
  unsigned long int *rank_offsets =
      (unsigned long int *)calloc((size_t)nproc + 1U,
                                  sizeof(*rank_offsets));
  if (SumMPI_i(rank_offsets == NULL ? 1 : 0) != 0) {
    fprintf(stderr, "%s: zero-dimension offset allocation failed\n", label);
    exit(1);
  }
  if (BuildSymmetryRepresentativeDirectory(
          NULL, 0UL, 0UL, 0UL, 0UL, rank_offsets,
          myrank, nproc, &directory) != 0 ||
      SymmetryRepresentativeDirectoryReady(directory) == 0) {
    fprintf(stderr, "%s: zero-dimension directory build failed\n", label);
    exit(1);
  }
  free(rank_offsets);
  return directory;
}

static void assert_c5_distributed_layout_model(
    enum C5ReferenceModel model,
    const char *label)
{
  struct BindStruct X;
  struct SymmetryBasisDigest replicated_digest;
  struct SymmetryBasisDigest distributed_digest;
  struct SymmetryBasisDigest reference_digest;
  struct SymmetryBasisRuntime reference_sym;
  struct SymmetryCanonicalResult canonical;
  struct SymmetryRepresentativeResult representative;
  struct SymmetryBasisVector *reference_basis;
  struct C5DirectoryTarget *directory_targets;
  size_t *reference_row_ptr;
  unsigned long int *reference_columns;
  double complex *reference_values;
  unsigned long int reference_local_offset;
  unsigned long int reference_local_dim;
  size_t reference_nnz;
  unsigned long int raw_dim;
  unsigned long int dim;
  unsigned long int directory_target_count;
  unsigned long int local_index;
  unsigned long long global_digest_count;
  unsigned long long global_digest_xor;
  unsigned long long global_digest_sum;
  double diagonal;
  int transition_found_count;
  int transition_interior_miss_count;
  int transition_range_miss_count;

  setup_c5_reference_bind(&X, model);
  raw_dim = X.Check.idim_max;
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: replicated reference build failed\n", label);
    exit(1);
  }
  assert_int_eq(X.Sym->basis_layout, SYMMETRY_BASIS_REPLICATED, label);
  assert_int_eq(X.Sym->representative_directory == NULL, 1, label);
  assert_int_eq(
      SymmetryBasisRepresentativeDirectoryReady(X.Sym), FALSE, label);
  assert_int_eq(
      ComputeSymmetryBasisDigest(X.Sym, &replicated_digest), 0, label);
  assert_int_eq(
      replicated_digest.algorithm ==
          SYMMETRY_BASIS_DIGEST_REPLICATED_FNV1A64,
      1, label);
  assert_int_eq(replicated_digest.fnv1a64 != 0U, 1, label);
  dim = X.Sym->dim;
  directory_targets = collect_c5_directory_targets(
      &X, model, &directory_target_count,
      &transition_found_count, &transition_interior_miss_count,
      &transition_range_miss_count, label);
  reference_basis = (struct SymmetryBasisVector *)malloc(
      ((size_t)dim + 1U) * sizeof(*reference_basis));
  if (reference_basis == NULL) {
    fprintf(stderr, "%s: reference basis allocation failed\n", label);
    exit(1);
  }
  memcpy(reference_basis, X.Sym->basis,
         ((size_t)dim + 1U) * sizeof(*reference_basis));
  if (setenv(
          "HPHI_SYMMETRY_VECTOR_EXCHANGE", "allgather", 1) != 0 ||
      ActivateSymmetryBasisDimension(&X) != 0 ||
      BuildSymmetryMatvecPlan(&X) != 0) {
    fprintf(stderr, "%s: replicated plan oracle build failed\n", label);
    exit(1);
  }
  unsetenv("HPHI_SYMMETRY_VECTOR_EXCHANGE");
  reference_local_offset = X.Sym->matvec_plan->local_offset;
  reference_local_dim = X.Sym->matvec_plan->local_dim;
  reference_nnz = X.Sym->matvec_plan->nnz;
  reference_row_ptr = (size_t *)malloc(
      ((size_t)reference_local_dim + 1U) *
      sizeof(*reference_row_ptr));
  reference_columns =
      reference_nnz == 0U
          ? NULL
          : (unsigned long int *)malloc(
                reference_nnz * sizeof(*reference_columns));
  reference_values =
      reference_nnz == 0U
          ? NULL
          : (double complex *)malloc(
                reference_nnz * sizeof(*reference_values));
  if (reference_row_ptr == NULL ||
      (reference_nnz > 0U &&
       (reference_columns == NULL || reference_values == NULL))) {
    fprintf(stderr, "%s: replicated plan oracle allocation failed\n",
            label);
    exit(1);
  }
  memcpy(
      reference_row_ptr, X.Sym->matvec_plan->row_ptr,
      ((size_t)reference_local_dim + 1U) *
          sizeof(*reference_row_ptr));
  if (reference_nnz > 0U) {
    memcpy(
        reference_columns, X.Sym->matvec_plan->col_index,
        reference_nnz * sizeof(*reference_columns));
    memcpy(
        reference_values, X.Sym->matvec_plan->values,
        reference_nnz * sizeof(*reference_values));
  }
  FreeSymmetryBasis(X.Sym);
  X.Sym = NULL;
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;

  setup_c5_reference_bind(&X, model);
  assert_ulong_eq(X.Check.idim_max, raw_dim, label);
  if (BuildSymmetryBasisForLayout(
          &X, SYMMETRY_BASIS_DISTRIBUTED) != 0) {
    fprintf(stderr, "%s: distributed layout build failed\n", label);
    exit(1);
  }
  assert_int_eq(X.Sym->basis_layout, SYMMETRY_BASIS_DISTRIBUTED, label);
  assert_ulong_eq(X.Sym->dim, dim, label);
  assert_int_eq(X.Sym->basis == NULL, 1, label);
  assert_ulong_eq(X.Sym->capacity, 0UL, label);
  assert_int_eq(X.Sym->rep_hash_size == 0UL &&
                    X.Sym->rep_hash_keys == NULL &&
                    X.Sym->rep_hash_values == NULL,
                1, label);
  assert_int_eq(X.Sym->rank_offsets != NULL, 1, label);
  assert_int_eq(X.Sym->representative_directory != NULL, 1, label);
  assert_int_eq(
      SymmetryRepresentativeDirectoryReady(
          X.Sym->representative_directory),
      1, label);
  assert_int_eq(
      SymmetryBasisRepresentativeDirectoryReady(X.Sym), TRUE, label);
  assert_ulong_eq(X.Sym->rank_offsets[0], 0UL, label);
  assert_ulong_eq(X.Sym->rank_offsets[nproc], dim, label);
  assert_ulong_eq(X.Sym->rank_offsets[myrank],
                  X.Sym->local_offset, label);
  assert_ulong_eq(X.Sym->rank_offsets[myrank + 1] -
                      X.Sym->rank_offsets[myrank],
                  X.Sym->local_dim, label);
  assert_int_eq(
      SymmetryBasisOwnedStorageReady(X.Sym, X.Sym->local_dim),
      TRUE, label);
  assert_int_eq(
      X.Sym->local_dim == 0UL ||
          (X.Sym->local_basis != NULL &&
           X.Sym->local_capacity >= X.Sym->local_dim + 1UL),
      1, label);
  assert_ulong_eq(
      (unsigned long int)X.Sym->distribution_stats.global_entries,
      dim, label);
  assert_ulong_eq(
      (unsigned long int)X.Sym->distribution_stats.rebalance_recv_entries,
      X.Sym->local_dim, label);
  assert_int_eq(
      (unsigned long long)X.Sym->local_capacity *
              (unsigned long long)sizeof(*X.Sym->local_basis) >=
          (unsigned long long)(X.Sym->local_dim + 1UL) *
              (unsigned long long)sizeof(*X.Sym->local_basis),
      1, label);

  for (local_index = 1UL;
       local_index <= X.Sym->local_dim;
       local_index++) {
    const struct SymmetryBasisVector *entry =
        SymmetryBasisLocalEntry(X.Sym, local_index);
    unsigned long int global_beta = X.Sym->local_offset + local_index;
    unsigned long int lookup_local_index = ULONG_MAX;
    unsigned long int lookup_global_beta = ULONG_MAX;
    double lookup_norm = -1.0;
    assert_int_eq(entry != NULL, 1, label);
    assert_int_eq(c1_basis_vector_fields_equal(
                      entry, &reference_basis[global_beta]),
                  1, label);
    assert_int_eq(
        GetOwnedHamiltonianDiagonal(&X, local_index, &diagonal), 0, label);
    assert_complex_close(
        diagonal, reference_basis[global_beta].diagonal, 0.0, label);
    assert_int_eq(
        SymmetryLookupDirectoryLocalRepresentative(
            X.Sym->representative_directory, entry->rep_state,
            &lookup_local_index, &lookup_global_beta,
            &lookup_norm, NULL),
        0, label);
    assert_ulong_eq(lookup_local_index, local_index, label);
    assert_ulong_eq(lookup_global_beta, global_beta, label);
    assert_int_eq(
        memcmp(&lookup_norm, &reference_basis[global_beta].norm,
               sizeof(lookup_norm)) == 0,
        1, label);
  }
  if (X.Sym->local_dim < dim) {
    unsigned long int remote_beta =
        X.Sym->local_offset > 0UL
            ? 1UL
            : X.Sym->local_offset + X.Sym->local_dim + 1UL;
    unsigned long int lookup_local_index = ULONG_MAX;
    unsigned long int lookup_global_beta = ULONG_MAX;
    double lookup_norm = -1.0;
    const double positive_zero = 0.0;
    assert_int_eq(remote_beta >= 1UL && remote_beta <= dim, 1, label);
    assert_int_eq(
        SymmetryLookupDirectoryLocalRepresentative(
            X.Sym->representative_directory,
            reference_basis[remote_beta].rep_state,
            &lookup_local_index, &lookup_global_beta,
            &lookup_norm, NULL),
        0, label);
    assert_ulong_eq(lookup_local_index, 0UL, label);
    assert_ulong_eq(lookup_global_beta, 0UL, label);
    assert_int_eq(
        memcmp(&lookup_norm, &positive_zero, sizeof(lookup_norm)) == 0,
        1, label);
  }
  assert_int_eq(SymmetryBasisLocalEntry(X.Sym, 0UL) == NULL, 1, label);
  assert_int_eq(
      SymmetryBasisLocalEntry(
          X.Sym, X.Sym->local_dim + 1UL) == NULL,
      1, label);
  assert_int_eq(
      SymmetryBasisReplicatedGlobalEntry(X.Sym, 1UL) == NULL,
      1, label);
  assert_int_eq(
      SymmetryFindRepresentative(&X, list_1[1], &representative), 0,
      label);
  assert_int_eq(
      representative.op_rep_to_state < X.Def.NSymTrans, 1, label);
  assert_c6_full_basis_directory_batch(X.Sym, reference_basis, label);
  assert_c5_directory_transition_batch(
      &X, model, directory_targets, directory_target_count, label);
  {
    struct SymmetryRepresentativeBatchStats before;
    struct SymmetryRepresentativeBatchStats after;
    assert_int_eq(
        GetSymmetryRepresentativeDirectoryBatchStats(
            X.Sym->representative_directory, &before),
        0, label);
    assert_c2_unresolved_blocks(&X, model, label);
    assert_int_eq(
        GetSymmetryRepresentativeDirectoryBatchStats(
            X.Sym->representative_directory, &after),
        0, label);
    assert_int_eq(
        representative_batch_stats_equal(&after, &before),
        1, "unresolved block build performs no directory batch");
  }
  assert_c3_plan_matches_replicated(
      &X, reference_row_ptr, reference_columns,
      reference_values, reference_local_offset,
      reference_local_dim, reference_nnz, label);

  assert_int_eq(
      ComputeSymmetryBasisDigest(X.Sym, &distributed_digest), 0, label);
  memset(&reference_sym, 0, sizeof(reference_sym));
  reference_sym.enabled = TRUE;
  reference_sym.basis_layout = SYMMETRY_BASIS_DISTRIBUTED;
  reference_sym.dim = dim;
  reference_sym.local_offset = X.Sym->local_offset;
  reference_sym.local_dim = X.Sym->local_dim;
  reference_sym.local_capacity = X.Sym->local_dim + 1UL;
  reference_sym.local_basis = reference_basis + X.Sym->local_offset;
  reference_sym.rank_offsets = X.Sym->rank_offsets;
  assert_int_eq(
      ComputeSymmetryBasisDigest(&reference_sym, &reference_digest),
      0, label);
  assert_int_eq(
      distributed_digest.algorithm ==
          SYMMETRY_BASIS_DIGEST_DISTRIBUTED_GLOBAL_BETA,
      1, label);
  assert_ulong_eq(
      (unsigned long int)distributed_digest.count,
      X.Sym->local_dim, label);
  assert_int_eq(
      distributed_digest.xor_hash == reference_digest.xor_hash &&
          distributed_digest.sum_hash == reference_digest.sum_hash,
      1, label);
  global_digest_count = distributed_digest.count;
  global_digest_xor = distributed_digest.xor_hash;
  global_digest_sum = distributed_digest.sum_hash;
#ifdef MPI
  if (nproc > 1) {
    unsigned long long reduced_count;
    unsigned long long reduced_xor;
    unsigned long long reduced_sum;
    if (MPI_Allreduce(&global_digest_count, &reduced_count, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&global_digest_xor, &reduced_xor, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_BXOR,
                      MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&global_digest_sum, &reduced_sum, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      fprintf(stderr, "%s: distributed digest reduction failed\n", label);
      exit(1);
    }
    global_digest_count = reduced_count;
    global_digest_xor = reduced_xor;
    global_digest_sum = reduced_sum;
  }
#endif
  assert_ulong_eq((unsigned long int)global_digest_count, dim, label);
  assert_int_eq(dim == 0UL ||
                    global_digest_xor != 0ULL ||
                    global_digest_sum != 0ULL,
                1, label);
  if (model == C5_REFERENCE_SPIN && sizeof(unsigned long int) == 8U) {
    assert_int_eq(
        global_digest_xor == UINT64_C(0x2f0fbd0254a21ba5),
        1, "C5 distributed Spin XOR literal");
    assert_int_eq(
        global_digest_sum == UINT64_C(0x10f042f5739b5c59),
        1, "C5 distributed Spin SUM literal");
  }

  {
    struct SymmetryBasisVector *saved_local_basis = X.Sym->local_basis;
    unsigned long int saved_local_capacity = X.Sym->local_capacity;
    unsigned long int saved_rank_end = X.Sym->rank_offsets[nproc];
    struct SymmetryRepresentativeDirectory *saved_directory =
        X.Sym->representative_directory;
    struct SymmetryRepresentativeDirectory *mismatched_directory;
    struct SymmetryBasisRuntime *partial_runtime;
    X.Sym->local_basis = NULL;
    assert_int_eq(
        SymmetryBasisOwnedStorageReady(X.Sym, X.Sym->local_dim),
        FALSE, label);
    assert_int_eq(
        ComputeSymmetryBasisDigest(X.Sym, &distributed_digest), -1, label);
    X.Sym->local_basis = saved_local_basis;
    X.Sym->local_capacity = X.Sym->local_dim;
    assert_int_eq(
        SymmetryBasisOwnedStorageReady(X.Sym, X.Sym->local_dim),
        FALSE, label);
    assert_int_eq(
        ComputeSymmetryBasisDigest(X.Sym, &distributed_digest), -1, label);
    assert_int_eq(ActivateSymmetryBasisDimension(&X), -1, label);
    assert_ulong_eq(X.Check.idim_max, raw_dim, label);
    X.Sym->local_capacity = saved_local_capacity;
    X.Sym->rank_offsets[nproc] =
        saved_rank_end == 0UL ? 1UL : saved_rank_end - 1UL;
    assert_int_eq(
        SymmetryBasisOwnedStorageReady(X.Sym, X.Sym->local_dim),
        FALSE, label);
    assert_int_eq(
        ComputeSymmetryBasisDigest(X.Sym, &distributed_digest), -1, label);
    X.Sym->rank_offsets[nproc] = saved_rank_end;
    assert_int_eq(
        SymmetryBasisOwnedStorageReady(X.Sym, X.Sym->local_dim),
        TRUE, label);
    X.Sym->representative_directory = NULL;
    assert_int_eq(
        SymmetryBasisRepresentativeDirectoryReady(X.Sym), FALSE, label);
    assert_int_eq(ActivateSymmetryBasisDimension(&X), -1, label);
    assert_ulong_eq(X.Check.idim_max, raw_dim, label);
    X.Sym->representative_directory = saved_directory;
    mismatched_directory = build_c6_zero_dimension_directory(label);
    X.Sym->representative_directory = mismatched_directory;
    assert_int_eq(
        SymmetryBasisRepresentativeDirectoryReady(X.Sym), FALSE, label);
    assert_int_eq(ActivateSymmetryBasisDimension(&X), -1, label);
    assert_ulong_eq(X.Check.idim_max, raw_dim, label);
    X.Sym->representative_directory = saved_directory;
    partial_runtime =
        (struct SymmetryBasisRuntime *)calloc(1U,
                                              sizeof(*partial_runtime));
    if (partial_runtime == NULL) {
      fprintf(stderr, "%s: partial runtime allocation failed\n", label);
      exit(1);
    }
    partial_runtime->representative_directory = mismatched_directory;
    FreeSymmetryBasis(partial_runtime);
    assert_int_eq(
        SymmetryBasisRepresentativeDirectoryReady(X.Sym), TRUE, label);
  }

  {
    struct SymmetryRepresentativeDirectory *saved_directory =
        X.Sym->representative_directory;
    assert_int_eq(ActivateSymmetryBasisDimension(&X), 0, label);
    assert_int_eq(
        X.Sym->representative_directory == saved_directory, 1, label);
  }
  assert_ulong_eq(X.Check.idim_max, X.Sym->local_dim, label);
  assert_ulong_eq(X.Check.idim_maxMPI, dim, label);
  assert_int_eq(
      SymmetryBasisOwnedStorageReady(X.Sym, X.Check.idim_max),
      TRUE, label);
  assert_int_eq(
      SymmetryBasisRepresentativeDirectoryReady(X.Sym), TRUE, label);
  assert_c4_distributed_solver_plan(
      &X, reference_row_ptr, reference_columns,
      reference_values, reference_local_offset,
      reference_local_dim, label);
  assert_int_eq(
      SymmetryCanonicalizeState(&X, 0UL, &canonical), -1, label);
  assert_int_eq(
      SymmetryEnumerateColumn(&X, 1UL, c5_noop_entry, NULL), -1, label);
  assert_int_eq(BuildSymmetryMatvecPlan(&X), -1, label);
  assert_int_eq(mltplySpinSym(&X, NULL, NULL), -1, label);
  assert_int_eq(X.Sym->matvec_plan == NULL &&
                    X.Sym->mpi_full_v1 == NULL &&
                    X.Sym->mpi_recvcounts == NULL &&
                    X.Sym->mpi_displs == NULL,
                1, label);

  if (model == C5_REFERENCE_SPIN && sizeof(unsigned long int) == 8U) {
    assert_ulong_eq(
        (unsigned long int)replicated_digest.fnv1a64,
        (unsigned long int)UINT64_C(0x14692a2818afe9ba),
        "C5 replicated FNV literal");
  }
  FreeSymmetryBasis(X.Sym);
  X.Sym = NULL;
  free(reference_basis);
  free(reference_row_ptr);
  free(reference_columns);
  free(reference_values);
  free(directory_targets);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

static void assert_c5_spin_interior_transition_batch(const char *label)
{
  struct BindStruct X;
  struct C5DirectoryTarget *targets;
  unsigned long int target_count;
  int found_count;
  int interior_miss_count;
  int range_miss_count;

  setup_bind(&X, 8U, 4U, 1U);
  if (BuildSymmetryBasis(&X) != 0) {
    fprintf(stderr, "%s: replicated Spin interior fixture failed\n", label);
    exit(1);
  }
  targets = collect_c5_directory_targets(
      &X, C5_REFERENCE_SPIN, &target_count,
      &found_count, &interior_miss_count, &range_miss_count, label);
  if (found_count == 0 || interior_miss_count == 0 ||
      range_miss_count == 0) {
    fprintf(stderr,
            "%s: Spin transition categories found=%d interior=%d range=%d\n",
            label, found_count, interior_miss_count, range_miss_count);
    exit(1);
  }
  FreeSymmetryBasis(X.Sym);
  X.Sym = NULL;
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;

  setup_bind(&X, 8U, 4U, 1U);
  if (BuildSymmetryBasisForLayout(
          &X, SYMMETRY_BASIS_DISTRIBUTED) != 0) {
    fprintf(stderr, "%s: distributed Spin interior fixture failed\n", label);
    exit(1);
  }
  assert_c5_directory_transition_batch(
      &X, C5_REFERENCE_SPIN, targets, target_count, label);
  FreeSymmetryBasis(X.Sym);
  X.Sym = NULL;
  free(targets);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
}

static void assert_c5_distributed_layout(const char *label)
{
  assert_c5_distributed_layout_model(C5_REFERENCE_SPIN, label);
  assert_c5_distributed_layout_model(C5_REFERENCE_SPINLESS, label);
  assert_c5_distributed_layout_model(C5_REFERENCE_HUBBARD, label);
  assert_c5_spin_interior_transition_batch(label);
}

static void assert_rank_local_basis_run_contract(const char *label)
{
  struct BindStruct X;
  struct SymmetryBasisRuntime local_sym;
  struct SymmetryBasisRun run = {NULL, 0UL, 0UL};
  struct SymmetryBasisOwnership ownership;
  struct SymmetryBasisDistributionStats distribution_stats;
  struct SymmetryBasisRun invalid_run;
  struct SymmetryBasisVector *invalid_entries;
  unsigned long int *presence;
  unsigned long int global_count;
  unsigned long long global_raw_states;
  unsigned long int index;
  int empty_rank_count;

  setup_spinless_bind(&X, 4U, 2U, 1U);
  memset(&local_sym, 0, sizeof(local_sym));
  memset(&ownership, 0, sizeof(ownership));
  memset(&distribution_stats, 0, sizeof(distribution_stats));
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

  memset(presence, 0,
         ((size_t)X.Sym->dim + 1U) * sizeof(*presence));
  if (SymmetrySampleSortBasisRun(
          &run, myrank, nproc, &distribution_stats) != 0) {
    fprintf(stderr, "%s: distributed sample sort failed\n", label);
    exit(1);
  }
  assert_ulong_eq(
      (unsigned long int)distribution_stats.global_entries,
      X.Sym->dim, label);
  assert_ulong_eq(
      (unsigned long int)distribution_stats.range_entries,
      run.count, label);
  assert_int_eq(
      distribution_stats.range_entries <=
          distribution_stats.bucket_entry_upper_bound,
      1, label);
  for (index = 1UL; index <= run.count; index++) {
    unsigned long int beta;
    int found = FALSE;
    if (index > 1UL) {
      assert_int_eq(run.entries[index - 1UL].rep_state <
                        run.entries[index].rep_state,
                    1, label);
    }
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
    if (MPI_Allreduce(MPI_IN_PLACE, presence, (int)(X.Sym->dim + 1UL),
                      MPI_UNSIGNED_LONG, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      fprintf(stderr, "%s: sample-sort presence reduction failed\n", label);
      exit(1);
    }
  }
#endif
  for (index = 1UL; index <= X.Sym->dim; index++) {
    assert_ulong_eq(presence[index], 1UL, label);
  }

  memset(presence, 0,
         ((size_t)X.Sym->dim + 1U) * sizeof(*presence));
  if (SymmetryExactRebalanceBasisRun(
          &run, myrank, nproc, &ownership,
          &distribution_stats) != 0) {
    fprintf(stderr, "%s: exact block rebalance failed\n", label);
    exit(1);
  }
  {
    unsigned long int expected_offset;
    unsigned long int expected_count;
    if (SymmetryBlockRange(X.Sym->dim, myrank, nproc,
                           &expected_offset, &expected_count) != 0) {
      fprintf(stderr, "%s: local exact block reference failed\n", label);
      exit(1);
    }
    assert_ulong_eq(ownership.dim, X.Sym->dim, label);
    assert_ulong_eq(ownership.local_offset, expected_offset, label);
    assert_ulong_eq(ownership.local_dim, expected_count, label);
    assert_ulong_eq(run.count, expected_count, label);
  }
  assert_ulong_eq(run.capacity, run.count + 1UL, label);
  for (index = 0UL; index <= (unsigned long int)nproc; index++) {
    unsigned long int expected_offset;
    unsigned long int expected_count;
    if (index == (unsigned long int)nproc) {
      assert_ulong_eq(ownership.rank_offsets[index],
                      X.Sym->dim, label);
      continue;
    }
    if (SymmetryBlockRange(X.Sym->dim, (int)index, nproc,
                           &expected_offset, &expected_count) != 0) {
      fprintf(stderr, "%s: exact ownership reference failed\n", label);
      exit(1);
    }
    assert_ulong_eq(ownership.rank_offsets[index],
                    expected_offset, label);
    assert_ulong_eq(ownership.rank_offsets[index + 1UL],
                    expected_offset + expected_count, label);
  }
  for (index = 1UL; index <= run.count; index++) {
    unsigned long int beta = ownership.local_offset + index;
    assert_int_eq(c1_basis_vector_fields_equal(
                      &run.entries[index], &X.Sym->basis[beta]),
                  1, label);
    presence[beta]++;
  }
#ifdef MPI
  if (nproc > 1) {
    if (MPI_Allreduce(MPI_IN_PLACE, presence, (int)(X.Sym->dim + 1UL),
                      MPI_UNSIGNED_LONG, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      fprintf(stderr, "%s: exact-rebalance presence reduction failed\n",
              label);
      exit(1);
    }
  }
#endif
  for (index = 1UL; index <= X.Sym->dim; index++) {
    assert_ulong_eq(presence[index], 1UL, label);
  }

  free(presence);
  FreeSymmetryBasisOwnership(&ownership);
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
  struct SymmetryBasisOwnership ownership;
  struct SymmetryBasisDistributionStats distribution_stats;
  unsigned long int global_count;

  setup_bind(&X, 4U, 0U, 1U);
  memset(&local_sym, 0, sizeof(local_sym));
  memset(&ownership, 0, sizeof(ownership));
  memset(&distribution_stats, 0, sizeof(distribution_stats));
  if (BuildRankLocalSymmetryBasisRun(&X, &local_sym, &run) != 0) {
    fprintf(stderr, "%s: empty rank-local run build failed\n", label);
    exit(1);
  }
  assert_ulong_eq(run.count, 0UL, label);
  assert_ulong_eq(run.capacity, 1UL, label);
  assert_int_eq(run.entries != NULL, 1, label);
  if (SymmetrySampleSortBasisRun(
          &run, myrank, nproc, &distribution_stats) != 0) {
    fprintf(stderr, "%s: empty sample sort failed\n", label);
    exit(1);
  }
  assert_ulong_eq(run.count, 0UL, label);
  assert_ulong_eq(run.capacity, 1UL, label);
  assert_int_eq(run.entries != NULL, 1, label);
  assert_ulong_eq(
      (unsigned long int)distribution_stats.global_entries,
      0UL, label);
  if (SymmetryExactRebalanceBasisRun(
          &run, myrank, nproc, &ownership,
          &distribution_stats) != 0) {
    fprintf(stderr, "%s: empty exact rebalance failed\n", label);
    exit(1);
  }
  assert_ulong_eq(run.count, 0UL, label);
  assert_ulong_eq(run.capacity, 1UL, label);
  assert_int_eq(run.entries != NULL, 1, label);
  assert_ulong_eq(ownership.dim, 0UL, label);
  assert_ulong_eq(ownership.local_offset, 0UL, label);
  assert_ulong_eq(ownership.local_dim, 0UL, label);
  assert_int_eq(ownership.rank_offsets != NULL, 1, label);
  {
    int peer;
    for (peer = 0; peer <= nproc; peer++) {
      assert_ulong_eq(ownership.rank_offsets[peer], 0UL, label);
    }
  }
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
  assert_int_eq(
      BuildSymmetryBasisForLayout(
          &X, SYMMETRY_BASIS_DISTRIBUTED),
      -1, "distributed builder rejects zero-dimensional sector");
  assert_int_eq(X.Sym == NULL, 1, label);
  FreeSymmetryBasisOwnership(&ownership);
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
  } else if (argc == 2 &&
             strcmp(argv[1], "--mpi-distributed-basis") == 0) {
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS ||
        MPI_Comm_size(MPI_COMM_WORLD, &nproc) != MPI_SUCCESS ||
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank) != MPI_SUCCESS) {
      fprintf(stderr, "MPI distributed-basis test initialization failed\n");
      return 1;
    }
    stdoutMPI = stderr;
    assert_c5_distributed_layout(
        "staged distributed basis matches replicated reference");
    if (myrank == 0) {
      fprintf(stdout,
              "staged distributed symmetry basis gate: PASS "
              "(%d MPI ranks)\n",
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
  assert_representative_discovery_model(
      C5_REFERENCE_SPIN,
      "Spin representative discovery is basis-layout independent");
  assert_representative_discovery_model(
      C5_REFERENCE_SPINLESS,
      "SpinlessFermion representative discovery preserves phase bits");
  assert_representative_discovery_model(
      C5_REFERENCE_HUBBARD,
      "Hubbard representative discovery preserves phase bits");
  assert_c5_distributed_layout(
      "serial staged distributed basis matches replicated reference");
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
  assert_owned_multi_block_plan(
      "owned multi-block plan preserves row order, views, remap, and apply");
  assert_zero_row_plan("local-row plan supports zero-row rank");
  assert_representative_hash_matches_basis(6, 3, 1,
                                           "C6 k=pi/3 representative hash matches basis");
  assert_hash_probe_lookup_handles_collision("representative hash probing handles collisions");
  return 0;
}
