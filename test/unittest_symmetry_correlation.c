#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "DefCommon.h"
#include "mltplySpinSym.h"
#include "makeHamSym.h"
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
double complex **Ham = NULL, *Ham_local = NULL;
int iHamPanelActive = 0, iHamSinkMode = 0;
long HamColBegin, HamColEnd, HamPanelLd;
void (*hamCollectSink)(long, long, double complex) = NULL;
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

static unsigned long int setup_full_spin_basis(unsigned int nsite)
{
  unsigned long int state, dim = 1UL << nsite;
  free(list_1);
  free(list_Diagonal);
  list_1 = (unsigned long int *)calloc(dim + 1UL, sizeof(*list_1));
  list_Diagonal = (double *)calloc(dim + 1UL, sizeof(*list_Diagonal));
  if (list_1 == NULL || list_Diagonal == NULL) {
    fprintf(stderr, "failed to allocate SpinGC basis\n");
    exit(1);
  }
  for (state = 0; state < dim; ++state) list_1[state + 1UL] = state;
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

static void setup_spingc_bind(struct BindStruct *X,
                              unsigned int nsite,
                              unsigned int momentum_index)
{
  memset(X, 0, sizeof(*X));
  setup_cyclic_def(&X->Def, nsite, momentum_index);
  X->Def.iCalcModel = SpinGC;
  X->Def.iFlgSzConserved = FALSE;
  setup_exchange_ring(&X->Def, nsite);
  X->Check.idim_max = setup_full_spin_basis(nsite);
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
#include "symmetry_terms.h"
#include "symmetry_correlation.h"

static unsigned long lcg_state = 12345UL;

static unsigned long lcg_next(void)
{
  lcg_state = lcg_state * 6364136223846793005UL + 1442695040888963407UL;
  return lcg_state >> 11;
}

static void fail(const char *label)
{
  fprintf(stderr, "unittest_symmetry_correlation: %s\n", label);
  exit(1);
}

static unsigned int spin_max(const struct DefineList *def)
{
  return def->iCalcModel == SpinlessFermion ? 0U : 1U;
}

static void test_matches_term(int model, unsigned int nsite, const char *label)
{
  struct DefineList def;
  unsigned int trial;
  unsigned int width = (model == Hubbard || model == tJ) ? 2U : 1U;
  memset(&def, 0, sizeof(def));
  def.iCalcModel = model;
  def.Nsite = nsite;
  for (trial = 0; trial < 1000U; ++trial) {
    struct SymmetryTerm term;
    unsigned long limit = 1UL << (nsite * width);
    unsigned long state = lcg_next() & (limit - 1UL);
    unsigned long out_term = 0UL, out_factors = 0UL;
    double complex value = 0.0;
    double sign = 0.0;
    unsigned int f;
    int status_term, status_factors;
    memset(&term, 0, sizeof(term));
    term.factors = 1U + (unsigned int)(lcg_next() % 2U);
    term.value = 0.37 - 0.21 * I;
    for (f = 0; f < term.factors; ++f) {
      int site_out = (int)(lcg_next() % nsite);
      term.index[4U*f] = site_out;
      term.index[4U*f+1U] = (int)(lcg_next() % (spin_max(&def) + 1U));
      term.index[4U*f+2U] = model == Spin ? site_out : (int)(lcg_next() % nsite);
      term.index[4U*f+3U] = (int)(lcg_next() % (spin_max(&def) + 1U));
    }
    status_term = ApplySymmetryTerm(&def, &term, state, &out_term, &value);
    status_factors = ApplySymmetryFactors(&def, term.factors, term.index,
                                          state, &out_factors, &sign);
    if (status_term != status_factors) fail(label);
    if (status_term == 1 &&
        (out_term != out_factors || value != sign * term.value)) fail(label);
  }
}

static void test_three_factors(void)
{
  struct DefineList def;
  int index[12] = {1,0,1,1, 0,1,0,0, 1,1,1,1};
  unsigned long out = 0UL;
  double sign = 0.0;
  memset(&def, 0, sizeof(def));
  def.iCalcModel = Hubbard;
  def.Nsite = 2U;
  if (ApplySymmetryFactors(&def, 3U, index, 9UL, &out, &sign) != 1 ||
      out != 6UL || sign != 1.0) fail("three-factor Hubbard product");
  {
    int reordered[12] = {1,1,1,1, 1,0,1,1, 0,1,0,0};
    if (ApplySymmetryFactors(&def, 3U, reordered, 9UL, &out, &sign) != 0)
      fail("factor order");
  }
  {
    int tj[8] = {1,0,0,0, 0,1,1,1};
    int final_double[4] = {0,1,1,1};
    def.iCalcModel = tJ;
    if (ApplySymmetryFactors(&def, 2U, tj, 9UL, &out, &sign) != 1 || out != 6UL)
      fail("tJ intermediate double occupancy");
    if (ApplySymmetryFactors(&def, 1U, final_double, 9UL, &out, &sign) != 0)
      fail("tJ final projection");
  }
  {
    int repeated[4U * 65U];
    int once[4] = {0,0,0,0};
    unsigned long out_once = 0UL;
    double sign_once = 0.0;
    unsigned int k;
    def.iCalcModel = Hubbard;
    for (k = 0; k < 65U; ++k) {
      repeated[4U*k] = 0;
      repeated[4U*k+1U] = 0;
      repeated[4U*k+2U] = 0;
      repeated[4U*k+3U] = 0;
    }
    if (ApplySymmetryFactors(&def, 65U, repeated, 9UL, &out, &sign) !=
        ApplySymmetryFactors(&def, 1U, once, 9UL, &out_once, &sign_once) ||
        out != out_once || sign != sign_once) fail("65-factor density");
    if (ApplySymmetryFactors(&def, UINT_MAX / 4U + 1U, repeated, 9UL,
                             &out, &sign) != -1) fail("factor overflow");
  }
}

static void assert_close(double complex got, double complex expected,
                         double tolerance, const char *label)
{
  if (cabs(got - expected) > tolerance) {
    fprintf(stderr, "%s: got (%.15g,%.15g), expected (%.15g,%.15g)\n",
            label, creal(got), cimag(got), creal(expected), cimag(expected));
    exit(1);
  }
}

static unsigned int def_spin_max(const struct DefineList *def)
{
  return def->iCalcModel == SpinlessFermion ? 0U : 1U;
}

static unsigned long int setup_tj_basis(unsigned int nsite,
                                        unsigned int nup,
                                        unsigned int ndown)
{
  unsigned long int state, dim = 0;
  unsigned long int limit = 1UL << (2U * nsite);
  free(list_1);
  free(list_Diagonal);
  list_1 = NULL;
  list_Diagonal = NULL;
  for (state = 0; state < limit; ++state)
    if ((state & (state >> 1U) & (ULONG_MAX / 3UL)) == 0UL &&
        (unsigned int)count_hubbard_spin(state, nsite, 0) == nup &&
        (unsigned int)count_hubbard_spin(state, nsite, 1) == ndown) ++dim;
  list_1 = (unsigned long int *)calloc(dim + 1UL, sizeof(*list_1));
  list_Diagonal = (double *)calloc(dim + 1UL, sizeof(*list_Diagonal));
  if (list_1 == NULL || list_Diagonal == NULL) {
    fprintf(stderr, "failed to allocate tJ basis\n");
    exit(1);
  }
  dim = 0;
  for (state = 0; state < limit; ++state)
    if ((state & (state >> 1U) & (ULONG_MAX / 3UL)) == 0UL &&
        (unsigned int)count_hubbard_spin(state, nsite, 0) == nup &&
        (unsigned int)count_hubbard_spin(state, nsite, 1) == ndown)
      list_1[++dim] = state;
  test_raw_dim = dim;
  return dim;
}

static void setup_tj_bind(struct BindStruct *X,
                          unsigned int nsite,
                          unsigned int nup,
                          unsigned int ndown,
                          unsigned int momentum_index)
{
  setup_hubbard_bind(X, nsite, nup, ndown, momentum_index);
  X->Def.iCalcModel = tJ;
  X->Check.idim_max = setup_tj_basis(nsite, nup, ndown);
}
/* 参照: sector vector を raw 基底へ展開し、applier で直接 <O> を取る。
 * 軌道グループ化・canonicalize・directory・halo を一切使わないので core の独立検証になる
 * （applier 自体は共有する。符号の独立検証は python 回帰の物理基底 builder が担う）。 */
static double complex reference_expectation(const struct BindStruct *X,
                                            const double complex *a,      /* a[1..dim] */
                                            unsigned int factors,
                                            const int *index)
{
  unsigned long int raw, r;
  unsigned long int dim = X->Sym->dim;
  double complex *c = (double complex *)calloc(test_raw_dim + 1UL, sizeof(*c));
  double complex value = 0.0;
  if (c == NULL) exit(1);
  for (raw = 1; raw <= test_raw_dim; raw++)
    for (r = 1; r <= dim; r++)
      c[raw] += a[r] * projected_coeff_for_state(X, r, list_1[raw]);
  for (raw = 1; raw <= test_raw_dim; raw++) {
    unsigned long int out = 0UL, out_index;
    double sign = 0.0;
    if (ApplySymmetryFactors(&X->Def, factors, index, list_1[raw], &out, &sign) != 1) continue;
    out_index = raw_index_for_state(out);           /* 0 if out is outside the raw list */
    if (out_index == 0UL) continue;
    value += c[raw] * sign * conj(c[out_index]);
  }
  free(c);
  return value;
}

static void random_sector_vector(double complex *a, unsigned long int dim)
{
  unsigned long int r;
  double norm = 0.0;
  for (r = 1; r <= dim; r++) {
    a[r] = ((double)(lcg_next() % 2001U) - 1000.0) / 1000.0 + I * (((double)(lcg_next() % 2001U) - 1000.0) / 1000.0);
    norm += creal(a[r]) * creal(a[r]) + cimag(a[r]) * cimag(a[r]);
  }
  norm = sqrt(norm);
  for (r = 1; r <= dim; r++) a[r] /= norm;
}

static void activate(struct BindStruct *X, const char *label)
{
  if (BuildSymmetryBasis(X) != 0 || ActivateSymmetryBasisDimension(X) != 0) {
    fprintf(stderr, "%s: basis build failed\n", label);
    exit(1);
  }
  X->Check.idim_max = X->Sym->local_dim;
}

/* 演算子集合: OneBodyG 全 pair、密度・交換型 TwoBodyG、3 因子 1 行、spinful では sector 外へ出る行 1 行。 */
static size_t build_request(const struct DefineList *def, int *index, struct SymmetryCorrelationOperator *ops)
{
  size_t n = 0; unsigned int i, j;
  unsigned int spin_max = def_spin_max(def);
  int *p = index;
  for (i = 0; i < def->Nsite; i++) for (j = 0; j < def->Nsite; j++) {
    unsigned int s;
    if ((def->iCalcModel == Spin || def->iCalcModel == SpinGC) && i != j) continue;
    for (s = 0; s <= spin_max; s++) {
      p[0] = (int)i; p[1] = (int)s; p[2] = (int)j; p[3] = (int)s;
      ops[n].factors = 1U; ops[n].index = p; n++; p += 4;
    }
  }
  for (i = 0; i < def->Nsite; i++) {
    j = (i + 1U) % def->Nsite;
    p[0] = (int)i; p[1] = 0; p[2] = (int)i; p[3] = 0; p[4] = (int)j; p[5] = 0; p[6] = (int)j; p[7] = 0;
    ops[n].factors = 2U; ops[n].index = p; n++; p += 8;
    if (spin_max == 1U) {
      p[0] = (int)i; p[1] = 0; p[2] = (int)i; p[3] = 1; p[4] = (int)j; p[5] = 1; p[6] = (int)j; p[7] = 0;
      ops[n].factors = 2U; ops[n].index = p; n++; p += 8;
    }
  }
  p[0] = 0; p[1] = 0; p[2] = 0; p[3] = 0;
  p[4] = 1; p[5] = 0;
  p[6] = (def->iCalcModel == Spin || def->iCalcModel == SpinGC) ? 1 : 2; p[7] = 0;
  p[8] = 2; p[9] = 0; p[10] = 2; p[11] = 0;
  ops[n].factors = 3U; ops[n].index = p; n++; p += 12;
  if (spin_max == 1U) {
    p[0] = 0; p[1] = 0; p[2] = 0; p[3] = 1;   /* S+_0 (Spin) or c^dag_{0,up} c_{0,down}: leaves the sector */
    ops[n].factors = 1U; ops[n].index = p; n++; p += 4;
  }
  if (def->iCalcModel == SpinGC) {
    static const int four[] = {0,1,0,0, 0,0,0,1, 1,1,1,0, 2,0,2,1};
    static const int six[] = {0,1,0,0, 0,0,0,1, 1,1,1,0,
                              1,0,1,1, 2,1,2,0, 3,0,3,1};
    memcpy(p, four, sizeof(four));
    ops[n].factors = 4U; ops[n].index = p; n++; p += 16;
    memcpy(p, six, sizeof(six));
    ops[n].factors = 6U; ops[n].index = p; n++;
  }
  return n;
}

static void test_expectation(struct BindStruct *X, size_t expected_orbits, unsigned long long min_waves,
                             unsigned long long exact_waves, unsigned int min_threads,
                             const char *label)
{
  int index[4096];
  struct SymmetryCorrelationOperator ops[512];
  double complex values[512];
  struct SymmetryCorrelationStats stats;
  double complex *a;
  size_t n, t, orbits = 0, members = 0;
  activate(X, label);
  a = (double complex *)calloc(X->Sym->dim + 1UL, sizeof(*a));
  if (a == NULL) exit(1);
  random_sector_vector(a, X->Sym->dim);
  n = build_request(&X->Def, index, ops);
  if (SymmetryCorrelationOrbitCount(&X->Def, ops, n, &orbits, &members) != 0 ||
      (expected_orbits != (size_t)-1 && orbits != expected_orbits)) {
    fprintf(stderr, "%s: orbit count %zu, expected %zu\n", label, orbits, expected_orbits);
    exit(1);
  }
  memset(&stats, 0, sizeof(stats));
  if (SymmetryCorrelationExpectationWithStats(X, a, ops, n, values, &stats) != 0) {
    fprintf(stderr, "%s: expectation failed\n", label);
    exit(1);
  }
  for (t = 0; t < n; t++)
    assert_close(values[t], reference_expectation(X, a, ops[t].factors, ops[t].index), 1.0e-12, label);
  if (def_spin_max(&X->Def) == 1U && X->Def.iCalcModel != SpinGC)
    assert_close(values[n-1], 0.0, 1.0e-15, label);
  if (stats.waves < min_waves || (exact_waves != 0ULL && stats.waves != exact_waves)) {
    fprintf(stderr, "%s: waves=%llu (min %llu, exact %llu)\n", label, stats.waves, min_waves, exact_waves);
    exit(1);
  }
  if (stats.threads < min_threads) { fprintf(stderr, "%s: threads=%u\n", label, stats.threads); exit(1); }
  if (X->Def.iCalcModel == SpinGC) {
    int offsite[4] = {0, 1, 1, 0};
    struct SymmetryCorrelationOperator invalid = {1U, offsite};
    size_t invalid_orbits = 0, invalid_members = 0;
    if (SymmetryCorrelationOrbitCount(&X->Def, &invalid, 1U,
                                      &invalid_orbits, &invalid_members) != -1)
      fail("SpinGC offsite OrbitCount must fail");
    if (SymmetryCorrelationExpectation(X, a, &invalid, 1U, values) != -1)
      fail("SpinGC offsite expectation must fail");
  }
  /* count == 0 は何もしない */
  values[0] = 42.0;
  if (SymmetryCorrelationExpectation(X, a, ops, 0, values) != 0 || values[0] != 42.0) {
    fprintf(stderr, "%s: count 0 must be a no-op\n", label);
    exit(1);
  }
  /* 有限上限: 軌道表すら入らない 256 byte では全 rank 一致で -1、0 で無制限に戻る。 */
  setenv("HPHI_SYMMETRY_MEMORY_LIMIT_BYTES", "256", 1);
  if (SymmetryCorrelationExpectation(X, a, ops, n, values) != -1) {
    fprintf(stderr, "%s: a 256-byte hard limit must fail\n", label);
    exit(1);
  }
  /* 対角専用の要求も policy を迂回しない。単一 orbit が 256 byte
   * 未満になる ABI もあるので、ここでは確実に不足する 1 byte を使う。 */
  setenv("HPHI_SYMMETRY_MEMORY_LIMIT_BYTES", "1", 1);
  if (SymmetryCorrelationExpectation(X, a, ops, 1, values) != -1) {
    fprintf(stderr, "%s: the diagonal-only path must honour the hard limit\n", label);
    exit(1);
  }
  setenv("HPHI_SYMMETRY_MEMORY_LIMIT_BYTES", "0", 1);
  if (SymmetryCorrelationExpectation(X, a, ops, 1, values) != 0) {
    fprintf(stderr, "%s: unlimited policy must succeed\n", label);
    exit(1);
  }
  free(a);
  FreeSymmetryBasis(X->Sym);
  X->Sym = NULL;
  free(list_1); list_1 = NULL;
  free(list_Diagonal); list_Diagonal = NULL;
}

int main(void)
{
  struct BindStruct X;
  unsigned int spingc_threads = 1U;
  stdoutMPI = stdout;
  test_matches_term(Spin, 6U, "spin wrapper equivalence");
  test_matches_term(SpinlessFermion, 6U, "spinless wrapper equivalence");
  test_matches_term(Hubbard, 4U, "Hubbard wrapper equivalence");
  test_matches_term(tJ, 4U, "tJ wrapper equivalence");
  test_three_factors();

  setup_bind(&X, 6U, 3U, 0U);
  test_expectation(&X, 6U, 1ULL, 0ULL, 1U, "spin k=0");
  setup_bind(&X, 6U, 3U, 1U);
  test_expectation(&X, 6U, 1ULL, 0ULL, 1U, "spin k=1 (complex character)");
  setup_spinless_bind(&X, 6U, 3U, 1U);
  test_expectation(&X, 8U, 1ULL, 0ULL, 1U, "spinless k=1");
  setup_hubbard_bind(&X, 4U, 2U, 2U, 2U);
  test_expectation(&X, 12U, 2ULL, 0ULL, 1U, "hubbard k=pi");
  setup_tj_bind(&X, 6U, 1U, 1U, 0U);
  test_expectation(&X, 16U, 5ULL, 5ULL, 1U, "tJ k=0");
  setup_spingc_bind(&X, 6U, 1U);
#ifdef _OPENMP
  omp_set_num_threads(3);
  spingc_threads = 3U;
#endif
  test_expectation(&X, (size_t)-1, 2ULL, 0ULL, spingc_threads,
                   "SpinGC k=1 OMP3");
  puts("unittest_symmetry_correlation: passed");
  return 0;
}
