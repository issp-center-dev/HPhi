#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "DefCommon.h"
#include "global.h"
struct BindStruct;
#include "struct.h"
#include "symmetry_kondo.h"
#include "symmetry_state_enumerator.h"
#include "symmetry_basis.h"

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


static unsigned int failures = 0U;

static void check_int(const char *label, int actual, int expected)
{
  if (actual != expected) {
    fprintf(stderr, "FAIL: %s: got %d, expected %d\n",
            label, actual, expected);
    failures++;
  }
}

static void check_ulong(const char *label, unsigned long actual,
                        unsigned long expected)
{
  if (actual != expected) {
    fprintf(stderr, "FAIL: %s: got %lu, expected %lu\n",
            label, actual, expected);
    failures++;
  }
}

static void check_unchanged(const char *label, const struct DefineList *actual,
                            const struct DefineList *expected)
{
  if (memcmp(actual, expected, sizeof(*actual)) != 0) {
    fprintf(stderr, "FAIL: %s modified the definition on failure\n", label);
    failures++;
  }
}

static struct DefineList base_normalize(int model, unsigned int nsite,
                                        unsigned int nlocal)
{
  struct DefineList def;
  memset(&def, 0, sizeof(def));
  def.iCalcModel = model;
  def.Nsite = nsite;
  def.NLocSpn = nlocal;
  return def;
}

static void test_model_classification(void)
{
  check_int("Kondo model", IsSymmetryKondoModel(Kondo), 1);
  check_int("KondoNConserved model",
            IsSymmetryKondoModel(KondoNConserved), 1);
  check_int("KondoGC model", IsSymmetryKondoModel(KondoGC), 1);
  check_int("Hubbard is not Kondo", IsSymmetryKondoModel(Hubbard), 0);
}

static void test_normalization_success(void)
{
  struct NormalizeCase {
    const char *label;
    struct DefineList input;
    int has_ncond;
    int has_sz;
    int has_nup;
    int has_ndown;
    int model;
    unsigned int ncond;
    unsigned int nup;
    unsigned int ndown;
    unsigned int ne;
    int total2sz;
    int sz_conserved;
  } cases[6];
  unsigned int i;

  memset(cases, 0, sizeof(cases));

  cases[0].label = "derive Ncond and 2Sz from Nup/Ndown";
  cases[0].input = base_normalize(Kondo, 6U, 3U);
  cases[0].input.Nup = 3U;
  cases[0].input.Ndown = 2U;
  cases[0].has_nup = cases[0].has_ndown = TRUE;
  cases[0].model = Kondo;
  cases[0].ncond = 2U;
  cases[0].nup = 3U;
  cases[0].ndown = 2U;
  cases[0].ne = 5U;
  cases[0].total2sz = 1;
  cases[0].sz_conserved = TRUE;

  cases[1].label = "derive Nup/Ndown from Ncond and 2Sz";
  cases[1].input = base_normalize(Kondo, 6U, 3U);
  cases[1].input.NCond = 4U;
  cases[1].input.Total2Sz = 1;
  cases[1].input.iFlgSzConserved = TRUE;
  cases[1].has_ncond = cases[1].has_sz = TRUE;
  cases[1].model = Kondo;
  cases[1].ncond = 4U;
  cases[1].nup = 4U;
  cases[1].ndown = 3U;
  cases[1].ne = 7U;
  cases[1].total2sz = 1;
  cases[1].sz_conserved = TRUE;

  cases[2].label = "Ncond-only vacuum with no local spins";
  cases[2].input = base_normalize(Kondo, 3U, 0U);
  cases[2].input.NCond = 0U;
  cases[2].input.Nup = 99U;
  cases[2].input.Ndown = 98U;
  cases[2].input.Total2Sz = -7;
  cases[2].input.iFlgSzConserved = TRUE;
  cases[2].has_ncond = TRUE;
  cases[2].model = KondoNConserved;
  cases[2].ncond = 0U;
  cases[2].ne = 0U;
  cases[2].sz_conserved = FALSE;

  cases[3].label = "explicit zero Nup/Ndown is a canonical sector";
  cases[3].input = base_normalize(Kondo, 3U, 0U);
  cases[3].input.Nup = 0U;
  cases[3].input.Ndown = 0U;
  cases[3].has_nup = cases[3].has_ndown = TRUE;
  cases[3].model = Kondo;
  cases[3].ncond = 0U;
  cases[3].ne = 0U;
  cases[3].sz_conserved = TRUE;

  cases[4].label = "spin-only space with no conduction sites";
  cases[4].input = base_normalize(Kondo, 3U, 3U);
  cases[4].input.Nup = 2U;
  cases[4].input.Ndown = 1U;
  cases[4].has_nup = cases[4].has_ndown = TRUE;
  cases[4].model = Kondo;
  cases[4].ncond = 0U;
  cases[4].nup = 2U;
  cases[4].ndown = 1U;
  cases[4].ne = 3U;
  cases[4].total2sz = 1;
  cases[4].sz_conserved = TRUE;

  cases[5].label = "KondoGC clears all fixed quantities";
  cases[5].input = base_normalize(KondoGC, 4U, 2U);
  cases[5].input.NCond = 3U;
  cases[5].input.Nup = 3U;
  cases[5].input.Ndown = 2U;
  cases[5].input.Ne = 5U;
  cases[5].input.Total2Sz = 1;
  cases[5].input.iFlgSzConserved = TRUE;
  cases[5].model = KondoGC;
  cases[5].sz_conserved = FALSE;

  for (i = 0U; i < sizeof(cases) / sizeof(cases[0]); ++i) {
    struct DefineList actual = cases[i].input;
    check_int(cases[i].label,
              NormalizeSymmetryKondoQuantumNumbers(
                  &actual, cases[i].has_ncond, cases[i].has_sz,
                  cases[i].has_nup, cases[i].has_ndown),
              0);
    check_int(cases[i].label, actual.iCalcModel, cases[i].model);
    check_ulong(cases[i].label, actual.NCond, cases[i].ncond);
    check_ulong(cases[i].label, actual.Nup, cases[i].nup);
    check_ulong(cases[i].label, actual.Ndown, cases[i].ndown);
    check_ulong(cases[i].label, actual.Ne, cases[i].ne);
    check_int(cases[i].label, actual.Total2Sz, cases[i].total2sz);
    check_int(cases[i].label, actual.iFlgSzConserved,
              cases[i].sz_conserved);
  }
}

static void expect_normalize_failure(const char *label,
                                     struct DefineList input,
                                     int has_ncond, int has_sz,
                                     int has_nup, int has_ndown)
{
  struct DefineList before = input;
  check_int(label, NormalizeSymmetryKondoQuantumNumbers(
                       &input, has_ncond, has_sz, has_nup, has_ndown),
            -1);
  check_unchanged(label, &input, &before);
}

static void test_normalization_failures(void)
{
  struct DefineList def;
  unsigned int flag;

  def = base_normalize(Kondo, 6U, 3U);
  def.Nup = 3U;
  def.Ndown = 2U;
  def.NCond = 4U;
  expect_normalize_failure("explicit Ncond mismatch", def,
                           TRUE, FALSE, TRUE, TRUE);

  def = base_normalize(Kondo, 6U, 3U);
  def.Nup = 3U;
  def.Ndown = 2U;
  def.Total2Sz = 0;
  def.iFlgSzConserved = TRUE;
  expect_normalize_failure("explicit 2Sz mismatch", def,
                           FALSE, TRUE, TRUE, TRUE);

  def = base_normalize(Kondo, 6U, 3U);
  def.NCond = 2U;
  def.Total2Sz = 0;
  def.iFlgSzConserved = TRUE;
  expect_normalize_failure("Ncond/2Sz parity mismatch", def,
                           TRUE, TRUE, FALSE, FALSE);

  def = base_normalize(Kondo, 6U, 3U);
  def.Nup = 3U;
  expect_normalize_failure("only Nup specified", def,
                           FALSE, FALSE, TRUE, FALSE);

  def = base_normalize(Kondo, 6U, 3U);
  expect_normalize_failure("no canonical quantities specified", def,
                           FALSE, FALSE, FALSE, FALSE);

  def = base_normalize(Kondo, 6U, 3U);
  def.Total2Sz = -1;
  def.iFlgSzConserved = TRUE;
  expect_normalize_failure("2Sz without Ncond", def,
                           FALSE, TRUE, FALSE, FALSE);

  def = base_normalize(Kondo, 3U, 2U);
  def.NCond = 3U;
  expect_normalize_failure("Ncond exceeds twice conduction sites", def,
                           TRUE, FALSE, FALSE, FALSE);

  def = base_normalize(Kondo, 4U, 3U);
  def.Nup = 0U;
  def.Ndown = 0U;
  expect_normalize_failure("inferred negative Ncond and L greater than Ne", def,
                           FALSE, FALSE, TRUE, TRUE);

  def = base_normalize(Kondo, 3U, 0U);
  def.NCond = UINT_MAX;
  expect_normalize_failure("wrapped negative Ncond is rejected", def,
                           TRUE, FALSE, FALSE, FALSE);

  def = base_normalize(Kondo, UINT_MAX, 0U);
  def.Nup = UINT_MAX;
  def.Ndown = 0U;
  expect_normalize_failure("derived 2Sz must fit its signed field", def,
                           FALSE, FALSE, TRUE, TRUE);

  def = base_normalize(Kondo, 2U, 3U);
  def.NCond = 0U;
  expect_normalize_failure("more local spins than sites", def,
                           TRUE, FALSE, FALSE, FALSE);

  def = base_normalize(Kondo, 0U, 0U);
  def.NCond = 0U;
  expect_normalize_failure("zero sites", def,
                           TRUE, FALSE, FALSE, FALSE);

  def = base_normalize(Hubbard, 4U, 0U);
  def.NCond = 2U;
  expect_normalize_failure("non-Kondo model", def,
                           TRUE, FALSE, FALSE, FALSE);

  for (flag = 0U; flag < 4U; ++flag) {
    def = base_normalize(KondoGC, 4U, 2U);
    def.NCond = 0U;
    def.Nup = 0U;
    def.Ndown = 0U;
    def.Total2Sz = 0;
    expect_normalize_failure("KondoGC rejects explicit zero", def,
                             flag == 0U, flag == 1U,
                             flag == 2U, flag == 3U);
  }
}

static struct DefineList base_physical(int model, int *locspn)
{
  struct DefineList def = base_normalize(model, 4U, 2U);
  def.LocSpn = locspn;
  def.NCond = 2U;
  def.Ne = 4U;
  if (model == Kondo) {
    def.Nup = 2U;
    def.Ndown = 2U;
    def.iFlgSzConserved = TRUE;
  }
  return def;
}

static void test_physical_space_and_mask(void)
{
  int locspn[4] = {LOCSPIN, ITINERANT, LOCSPIN, ITINERANT};
  struct DefineList def = base_physical(Kondo, locspn);
  unsigned long mask = ULONG_MAX;
  unsigned long physical = 1UL | (3UL << 2U) | (2UL << 4U);

  check_int("validate canonical Kondo space",
            ValidateSymmetryKondoSpace(&def), 0);
  check_int("local mask", SymmetryKondoLocalMask(&def, &mask), 0);
  check_ulong("local mask bits", mask, 0x5UL);
  check_int("physical mixed Kondo state",
            SymmetryKondoStateIsPhysical(&def, physical), 1);
  check_int("empty local digit is unphysical",
            SymmetryKondoStateIsPhysical(&def, physical & ~3UL), 0);
  check_int("double local digit is unphysical",
            SymmetryKondoStateIsPhysical(&def, physical | 2UL), 0);
  check_int("bits beyond the lattice are unphysical",
            SymmetryKondoStateIsPhysical(&def,
                                         physical | (1UL << 12U)), 0);

  def = base_physical(KondoNConserved, locspn);
  check_int("validate KondoNConserved space",
            ValidateSymmetryKondoSpace(&def), 0);

  def = base_physical(KondoGC, locspn);
  def.NCond = 0U;
  def.Ne = 0U;
  check_int("validate KondoGC space", ValidateSymmetryKondoSpace(&def), 0);
}

static void test_physical_space_failures(void)
{
  int locspn[4] = {LOCSPIN, ITINERANT, LOCSPIN, ITINERANT};
  int negative[4] = {LOCSPIN, -1, LOCSPIN, ITINERANT};
  int general_spin[4] = {2, ITINERANT, LOCSPIN, ITINERANT};
  int wide_locspn[sizeof(unsigned long) * CHAR_BIT / 2U + 1U];
  struct DefineList def = base_physical(Kondo, locspn);
  unsigned int i;
  unsigned long mask = 17UL;

  def.iCalcModel = Hubbard;
  check_int("reject non-Kondo physical space",
            ValidateSymmetryKondoSpace(&def), -1);

  def = base_physical(Kondo, NULL);
  check_int("reject missing LocSpn", ValidateSymmetryKondoSpace(&def), -1);

  def = base_physical(Kondo, locspn);
  def.NLocSpn = 1U;
  check_int("reject mismatched NLocSpn",
            ValidateSymmetryKondoSpace(&def), -1);
  check_int("mask rejects invalid space", SymmetryKondoLocalMask(&def, &mask),
            -1);
  check_ulong("mask unchanged on failure", mask, 17UL);

  def = base_physical(Kondo, negative);
  check_int("reject negative LocSpn", ValidateSymmetryKondoSpace(&def), -1);

  def = base_physical(Kondo, general_spin);
  check_int("reject general local spin",
            ValidateSymmetryKondoSpace(&def), -1);

  for (i = 0U; i < sizeof(wide_locspn) / sizeof(wide_locspn[0]); ++i)
    wide_locspn[i] = ITINERANT;
  def = base_normalize(Kondo, (unsigned int)(sizeof(unsigned long) * CHAR_BIT / 2U), 0U);
  def.LocSpn = wide_locspn;
  def.NCond = 0U;
  def.Ne = 0U;
  def.iFlgSzConserved = TRUE;
  check_int("canonical Kondo accepts the top word bit",
            ValidateSymmetryKondoSpace(&def), 0);

  def.iCalcModel = KondoGC;
  def.iFlgSzConserved = FALSE;
  check_int("KondoGC reserves the top word bit",
            ValidateSymmetryKondoSpace(&def), -1);

  def.iCalcModel = Kondo;
  def.Nsite++;
  check_int("canonical Kondo rejects more than one word",
            ValidateSymmetryKondoSpace(&def), -1);
}

static void test_permutation_sign(void)
{
  int locspn[4] = {LOCSPIN, ITINERANT, LOCSPIN, ITINERANT};
  int identity[4] = {0, 1, 2, 3};
  int swap_pairs[4] = {2, 3, 0, 1};
  int duplicate[4] = {0, 1, 0, 3};
  int out_of_range[4] = {0, 1, 2, 4};
  int changes_kind[4] = {1, 0, 2, 3};
  struct DefineList def = base_physical(Kondo, locspn);
  int sign = 0;

  check_int("identity permutation",
            SymmetryKondoPermutationSign(&def, identity, &sign), 0);
  check_int("identity local sign", sign, 1);

  check_int("local transposition",
            SymmetryKondoPermutationSign(&def, swap_pairs, &sign), 0);
  check_int("odd local sign", sign, -1);

  sign = 7;
  check_int("reject duplicate permutation image",
            SymmetryKondoPermutationSign(&def, duplicate, &sign), -1);
  check_int("sign unchanged for duplicate image", sign, 7);
  check_int("reject out-of-range permutation image",
            SymmetryKondoPermutationSign(&def, out_of_range, &sign), -1);
  check_int("reject local/conduction mixing",
            SymmetryKondoPermutationSign(&def, changes_kind, &sign), -1);
}


/* Independent word filter: deleting the local occupancy restriction or changing
 * the unranking order must fail the complete word-by-word comparison. */
static int brute_accept(const struct DefineList *d, unsigned long word)
{
  unsigned int up = 0, down = 0, nc = 0, i;
  for (i = 0; i < d->Nsite; ++i) {
    unsigned int digit = (word >> (2 * i)) & 3UL;
    if (d->LocSpn[i] && digit != 1 && digit != 2)
      return 0;
    up += digit & 1U;
    down += (digit >> 1) & 1U;
    if (!d->LocSpn[i])
      nc += (digit & 1U) + ((digit >> 1) & 1U);
  }
  return d->iCalcModel == KondoGC ||
         (d->iCalcModel == KondoNConserved ? nc == d->NCond
                                           : up == d->Nup && down == d->Ndown);
}

static struct DefineList enumeration_case(int model, unsigned int n, unsigned int l,
                                          int *loc, unsigned int nc, int sz)
{
  struct DefineList d = base_normalize(model, n, l);
  d.LocSpn = loc;
  if (model != KondoGC) {
    d.NCond = nc;
    d.Ne = l + nc;
  }
  if (model == Kondo) {
    d.Nup = (d.Ne + sz) / 2;
    d.Ndown = d.Ne - d.Nup;
    d.Total2Sz = sz;
    d.iFlgSzConserved = TRUE;
  }
  return d;
}

static void compare_enumeration(struct DefineList *d, unsigned long want)
{
  struct SymmetryStateEnumerator e;
  unsigned long dim = 0, word, got = 0, rank = 0;
  check_int("checked dimension", ComputeSymmetryKondoDimension(d, &dim), 0);
  check_ulong("dimension fixture", dim, want);
  if (InitSymmetryStateEnumerator(d, want, &e)) {
    check_int("initialize Kondo enumerator", -1, 0);
    return;
  }
  for (word = 0; word < (1UL << (2 * d->Nsite)); ++word)
    if (brute_accept(d, word)) {
      ++rank;
      check_int("direct unranking", SymmetryStateEnumeratorStateAt(&e, rank, &got), 0);
      check_ulong("word order", got, word);
    }
  check_ulong("brute-force dimension", rank, want);
  check_int("rank zero", SymmetryStateEnumeratorStateAt(&e, 0, &got), -1);
  check_int("rank beyond dimension", SymmetryStateEnumeratorStateAt(&e, want + 1, &got),
            -1);
  check_int("expected dimension mismatch", InitSymmetryStateEnumerator(d, want + 1, &e),
            -1);
}

static void test_enumeration(void)
{
  unsigned int p, layout, m, i;
  int models[] = {Kondo, KondoNConserved, KondoGC};
  unsigned long dims[2][3] = {{39, 120, 512}, {144, 448, 4096}};
  for (p = 3; p <= 4; ++p)
    for (layout = 0; layout < 2; ++layout) {
      int loc[8];
      for (i = 0; i < 2 * p; ++i)
        loc[i] = layout ? i % 2 == 0 : i < p;
      for (m = 0; m < 3; ++m) {
        struct DefineList d =
            enumeration_case(models[m], 2 * p, p, loc, 2, p == 3 ? 1 : 0);
        compare_enumeration(&d, dims[p - 3][m]);
      }
    }
}

static void test_completion_counts(void)
{
  unsigned int l, c, m, i;
  int models[] = {Kondo, KondoNConserved, KondoGC};
  for (l = 0; l <= 3; ++l)
    for (c = 0; c <= 3; ++c)
      for (m = 0; m < 3; ++m) {
        int up, down, nc, loc[6];
        for (i = 0; i < l + c; ++i)
          loc[i] = i < l;
        for (up = -1; up <= (int)(l + c) + 1; ++up)
          for (down = -1; down <= (int)(l + c) + 1; ++down) {
            unsigned long count = 99, want = 0, w;
            nc = up + down - (int)l;
            struct DefineList d = enumeration_case(models[m], l + c, l, loc,
                                                   nc < 0 ? 0 : (unsigned int)nc, 0);
            d.Nup = up;
            d.Ndown = down;
            d.NCond = nc;
            for (w = 0; w < (1UL << (2 * (l + c))); ++w)
              want += brute_accept(&d, w);
            check_int(
                "completion status",
                CountSymmetryKondoCompletions(models[m], l, c, up, down, nc, &count), 0);
            check_ulong("completion vs independent filter", count, want);
            if (l + c && want) {
              d.Ne = l + nc;
              d.Total2Sz = up - down;
              if (models[m] != Kondo)
                d.Nup = d.Ndown = d.Total2Sz = 0;
              if (models[m] == KondoGC)
                d.Ne = d.NCond = 0;
              compare_enumeration(&d, want);
            }
            if (models[m] == Kondo && up >= 0 && down >= 0 && l + c) {
              struct DefineList normalized = base_normalize(Kondo, l + c, l);
              normalized.Nup = up;
              normalized.Ndown = down;
              check_int("existence inequality vs completion",
                        NormalizeSymmetryKondoQuantumNumbers(&normalized, 0, 0, 1, 1) ==
                            0,
                        want != 0);
            }
          }
      }
}

static void test_word_boundaries(void)
{
  const unsigned int n = HPHI_SYMMETRY_STATE_WORD_BITS / 2;
  int loc[HPHI_SYMMETRY_STATE_WORD_BITS / 2];
  unsigned int i;
  unsigned long dim = 0, word = 0, want = 0;
  struct SymmetryStateEnumerator e;
  for (i = 0; i < n; ++i) {
    loc[i] = 1;
    want |= 1UL << (2 * i);
  }
  struct DefineList d = enumeration_case(Kondo, n, n, loc, 0, n);
  check_int("top-bit canonical init", InitSymmetryStateEnumerator(&d, 1, &e), 0);
  check_int("top-bit canonical state", SymmetryStateEnumeratorStateAt(&e, 1, &word), 0);
  check_ulong("all-up canonical word", word, want);
  for (i = 0; i < n; ++i)
    loc[i] = 0;
  memset(&d, 0, sizeof(d));
  check_int("enumerator owns its layout", SymmetryStateEnumeratorStateAt(&e, 1, &word),
            0);
  check_ulong("state independent of original def and LocSpn", word, want);
  d = enumeration_case(KondoNConserved, n, 0, loc, 2 * n, 0);
  check_int("full conduction init", InitSymmetryStateEnumerator(&d, 1, &e), 0);
  check_int("full conduction state", SymmetryStateEnumeratorStateAt(&e, 1, &word), 0);
  check_ulong("full conduction word", word, ULONG_MAX);
  check_int("completion product overflow",
            CountSymmetryKondoCompletions(KondoNConserved, n, n, 0, 0, n, &dim), -1);
  check_int(
      "completion sum overflow",
      CountSymmetryKondoCompletions(Kondo, n / 16 + 3, n, n / 2 + 2, n / 2 + 2, 0, &dim),
      -1);
  check_int("completion power overflow",
            CountSymmetryKondoCompletions(KondoGC, n, n, 0, 0, 0, &dim), -1);
  for (i = 0; i < n; ++i)
    loc[i] = 1;
  d = enumeration_case(KondoGC, n - 1, n - 1, loc, 0, 0);
  dim = 1UL << (n - 1);
  want &= ~(1UL << (2 * (n - 1)));
  check_int("large GC init without enumeration", InitSymmetryStateEnumerator(&d, dim, &e),
            0);
  check_int("large GC first", SymmetryStateEnumeratorStateAt(&e, 1, &word), 0);
  check_ulong("large GC first word", word, want);
  check_int("large GC last", SymmetryStateEnumeratorStateAt(&e, dim, &word), 0);
  check_ulong("large GC last word", word, want << 1);
  d.Nsite = d.NLocSpn = n;
  check_int("GC strict boundary", ComputeSymmetryKondoDimension(&d, &dim), -1);
}

static void test_physical_group_action(void)
{
  unsigned int p, layout, g, h, i;
  for (p = 3; p <= 4; ++p)
    for (layout = 0; layout < 2; ++layout) {
      int loc[8], perm[4][8], *rows[4];
      unsigned long w, vacuum = 0;
      for (i = 0; i < 2 * p; ++i) {
        loc[i] = layout ? i % 2 == 0 : i < p;
        if (loc[i])
          vacuum |= 1UL << (2 * i);
      }
      for (g = 0; g < p; ++g) {
        rows[g] = perm[g];
        for (i = 0; i < 2 * p; ++i)
          perm[g][i] =
              layout ? 2 * ((i / 2 + g) % p) + i % 2 : ((i % p + g) % p) + (i / p) * p;
      }
      struct DefineList d = enumeration_case(KondoGC, 2 * p, p, loc, 0, 0);
      d.NSymTrans = p;
      d.SymTrans = rows;
      for (w = 0; w < (1UL << (4 * p)); ++w)
        if (brute_accept(&d, w)) {
          for (g = 0; g < p; ++g)
            for (h = 0; h < p; ++h) {
              struct SymmetryTransformResult a, b, c;
              if (SymmetryApplyToState(&d, w, h, &a) ||
                  SymmetryApplyToState(&d, a.state, g, &b) ||
                  SymmetryApplyToState(&d, w, (g + h) % p, &c)) {
                check_int("Kondo physical group action supported", -1, 0);
                return;
              }
              check_ulong("group composition state", b.state, c.state);
              check_int("group composition sign", a.amplitude * b.amplitude, c.amplitude);
              if (w == vacuum)
                check_int("polarized vacuum momentum zero", a.amplitude, 1);
            }
        }
    }
}

/* Character projection trace is an independent sector-dimension calculation;
 * it detects a wrong local parity even when a wrong action is still a group. */
static void test_sector_dimensions(void)
{
  const int models[3] = {Kondo, KondoNConserved, KondoGC};
  const unsigned long raw[2][3] = {{39, 120, 512}, {144, 448, 4096}};
  const unsigned long expected[2][3][4] = {
      {{13, 13, 13, 0}, {40, 40, 40, 0}, {176, 168, 168, 0}},
      {{36, 36, 36, 36}, {108, 116, 108, 116}, {1024, 1024, 1024, 1024}}};
  unsigned int p, layout, m, g, k, i;
  for (p = 3; p <= 4; ++p)
    for (layout = 0; layout < 2; ++layout) {
      int loc[8], storage[4][8], *rows[4];
      for (i = 0; i < 2 * p; ++i)
        loc[i] = layout ? i % 2 == 0 : i < p;
      for (g = 0; g < p; ++g) {
        rows[g] = storage[g];
        for (i = 0; i < 2 * p; ++i)
          storage[g][i] =
              layout ? 2 * ((i / 2 + g) % p) + i % 2 : ((i % p + g) % p) + (i / p) * p;
      }
      for (m = 0; m < 3; ++m) {
        struct DefineList d =
            enumeration_case(models[m], 2 * p, p, loc, 2, p == 3 ? 1 : 0);
        struct SymmetryStateEnumerator e;
        double traces[4] = {0};
        unsigned long rank, word;
        d.NSymTrans = p;
        d.SymTrans = rows;
        if (InitSymmetryStateEnumerator(&d, raw[p - 3][m], &e)) {
          check_int("sector trace enumeration", -1, 0);
          continue;
        }
        for (rank = 1; rank <= e.raw_dim; ++rank) {
          check_int("trace unranking", SymmetryStateEnumeratorStateAt(&e, rank, &word),
                    0);
          for (g = 0; g < p; ++g) {
            struct SymmetryTransformResult moved;
            check_int("trace physical translation",
                      SymmetryApplyToState(&d, word, g, &moved), 0);
            if (moved.state == word)
              traces[g] += moved.amplitude;
          }
        }
        for (k = 0; k < p; ++k) {
          double complex dimension = 0;
          for (g = 0; g < p; ++g)
            dimension += traces[g] * cexp(2 * I * acos(-1) * k * g / p) / p;
          check_int("literal sector dimension",
                    cabs(dimension - expected[p - 3][m][k]) < 1e-10, 1);
        }
      }
    }
}

int main(void)
{
  stdoutMPI = tmpfile();
  if (stdoutMPI == NULL) stdoutMPI = stderr;

  test_sector_dimensions();
  test_enumeration();
  test_completion_counts();
  test_word_boundaries();
  test_physical_group_action();
  test_model_classification();
  test_normalization_success();
  test_normalization_failures();
  test_physical_space_and_mask();
  test_physical_space_failures();
  test_permutation_sign();

  if (stdoutMPI != stderr) fclose(stdoutMPI);
  if (failures != 0U) {
    fprintf(stderr, "%u symmetry Kondo unit checks failed\n", failures);
    return 1;
  }
  puts("symmetry Kondo unit checks passed");
  return 0;
}
