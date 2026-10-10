#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "DefCommon.h"
#include "global.h"
struct BindStruct;
#include "struct.h"
#include "symmetry_kondo.h"

FILE *stdoutMPI = NULL;

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

int main(void)
{
  stdoutMPI = tmpfile();
  if (stdoutMPI == NULL) stdoutMPI = stderr;

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
