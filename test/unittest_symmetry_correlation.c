#include <complex.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Common.h"
#include "symmetry_terms.h"

FILE *stdoutMPI;

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

int main(void)
{
  stdoutMPI = stdout;
  test_matches_term(Spin, 6U, "spin wrapper equivalence");
  test_matches_term(SpinlessFermion, 6U, "spinless wrapper equivalence");
  test_matches_term(Hubbard, 4U, "Hubbard wrapper equivalence");
  test_matches_term(tJ, 4U, "tJ wrapper equivalence");
  test_three_factors();
  puts("unittest_symmetry_correlation: task 1 passed");
  return 0;
}
