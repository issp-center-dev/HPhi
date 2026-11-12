#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "DefCommon.h"
struct BindStruct;
#include "struct.h"
#include "symmetry_diagonal.h"
#include "symmetry_terms.h"

FILE *stdoutMPI = NULL;
int nproc = 1;
int myrank = 0;

struct MatrixContext {
  const struct DefineList *def;
  unsigned long dimension;
  double complex matrix[8][8];
  unsigned int term_count;
};

static void fail(const char *label)
{
  fprintf(stderr, "FAIL: %s\n", label);
  exit(1);
}

static void expect_close(double complex actual, double complex expected,
                         const char *label)
{
  if (cabs(actual - expected) > 1e-14) fail(label);
}

static int accumulate_matrix(const struct SymmetryTerm *term, void *context)
{
  struct MatrixContext *matrix = context;
  unsigned long state;
  matrix->term_count++;
  for (state = 0; state < matrix->dimension; ++state) {
    unsigned long out = 99UL;
    double complex value = 0.0;
    int status = ApplySymmetryTerm(matrix->def, term, state, &out, &value);
    if (status < 0) return -1;
    if (status == 1) {
      if (out >= matrix->dimension) return -1;
      matrix->matrix[out][state] += value;
    }
  }
  return 0;
}

static int count_terms(const struct SymmetryTerm *term, void *context)
{
  unsigned int *count = context;
  (void)term;
  (*count)++;
  return 0;
}

static void expect_pair_matrix(const struct MatrixContext *matrix,
                               const int (*pairs)[2], const double *couplings,
                               unsigned int pair_count, const char *label)
{
  double complex expected[8][8] = {{0.0}};
  unsigned int p;
  unsigned long state;
  for (p = 0; p < pair_count; ++p) {
    unsigned long mask_a = 1UL << (unsigned int)pairs[p][0];
    unsigned long mask_b = 1UL << (unsigned int)pairs[p][1];
    if (pairs[p][0] == pairs[p][1]) continue;
    for (state = 0; state < matrix->dimension; ++state) {
      int bit_a = (state & mask_a) != 0UL;
      int bit_b = (state & mask_b) != 0UL;
      if (bit_a == bit_b) {
        unsigned long out = state ^ mask_a ^ mask_b;
        expected[out][state] += couplings[p];
      }
    }
  }
  for (state = 0; state < matrix->dimension; ++state) {
    unsigned long out;
    for (out = 0; out < matrix->dimension; ++out)
      if (cabs(matrix->matrix[out][state] - expected[out][state]) > 1e-14)
        fail(label);
  }
}

static void enumerate_pair_matrix(unsigned int nsite, int (*rows)[2],
                                  double *couplings, unsigned int row_count,
                                  struct MatrixContext *matrix)
{
  struct DefineList def;
  int *row_ptrs[4];
  unsigned int p;
  memset(&def, 0, sizeof(def));
  memset(matrix, 0, sizeof(*matrix));
  for (p = 0; p < row_count; ++p) row_ptrs[p] = rows[p];
  def.iCalcModel = SpinGC;
  def.Nsite = nsite;
  def.NPairLiftCoupling = row_count;
  def.PairLiftCoupling = row_ptrs;
  def.ParaPairLiftCoupling = couplings;
  matrix->def = &def;
  matrix->dimension = 1UL << nsite;
  if (EnumerateSymmetryTerms(&def, -1, accumulate_matrix, matrix) != 0)
    fail("enumerate PairLift");
  expect_pair_matrix(matrix, (const int (*)[2])rows, couplings, row_count,
                     "PairLift matrix");
}

static void test_local_algebra(void)
{
  struct DefineList def;
  struct SymmetryTerm plus_plus = {2, {0,1,0,0, 1,1,1,0}, 0.31};
  struct SymmetryTerm product = {2, {0,1,0,0, 0,0,0,1}, 1.0};
  struct SymmetryTerm offsite = {1, {0,1,1,0}, 1.0};
  struct SymmetryTerm bad_local = {1, {0,2,0,0}, 1.0};
  unsigned long out = 99;
  double complex value = 0;
  int rc;
  memset(&def, 0, sizeof(def));
  def.iCalcModel = SpinGC;
  def.Nsite = 2;

  rc = ApplySymmetryTerm(&def, &plus_plus, 0UL, &out, &value);
  if (!(rc == 1 && out == 3UL && cabs(value - 0.31) < 1e-14))
    fail("SpinGC E10 E10");
  plus_plus.index[4] = plus_plus.index[6] = 0;
  if (ApplySymmetryTerm(&def, &plus_plus, 0UL, &out, &value) != 0)
    fail("same-site E10 E10 is zero");

  if (ApplySymmetryTerm(&def, &product, 1UL, &out, &value) != 1 || out != 1UL)
    fail("same-site E10 E01");
  product.index[1] = 0;
  product.index[3] = 1;
  product.index[5] = 1;
  product.index[7] = 0;
  if (ApplySymmetryTerm(&def, &product, 0UL, &out, &value) != 1 || out != 0UL)
    fail("same-site E01 E10");
  if (ApplySymmetryTerm(&def, &offsite, 0UL, &out, &value) != -1)
    fail("off-site SpinGC Transfer rejected");
  if (ApplySymmetryTerm(&def, &bad_local, 0UL, &out, &value) != -1)
    fail("invalid SpinGC local index rejected");
}

static void test_pair_lift(void)
{
  struct MatrixContext matrix;
  int one[1][2] = {{0, 1}};
  int reverse[2][2] = {{0, 1}, {1, 0}};
  int duplicate[2][2] = {{0, 1}, {0, 1}};
  int onsite[1][2] = {{0, 0}};
  int triangle[3][2] = {{0, 1}, {1, 2}, {0, 2}};
  double j_one[1] = {0.31};
  double j_two[2] = {0.31, 0.31};
  double j_zero[1] = {0.31};
  double j_triangle[3] = {0.11, 0.17, 0.23};

  enumerate_pair_matrix(2, one, j_one, 1, &matrix);
  if (matrix.term_count != 2U) fail("one PairLift emits two terms");
  expect_close(matrix.matrix[3][0], 0.31, "PairLift 00 to 11");
  expect_close(matrix.matrix[0][3], 0.31, "PairLift 11 to 00");

  enumerate_pair_matrix(2, reverse, j_two, 2, &matrix);
  expect_close(matrix.matrix[3][0], 0.62, "reverse PairLift accumulates");
  expect_close(matrix.matrix[0][3], 0.62, "reverse PairLift adjoint accumulates");

  enumerate_pair_matrix(2, duplicate, j_two, 2, &matrix);
  expect_close(matrix.matrix[3][0], 0.62, "duplicate PairLift accumulates");
  expect_close(matrix.matrix[0][3], 0.62, "duplicate PairLift adjoint accumulates");

  enumerate_pair_matrix(2, onsite, j_zero, 1, &matrix);
  enumerate_pair_matrix(3, triangle, j_triangle, 3, &matrix);
}

static void test_iterator_families(void)
{
  struct DefineList def;
  struct MatrixContext matrix;
  unsigned int count = 99U;
  int exchange_row[1][2] = {{0, 1}};
  int *exchange_rows[1] = {exchange_row[0]};
  double exchange_values[1] = {0.5};
  int coulomb_row[1][2] = {{0, 1}};
  int *coulomb_rows[1] = {coulomb_row[0]};
  double coulomb_values[1] = {-0.25};
  int hund_row[1][2] = {{0, 1}};
  int *hund_rows[1] = {hund_row[0]};
  double hund_values[1] = {-0.5};
  double diagonal;

  memset(&def, 0, sizeof(def));
  def.iCalcModel = SpinGC;
  def.Nsite = 2;
  count = 0;
  if (!SymmetryUsesExtendedTerms(&def)) fail("empty SpinGC uses term iterator");
  if (EnumerateSymmetryTerms(&def, -1, count_terms, &count) != 0 || count != 0U)
    fail("empty SpinGC Hamiltonian");
  if (EvaluateSymmetryStateDiagonal(&def, 0UL, &diagonal) != 0 || diagonal != 0.0)
    fail("SpinGC diagonal model gate");

  def.NExchangeCoupling = 1;
  def.ExchangeCoupling = exchange_rows;
  def.ParaExchangeCoupling = exchange_values;
  memset(&matrix, 0, sizeof(matrix));
  matrix.def = &def;
  matrix.dimension = 4UL;
  if (EnumerateSymmetryTerms(&def, -1, accumulate_matrix, &matrix) != 0 ||
      matrix.term_count != 2U)
    fail("pure Exchange SpinGC iterator");
  expect_close(matrix.matrix[2][1], 0.5, "SpinGC Exchange 01 to 10");
  expect_close(matrix.matrix[1][2], 0.5, "SpinGC Exchange 10 to 01");

  def.NExchangeCoupling = 0;
  def.NIsingCoupling = 1;
  def.NCoulombInter = 1;
  def.CoulombInter = coulomb_rows;
  def.ParaCoulombInter = coulomb_values;
  def.NHundCoupling = 1;
  def.HundCoupling = hund_rows;
  def.ParaHundCoupling = hund_values;
  count = 0;
  if (EnumerateSymmetryTerms(&def, -1, count_terms, &count) != 0 || count != 6U)
    fail("pure Ising reader expansion only");
  if (EvaluateSymmetryStateDiagonal(&def, 0UL, &diagonal) != 0)
    fail("SpinGC Ising diagonal");
  expect_close(diagonal, 0.25, "SpinGC Ising parallel diagonal");
}

static void test_invalid_enumeration(void)
{
  struct DefineList def;
  int transfer[1][4] = {{0, 1, 0, 0}};
  int *transfer_rows[1] = {transfer[0]};
  double complex transfer_values[1] = {NAN};
  unsigned int count = 0;
  memset(&def, 0, sizeof(def));
  def.iCalcModel = SpinGC;
  def.Nsite = 2;
  def.EDNTransfer = 1;
  def.EDGeneralTransfer = transfer_rows;
  def.EDParaGeneralTransfer = transfer_values;
  if (EnumerateSymmetryTerms(&def, -1, count_terms, &count) != -1)
    fail("non-finite SpinGC coefficient rejected");
  transfer_values[0] = 1.0;
  transfer[0][1] = 2;
  if (EnumerateSymmetryTerms(&def, -1, count_terms, &count) != -1)
    fail("invalid enumerated local index rejected");
  transfer[0][1] = 1;
  transfer[0][2] = 1;
  if (EnumerateSymmetryTerms(&def, -1, count_terms, &count) != -1)
    fail("off-site enumerated SpinGC Transfer rejected");
}

static void test_validation(void)
{
  struct DefineList def;
  int transfer[2][4] = {{0, 1, 0, 0}, {1, 1, 1, 0}};
  int *transfer_rows[2] = {transfer[0], transfer[1]};
  double complex transfer_values[2] = {0.4, 0.4};
  int permutation[2] = {1, 0};
  int *permutations[1] = {permutation};
  memset(&def, 0, sizeof(def));
  def.iCalcModel = Spin;
  def.Nsite = 2;
  def.EDNTransfer = 1;
  def.EDGeneralTransfer = transfer_rows;
  def.EDParaGeneralTransfer = transfer_values;
  if (ValidateSymmetryTerms(&def) != -1)
    fail("canonical Spin transverse field fixed-Sz rejection");

  def.iCalcModel = SpinGC;
  def.EDNTransfer = 2;
  def.NSymTrans = 1;
  def.SymTrans = permutations;
  if (ValidateSymmetryTerms(&def) != 0)
    fail("SpinGC transverse field and permutation invariance");
  def.EDNTransfer = 1;
  if (ValidateSymmetryTerms(&def) != -1)
    fail("SpinGC still checks site permutation invariance");
}

/* Px maps E_st to E_(1-s)(1-t) without conjugating the value: Sx and
 * Sy_0 Sz_1 + Sz_0 Sy_1 (complex units) are invariant, Sy and Sz are not. */
static void test_spin_flip_validation(void)
{
  struct DefineList def;
  int identity[2] = {0, 1};
  int *permutations[2] = {identity, identity};
  int flips[2] = {0, 1};
  int transfer[2][4] = {{0, 1, 0, 0}, {0, 0, 0, 1}};
  int *transfer_rows[2] = {transfer[0], transfer[1]};
  double complex sx[2] = {-0.4, -0.4}, sy[2] = {0.5 * I, -0.5 * I};
  int diagonal[2][4] = {{0, 1, 0, 1}, {0, 0, 0, 0}};
  int *diagonal_rows[2] = {diagonal[0], diagonal[1]};
  double complex sz[2] = {-0.5, 0.5};
  int interall[8][8];
  int *interall_rows[8];
  double complex interall_values[8];
  unsigned int row = 0, c, site;
  for (site = 0; site < 2; ++site) {
    for (c = 0; c < 2; ++c) {
      double sz_c = c ? 0.5 : -0.5;
      /* Sy_site Sz_other: E10 carries -i/2, E01 carries +i/2. */
      int pattern[2][8] = {{(int)site, 1, (int)site, 0, (int)(1-site), (int)c, (int)(1-site), (int)c},
                           {(int)site, 0, (int)site, 1, (int)(1-site), (int)c, (int)(1-site), (int)c}};
      memcpy(interall[row], pattern[0], sizeof(pattern[0]));
      interall_values[row++] = -0.5 * I * sz_c;
      memcpy(interall[row], pattern[1], sizeof(pattern[1]));
      interall_values[row++] = 0.5 * I * sz_c;
    }
  }
  for (row = 0; row < 8; ++row) interall_rows[row] = interall[row];
  memset(&def, 0, sizeof(def));
  def.iFlgSymmetryBasis = TRUE;
  def.iCalcModel = SpinGC;
  def.Nsite = 1;
  def.iSpinFlipParity = -1;
  def.NSymTrans = 2;
  def.NSymTransInput = 1;
  def.SymTrans = permutations;
  def.SymTransFlip = flips;
  def.EDNTransfer = 2;
  def.EDGeneralTransfer = transfer_rows;
  def.EDParaGeneralTransfer = sx;
  if (ValidateSymmetryTerms(&def) != 0) fail("Px keeps the transverse field");
  def.EDParaGeneralTransfer = sy;
  if (ValidateSymmetryTerms(&def) != -1) fail("Px rejects a Sy field");
  def.EDGeneralTransfer = diagonal_rows;
  def.EDParaGeneralTransfer = sz;
  if (ValidateSymmetryTerms(&def) != -1) fail("Px rejects a Sz field");
  def.Nsite = 2;
  def.EDNTransfer = 0;
  def.NInterAll_OffDiagonal = 8;
  def.InterAll_OffDiagonal = interall_rows;
  def.ParaInterAll_OffDiagonal = interall_values;
  if (ValidateSymmetryTerms(&def) != 0) fail("Px keeps Sy_0 Sz_1 + Sz_0 Sy_1");
  /* Sy_0 (E00 + E11)_1 / 2: the same complex units, now Px odd. */
  def.NInterAll_OffDiagonal = 4;
  interall_values[0] = -interall_values[0];
  interall_values[1] = -interall_values[1];
  if (ValidateSymmetryTerms(&def) != -1) fail("Px rejects Sy_0 times the identity");
  interall_values[0] = -interall_values[0];
  interall_values[1] = -interall_values[1];
  def.NInterAll_OffDiagonal = 8;
  flips[1] = 0;
  if (ValidateSymmetryTerms(&def) != 0) fail("the flag selects Px per operation");
}

int main(void)
{
  stdoutMPI = stdout;
  test_local_algebra();
  test_pair_lift();
  test_iterator_families();
  test_invalid_enumeration();
  test_validation();
  test_spin_flip_validation();
  puts("unittest_symmetry_spingc_terms: PASS");
  return 0;
}
