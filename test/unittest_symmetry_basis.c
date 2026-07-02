#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

static int perm_storage[4][4];
static int anti_storage[4][4];
static int *perm_rows[4];
static int *anti_rows[4];
static double complex chars4[4];

static void setup_c4_def(struct DefineList *def)
{
  int g, s;
  memset(def, 0, sizeof(*def));
  for (g = 0; g < 4; g++) {
    perm_rows[g] = perm_storage[g];
    anti_rows[g] = anti_storage[g];
  }
  def->iFlgSymmetryBasis = TRUE;
  def->Nsite = 4;
  def->NSymTrans = 4;
  def->SymTrans = perm_rows;
  def->SymTransAnti = anti_rows;
  def->SymTransChar = chars4;
  for (g = 0; g < 4; g++) {
    chars4[g] = 1.0 + 0.0 * I;
    for (s = 0; s < 4; s++) {
      perm_storage[g][s] = (s + g) % 4;
      anti_storage[g][s] = 1;
    }
  }
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
  return 0;
}
