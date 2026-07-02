#include <stdio.h>
#include <stdlib.h>
#include "symmetry_basis.h"

FILE *stdoutMPI = NULL;
long unsigned int *list_1 = NULL;
long unsigned int *list_2_1 = NULL;
long unsigned int *list_2_2 = NULL;
double *list_Diagonal = NULL;
int g_tj_odd_split_guard_enabled = 0;
long unsigned int g_tj_odd_split_up_mask = 0;
long unsigned int g_tj_odd_split_down_mask = 0;

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
  return 0;
}
