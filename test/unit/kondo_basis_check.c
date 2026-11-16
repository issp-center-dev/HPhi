/* Compare the production canonical Kondo maps with exhaustive local states. */
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "sz.h"
#include "bitcalc.h"

static unsigned popcount(unsigned long word)
{
  unsigned count = 0;
  for (; word; word >>= 1) count += word & 1UL;
  return count;
}

static void setup(struct BindStruct *x, int sites, unsigned local,
                  int up, int down)
{
  memset(x, 0, sizeof(*x));
  x->Def.iCalcModel = Kondo;
  x->Def.Nsite = x->Def.NsiteMPI = sites;
  x->Def.Nup = up;
  x->Def.Ndown = down;
  x->Def.Ne = up + down;
  x->Def.NLocSpn = popcount(local);
  x->Check.sdim = 1UL << sites;
  x->Def.LocSpn = calloc(sites, sizeof(*x->Def.LocSpn));
  x->Def.Tpow = calloc(2 * sites + 1, sizeof(*x->Def.Tpow));
  if (!x->Def.LocSpn || !x->Def.Tpow) exit(2);
  for (int i = 0; i < sites; ++i) x->Def.LocSpn[i] = (local >> i) & 1U;
  for (int i = 0; i <= 2 * sites; ++i) x->Def.Tpow[i] = 1UL << i;
}

static unsigned long reference(const struct BindStruct *x, unsigned long *words)
{
  unsigned long count = 0;
  for (unsigned long w = 0; w < (1UL << (2 * x->Def.Nsite)); ++w) {
    int up = 0, down = 0, valid = 1;
    for (int i = 0; i < x->Def.Nsite; ++i) {
      unsigned digit = (w >> (2 * i)) & 3UL;
      if (x->Def.LocSpn[i] && digit != 1 && digit != 2) valid = 0;
      up += digit & 1U;
      down += (digit >> 1) & 1U;
    }
    if (valid && up == x->Def.Nup && down == x->Def.Ndown)
      words[++count] = w;
  }
  return count;
}

static void release(struct BindStruct *x)
{
  free(x->Def.LocSpn);
  free(x->Def.Tpow);
}

/* A valid lower-half vacuum must be visited without calling snoob(0). */
static int check_zero(void)
{
  struct BindStruct x;
  unsigned long words[14] = {0}, right[16] = {0}, left[16] = {0}, offsets[16];
  setup(&x, 4, 1U << 3, 2, 1);
  calculate_jb_Kondo(&x, offsets, 16);
  int count = omp_sz_Kondo_hacker(7, 16, &x, words, right, left, offsets);
  int ok = count == 1 && words[offsets[7] + 1] == 0x70;
  if (!ok) fprintf(stderr, "lower-half vacuum: count=%d (expected 1), state=0x%lx (expected 0x70)\n",
                   count, words[offsets[7] + 1]);
  release(&x);
  return ok;
}

static int check_case(int sites, unsigned local, int up, int down)
{
  struct BindStruct x;
  setup(&x, sites, local, up, down);
  unsigned long size = x.Check.sdim;
  unsigned long *expected = calloc(size * size + 1, sizeof(*expected));
  if (!expected) exit(2);
  unsigned long count = reference(&x, expected);
  if (count == 0) {
    free(expected);
    release(&x);
    return 2; /* Production sz() also skips zero-dimensional sectors. */
  }
  x.Check.idim_max = count;
  unsigned long *offsets = calloc(size, sizeof(*offsets));
  unsigned long *words = calloc(count + 2, sizeof(*words));
  unsigned long *right = calloc(size, sizeof(*right));
  unsigned long *left = calloc(size, sizeof(*left));
  if (!offsets || !words || !right || !left) exit(2);
  calculate_jb_Kondo(&x, offsets, size);
  unsigned long irght, ilft, ihfbit;
  if (GetSplitBitByModel(sites, Kondo, &irght, &ilft, &ihfbit) != 0) exit(2);
  int ok = 1;
  for (int fast = 0; fast <= 1 && ok; ++fast) {
    memset(words, 0xff, (count + 2) * sizeof(*words));
    memset(right, 0, size * sizeof(*right));
    memset(left, 0, size * sizeof(*left));
    unsigned long generated = 0;
#pragma omp parallel for reduction(+:generated)
    for (unsigned long ib = 0; ib < size; ++ib) {
      generated += fast ? omp_sz_Kondo_hacker(ib, size, &x, words, right, left, offsets)
                        : omp_sz_Kondo(ib, size, &x, words, right, left, offsets);
    }
    if (generated != count || words[0] != ULONG_MAX || words[count + 1] != ULONG_MAX) ok = 0;
    for (unsigned long i = 1; i <= count && ok; ++i) {
      unsigned long index = 0;
      if (words[i] != expected[i] ||
          !GetOffComp(right, left, expected[i], irght, ilft, ihfbit, &index) || index != i)
        ok = 0;
    }
    if (!ok) fprintf(stderr, "basis mismatch: sites=%d local=0x%x up=%d down=%d CalcHS=%d count=%lu expected=%lu\n",
                     sites, local, up, down, fast, generated, count);
  }
  free(expected); free(offsets); free(words); free(right); free(left);
  release(&x);
  return ok;
}

int main(int argc, char **argv)
{
  stdoutMPI = tmpfile(); /* Suppress the production prefix-count diagnostics. */
  if (!stdoutMPI) return 2;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  if (argc == 2 && strcmp(argv[1], "--underflow") == 0) {
    int ok = check_case(4, 1, 2, 1);
    fclose(stdoutMPI);
    return ok ? 0 : 1;
  }
  if (!check_zero()) return 1;
  if (argc == 2 && strcmp(argv[1], "--zero") == 0) {
    fclose(stdoutMPI);
    return 0;
  }
  if (argc != 1) return 2;
  unsigned cases = 0, empty = 0;
  /* All mixed local/itinerant layouts and quantum numbers, including empty
     blocks, fully polarized spins, odd central sites and full occupation. */
  for (int sites = 2; sites <= 5; ++sites)
    for (unsigned local = 1; local + 1 < (1U << sites); ++local)
      for (int up = 0; up <= sites; ++up)
        for (int down = 0; down <= sites; ++down) {
          int result = check_case(sites, local, up, down);
          if (!result) return 1;
          if (result == 2) ++empty;
          else ++cases;
        }
  fclose(stdoutMPI);
  printf("Kondo basis: %u nonempty sectors passed (CalcHS 0/1, %d threads); %u empty sectors skipped\n",
         cases, nthreads, empty);
  return 0;
}
