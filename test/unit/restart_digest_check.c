#include <complex.h>
#include <stdio.h>
#include "restart_checkpoint.h"
#include "restart_check.h"
int main(void)
{
  enum { N = 37 };
  double complex v[N + 1] = {0}, w[N + 1];
  for (int i = 1; i <= N; ++i) v[i] = (0.1 * i) + I * (1.0 / i);
  uint64_t x1 = 0, s1 = 0, x2 = 0, s2 = 0, x3 = 0, s3 = 0, lo, hi, word;
  RestartDigestAccumulate(0, v, 10, &x1, &s1);
  RestartDigestAccumulate(10, v + 10, N - 10, &x1, &s1);          /* partition A: [0,10) [10,37) */
  RestartDigestAccumulate(0, v, 20, &x2, &s2);
  RestartDigestAccumulate(20, v + 20, N - 20, &x2, &s2);          /* partition B: [0,20) [20,37) */
  CHECK(RW_FLAGS == 39);
  CHECK(x1 == x2 && s1 == s2);
  for (int i = 0; i <= N; ++i) w[i] = v[i];
  w[5] = v[6]; w[6] = v[5];
  RestartDigestAccumulate(0, w, N, &x3, &s3);
  CHECK(x3 != x1);                                                /* reordering is detected */
  CHECK(RestartIntersect(0, 8, 2, 2, &lo, &hi) == 1 && lo == 2 && hi == 4);   /* np4 -> np1 second file */
  CHECK(RestartIntersect(0, 8, 0, 2, &lo, &hi) == 1 && lo == 0 && hi == 2);
  CHECK(RestartIntersect(7, 14, 10, 2, &lo, &hi) == 1 && lo == 10 && hi == 12);
  CHECK(RestartIntersect(0, 4, 6, 2, &lo, &hi) == 0);                        /* disjoint */
  CHECK(RestartIntersect(0, 4, 2, 0, &lo, &hi) == 0);                        /* empty old rank */
  FILE *fp = tmpfile(); CHECK(fp != NULL);
  CHECK(RestartPutWord(fp, UINT64_C(0x0a31565349485048)) == 0); rewind(fp);
  unsigned char b[8]; CHECK(fread(b, 1, 8, fp) == 8 && b[0] == 0x48 && b[7] == 0x0a); rewind(fp);
  CHECK(RestartGetWord(fp, &word) == 0 && word == UINT64_C(0x0a31565349485048)); fclose(fp);
  CHECK(RestartWordDouble(RestartDoubleWord(-0.0)) == 0.0 && RestartDoubleWord(-0.0) != RestartDoubleWord(0.0));
  puts("restart_digest_check: ok");
  return 0;
}
