#include <string.h>
#include "restart_checkpoint.h"
#define RESTART_FNV_OFFSET UINT64_C(14695981039346656037)
#define RESTART_FNV_PRIME  UINT64_C(1099511628211)

static uint64_t hash_bytes(uint64_t hash, const unsigned char *bytes, size_t n)
{
  for (size_t i = 0; i < n; ++i) { hash ^= bytes[i]; hash *= RESTART_FNV_PRIME; }
  return hash;
}
uint64_t RestartHashTag(const char *tag) { return hash_bytes(RESTART_FNV_OFFSET, (const unsigned char *)tag, strlen(tag)); }
uint64_t RestartHashWord(uint64_t hash, uint64_t value)
{
  unsigned char bytes[8];
  for (int i = 0; i < 8; ++i) { bytes[i] = value & 255; value >>= 8; }
  return hash_bytes(hash, bytes, 8);
}
uint64_t RestartDoubleWord(double value) { uint64_t w; memcpy(&w, &value, 8); return w; }
double   RestartWordDouble(uint64_t word) { double v; memcpy(&v, &word, 8); return v; }
int RestartPutWord(FILE *fp, uint64_t word)
{
  unsigned char bytes[8];
  for (int i = 0; i < 8; ++i) { bytes[i] = word & 255; word >>= 8; }
  return fwrite(bytes, 1, 8, fp) == 8 ? 0 : -1;
}
int RestartGetWord(FILE *fp, uint64_t *word)
{
  unsigned char bytes[8];
  if (fread(bytes, 1, 8, fp) != 8) return -1;
  *word = 0;
  for (int i = 7; i >= 0; --i) *word = (*word << 8) | bytes[i];
  return 0;
}
void RestartDigestAccumulate(uint64_t offset, const double complex *vector, unsigned long count,
                             uint64_t *xor_hash, uint64_t *sum_hash)
{
  const uint64_t tag = RestartHashTag("hphi-ckpt-data-v3");
  for (unsigned long i = 1; i <= count; ++i) {
    uint64_t h = RestartHashWord(tag, offset + i);
    h = RestartHashWord(h, RestartDoubleWord(creal(vector[i])));
    h = RestartHashWord(h, RestartDoubleWord(cimag(vector[i])));
    *xor_hash ^= h; *sum_hash += h;
  }
}
int RestartIntersect(uint64_t a, uint64_t b, uint64_t offset, uint64_t count, uint64_t *lo, uint64_t *hi)
{
  uint64_t end = offset + count;   /* caller checked offset + count against overflow */
  *lo = a > offset ? a : offset;
  *hi = b < end ? b : end;
  return *hi > *lo;
}
