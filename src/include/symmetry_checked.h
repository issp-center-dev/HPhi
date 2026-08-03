#ifndef HPHI_SYMMETRY_CHECKED_H
#define HPHI_SYMMETRY_CHECKED_H

#include <limits.h>
#include <stddef.h>
#include <stdint.h>

#include "symmetry_distribution.h"

static inline int SymmetryCheckedU64Add(
    uint64_t lhs,
    uint64_t rhs,
    uint64_t *result)
{
  if (result == NULL || lhs > UINT64_MAX - rhs) return -1;
  *result = lhs + rhs;
  return 0;
}

static inline int SymmetryCheckedU64Mul(
    uint64_t lhs,
    uint64_t rhs,
    uint64_t *result)
{
  if (result == NULL || (lhs != 0U && rhs > UINT64_MAX / lhs)) return -1;
  *result = lhs * rhs;
  return 0;
}

static inline int SymmetryCheckedSizeMul(
    size_t lhs,
    size_t rhs,
    size_t *result)
{
  if (result == NULL || (lhs != 0U && rhs > SIZE_MAX / lhs)) return -1;
  *result = lhs * rhs;
  return 0;
}

static inline int SymmetryCheckedSizeAdd(
    size_t lhs,
    size_t rhs,
    size_t *result)
{
  if (result == NULL || lhs > SIZE_MAX - rhs) return -1;
  *result = lhs + rhs;
  return 0;
}

static inline int SymmetryCheckedU64ToSize(
    uint64_t value,
    size_t *result)
{
  if (result == NULL || value > (uint64_t)SIZE_MAX) return -1;
  *result = (size_t)value;
  return 0;
}

static inline int SymmetryCheckedUlongToU64(
    unsigned long int value,
    uint64_t *result)
{
  uint64_t converted;
  if (result == NULL) return -1;
  converted = (uint64_t)value;
  if ((unsigned long int)converted != value) return -1;
  *result = converted;
  return 0;
}

static inline int SymmetryBasisRunIsValid(
    const struct SymmetryBasisRun *run)
{
  if (run == NULL) return 0;
  if (run->entries == NULL) {
    return run->count == 0UL && run->capacity == 0UL;
  }
  if (run->count == ULONG_MAX || run->capacity == 0UL) return 0;
  return run->capacity >= run->count + 1UL;
}

#endif /* HPHI_SYMMETRY_CHECKED_H */
