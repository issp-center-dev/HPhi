#include <limits.h>
#include <stdint.h>
#include <stdio.h>

#include "DefCommon.h"
#include "global.h"
struct BindStruct;
#include "struct.h"
#include "symmetry_kondo.h"

static FILE *kondo_diagnostic_stream(void)
{
  return stdoutMPI != NULL ? stdoutMPI : stderr;
}

static int kondo_error(const char *message)
{
  fprintf(kondo_diagnostic_stream(), "Error: %s\n", message);
  return -1;
}

int IsSymmetryKondoModel(int model)
{
  return model == Kondo || model == KondoNConserved || model == KondoGC;
}

static int canonical_sector_exists(int64_t local_count,
                                   int64_t conduction_count,
                                   int64_t nup, int64_t ndown)
{
  int64_t lower = 0;
  int64_t upper = local_count;
  int64_t candidate;

  candidate = nup - conduction_count;
  if (candidate > lower) lower = candidate;
  candidate = local_count - ndown;
  if (candidate > lower) lower = candidate;

  if (nup < upper) upper = nup;
  candidate = local_count + conduction_count - ndown;
  if (candidate < upper) upper = candidate;
  return lower <= upper;
}

int NormalizeSymmetryKondoQuantumNumbers(struct DefineList *def,
                                         int has_ncond, int has_sz,
                                         int has_nup, int has_ndown)
{
  struct DefineList normalized;
  int64_t nsite;
  int64_t nlocal;
  int64_t nconduction;
  int64_t ncond;
  int64_t nup = 0;
  int64_t ndown = 0;
  int64_t total2sz = 0;
  int64_t ne;

  if (def == NULL) return kondo_error("missing Kondo definition");
  if (def->iCalcModel != Kondo && def->iCalcModel != KondoGC)
    return kondo_error("Kondo quantum-number normalization used for another model");
  nsite = (int64_t)def->Nsite;
  nlocal = (int64_t)def->NLocSpn;
  if (nsite <= 0) return kondo_error("Kondo TransSym requires Nsite > 0");
  if (nlocal > nsite)
    return kondo_error("NLocalSpin must not exceed Nsite");
  nconduction = nsite - nlocal;
  normalized = *def;

  if (def->iCalcModel == KondoGC) {
    if (has_ncond || has_sz || has_nup || has_ndown)
      return kondo_error(
          "KondoGC TransSym does not accept explicit 2Sz/Nup/Ndown/Ncond");
    normalized.NCond = 0U;
    normalized.Nup = 0U;
    normalized.Ndown = 0U;
    normalized.Ne = 0U;
    normalized.Total2Sz = 0;
    normalized.iFlgSzConserved = FALSE;
    *def = normalized;
    return 0;
  }

  if ((has_nup != FALSE) != (has_ndown != FALSE))
    return kondo_error("Nup and Ndown must be specified together");

  if (has_nup && has_ndown) {
    nup = (int64_t)def->Nup;
    ndown = (int64_t)def->Ndown;
    ncond = nup + ndown - nlocal;
    total2sz = nup - ndown;
    if (has_ncond && (int64_t)def->NCond != ncond)
      return kondo_error("Ncond conflicts with Nup/Ndown");
    if (has_sz && (int64_t)def->Total2Sz != total2sz)
      return kondo_error("2Sz conflicts with Nup/Ndown");
  } else {
    if (!has_ncond)
      return kondo_error("Kondo TransSym requires Ncond or Nup/Ndown");
    ncond = (int64_t)def->NCond;
    if (!has_sz) {
      ne = nlocal + ncond;
      if (ncond < 0 || ncond > 2 * nconduction || ne < nlocal ||
          ne > UINT_MAX)
        return kondo_error("Ncond is outside the Kondo physical space");
      normalized.iCalcModel = KondoNConserved;
      normalized.NCond = (unsigned int)ncond;
      normalized.Nup = 0U;
      normalized.Ndown = 0U;
      normalized.Ne = (unsigned int)ne;
      normalized.Total2Sz = 0;
      normalized.iFlgSzConserved = FALSE;
      *def = normalized;
      return 0;
    }
    total2sz = (int64_t)def->Total2Sz;
    ne = nlocal + ncond;
    if (((ne + total2sz) % 2) != 0 || ((ne - total2sz) % 2) != 0)
      return kondo_error("Ncond and 2Sz have incompatible parity");
    nup = (ne + total2sz) / 2;
    ndown = (ne - total2sz) / 2;
  }

  ne = nlocal + ncond;
  if (ncond < 0 || ncond > 2 * nconduction || ne < nlocal ||
      ne > UINT_MAX)
    return kondo_error("Ncond is outside the Kondo physical space");
  if (nup < 0 || ndown < 0 || nup > nsite || ndown > nsite ||
      nup + ndown != ne || nup - ndown != total2sz ||
      total2sz < INT_MIN || total2sz > INT_MAX)
    return kondo_error("Kondo fixed quantum numbers are inconsistent");
  if (!canonical_sector_exists(nlocal, nconduction, nup, ndown))
    return kondo_error("Kondo fixed sector has no physical states");

  normalized.iCalcModel = Kondo;
  normalized.NCond = (unsigned int)ncond;
  normalized.Nup = (unsigned int)nup;
  normalized.Ndown = (unsigned int)ndown;
  normalized.Ne = (unsigned int)ne;
  normalized.Total2Sz = (int)total2sz;
  normalized.iFlgSzConserved = TRUE;
  *def = normalized;
  return 0;
}

static int kondo_word_width_valid(const struct DefineList *def)
{
  const unsigned int word_bits =
      (unsigned int)(CHAR_BIT * sizeof(unsigned long));
  if (def->Nsite == 0U) return FALSE;
  if (def->iCalcModel == KondoGC)
    return def->Nsite < word_bits / 2U;
  return def->Nsite <= word_bits / 2U;
}

int ValidateSymmetryKondoSpace(const struct DefineList *def)
{
  unsigned int site;
  unsigned int local_count = 0U;
  int64_t nsite;
  int64_t nlocal;
  int64_t nconduction;
  int64_t ncond;
  int64_t ne;
  int64_t nup;
  int64_t ndown;

  if (def == NULL) return kondo_error("missing Kondo definition");
  if (!IsSymmetryKondoModel(def->iCalcModel))
    return kondo_error("Kondo physical-space validation used for another model");
  if (!kondo_word_width_valid(def))
    return kondo_error("Kondo TransSym site count exceeds the state word");
  if (def->LocSpn == NULL)
    return kondo_error("Kondo TransSym requires LocSpn for every site");

  for (site = 0U; site < def->Nsite; ++site) {
    if (def->LocSpn[site] > LOCSPIN || def->iFlgGeneralSpin == TRUE)
      return kondo_error("Kondo TransSym supports only Spin-1/2 local spins");
    if (def->LocSpn[site] < ITINERANT)
      return kondo_error("Kondo TransSym has an invalid LocSpn value");
    if (def->LocSpn[site] == LOCSPIN) local_count++;
  }
  if (local_count != def->NLocSpn)
    return kondo_error("NLocalSpin does not match the LocSpn array");

  nsite = (int64_t)def->Nsite;
  nlocal = (int64_t)def->NLocSpn;
  nconduction = nsite - nlocal;
  ncond = (int64_t)def->NCond;
  ne = (int64_t)def->Ne;
  if (def->iCalcModel == KondoGC) {
    if (def->NCond != 0U || def->Nup != 0U || def->Ndown != 0U ||
        def->Ne != 0U || def->Total2Sz != 0 ||
        def->iFlgSzConserved != FALSE)
      return kondo_error("KondoGC fixed quantities were not normalized");
    return 0;
  }
  if (ncond < 0 || ncond > 2 * nconduction || ne != nlocal + ncond ||
      ne < nlocal)
    return kondo_error("Kondo particle numbers are inconsistent");
  if (def->iCalcModel == KondoNConserved) {
    if (def->Nup != 0U || def->Ndown != 0U || def->Total2Sz != 0 ||
        def->iFlgSzConserved != FALSE)
      return kondo_error("KondoNConserved spin quantities were not normalized");
    return 0;
  }

  nup = (int64_t)def->Nup;
  ndown = (int64_t)def->Ndown;
  if (def->iFlgSzConserved != TRUE || nup > nsite || ndown > nsite ||
      nup + ndown != ne || nup - ndown != (int64_t)def->Total2Sz)
    return kondo_error("Kondo fixed quantum numbers are inconsistent");
  if (!canonical_sector_exists(nlocal, nconduction, nup, ndown))
    return kondo_error("Kondo fixed sector has no physical states");
  return 0;
}

int SymmetryKondoLocalMask(const struct DefineList *def, unsigned long *mask)
{
  unsigned int site;
  unsigned long result = 0UL;
  if (mask == NULL) return kondo_error("missing Kondo local mask output");
  if (ValidateSymmetryKondoSpace(def) != 0) return -1;
  for (site = 0U; site < def->Nsite; ++site) {
    if (def->LocSpn[site] == LOCSPIN) result |= 1UL << site;
  }
  *mask = result;
  return 0;
}

int SymmetryKondoStateIsPhysical(const struct DefineList *def,
                                 unsigned long state)
{
  const unsigned int word_bits =
      (unsigned int)(CHAR_BIT * sizeof(unsigned long));
  unsigned int site;
  unsigned int used_bits;
  if (def == NULL || def->LocSpn == NULL ||
      !IsSymmetryKondoModel(def->iCalcModel) ||
      !kondo_word_width_valid(def))
    return 0;
  used_bits = 2U * def->Nsite;
  if (used_bits < word_bits && (state >> used_bits) != 0UL) return 0;
  for (site = 0U; site < def->Nsite; ++site) {
    unsigned int digit = (unsigned int)((state >> (2U * site)) & 3UL);
    if (def->LocSpn[site] == LOCSPIN && digit != 1U && digit != 2U)
      return 0;
  }
  return 1;
}

int SymmetryKondoPermutationSign(const struct DefineList *def,
                                 const int *permutation, int *sign)
{
  unsigned int source;
  unsigned int other;
  unsigned long seen = 0UL;
  unsigned int inversion_parity = 0U;
  int result;

  if (permutation == NULL || sign == NULL)
    return kondo_error("missing Kondo permutation or sign output");
  if (ValidateSymmetryKondoSpace(def) != 0) return -1;
  for (source = 0U; source < def->Nsite; ++source) {
    int target = permutation[source];
    if (target < 0 || (unsigned int)target >= def->Nsite)
      return kondo_error("Kondo site permutation is out of range");
    if ((seen & (1UL << (unsigned int)target)) != 0UL)
      return kondo_error("Kondo site permutation is not bijective");
    seen |= 1UL << (unsigned int)target;
    if (def->LocSpn[source] != def->LocSpn[target])
      return kondo_error("Kondo site permutation changes site type");
  }
  for (source = 0U; source < def->Nsite; ++source) {
    if (def->LocSpn[source] != LOCSPIN) continue;
    for (other = source + 1U; other < def->Nsite; ++other) {
      if (def->LocSpn[other] == LOCSPIN &&
          permutation[source] > permutation[other])
        inversion_parity ^= 1U;
    }
  }
  result = inversion_parity == 0U ? 1 : -1;
  *sign = result;
  return 0;
}
