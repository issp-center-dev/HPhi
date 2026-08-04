#include <limits.h>

#include "DefCommon.h"
#include "symmetry_diagonal.h"
#include "struct.h"

static int state_fits_width(unsigned long int state, unsigned int bit_count)
{
  const unsigned int word_bits =
      (unsigned int)(CHAR_BIT * sizeof(unsigned long int));
  if (bit_count == 0U || bit_count > word_bits) return FALSE;
  if (bit_count == word_bits) return TRUE;
  return (state >> bit_count) == 0UL ? TRUE : FALSE;
}

static int valid_site_pair_storage(unsigned int count,
                                   int **sites,
                                   const double *parameters,
                                   unsigned int nsite)
{
  unsigned int index;
  if (count == 0U) return TRUE;
  if (sites == NULL || parameters == NULL) return FALSE;
  for (index = 0U; index < count; index++) {
    if (sites[index] == NULL ||
        sites[index][0] < 0 || sites[index][1] < 0 ||
        (unsigned int)sites[index][0] >= nsite ||
        (unsigned int)sites[index][1] >= nsite) {
      return FALSE;
    }
  }
  return TRUE;
}

static int valid_site_storage(unsigned int count,
                              int **sites,
                              const double *parameters,
                              unsigned int nsite)
{
  unsigned int index;
  if (count == 0U) return TRUE;
  if (sites == NULL || parameters == NULL) return FALSE;
  for (index = 0U; index < count; index++) {
    if (sites[index] == NULL ||
        sites[index][0] < 0 ||
        (unsigned int)sites[index][0] >= nsite) {
      return FALSE;
    }
  }
  return TRUE;
}

static int has_unsupported_general_diagonal(const struct DefineList *def)
{
  return def->EDNChemi > 0U ||
         def->NInterAll > 0U ||
         def->NInterAll_Diagonal > 0U ||
         def->NNBodyInterAll > 0U ||
         def->NNBodyInterAll_Diagonal > 0U;
}

static int evaluate_spin_diagonal(const struct DefineList *def,
                                  unsigned long int state,
                                  double *value)
{
  unsigned int index;
  if (def->iFlgGeneralSpin != FALSE ||
      has_unsupported_general_diagonal(def) ||
      def->NCoulombIntra > 0U ||
      def->NCoulombInter != def->NIsingCoupling ||
      def->NHundCoupling != def->NIsingCoupling ||
      !valid_site_pair_storage(
          def->NCoulombInter, def->CoulombInter,
          def->ParaCoulombInter, def->Nsite) ||
      !valid_site_pair_storage(
          def->NHundCoupling, def->HundCoupling,
          def->ParaHundCoupling, def->Nsite)) {
    return -1;
  }
  for (index = 0U; index < def->NCoulombInter; index++) {
    *value += def->ParaCoulombInter[index];
  }
  for (index = 0U; index < def->NHundCoupling; index++) {
    unsigned int site0 = (unsigned int)def->HundCoupling[index][0];
    unsigned int site1 = (unsigned int)def->HundCoupling[index][1];
    unsigned long int bit0 = (state >> site0) & 1UL;
    unsigned long int bit1 = (state >> site1) & 1UL;
    if (bit0 == bit1) *value += -def->ParaHundCoupling[index];
  }
  return 0;
}

static int evaluate_spinless_diagonal(const struct DefineList *def,
                                      unsigned long int state,
                                      double *value)
{
  unsigned int index;
  if (has_unsupported_general_diagonal(def) ||
      def->NCoulombIntra > 0U ||
      def->NHundCoupling > 0U ||
      def->NIsingCoupling > 0U ||
      !valid_site_pair_storage(
          def->NCoulombInter, def->CoulombInter,
          def->ParaCoulombInter, def->Nsite)) {
    return -1;
  }
  for (index = 0U; index < def->NCoulombInter; index++) {
    unsigned int site0 = (unsigned int)def->CoulombInter[index][0];
    unsigned int site1 = (unsigned int)def->CoulombInter[index][1];
    if (((state >> site0) & 1UL) != 0UL &&
        ((state >> site1) & 1UL) != 0UL) {
      *value += def->ParaCoulombInter[index];
    }
  }
  return 0;
}

static int evaluate_hubbard_diagonal(const struct DefineList *def,
                                     unsigned long int state,
                                     double *value)
{
  unsigned int index;
  if (has_unsupported_general_diagonal(def) ||
      def->NCoulombInter > 0U ||
      def->NHundCoupling > 0U ||
      def->NIsingCoupling > 0U ||
      !valid_site_storage(
          def->NCoulombIntra, def->CoulombIntra,
          def->ParaCoulombIntra, def->Nsite)) {
    return -1;
  }
  for (index = 0U; index < def->NCoulombIntra; index++) {
    unsigned int site = (unsigned int)def->CoulombIntra[index][0];
    if (((state >> (2U * site)) & 1UL) != 0UL &&
        ((state >> (2U * site + 1U)) & 1UL) != 0UL) {
      *value += def->ParaCoulombIntra[index];
    }
  }
  return 0;
}

int EvaluateSymmetryStateDiagonal(
    const struct DefineList *def,
    unsigned long int state,
    double *diagonal)
{
  const unsigned int word_bits =
      (unsigned int)(CHAR_BIT * sizeof(unsigned long int));
  unsigned int bit_count;
  double value = 0.0;
  int status;
  if (def == NULL || diagonal == NULL || def->Nsite == 0U) return -1;
  if (def->iCalcModel == Hubbard) {
    if (def->Nsite > word_bits / 2U) return -1;
    bit_count = 2U * def->Nsite;
  } else if (def->iCalcModel == Spin ||
             def->iCalcModel == SpinlessFermion) {
    if (def->Nsite > word_bits) return -1;
    bit_count = def->Nsite;
  } else {
    return -1;
  }
  if (!state_fits_width(state, bit_count)) return -1;

  if (def->iCalcModel == Spin) {
    status = evaluate_spin_diagonal(def, state, &value);
  } else if (def->iCalcModel == SpinlessFermion) {
    status = evaluate_spinless_diagonal(def, state, &value);
  } else {
    status = evaluate_hubbard_diagonal(def, state, &value);
  }
  if (status != 0) return -1;
  *diagonal = value;
  return 0;
}
