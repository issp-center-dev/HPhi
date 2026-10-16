#include <limits.h>
#include <math.h>
#include "symmetry_terms.h"
#include "symmetry_kondo.h"

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

struct DiagonalContext {
  const struct DefineList *def;
  unsigned long state;
  double value;
};

static int accumulate_diagonal(const struct SymmetryTerm *term, void *context)
{
  struct DiagonalContext *diagonal = context;
  unsigned long out;
  double complex value;
  int status = ApplySymmetryTerm(diagonal->def, term, diagonal->state, &out, &value);
  if (status < 0) return -1;
  if (status) {
    if (out != diagonal->state || fabs(cimag(value)) > 1e-10) return -1;
    diagonal->value += creal(value);
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
  struct DiagonalContext context = {def, state, 0.0};
  if (def == NULL || diagonal == NULL || def->Nsite == 0U) return -1;
  if (def->iCalcModel == Hubbard || def->iCalcModel == tJ ||
      IsSymmetryKondoModel(def->iCalcModel)) {
    if (def->Nsite > word_bits / 2U) return -1;
    bit_count = 2U * def->Nsite;
  } else if (def->iCalcModel == Spin || def->iCalcModel == SpinGC ||
             def->iCalcModel == SpinlessFermion) {
    if (def->Nsite > word_bits) return -1;
    bit_count = def->Nsite;
  } else {
    return -1;
  }
  if (!state_fits_width(state, bit_count)) return -1;
  if (IsSymmetryKondoModel(def->iCalcModel) &&
      !SymmetryKondoStateIsPhysical(def, state)) return -1;

  if (EnumerateSymmetryTerms(def, 0, accumulate_diagonal, &context)) return -1;
  *diagonal = context.value;
  return 0;
}
