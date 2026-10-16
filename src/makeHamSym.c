#include <limits.h>
#include <math.h>
#include <stdint.h>
#include "Common.h"
#include "hamstore.h"
#include "makeHamSym.h"
#include "symmetry_basis.h"
#include "symmetry_matvec_plan.h"

struct DenseSymmetryColumn {
  unsigned long dimension;
  long column;
};

static int store_entry(unsigned long row, double complex value, void *opaque)
{
  const struct DenseSymmetryColumn *column = opaque;
  if (row == 0 || row > column->dimension ||
      !isfinite(creal(value)) || !isfinite(cimag(value))) return -1;
  AddHamElem(row, column->column, value);
  return 0;
}

int makeHamSym(const struct BindStruct *X)
{
  unsigned long dimension, row;
  long first, last, column;
  struct DenseSymmetryColumn sink;
  if (X == NULL || X->Sym == NULL || X->Sym->enabled != TRUE ||
      X->Sym->basis_layout != SYMMETRY_BASIS_REPLICATED ||
      X->Sym->basis == NULL || X->Sym->dim == 0 || X->Sym->dim >= LONG_MAX ||
      iHamSinkMode == HAM_SINK_TRACE_COLLECT) return -1;
  dimension = X->Sym->dim;
  first = 1;
  last = (long)dimension;
  if (iHamPanelActive) {
    size_t elements;
    if (Ham_local == NULL || HamPanelLd != (long)dimension) return -1;
    first = HamColBegin;
    last = HamColEnd;
    /* Empty panels use the same canonical range as setmem_large(). */
    if (first == 1 && last == 0) return 0;
    if (first < 1 || last < first || (unsigned long)last > dimension ||
        (size_t)(last-first+1) > SIZE_MAX / dimension / sizeof(*Ham_local)) return -1;
    elements = (size_t)dimension * (size_t)(last-first+1);
    memset(Ham_local, 0, elements * sizeof(*Ham_local));
  } else {
    if (Ham == NULL || dimension > SIZE_MAX / sizeof(**Ham) - 1) return -1;
    for (row = 0; row <= dimension; ++row)
      if (Ham[row] == NULL) return -1;
    for (row = 0; row <= dimension; ++row)
      memset(Ham[row], 0, (dimension+1) * sizeof(**Ham));
  }
  sink.dimension = dimension;
  for (column = first; column <= last; ++column) {
    sink.column = column;
    if (SymmetryEnumerateColumn(X, (unsigned long)column, store_entry, &sink) != 0)
      return -1;
  }
  return 0;
}
