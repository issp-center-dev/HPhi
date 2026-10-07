#include <math.h>
#include "struct.h"
#include "wrapperMPI.h"
#include "symmetry_basis.h"
#include "symmetry_observables.h"

int EvaluateSymmetrySpinGCMoments(struct BindStruct *X,
                                  const double complex *vec)
{
  double local_sz = 0.0;
  double local_sz2 = 0.0;
  int invalid = X == NULL || vec == NULL;

  if (!invalid) {
    invalid = X->Def.iCalcModel != SpinGC ||
              X->Def.iFlgGeneralSpin != FALSE ||
              X->Def.iFlgSymmetryBasis != TRUE ||
              X->Sym == NULL || X->Sym->enabled != TRUE ||
              X->Sym->nsite != (unsigned int)X->Def.Nsite ||
              X->Check.idim_max != X->Sym->local_dim ||
              SymmetryBasisOwnedStorageReady(
                  X->Sym, X->Sym->local_dim) != TRUE;
  }
  if (SumMPI_i(invalid) != 0) return -1;

  invalid = 0;
  for (unsigned long int j = 1; j <= X->Sym->local_dim; ++j) {
    const struct SymmetryBasisVector *entry =
        SymmetryBasisLocalEntry(X->Sym, j);
    double real = creal(vec[j]);
    double imag = cimag(vec[j]);
    if (entry == NULL || !isfinite(real) || !isfinite(imag)) {
      invalid = 1;
      continue;
    }
    unsigned long int bits = entry->rep_state;
    unsigned int count = 0;
    while (bits != 0UL) {
      ++count;
      bits &= bits - 1UL;
    }
    double magnetization = (double)count - 0.5*(double)X->Def.Nsite;
    double weight = real*real + imag*imag;
    local_sz += weight*magnetization;
    local_sz2 += weight*magnetization*magnetization;
  }
  invalid |= !isfinite(local_sz) || !isfinite(local_sz2);
  if (SumMPI_i(invalid) != 0) return -1;

  double sz = SumMPI_d(local_sz);
  double sz2 = SumMPI_d(local_sz2);
  if (SumMPI_i(!isfinite(sz) || !isfinite(sz2)) != 0) return -1;

  X->Phys.Sz = sz;
  X->Phys.Sz2 = sz2;
  X->Phys.num = (double)X->Def.Nsite;
  X->Phys.num2 = (double)X->Def.Nsite*(double)X->Def.Nsite;
  X->Phys.doublon = 0.0;
  X->Phys.doublon2 = 0.0;
  X->Phys.num_up = 0.5*(double)X->Def.Nsite + sz;
  X->Phys.num_down = 0.5*(double)X->Def.Nsite - sz;
  return 0;
}
