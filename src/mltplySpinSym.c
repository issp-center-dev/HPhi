#include "mltplySpinSym.h"
#ifdef MPI
#include <mpi.h>
#endif
#include <limits.h>
#include "DefCommon.h"
#include "bitcalc.h"
#include "global.h"
#include "symmetry_basis.h"
#include "struct.h"
#include "wrapperMPI.h"

static int apply_exchange_halfspin(unsigned long int state,
                                   int site0,
                                   int site1,
                                   unsigned long int *out_state)
{
  unsigned long int b0 = (state >> (unsigned int)site0) & 1UL;
  unsigned long int b1 = (state >> (unsigned int)site1) & 1UL;
  if (b0 == b1) return FALSE;
  *out_state = state ^ (1UL << (unsigned int)site0) ^ (1UL << (unsigned int)site1);
  return TRUE;
}

static unsigned long int mask_between_sites(unsigned int site0, unsigned int site1)
{
  unsigned int lo = site0 < site1 ? site0 : site1;
  unsigned int hi = site0 < site1 ? site1 : site0;
  if (hi <= lo + 1U) return 0UL;
  return (1UL << hi) - (1UL << (lo + 1U));
}

static int apply_spinless_hopping_hermite(unsigned long int state,
                                          unsigned int site1,
                                          unsigned int site2,
                                          double complex trans,
                                          unsigned long int *out_state,
                                          double complex *hval)
{
  unsigned long int mask1 = 1UL << site1;
  unsigned long int mask2 = 1UL << site2;
  unsigned long int occupied1 = state & mask1;
  unsigned long int occupied2 = state & mask2;
  int sgn = 1;
  if (site1 == site2) return FALSE;
  if ((occupied1 == 0UL && occupied2 == 0UL) ||
      (occupied1 != 0UL && occupied2 != 0UL)) {
    return FALSE;
  }
  SgnBit(state & mask_between_sites(site1, site2), &sgn);
  *out_state = state ^ mask1 ^ mask2;
  if (occupied1 != 0UL && occupied2 == 0UL) {
    *hval = (double)sgn * conj(trans);
  } else {
    *hval = (double)sgn * trans;
  }
  return TRUE;
}

static int apply_hubbard_hopping_hermite(unsigned long int state,
                                         unsigned int site1,
                                         unsigned int spin1,
                                         unsigned int site2,
                                         unsigned int spin2,
                                         double complex trans,
                                         unsigned long int *out_state,
                                         double complex *hval)
{
  const unsigned int max_bits = (unsigned int)(sizeof(unsigned long int) * CHAR_BIT);
  unsigned int orbital1;
  unsigned int orbital2;
  if (spin1 > 1U || spin2 > 1U) return FALSE;
  if (site1 > (max_bits - spin1) / 2U ||
      site2 > (max_bits - spin2) / 2U) {
    return FALSE;
  }
  orbital1 = 2U * site1 + spin1;
  orbital2 = 2U * site2 + spin2;
  if (orbital1 >= max_bits || orbital2 >= max_bits) return FALSE;
  return apply_spinless_hopping_hermite(state, orbital1, orbital2, trans,
                                        out_state, hval);
}

static int add_canonicalized_transition(struct BindStruct *X,
                                        double complex *tmp_v0,
                                        const double complex *full_v1,
                                        double complex *prdct,
                                        unsigned long int beta,
                                        unsigned long int to_state,
                                        double complex hval,
                                        double complex input_amp)
{
  unsigned long int alpha;
  unsigned long int local_alpha;
  double norm_factor;
  double complex contribution;
  struct SymmetryCanonicalResult result;
  if (SymmetryCanonicalizeState(X, to_state, &result) != 0) return -1;
  if (result.found != TRUE) return 0;
  alpha = result.basis_index;
  if (SymmetryBasisGlobalToLocal(X->Sym, alpha, &local_alpha) != TRUE) return 0;
  norm_factor = X->Sym->basis[alpha].norm / X->Sym->basis[beta].norm;
  contribution = hval * result.phase * norm_factor * input_amp;
  tmp_v0[local_alpha] += contribution;
  *prdct += conj(full_v1[alpha]) * contribution;
  return 0;
}

static const double complex *get_full_input_vector(struct BindStruct *X,
                                                   double complex *tmp_v1)
{
#ifdef MPI
  if (nproc > 1) {
    int ierr;
    if (X->Sym->mpi_full_v1 == NULL || X->Sym->mpi_recvcounts == NULL ||
        X->Sym->mpi_displs == NULL) {
      return NULL;
    }
    ierr = MPI_Allgatherv(&tmp_v1[1], (int)X->Sym->local_dim,
                          MPI_DOUBLE_COMPLEX,
                          &X->Sym->mpi_full_v1[1], X->Sym->mpi_recvcounts,
                          X->Sym->mpi_displs, MPI_DOUBLE_COMPLEX,
                          MPI_COMM_WORLD);
    if (ierr != 0) return NULL;
    return X->Sym->mpi_full_v1;
  }
#else
  (void)X;
#endif
  return tmp_v1;
}

int mltplySpinSym(struct BindStruct *X,
                  double complex *tmp_v0,
                  double complex *tmp_v1)
{
  unsigned long int beta;
  double complex prdct = 0.0;
  const double complex *full_v1 = get_full_input_vector(X, tmp_v1);
  if (full_v1 == NULL) return -1;

  for (beta = 1; beta <= X->Sym->dim; beta++) {
    unsigned int p;
    double complex vin = full_v1[beta];
    if (cabs(vin) == 0.0) continue;

    {
      unsigned long int local_beta;
      if (SymmetryBasisGlobalToLocal(X->Sym, beta, &local_beta) == TRUE) {
        tmp_v0[local_beta] += X->Sym->sym_diagonal[beta] * vin;
        prdct += X->Sym->sym_diagonal[beta] * conj(vin) * vin;
      }
    }

    if (X->Def.iCalcModel == Spin) {
      for (p = 0; p < X->Def.NExchangeCoupling; p++) {
        unsigned long int out_state;
        if (apply_exchange_halfspin(X->Sym->basis[beta].rep_state,
                                    X->Def.ExchangeCoupling[p][0],
                                    X->Def.ExchangeCoupling[p][1],
                                    &out_state) == TRUE) {
          if (add_canonicalized_transition(X, tmp_v0, full_v1, &prdct, beta, out_state,
                                           X->Def.ParaExchangeCoupling[p], vin) != 0) {
            return -1;
          }
        }
      }
    } else if (X->Def.iCalcModel == SpinlessFermion) {
      for (p = 0; p < X->Def.EDNTransfer; p += 2U) {
        unsigned long int out_state;
        double complex hval;
        double complex trans = -X->Def.EDParaGeneralTransfer[p];
        unsigned int site1 = (unsigned int)X->Def.EDGeneralTransfer[p][0];
        unsigned int site2 = (unsigned int)X->Def.EDGeneralTransfer[p][2];
        if (apply_spinless_hopping_hermite(X->Sym->basis[beta].rep_state,
                                           site1, site2, trans,
                                           &out_state, &hval) == TRUE) {
          if (add_canonicalized_transition(X, tmp_v0, full_v1, &prdct, beta, out_state,
                                           hval, vin) != 0) {
            return -1;
          }
        }
      }
    } else if (X->Def.iCalcModel == Hubbard) {
      for (p = 0; p < X->Def.EDNTransfer; p += 2U) {
        unsigned long int out_state;
        double complex hval;
        double complex trans = -X->Def.EDParaGeneralTransfer[p];
        unsigned int site1 = (unsigned int)X->Def.EDGeneralTransfer[p][0];
        unsigned int spin1 = (unsigned int)X->Def.EDGeneralTransfer[p][1];
        unsigned int site2 = (unsigned int)X->Def.EDGeneralTransfer[p][2];
        unsigned int spin2 = (unsigned int)X->Def.EDGeneralTransfer[p][3];
        if (apply_hubbard_hopping_hermite(X->Sym->basis[beta].rep_state,
                                          site1, spin1, site2, spin2, trans,
                                          &out_state, &hval) == TRUE) {
          if (add_canonicalized_transition(X, tmp_v0, full_v1, &prdct, beta, out_state,
                                           hval, vin) != 0) {
            return -1;
          }
        }
      }
    }
  }
  X->Large.prdct += prdct;
  return 0;
}
