#include <bitcalc.h>
#include "mltplySpinSym.h"
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

static int raw_index_from_state_sym(const struct BindStruct *X,
                                    unsigned long int state,
                                    unsigned long int *raw_index)
{
  return GetOffComp(list_2_1, list_2_2, state,
                    X->Large.irght, X->Large.ilft, X->Large.ihfbit,
                    raw_index);
}

static void add_raw_transition(struct BindStruct *X,
                               double complex *tmp_v0,
                               double complex *tmp_v1,
                               double complex *prdct,
                               unsigned long int beta,
                               double complex beta_amp,
                               unsigned long int to_state,
                               double complex hval,
                               double complex input_amp)
{
  unsigned long int to_raw = 0;
  unsigned long int alpha;
  double complex alpha_coeff;
  double complex contribution;
  if (raw_index_from_state_sym(X, to_state, &to_raw) != TRUE) return;
  alpha = X->Sym->raw_to_sym[to_raw];
  if (alpha == 0) return;
  alpha_coeff = X->Sym->raw_to_coeff[to_raw];
  contribution = hval * beta_amp * conj(alpha_coeff) * input_amp;
  tmp_v0[alpha] += contribution;
  *prdct += conj(tmp_v1[alpha]) * contribution;
  (void)beta;
}

int mltplySpinSym(struct BindStruct *X,
                  double complex *tmp_v0,
                  double complex *tmp_v1)
{
  unsigned long int beta;
  double complex prdct = 0.0;

  for (beta = 1; beta <= X->Sym->dim; beta++) {
    unsigned int p;
    double complex vin = tmp_v1[beta];
    if (cabs(vin) == 0.0) continue;

    tmp_v0[beta] += X->Sym->sym_diagonal[beta] * vin;
    prdct += X->Sym->sym_diagonal[beta] * conj(vin) * vin;

    for (p = 0; p < X->Sym->basis[beta].count; p++) {
      unsigned long int raw = X->Sym->basis[beta].raw_index[p];
      unsigned long int state = list_1[raw];
      double complex basis_coeff = X->Sym->basis[beta].coeff[p];
      unsigned int term;
      for (term = 0; term < X->Def.NExchangeCoupling; term++) {
        unsigned long int out_state;
        if (apply_exchange_halfspin(state,
                                    X->Def.ExchangeCoupling[term][0],
                                    X->Def.ExchangeCoupling[term][1],
                                    &out_state) == TRUE) {
          add_raw_transition(X, tmp_v0, tmp_v1, &prdct, beta, basis_coeff, out_state,
                             X->Def.ParaExchangeCoupling[term], vin);
        }
      }
    }
  }
  X->Large.prdct += prdct;
  return 0;
}
