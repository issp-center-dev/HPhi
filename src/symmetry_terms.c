#include <limits.h>
#include <math.h>
#include <stdint.h>
#include "Common.h"
#include "symmetry_terms.h"

int SymmetryUsesExtendedTerms(const struct DefineList *def)
{
  return def->iCalcModel == tJ || def->EDNChemi || def->NInterAll || def->NInterAll_Diagonal ||
      def->NInterAll_OffDiagonal || def->NPairHopping ||
      (def->iCalcModel == Hubbard && (def->NCoulombInter || def->NHundCoupling ||
       def->NExchangeCoupling || def->NIsingCoupling)) ||
      (def->iCalcModel == Spin && (def->NTransfer ||
       def->NCoulombInter != def->NIsingCoupling ||
       def->NHundCoupling != def->NIsingCoupling));
}

static int emit(const struct DefineList *def, int kind, unsigned int factors,
                const int *index, double complex value,
                SymmetryTermCallback callback, void *context)
{
  struct SymmetryTerm term = {0};
  unsigned int i, j;
  int diagonal = 1;
  term.factors = factors;
  term.value = value;
  if (!isfinite(creal(value)) || !isfinite(cimag(value))) return -1;
  for (i = 0; i < 4 * factors; i += 2) {
    if (index[i] < 0 || (unsigned int)index[i] >= def->Nsite ||
        index[i+1] < 0 || index[i+1] > (def->iCalcModel == SpinlessFermion ? 0 : 1))
      return -1;
    term.index[i] = index[i];
    term.index[i+1] = index[i+1];
  }
  if (def->iCalcModel == Spin) {
    for (i = 0; i < factors; ++i)
      if (index[4*i] != index[4*i+2]) return -1;
  }
  /* A product is diagonal precisely when the occupation change at each
   * orbital (or the local Spin state) vanishes. */
  for (i = 0; i < factors; ++i) {
    int balance = 0;
    for (j = 0; j < factors; ++j) {
      if (index[4*i] == index[4*j] && index[4*i+1] == index[4*j+1]) ++balance;
      if (index[4*i] == index[4*j+2] && index[4*i+1] == index[4*j+3]) --balance;
    }
    if (balance != 0) diagonal = 0;
  }
  if (kind >= 0 && kind != !diagonal) return 0;
  return callback(&term, context);
}

int EnumerateSymmetryTerms(const struct DefineList *def, int kind,
                          SymmetryTermCallback callback, void *context)
{
  unsigned int p;
  int a, b, s, t;
  if (!def || !callback || kind < -1 || kind > 1 || def->Nsite == 0 ||
      (def->iCalcModel != Spin && def->iCalcModel != SpinlessFermion &&
       def->iCalcModel != Hubbard && def->iCalcModel != tJ) || def->iFlgGeneralSpin ||
      def->NNBodyInterAll || def->NAnomalousTerm || def->NPairLiftCoupling ||
      def->NIsingCoupling > def->NCoulombInter ||
      def->NIsingCoupling > def->NHundCoupling) return -1;
#define EMIT(n, v, ...) do { int ix[] = {__VA_ARGS__}; \
  if (emit(def, kind, n, ix, v, callback, context)) return -1; } while (0)
#define CHECK_STORAGE(n, rows, values) do { \
  if ((n) && (!(rows) || !(values))) return -1; \
  for (p = 0; p < (n); ++p) if (!(rows)[p]) return -1; } while (0)
  CHECK_STORAGE(def->EDNTransfer, def->EDGeneralTransfer, def->EDParaGeneralTransfer);
  for (p = 0; p < def->EDNTransfer; ++p)
    if (emit(def, kind, 1, def->EDGeneralTransfer[p],
             -def->EDParaGeneralTransfer[p], callback, context)) return -1;
  if (def->EDNChemi && (!def->EDChemi || !def->EDSpinChemi || !def->EDParaChemi))
    return -1;
  for (p = 0; p < def->EDNChemi; ++p) {
    a = def->EDChemi[p]; s = def->EDSpinChemi[p];
    EMIT(1, -def->EDParaChemi[p], a,s,a,s);
  }
  CHECK_STORAGE(def->NInterAll_Diagonal, def->InterAll_Diagonal, def->ParaInterAll_Diagonal);
  for (p = 0; p < def->NInterAll_Diagonal; ++p) {
    a = def->InterAll_Diagonal[p][0]; s = def->InterAll_Diagonal[p][1];
    b = def->InterAll_Diagonal[p][2]; t = def->InterAll_Diagonal[p][3];
    EMIT(2, def->ParaInterAll_Diagonal[p], a,s,a,s,b,t,b,t);
  }
  CHECK_STORAGE(def->NInterAll_OffDiagonal, def->InterAll_OffDiagonal, def->ParaInterAll_OffDiagonal);
  for (p = 0; p < def->NInterAll_OffDiagonal; ++p)
    if (emit(def, kind, 2, def->InterAll_OffDiagonal[p],
             def->ParaInterAll_OffDiagonal[p], callback, context)) return -1;
  CHECK_STORAGE(def->NCoulombIntra, def->CoulombIntra, def->ParaCoulombIntra);
  if (def->NCoulombIntra && def->iCalcModel != Hubbard && def->iCalcModel != tJ) return -1;
  for (p = 0; p < def->NCoulombIntra; ++p) {
    a = def->CoulombIntra[p][0];
    EMIT(2, def->ParaCoulombIntra[p], a,0,a,0,a,1,a,1);
  }
  CHECK_STORAGE(def->NCoulombInter, def->CoulombInter, def->ParaCoulombInter);
  for (p = 0; p < def->NCoulombInter; ++p) {
    a = def->CoulombInter[p][0]; b = def->CoulombInter[p][1];
    for (s = 0; s <= (def->iCalcModel == SpinlessFermion ? 0 : 1); ++s)
      for (t = 0; t <= (def->iCalcModel == SpinlessFermion ? 0 : 1); ++t)
        EMIT(2, def->ParaCoulombInter[p], a,s,a,s,b,t,b,t);
  }
  CHECK_STORAGE(def->NHundCoupling, def->HundCoupling, def->ParaHundCoupling);
  if (def->NHundCoupling && def->iCalcModel == SpinlessFermion) return -1;
  for (p = 0; p < def->NHundCoupling; ++p) {
    a = def->HundCoupling[p][0]; b = def->HundCoupling[p][1];
    for (s = 0; s < 2; ++s) EMIT(2, -def->ParaHundCoupling[p], a,s,a,s,b,s,b,s);
  }
  CHECK_STORAGE(def->NExchangeCoupling, def->ExchangeCoupling, def->ParaExchangeCoupling);
  if (def->NExchangeCoupling && def->iCalcModel == SpinlessFermion) return -1;
  for (p = 0; p < def->NExchangeCoupling; ++p) {
    a = def->ExchangeCoupling[p][0]; b = def->ExchangeCoupling[p][1];
    if (a == b) continue; /* Existing Exchange semantics. */
    for (s = 0; s < 2; ++s) {
      if (def->iCalcModel == Spin) {
        EMIT(2, def->ParaExchangeCoupling[p], a,s,a,1-s,b,1-s,b,s);
      } else {
        EMIT(2, def->ParaExchangeCoupling[p], a,s,b,s,b,1-s,a,1-s);
      }
    }
  }
  CHECK_STORAGE(def->NPairHopping, def->PairHopping, def->ParaPairHopping);
  if (def->NPairHopping && def->iCalcModel != Hubbard && def->iCalcModel != tJ) return -1;
  for (p = 0; p < def->NPairHopping; ++p) {
    a = def->PairHopping[p][0]; b = def->PairHopping[p][1];
    EMIT(2, def->ParaPairHopping[p], a,0,b,0,a,1,b,1);
  }
#undef CHECK_STORAGE
#undef EMIT
  return 0;
}

int ApplySymmetryTerm(const struct DefineList *def,
                      const struct SymmetryTerm *term, unsigned long state,
                      unsigned long *out, double complex *value)
{
  int f;
  double sign = 1;
  unsigned int width = (def && (def->iCalcModel == Hubbard || def->iCalcModel == tJ)) ? 2U : 1U;
  if (!def || !term || !out || !value || term->factors < 1 || term->factors > 2 ||
      def->Nsite == 0 || def->Nsite > CHAR_BIT * sizeof(state) / width ||
      (def->iCalcModel != Spin && def->iCalcModel != SpinlessFermion &&
       def->iCalcModel != Hubbard && def->iCalcModel != tJ)) return -1;
  for (f = 0; f < (int)(4*term->factors); f += 2)
    if (term->index[f] < 0 || (unsigned int)term->index[f] >= def->Nsite ||
        term->index[f+1] < 0 || term->index[f+1] > (def->iCalcModel == SpinlessFermion ? 0 : 1))
      return -1;
  if (def->iCalcModel == Spin)
    for (f = 0; f < (int)term->factors; ++f)
      if (term->index[4*f] != term->index[4*f+2]) return -1;
  for (f = (int)term->factors - 1; f >= 0; --f) {
    const int *x = term->index + 4*f;
    if (def->iCalcModel == Spin) {
      unsigned long mask = 1UL << x[0];
      if ((int)((state >> x[0]) & 1UL) != x[3]) return 0;
      state = (state & ~mask) | ((unsigned long)x[1] << x[0]);
    } else {
      int op;
      for (op = 1; op >= 0; --op) { /* annihilation, then creation */
        unsigned int orbital = width * (unsigned int)x[2*op] + (unsigned int)x[2*op+1];
        unsigned long mask = 1UL << orbital, below = state & (mask - 1UL);
        if (((state & mask) != 0) != (op == 1)) return 0;
        while (below) { sign = -sign; below &= below - 1UL; }
        state ^= mask;
      }
    }
  }
  /* Project the final state onto the physical tJ Hilbert space. */
  if (def->iCalcModel == tJ && (state & (state >> 1U) & (ULONG_MAX / 3UL))) return 0;
  *out = state;
  *value = sign * term->value;
  return 1;
}

/* Canonical polynomials are used only at validation time. Fermions use
 * normal-ordered creation/annihilation strings; Spin uses the independent
 * local basis {I, E01, E10, E11}, with E00 = I - E11. */
struct Monomial { int key[7]; double complex value; };
struct Polynomial {
  struct Monomial *terms;
  size_t count, capacity;
  const struct DefineList *def;
  const int *permutation;
};

static int append_monomial(struct Polynomial *poly, const int key[7], double complex value)
{
  if (poly->count == poly->capacity) {
    size_t cap = poly->capacity ? poly->capacity * 2 : 32;
    struct Monomial *next;
    if (cap < poly->capacity || cap > SIZE_MAX / sizeof(*next)) return -1;
    next = realloc(poly->terms, cap * sizeof(*next));
    if (!next) return -1;
    poly->terms = next; poly->capacity = cap;
  }
  memcpy(poly->terms[poly->count].key, key, 7 * sizeof(int));
  poly->terms[poly->count++].value = value;
  return 0;
}

static int canonical_term(const struct SymmetryTerm *input, void *context)
{
  struct Polynomial *poly = context;
  int x[8], key[7] = {0}, n = (int)input->factors, i;
  double complex value = input->value;
  memcpy(x, input->index, sizeof(x));
  if (poly->permutation)
    for (i = 0; i < 4*n; i += 2) x[i] = poly->permutation[x[i]];
  if (poly->def->iCalcModel == Spin) {
    if (n == 2 && x[0] == x[4]) {
      if (x[3] != x[5]) return 0;
      x[3] = x[7]; n = 1;
    }
    if (n == 2 && x[0] > x[4])
      for (i = 0; i < 4; ++i) { int tmp = x[i]; x[i] = x[i+4]; x[i+4] = tmp; }
    /* At most two E00 factors: expand their identity/E11 choices. */
    for (i = 0; i < (1 << n); ++i) {
      int f, count = 0, valid = 1;
      double complex coefficient = value;
      memset(key, 0, sizeof(key));
      for (f = 0; f < n; ++f) {
        int zero = x[4*f+1] == 0 && x[4*f+3] == 0;
        int choose = (i >> f) & 1;
        if (!zero && choose) { valid = 0; break; }
        if (zero && !choose) continue;
        key[1+3*count] = x[4*f];
        key[2+3*count] = zero ? 1 : x[4*f+1];
        key[3+3*count] = zero ? 1 : x[4*f+3];
        if (zero) coefficient = -coefficient;
        ++count;
      }
      key[0] = count;
      if (valid && append_monomial(poly, key, coefficient)) return -1;
    }
    return 0;
  } else {
    int width = (poly->def->iCalcModel == Hubbard || poly->def->iCalcModel == tJ) ? 2 : 1;
    int a = width*x[0]+x[1], b = width*x[2]+x[3];
    if (n == 1) {
      key[0] = 1; key[1] = a; key[2] = b;
      return append_monomial(poly, key, value);
    } else {
      int c = width*x[4]+x[5], d = width*x[6]+x[7], tmp;
      if (b == c) {
        key[0] = 1; key[1] = a; key[2] = d;
        if (append_monomial(poly, key, value)) return -1;
      }
      if (a == c || b == d) return 0;
      /* P O P vanishes if a normal-ordered string creates or annihilates
       * a doublon. Retain any contraction emitted above before discarding it. */
      if (poly->def->iCalcModel == tJ && (a / 2 == c / 2 || b / 2 == d / 2)) return 0;
      value = -value;
      if (a > c) { tmp = a; a = c; c = tmp; value = -value; }
      if (b > d) { tmp = b; b = d; d = tmp; value = -value; }
      memset(key, 0, sizeof(key));
      key[0] = 2; key[1] = a; key[2] = c; key[3] = b; key[4] = d;
      return append_monomial(poly, key, value);
    }
  }
}

static int compare_monomials(const void *a, const void *b)
{
  const struct Monomial *x = a, *y = b;
  int i;
  for (i = 0; i < 7; ++i) {
    if (x->key[i] < y->key[i]) return -1;
    if (x->key[i] > y->key[i]) return 1;
  }
  return 0;
}

static int collect_polynomial(struct Polynomial *poly)
{
  size_t i, next = 0;
  if (EnumerateSymmetryTerms(poly->def, -1, canonical_term, poly)) return -1;
  if (poly->count) qsort(poly->terms, poly->count, sizeof(*poly->terms), compare_monomials);
  for (i = 0; i < poly->count; ++i) {
    if (next && compare_monomials(&poly->terms[next-1], &poly->terms[i]) == 0)
      poly->terms[next-1].value += poly->terms[i].value;
    else poly->terms[next++] = poly->terms[i];
  }
  poly->count = 0;
  for (i = 0; i < next; ++i) {
    if (!isfinite(creal(poly->terms[i].value)) || !isfinite(cimag(poly->terms[i].value))) return -1;
    if (cabs(poly->terms[i].value) > 1e-10) poly->terms[poly->count++] = poly->terms[i];
  }
  return 0;
}

int ValidateSymmetryTerms(const struct DefineList *def)
{
  struct Polynomial base = {0}, mapped = {0};
  size_t i;
  unsigned int g;
  int error = 1;
  base.def = mapped.def = def;
  if (collect_polynomial(&base)) goto done;
  for (i = 0; i < base.count; ++i) {
    const int *key = base.terms[i].key;
    int delta = 0, f;
    if (def->iCalcModel == Spin)
      for (f = 0; f < key[0]; ++f) delta += key[2+3*f] - key[3+3*f];
    else if (def->iCalcModel == Hubbard || def->iCalcModel == tJ)
      for (f = 0; f < key[0]; ++f) delta += key[1+f]%2 - key[1+key[0]+f]%2;
    if (delta) {
      fprintf(stdoutMPI, "Error: TransSym Hamiltonian does not conserve fixed Sz/Nup/Ndown.\n");
      goto done;
    }
  }
  for (g = 0; g < def->NSymTrans; ++g) {
    mapped.count = 0;
    mapped.permutation = def->SymTrans[g];
    if (collect_polynomial(&mapped)) goto done;
    if (base.count != mapped.count) goto invariant_error;
    for (i = 0; i < base.count; ++i)
      if (compare_monomials(&base.terms[i], &mapped.terms[i]) ||
          cabs(base.terms[i].value - mapped.terms[i].value) > 1e-10)
        goto invariant_error;
  }
  error = 0;
  goto done;
invariant_error:
  fprintf(stdoutMPI, "Error: TransSym Hamiltonian invariance failed for normalized terms under op %u.\n", g);
done:
  free(base.terms); free(mapped.terms);
  if (error) fprintf(stdoutMPI, "Error: invalid or unsupported TransSym Hamiltonian terms.\n");
  return error ? -1 : 0;
}
