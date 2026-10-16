#include <math.h>
#include <stdint.h>
#include "Common.h"
#include "symmetry_kondo.h"
#include "symmetry_kondo_terms.h"
#include "symmetry_terms.h"

/* Only a single Hamiltonian term is expanded here: at most four fermions.
 * No object in the normalizer scales with the Hilbert-space dimension. */
struct Fermion { int orbital, create; };
struct LocalBlock { int site; double coefficient[4]; };
struct Expansion {
  struct LocalBlock local[2];
  unsigned int nlocal;
  SymmetryKondoMonomialCallback callback;
  void *context;
};

static int expand_local(const struct Expansion *expansion, unsigned int site,
                        struct SymmetryKondoMonomial monomial)
{
  unsigned int choice;
  if (site == expansion->nlocal) {
    if (monomial.value == 0.0) return 0;
    return expansion->callback(&monomial, expansion->context) ? -1 : 0;
  }
  for (choice = 0; choice < 4; ++choice) {
    struct SymmetryKondoMonomial next = monomial;
    double coefficient = expansion->local[site].coefficient[choice];
    if (coefficient == 0.0) continue;
    next.value *= coefficient;
    if (choice) {
      unsigned int n = next.nlocal++;
      next.local_site[n] = expansion->local[site].site;
      next.local_out[n] = choice == 1 ? 0 : 1;
      next.local_in[n] = choice == 2 ? 0 : 1;
    }
    if (expand_local(expansion, site + 1, next)) return -1;
  }
  return 0;
}

static int normal_order(const struct Expansion *expansion,
                        const struct Fermion *ops, unsigned int count,
                        double complex value)
{
  struct SymmetryKondoMonomial monomial = {0};
  unsigned int i, j;
  /* a_p c^dagger_q = delta_pq - c^dagger_q a_p. Recursion strictly
   * decreases either the string length or its number of inversions. */
  for (i = 0; i + 1 < count; ++i) {
    if (!ops[i].create && ops[i+1].create) {
      struct Fermion next[4];
      memcpy(next, ops, count * sizeof(*ops));
      if (ops[i].orbital == ops[i+1].orbital) {
        for (j = i; j + 2 < count; ++j) next[j] = ops[j+2];
        if (normal_order(expansion, next, count - 2, value)) return -1;
      }
      memcpy(next, ops, count * sizeof(*ops));
      next[i] = ops[i+1]; next[i+1] = ops[i];
      return normal_order(expansion, next, count, -value);
    }
  }
  for (i = 0; i < count; ++i) {
    unsigned int *n = ops[i].create ? &monomial.ncreate : &monomial.nannihilate;
    int *orbitals = ops[i].create ? monomial.create : monomial.annihilate;
    if (*n >= 2) return -1;
    for (j = 0; j < *n; ++j) if (orbitals[j] == ops[i].orbital) return 0;
    j = (*n)++;
    while (j && orbitals[j-1] > ops[i].orbital) {
      orbitals[j] = orbitals[j-1]; --j; value = -value;
    }
    orbitals[j] = ops[i].orbital;
  }
  monomial.value = value;
  return expand_local(expansion, 0, monomial);
}

/* Evaluate a local product in the four-state Fock space, retaining only its
 * final singly occupied block. In particular, intermediate empty/double
 * states are essential for the crossed Exchange encoding. */
static void local_block(const struct Fermion *ops, unsigned int count,
                        double coefficient[4])
{
  double block[2][2] = {{0}};
  unsigned int input;
  for (input = 0; input < 2; ++input) {
    unsigned int digit = 1U << input, i = count;
    double sign = 1.0;
    while (i-- > 0) {
      unsigned int spin = (unsigned int)ops[i].orbital % 2;
      unsigned int mask = 1U << spin;
      if (((digit & mask) != 0) == ops[i].create) { sign = 0; break; }
      if (spin && (digit & 1U)) sign = -sign;
      digit ^= mask;
    }
    if (digit == 1U || digit == 2U) block[digit == 2U][input] = sign;
  }
  coefficient[0] = block[0][0];
  coefficient[1] = block[0][1];
  coefficient[2] = block[1][0];
  coefficient[3] = block[1][1] - block[0][0];
}

int CanonicalizeSymmetryKondoTerm(const struct DefineList *def,
    const struct SymmetryTerm *term, const int *permutation,
    SymmetryKondoMonomialCallback callback, void *context)
{
  struct Fermion ops[4];
  struct Expansion expansion = {0};
  unsigned int count, i, begin;
  double complex value;
  int permutation_sign;
  if (!def || !term || !callback || term->factors < 1 || term->factors > 2 ||
      !isfinite(creal(term->value)) || !isfinite(cimag(term->value)) ||
      ValidateSymmetryKondoSpace(def)) return -1;
  if (permutation && SymmetryKondoPermutationSign(def, permutation, &permutation_sign)) return -1;
  count = 2 * term->factors;
  value = term->value;
  expansion.callback = callback; expansion.context = context;
  for (i = 0; i < count; ++i) {
    int site = term->index[2*i], spin = term->index[2*i+1];
    if (site < 0 || (unsigned int)site >= def->Nsite || spin < 0 || spin > 1) return -1;
    if (permutation) site = permutation[site];
    ops[i].orbital = 2 * site + spin;
    ops[i].create = i % 2 == 0;
  }
  /* Stable insertion sort: local sites first in ascending order; all
   * conduction operators share a key and keep their original order. */
  for (i = 1; i < count; ++i) {
    struct Fermion op = ops[i];
    unsigned int j = i;
    int site = op.orbital / 2;
    unsigned int key = def->LocSpn[site] == LOCSPIN ? (unsigned int)site : def->Nsite;
    while (j) {
      int previous = ops[j-1].orbital / 2;
      unsigned int previous_key = def->LocSpn[previous] == LOCSPIN ? (unsigned int)previous : def->Nsite;
      if (previous_key <= key) break;
      ops[j] = ops[j-1]; --j; value = -value;
    }
    ops[j] = op;
  }
  begin = 0;
  while (begin < count && def->LocSpn[ops[begin].orbital / 2] == LOCSPIN) {
    unsigned int end = begin + 1;
    struct LocalBlock *block;
    while (end < count && ops[end].orbital / 2 == ops[begin].orbital / 2) ++end;
    /* Odd local products have a zero singly occupied block. This also
     * bounds the number of nonzero local blocks to two. */
    if ((end - begin) % 2) return 0;
    block = &expansion.local[expansion.nlocal++];
    block->site = ops[begin].orbital / 2;
    local_block(ops + begin, end - begin, block->coefficient);
    begin = end;
  }
  return normal_order(&expansion, ops + begin, count - begin, value);
}

struct Polynomial {
  struct SymmetryKondoMonomial *terms;
  size_t count, capacity;
  const struct DefineList *def;
  const int *permutation;
};
static int append(const struct SymmetryKondoMonomial *term, void *context)
{
  struct Polynomial *p = context;
  if (p->count == p->capacity) {
    size_t capacity = p->capacity ? 2 * p->capacity : 32;
    struct SymmetryKondoMonomial *next;
    if (capacity < p->capacity || capacity > SIZE_MAX / sizeof(*next)) return -1;
    next = realloc(p->terms, capacity * sizeof(*next));
    if (!next) return -1;
    p->terms = next; p->capacity = capacity;
  }
  p->terms[p->count++] = *term;
  return 0;
}
static int canonical_term(const struct SymmetryTerm *term, void *context)
{
  struct Polynomial *p = context;
  return CanonicalizeSymmetryKondoTerm(p->def, term, p->permutation, append, p);
}
static int compare(const void *a, const void *b)
{
  const struct SymmetryKondoMonomial *x = a, *y = b;
  unsigned int i;
#define CMP(a, b) do { if ((a) < (b)) return -1; if ((a) > (b)) return 1; } while (0)
  CMP(x->nlocal, y->nlocal);
  for (i = 0; i < x->nlocal; ++i) {
    CMP(x->local_site[i], y->local_site[i]);
    CMP(x->local_out[i], y->local_out[i]); CMP(x->local_in[i], y->local_in[i]);
  }
  CMP(x->ncreate, y->ncreate); CMP(x->nannihilate, y->nannihilate);
  for (i = 0; i < x->ncreate; ++i) CMP(x->create[i], y->create[i]);
  for (i = 0; i < x->nannihilate; ++i) CMP(x->annihilate[i], y->annihilate[i]);
#undef CMP
  return 0;
}
static int collect(struct Polynomial *p)
{
  size_t i, count = 0;
  p->count = 0;
  if (EnumerateSymmetryTerms(p->def, -1, canonical_term, p)) return -1;
  if (p->count) qsort(p->terms, p->count, sizeof(*p->terms), compare);
  for (i = 0; i < p->count; ++i) {
    if (count && compare(&p->terms[count-1], &p->terms[i]) == 0)
      p->terms[count-1].value += p->terms[i].value;
    else p->terms[count++] = p->terms[i];
  }
  p->count = 0;
  for (i = 0; i < count; ++i) {
    double complex value = p->terms[i].value;
    if (!isfinite(creal(value)) || !isfinite(cimag(value))) return -1;
    if (cabs(value) > 1e-10) p->terms[p->count++] = p->terms[i];
  }
  return 0;
}
int ValidateSymmetryKondoTerms(const struct DefineList *def)
{
  struct Polynomial base = {0}, mapped = {0};
  size_t i;
  unsigned int g, f;
  int error = 1;
  base.def = mapped.def = def;
  if (ValidateSymmetryKondoSpace(def) || (def->NSymTrans && !def->SymTrans) ||
      collect(&base)) goto done;
  for (i = 0; i < base.count; ++i) {
    const struct SymmetryKondoMonomial *m = &base.terms[i];
    int delta_sz = 0;
    for (f = 0; f < m->nlocal; ++f) delta_sz += 2 * (m->local_in[f] - m->local_out[f]);
    for (f = 0; f < m->ncreate; ++f) delta_sz += 1 - 2 * (m->create[f] % 2);
    for (f = 0; f < m->nannihilate; ++f) delta_sz -= 1 - 2 * (m->annihilate[f] % 2);
    if ((def->iCalcModel == Kondo && delta_sz) ||
        (def->iCalcModel == KondoNConserved && m->ncreate != m->nannihilate)) {
      fprintf(stdoutMPI, "Error: TransSym Kondo Hamiltonian does not conserve fixed quantum numbers.\n");
      goto done;
    }
  }
  for (g = 0; g < def->NSymTrans; ++g) {
    int sign;
    mapped.permutation = def->SymTrans[g];
    if (SymmetryKondoPermutationSign(def, mapped.permutation, &sign) || collect(&mapped)) goto done;
    if (base.count != mapped.count) goto invariant_error;
    for (i = 0; i < base.count; ++i)
      if (compare(&base.terms[i], &mapped.terms[i]) ||
          cabs(base.terms[i].value - mapped.terms[i].value) > 1e-10)
        goto invariant_error;
  }
  error = 0;
  goto done;
invariant_error:
  fprintf(stdoutMPI, "Error: TransSym Kondo Hamiltonian invariance failed for normalized terms under op %u.\n", g);
done:
  free(base.terms); free(mapped.terms);
  if (error) fprintf(stdoutMPI, "Error: invalid or unsupported TransSym Kondo Hamiltonian terms.\n");
  return error ? -1 : 0;
}
