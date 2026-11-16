#!/usr/bin/env python3
"""Independent oracle for HubbardGC/SpinlessFermionGC sectors with optional N parity.

Fock words, translations and pair operators are built here from explicit
creation/annihilation strings; no HPhi enumerator or term builder is used.
"""
from dataclasses import dataclass
import argparse
import numpy as np
from symmetry_kondo_reference import _car, _parity
from symmetry_nconserved_reference import (DEFAULT as HUBBARD_DEFAULT, hopping,
                                           operator_matrix, orbital,
                                           physical_terms as hubbard_terms)

MODELS = ('HubbardGC', 'SpinlessFermionGC')
# ParticleParity: 0 all N, 1 even N, -1 odd N.
PARITIES = (0, 1, -1)
SPINLESS_DEFAULT = dict(t=.61*np.exp(.23j), mu=.17, v=.29)
# Translation-invariant pair amplitudes; each term is added with its adjoint.
PAIR_DEFAULT = dict(onsite=.23+.11j, singlet=.17-.06j, triplet=.09+.05j, spinless=.31-.12j)


@dataclass(frozen=True)
class Case:
    model: str
    sites: int
    parity: int
    momentum: int


@dataclass
class Reference:
    words: np.ndarray
    basis: np.ndarray
    representatives: np.ndarray
    hamiltonian: np.ndarray
    moments: dict[str, np.ndarray]


def width(model):
    return 2 if model == 'HubbardGC' else 1


def _validate(case):
    if case.model not in MODELS:
        raise ValueError('unsupported model')
    if not 1 <= case.sites <= (6 if case.model == 'HubbardGC' else 12):
        raise ValueError('reference words must stay small')
    if case.parity not in PARITIES:
        raise ValueError('parity must be 0, 1 or -1')
    if not 0 <= case.momentum < case.sites:
        raise ValueError('momentum index outside the translation group')


def raw_words(case: Case) -> np.ndarray:
    """All Fock words of the GC space, or those with (-1)^N equal to the parity."""
    _validate(case)
    words = range(1 << (width(case.model)*case.sites))
    if case.parity:
        words = [w for w in words if (-1)**bin(w).count('1') == case.parity]
    return np.array(list(words), dtype=np.uint64)


def translate_word(case: Case, word, g):
    """T_g c^dagger_{j,s} T_g^dagger = c^dagger_{j+g,s}, with its reordering sign."""
    w = width(case.model)
    occupied = [o for o in range(w*case.sites) if (word >> o) & 1]
    mapped = [w*((o//w+g) % case.sites)+o % w for o in occupied]
    return sum(1 << o for o in mapped), _parity(mapped)


def translation_matrix(case: Case, words: np.ndarray, g: int) -> np.ndarray:
    index = {int(w): i for i, w in enumerate(words)}
    out = np.zeros((len(words), len(words)))
    for source, word in enumerate(words):
        target, sign = translate_word(case, int(word), g)
        out[index[target], source] = sign
    return out


def sector_basis(case: Case, words: np.ndarray):
    """Orthonormal columns sum_g conj(chi(g)) T_g|r>, chi(g)=exp(-2 pi i k g/L)."""
    projector = sum(np.exp(2j*np.pi*case.momentum*g/case.sites) *
                    translation_matrix(case, words, g) for g in range(case.sites))/case.sites
    columns, representatives, covered = [], [], np.zeros(len(words), bool)
    for r in range(len(words)):
        if covered[r]:
            continue
        column = projector[:, r].copy()
        covered |= np.abs(column) > 1e-12
        norm = np.linalg.norm(column)
        if norm < 1e-12:
            continue
        column *= np.conj(column[r])/abs(column[r])/norm
        columns.append(column)
        representatives.append(words[r])
    basis = np.column_stack(columns) if columns else np.zeros((len(words), 0), complex)
    return basis, np.array(representatives, dtype=np.uint64)


def spinless_terms(sites, parameters=None):
    """-t c+_{j+1} c_j + h.c. - mu n_j + V n_j n_{j+1}, as (value, string) pairs."""
    p = parameters or SPINLESS_DEFAULT
    terms = []
    for j in range(sites):
        nxt = (j+1) % sites
        terms += [(-p['t'], [(1, nxt), (0, j)]), (-np.conj(p['t']), [(1, j), (0, nxt)]),
                  (-p['mu'], [(1, j), (0, j)]),
                  (p['v'], [(1, j), (0, j), (1, nxt), (0, nxt)])]
    return terms


def physical_terms(case: Case, parameters=None):
    if case.model == 'HubbardGC':
        return hubbard_terms(case.sites, parameters or HUBBARD_DEFAULT)
    return spinless_terms(case.sites, parameters)


def pair_operator(model, kind, i, si, j, sj):
    """HPhi AnomalousTerm/AnomalousG order: type 0 is c_i c_j, type 1 is c+_i c+_j."""
    a, b = (orbital(i, si), orbital(j, sj)) if width(model) == 2 else (i, j)
    return [(kind, a), (kind, b)]


def pair_rows(model, sites, parameters=None):
    """(type, i, si, j, sj, value) rows of a translation-invariant pair H.

    Each V c+_i c+_j is followed by its adjoint conj(V) c_j c_i, which is
    HPhi's adjacent Hermitian pair convention.
    """
    p = parameters or PAIR_DEFAULT
    rows = []

    def add(i, si, j, sj, value):
        rows.append((1, i, si, j, sj, complex(value)))
        rows.append((0, j, sj, i, si, complex(value).conjugate()))
    for j in range(sites):
        nxt = (j+1) % sites
        if model == 'HubbardGC':
            add(j, 0, j, 1, p['onsite'])
            add(j, 0, nxt, 1, p['singlet'])
            add(j, 1, nxt, 0, -p['singlet'])
            add(j, 0, nxt, 0, p['triplet'])
        else:
            add(j, 0, nxt, 0, p['spinless'])
    return rows


def equivalent_pair_rows(rows):
    """The same pair H as another valid input, no longer row-wise translated.

    The first Hermitian pair is split in halves, some swap their operators
    (V c+_i c+_j = -V c+_j c+_i) or put the adjoint row first, and the pair
    order is reversed.
    """
    def block(kind, i, si, j, sj, v, swap, adjoint_first):
        if swap:
            i, si, j, sj, v = j, sj, i, si, -v
        pair = [(kind, i, si, j, sj, v), (1-kind, j, sj, i, si, v.conjugate())]
        return pair[::-1] if adjoint_first else pair
    blocks = []
    for n, (kind, i, si, j, sj, v) in enumerate(rows[::2]):
        if n == 0:
            blocks += [block(kind, i, si, j, sj, v/2, False, False),
                       block(kind, i, si, j, sj, v/2, True, True)]
        else:
            blocks.append(block(kind, i, si, j, sj, v, n % 3 == 1, n % 2 == 0))
    return [row for pair in reversed(blocks) for row in pair]


def pair_terms(model, rows):
    return [(value, pair_operator(model, kind, i, si, j, sj))
            for kind, i, si, j, sj, value in rows]


def moments(case: Case, words: np.ndarray) -> dict[str, np.ndarray]:
    n = np.array([bin(int(w)).count('1') for w in words], float)
    if case.model == 'HubbardGC':
        up = np.array([sum((int(w) >> (2*j)) & 1 for j in range(case.sites)) for w in words], float)
        doublons = np.array([sum(((int(w) >> (2*j)) & 3) == 3 for j in range(case.sites))
                             for w in words], float)
        result = dict(N=n, D=doublons, Sz=.5*(2*up-n))
    else:
        result = dict(N=n, D=0*n, Sz=0*n)
    result.update({key+'2': value**2 for key, value in list(result.items())})
    return result


def make_reference(case: Case, parameters=None, *, pairs=None, empty=False) -> Reference:
    """pairs: None for no pair H, or rows as returned by pair_rows."""
    words = raw_words(case)
    basis, representatives = sector_basis(case, words)
    terms = [] if empty else list(physical_terms(case, parameters))
    if pairs:
        terms += pair_terms(case.model, pairs)
    h = operator_matrix(case, words, terms)
    return Reference(words, basis, representatives, h, moments(case, words))


def sector_matrix(reference: Reference) -> np.ndarray:
    return reference.basis.conj().T @ reference.hamiltonian @ reference.basis


def direct_unrank(rank, parity):
    """The production formula: append the lowest bit that sets the popcount parity."""
    if not parity:
        return rank
    return (rank << 1) | ((bin(rank).count('1') & 1) ^ (0 if parity > 0 else 1))


def self_test():
    # The direct unrank formula enumerates each parity class in increasing order.
    for m in range(1, 13):
        for parity in (1, -1):
            case_words = [w for w in range(1 << m) if (-1)**bin(w).count('1') == parity]
            assert [direct_unrank(r, parity) for r in range(1 << (m-1))] == case_words
    for model in MODELS:
        for sites in ((1, 2, 3, 4) if model == 'HubbardGC' else (1, 2, 3, 4, 5)):
            full = make_reference(Case(model, sites, 0, 0), pairs=pair_rows(model, sites))
            h_full = full.hamiltonian
            rewritten = make_reference(Case(model, sites, 0, 0),
                                       pairs=equivalent_pair_rows(pair_rows(model, sites)))
            assert np.allclose(rewritten.hamiltonian, h_full, rtol=0, atol=1e-14)
            n = full.moments['N']
            # Pair terms change N by 2 and keep (-1)^N.
            odd_block = h_full[np.ix_(n % 2 == 0, n % 2 == 1)]
            assert np.max(abs(odd_block), initial=0) < 1e-14
            # A spinless p-wave pair cancels on rings shorter than three sites.
            if model == 'HubbardGC' or sites >= 3:
                assert np.count_nonzero(abs(h_full[np.ix_(n == 2, n == 0)]) > 1e-12)
            spectra = {}
            for parity in PARITIES:
                total, spectrum = 0, []
                for k in range(sites):
                    case = Case(model, sites, parity, k)
                    ref = make_reference(case, pairs=pair_rows(model, sites))
                    assert len(ref.words) == (1 << (width(model)*sites)) >> (parity != 0)
                    assert np.all(np.diff(ref.words.astype(object)) > 0)
                    h, b = ref.hamiltonian, ref.basis
                    assert np.allclose(h, h.conj().T, atol=1e-13)
                    assert np.allclose(b.conj().T @ b, np.eye(b.shape[1]), atol=1e-12)
                    t1 = translation_matrix(case, ref.words, 1)
                    assert np.allclose(t1 @ b, np.exp(-2j*np.pi*k/sites)*b, atol=1e-12)
                    assert np.allclose(np.linalg.matrix_power(t1, sites), np.eye(len(ref.words)))
                    assert np.allclose(t1 @ h, h @ t1, atol=1e-12)
                    rows = np.searchsorted(ref.words, ref.representatives)
                    anchors = b[rows, np.arange(len(rows))]
                    assert np.all(anchors.real > 0) and np.max(abs(anchors.imag), initial=0) < 1e-12
                    total += b.shape[1]
                    spectrum += list(np.linalg.eigvalsh(sector_matrix(ref)))
                assert total == len(raw_words(Case(model, sites, parity, 0)))
                spectra[parity] = np.sort(spectrum)
            assert np.allclose(spectra[0], np.linalg.eigvalsh(h_full), atol=1e-10)
            assert np.allclose(np.sort(np.r_[spectra[1], spectra[-1]]), spectra[0], atol=1e-10)
        # F_ij = <c_i c_j> = -F_ji and <c+_j c+_i> = conj(F_ij) for any state.
        case = Case(model, 3, 0, 0)
        words = raw_words(case)
        rng = np.random.default_rng(20261019)
        psi = rng.normal(size=len(words)) + 1j*rng.normal(size=len(words))
        psi /= np.linalg.norm(psi)

        def expect(kind, i, si, j, sj):
            op = operator_matrix(case, words, [(1, pair_operator(model, kind, i, si, j, sj))])
            return np.vdot(psi, op @ psi)
        spin = 1 if model == 'HubbardGC' else 0
        f = expect(0, 0, 0, 1, spin)
        assert abs(f) > 1e-3
        assert abs(f + expect(0, 1, spin, 0, 0)) < 1e-13
        assert abs(expect(1, 1, spin, 0, 0) - np.conj(f)) < 1e-13
    # c+_{0up} c+_{0dn}|0> = |0up,0dn>, i.e. the bit-0 creator acts last.
    case = Case('HubbardGC', 1, 0, 0)
    op = operator_matrix(case, raw_words(case), [(1, pair_operator('HubbardGC', 1, 0, 0, 0, 1))])
    assert op[3, 0] == 1 and np.count_nonzero(op) == 1
    print('fermion GC reference self-test passed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true')
    if parser.parse_args().self_test:
        self_test()
