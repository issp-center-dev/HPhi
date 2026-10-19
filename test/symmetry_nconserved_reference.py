#!/usr/bin/env python3
"""Independent fixed-N Fock reference for HubbardNConserved/tJNConserved sectors.

Raw words come from orbital combinations, never from HPhi's unrank. Operators
act through the CAR string and permutation parity of the Kondo reference with
no local spins: orbital 2*site+spin, spin 0 = up, HPhi's word order.
"""
from dataclasses import dataclass
from itertools import combinations
import argparse
import numpy as np
from symmetry_kondo_reference import _car, _parity

MODELS = ('HubbardNConserved', 'tJNConserved')


@dataclass(frozen=True)
class Case:
    model: str
    sites: int
    n: int
    momentum: int


@dataclass
class Reference:
    words: np.ndarray
    basis: np.ndarray
    representatives: np.ndarray
    hamiltonian: np.ndarray
    moments: dict[str, np.ndarray]


def _validate(case):
    if case.model not in MODELS:
        raise ValueError('unsupported model')
    if not 1 <= case.sites <= 16:
        raise ValueError('reference words must stay small')
    # N=0 is a valid physical block for references, not an HPhi input.
    if not 0 <= case.n <= (case.sites if case.model == 'tJNConserved' else 2*case.sites):
        raise ValueError('particle number outside the space')
    if not 0 <= case.momentum < case.sites:
        raise ValueError('momentum index outside the translation group')


def _doubly_occupied(word, sites):
    return any(((word >> (2*j)) & 3) == 3 for j in range(sites))


def raw_words(case: Case) -> np.ndarray:
    _validate(case)
    words = [sum(1 << orbital for orbital in occupied)
             for occupied in combinations(range(2*case.sites), case.n)]
    if case.model == 'tJNConserved':
        words = [w for w in words if not _doubly_occupied(w, case.sites)]
    return np.array(sorted(words), dtype=np.uint64)


def translate_word(sites, word, g):
    """T_g c^dagger_{j,s} T_g^dagger = c^dagger_{j+g,s}, with its reordering sign."""
    occupied = [orbital for orbital in range(2*sites) if (word >> orbital) & 1]
    mapped = [2*((orbital//2+g) % sites)+orbital % 2 for orbital in occupied]
    return sum(1 << orbital for orbital in mapped), _parity(mapped)


def translation_matrix(case: Case, words: np.ndarray, g: int) -> np.ndarray:
    index = {int(w): i for i, w in enumerate(words)}
    out = np.zeros((len(words), len(words)))
    for source, word in enumerate(words):
        target, sign = translate_word(case.sites, int(word), g)
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


def operator_matrix(case: Case, words: np.ndarray, terms) -> np.ndarray:
    """Sum of coefficient * ordered string, each (create, orbital) acting from the right.

    The whole string acts in the full Fock space; tJ projects only its result.
    """
    index = {int(w): i for i, w in enumerate(words)}
    out = np.zeros((len(words), len(words)), complex)
    for source, word in enumerate(words):
        for coefficient, string in terms:
            target, amplitude = _car(int(word), string)
            if amplitude and target in index:
                out[index[target], source] += coefficient*amplitude
    return out


def orbital(site, spin):
    return 2*site+spin


def hopping(i, s, j, t, value):
    """value c^dagger_{i,s} c_{j,t}."""
    return (value, [(1, orbital(i, s)), (0, orbital(j, t))])


def density_density(i, s, j, t, value):
    return (value, [(1, orbital(i, s)), (0, orbital(i, s)),
                    (1, orbital(j, t)), (0, orbital(j, t))])


def physical_terms(sites, parameters):
    """Translation-invariant spin-mixing chain; every off-diagonal term with its adjoint."""
    p = parameters
    terms = []
    for j in range(sites):
        nxt = (j+1) % sites
        for s in (0, 1):
            terms += [hopping(nxt, s, j, s, -p['t']), hopping(j, s, nxt, s, -np.conj(p['t']))]
            terms.append(hopping(j, s, j, s, -p['mu'] + (-.5 if s else .5)*p['hz']))
        terms += [hopping(nxt, 0, j, 1, -p['soc']), hopping(j, 1, nxt, 0, -np.conj(p['soc'])),
                  hopping(nxt, 1, j, 0, -p['soc2']), hopping(j, 0, nxt, 1, -np.conj(p['soc2'])),
                  hopping(j, 0, j, 1, -p['hx']), hopping(j, 1, j, 0, -np.conj(p['hx']))]
        terms.append(density_density(j, 0, j, 1, p['u']))
        for s in (0, 1):
            for t in (0, 1):
                terms.append(density_density(j, s, nxt, t, p['v']))
        # S+_j S+_{j+1} and its adjoint change Sz by +-2 and keep N.
        terms += [(p['w'], [(1, orbital(j, 0)), (0, orbital(j, 1)),
                            (1, orbital(nxt, 0)), (0, orbital(nxt, 1))]),
                  (np.conj(p['w']), [(1, orbital(nxt, 1)), (0, orbital(nxt, 0)),
                                     (1, orbital(j, 1)), (0, orbital(j, 0))])]
    return terms


DEFAULT = dict(t=.73*np.exp(.17j), mu=.23, hz=.11, soc=.21+.07j, soc2=-.09+.04j,
               hx=.27-.05j, u=.41, v=.13, w=.08-.03j)


def make_reference(case: Case, parameters=None, *, empty=False) -> Reference:
    words = raw_words(case)
    basis, representatives = sector_basis(case, words)
    terms = [] if empty else physical_terms(case.sites, parameters or DEFAULT)
    h = operator_matrix(case, words, terms)
    n_up = np.array([sum((int(w) >> (2*j)) & 1 for j in range(case.sites)) for w in words], float)
    n_down = np.array([sum((int(w) >> (2*j+1)) & 1 for j in range(case.sites)) for w in words], float)
    doublons = np.array([sum(((int(w) >> (2*j)) & 3) == 3 for j in range(case.sites))
                         for w in words], float)
    moments = dict(N=n_up+n_down, D=doublons, Sz=.5*(n_up-n_down))
    moments.update({key+'2': value**2 for key, value in list(moments.items())})
    return Reference(words, basis, representatives, h, moments)


def sector_matrix(reference: Reference) -> np.ndarray:
    return reference.basis.conj().T @ reference.hamiltonian @ reference.basis


def raw_dimension(model, sites, n):
    from math import comb
    return comb(2*sites, n) if model == 'HubbardNConserved' else comb(sites, n)*2**n


def self_test():
    rng = np.random.default_rng(20261019)
    for model in MODELS:
        for sites in (2, 3, 4):
            top = sites if model == 'tJNConserved' else 2*sites
            for n in range(top+1):
                total = 0
                spectrum = []
                for k in range(sites):
                    case = Case(model, sites, n, k)
                    ref = make_reference(case)
                    assert len(ref.words) == raw_dimension(model, sites, n)
                    assert np.all(np.diff(ref.words.astype(object)) > 0)
                    h, b = ref.hamiltonian, ref.basis
                    assert np.allclose(h, h.conj().T, atol=1e-13)
                    assert np.allclose(b.conj().T @ b, np.eye(b.shape[1]), atol=1e-12)
                    # Character, closure and H symmetry of the group action.
                    t1 = translation_matrix(case, ref.words, 1)
                    assert np.allclose(t1 @ b, np.exp(-2j*np.pi*k/sites)*b, atol=1e-12)
                    assert np.allclose(np.linalg.matrix_power(t1, sites), np.eye(len(ref.words)))
                    assert np.allclose(t1 @ h, h @ t1, atol=1e-12)
                    rows = np.searchsorted(ref.words, ref.representatives)
                    anchors = b[rows, np.arange(len(rows))]
                    assert np.all(anchors.real > 0) and np.max(abs(anchors.imag), initial=0) < 1e-12
                    total += b.shape[1]
                    spectrum += list(np.linalg.eigvalsh(sector_matrix(ref)))
                assert total == raw_dimension(model, sites, n)
                full = np.linalg.eigvalsh(make_reference(Case(model, sites, n, 0)).hamiltonian)
                assert np.allclose(np.sort(spectrum), full, atol=1e-10)
    # The spin-mixing terms connect different Sz in the same N block.
    ref = make_reference(Case('HubbardNConserved', 3, 2, 0))
    sz = ref.moments['Sz']
    assert np.count_nonzero(abs(ref.hamiltonian[np.ix_(sz == 1, sz == -1)]) > 1e-12)
    # tJ projects only the result: c+_{1dn} c_{0dn} c+_{0up} c_{1up} passes through
    # the doubly occupied site 0 and maps |0dn,1up> to -|0up,1dn>.
    case = Case('tJNConserved', 2, 2, 0)
    words = raw_words(case)
    assert not any(_doubly_occupied(int(w), 2) for w in words)
    product = operator_matrix(case, words, [(1, [(1, 3), (0, 1), (1, 0), (0, 2)])])
    source, target = list(words).index(6), list(words).index(9)
    assert product[target, source] == -1 and np.count_nonzero(product) == 1
    print('nconserved reference self-test passed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true')
    if parser.parse_args().self_test:
        self_test()
