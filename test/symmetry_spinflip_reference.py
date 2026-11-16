#!/usr/bin/env python3
"""Independent oracle for Spin(Sz=0)/SpinGC sectors with the spin flip Px.

Px = prod_i sigma^x_i exchanges up and down at every site with amplitude +1.
A word holds the local state of site i in bit i (1 = up, Sz = +1/2). Site
permutations and Px are explicit matrices built here from Kronecker products
of 2x2 Pauli matrices and from the images of basis words; the sector is the
range of Pi_k (1 + p Px)/2. No HPhi enumerator, canonicalizer or term builder
is used.
"""
from dataclasses import dataclass
import argparse
import numpy as np

MODELS = ('SpinGC', 'Spin')
# SpinFlipParity: 0 leaves Px out (control), +1/-1 select its eigenvalue.
PARITIES = (1, -1)
GROUPS = ('translation', 'identity', 'reflection', 'translation_reflection')

SIGMA_X = np.array([[0, 1], [1, 0]], dtype=complex)
E10 = np.array([[0, 0], [1, 0]], dtype=complex)   # |1><0|, raises Sz
E01 = E10.conj().T
SX, SY = (E10 + E01)/2, (E10 - E01)/(2j)
SZ = np.diag([-.5, .5]).astype(complex)

# Px-invariant couplings. Sy_i Sz_j has complex local matrix elements, so a
# conjugated Px image would break its invariance; PairLift is real.
SPINGC_DEFAULT = dict(jz=.37, jxy=.23, pair_lift=.17, hx=.83, yz=.19, j2=.11)
SPIN_DEFAULT = dict(jz=1., jxy=1., j2=.27)


@dataclass(frozen=True)
class Case:
    model: str
    sites: int
    parity: int
    momentum: int = 0
    group: str = 'translation'
    reflection: int = 1


@dataclass
class Reference:
    words: np.ndarray            # raw words of the model (Sz=0 only for Spin)
    basis: np.ndarray            # columns of the sector in the raw words
    representatives: np.ndarray  # smallest word of each G x Z2 orbit
    hamiltonian: np.ndarray      # raw H
    moments: dict[str, np.ndarray]
    operations: list             # (permutation matrix, flip, character) of G x Z2


def _validate(case):
    if case.model not in MODELS or case.group not in GROUPS:
        raise ValueError('unsupported model or group')
    if not 1 <= case.sites <= 10:
        raise ValueError('reference matrices must stay small')
    if case.model == 'Spin' and case.sites % 2:
        raise ValueError('a Spin Px sector needs Sz=0, so an even L')
    if case.parity not in (0, 1, -1) or case.reflection not in (1, -1):
        raise ValueError('parities must be 0, 1 or -1')
    if case.group in ('identity', 'reflection') and case.momentum != 0:
        raise ValueError('only translations carry a momentum')
    if case.group == 'translation_reflection' and 2*case.momentum % case.sites:
        raise ValueError('reflections keep only k=0 and k=pi one dimensional')
    if not 0 <= case.momentum < case.sites:
        raise ValueError('momentum index outside the translation group')


def tensor_at(nsite, site, matrix):
    """matrix at site, identities elsewhere; site 0 is the least significant bit."""
    result = np.ones((1, 1), complex)
    for pos in reversed(range(nsite)):
        result = np.kron(result, matrix if pos == site else np.eye(2))
    return result


def spin_operators(nsite):
    return {name: [tensor_at(nsite, i, m) for i in range(nsite)]
            for name, m in dict(plus=E10, minus=E01, x=SX, y=SY, z=SZ).items()}


def flip_matrix(nsite):
    """Px = sigma^x (x) ... (x) sigma^x on all 2^L words."""
    result = np.ones((1, 1), complex)
    for _ in range(nsite):
        result = np.kron(result, SIGMA_X)
    return result


def permutation_matrix(nsite, perm):
    """U|w> moves the state of site i to site perm[i]."""
    dim = 1 << nsite
    out = np.zeros((dim, dim))
    for w in range(dim):
        out[sum(((w >> i) & 1) << perm[i] for i in range(nsite)), w] = 1
    return out


def input_group(case: Case):
    """(permutation, character) rows of the TransSym file, before Px."""
    _validate(case)
    n = case.sites
    identity = list(range(n))
    mirror = [(-i) % n for i in range(n)]
    if case.group == 'identity':
        return [(identity, 1.)]
    if case.group == 'reflection':
        return [(identity, 1.), (mirror, float(case.reflection))]
    rows = [([(i+g) % n for i in range(n)], np.exp(-2j*np.pi*case.momentum*g/n))
            for g in range(n)]
    if case.group == 'translation_reflection':
        # R T^g, with chi(R T^g) = r chi(T^g); k=0 or pi keeps it one dimensional.
        rows += [([mirror[(i+g) % n] for i in range(n)], case.reflection*chi)
                 for g, (_, chi) in enumerate(rows[:n])]
    return rows


def raw_words(case: Case) -> np.ndarray:
    _validate(case)
    words = range(1 << case.sites)
    if case.model == 'Spin':
        words = [w for w in words if 2*bin(w).count('1') == case.sites]
    return np.array(list(words), dtype=np.uint64)


def operations(case: Case):
    """Explicit matrices of G, and of G x Z2 with chi(g,1) = chi(g) p."""
    n = case.sites
    ops = [(permutation_matrix(n, perm), 0, complex(chi)) for perm, chi in input_group(case)]
    if case.parity:
        px = flip_matrix(n)
        ops += [(u @ px, 1, case.parity*chi) for u, _, chi in ops]
    return ops


def sector_basis(case: Case, words: np.ndarray):
    """Columns P|r>/|P|r>|, P = sum_g conj(chi(g)) U_g/|G|, positive at r = min orbit."""
    ops = operations(case)
    dim = 1 << case.sites
    projector = sum(np.conj(chi)*u for u, _, chi in ops)/len(ops)
    index = np.asarray(words, dtype=np.int64)
    columns, representatives = [], []
    for r in index:
        # The orbit is read off the explicit operator matrices.
        orbit = {int(np.flatnonzero(abs(u[:, r]) > .5)[0]) for u, _, _ in ops}
        if min(orbit) != r:
            continue
        column = projector[:, r].copy()
        norm = np.linalg.norm(column)
        if norm < 1e-12:
            continue
        column *= np.conj(column[r])/abs(column[r])/norm
        columns.append(column[index])
        representatives.append(r)
    outside = np.setdiff1d(np.arange(dim), index)
    assert np.max(abs(projector[np.ix_(outside, index)]), initial=0) < 1e-14, \
        'the group leaves the raw space'
    basis = np.column_stack(columns) if columns else np.zeros((len(words), 0), complex)
    return basis, np.array(representatives, dtype=np.uint64), projector[np.ix_(index, index)]


def nearest_bonds(n):
    return [(i, (i+1) % n) for i in range(n)] if n > 2 else ([(0, 1)] if n == 2 else [])


def second_bonds(n):
    return [(i, (i+2) % n) for i in range(n)] if n > 4 else []


def physical_hamiltonian(case: Case, parameters=None, *, extra=()):
    """Translation/reflection invariant H of the model, as a 2^L matrix.

    extra: names of added terms that break Px ('hz', 'sy', 'dm', 'complex_pair').
    """
    n = case.sites
    s = spin_operators(n)
    p, m, x, y, z = (s[key] for key in ('plus', 'minus', 'x', 'y', 'z'))
    h = np.zeros((1 << n, 1 << n), complex)
    if case.model == 'SpinGC':
        q = dict(SPINGC_DEFAULT, **(parameters or {}))
        for i, j in nearest_bonds(n):
            h += (q['jz']*z[i]@z[j] + q['jxy']/2*(p[i]@m[j] + m[i]@p[j])
                  + q['pair_lift']*(p[i]@p[j] + m[i]@m[j])
                  + q['yz']*(y[i]@z[j] + z[i]@y[j]))
        for i in range(n):
            h -= q['hx']*x[i]
    else:
        q = dict(SPIN_DEFAULT, **(parameters or {}))
        for i, j in nearest_bonds(n):
            h += q['jz']*z[i]@z[j] + q['jxy']/2*(p[i]@m[j] + m[i]@p[j])
        for i, j in second_bonds(n):
            h += q['j2']*(z[i]@z[j] + (p[i]@m[j] + m[i]@p[j])/2)
    for name in extra:
        for i, j in nearest_bonds(n):
            if name == 'dm':
                h += .13j*(p[i]@m[j] - m[i]@p[j])
            elif name == 'complex_pair':
                h += .13j*(p[i]@p[j] - m[i]@m[j])
        for i in range(n):
            if name == 'hz':
                h += .21*z[i]
            elif name == 'sy':
                h += .21*y[i]
    return h


def moments(case: Case, words: np.ndarray) -> dict[str, np.ndarray]:
    sz = np.array([bin(int(w)).count('1') - case.sites/2 for w in words])
    return dict(Sz=sz, Sz2=sz**2)


def make_reference(case: Case, parameters=None, *, extra=(), empty=False) -> Reference:
    words = raw_words(case)
    basis, representatives, _ = sector_basis(case, words)
    index = np.asarray(words, dtype=np.int64)
    h = np.zeros((len(words), len(words)), complex) if empty else \
        physical_hamiltonian(case, parameters, extra=extra)[np.ix_(index, index)]
    return Reference(words, basis, representatives, h, moments(case, words), operations(case))


def sector_matrix(reference: Reference) -> np.ndarray:
    return reference.basis.conj().T @ reference.hamiltonian @ reference.basis


def raw_operator(case: Case, words: np.ndarray, factors) -> np.ndarray:
    """Ordered product of local units E_st = |s><t|, factors (site, s, site, t)."""
    op = np.eye(1 << case.sites, dtype=complex)
    for i, a, j, b in factors:
        assert i == j
        unit = np.zeros((2, 2), complex)
        unit[a, b] = 1
        op = op @ tensor_at(case.sites, i, unit)
    index = np.asarray(words, dtype=np.int64)
    return op[np.ix_(index, index)]


def check_reference(case: Case, ref: Reference, invariant=True):
    """Projector algebra of one sector: Q^dagger Q = 1, U_g Q = chi(g) Q, Px Q = p Q."""
    b, words = ref.basis, np.asarray(ref.words, dtype=np.int64)
    assert np.allclose(b.conj().T @ b, np.eye(b.shape[1]), atol=1e-12)
    _, _, projector = sector_basis(case, ref.words)
    assert np.allclose(projector @ projector, projector, atol=1e-12)
    assert np.allclose(projector, projector.conj().T, atol=1e-12)
    assert abs(np.trace(projector).real - b.shape[1]) < 1e-9
    assert np.allclose(projector @ b, b, atol=1e-12)
    for u, flip, chi in ref.operations:
        moved = u[np.ix_(words, words)] @ b
        assert np.allclose(moved, chi*b, atol=1e-12)
    if case.parity:
        px = flip_matrix(case.sites)[np.ix_(words, words)]
        assert np.allclose(px @ b, case.parity*b, atol=1e-12)
    rows = np.searchsorted(ref.words, ref.representatives)
    anchors = b[rows, np.arange(len(rows))]
    assert np.all(anchors.real > 0) and np.max(abs(anchors.imag), initial=0) < 1e-12
    h = ref.hamiltonian
    assert np.allclose(h, h.conj().T, atol=1e-13)
    if invariant:
        for u, _, _ in ref.operations:
            uu = u[np.ix_(words, words)]
            assert np.allclose(uu @ h, h @ uu, atol=1e-12)


def self_test():
    for model, sizes in (('SpinGC', (1, 2, 3, 4, 5)), ('Spin', (2, 4, 6))):
        for n in sizes:
            words = raw_words(Case(model, n, 0))
            full = np.linalg.eigvalsh(physical_hamiltonian(Case(model, n, 0))[
                np.ix_(words.astype(np.int64), words.astype(np.int64))])
            spectrum, total = [], 0
            for k in range(n):
                block = []
                for parity in PARITIES:
                    case = Case(model, n, parity, k)
                    ref = make_reference(case)
                    check_reference(case, ref)
                    total += ref.basis.shape[1]
                    block += list(np.linalg.eigvalsh(sector_matrix(ref)))
                control = make_reference(Case(model, n, 0, k))
                # D(k,+) + D(k,-) = D(k), and the two spectra make up the k block.
                assert len(block) == control.basis.shape[1]
                assert np.allclose(np.sort(block), np.linalg.eigvalsh(sector_matrix(control)), atol=1e-10)
                spectrum += block
            assert total == len(words)
            assert np.allclose(np.sort(spectrum), full, atol=1e-10)
    # Px commutes with every permutation, maps Sz to -Sz and keeps Sx.
    n = 4
    px, s = flip_matrix(n), spin_operators(n)
    for perm in ([1, 2, 3, 0], [0, 3, 2, 1]):
        u = permutation_matrix(n, perm)
        assert np.allclose(u @ px, px @ u)
    assert np.allclose(px @ s['z'][1] @ px, -s['z'][1])
    assert np.allclose(px @ s['y'][1] @ px, -s['y'][1])
    assert np.allclose(px @ s['x'][1] @ px, s['x'][1])
    # The added terms break Px but keep the translations.
    for model, extra in (('SpinGC', 'hz'), ('SpinGC', 'sy'), ('SpinGC', 'complex_pair'), ('Spin', 'dm'),
                         ('Spin', 'hz')):
        h = physical_hamiltonian(Case(model, n, 0), extra=(extra,))
        assert np.linalg.norm(px @ h - h @ px) > 1e-3, (model, extra)
    # The Neel word is fixed by T Px: at L=4 it survives only for (k=0,p=+1) or (k=2,p=-1).
    for k, parity, present in ((0, 1, True), (0, -1, False), (2, 1, False), (2, -1, True)):
        ref = make_reference(Case('SpinGC', 4, parity, k))
        assert (0b0101 in ref.representatives) == present, (k, parity)
    # Px alone and with reflections.
    for case in (Case('SpinGC', 3, 1, 0, 'identity'), Case('Spin', 4, -1, 0, 'reflection', -1),
                 Case('SpinGC', 4, -1, 2, 'translation_reflection', -1)):
        check_reference(case, make_reference(case))
    print('spin-flip reference self-test passed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true')
    if parser.parse_args().self_test:
        self_test()
