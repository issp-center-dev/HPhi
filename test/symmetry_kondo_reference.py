#!/usr/bin/env python3
"""Independent spin tensor / conduction Fock reference for Kondo sectors."""
from dataclasses import dataclass
from typing import Callable
import argparse
import numpy as np


@dataclass(frozen=True)
class Case:
    model: str
    cells: int
    layout: str
    ncond: int | None
    sz2: int | None
    momentum: int


@dataclass
class Reference:
    words: np.ndarray
    basis: np.ndarray
    representatives: np.ndarray
    apply_h: Callable[[np.ndarray], np.ndarray]
    moments: dict[str, np.ndarray]


# Local operators act on |up>, |down>; only the conduction space carries CAR.
_SP = np.array([[0, 1], [0, 0]], dtype=complex)
_SM = _SP.T


def _sites(case: Case):
    if case.cells < 1 or case.layout not in ('block', 'alternating'):
        raise ValueError('positive cell count and block/alternating layout required')
    p = case.cells
    if case.layout == 'block':
        return tuple(range(p)), tuple(range(p, 2*p))
    return tuple(range(0, 2*p, 2)), tuple(range(1, 2*p, 2))


def _validate_case(case):
    _sites(case)
    if case.model not in ('Kondo', 'KondoNConserved', 'KondoGC'):
        raise ValueError('unsupported Kondo model')
    if 4*case.cells > 64:
        raise ValueError('reference words must fit uint64')
    if case.model == 'KondoGC':
        if case.ncond is not None or case.sz2 is not None:
            raise ValueError('GC has no fixed particle or spin quantities')
    elif case.ncond is None or not 0 <= case.ncond <= 2*case.cells:
        raise ValueError('canonical reference requires a valid ncond')
    if (case.model == 'Kondo') != (case.sz2 is not None):
        raise ValueError('only Kondo fixes sz2')


@dataclass
class _Space:
    words: np.ndarray
    spins: list[tuple[int, ...]]
    focks: list[int]
    etas: np.ndarray
    lookup: dict[tuple[tuple[int, ...], int], int]


def _embedding(case, spins, fock):
    local, conduction = _sites(case)
    word = sum((1 << spin) << (2*site) for site, spin in zip(local, spins))
    parity = 0
    for j, site in enumerate(conduction):
        digit = (fock >> (2*j)) & 3
        word |= digit << (2*site)
        parity += digit.bit_count() * sum(site < l for l in local)
    return word, (-1 if parity % 2 else 1)


def _space(case):
    # Enumerate the tensor product directly, never filter HPhi raw words.
    from itertools import product
    _validate_case(case)
    entries = []
    for spins in product((0, 1), repeat=case.cells):
        for fock in range(1 << (2*case.cells)):
            if case.ncond is not None and fock.bit_count() != case.ncond:
                continue
            if case.sz2 is not None:
                sz2 = sum(1-2*s for s in spins)
                sz2 += sum(((fock >> (2*j)) & 1)-((fock >> (2*j+1)) & 1)
                           for j in range(case.cells))
                if sz2 != case.sz2:
                    continue
            word, eta = _embedding(case, spins, fock)
            entries.append((word, spins, fock, eta))
    entries.sort()
    return _Space(np.array([e[0] for e in entries], dtype=np.uint64),
                  [e[1] for e in entries], [e[2] for e in entries],
                  np.array([e[3] for e in entries], dtype=np.int8),
                  {(e[1], e[2]): i for i, e in enumerate(entries)})


def _parity(mapped):
    return (-1)**sum(a > b for i, a in enumerate(mapped) for b in mapped[i+1:])


def _translation_data(case, space, g):
    p = case.cells
    destinations = np.empty(len(space.words), dtype=np.intp)
    amplitudes = np.empty(len(space.words), dtype=np.int8)
    for source, (spins, fock) in enumerate(zip(space.spins, space.focks)):
        rotated = tuple(spins[(j-g) % p] for j in range(p))
        occupied = [orb for orb in range(2*p) if (fock >> orb) & 1]
        mapped = [2*((orb//2+g) % p)+orb % 2 for orb in occupied]
        target = space.lookup[rotated, sum(1 << orb for orb in mapped)]
        destinations[source] = target
        amplitudes[source] = _parity(mapped)*int(space.etas[source])*int(space.etas[target])
    return destinations, amplitudes


def _translate(case, vector, g):
    space = _space(case)
    dest, amp = _translation_data(case, space, g)
    out = np.empty_like(vector, dtype=complex)
    out[dest] = amp.reshape((-1,) + (1,)*(vector.ndim-1))*vector
    return out


def _check_eta(case):
    # Independent full-word permutation is only a test oracle for U T U†,
    # never the reference's physical translation implementation.
    space = _space(case)
    local, conduction = _sites(case)
    for g in range(case.cells):
        perm = {sites[j]: sites[(j+g) % case.cells]
                for sites in (local, conduction) for j in range(case.cells)}
        local_sign = _parity([perm[j] for j in local])
        dest, amp = _translation_data(case, space, g)
        for source, word in enumerate(space.words):
            occupied = [orb for orb in range(4*case.cells) if (int(word) >> orb) & 1]
            mapped = [2*perm[orb//2]+orb % 2 for orb in occupied]
            assert int(space.words[dest[source]]) == sum(1 << orb for orb in mapped)
            assert amp[source] == local_sign*_parity(mapped)


def _car(fock, operators):
    """Apply a left-to-right conduction creation/annihilation string."""
    amplitude = 1
    for create, orbital in reversed(operators):
        occupied = (fock >> orbital) & 1
        if occupied == create:
            return fock, 0
        amplitude *= (-1)**((fock & ((1 << orbital)-1)).bit_count())
        fock ^= 1 << orbital
    return fock, amplitude


def _tensor_action(spins, fock, spin_factors, fermions):
    # Each matrix used here has at most one nonzero entry per input column.
    current = list(spins)
    amplitude = 1+0j
    for site, matrix in reversed(spin_factors):
        column = matrix[:, current[site]]
        targets = np.flatnonzero(column)
        if not len(targets):
            return spins, fock, 0j
        if len(targets) != 1:
            raise ValueError('split sums into elementary spin matrices')
        target = int(targets[0])
        amplitude *= column[target]
        current[site] = target
    fock, sign = _car(fock, fermions)
    return tuple(current), fock, amplitude*sign


def _transition_action(space, terms, diagonal=None):
    sources, targets, amplitudes = [], [], []
    for source, (spins, fock) in enumerate(zip(space.spins, space.focks)):
        for coefficient, spin_factors, fermions in terms:
            out_spins, out_fock, amplitude = _tensor_action(spins, fock, spin_factors, fermions)
            target = space.lookup.get((out_spins, out_fock))
            if amplitude and target is not None:
                sources.append(source)
                targets.append(target)
                amplitudes.append(coefficient*amplitude*int(space.etas[source])*int(space.etas[target]))
    sources = np.asarray(sources, dtype=np.intp)
    targets = np.asarray(targets, dtype=np.intp)
    amplitudes = np.asarray(amplitudes, dtype=complex)
    if diagonal is None:
        diagonal = np.zeros(len(space.words))

    def apply(vector):
        vector = np.asarray(vector, dtype=complex)
        if vector.shape != (len(space.words),):
            raise ValueError('expected one raw-space vector')
        out = diagonal*vector
        np.add.at(out, targets, amplitudes*vector[sources])
        return out
    return apply


def apply_operator(case: Case, factors: tuple[tuple[int, int, int, int], ...],
                   vector: np.ndarray) -> np.ndarray:
    """Apply P O P for ordered bilinears, in sorted raw-word coordinates.

    Each local bilinear must be onsite E_ab; all other bilinears must join
    conduction sites. Crossed local bilinears (an internal HPhi encoding of
    Exchange) are deliberately rejected: use their physical spin equivalent.
    Intermediate spin/particle sectors are retained until the whole product
    has acted, so e.g. S+ S- is evaluated without intermediate projection.
    """
    space = _space(case)
    local, conduction = _sites(case)
    spin_factors, fermions = [], []
    for i, a, j, b in factors:
        if a not in (0, 1) or b not in (0, 1):
            raise ValueError('spin index must be 0 or 1')
        if i in local and j == i:
            matrix = np.zeros((2, 2), complex)
            matrix[a, b] = 1
            spin_factors.append((local.index(i), matrix))
        elif i in conduction and j in conduction:
            fermions.extend(((1, 2*conduction.index(i)+a),
                             (0, 2*conduction.index(j)+b)))
        else:
            raise ValueError('only onsite local or conduction bilinears are supported')
    return _transition_action(space, [(1, spin_factors, fermions)])(vector)


def _hamiltonian(case, space, empty):
    p = case.cells
    diagonal = np.zeros(len(space.words))
    if empty:
        return _transition_action(space, [], diagonal)
    t, phi, u, v, jz, jp, coupling, hz, mu = .73, .17, .41, .13, 1.10, .83, .19, .11, .23
    for row, (spins, fock) in enumerate(zip(space.spins, space.focks)):
        sz = [.5-s for s in spins]
        digits = [(fock >> (2*j)) & 3 for j in range(p)]
        n = [d.bit_count() for d in digits]
        sc = [.5*((d & 1)-((d >> 1) & 1)) for d in digits]
        diagonal[row] = sum(u*(digits[j] == 3) + v*n[j]*n[(j+1) % p]
                            + jz*sz[j]*sc[j] + coupling*sz[j]*sz[(j+1) % p]
                            - hz*(sz[j]+.7*sc[j])-mu*n[j] for j in range(p))
    terms = []
    for j in range(p):
        nxt = (j+1) % p
        for sigma in (0, 1):
            terms.append((-t*np.exp(1j*phi), [], [(1, 2*nxt+sigma), (0, 2*j+sigma)]))
            terms.append((-t*np.exp(-1j*phi), [], [(1, 2*j+sigma), (0, 2*nxt+sigma)]))
        terms.append((jp/2, [(j, _SP)], [(1, 2*j+1), (0, 2*j)]))
        terms.append((jp/2, [(j, _SM)], [(1, 2*j), (0, 2*j+1)]))
        if case.model != 'Kondo':
            terms.extend([(-.27/2, [(j, _SP)], []), (-.27/2, [(j, _SM)], []),
                          (.08j, [], [(1, 2*j), (0, 2*j+1)]),
                          (-.08j, [], [(1, 2*j+1), (0, 2*j)]),
                          (.09, [(j, _SP), (nxt, _SP)], []),
                          (.09, [(j, _SM), (nxt, _SM)], [])])
    return _transition_action(space, terms, diagonal)


def make_reference(case: Case, *, empty: bool = False) -> Reference:
    space = _space(case)
    p = case.cells
    translations = [_translation_data(case, space, g) for g in range(p)]
    visited = set()
    columns, representatives = [], []
    for source, word in enumerate(space.words):
        if source in visited:
            continue
        column = {}
        for g, (dest, amp) in enumerate(translations):
            target = int(dest[source])
            visited.add(target)
            column[target] = column.get(target, 0j) + np.exp(2j*np.pi*case.momentum*g/p)*amp[source]/p
        norm = np.sqrt(sum(abs(value)**2 for value in column.values()))
        if norm < 1e-12:
            continue
        phase = column[source]/abs(column[source])
        columns.append({i: value/(norm*phase) for i, value in column.items()})
        representatives.append(word)
    basis = np.zeros((len(space.words), len(columns)), dtype=np.complex128)
    for j, column in enumerate(columns):
        for i, value in column.items():
            basis[i, j] = value
    nc = np.array([f.bit_count() for f in space.focks], dtype=float)
    doublons = np.array([sum(((f >> (2*j)) & 3) == 3 for j in range(p))
                        for f in space.focks], dtype=float)
    sz = np.array([sum(.5-s for s in spins) +
                   .5*sum(((f >> (2*j)) & 1)-((f >> (2*j+1)) & 1) for j in range(p))
                   for spins, f in zip(space.spins, space.focks)])
    moments = {'N': nc+p, 'D': doublons, 'Sz': sz, 'Ncond': nc}
    moments.update({key+'2': value**2 for key, value in list(moments.items())})
    return Reference(space.words, basis, np.asarray(representatives, dtype=np.uint64),
                     _hamiltonian(case, space, empty), moments)


def sector_hamiltonian(reference: Reference) -> np.ndarray:
    basis = reference.basis
    out = np.empty((basis.shape[1], basis.shape[1]), dtype=complex)
    for column in range(basis.shape[1]):
        out[:, column] = basis.conj().T @ reference.apply_h(basis[:, column])
    return out


def _self_test():
    # Literal counts catch wrong fermion signs, parity correction or character.
    dimensions = {3: ([13]*3, [40]*3, [176, 168, 168]),
                  4: ([36]*4, [108, 116, 108, 116], [1024]*4)}
    raw_dimensions = {3: [39, 120, 512], 4: [144, 448, 4096]}
    models = ('Kondo', 'KondoNConserved', 'KondoGC')
    rng = np.random.default_rng(20261010)
    errors = dict(norm=0., hermitian=0., commutator=0., character=0.)
    for p in (3, 4):
        for layout in ('block', 'alternating'):
            for m, model in enumerate(models):
                for k in range(p):
                    case = Case(model, p, layout,
                                None if m == 2 else 2,
                                (1 if p == 3 else 0) if m == 0 else None, k)
                    ref = make_reference(case)
                    assert len(ref.words) == raw_dimensions[p][m]
                    assert ref.basis.shape[1] == dimensions[p][m][k]
                    assert ref.words.dtype == np.uint64
                    assert ref.basis.dtype == np.complex128
                    assert np.all(np.diff(ref.words) > 0)
                    assert np.all(np.diff(ref.representatives) > 0)
                    rows = np.searchsorted(ref.words, ref.representatives)
                    anchors = ref.basis[rows, np.arange(len(rows))]
                    assert np.max(abs(anchors.imag), initial=0) < 1e-13
                    assert np.all(anchors.real > 0)
                    # Orbit columns have disjoint support, so their normalized
                    # norms and absence of overlapping rows prove B^dagger B=I.
                    norm_error = np.max(abs(np.sum(abs(ref.basis)**2, axis=0)-1), initial=0)
                    errors['norm'] = max(errors['norm'], norm_error)
                    assert norm_error < 1e-12
                    assert np.all(np.count_nonzero(abs(ref.basis) > 1e-12, axis=1) <= 1)
                    v = rng.normal(size=len(ref.words)) + 1j*rng.normal(size=len(ref.words))
                    w = rng.normal(size=len(ref.words)) + 1j*rng.normal(size=len(ref.words))
                    hermitian_error = abs(np.vdot(w, ref.apply_h(v)) - np.vdot(ref.apply_h(w), v))
                    errors['hermitian'] = max(errors['hermitian'], hermitian_error)
                    assert hermitian_error < 1e-10
                    project = lambda x: ref.basis @ (ref.basis.conj().T @ x)
                    commutator_error = np.linalg.norm(ref.apply_h(project(v))-project(ref.apply_h(v)))
                    errors['commutator'] = max(errors['commutator'], commutator_error)
                    assert commutator_error < 1e-10
                    for g in range(p):
                        character_error = np.max(abs(_translate(case, ref.basis, g)-
                                           np.exp(-2j*np.pi*k*g/p)*ref.basis), initial=0)
                        errors['character'] = max(errors['character'], character_error)
                        assert character_error < 1e-12
                        for h in range(p):
                            assert np.allclose(_translate(case, _translate(case, v, h), g),
                                               _translate(case, v, (g+h)%p))
                    _check_eta(case)
                    n = ref.moments['Ncond']
                    assert np.array_equal(ref.moments['N'], n+p)
                    for key in ('N', 'D', 'Sz', 'Ncond'):
                        assert np.array_equal(ref.moments[key+'2'], ref.moments[key]**2)
                    if p == 3 and k == 1:
                        hs = sector_hamiltonian(ref)
                        assert np.allclose(hs, hs.conj().T, atol=1e-12)
                print(f'P={p} {layout} {model}: raw/sector dimensions, eta, group, H passed')
    for layout in ('block', 'alternating'):
        for model in models:
            # All-up local spins and the conduction vacuum are invariant,
            # including inside the larger NConserved and GC spaces.
            case = Case(model, 4, layout, None if model == 'KondoGC' else 0,
                        4 if model == 'Kondo' else None, 0)
            ref = make_reference(case, empty=True)
            local, _ = _sites(case)
            word = sum(1 << (2*i) for i in local)
            v = np.zeros(len(ref.words), complex)
            v[np.searchsorted(ref.words, word)] = 1
            assert np.allclose(_translate(case, v, 1), v)
            assert np.allclose(ref.basis @ (ref.basis.conj().T @ v), v)
            assert not np.any(ref.apply_h(v))
        for k in (1, 2):
            case = Case('Kondo', 3, layout, 0, 1, k)
            ref = make_reference(case, empty=True)
            local, _ = _sites(case)
            base = sum(1 << (2*i) for i in local)
            fourier = np.zeros(3, complex)
            for j, site in enumerate(local):
                word = base ^ (3 << (2*site))
                fourier[np.searchsorted(ref.words, word)] = np.exp(2j*np.pi*k*j/3)/np.sqrt(3)
            assert np.allclose(_translate(case, fourier, 1), np.exp(-2j*np.pi*k/3)*fourier)
            assert np.allclose(ref.basis @ (ref.basis.conj().T @ fourier), fourier)
        # Product ordering matters: S+ S- is |up><up|, S- S+ is |down><down|.
        case = Case('KondoGC', 1, layout, None, None, 0)
        ref = make_reference(case, empty=True)
        l, c = (a[0] for a in _sites(case))
        v = np.ones(8, complex)
        plus = (l, 0, l, 1)
        minus = (l, 1, l, 0)
        up = apply_operator(case, (plus, minus), v)
        down = apply_operator(case, (minus, plus), v)
        assert np.array_equal(up + down, v)
        assert np.vdot(up, down) == 0
        assert not np.any(apply_operator(case, ((l, 0, l, 0), (l, 1, l, 1)), v))
        hop = (c, 0, c, 1)
        assert np.array_equal(apply_operator(case, (hop,), apply_operator(case, ((c,1,c,0),), v)),
                              apply_operator(case, ((c,0,c,0),), v) -
                              apply_operator(case, ((c,0,c,0),(c,1,c,1)), v))
    case = Case('KondoGC', 1, 'block', None, None, 0)
    ref = make_reference(case)
    h = np.column_stack([ref.apply_h(v) for v in np.eye(8)])
    hop = -1.46*np.cos(.17)
    want = np.diag([-.0075, .1025, hop+.129, hop-.311,
                    hop-.344, hop+.316, 2*hop+.4625, 2*hop+.5725]).astype(complex)
    for i in (0, 2, 4, 6):
        want[i, i+1] = want[i+1, i] = -.135
    for i in (2, 3):
        want[i, i+2] = .08j
        want[i+2, i] = -.08j
    want[3, 4] = want[4, 3] = .415
    assert np.allclose(h, want, atol=1e-13)
    try:
        apply_operator(case, ((0, 0, 1, 0),), np.ones(8))
    except ValueError:
        pass
    else:
        raise AssertionError('crossed local bilinear accepted')
    # An intermediate Sz-changing local action must not be projected away.
    case = Case('Kondo', 1, 'block', 1, 0, 0)
    ref = make_reference(case, empty=True)
    assert np.array_equal(apply_operator(case, ((0,0,0,1),(0,1,0,0)), np.ones(2)), [0, 1])
    print('maximum absolute errors:', {k: f'{v:.3e}' for k, v in errors.items()})
    print('reference self-test passed')


def _drive_action(case: Case, kind: str, amplitude: float):
    """Physical tensor/CAR drive; amplitude is A(t) for a Peierls drive."""
    space = _space(case)
    diagonal = np.zeros(len(space.words))
    terms = []
    if kind == 'onebody' and case.model == 'Kondo':
        diagonal = np.array([.21*amplitude*sum(.5-s for s in spins)
                             for spins in space.spins])
    elif kind == 'twobody' and case.model == 'Kondo':
        diagonal = np.array([.17*amplitude*sum((.5-spins[j]) *
                              .5*(((f >> (2*j)) & 1)-((f >> (2*j+1)) & 1))
                              for j in range(case.cells))
                             for spins, f in zip(space.spins, space.focks)])
    elif kind in ('onebody', 'twobody'):
        for j in range(case.cells):
            if kind == 'onebody':
                terms.extend([(-.105j*amplitude, [(j, _SP)], []),
                              (.105j*amplitude, [(j, _SM)], [])])
            else:
                for s, sz in ((0, .5), (1, -.5)):
                    for spin in (_SP, _SM):
                        terms.append((.085*sz*amplitude, [(j, spin)],
                                      [(1, 2*j+s), (0, 2*j+s)]))
    elif kind == 'laser':
        forward = -.73*np.exp(.17j)*(np.exp(-1j*amplitude)-1)
        for j in range(case.cells):
            nxt = (j+1) % case.cells
            for s in (0, 1):
                terms.extend([(forward, [], [(1, 2*nxt+s), (0, 2*j+s)]),
                              (forward.conjugate(), [], [(1, 2*j+s), (0, 2*nxt+s)])])
    else:
        raise ValueError(kind)
    static = _hamiltonian(case, space, False)
    drive = _transition_action(space, terms, diagonal)
    return lambda vector: static(vector)+drive(vector)


def drive_hamiltonian(case: Case, kind: str, amplitude: float) -> np.ndarray:
    """Return static H plus the physical onebody/twobody/laser sector drive."""
    reference = make_reference(case)
    reference.apply_h = _drive_action(case, kind, amplitude)
    return sector_hamiltonian(reference)


def make_empty_rank_reference(model: str) -> Reference:
    """Odd sector: P2/Nc2, or two local spins and one fixed conduction site.

    The GC fixture uses an explicit local mask (sites 0, 1); it is not a
    cell chain. Construct its tensor Hamiltonian before embedding raw words.
    """
    if model != 'KondoGC':
        return make_reference(Case(model, 2, 'block', 2,
                                   0 if model == 'Kondo' else None, 1))
    from itertools import product
    entries = sorted((sum((1 << s) << (2*l) for l, s in enumerate(spins)) | (f << 4), spins, f)
                     for spins in product((0, 1), repeat=2) for f in range(4))
    space = _Space(np.array([e[0] for e in entries], dtype=np.uint64),
                   [e[1] for e in entries], [e[2] for e in entries],
                   np.ones(16, dtype=np.int8),
                   {(s, f): i for i, (_, s, f) in enumerate(entries)})
    columns, representatives = [], []
    for i, (word, spins, fock) in enumerate(entries):
        j = space.lookup[spins[::-1], fock]
        if i >= j:
            continue
        vector = np.zeros(16, complex)
        vector[i], vector[j] = 1/np.sqrt(2), -1/np.sqrt(2)
        columns.append(vector)
        representatives.append(word)
    diagonal = np.array([.6*sum(.5-s for s in spins)*.5*((f&1)-((f>>1)&1))
                          + .2*(.5-spins[0])*(.5-spins[1])
                          + .41*(f == 3)-.23*f.bit_count() for _, spins, f in entries])
    terms = [(.1, [(0, _SP), (1, _SM)], []), (.1, [(0, _SM), (1, _SP)], []),
             (.085j, [], [(1, 0), (0, 1)]), (-.085j, [], [(1, 1), (0, 0)])]
    for l in (0, 1):
        terms.extend([(.3, [(l, _SP)], [(1, 1), (0, 0)]),
                      (.3, [(l, _SM)], [(1, 0), (0, 1)]),
                      (-.065, [(l, _SP)], []), (-.065, [(l, _SM)], [])])
    nc = np.array([f.bit_count() for f in space.focks], float)
    moments = dict(N=nc+2, Ncond=nc,
                   D=np.array([f == 3 for f in space.focks], float),
                   Sz=np.array([sum(.5-s for s in spins)+.5*((f&1)-((f>>1)&1))
                                for spins, f in zip(space.spins, space.focks)]))
    moments.update({key+'2': value**2 for key, value in list(moments.items())})
    return Reference(space.words, np.column_stack(columns),
                     np.array(representatives, dtype=np.uint64),
                     _transition_action(space, terms, diagonal), moments)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        _self_test()
