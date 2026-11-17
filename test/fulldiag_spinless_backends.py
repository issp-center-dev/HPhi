"""Raw spinless FullDiag states/observables against an independent Fock oracle.

The existing fulldiag_spinless test retains its all-up Hubbard/sector-union
comparison. Here every vector and complex correlation is checked independently
for LAPACK and the configured distributed backends, including ExpecMode 0/1/2.
"""
from pathlib import Path
import os
import sys

import numpy as np
from symmetry_spingc_common import definition, run_case
from symmetry_fulldiag_common import (launcher_from_env, solver_backends, suite_root,
                                      require_available_backends, require_backends)
from verify_fulldiag_eigenvectors import read_vectors

L = 5
HOPS = [(0, 1, 1.0+.3j), (1, 2, .8), (2, 3, 1.1-.2j), (3, 4, .9+.1j)]
TRANS = [(i, j, z) for a, b, t in HOPS for i, j, z in ((a, b, t), (b, a, np.conj(t)))]
TRANS += [(2, 2, .35)]
COULOMB = [(i, i+1, v) for i, v in enumerate((.7, .4, .9, .6))]
PAIRS = [(0, 2, 1, 3, .45+.15j), (3, 1, 2, 0, .45-.15j)]
ONE = [(i, 0, j, 0) for i in range(L) for j in range(L)]
TWO = [(i, 0, j, 0, k, 0, m, 0) for i, j, k, m, _ in PAIRS]


def operator(words, orbitals):
    """Ordered c^dagger_i c_j ... via CAR, with ascending bit-word basis."""
    index = {word: j for j, word in enumerate(words)}
    out = np.zeros((len(words), len(words)), complex)
    actions = [(pos % 2 == 0, orbital) for pos, orbital in enumerate(orbitals)]
    for col, word in enumerate(words):
        sign = 1
        for create, orbital in reversed(actions):
            mask = 1 << orbital
            if bool(word & mask) == create:
                break
            sign *= (-1)**((word & (mask-1)).bit_count())
            word ^= mask
        else:
            if word in index:
                out[index[word], col] = sign
    return out


def reference(n):
    words = [w for w in range(1 << L) if n is None or w.bit_count() == n]
    one = np.array([operator(words, row[::2]) for row in ONE])
    two = np.array([operator(words, row[::2]) for row in TWO])
    h = sum(-t*operator(words, (i, j)) for i, j, t in TRANS)
    h += sum(v*operator(words, (i, i, j, j)) for i, j, v in COULOMB)
    h += sum(t*op for (*_, t), op in zip(PAIRS, two))
    np.testing.assert_allclose(h, h.conj().T, atol=1e-14, rtol=0)
    values, vectors = np.linalg.eigh(h)
    # This fixture is deliberately nondegenerate: per-state references are
    # unambiguous across solvers, while a vector's overall phase is irrelevant.
    if len(values) > 1:
        assert np.min(np.diff(values)) > 1e-5
    return h, values, vectors, one, two, np.array([w.bit_count() for w in words])


def write_input(path, n, solver, mode):
    path.mkdir()
    calc = dict(CalcType=2, CalcModel=8 if n is None else 7, Solver=solver,
                NGPU=0, ExpecMode=mode, OutputEigenVec=1)
    (path/'calc.def').write_text(''.join(f'{k} {v}\n' for k, v in calc.items()))
    mod = dict(Nsite=L, Lanczos_max=100, initial_iv=-1, exct=1, LanczosEps=14,
               LanczosTarget=2, LargeValue=12., NumAve=5, ExpecInterval=20)
    if n is not None:
        mod['Ncond'] = n
    (path/'mod.def').write_text('====\nModel_Parameters 0\n====\n====\n====\n'
        'CDataFileHead zvo\nCParaFileHead zqp\n====\n'
        + ''.join(f'{k} {v}\n' for k, v in mod.items()))
    definition(path, 'loc.def', [(i, 0) for i in range(L)], count=0)
    families = {
        'Trans': [(i, 0, j, 0, complex(t).real, complex(t).imag) for i, j, t in TRANS],
        'CoulombInter': COULOMB,
        'NBodyInterAll': [(2, i, 0, j, 0, k, 0, m, 0, complex(t).real, complex(t).imag)
                         for i, j, k, m, t in PAIRS],
        'OneBodyG': ONE, 'TwoBodyG': TWO,
    }
    names = ['CalcMod calc.def', 'ModPara mod.def', 'LocSpin loc.def']
    for key, rows in families.items():
        definition(path, key+'.def', rows)
        names.append(f'{key} {key}.def')
    (path/'sym.def').write_text('\n'.join(names)+'\n')


def check_output(path, ref):
    h, expected, oracle_vectors, one, two, counts = ref
    dim = len(h)
    output = path/'output'
    spectrum = np.loadtxt(output/'Eigenvalue.dat', ndmin=2)
    assert spectrum.shape == (dim, 2) and np.isfinite(spectrum).all(), path
    np.testing.assert_array_equal(spectrum[:, 0], np.arange(dim))
    np.testing.assert_allclose(spectrum[:, 1], expected, atol=1e-9, rtol=0)
    vectors = np.asarray(read_vectors(output, 'zvo')).T
    assert vectors.shape == (dim, dim) and np.isfinite(vectors).all(), path
    np.testing.assert_allclose(vectors.conj().T @ vectors, np.eye(dim), atol=1e-10, rtol=0)
    residual = np.max(np.linalg.norm(h @ vectors-vectors*expected, axis=0))
    assert residual < 1e-10*max(1., np.linalg.norm(h, 2)), (path, residual)
    expected_n = counts @ abs(oracle_vectors)**2
    np.testing.assert_allclose(counts @ abs(vectors)**2, expected_n, atol=1e-10, rtol=0)
    phys_files = list(output.glob('zvo_phys*.dat'))
    assert len(phys_files) == 1, phys_files
    phys = np.loadtxt(phys_files[0], skiprows=1, ndmin=2)
    assert phys.shape == (dim, 5) and np.isfinite(phys).all(), path
    # phys uses six decimals; vectors and complex Green functions have higher precision.
    np.testing.assert_allclose(phys[:, 0], expected, atol=5.1e-7, rtol=0)
    np.testing.assert_allclose(phys[:, 1], expected_n, atol=5.1e-7, rtol=0)
    np.testing.assert_allclose(phys[:, 2:], 0., atol=1e-12, rtol=0)
    max_green = 0.
    max_imaginary = 0.
    for state in range(dim):
        psi = oracle_vectors[:, state]
        for suffix, rows, ops in (('cisajs', ONE, one), ('cisajscktalt', TWO, two)):
            data = np.loadtxt(output/f'zvo_{suffix}_eigen{state}.dat', ndmin=2)
            assert data.shape == (len(rows), len(rows[0])+2) and np.isfinite(data).all()
            np.testing.assert_array_equal(data[:, :-2], rows)
            actual = data[:, -2] + 1j*data[:, -1]
            want = np.array([np.vdot(psi, op @ psi) for op in ops])
            max_green = max(max_green, float(np.max(abs(actual-want))))
            max_imaginary = max(max_imaginary, float(np.max(abs(want.imag))))
            np.testing.assert_allclose(actual, want, atol=1e-9, rtol=0)
            if suffix == 'cisajs':
                np.testing.assert_allclose(actual.reshape(L, L).trace(), expected_n[state],
                                           atol=1e-9, rtol=0)
    if dim > 1:
        assert max_imaginary > 1e-3, 'fixture must exercise complex correlations'
    print(f'{path.name}: {dim} states, residual={residual:.3e}, '
          f'max Green error={max_green:.3e}; norm/N/finite checks passed', flush=True)


def main():
    hphi = Path(sys.argv[1]).resolve()
    launcher, ranks = launcher_from_env()
    backends = solver_backends(launcher, ranks)
    require_available_backends(backends)
    root = suite_root('spinless-backends')
    env = dict(os.environ)
    env.pop('HPHI_SYMMETRY_BASIS_LAYOUT', None)
    exercised = {}
    for n in (2, None, L):
        name = 'gc' if n is None else f'n{n}'
        ref = reference(n)
        for solver, label, command in backends:
            for mode in ((0,) if solver == 0 else (0, 1, 2)):
                path = root/f'{name}_{label}_mode{mode}'
                write_input(path, n, solver, mode)
                if solver == 3 and len(ref[0]) == 1:
                    # ELPA explicitly rejects a matrix smaller than its process
                    # grid. This expected rejection never counts as execution.
                    run_case(path, hphi, 'grid-rejection', command, env,
                             expected_error='smaller than the process grid')
                    print(f'{path.name}: ELPA dimension-one guard passed', flush=True)
                    continue
                text = run_case(path, hphi, 'run', command, env)
                assert f'FullDiag solver: {label}' in text, text
                check_output(path, ref)
                exercised.setdefault(label, []).append((name, len(ref[0]), ranks if command else 1))
    require_backends(exercised)
    print('spinless FullDiag backends passed: '+', '.join(exercised), flush=True)


if __name__ == '__main__':
    main()
