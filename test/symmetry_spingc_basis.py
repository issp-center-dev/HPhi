"""SpinGC input boundaries, independent matrices and LAPACK/ScaLAPACK/ELPA spectra."""
import os
from pathlib import Path
import re
import shlex
import struct
import sys
import tempfile

import numpy as np
import symmetry_spingc_common as c


def translation(n, k):
    return [[(i+g) % n for i in range(n)] for g in range(n)], np.exp(2j*np.pi*k*np.arange(n)/n)


def close_matrix(actual, expected, tolerance=1e-11):
    assert np.isfinite(actual).all() and np.isfinite(expected).all()
    error = np.linalg.norm(actual-expected, 2)
    assert error <= tolerance*max(1, np.linalg.norm(expected, 2)), error


def fixtures():
    n = 8
    p, chi = translation(n, 0)
    field = {'Trans': [(i, a, i, b, .5, 0) for i in range(n) for a, b in [(1, 0), (0, 1)]]}
    yield 'A', n, p, chi, field, -sum(c.spin_operators(n)['x'])
    for n, prefix, momenta in [(8, 'B', range(8)), (6, 'D', [1])]:
        h = c.mixed_hamiltonian(n)
        families = c.mixed_families(n)
        close_matrix(c.families_matrix(n, families), h, 1e-12)
        control = c.mixed_families(n, all_interall=True)
        close_matrix(c.families_matrix(n, control), h, 1e-12)
        for k in momenta:
            p, chi = translation(n, k)
            yield '{}k{}'.format(prefix, k), n, p, chi, families, h
            if k == 0:
                yield prefix+'interall', n, p, chi, control, h
    for label, rows in [('one', [(0, 1, .31)]), ('reverse', [(0, 1, .31), (1, 0, .31)]),
                        ('duplicate', [(0, 1, .31)]*2), ('onsite', [(0, 0, .31)]),
                        ('reverseonly', [(1, 0, .31)]),
                        ('onsite_added', [(0, 1, .31), (0, 0, .31)])]:
        f = {'PairLift': rows}
        yield 'C'+label, 2, [[0, 1]], [1], f, c.families_matrix(2, f)
    f = {'PairLift': [(i, (i+1) % 3, .31) for i in range(3)]}
    h = c.families_matrix(3, f)
    for k in range(3):
        p, chi = translation(3, k)
        yield 'C3k'+str(k), 3, p, chi, f, h
    for parity in [1, -1]:
        yield 'E'+str(parity), 4, [list(range(4)), [(-i) % 4 for i in range(4)]], [1, parity], c.mixed_families(4, 0), c.mixed_hamiltonian(4, 0)
    yield 'F', 2, [[0, 1]], [1], {}, np.zeros((4, 4), complex)


def reference(n, perms, chars, h):
    b, reps, projector = c.sector_basis(n, perms, chars)
    close_matrix(b.conj().T@b, np.eye(len(reps)), 1e-12)
    close_matrix(h, h.conj().T, 1e-12)
    assert np.linalg.norm(h@projector-projector@h, 2) <= 1e-12*max(1, np.linalg.norm(h, 2))
    for perm, chi in zip(perms, chars):
        moved = [sum(((r >> i) & 1) << perm[i] for i in range(n)) for r in range(1 << n)]
        transformed = np.empty_like(b)
        transformed[moved] = b
        assert np.linalg.norm(transformed-chi*b, 2) <= 1e-12
    return b.conj().T@h@b, reps


def manifest(path, dim, n, families):
    data = dict(line.split('=', 1) for line in (path/'output/zvo_symmetry_sector.dat').read_text().splitlines())
    assert data['model'] == 'SpinGC' and data['fixed_quantities'] == 'none'
    assert not any(k in data for k in ['fixed_2sz', 'fixed_ne', 'fixed_nup', 'fixed_ndown'])
    assert int(data['sector_dim']) == dim and int(data['full_dim']) == 1 << n
    assert data['hamiltonian_digest'].startswith('hphi-parsed-hamiltonian-fnv1a64-v3:')
    assert int(data['pair_lift']) == len(families.get('PairLift', []))
    return data


def probe_matrix(path, dim, reps, ranks):
    files = sorted(path.glob('spingc_probe_rank_*.dat'))
    assert len(files) == ranks
    matrix = np.empty((dim, dim), complex)
    seen = np.zeros((dim, dim), bool)
    coverage = np.zeros(dim, bool)
    empty = 0
    for file in files:
        info = dict(line.split('=', 1) for line in file.with_suffix('.info').read_text().splitlines())
        assert int(info['dim']) == dim
        offset, size = int(info['offset']), int(info['local_dim'])
        assert not coverage[offset:offset+size].any()
        coverage[offset:offset+size] = True
        assert int(info['raw_basis_elements']) == int(info['raw_diagonal_elements']) == 0
        assert int(info['initial_vector_elements']) == 3*(size+1)
        assert int(info['probe_vector_elements']) == 2*(size+1)
        lines = file.read_text().splitlines()
        assert len(lines) == dim*size
        empty += size == 0
        for line in lines:
            row, col, rep, real, imag = line.split()
            row, col = int(row)-1, int(col)-1
            assert offset <= row < offset+size and 0 <= col < dim and not seen[row, col]
            assert int(rep) == reps[row]
            matrix[row, col] = complex(float(real), float(imag))
            seen[row, col] = True
    assert coverage.all() and seen.all()
    assert empty == max(0, ranks-dim)
    return matrix


def positive(root, hphi, probe, mode):
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = c.mpi_size(launcher)
    assert ranks is not None, 'Specify MPI process count explicitly for probe rank checks'
    dimensions, spectra, identities = {}, [], {}
    for name, n, perms, chars, families, h in fixtures():
        expected, reps = reference(n, perms, chars, h)
        dim = len(reps)
        dimensions[name] = dim
        if name == 'Dk1':
            assert dim == 9
        if mode == 'basis':
            for layout in ['replicated', 'distributed']:
                path = root/(name+'_'+layout)
                c.write_case(path, n, perms, chars, families, 3, {})
                env = dict(os.environ, HPHI_TEST_SPINGC_ACTION='matvec', HPHI_SYMMETRY_BASIS_LAYOUT=layout)
                c.run_case(path, probe, 'probe', launcher, env)
                actual = probe_matrix(path, dim, reps, ranks)
                close_matrix(actual, expected)
                if name == 'F':
                    assert np.count_nonzero(actual) == 0
                identity = manifest(path, dim, n, families)
                keys = ('group_digest', 'sector_digest', 'hamiltonian_digest')
                digests = tuple(identity[key] for key in keys)
                if name in identities:
                    assert identities[name] == digests
                identities[name] = digests
        else:
            path = root/name
            c.write_case(path, n, perms, chars, families, 2, {'CalcMod': {'Solver': 0}})
            env = dict(os.environ)
            env.pop('HPHI_SYMMETRY_BASIS_LAYOUT', None)
            text = c.run_case(path, hphi, 'fulldiag', [], env)
            values = np.loadtxt(path/'output/zvo_energy_sector.dat', ndmin=2)
            np.testing.assert_array_equal(values[:, 0], np.arange(dim))
            assert np.isfinite(values).all()
            np.testing.assert_allclose(values[:, 1], np.linalg.eigvalsh(expected), atol=3e-8, rtol=0)
            if name == 'F':
                assert np.count_nonzero(values[:, 1]) == 0
            assert 'raw_basis_elements=0' in text
            manifest(path, dim, n, families)
            if re.fullmatch('Bk[0-7]', name):
                spectra.extend(values[:, 1])
            if name == 'A':
                assert abs(values[0, 1]+4) < 3e-8
        print('{} {} dim={} passed'.format(mode, name, dim), flush=True)
    assert sum(dimensions['Bk'+str(k)] for k in range(8)) == 256
    assert dimensions['E1']+dimensions['E-1'] == 16
    if mode == 'basis':
        assert len({identities[name][2] for name in
                    ['Cone', 'Creverse', 'Cduplicate', 'Consite', 'Creverseonly', 'Consite_added']}) == 6
        # A failure opening a single rank file must stop all participants,
        # including ranks with no matrix rows (D at np16).
        p, chi = translation(6, 1)
        for layout in ['replicated', 'distributed']:
            path = root/('probe_open_failure_'+layout)
            c.write_case(path, 6, p, chi, c.mixed_families(6), 3, {})
            (path/('spingc_probe_rank_{}.dat'.format(ranks-1))).mkdir()
            c.run_case(path, probe, 'reject', launcher,
                       dict(os.environ, HPHI_TEST_SPINGC_ACTION='matvec',
                            HPHI_SYMMETRY_BASIS_LAYOUT=layout),
                       # Non-MPI exitMPI has no MPI abort banner. Require
                       # completed setup plus nonzero exit in that build.
                       'MPI-ERROR MESSAGE' if ranks > 1 else
                       ('Symmetry matvec:' if layout == 'replicated' else
                        'Symmetry distributed matvec:'))
    if mode == 'fulldiag':
        np.testing.assert_allclose(sorted(spectra), np.linalg.eigvalsh(c.mixed_hamiltonian(8)), atol=3e-8, rtol=0)
        backend_full_diag(root, hphi, launcher, ranks)


def backend_full_diag(root, hphi, launcher, ranks):
    """Only B k0/1 need distributed solvers; LAPACK covers the full spectrum."""
    env = dict(os.environ)
    env.pop('HPHI_SYMMETRY_BASIS_LAYOUT', None)
    backends = []
    for solver, flag, label in [(1, 'HPHI_HAS_SCALAPACK', 'ScaLAPACK'),
                                 (3, 'HPHI_HAS_ELPA', 'ELPA')]:
        if os.environ.get(flag) == '1' and ranks > 1:
            backends.append((solver, label))
        else:
            print('{} not exercised: backend unavailable or launcher has one rank'.format(label), flush=True)
    for solver, label in backends:
        for k in (0, 1):
            p, chi = translation(8, k)
            families = c.mixed_families(8)
            expected, reps = reference(8, p, chi, c.mixed_hamiltonian(8))
            path = root/'Bk{}_solver{}'.format(k, solver)
            c.write_case(path, 8, p, chi, families, 2,
                         {'CalcMod': {'Solver': solver, 'NGPU': 0}})
            text = c.run_case(path, hphi, 'fulldiag', launcher, env)
            values = np.loadtxt(path/'output/zvo_energy_sector.dat', ndmin=2)
            assert np.isfinite(values).all()
            np.testing.assert_array_equal(values[:, 0], np.arange(len(reps)))
            np.testing.assert_allclose(values[:, 1], np.linalg.eigvalsh(expected), atol=3e-8, rtol=0)
            data = manifest(path, len(reps), 8, families)
            assert int(data['solver_id']) == solver and data['output_scope'] == 'eigenvalues'
            assert data['basis_layout'] == 'replicated'
            storage = 'replicated' if solver == 1 else 'column_panel'
            assert 'matrix={} eigenvectors=distributed'.format(storage) in text
            assert not list((path/'output').glob('*eigenvec*'))
            print('{} Solver{} B k{} np{} dim={} passed'.format(label, solver, k, ranks, len(reps)), flush=True)
        # Nonserial ExpecMode must be rejected on actual multiple ranks;
        # one-rank readdef intentionally demotes it before the sector gate.
        path = root/'expecmode_solver{}'.format(solver)
        c.write_case(path, 8, p, chi, families, 2,
                     {'CalcMod': {'Solver': solver, 'ExpecMode': 1}})
        c.run_case(path, hphi, 'reject', launcher, env, 'MAGMA and ExpecMode are not supported')
    if any(solver == 3 for solver, _ in backends):
        # Two-site odd translation sector has dimension one: too small
        # for any multiple-rank ELPA process grid, including required np2.
        p, chi = translation(2, 1)
        path = root/'elpa_small_sector'
        c.write_case(path, 2, p, chi, {}, 2, {'CalcMod': {'Solver': 3}})
        _, reps, _ = c.sector_basis(2, p, chi)
        assert len(reps) == 1
        c.run_case(path, hphi, 'reject', launcher, env, 'smaller than the process grid')
        assert not (path/'output/zvo_energy_sector.dat').exists()
        print('ELPA Solver3 process-grid rejection np{} dim=1 passed'.format(ranks), flush=True)


def validation(root, hphi):
    p, chi = translation(4, 0)
    base = root/'valid_input'
    c.write_case(base, 4, p, chi, c.mixed_families(4), 3, {})
    def negative(name, error, mutate, template=base):
        import shutil
        path = root/name
        shutil.copytree(template, path)
        mutate(path)
        c.run_case(path, hphi, 'reject', [], dict(os.environ), error)
        assert not list((path/'output').glob('*energy*'))
        assert not list((path/'output').glob('*cisajs*'))
        print('validation {} passed'.format(name), flush=True)
    def append(path, file, text):
        with (path/file).open('a') as stream:
            stream.write(text)
    fixed = [{'2Sz': 0}, {'2Sz': 2}, {'Ncond': 0}, {'Ncond': 2},
             {'Nup': 0, 'Ndown': 0}, {'Nup': 1, 'Ndown': 0}, {'Nup': 0, 'Ndown': 1}]
    for i, options in enumerate(fixed):
        negative('fixed'+str(i), 'SpinGC TransSym does not accept explicit 2Sz/Nup/Ndown/Ncond.',
                 lambda path, options=options: append(path, 'mod.def', ''.join('{} {}\n'.format(k, v) for k, v in options.items())))
    for key in ['Nup', 'Ndown']:
        negative('unpaired_'+key, 'must be specified together', lambda path, key=key: append(path, 'mod.def', key+' 0\n'))
    def replace(path, file, old, new):
        (path/file).write_text((path/file).read_text().replace(old, new))
    empty = root/'valid_empty'
    c.write_case(empty, 4, p, chi, {}, 3, {})
    negative('general_spin', 'supports only Spin-1/2',
             lambda path: replace(path, 'loc.def', '0 1\n', '0 2\n'), empty)
    def extra(path, key, rows):
        c.definition(path, key+'.def', rows)
        append(path, 'sym.def', '{} {}.def\n'.format(key, key))
    negative('onebody_offsite', 'SpinGC TransSym OneBodyG requires onsite operators.',
             lambda path: extra(path, 'OneBodyG', [(0, 1, 1, 0)]))
    negative('onebody_offsite_second', 'SpinGC TransSym OneBodyG requires onsite operators.',
             lambda path: extra(path, 'OneBodyG', [(0, 1, 0, 0), (0, 1, 1, 0)]))
    negative('laser', 'TEOneBody/TETwoBody', lambda path: extra(path, 'Laser', [('amplitude', .1)]))
    def boost(path):
        (path/'boost.def').write_text('0 0 0\n0\n0 0 0 0\n')
        append(path, 'sym.def', 'Boost boost.def\n')
    negative('boost', 'TransSym does not support Boost', boost)
    negative('unsupported_family', 'unsupported term families', lambda path: extra(path, 'CoulombIntra', [(0, .2)]))
    negative('nbody_family', 'unsupported term families',
             lambda path: extra(path, 'NBodyInterAll', [(3, 0, 1, 0, 1, 1, 0, 1, 0, 2, 1, 2, 1, .3, 0)]))
    negative('pairhopping_family', 'PairHop is not active in Spin and SpinGC.',
             lambda path: extra(path, 'PairHop', [(0, 1, .2)]))
    negative('unsupported_method', 'CalcType:',
             lambda path: replace(path, 'calc.def', 'CalcType 3', 'CalcType 6'))
    negative('unsupported_model', 'supports only SpinGC, Spin',
             lambda path: replace(path, 'calc.def', 'CalcModel 4', 'CalcModel 5'), empty)
    negative('restart', 'does not support ReStart', lambda path: append(path, 'calc.def', 'ReStart 1\n'))
    negative('spectrum', 'spectrum calculations', lambda path: append(path, 'calc.def', 'CalcSpec 1\n'))
    negative('fulldiag_vectors', 'does not support OutputEigenVec',
             lambda path: (replace(path, 'calc.def', 'CalcType 3', 'CalcType 2'), append(path, 'calc.def', 'OutputEigenVec 1\n')))
    full = root/'valid_fulldiag'
    c.write_case(full, 4, p, chi, c.mixed_families(4), 2,
                 {'CalcMod': {'Solver': 0}})
    negative('fulldiag_correlations', 'sector FullDiag outputs eigenvalues only',
             lambda path: extra(path, 'OneBodyG', [(0, 1, 0, 0)]), full)
    negative('fulldiag_magma', 'MAGMA',
             lambda path: replace(path, 'calc.def', 'Solver 0', 'Solver 2'), full)
    negative('fulldiag_expecmode', 'ExpecMode',
             lambda path: append(path, 'calc.def', 'ExpecMode 1\n'), full)
    c.run_case(full, hphi, 'reject_distributed', [],
               dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='distributed'),
               'distributed symmetry basis is supported for TransSym Lanczos, TPQ, CG, TimeEvolution and cTPQ runs only')
    negative('invalid_character', 'character', lambda path: replace(path, 'group.def', '1 1.0 0.0', '1 0.5 0.0'))
    negative('nonfinite_character', 'character must be finite',
             lambda path: replace(path, 'group.def', '1 1.0 0.0', '1 nan 0.0'))
    negative('invalid_site', 'TransSym', lambda path: replace(path, 'group.def', '0 0 0 1', '0 0 4 1'))
    negative('noninvariant', 'invariance', lambda path: replace(path, 'PairLift.def', '0 1 0.17', '0 1 0.18'))
    for n in [0, 8*struct.calcsize("@L")]:
        def size(path, n=n):
            # Use empty H and identity: only the site-count boundary is invalid.
            c.write_case(path, n, [list(range(n))], [1], {}, 3, {})
        negative('nsite'+str(n), 'Nsite' if n == 0 else 'Nsite must satisfy', size)


def main():
    mode, hphi, probe = sys.argv[1:]
    hphi, probe = Path(hphi).resolve(), Path(probe).resolve()
    root = Path(tempfile.mkdtemp(prefix='symmetry_spingc_'+mode+'_', dir='.'))
    print('artifacts: {}'.format(root.resolve()), flush=True)
    if mode == 'validation':
        validation(root, hphi)
    elif mode in ['basis', 'fulldiag']:
        positive(root, hphi, probe, mode)
    else:
        raise ValueError(mode)


if __name__ == '__main__':
    main()
