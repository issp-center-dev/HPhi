#!/usr/bin/env python3
"""GC fermion sectors with pair H: moments, CG/Lanczos/FullDiag, normal and anomalous correlations."""
from pathlib import Path
import os
import re
import shlex
import shutil
import sys
import numpy as np
from symmetry_kondo_common import read_manifest, read_vector_parts, run_probe, selected_layouts
from symmetry_kondo_observables import (add_correlation_requests, arbitrary_vector,
                                        check_correlation_outputs, parse_energy, parse_moments)
from symmetry_spingc_common import definition, mpi_size, run_case
from symmetry_fulldiag_common import require_backends
from symmetry_fermion_gc_common import expert_case
from symmetry_fermion_gc_reference import (Case, MODELS, PARITIES, make_reference,
                                           operator_matrix, pair_operator, pair_rows,
                                           sector_matrix, width)

LAYOUT = os.environ.get('HPHI_SYMMETRY_BASIS_LAYOUT', 'distributed')
REFERENCE_TOL = 1e-10
# LOBCG uses 10**(-LanczosEps/2); leave margin below the accepted residual.
CG_LANCZOS_EPS = 22
SITES = {'HubbardGC': 4, 'SpinlessFermionGC': 5}


def reference_tolerance(expected):
    return REFERENCE_TOL * max(1., float(np.max(np.abs(expected))))


def check_ground_state(vector, h, expected_energy, label, *, require_residual=True):
    np.testing.assert_allclose(np.linalg.norm(vector), 1., atol=1e-12, rtol=0)
    hv = h @ vector
    energy = float(np.vdot(vector, hv).real)
    assert abs(energy-expected_energy) <= reference_tolerance(expected_energy), label
    residual = np.linalg.norm(hv-energy*vector) / max(1., np.linalg.norm(h, 2))
    if require_residual:
        assert residual <= REFERENCE_TOL, (label, residual)
    residual_check = 'passed' if require_residual else 'diagnostic only'
    print(f'{label} ground-state energy passed; relative residual={residual:.3e} '
          f'({residual_check})', flush=True)
    return energy


def gc_case(root, name, case, method, *, pairs=True, options=None, parameters=None):
    rows = pair_rows(case.model, case.sites) if pairs else None
    return expert_case(root/name, case, method=method, options=options,
                       parameters=parameters, pairs=rows)


def check_moments(root, probe):
    """Arbitrary vectors mix N inside the sector; every moment comes from the state."""
    # Spinless L=5 has N and N+2 at k=1 for both parities; L=3 or 4 does not.
    for model, sites in (('HubbardGC', 3), ('SpinlessFermionGC', 5)):
        for parity in PARITIES:
            case = Case(model, sites, parity, 1)
            ref = make_reference(case)
            vector = arbitrary_vector(ref.basis.shape[1], 3)
            weights = abs(ref.basis @ vector)**2
            expected = {key: float(weights @ value) for key, value in ref.moments.items()}
            assert expected['N2']-expected['N']**2 > 1e-3, (model, parity, expected)
            for layout in selected_layouts():
                path = gc_case(root, f'moments_{model}_{parity}_{layout}', case, 3)
                actual = parse_moments(run_probe(path, case, probe, action='moments',
                                                 layout=layout, vector=vector)['text'])
                for key, value in actual.items():
                    assert abs(value-expected[key]) <= 1e-10, (model, parity, key, value, expected[key])
    print('GC moments passed', flush=True)


def check_solvers(root, probe, launcher):
    """CG and the compatibility Lanczos path reach the pair-H ground state of a complex k."""
    for model in MODELS:
        for parity, method in ((0, 3), (-1, 3), (0, 0)):
            case = Case(model, SITES[model], parity, 1)
            ref = make_reference(case, pairs=pair_rows(model, case.sites))
            h = sector_matrix(ref)
            name = {3: 'cg', 0: 'lanczos'}[method]
            path = gc_case(root, f'{model}_{parity}_{name}', case, method,
                           options={'CalcMod': {'CalcEigenVec': 1},
                                    'ModPara': {'LanczosEps': CG_LANCZOS_EPS if method == 3 else 18}})
            env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=LAYOUT, HPHI_TEST_SYMMETRY_CAPTURE='1')
            if method == 0:
                env['HPHI_TEST_SYMMETRY_ACTION'] = 'lanczos'
            run_case(path, probe, name, launcher, env)
            reported, text = parse_energy(path)
            vector = read_vector_parts(path, 'sector_final_sample0_step0.rank', len(h))
            # Lanczos compatibility is checked against the reference ground energy:
            # its energy-change stopping rule does not enforce the CG residual bound.
            energy = check_ground_state(vector, h, np.linalg.eigvalsh(h)[0], path.name,
                                        require_residual=(method == 3))
            assert abs(reported-energy) <= reference_tolerance(energy)
            weights = abs(ref.basis @ vector)**2
            for key, label in (('Sz', 'Sz'), ('D', 'Doublon')):
                actual = float(re.search(rf'^\s*{label}\s+([^\s]+)', text, re.M).group(1))
                expected = weights @ ref.moments[key]
                assert abs(actual-expected) <= reference_tolerance(expected), (model, key)


def check_parity_ground_states(root, hphi, launcher):
    """The lowest state of each parity is found directly; full GC gives the smaller one."""
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=LAYOUT)
    for model in MODELS:
        energies = {}
        for parity in PARITIES:
            case = Case(model, SITES[model], parity, 0)
            expected = np.linalg.eigvalsh(sector_matrix(
                make_reference(case, pairs=pair_rows(model, case.sites))))[0]
            path = gc_case(root, f'parity_{model}_{parity}', case, 3,
                           options={'ModPara': {'LanczosEps': CG_LANCZOS_EPS}})
            run_case(path, hphi, 'cg', launcher, env)
            energies[parity] = parse_energy(path)[0]
            assert abs(energies[parity]-expected) <= reference_tolerance(expected), (model, parity)
        assert abs(energies[1]-energies[-1]) > 1e-6, energies
        assert abs(energies[0]-min(energies[1], energies[-1])) <= reference_tolerance(energies[0])
        print(f'{model} even/odd ground states {energies[1]:.10f}/{energies[-1]:.10f} passed', flush=True)


def normal_requests(model):
    if model == 'HubbardGC':
        hop, flip, mixed, density = [1, 0, 0, 0], [0, 0, 0, 1], [1, 0, 0, 1], [2, 1, 2, 1]
        return {'one': [hop, flip, mixed, hop],
                'two': [hop + density, flip + [1, 0, 1, 1], flip + [1, 1, 1, 0], hop + density],
                'nbody': [[3] + hop + density + flip, [3] + hop + density + flip]}
    hop, back, density = [1, 0, 0, 0], [0, 0, 2, 0], [2, 0, 2, 0]
    return {'one': [hop, back, density, hop],
            'two': [hop + density, back + [1, 0, 1, 0], hop + density],
            'nbody': [[2] + hop + density, [2] + hop + density]}


def anomalous_requests(model):
    """F_ij, F_ji, the adjoint row and a repeated first row."""
    if model == 'HubbardGC':
        return [[0, 0, 0, 0, 1], [0, 0, 1, 0, 0], [1, 0, 1, 0, 0], [0, 0, 0, 1, 1],
                [0, 1, 1, 0, 0], [1, 0, 0, 1, 0], [0, 0, 0, 0, 1]]
    return [[0, 0, 0, 1, 0], [0, 1, 0, 0, 0], [1, 1, 0, 0, 0], [0, 0, 0, 2, 0],
            [1, 0, 0, 2, 0], [0, 0, 0, 1, 0]]


def normal_reference(case, ref, psi, requests):
    w = width(case.model)
    values = {}
    for kind, rows in requests.items():
        result = []
        for row in rows:
            data = row[1:] if kind == 'nbody' else row
            string = []
            for f in range(0, len(data), 4):
                i, s, j, t = data[f:f+4]
                string += [(1, w*i+s), (0, w*j+t)]
            result.append(np.vdot(psi, operator_matrix(case, ref.words, [(1, string)]) @ psi))
        values[kind] = np.asarray(result)
        assert abs(values[kind][0]-values[kind][-1]) < 1e-12
    return values


def anomalous_reference(case, ref, psi, rows):
    values = np.array([np.vdot(psi, operator_matrix(
        case, ref.words, [(1, pair_operator(case.model, *row))]) @ psi) for row in rows])
    # Antisymmetry, the adjoint relation and a nonreal value are part of the fixture.
    assert abs(values[0]) > 1e-3 and abs(values[0]+values[1]) < 1e-12
    assert abs(values[2]-np.conj(values[0])) < 1e-12 and abs(values[0].imag) > 1e-4
    assert abs(values[0]-values[-1]) < 1e-12
    return values


def read_anomalous(path, aggregate, state=None):
    name = ('zvo_AnomalousG_eigen.dat' if aggregate else 'zvo_AnomalousG.dat') \
        if state is None else f'zvo_AnomalousG_eigen{state}.dat'
    data = np.loadtxt(path/'output'/name, ndmin=2)
    prefix = 1 if aggregate else 0
    if aggregate:
        np.testing.assert_array_equal(data[:, 0], 0)
    return data[:, prefix:prefix+5].astype(int), data[:, -2]+1j*data[:, -1]


def check_correlations(root, probe, launcher):
    for model in MODELS:
        case = Case(model, SITES[model], 1, 1)
        ref = make_reference(case, pairs=pair_rows(model, case.sites))
        h = sector_matrix(ref)
        requests, rows = normal_requests(model), anomalous_requests(model)
        for aggregate in (False, True):
            form = 'aggregate' if aggregate else 'legacy'
            path = gc_case(root, f'correlation_{model}_{form}', case, 3,
                           options={'CalcMod': {'CalcEigenVec': 1},
                                    'ModPara': {'LanczosEps': CG_LANCZOS_EPS}})
            add_correlation_requests(path, requests, aggregate)
            definition(path, 'anomalousg.def', rows)
            with (path/'sym.def').open('a') as f:
                f.write('AnomalousG anomalousg.def\n')
            env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=LAYOUT, HPHI_TEST_SYMMETRY_CAPTURE='1')
            run_case(path, probe, 'correlation', launcher, env)
            vector = read_vector_parts(path, 'sector_final_sample0_step0.rank', len(h))
            check_ground_state(vector, h, np.linalg.eigvalsh(h)[0], path.name)
            psi = ref.basis @ vector
            check_correlation_outputs(path, aggregate, requests,
                                      normal_reference(case, ref, psi, requests), atol=REFERENCE_TOL)
            indices, actual = read_anomalous(path, aggregate)
            np.testing.assert_array_equal(indices, np.asarray(rows))
            np.testing.assert_allclose(actual, anomalous_reference(case, ref, psi, rows),
                                       atol=REFERENCE_TOL, rtol=0)
        # Without pair terms the nondegenerate ground state has a fixed N.
        path = gc_case(root, f'correlation_{model}_no_pairs', case, 3, pairs=False,
                       options={'ModPara': {'LanczosEps': CG_LANCZOS_EPS}})
        definition(path, 'anomalousg.def', rows)
        with (path/'sym.def').open('a') as f:
            f.write('AnomalousG anomalousg.def\n')
        run_case(path, probe, 'correlation', launcher, dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=LAYOUT))
        assert np.max(abs(read_anomalous(path, False)[1])) <= REFERENCE_TOL
        print(f'{model} normal and anomalous correlations passed', flush=True)


def fulldiag_spectrum(path):
    data = np.loadtxt(path/'output/zvo_energy_sector.dat', ndmin=2)
    np.testing.assert_array_equal(data[:, 0], np.arange(len(data)))
    assert np.isfinite(data).all()
    return data[:, 1]


def check_fulldiag(root, hphi, launcher):
    ranks = mpi_size(launcher)
    backends = [(0, 'LAPACK', [])]
    for solver, flag, label in [(1, 'HPHI_HAS_SCALAPACK', 'ScaLAPACK'), (3, 'HPHI_HAS_ELPA', 'ELPA')]:
        if os.environ.get(flag) == '1' and ranks > 1:
            backends.append((solver, label, launcher))
        else:
            print(f'{label} not exercised: backend unavailable or launcher has one rank', flush=True)
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated')
    exercised = {}
    for solver, label, backend_launcher in backends:
        for model in MODELS:
            case = Case(model, SITES[model]-1, -1, 1)
            ref = make_reference(case, pairs=pair_rows(model, case.sites))
            rows = anomalous_requests(model)
            path = gc_case(root, f'fulldiag_{model}_solver{solver}', case, 2,
                           options={'CalcMod': {'Solver': solver, 'NGPU': 0}})
            definition(path, 'anomalousg.def', rows)
            with (path/'sym.def').open('a') as f:
                f.write('AnomalousG anomalousg.def\n')
            text = run_case(path, hphi, 'fulldiag', backend_launcher, env)
            assert f'FullDiag solver: {label}' in text, text
            energies, vectors = np.linalg.eigh(sector_matrix(ref))
            error = np.max(abs(fulldiag_spectrum(path)-energies))
            assert error <= reference_tolerance(energies), (model, label, error)
            assert read_manifest(path)['model'] == model
            # Anomalous correlations of the eigenstates without degenerate partners.
            gaps = np.diff(energies)
            checked = 0
            for n in range(len(energies)):
                if (n and gaps[n-1] < 1e-8) or (n+1 < len(energies) and gaps[n] < 1e-8):
                    continue
                psi = ref.basis @ vectors[:, n]
                expected = np.array([np.vdot(psi, operator_matrix(
                    case, ref.words, [(1, pair_operator(model, *row))]) @ psi) for row in rows])
                indices, actual = read_anomalous(path, False, n)
                np.testing.assert_array_equal(indices, np.asarray(rows))
                np.testing.assert_allclose(actual, expected, atol=REFERENCE_TOL, rtol=0)
                checked += 1
            assert checked >= len(energies)//2, (model, checked)
            exercised.setdefault(label, []).append((model, len(energies), ranks))
            print(f'{model} FullDiag {label} error={error:.3e}, {checked} AnomalousG states passed',
                  flush=True)
    # A CI job that requires a backend must not pass on its skip message.
    require_backends(exercised)


def check_pair_free_limit(root, hphi):
    """Without pairs, one momentum of the GC space is the sum of its fixed-N sectors."""
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated')
    for model, fixed_model, sites in (('HubbardGC', 0, 3), ('SpinlessFermionGC', 7, 4)):
        case = Case(model, sites, 0, 0)
        path = gc_case(root, f'limit_{model}', case, 2, pairs=False)
        run_case(path, hphi, 'gc', [], env)
        # Only the vacuum has N=0; with no constant term it has energy zero.
        combined = [0.]
        ref = make_reference(case)
        # Translations keep N, so every GC basis column has one N.
        n_of_column = np.rint(ref.moments['N'] @ abs(ref.basis)**2).astype(int)
        for n in range(1, width(model)*sites+1):
            if not np.any(n_of_column == n):
                continue
            path = gc_case(root, f'limit_{model}_N{n}', case, 2, pairs=False,
                           options={'CalcMod': {'CalcModel': fixed_model}, 'ModPara': {'Ncond': n}})
            run_case(path, hphi, 'fixed_n', [], env)
            assert read_manifest(path)['model'] == ('HubbardNConserved' if fixed_model == 0
                                                    else 'SpinlessFermion')
            combined += list(fulldiag_spectrum(path))
        np.testing.assert_allclose(np.sort(combined), fulldiag_spectrum(root/f'limit_{model}'),
                                   rtol=0, atol=1e-10)
        print(f'{model} pair-free limit matches the fixed-N sectors', flush=True)


def check_raw(root, hphi):
    """Raw HubbardGC with the same pairs reaches the lowest sector energy and its F."""
    model = 'HubbardGC'
    energies = [np.linalg.eigvalsh(sector_matrix(make_reference(
        Case(model, SITES[model], 0, k), pairs=pair_rows(model, SITES[model]))))[:2]
        for k in range(SITES[model])]
    k = int(np.argmin([e[0] for e in energies]))
    lowest = sorted(e for pair in energies for e in pair)
    assert lowest[1]-lowest[0] > 1e-6, lowest[:2]
    rows = anomalous_requests(model)
    values = {}
    for label, transsym in (('raw', False), ('sector', True)):
        # An explicit ParticleParity 0 is the full space in raw input and sectors alike.
        path = gc_case(root, f'raw_{model}_{label}', Case(model, SITES[model], 0, k), 3,
                       options={'ModPara': {'LanczosEps': CG_LANCZOS_EPS, 'ParticleParity': 0}})
        definition(path, 'anomalousg.def', rows)
        names = (path/'sym.def').read_text() + 'AnomalousG anomalousg.def\n'
        if not transsym:
            names = names.replace('TransSym group.def\n', '')
        (path/'sym.def').write_text(names)
        # Raw MPI site splitting is not part of this suite; run it serially.
        run_case(path, hphi, label, [], dict(os.environ))
        assert (path/'output/zvo_symmetry_sector.dat').exists() == transsym
        assert abs(parse_energy(path)[0]-lowest[0]) <= reference_tolerance(lowest[0])
        values[label] = read_anomalous(path, False)[1]
    np.testing.assert_allclose(values['sector'], values['raw'], atol=1e-8, rtol=0)
    print('raw HubbardGC pair ground state and AnomalousG agree with the sector', flush=True)


def check_rejections(root, hphi, launcher):
    case = Case('HubbardGC', 3, 0, 1)
    path = gc_case(root, 'reject_calcspec', case, 3, options={'CalcMod': {'CalcSpec': 1}})
    run_case(path, hphi, 'reject', launcher, dict(os.environ),
             expected_error='AnomalousTerm does not support CalcSpec.')


def main():
    hphi, probe = (Path(p).resolve() for p in sys.argv[1:3])
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    assert mpi_size(launcher) is not None
    root = Path('symmetry_fermion_gc_core')
    # Captured vectors of a run with another rank count must not be reread.
    if root.exists():
        shutil.rmtree(root)
    # The ELPA job selects only the dense backends with scope=fulldiag.
    scope = os.environ.get('HPHI_TEST_FERMION_GC_CORE_SCOPE', 'all')
    assert scope in ('all', 'fulldiag'), scope
    if scope == 'all':
        check_rejections(root, hphi, launcher)
        check_moments(root, probe)
        check_solvers(root, probe, launcher)
        check_parity_ground_states(root, hphi, launcher)
        check_correlations(root, probe, launcher)
    check_fulldiag(root, hphi, launcher)
    if scope == 'all':
        check_pair_free_limit(root, hphi)
        check_raw(root, hphi)
    print(f'Fermion GC core ({scope}) passed')


if __name__ == '__main__':
    main()
