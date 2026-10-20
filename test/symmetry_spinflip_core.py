#!/usr/bin/env python3
"""Px sectors: moments, CG/Lanczos/FullDiag and correlations (Sx, Sy=Sz=0, densities, products)."""
from pathlib import Path
import os
import re
import shlex
import shutil
import sys
import numpy as np
import symmetry_fulldiag_common as fd
from symmetry_kondo_common import read_manifest, read_vector_parts, run_probe, selected_layouts
from symmetry_kondo_observables import (add_correlation_requests, arbitrary_vector,
                                        check_correlation_outputs, parse_energy, parse_moments,
                                        write_probe_inputs)
from symmetry_fermion_gc_core import check_ground_state
from symmetry_spingc_common import mpi_size, run_case
from symmetry_spinflip_common import check_manifest, expert_case
from symmetry_spinflip_reference import (Case, PARITIES, flip_matrix, make_reference, raw_operator,
                                         sector_matrix)

LAYOUT = os.environ.get('HPHI_SYMMETRY_BASIS_LAYOUT', 'distributed')
REFERENCE_TOL = 1e-10
# LOBCG uses 10**(-LanczosEps/2); leave margin below the accepted residual.
CG_LANCZOS_EPS = 22
SOLVER_CASES = (Case('SpinGC', 8, 1, 1), Case('SpinGC', 8, -1, 1),
                Case('Spin', 8, 1, 1), Case('Spin', 8, -1, 1))


def reference_tolerance(expected):
    return REFERENCE_TOL * max(1., float(np.max(np.abs(expected))))


def sector_label(case):
    return f'{case.model}_L{case.sites}_k{case.momentum}_p{case.parity}'


def other_parity_energy(case):
    other = Case(case.model, case.sites, -case.parity, case.momentum, case.group, case.reflection)
    return np.linalg.eigvalsh(sector_matrix(make_reference(other)))[0]


def check_moments(root, probe):
    """<Sz>=0 for every vector of a Px sector; Sz^2, N and N^2 from the state."""
    for case in (Case('SpinGC', 5, 1, 1), Case('SpinGC', 5, -1, 1), Case('SpinGC', 4, 1, 0),
                 Case('SpinGC', 4, -1, 2), Case('SpinGC', 4, 0, 1)):
        ref = make_reference(case)
        vector = arbitrary_vector(ref.basis.shape[1], 5)
        weights = abs(ref.basis @ vector)**2
        expected = {key: float(weights @ value) for key, value in ref.moments.items()}
        if case.parity:
            assert abs(expected['Sz']) < 1e-14 and expected['Sz2'] > .1, (case, expected)
        else:
            # The Px-off control shows that the probe sees a magnetized vector.
            assert abs(expected['Sz']) > 1e-3, expected
        for layout in selected_layouts():
            path = expert_case(root/f'moments_{sector_label(case)}_{layout}', case, method=3)
            actual = parse_moments(run_probe(path, case, probe, action='moments', layout=layout,
                                             vector=vector)['text'])
            assert actual['N'] == case.sites and actual['N2'] == case.sites**2, actual
            assert actual['D'] == actual['D2'] == 0, actual
            for key in ('Sz', 'Sz2'):
                assert abs(actual[key]-expected[key]) <= 1e-12, (case, key, actual[key], expected[key])
            if case.parity:
                assert actual['Sz'] == 0, actual
    print('Px moments passed: Sz=0, Sz^2 from the representatives', flush=True)


def check_moment_failures(root, probe, launcher):
    """A NaN or corrupt storage on one rank fails the Px moments on every rank."""
    case = Case('SpinGC', 5, -1, 1)
    vector = arbitrary_vector(make_reference(case).basis.shape[1], 9)
    for name, key in (('nan', 'HPHI_TEST_SYMMETRY_NAN_RANK'),
                      ('storage', 'HPHI_TEST_SYMMETRY_INVALID_STORAGE_RANK')):
        path = expert_case(root/f'moments_failure_{name}', case, method=3)
        write_probe_inputs(path, vector, mpi_size(launcher))
        env = dict(os.environ, HPHI_TEST_SYMMETRY_ACTION='moments', HPHI_SYMMETRY_BASIS_LAYOUT=LAYOUT,
                   **{key: '0'})
        run_case(path, probe, name, launcher, env, 'Error: sector moments probe failed.')
    print('Px moments fail collectively on a bad rank', flush=True)


def check_solvers(root, probe, launcher):
    """CG reaches the residual bound; the compatibility Lanczos path the reference energy."""
    for case in SOLVER_CASES:
        ref = make_reference(case)
        h = sector_matrix(ref)
        expected = np.linalg.eigvalsh(h)[0]
        # A run of the other parity of equal momentum cannot pass on its energy.
        assert abs(expected-other_parity_energy(case)) > 1e-6, case
        for method in (3, 0):
            name = {3: 'cg', 0: 'lanczos'}[method]
            path = expert_case(root/f'{sector_label(case)}_{name}', case, method=method,
                               options={'CalcMod': {'CalcEigenVec': 1},
                                        'ModPara': {'LanczosEps': CG_LANCZOS_EPS if method == 3 else 18}})
            env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=LAYOUT, HPHI_TEST_SYMMETRY_CAPTURE='1')
            if method == 0:
                env['HPHI_TEST_SYMMETRY_ACTION'] = 'lanczos'
            run_case(path, probe, name, launcher, env)
            check_manifest(read_manifest(path), case, len(h))
            reported, text = parse_energy(path)
            vector = read_vector_parts(path, 'sector_final_sample0_step0.rank', len(h))
            # Lanczos compatibility is checked against the reference ground energy:
            # its energy-change stopping rule does not enforce the CG residual bound.
            energy = check_ground_state(vector, h, expected, path.name, require_residual=(method == 3))
            assert abs(reported-energy) <= reference_tolerance(energy)
            index = ref.words.astype(np.int64)
            psi = ref.basis @ vector
            np.testing.assert_allclose(flip_matrix(case.sites)[np.ix_(index, index)] @ psi,
                                       case.parity*psi, rtol=0, atol=1e-12)
            assert float(re.search(r'^\s*Sz\s+([^\s]+)', text, re.M).group(1)) == 0., text
    print('Px CG/Lanczos ground states passed', flush=True)


def check_parity_ground_states(root, hphi, launcher):
    """Each parity reaches its own lowest state; the Px-off sector gives the smaller one."""
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=LAYOUT)
    for model in ('SpinGC', 'Spin'):
        energies = {}
        for parity in (0, *PARITIES):
            case = Case(model, 8, parity, 0)
            expected = np.linalg.eigvalsh(sector_matrix(make_reference(case)))[0]
            path = expert_case(root/f'ground_{sector_label(case)}', case, method=3,
                               options={'ModPara': {'LanczosEps': CG_LANCZOS_EPS}})
            run_case(path, hphi, 'cg', launcher, env)
            energies[parity] = parse_energy(path)[0]
            assert abs(energies[parity]-expected) <= reference_tolerance(expected), (model, parity)
        assert abs(energies[1]-energies[-1]) > 1e-6, energies
        assert abs(energies[0]-min(energies[1], energies[-1])) <= reference_tolerance(energies[0])
        print(f'{model} k=0 Px=+1/-1 ground states {energies[1]:.10f}/{energies[-1]:.10f} passed',
              flush=True)


def correlation_requests(model):
    """Rows (site, s, site, t) of local units E_st; Spin keeps one, two and NBodyG."""
    up, down, plus, minus = [0, 1, 0, 1], [0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]
    up1, plus1, minus1, up2 = [1, 1, 1, 1], [1, 1, 1, 0], [1, 0, 1, 1], [2, 1, 2, 1]
    requests = {
        'one': [plus, minus, up, down, [3, 1, 3, 0], plus],
        'two': [plus+minus1, minus+plus1, up+up1, up+[2, 1, 2, 1], plus+up1, up+[0, 1, 0, 1],
                plus+minus1],
        'nbody': [[3]+plus+minus1+up2, [3]+up+up1+up2, [3]+plus+minus1+up2]}
    if model == 'SpinGC':
        requests['two'][-1:-1] = [plus+plus1, minus+minus1]
        requests['three'] = [plus+up1+up2, up+up1+up2, plus+plus1+up2, plus+up1+up2]
    return requests


def correlation_reference(case, ref, psi, requests):
    values = {}
    for kind, rows in requests.items():
        result = []
        for row in rows:
            data = row[1:] if kind == 'nbody' else row
            factors = [tuple(data[f:f+4]) for f in range(0, len(data), 4)]
            result.append(np.vdot(psi, raw_operator(case, ref.words, factors) @ psi))
        values[kind] = np.asarray(result)
        assert abs(values[kind][0]-values[kind][-1]) < 1e-12
    one = values['one']
    # Px: <E00> = <E11> = 1/2, <Sz> = 0 and <Sy> = 0, so <S+> = <S-> is real.
    assert abs(one[2]-.5) < 1e-12 and abs(one[3]-.5) < 1e-12, one
    assert abs(one[0]-one[1]) < 1e-12 and abs(one[0].imag) < 1e-12, one
    if case.model == 'SpinGC':
        assert abs(one[0]) > 1e-3, one        # <Sx> survives the transverse field
    else:
        assert abs(one[0]) < 1e-14, one       # S+ leaves the Sz=0 sector
    # Sz^2-type products stay finite; S+_0 S-_1 carries the exchange.
    assert abs(values['two'][0]) > 1e-3 and abs(values['two'][2]) > 1e-3, values['two']
    return values


def check_correlations(root, probe, launcher):
    for case in (Case('SpinGC', 8, -1, 1), Case('Spin', 8, 1, 1)):
        ref = make_reference(case)
        h = sector_matrix(ref)
        requests = correlation_requests(case.model)
        for aggregate in (False, True):
            form = 'aggregate' if aggregate else 'legacy'
            path = expert_case(root/f'correlation_{sector_label(case)}_{form}', case, method=3,
                               options={'CalcMod': {'CalcEigenVec': 1},
                                        'ModPara': {'LanczosEps': CG_LANCZOS_EPS}})
            add_correlation_requests(path, requests, aggregate)
            env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=LAYOUT, HPHI_TEST_SYMMETRY_CAPTURE='1')
            run_case(path, probe, 'correlation', launcher, env)
            vector = read_vector_parts(path, 'sector_final_sample0_step0.rank', len(h))
            check_ground_state(vector, h, np.linalg.eigvalsh(h)[0], path.name)
            psi = ref.basis @ vector
            check_correlation_outputs(path, aggregate, requests,
                                      correlation_reference(case, ref, psi, requests),
                                      atol=REFERENCE_TOL)
        print(f'{sector_label(case)} correlations passed', flush=True)


def check_fulldiag(root, hphi, launcher):
    """Every eigenstate of LAPACK/ScaLAPACK/ELPA FullDiag and its correlations."""
    ranks = mpi_size(launcher)
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated')
    exercised = {}
    for solver, label, backend_launcher in fd.solver_backends(launcher, ranks):
        backend_ranks = mpi_size(backend_launcher)
        for case in (Case('SpinGC', 6, -1, 1), Case('Spin', 8, 1, 0)):
            ref = make_reference(case)
            h = sector_matrix(ref)
            if solver == 3 and len(h) < fd.grid_max(backend_ranks):
                print(f'{label} skipped for D={len(h)} < process grid of np={backend_ranks}', flush=True)
                continue
            requests = correlation_requests(case.model)
            requests.pop('nbody')
            path = expert_case(root/f'fulldiag_{sector_label(case)}_solver{solver}', case, method=2,
                               options={'CalcMod': {'Solver': solver, 'NGPU': 0, 'OutputEigenVec': 1}})
            add_correlation_requests(path, requests, False)
            text = run_case(path, hphi, 'fulldiag', backend_launcher, env)
            assert f'FullDiag solver: {label}' in text, text
            check_manifest(read_manifest(path), case, len(h))
            energies = np.loadtxt(path/'output/zvo_energy_sector.dat', ndmin=2)[:, 1]
            vectors, _ = fd.read_states(path/'output', len(h), backend_ranks)
            fd.check_eigensystem(h, energies, vectors)
            for kind, rows in requests.items():
                actual = fd.read_green(path/'output', kind, len(h), False, rows)
                for state in range(len(h)):
                    psi = ref.basis @ vectors[:, state]
                    expected = [np.vdot(psi, raw_operator(case, ref.words,
                                [tuple(row[f:f+4]) for f in range(0, len(row), 4)]) @ psi)
                                for row in rows]
                    np.testing.assert_allclose(actual[state], expected, rtol=0, atol=REFERENCE_TOL,
                                               err_msg=f'{path.name} {kind} state {state}')
            exercised.setdefault(label, []).append((sector_label(case), len(h), backend_ranks))
            print(f'{sector_label(case)} FullDiag {label}: {len(h)} states and correlations passed',
                  flush=True)
    # A CI job that requires a backend must not pass on its skip message.
    fd.require_backends(exercised)


def check_raw(root, hphi):
    """The raw space (no TransSym) reaches the lowest of all (k, p) sector energies."""
    env = dict(os.environ)
    for model in ('SpinGC', 'Spin'):
        sites = 8
        lowest = min(np.linalg.eigvalsh(sector_matrix(make_reference(Case(model, sites, p, k))))[0]
                     for k in range(sites) for p in PARITIES)
        path = expert_case(root/f'raw_{model}', Case(model, sites, 0, 0), method=3,
                           options={'ModPara': {'LanczosEps': CG_LANCZOS_EPS}})
        names = (path/'sym.def').read_text().replace('TransSym group.def\n', '')
        (path/'sym.def').write_text(names)
        # Raw MPI site splitting is not part of this suite; run it serially.
        run_case(path, hphi, 'raw', [], env)
        assert not (path/'output/zvo_symmetry_sector.dat').exists()
        assert abs(parse_energy(path)[0]-lowest) <= 1e-8*max(1, abs(lowest)), (model, lowest)
    print('raw Spin/SpinGC ground states agree with the lowest Px sector', flush=True)


def main():
    hphi, probe = (Path(p).resolve() for p in sys.argv[1:3])
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    assert mpi_size(launcher) is not None
    root = Path('symmetry_spinflip_core')
    # Captured vectors of a run with another rank count must not be reread.
    if root.exists():
        shutil.rmtree(root)
    # The ELPA job selects only the dense backends with scope=fulldiag.
    scope = os.environ.get('HPHI_TEST_SPINFLIP_CORE_SCOPE', 'all')
    assert scope in ('all', 'fulldiag'), scope
    if scope == 'all':
        check_moments(root, probe)
        check_moment_failures(root, probe, launcher)
        check_solvers(root, probe, launcher)
        check_parity_ground_states(root, hphi, launcher)
        check_correlations(root, probe, launcher)
    check_fulldiag(root, hphi, launcher)
    if scope == 'all':
        check_raw(root, hphi)
    print(f'Spin-flip core ({scope}) passed')


if __name__ == '__main__':
    main()
