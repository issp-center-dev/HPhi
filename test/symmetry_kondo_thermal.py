#!/usr/bin/env python3
"""Finite-sample Kondo TPQ/cTPQ trajectories from captured initial states."""
import os
from pathlib import Path
import shlex
import shutil
import sys
import tempfile

import numpy as np

from symmetry_kondo_common import (expert_case, read_manifest, read_vector_parts)
from symmetry_kondo_observables import (CORRELATION_FILES,
                                        add_correlation_requests,
                                        correlation_requests)
from symmetry_kondo_reference import (Case, apply_operator, make_reference,
                                      sector_hamiltonian)
from symmetry_spingc_common import mpi_size, run_case


ATOL_DATA = 2e-11
ATOL_CORRELATION = 1e-8
NUM_AVE = 2
NSTEPS = 5
UNIFORM_BETAS = np.array([0., .02, .04, .06, .08])
EXPLICIT_BETAS = np.array([0., .02, .05, .05, .08])
EXPLICIT_FLAGS = (0, 1, 0, 1, 0)


def captured_initials(path: Path, dimension: int, representatives: np.ndarray,
                      ranks: int) -> tuple[list[np.ndarray], list[float]]:
    vectors, prenorms = [], []
    for sample in range(NUM_AVE):
        stem = f'sector_initial_sample{sample}_step0.rank'
        files = sorted(path.glob(stem + '*.dat'))
        assert len(files) == ranks, (stem, files)
        infos = [read_manifest(file.with_suffix('.info')) for file in files]
        assert len({info['prenorm'] for info in infos}) == 1, infos
        assert all(int(info['dim']) == dimension for info in infos)
        assert sum(int(info['local_dim']) for info in infos) == dimension
        assert all(int(info['probe_vector_elements']) == 0 for info in infos)
        vector = read_vector_parts(path, stem, dimension)
        actual_representatives = np.zeros(dimension, dtype=np.uint64)
        for file in files:
            for line in file.read_text().splitlines():
                index, representative, *_ = line.split()
                actual_representatives[int(index) - 1] = int(representative)
        np.testing.assert_array_equal(actual_representatives, representatives)
        np.testing.assert_allclose(np.vdot(vector, vector), 1., atol=1e-12, rtol=0)
        prenorm = float(infos[0]['prenorm'])
        assert np.isfinite(prenorm) and prenorm > 0
        vectors.append(vector)
        prenorms.append(prenorm)
    return vectors, prenorms


def tpq_step(hamiltonian: np.ndarray, vector: np.ndarray,
             nsite: int) -> tuple[np.ndarray, float]:
    work = 4 * vector - hamiltonian @ vector / nsite
    norm = np.linalg.norm(work)
    assert np.isfinite(norm) and norm > 0
    return work / norm, norm


def ctpq_step(hamiltonian: np.ndarray, vector: np.ndarray,
              delta_beta: float) -> tuple[np.ndarray, float]:
    work = vector.copy()
    term = vector.copy()
    for order in range(1, 13):
        term = (-delta_beta / (2 * order)) * (hamiltonian @ term)
        work += term
    norm = np.linalg.norm(work)
    assert np.isfinite(norm) and norm > 0
    return work / norm, norm


def expected_history(kind: str, hamiltonian: np.ndarray, initial: np.ndarray,
                     prenorm: float, case: Case, reference,
                     betas: np.ndarray | None) -> dict:
    vector, norm = initial.copy(), prenorm
    result = {'SS': [], 'Norm': [], 'Flct': [], 'vectors': []}
    for step in range(NSTEPS):
        result['vectors'].append(vector.copy())
        psi = reference.basis @ vector
        weights = np.abs(psi) ** 2
        moments = {key: float(weights @ reference.moments[key])
                   for key in ('N', 'N2', 'D', 'D2', 'Sz', 'Sz2')}
        hv = hamiltonian @ vector
        energy = float(np.vdot(vector, hv).real)
        energy2 = float(np.vdot(hv, hv).real)
        beta = (2 * step / (2 * case.cells * 4 - energy)
                if kind == 'tpq' else float(betas[step]))
        result['SS'].append([beta, energy, energy2,
                             moments['D'], moments['N'], step])
        result['Norm'].append([beta, norm, prenorm, step])
        if kind != 'tpq' or step:
            result['Flct'].append([beta, moments['N'], moments['N2'],
                                   moments['D'], moments['D2'], moments['Sz'],
                                   moments['Sz2'], step])
        if step + 1 < NSTEPS:
            if kind == 'tpq':
                vector, norm = tpq_step(hamiltonian, vector, 2 * case.cells)
            else:
                vector, norm = ctpq_step(
                    hamiltonian, vector, float(betas[step + 1] - betas[step]))
                if betas[step + 1] == betas[step]:
                    np.testing.assert_allclose(vector, result['vectors'][-1],
                                               atol=5e-15, rtol=0)
                    np.testing.assert_allclose(norm, 1., atol=5e-15, rtol=0)
    return {key: np.asarray(value) if key != 'vectors' else value
            for key, value in result.items()}


def finite_data(path: Path) -> np.ndarray:
    data = np.loadtxt(path, ndmin=2)
    assert np.isfinite(data).all(), path
    return data


def check_data(path: Path, aggregate: bool, expected: list[dict]) -> None:
    for family in ('SS', 'Norm', 'Flct'):
        if aggregate:
            data = finite_data(path / 'output' / f'zvo_{family}_tpq.dat')
            prefix = [(sample, int(row[-1])) for sample in range(NUM_AVE)
                      for row in expected[sample][family]]
            np.testing.assert_array_equal(data[:, :2], prefix)
            want = np.vstack([item[family][:, :-1] for item in expected])
            np.testing.assert_allclose(data[:, 2:], want, atol=ATOL_DATA,
                                       rtol=ATOL_DATA, err_msg=family)
            assert not list((path / 'output').glob(f'zvo_{family}_rand*.dat'))
        else:
            assert not (path / 'output' / f'zvo_{family}_tpq.dat').exists()
            for sample in range(NUM_AVE):
                data = finite_data(path / 'output' / f'zvo_{family}_rand{sample}.dat')
                np.testing.assert_array_equal(data[:, -1], expected[sample][family][:, -1])
                np.testing.assert_allclose(data[:, :-1], expected[sample][family][:, :-1],
                                           atol=ATOL_DATA, rtol=ATOL_DATA,
                                           err_msg=f'{family} sample{sample}')


def check_correlation_block(data: np.ndarray, rows: list[list[int]],
                            expected: np.ndarray, prefix: int, label: str) -> None:
    assert data.shape == (len(rows), prefix + len(rows[0]) + 2), (label, data.shape)
    np.testing.assert_array_equal(data[:, prefix:-2], rows)
    actual = data[:, -2] + 1j * data[:, -1]
    assert np.isfinite(actual).all() and np.isfinite(expected).all()
    np.testing.assert_allclose(actual, expected, atol=ATOL_CORRELATION, rtol=0,
                               err_msg=label)


def thermal_correlation_reference(case: Case, reference, vector: np.ndarray,
                                  requests: dict[str, list[list[int]]]) -> dict:
    psi = reference.basis @ vector
    values = {}
    for kind, rows in requests.items():
        result = []
        for row in rows:
            data = row[1:] if kind == 'nbody' else row
            factors = tuple(tuple(data[i:i + 4]) for i in range(0, len(data), 4))
            result.append(np.vdot(psi, apply_operator(case, factors, psi)))
        result = np.asarray(result)
        assert np.isfinite(result).all() and np.max(abs(result)) > 1e-10, (kind, result)
        values[kind] = result
    return values


def check_correlations(path: Path, aggregate: bool, emitted_steps: list[int],
                       expected: list[dict], case: Case, reference,
                       requests: dict[str, list[list[int]]]) -> None:
    expected_values = [[thermal_correlation_reference(case, reference, vector, requests)
                        for vector in item['vectors']] for item in expected]
    for kind, rows in requests.items():
        stem = CORRELATION_FILES[kind]
        if aggregate:
            data = finite_data(path / 'output' / f'zvo_{stem}_tpq.dat')
            prefixes = [(sample, step) for sample in range(NUM_AVE)
                        for step in emitted_steps for _ in rows]
            np.testing.assert_array_equal(data[:, :2], prefixes)
            offset = 0
            for sample in range(NUM_AVE):
                for step in emitted_steps:
                    block = data[offset:offset + len(rows)]
                    check_correlation_block(block, rows,
                                            expected_values[sample][step][kind], 2,
                                            f'{kind} sample{sample} step{step}')
                    offset += len(rows)
        else:
            for sample in range(NUM_AVE):
                for step in range(NSTEPS):
                    file = path / 'output' / f'zvo_{stem}_set{sample}step{step}.dat'
                    if step not in emitted_steps:
                        assert not file.exists(), file
                    else:
                        check_correlation_block(
                            finite_data(file), rows, expected_values[sample][step][kind], 0,
                            f'{kind} sample{sample} step{step}')


def check_metadata(path: Path, case: Case, kind: str, ranks: int,
                   explicit: bool, dimension: int) -> None:
    values = read_manifest(path)
    assert values['ensemble'] == 'single_symmetry_sector'
    assert values['calc_type'] == ('TPQ' if kind == 'tpq' else 'cTPQ')
    assert int(values['sector_dim']) == dimension
    assert int(values['num_ave']) == NUM_AVE
    assert int(values['initial_vec_type']) == (1 if case.momentum == 0 else 0)
    assert int(values['mpi_ranks']) == ranks
    if kind == 'tpq':
        assert float(values['large_value']) == 4
    else:
        assert int(values['canonical_tpq_steps']) == NSTEPS
        assert values['beta_schedule'] == ('explicit' if explicit else 'uniform')
        if explicit:
            for step, beta in enumerate(EXPLICIT_BETAS):
                np.testing.assert_allclose(
                    list(map(float, values[f'invtemp_row_{step}'].split())),
                    [beta, 12], atol=0, rtol=0)
        else:
            assert float(values['large_value']) == 50
            assert float(values['beta_step']) == .02
            assert int(values['expand_coef']) == 12


def run_thermal(root: Path, executable: Path, launcher: list[str], case: Case,
                kind: str, explicit: bool, aggregate: bool) -> Path:
    reference = make_reference(case)
    hamiltonian = sector_hamiltonian(reference)
    schedule = 'explicit' if explicit else 'uniform'
    path = root / f'{case.model}_k{case.momentum}_{kind}_{schedule}'
    options = {
        'CalcMod': {'InitialVecType': 1 if case.momentum == 0 else 0,
                    'OutputGreenFormat': int(aggregate)},
        'ModPara': {'Lanczos_max': NSTEPS, 'initial_iv': 7,
                    'NumAve': NUM_AVE, 'ExpecInterval': 2,
                    'LargeValue': 4 if kind == 'tpq' else 50,
                    **({'ExpandCoef': 12} if kind == 'ctpq' else {})},
    }
    expert_case(path, case, method=1 if kind == 'tpq' else 5, options=options)
    requests = correlation_requests()
    add_correlation_requests(path, requests, aggregate)
    if explicit:
        (path / 'beta.def').write_text(''.join(
            f'{beta} 12 {flag} 0\n'
            for beta, flag in zip(EXPLICIT_BETAS, EXPLICIT_FLAGS)))
        with (path / 'sym.def').open('a') as handle:
            handle.write('InvTemp beta.def\n')
        betas, emitted_steps = EXPLICIT_BETAS, [0, 1, 3]
    elif kind == 'ctpq':
        betas, emitted_steps = UNIFORM_BETAS, [0, 2, 4]
    else:
        betas, emitted_steps = None, [0, 2, 4]
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='distributed',
               HPHI_TEST_SYMMETRY_CAPTURE='1')
    text = run_case(path, executable, 'run', launcher, env)
    assert 'Symmetry basis layout: distributed' in text
    vectors, prenorms = captured_initials(
        path, len(reference.representatives), reference.representatives,
        mpi_size(launcher))
    expected = [expected_history(kind, hamiltonian, vector, prenorm,
                                 case, reference, betas)
                for vector, prenorm in zip(vectors, prenorms)]
    check_metadata(path, case, kind, mpi_size(launcher), explicit, len(hamiltonian))
    check_data(path, aggregate, expected)
    check_correlations(path, aggregate, emitted_steps, expected,
                       case, reference, requests)
    print(f'{case.model} k={case.momentum} {kind} {schedule} dim={len(hamiltonian)} '
          f'prenorms={prenorms} passed', flush=True)
    return path


def reject_invtemp_eigen(root: Path, source: Path, executable: Path) -> None:
    path = root / 'reject_invtemp_eigen'
    shutil.copytree(source, path)
    rows = (path / 'beta.def').read_text().splitlines()
    fields = rows[3].split()
    fields[-1] = '1'
    rows[3] = ' '.join(fields)
    (path / 'beta.def').write_text('\n'.join(rows) + '\n')
    shutil.rmtree(path / 'output')
    run_case(path, executable, 'reject-eigen', [],
             dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='distributed'),
             'does not support InvTemp eigenvector output')
    assert not list((path / 'output').glob('*eigen*'))


def main() -> None:
    executable = Path(sys.argv[1]).resolve()
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = mpi_size(launcher)
    assert ranks in (1, 2, 4, 16), 'Kondo thermal suite requires np1, np2, np4, or np16'
    root = Path(tempfile.mkdtemp(prefix='symmetry_kondo_thermal_', dir='.'))
    print(f'artifacts: {root.resolve()}', flush=True)
    rejection_source = None
    for model in ('Kondo', 'KondoNConserved', 'KondoGC'):
        for momentum in (0, 1):
            case = Case(model, 3, 'block', None if model == 'KondoGC' else 2,
                        1 if model == 'Kondo' else None, momentum)
            aggregate = momentum == 1
            run_thermal(root, executable, launcher, case, 'tpq', False, aggregate)
            run_thermal(root, executable, launcher, case, 'ctpq', False, aggregate)
            explicit_path = run_thermal(root, executable, launcher, case,
                                        'ctpq', True, aggregate)
            if ranks == 1 and model == 'Kondo' and momentum == 0:
                rejection_source = explicit_path
    if ranks == 1:
        assert rejection_source is not None
        reject_invtemp_eigen(root, rejection_source, executable)


if __name__ == '__main__':
    main()
