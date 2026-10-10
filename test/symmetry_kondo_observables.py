#!/usr/bin/env python3
"""Kondo symmetry-sector moments, iterative solvers, and correlations."""
from pathlib import Path
import os
import re
import shlex
import sys
import tempfile

import numpy as np

from symmetry_kondo_common import expert_case, read_vector_parts, run_probe
from symmetry_kondo_reference import Case, apply_operator, make_reference, sector_hamiltonian
from symmetry_spingc_common import definition, mpi_size, run_case


MOMENT_KEYS = ('N', 'N2', 'D', 'D2', 'Sz', 'Sz2')
CORRELATION_FILES = {
    'one': 'cisajs', 'two': 'cisajscktalt', 'three': 'ThreeBody',
    'four': 'FourBody', 'six': 'SixBody', 'nbody': 'NBodyG',
}
CORRELATION_KEYWORDS = {
    'one': 'OneBodyG', 'two': 'TwoBodyG', 'three': 'ThreeBodyG',
    'four': 'FourBodyG', 'six': 'SixBodyG', 'nbody': 'NBodyG',
}


def arbitrary_vector(dimension: int, seed: int) -> np.ndarray:
    index = np.arange(1, dimension + 1, dtype=float)
    vector = np.cos(index * (.17 + seed * .01)) + 1j * np.sin(index * (.29 + seed * .02))
    return vector / np.linalg.norm(vector)


def write_probe_inputs(path: Path, vector: np.ndarray, ranks: int) -> None:
    for rank in range(ranks):
        base, remainder = divmod(len(vector), ranks)
        start = base * rank + min(rank, remainder)
        stop = start + base + (rank < remainder)
        np.savetxt(path / f'probe-input.rank{rank}.dat',
                   np.column_stack((vector[start:stop].real, vector[start:stop].imag)),
                   fmt='%.17g')


def parse_moments(text: str) -> dict[str, float]:
    result = {}
    for key in MOMENT_KEYS:
        match = re.search(rf'^SectorProbe {key} ([^\n]+)$', text, re.M)
        assert match, (key, text)
        result[key] = float(match.group(1))
    assert all(np.isfinite(value) for value in result.values())
    return result


def check_arbitrary_moments(root: Path, probe: Path, layout: str) -> None:
    for momentum in (0, 1):
        case = Case('KondoGC', 3, 'block', None, None, momentum)
        ref = make_reference(case)
        vector = arbitrary_vector(ref.basis.shape[1], 7 + momentum)
        path = expert_case(root / f'moments_k{momentum}_{layout}', case,
                           method=3, options={})
        actual = parse_moments(run_probe(path, case, probe, action='moments',
                                         layout=layout, vector=vector)['text'])
        psi = ref.basis @ vector
        weights = np.abs(psi) ** 2
        expected = {name: float(weights @ values) for name, values in ref.moments.items()}
        assert expected['Ncond2'] - expected['Ncond'] ** 2 > 1e-6
        assert abs(expected['Ncond'] - (expected['N'] - case.cells)) < 1e-12
        assert abs(expected['Ncond2'] -
                   (expected['N2'] - 2 * case.cells * expected['N'] + case.cells ** 2)) < 1e-12
        for name in MOMENT_KEYS:
            assert abs(actual[name] - expected[name]) <= 1e-8, (name, actual[name], expected[name])


def parse_energy(path: Path) -> tuple[float, str]:
    text = (path / 'output/zvo_energy.dat').read_text()
    match = re.search(r'^\s*Energy\s+([^\s]+)', text, re.M)
    assert match, text
    return float(match.group(1)), text


def check_eigenvector(vector: np.ndarray, hamiltonian: np.ndarray,
                      expected_energy: float) -> float:
    energy = float(np.vdot(vector, hamiltonian @ vector).real)
    scale = max(1.0, np.linalg.norm(hamiltonian, 2) * np.linalg.norm(vector))
    assert abs(energy - expected_energy) <= 3e-8
    residual = np.linalg.norm(hamiltonian @ vector - energy * vector) / scale
    assert residual <= 1e-8, residual
    return energy


def solver_case(root: Path, executable: Path, launcher: list[str], case: Case,
                method: int) -> None:
    ref = make_reference(case)
    hamiltonian = sector_hamiltonian(ref)
    name = 'cg' if method == 3 else 'lanczos'
    path = expert_case(root / f'{case.model}_P{case.cells}_{case.layout}_k{case.momentum}_{name}',
                       case, method=method,
                       options={'CalcMod': {'CalcEigenVec': 1},
                                'ModPara': {'LanczosEps': 18}})
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='distributed',
               HPHI_TEST_SYMMETRY_CAPTURE='1')
    if method == 0:
        env['HPHI_TEST_SYMMETRY_ACTION'] = 'lanczos'
    run_case(path, executable, name, launcher, env)
    reported, energy_text = parse_energy(path)
    vector = read_vector_parts(path, 'sector_final_sample0_step0.rank', len(hamiltonian))
    energy = check_eigenvector(vector, hamiltonian, np.linalg.eigvalsh(hamiltonian)[0])
    assert abs(reported - energy) <= 3e-8
    psi = ref.basis @ vector
    expected_sz = float(abs(psi) ** 2 @ ref.moments['Sz'])
    actual_sz = float(re.search(r'^\s*Sz\s+([^\s]+)', energy_text, re.M).group(1))
    assert abs(actual_sz - expected_sz) <= 1e-8


def check_solvers(root: Path, executable: Path, launcher: list[str]) -> None:
    for model in ('Kondo', 'KondoNConserved', 'KondoGC'):
        for momentum in (0, 1):
            case = Case(model, 3, 'block', None if model == 'KondoGC' else 2,
                        1 if model == 'Kondo' else None, momentum)
            solver_case(root, executable, launcher, case, 3)
            solver_case(root, executable, launcher, case, 0)
        for layout in ('block', 'alternating'):
            case = Case(model, 4, layout, None if model == 'KondoGC' else 2,
                        0 if model == 'Kondo' else None, 1)
            solver_case(root, executable, launcher, case, 3)


def correlation_requests() -> dict[str, list[list[int]]]:
    hop = [3, 0, 4, 0]
    local_flip = [0, 0, 0, 1]
    flip_pair = local_flip + [3, 1, 3, 0]
    local_up = [1, 0, 1, 0]
    cond_up = [5, 0, 5, 0]
    plus_minus = [0, 0, 0, 1, 0, 1, 0, 0]
    minus_plus = [0, 1, 0, 0, 0, 0, 0, 1]
    requests = {
        'one': [hop, local_flip, hop],
        'two': [hop + local_up, flip_pair, plus_minus, minus_plus, hop + local_up],
        'three': [hop + local_up + cond_up, flip_pair + local_up,
                  local_flip + [1, 1, 1, 0] + [4, 0, 4, 1],
                  hop + local_up + cond_up],
        'four': [hop + local_up + cond_up + [2, 1, 2, 1],
                 flip_pair + local_up + cond_up,
                 plus_minus + cond_up + [2, 1, 2, 1],
                 hop + local_up + cond_up + [2, 1, 2, 1]],
        'six': [hop + local_up * 5,
                flip_pair + local_up * 4,
                plus_minus + local_up * 4],
    }
    requests['six'].append(requests['six'][0])
    five_factor = [row + local_up for row in requests['four']]
    requests['nbody'] = [[5] + row for row in five_factor]
    return requests


def add_correlation_requests(path: Path, requests: dict[str, list[list[int]]],
                             aggregate: bool) -> None:
    namelist = (path / 'sym.def').read_text()
    for kind, rows in requests.items():
        definition(path, CORRELATION_FILES[kind] + '.def', rows)
        namelist += f'{CORRELATION_KEYWORDS[kind]} {CORRELATION_FILES[kind]}.def\n'
    (path / 'sym.def').write_text(namelist)
    if aggregate:
        with (path / 'calc.def').open('a') as handle:
            handle.write('OutputGreenFormat 1\n')


def correlation_reference(case: Case, ref, vector: np.ndarray,
                          requests: dict[str, list[list[int]]]) -> dict[str, np.ndarray]:
    psi = ref.basis @ vector
    values = {}
    for kind, rows in requests.items():
        result = []
        for row in rows:
            data = row[1:] if kind == 'nbody' else row
            factors = tuple(tuple(data[i:i + 4]) for i in range(0, len(data), 4))
            result.append(np.vdot(psi, apply_operator(case, factors, psi)))
        result = np.asarray(result)
        assert np.isfinite(result).all() and np.max(abs(result)) > 1e-7, (kind, result)
        assert np.max(abs(result.imag)) > 1e-8, (kind, result)
        assert abs(result[0] - result[-1]) < 1e-12, (kind, result)
        values[kind] = result
    assert abs(values['two'][2] - values['two'][3]) > 1e-6
    return values


def check_correlation_outputs(path: Path, aggregate: bool,
                              requests: dict[str, list[list[int]]],
                              expected: dict[str, np.ndarray]) -> None:
    for kind, rows in requests.items():
        stem = CORRELATION_FILES[kind]
        filename = f'zvo_{stem}_eigen.dat' if aggregate else f'zvo_{stem}_eigen0.dat'
        prefix = 1 if aggregate else 0
        data = np.loadtxt(path / 'output' / filename, ndmin=2)
        assert data.shape == (len(rows), len(rows[0]) + 2 + prefix), (filename, data.shape)
        if aggregate:
            np.testing.assert_array_equal(data[:, 0], 0)
        np.testing.assert_array_equal(data[:, prefix:-2], np.asarray(rows))
        actual = data[:, -2] + 1j * data[:, -1]
        np.testing.assert_allclose(actual, expected[kind], atol=1e-8, rtol=0, err_msg=filename)


def check_correlations(root: Path, executable: Path, launcher: list[str]) -> None:
    requests = correlation_requests()
    for model in ('Kondo', 'KondoGC'):
        case = Case(model, 3, 'block', None if model == 'KondoGC' else 2,
                    1 if model == 'Kondo' else None, 1)
        ref = make_reference(case)
        hamiltonian = sector_hamiltonian(ref)
        for aggregate in ((False, True) if model == 'KondoGC' else (False,)):
            form = 'aggregate' if aggregate else 'legacy'
            path = expert_case(root / f'correlation_{model}_{form}', case, method=3,
                               options={'CalcMod': {'CalcEigenVec': 1},
                                        'ModPara': {'LanczosEps': 18}})
            add_correlation_requests(path, requests, aggregate)
            env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='distributed',
                       HPHI_TEST_SYMMETRY_CAPTURE='1')
            run_case(path, executable, 'correlation', launcher, env)
            vector = read_vector_parts(path, 'sector_final_sample0_step0.rank', len(hamiltonian))
            check_eigenvector(vector, hamiltonian, np.linalg.eigvalsh(hamiltonian)[0])
            check_correlation_outputs(path, aggregate, requests,
                                      correlation_reference(case, ref, vector, requests))

    invalid_case = Case('KondoGC', 3, 'block', None, None, 1)
    invalid = expert_case(root / 'invalid_local_offsite', invalid_case, method=3, options={})
    run_case(invalid, executable, 'invalid-correlation', launcher,
             dict(os.environ, HPHI_TEST_SYMMETRY_ACTION='invalid-correlation',
                  HPHI_SYMMETRY_BASIS_LAYOUT='distributed'))


def check_collective_failures(root: Path, probe: Path, launcher: list[str]) -> None:
    ranks = mpi_size(launcher)
    case = Case('KondoGC', 3, 'block', None, None, 1)
    ref = make_reference(case)
    vector = arbitrary_vector(ref.basis.shape[1], 19)
    for name in ('nan', 'storage'):
        path = expert_case(root / f'failure_{name}', case, method=3, options={})
        write_probe_inputs(path, vector, ranks)
        env = dict(os.environ, HPHI_TEST_SYMMETRY_ACTION='moments',
                   HPHI_SYMMETRY_BASIS_LAYOUT='distributed')
        if name == 'nan':
            env['HPHI_TEST_SYMMETRY_NAN_RANK'] = '0'
        else:
            env['HPHI_TEST_SYMMETRY_INVALID_STORAGE_RANK'] = '0'
        run_case(path, probe, name, launcher, env, 'Error: sector moments probe failed.')
    if ranks > 1:
        empty_case = Case('Kondo', 1, 'block', 0, 1, 0)
        empty_ref = make_reference(empty_case, empty=True)
        assert empty_ref.basis.shape[1] == 1 < ranks
        path = expert_case(root / 'empty_rank_moments', empty_case, method=3,
                           options={}, empty=True)
        result = run_probe(path, empty_case, probe, action='moments',
                           layout='distributed', vector=np.ones(1, complex))
        actual = parse_moments(result['text'])
        assert actual == {'N': 1., 'N2': 1., 'D': 0., 'D2': 0., 'Sz': .5, 'Sz2': .25}
        for name in ('nan', 'storage'):
            failure = expert_case(root / f'empty_rank_failure_{name}', empty_case,
                                  method=3, options={}, empty=True)
            write_probe_inputs(failure, np.ones(1, complex), ranks)
            env = dict(os.environ, HPHI_TEST_SYMMETRY_ACTION='moments',
                       HPHI_SYMMETRY_BASIS_LAYOUT='distributed')
            if name == 'nan':
                env['HPHI_TEST_SYMMETRY_NAN_RANK'] = '0'
            else:
                env['HPHI_TEST_SYMMETRY_INVALID_STORAGE_RANK'] = '0'
            run_case(failure, probe, f'empty_rank_{name}', launcher, env,
                     'Error: sector moments probe failed.')


def main() -> None:
    hphi, probe = [Path(arg).resolve() for arg in sys.argv[1:]]
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = mpi_size(launcher)
    assert ranks in (1, 2, 4, 16), 'Kondo observables suite requires serial, np2, np4, or np16'
    root = Path(tempfile.mkdtemp(prefix='symmetry_kondo_observables_', dir='.'))
    print(f'artifacts: {root.resolve()}', flush=True)
    scope = os.environ.get('HPHI_TEST_KONDO_OBSERVABLES_SCOPE', 'all')
    assert scope in ('all', 'moments', 'solvers', 'correlation', 'failures')
    if scope in ('all', 'moments'):
        for layout in ('replicated', 'distributed'):
            check_arbitrary_moments(root, probe, layout)
    if scope in ('all', 'solvers'):
        check_solvers(root, probe, launcher)
    if scope in ('all', 'correlation'):
        check_correlations(root, probe, launcher)
    if scope in ('all', 'failures'):
        check_collective_failures(root, probe, launcher)


if __name__ == '__main__':
    main()
