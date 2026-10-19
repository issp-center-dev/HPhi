#!/usr/bin/env python3
"""Fixed-N spin-mixing sectors: moments, CG/Lanczos/FullDiag and correlations."""
from pathlib import Path
import os
import re
import shlex
import shutil
import sys
import numpy as np
from symmetry_kondo_common import read_manifest, read_vector_parts, run_probe, selected_layouts
from symmetry_kondo_observables import (add_correlation_requests, arbitrary_vector,
                                        check_correlation_outputs, check_eigenvector,
                                        parse_energy, parse_moments)
from symmetry_spingc_common import mpi_size, run_case
from symmetry_nconserved_common import expert_case
from symmetry_nconserved_reference import (Case, DEFAULT, MODELS, make_reference,
                                           operator_matrix, sector_matrix)

LAYOUT = os.environ.get('HPHI_SYMMETRY_BASIS_LAYOUT', 'distributed')


def check_moments(root, probe):
    """Arbitrary Sz-mixing sector vectors: N is fixed, Sz and Sz^2 are not."""
    for model in MODELS:
        case = Case(model, 4, 3, 1)
        ref = make_reference(case)
        vector = arbitrary_vector(ref.basis.shape[1], 3)
        for layout in selected_layouts():
            path = expert_case(root/f'moments_{model}_{layout}', case, method=3)
            actual = parse_moments(run_probe(path, case, probe, action='moments', layout=layout,
                                             vector=vector)['text'])
            weights = abs(ref.basis @ vector)**2
            expected = {key: float(weights @ value) for key, value in ref.moments.items()}
            assert expected['Sz2']-expected['Sz']**2 > 1e-3, expected
            assert expected['N'] == 3 and expected['N2'] == 9
            if model == 'tJNConserved':
                assert expected['D'] == 0
            for key, value in actual.items():
                assert abs(value-expected[key]) <= 1e-10, (model, key, value, expected[key])


def check_solvers(root, executable, launcher):
    for model in MODELS:
        for k, method in ((1, 3), (2, 3), (1, 0)):
            case = Case(model, 4, 3, k)
            ref = make_reference(case)
            h = sector_matrix(ref)
            name = {3: 'cg', 0: 'lanczos'}[method]
            path = expert_case(root/f'{model}_k{k}_{name}', case, method=method,
                               options={'CalcMod': {'CalcEigenVec': 1},
                                        'ModPara': {'LanczosEps': 18}})
            env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=LAYOUT, HPHI_TEST_SYMMETRY_CAPTURE='1')
            if method == 0:
                env['HPHI_TEST_SYMMETRY_ACTION'] = 'lanczos'
            run_case(path, executable, name, launcher, env)
            reported, text = parse_energy(path)
            vector = read_vector_parts(path, 'sector_final_sample0_step0.rank', len(h))
            energy = check_eigenvector(vector, h, np.linalg.eigvalsh(h)[0])
            assert abs(reported-energy) <= 3e-8
            weights = abs(ref.basis @ vector)**2
            for key, label in (('Sz', 'Sz'), ('D', 'Doublon')):
                actual = float(re.search(rf'^\s*{label}\s+([^\s]+)', text, re.M).group(1))
                assert abs(actual-weights @ ref.moments[key]) <= 1e-8, (model, k, key)
            # Sz is evaluated from the state, never copied from an input value.
            assert abs(weights @ ref.moments['Sz']) > 1e-3


def orbital_factor(i, s, j, t):
    return [i, s, j, t]


def correlation_requests():
    hop = orbital_factor(1, 0, 0, 0)
    flip = orbital_factor(0, 0, 0, 1)
    mixed_hop = orbital_factor(1, 0, 0, 1)
    density = orbital_factor(2, 1, 2, 1)
    double_flip = flip + orbital_factor(1, 0, 1, 1)
    exchange = flip + orbital_factor(1, 1, 1, 0)
    return {
        'one': [hop, flip, mixed_hop, hop],
        'two': [hop + density, double_flip, exchange, mixed_hop + density, hop + density],
        'three': [hop + density + flip, double_flip + density, hop + density + flip],
        'nbody': [[3] + hop + density + flip, [3] + double_flip + density,
                  [3] + exchange + density, [3] + hop + density + flip],
    }


def correlation_operators(case, ref, requests):
    """Raw-space matrices of each requested ordered c^dagger c product."""
    operators = {}
    for kind, rows in requests.items():
        operators[kind] = []
        for row in rows:
            data = row[1:] if kind == 'nbody' else row
            string = []
            for f in range(0, len(data), 4):
                i, s, j, t = data[f:f+4]
                string += [(1, 2*i+s), (0, 2*j+t)]
            operators[kind].append(operator_matrix(case, ref.words, [(1, string)]))
    return operators


def expectations(ref, vector, operators):
    psi = ref.basis @ vector
    return {kind: np.array([np.vdot(psi, op @ psi) for op in ops])
            for kind, ops in operators.items()}


def correlation_reference(case, ref, vector, requests):
    values = expectations(ref, vector, correlation_operators(case, ref, requests))
    for kind, result in values.items():
        assert np.isfinite(result).all() and abs(result[0]-result[-1]) < 1e-12
        # Spin-flip rows are nonzero only because Sz is not fixed.
        assert abs(result[1]) > 1e-6, (kind, result)
    return values


def check_correlations(root, executable, launcher):
    requests = correlation_requests()
    for model in MODELS:
        case = Case(model, 4, 3, 1)
        ref = make_reference(case)
        h = sector_matrix(ref)
        for aggregate in ((False, True) if model == 'HubbardNConserved' else (False,)):
            form = 'aggregate' if aggregate else 'legacy'
            path = expert_case(root/f'correlation_{model}_{form}', case, method=3,
                               options={'CalcMod': {'CalcEigenVec': 1},
                                        'ModPara': {'LanczosEps': 18}})
            add_correlation_requests(path, requests, aggregate)
            env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=LAYOUT, HPHI_TEST_SYMMETRY_CAPTURE='1')
            run_case(path, executable, 'correlation', launcher, env)
            vector = read_vector_parts(path, 'sector_final_sample0_step0.rank', len(h))
            check_eigenvector(vector, h, np.linalg.eigvalsh(h)[0])
            check_correlation_outputs(path, aggregate, requests,
                                      correlation_reference(case, ref, vector, requests))


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
    for solver, label, backend_launcher in backends:
        for model in MODELS:
            case = Case(model, 4, 3, 1)
            ref = make_reference(case)
            path = expert_case(root/f'fulldiag_{model}_solver{solver}', case, method=2,
                               options={'CalcMod': {'Solver': solver, 'NGPU': 0}})
            text = run_case(path, hphi, 'fulldiag', backend_launcher, env)
            assert f'FullDiag solver: {label}' in text, text
            error = np.max(abs(fulldiag_spectrum(path)-np.linalg.eigvalsh(sector_matrix(ref))))
            assert error <= 3e-8, (model, label, error)
            assert read_manifest(path)['model'] == model
            print(f'{model} FullDiag {label} error={error:.3e} passed', flush=True)


def check_sz_conserving_limit(root, hphi):
    """With spin mixing switched off, every fixed-Sz sector of the same k adds up."""
    parameters = dict(DEFAULT, soc=0, soc2=0, hx=0, w=0)
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated')
    for model in MODELS:
        case = Case(model, 4, 3, 1)
        path = expert_case(root/f'limit_{model}', case, method=2, parameters=parameters)
        run_case(path, hphi, 'nconserved', [], env)
        combined = []
        for sz2 in (-3, -1, 1, 3):
            path = expert_case(root/f'limit_{model}_2Sz{sz2}', case, method=2,
                               parameters=parameters, options={'ModPara': {'2Sz': sz2}})
            run_case(path, hphi, 'fixed_sz', [], env)
            assert read_manifest(path)['model'] == model.replace('NConserved', '')
            combined += list(fulldiag_spectrum(path))
        np.testing.assert_allclose(np.sort(combined), fulldiag_spectrum(root/f'limit_{model}'),
                                   rtol=0, atol=1e-10)


def check_raw(root, hphi):
    """The raw Ncond-only path is unchanged and reaches the lowest sector energy."""
    for model in MODELS:
        energies = [np.linalg.eigvalsh(sector_matrix(make_reference(Case(model, 4, 3, k))))[0]
                    for k in range(4)]
        path = expert_case(root/f'raw_{model}', Case(model, 4, 3, 0), method=3,
                           options={'ModPara': {'LanczosEps': 18}})
        names = (path/'sym.def').read_text().replace('TransSym group.def\n', '')
        (path/'sym.def').write_text(names)
        # Raw MPI site splitting is not part of this suite; run it serially.
        run_case(path, hphi, 'raw', [], dict(os.environ))
        assert not (path/'output/zvo_symmetry_sector.dat').exists()
        assert abs(parse_energy(path)[0]-min(energies)) <= 3e-8


def check_rejections(root, hphi, launcher):
    case = Case('HubbardNConserved', 4, 3, 1)
    path = expert_case(root/'reject_calcspec', case, method=3,
                       options={'CalcMod': {'CalcSpec': 1}})
    run_case(path, hphi, 'reject', launcher, dict(os.environ),
             expected_error='spectrum calculations')


def main():
    hphi, probe = (Path(p).resolve() for p in sys.argv[1:3])
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    assert mpi_size(launcher) is not None
    root = Path('symmetry_nconserved_core')
    # Captured vectors of a run with another rank count must not be reread.
    if root.exists():
        shutil.rmtree(root)
    # The ELPA job selects only the dense backends with scope=fulldiag.
    scope = os.environ.get('HPHI_TEST_NCONSERVED_CORE_SCOPE', 'all')
    assert scope in ('all', 'fulldiag'), scope
    if scope == 'all':
        check_rejections(root, hphi, launcher)
        check_moments(root, probe)
        check_solvers(root, probe, launcher)
        check_correlations(root, probe, launcher)
    check_fulldiag(root, hphi, launcher)
    if scope == 'all':
        check_sz_conserving_limit(root, hphi)
        check_raw(root, hphi)
    print(f'NConserved core ({scope}) passed')


if __name__ == '__main__':
    main()
