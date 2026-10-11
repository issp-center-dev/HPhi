"""FullDiag sector states as TE initial states and CG re-evaluation inputs."""
from pathlib import Path
import os
import shutil
import struct
import sys

import numpy as np

import symmetry_fulldiag_common as common
import symmetry_spingc_common as spingc
from symmetry_fulldiag_common import W

# Sector TE needs an explicit time grid; zero TEOneBody rows keep H static.
TIMES = (0., .001, .002, .003, .004)
QUENCH = .17
ATOL = 1e-8
FULLDIAG, CG, TE = 2, 3, 4
KIND_LINE = 'source_hamiltonian_kind=sector_matrix_v1 hamiltonian_changed=unknown'


def exact_state(h, initial, time):
    energy, rotation = np.linalg.eigh(h)
    return rotation @ (np.exp(-1j*time*energy) * (rotation.conj().T @ initial))


def state_distance(actual, expected):
    overlap = np.vdot(expected, actual)
    assert abs(overlap) > 0
    return np.linalg.norm(actual*np.conj(overlap)/abs(overlap)-expected)


def environment(layout='replicated'):
    env = dict(os.environ)
    env.pop('HPHI_SYMMETRY_BASIS_LAYOUT', None)
    if layout is not None:
        env['HPHI_SYMMETRY_BASIS_LAYOUT'] = layout
    return env


def consumer(case, name):
    """A copy of the producer inputs without outputs."""
    target = case.path.parent/name
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(case.path, target, ignore=shutil.ignore_patterns('output', '*.log'))
    (target/'output').mkdir()
    return target


def produce(case, hphi, solver, launcher, label):
    """FullDiag with every state saved; the states are kept in a store."""
    if (case.path/'output').exists():
        shutil.rmtree(case.path/'output')
    common.run_fd(case, hphi, solver, {'OutputEigenVec': 1, 'OutputGreenFormat': 0},
                  launcher, f'produce_{label}')
    store = case.path.parent/f'{case.name}_{label}_states'
    if store.exists():
        shutil.rmtree(store)
    shutil.copytree(case.path/'output', store)
    return store


def copy_state(store, state, target, name, ranks):
    for rank in range(ranks):
        shutil.copy2(store/f'zvo_eigenvec_{state}_rank_{rank}.dat',
                     target/'output'/f'{name}_rank_{rank}.dat')


def set_te_inputs(path, quench_rows=None):
    common.set_calc_file(path, {'CalcType': TE, 'InputEigenVec': 1, 'OutputEigenVec': 1,
                                'OutputGreenFormat': 0})
    common.set_mod_file(path, {'Lanczos_max': len(TIMES), 'ExpandCoef': 12,
                               'OutputInterval': 1, 'ExpecInterval': 1})
    (path/'TEOneBody.def').write_text(f'====\nNTimeSteps {len(TIMES)}\n====\n====\n====\n'
                                      + ''.join(f'{t} 0\n' for t in TIMES))
    names = (path/'sym.def').read_text()
    for line in ('SpectrumVec seed\n', 'TEOneBody TEOneBody.def\n'):
        if line not in names:
            names += line
    (path/'sym.def').write_text(names)
    if quench_rows:
        lines = (path/'InterAll.def').read_text().splitlines()
        rows = [line.split() for line in lines[5:] if line.strip()] + \
            [list(map(str, row)) for row in quench_rows]
        spingc.definition(path, 'InterAll.def', rows)


def te_states(path, ranks, steps):
    states, headers = [], []
    for step in range(steps):
        vector, parts = common.read_fd_state(path/'output', step, ranks, expect_state=None)
        states.append(vector)
        headers.append(parts[0])
    return states, headers


def word_time(word):
    return struct.unpack('<d', struct.pack('<Q', int(word)))[0]


def check_te(case, path, hphi, launcher, ranks, label, seed, state, h, source_digest):
    text = spingc.run_case(path, hphi, label, launcher, environment())
    assert KIND_LINE in text, text
    manifest = common.read_manifest(path/'output')
    assert int(manifest['source_method']) == FULLDIAG and int(manifest['source_state']) == state
    assert manifest['source_hamiltonian_digest'] == f'{source_digest:016x}'
    assert manifest['source_hamiltonian_kind'] == 'sector_matrix_v1', manifest
    parsed = int(manifest['hamiltonian_digest'].rsplit(':', 1)[1], 16)
    states, headers = te_states(path, ranks, len(TIMES))
    worst = 0.
    for step, (vector, header) in enumerate(zip(states, headers)):
        # TE writes its own H: the parsed digest of the consumer input.
        assert int(header[W.METHOD]) == TE and int(header[W.STATE]) == state
        assert int(header[W.STEP]) == step and int(header[W.HAMILTONIAN]) == parsed
        time = word_time(header[W.TIME])
        assert abs(time - TIMES[step]) < 1e-15, (time, TIMES[step])
        np.testing.assert_allclose(np.linalg.norm(vector), 1, rtol=0, atol=1e-10)
        expected = exact_state(h, seed, time)
        distance = state_distance(vector, expected)
        assert distance <= ATOL, (case.name, label, step, distance)
        worst = max(worst, distance)
        for family, requests in case.operators.items():
            actual = common.read_green_file(path/'output'/f'zvo_{common.GREEN_STEMS[family]}_step{step}.dat',
                                            [row for row, _ in requests])
            reference = np.array([np.vdot(expected, oq @ expected) for _, oq in requests])
            np.testing.assert_allclose(actual.real, reference.real, rtol=0, atol=ATOL)
            np.testing.assert_allclose(actual.imag, reference.imag, rtol=0, atol=ATOL)
    return worst


def quench_case(case, states, required):
    """Uniform added term and a seed for which the quench is not trivial.

    A required quench (every model of the all selection) fails instead of
    skipping; a boundary fixture may be unsuitable (D=1, no oracle)."""
    if case.model in common.KONDO_MODELS and case.kondo is None:
        assert not required, case.name
        return None
    rows, delta = common.uniform_quench(case, QUENCH)
    after = case.hq + delta
    commutator = np.linalg.norm(case.hq @ after - after @ case.hq)
    if commutator <= 1e-8:
        assert not required, (case.name, commutator)
        return None
    for state in range(len(case.hq)):
        seed = states[:, state]
        variance = np.vdot(seed, after @ after @ seed).real - np.vdot(seed, after @ seed).real**2
        if variance > 1e-8:
            return rows, after, state
    assert not required, f'{case.name}: no state has a nonzero quench variance'
    return None


def check_cg(case, path, hphi, launcher, ranks, label, vectors, energies):
    count = min(2, len(case.hq))
    common.set_calc_file(path, {'CalcType': CG, 'InputEigenVec': 1, 'OutputEigenVec': 0,
                                'OutputGreenFormat': 0})
    common.set_mod_file(path, {'exct': count})
    text = spingc.run_case(path, hphi, label, launcher, environment())
    assert 'An Eigenvector is inputted.' in text and 'Residual-2-norm' not in text, text
    assert KIND_LINE in text
    energy = [float(line.split()[1]) for line in (path/'output'/'zvo_energy.dat').read_text().splitlines()
              if line.strip().startswith('Energy')]
    assert len(energy) == count
    for n in range(count):
        v = vectors[:, n]
        np.testing.assert_allclose(energy[n], np.vdot(v, case.hq @ v).real, rtol=0, atol=ATOL)
        np.testing.assert_allclose(energy[n], energies[n], rtol=0, atol=ATOL)
    for family, requests in case.operators.items():
        actual = common.read_green(path/'output', family, count, False, [row for row, _ in requests])
        for n in range(count):
            v = vectors[:, n]
            reference = np.array([np.vdot(v, oq @ v) for _, oq in requests])
            np.testing.assert_allclose(actual[n].real, reference.real, rtol=0, atol=ATOL)
            np.testing.assert_allclose(actual[n].imag, reference.imag, rtol=0, atol=ATOL)
    assert not list((path/'output').glob('zvo_eigenvec_*_rank_*.dat.part'))


def negatives(cases, case, store, hphi, launcher, ranks, label):
    """Mismatched state imports are rejected before any evolution."""
    header_error = 'checkpoint header / sector / layout validation failed'
    path = consumer(case, f'{case.name}_{label}_default_layout')
    copy_state(store, 0, path, 'seed', ranks)
    set_te_inputs(path)
    spingc.run_case(path, hphi, 'te_default_layout', launcher, environment(None), header_error)
    if ranks > 1:
        path = consumer(case, f'{case.name}_{label}_rank_count')
        copy_state(store, 0, path, 'seed', ranks)
        set_te_inputs(path)
        spingc.run_case(path, hphi, 'te_rank_count', [], environment(), header_error)
    same = next(other for other in cases if other.model == case.model and other is not case
                and len(other.hq) == len(case.hq))
    path = consumer(same, f'{same.name}_{label}_other_sector')
    copy_state(store, 0, path, 'zvo_eigenvec_0', ranks)
    common.set_calc_file(path, {'CalcType': CG, 'InputEigenVec': 1, 'OutputEigenVec': 0})
    common.set_mod_file(path, {'exct': 1})
    spingc.run_case(path, hphi, 'cg_other_sector', launcher, environment(), header_error)
    if len(case.hq) > 1:
        path = consumer(case, f'{case.name}_{label}_state_number')
        copy_state(store, 1, path, 'zvo_eigenvec_0', ranks)
        common.set_calc_file(path, {'CalcType': CG, 'InputEigenVec': 1, 'OutputEigenVec': 0})
        common.set_mod_file(path, {'exct': 1})
        spingc.run_case(path, hphi, 'cg_state_number', launcher, environment(),
                        'symmetry checkpoint state index or step is invalid for CG input')
    print(f'{case.name} {label}: layout/rank/sector/state rejections PASS', flush=True)


def kondo_layout_negative(root, case, store, hphi, launcher, ranks, label):
    from symmetry_kondo_common import expert_case
    from symmetry_kondo_reference import Case
    other = Case(case.kondo.model, case.kondo.cells, 'alternating', case.kondo.ncond,
                 case.kondo.sz2, case.kondo.momentum)
    path = expert_case(root/f'{case.name}_{label}_alternating', other, method=CG,
                       options={'CalcMod': {'InputEigenVec': 1}, 'ModPara': {'exct': 1}})
    (path/'output').mkdir(exist_ok=True)
    copy_state(store, 0, path, 'zvo_eigenvec_0', ranks)
    spingc.run_case(path, hphi, 'cg_kondo_layout', launcher, environment(),
                    'checkpoint header / sector / layout validation failed')
    print(f'{case.name} {label}: Kondo local-site layout rejection PASS', flush=True)


def main():
    hphi = Path(sys.argv[1]).resolve()
    mode = sys.argv[2] if len(sys.argv) > 2 else 'all'
    assert mode in ('all', 'generated'), mode
    launcher, ranks = common.launcher_from_env()
    selection = common.selected_cases()
    root = common.suite_root(f'reuse-{selection}')
    cases = common.cases(root, selection)
    chosen, seen = [], set()
    for case in cases:
        if case.model not in seen and (selection == 'boundary' or len(case.hq) > 1):
            chosen.append(case)
            seen.add(case.model)
    if selection == 'boundary':
        chosen = cases
    worst, negatives_done, quenches = 0., set(), 0
    for case in chosen:
        if not common.attach_correlations(case, False):
            case.operators = {}
        for solver, label, backend_launcher in common.solver_backends(launcher, ranks):
            n = common.spingc.mpi_size(backend_launcher)
            if solver == 3 and len(case.hq) < common.grid_max(n):
                continue
            store = produce(case, hphi, solver, backend_launcher, label)
            dim = len(case.hq)
            energies = np.loadtxt(store/'zvo_energy_sector.dat', ndmin=2)[:, 1]
            vectors, headers = common.read_states(store, dim, n)
            source_digest = int(headers[0][0][W.HAMILTONIAN])
            state = min(1, dim-1)
            # Static H: the FullDiag state only acquires a phase.
            path = consumer(case, f'{case.name}_{label}_te')
            copy_state(store, state, path, 'seed', n)
            set_te_inputs(path)
            worst = max(worst, check_te(case, path, hphi, backend_launcher, n, 'te_static',
                                        vectors[:, state], state, case.hq, source_digest))
            quench = quench_case(case, vectors, selection == 'all') if dim > 1 else None
            if quench is not None:
                quenches += 1
                rows, after, quench_state = quench
                path = consumer(case, f'{case.name}_{label}_quench')
                copy_state(store, quench_state, path, 'seed', n)
                set_te_inputs(path, rows)
                worst = max(worst, check_te(case, path, hphi, backend_launcher, n, 'te_quench',
                                            vectors[:, quench_state], quench_state, after,
                                            source_digest))
            path = consumer(case, f'{case.name}_{label}_cg')
            for k in range(min(2, dim)):
                copy_state(store, k, path, f'zvo_eigenvec_{k}', n)
            check_cg(case, path, hphi, backend_launcher, n, 'cg', vectors, energies)
            print(f'{case.name} {label} np={n}: TE static/quench and CG re-evaluation PASS',
                  flush=True)
            if label not in negatives_done and selection == 'all':
                negatives(cases, case, store, hphi, backend_launcher, n, label)
                negatives_done.add(label)
            if case.model == 'Kondo' and selection == 'all' and f'kondo_{label}' not in negatives_done:
                kondo_layout_negative(root, case, store, hphi, backend_launcher, n, label)
                negatives_done.add(f'kondo_{label}')
    assert quenches > 0, 'no noncommuting quench was exercised'
    print(f'sector FullDiag reuse ({selection}): {len(seen)} models, {quenches} quenches, '
          f'worst state error {worst:.2e} PASS', flush=True)


if __name__ == '__main__':
    main()
