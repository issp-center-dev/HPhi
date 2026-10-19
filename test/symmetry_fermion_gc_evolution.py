#!/usr/bin/env python3
"""GC fermion sectors with pair H: TPQ/cTPQ samples and TE trajectories vs explicit updates."""
from pathlib import Path
import os
import shlex
import shutil
import sys
import numpy as np
from symmetry_kondo_common import merged_rows, read_manifest, read_vector_parts
from symmetry_kondo_observables import CORRELATION_FILES, add_correlation_requests
from symmetry_kondo_te import attach_seed, exact_endpoint_step, te_step
from symmetry_kondo_thermal import (EXPLICIT_BETAS, EXPLICIT_FLAGS, NSTEPS, NUM_AVE,
                                    UNIFORM_BETAS, captured_initials, check_correlation_block,
                                    check_data, check_metadata, ctpq_step, finite_data, tpq_step)
from symmetry_spingc_common import checkpoint_digest, definition, mpi_size, run_case
from symmetry_spingc_te import AMPLITUDES, TIMES, checkpoints, schedule
from symmetry_fermion_gc_common import (INPUT_MODEL, PARITY_FLAGS, expert_case, join_gc_vectors,
                                        read_gc_checkpoint)
from symmetry_fermion_gc_core import REFERENCE_TOL, anomalous_requests, normal_requests
from symmetry_fermion_gc_reference import (Case, MODELS, make_reference, operator_matrix,
                                           pair_operator, pair_rows, sector_matrix, width)
from symmetry_nconserved_reference import drive_terms as hubbard_drive_terms
from symmetry_nconserved_common import drive_rows as hubbard_drive_rows

MOMENTS = ('N', 'N2', 'D', 'D2', 'Sz', 'Sz2')
ERRORS = dict(vector=0., endpoint=0., data=0.)
CASES = {'HubbardGC': Case('HubbardGC', 4, 1, 1), 'SpinlessFermionGC': Case('SpinlessFermionGC', 6, 0, 1)}
SPINLESS_DRIVE = .07+.04j


def layout_env(layout, **extra):
    return dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout, HPHI_TEST_SYMMETRY_CAPTURE='1',
                **extra)


def drive_terms(case, kind, amplitude=1.):
    if case.model == 'HubbardGC':
        return hubbard_drive_terms(case.sites, kind, amplitude)
    assert kind == 'onebody', kind
    f = amplitude*SPINLESS_DRIVE
    terms = []
    for j in range(case.sites):
        nxt = (j+1) % case.sites
        terms += [(-f, [(1, nxt), (0, j)]), (-np.conj(f), [(1, j), (0, nxt)])]
    return terms


def drive_rows(case, kind):
    """TEOneBody (Transfer sign) or TETwoBody (InterAll) rows of drive_terms."""
    if case.model == 'HubbardGC':
        return hubbard_drive_rows(case.sites, kind)
    rows = []
    for j in range(case.sites):
        nxt = (j+1) % case.sites
        rows += [((nxt, 0, j, 0), SPINLESS_DRIVE), ((j, 0, nxt, 0), np.conj(SPINLESS_DRIVE))]
    return [(*indices, z.real, z.imag) for indices, z in merged_rows(rows)]


def drive_matrix(ref, case, kind, amplitude=1.):
    raw = operator_matrix(case, ref.words, drive_terms(case, kind, amplitude))
    return ref.basis.conj().T @ raw @ ref.basis


def moments(ref, vector):
    weights = abs(ref.basis @ vector)**2
    return {key: float(weights @ ref.moments[key]) for key in MOMENTS}


def requests_for(model):
    return normal_requests(model), anomalous_requests(model)


def operators_for(case, ref, requests, rows):
    """Raw-space matrices of the normal requests and the AnomalousG rows."""
    w = width(case.model)
    out = {}
    for kind, items in requests.items():
        out[kind] = []
        for row in items:
            data = row[1:] if kind == 'nbody' else row
            string = []
            for f in range(0, len(data), 4):
                i, s, j, t = data[f:f+4]
                string += [(1, w*i+s), (0, w*j+t)]
            out[kind].append(operator_matrix(case, ref.words, [(1, string)]))
    out['anomalous'] = [operator_matrix(case, ref.words, [(1, pair_operator(case.model, *row))])
                        for row in rows]
    return out


def expectations(ref, vector, operators):
    psi = ref.basis @ vector
    return {kind: np.array([np.vdot(psi, op @ psi) for op in ops]) for kind, ops in operators.items()}


def add_anomalous_request(path, rows):
    definition(path, 'anomalousg.def', rows)
    with (path/'sym.def').open('a') as f:
        f.write('AnomalousG anomalousg.def\n')


def check_anomalous_block(data, rows, expected, prefix, label):
    np.testing.assert_array_equal(data[:, prefix:prefix+5].astype(int), np.asarray(rows))
    np.testing.assert_allclose(data[:, -2]+1j*data[:, -1], expected, atol=REFERENCE_TOL, rtol=0,
                               err_msg=label)


def thermal_history(kind, h, initial, prenorm, ref, nsite, betas):
    """mTPQ (LargeValue 4) or cTPQ (12th-order exp(-dbeta H/2)) from the captured state."""
    vector, norm = initial.copy(), prenorm
    out = {'SS': [], 'Norm': [], 'Flct': [], 'vectors': []}
    for step in range(NSTEPS):
        out['vectors'].append(vector.copy())
        m = moments(ref, vector)
        hv = h @ vector
        energy = float(np.vdot(vector, hv).real)
        beta = 2*step/(nsite*4-energy) if kind == 'tpq' else float(betas[step])
        out['SS'].append([beta, energy, float(np.vdot(hv, hv).real), m['D'], m['N'], step])
        out['Norm'].append([beta, norm, prenorm, step])
        if kind != 'tpq' or step:
            out['Flct'].append([beta, m['N'], m['N2'], m['D'], m['D2'], m['Sz'], m['Sz2'], step])
        if step+1 < NSTEPS:
            if kind == 'tpq':
                vector, norm = tpq_step(h, vector, nsite)
            else:
                vector, norm = ctpq_step(h, vector, float(betas[step+1]-betas[step]))
    return {key: value if key == 'vectors' else np.asarray(value) for key, value in out.items()}


def check_thermal_correlations(path, aggregate, emitted, history, ref, operators, requests, rows):
    families = dict(requests, anomalous=rows)
    for kind, items in families.items():
        stem = 'AnomalousG' if kind == 'anomalous' else CORRELATION_FILES[kind]
        expected = [[expectations(ref, v, operators)[kind] for v in item['vectors']]
                    for item in history]

        def check(data, prefix, sample, step):
            label = f'{kind} {sample} {step}'
            if kind == 'anomalous':
                check_anomalous_block(data, items, expected[sample][step], prefix, label)
            else:
                check_correlation_block(data, items, expected[sample][step], prefix, label,
                                        atol=REFERENCE_TOL)
        if aggregate:
            data = finite_data(path/'output'/f'zvo_{stem}_tpq.dat')
            prefixes = [(s, step) for s in range(NUM_AVE) for step in emitted for _ in items]
            np.testing.assert_array_equal(data[:, :2], prefixes)
            offset = 0
            for sample in range(NUM_AVE):
                for step in emitted:
                    check(data[offset:offset+len(items)], 2, sample, step)
                    offset += len(items)
        else:
            for sample in range(NUM_AVE):
                for step in range(NSTEPS):
                    file = path/'output'/f'zvo_{stem}_set{sample}step{step}.dat'
                    if step not in emitted:
                        assert not file.exists(), file
                    else:
                        check(finite_data(file), 0, sample, step)


def run_thermal(root, hphi, launcher, case, kind, explicit, aggregate):
    pairs = pair_rows(case.model, case.sites)
    ref = make_reference(case, pairs=pairs)
    h = sector_matrix(ref)
    schedule_name = 'explicit' if explicit else 'uniform'
    path = expert_case(root/f'{case.model}_{kind}_{schedule_name}', case, pairs=pairs,
                       method=1 if kind == 'tpq' else 5, options={
        'CalcMod': {'InitialVecType': 1 if case.momentum == 0 else 0,
                    'OutputGreenFormat': int(aggregate)},
        'ModPara': {'Lanczos_max': NSTEPS, 'initial_iv': 7, 'NumAve': NUM_AVE,
                    'ExpecInterval': 2, 'LargeValue': 4 if kind == 'tpq' else 50,
                    **({'ExpandCoef': 12} if kind == 'ctpq' else {})}})
    requests, rows = requests_for(case.model)
    add_correlation_requests(path, requests, aggregate)
    add_anomalous_request(path, rows)
    if explicit:
        (path/'beta.def').write_text(''.join(f'{beta} 12 {flag} 0\n'
                                             for beta, flag in zip(EXPLICIT_BETAS, EXPLICIT_FLAGS)))
        with (path/'sym.def').open('a') as f:
            f.write('InvTemp beta.def\n')
        betas, emitted = EXPLICIT_BETAS, [0, 1, 3]
    else:
        betas, emitted = (UNIFORM_BETAS if kind == 'ctpq' else None), [0, 2, 4]
    layout = os.environ.get('HPHI_SYMMETRY_BASIS_LAYOUT', 'distributed')
    run_case(path, hphi, 'run', launcher, layout_env(layout))
    vectors, prenorms = captured_initials(path, len(h), ref.representatives, mpi_size(launcher))
    history = [thermal_history(kind, h, v, p, ref, case.sites, betas)
               for v, p in zip(vectors, prenorms)]
    check_metadata(path, case, kind, mpi_size(launcher), explicit, len(h))
    check_data(path, aggregate, history)
    check_thermal_correlations(path, aggregate, emitted, history, ref,
                               operators_for(case, ref, requests, rows), requests, rows)
    # N fluctuates: Flct reports the state value of a GC sample.
    flct = np.vstack([item['Flct'] for item in history])
    assert np.max(flct[:, 2]-flct[:, 1]**2) > 1e-4
    print(f'{case.model} {kind} {schedule_name} dim={len(h)} passed', flush=True)


def seed_cg(root, name, case, hphi, launcher, env, pairs=True):
    path = expert_case(root/name, case, method=3, pairs=pair_rows(case.model, case.sites) if pairs else None,
                       options={'CalcMod': {'OutputEigenVec': 1}, 'ModPara': {'LanczosEps': 18}})
    run_case(path, hphi, 'cg_seed', launcher, env)
    return path


def seed_fulldiag(root, name, case, hphi, env):
    """LAPACK FullDiag ground state on one rank, read by TE with any rank count."""
    path = expert_case(root/name, case, method=2, pairs=pair_rows(case.model, case.sites),
                       options={'CalcMod': {'OutputEigenVec': 1}})
    run_case(path, hphi, 'fd_seed', [], dict(env, HPHI_SYMMETRY_BASIS_LAYOUT='replicated'))
    (path/'output'/'zvo_eigenvec_0_rank_0.dat').rename(path/'output'/'fdseed_rank_0.dat')
    return path


def te_hamiltonians(case, drive):
    """free_quench evolves without pairs; every other drive keeps the static pairs."""
    pairs = None if drive == 'free_quench' else pair_rows(case.model, case.sites)
    ref = make_reference(case, pairs=pairs)
    base = sector_matrix(ref)
    if drive in ('static', 'pair_quench', 'free_quench'):
        return ref, [base]*len(TIMES)
    return ref, [base+drive_matrix(ref, case, drive, amplitude) for amplitude in AMPLITUDES]


def evolve(root, name, case, source, label, drive, hphi, launcher, env, interval=1):
    ref, hamiltonians = te_hamiltonians(case, drive)
    pairs = None if drive == 'free_quench' else pair_rows(case.model, case.sites)
    path = expert_case(root/name, case, method=4, pairs=pairs, options={
        'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1},
        'ModPara': {'Lanczos_max': len(TIMES), 'ExpandCoef': 8, 'ExpecInterval': interval,
                    'OutputInterval': interval}})
    requests, rows = requests_for(case.model)
    add_correlation_requests(path, requests, False)
    add_anomalous_request(path, rows)
    attach_seed(path, source, label)
    family = 'TETwoBody' if drive == 'twobody' else 'TEOneBody'
    schedule(path, family, drive_rows(case, drive) if drive in ('onebody', 'twobody') else [])
    text = run_case(path, hphi, 'evolve', launcher, env)
    ranks = mpi_size(launcher)
    dim = ref.basis.shape[1]
    initial = join_gc_vectors(checkpoints(source, label), dim, case)
    np.testing.assert_array_equal(read_vector_parts(path, 'sector_initial_sample0_step0.rank', dim),
                                  initial)
    meta = read_manifest(path)
    dynamic = drive in ('onebody', 'twobody')
    assert meta['te_hamiltonian'] == ('time_dependent' if dynamic else 'static')
    if dynamic:
        assert f'all {len(TIMES)} Hamiltonians passed preflight' in text
    if drive.endswith('quench'):
        assert 'hamiltonian_changed=yes' in text
    operators = operators_for(case, ref, requests, rows)
    vector, exact = initial.copy(), initial.copy()
    ss, norms, flct, anomalous = [], [], [], []
    for step, (time, h) in enumerate(zip(TIMES, hamiltonians)):
        np.testing.assert_allclose(h, h.conj().T, atol=1e-12, rtol=0)
        dt = time-(TIMES[step-1] if step else TIMES[0])
        work, norm = te_step(h, vector, dt)
        vector = work/norm
        exact = exact_endpoint_step(h, exact, dt)
        ERRORS['endpoint'] = max(ERRORS['endpoint'], np.linalg.norm(vector-exact))
        assert np.linalg.norm(vector-exact) <= 1e-8 and abs(norm-1) <= 1e-8
        hv = h @ vector
        m = moments(ref, vector)
        ss.append([time, np.vdot(vector, hv).real, np.vdot(hv, hv).real, m['D'], m['N'], step])
        norms.append([time, norm, step])
        flct.append([time, m['N'], m['N2'], m['D'], m['D2'], m['Sz'], m['Sz2'], step])
        files = checkpoints(path, f'zvo_eigenvec_{step}')
        if step % interval:
            assert not files
            continue
        assert len(files) == ranks
        error = np.linalg.norm(join_gc_vectors(files, dim, case)-vector)
        ERRORS['vector'] = max(ERRORS['vector'], error)
        assert error <= 3e-11, (name, step, error)
        for file in files:
            header, _ = read_gc_checkpoint(file, case)
            assert tuple(header[22:25]) == (4, 0, step)
        values = expectations(ref, vector, operators)
        for kind, items in requests.items():
            check_correlation_block(finite_data(path/'output'/f'zvo_{CORRELATION_FILES[kind]}_step{step}.dat'),
                                    items, values[kind], 0, f'{name} {kind} step{step}',
                                    atol=REFERENCE_TOL)
        check_anomalous_block(finite_data(path/'output'/f'zvo_AnomalousG_step{step}.dat'),
                              rows, values['anomalous'], 0, f'{name} AnomalousG step{step}')
        anomalous.append(values['anomalous'][0])
    for family_name, expected in (('SS', ss), ('Norm', norms), ('Flct', flct)):
        actual = finite_data(path/'output'/f'zvo_{family_name}.dat')
        expected = np.asarray(expected)
        np.testing.assert_array_equal(actual[:, -1], expected[:, -1])
        error = np.max(abs(actual-expected))
        ERRORS['data'] = max(ERRORS['data'], error)
        assert error <= 4e-11, (name, family_name, error)
    flct, anomalous = np.asarray(flct), np.asarray(anomalous)
    if drive == 'pair_quench':
        # From a fixed-N state, the pairs move weight between N and N +- 2. The
        # CG seed keeps other N only at its convergence level.
        assert np.ptp(flct[:, 1]) > 1e-6 and abs(anomalous[0]) < 1e-8, (name, flct[:, 1])
    if drive == 'free_quench':
        # An N-conserving H keeps N while F of the GC superposition evolves.
        assert np.ptp(flct[:, 1]) < 1e-12 and np.min(abs(anomalous)) > 1e-3, (name, anomalous)
        assert np.ptp(anomalous.real)+np.ptp(anomalous.imag) > 1e-6, anomalous
    finals = checkpoints(path, 'zvo_eigenvec_final')
    parts = [read_gc_checkpoint(file, case) for file in finals]
    assert len(parts) == ranks and all(tuple(h[26:28]) == checkpoint_digest(parts) for h, _ in parts)
    assert np.linalg.norm(join_gc_vectors(finals, dim, case)-vector) <= 3e-11
    print(f'{name} trajectory/outputs passed', flush=True)
    return path


def te_to_cg(root, case, source, hphi, launcher, env):
    """The final TE state is a CG input of the same space; the payload is unchanged."""
    path = expert_case(root/(source.name+'_to_cg'), case, method=3,
                       pairs=pair_rows(case.model, case.sites),
                       options={'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1}})
    (path/'output').mkdir()
    for rank, file in enumerate(checkpoints(source, 'zvo_eigenvec_final')):
        shutil.copy2(file, path/'output'/f'zvo_eigenvec_0_rank_{rank}.dat')
    run_case(path, hphi, 'te_to_cg', launcher, env)
    dim = make_reference(case).basis.shape[1]
    np.testing.assert_array_equal(join_gc_vectors(checkpoints(path, 'zvo_eigenvec_0'), dim, case),
                                  join_gc_vectors(checkpoints(source, 'zvo_eigenvec_final'), dim, case))
    for file in checkpoints(path, 'zvo_eigenvec_0'):
        assert read_gc_checkpoint(file, case)[0][22] == 3


def preflight_negatives(root, hphi, env):
    """A drive that breaks the translation is rejected before any TE output."""
    case = CASES['HubbardGC']
    seed = seed_cg(root, 'negative_seed', case, hphi, [], env)
    path = expert_case(root/'reject_noninvariant', case, method=4,
                       pairs=pair_rows(case.model, case.sites), options={
        'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1},
        'ModPara': {'Lanczos_max': len(TIMES), 'ExpandCoef': 8, 'ExpecInterval': 1,
                    'OutputInterval': 1}})
    attach_seed(path, seed, 'zvo_eigenvec_0')
    schedule(path, 'TEOneBody', [(0, 0, 0, 1, .3, 0.), (0, 1, 0, 0, .3, 0.)],
             amplitudes=[0, 0, 0, 1, 0])
    run_case(path, hphi, 'reject', [], env, 'preflight failed at step 3')
    assert not (path/'output/zvo_SS.dat').exists()
    assert not list((path/'output').glob('zvo_eigenvec_[0-9]*'))
    print('rejected noninvariant drive before TE output', flush=True)


def main():
    hphi = Path(sys.argv[1]).resolve()
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = mpi_size(launcher)
    assert ranks is not None
    root = Path('symmetry_fermion_gc_evolution')
    if root.exists():
        shutil.rmtree(root)
    root.mkdir()
    for model in MODELS:
        aggregate = model == 'SpinlessFermionGC'
        run_thermal(root, hphi, launcher, CASES[model], 'tpq', False, aggregate)
        run_thermal(root, hphi, launcher, CASES[model], 'ctpq', False, aggregate)
        run_thermal(root, hphi, launcher, CASES[model], 'ctpq', True, aggregate)
    layout = os.environ.get('HPHI_SYMMETRY_BASIS_LAYOUT', 'distributed')
    env = layout_env(layout)
    if ranks == 1:
        preflight_negatives(root, hphi, env)
    for model in MODELS:
        case = CASES[model]
        seed = seed_cg(root, model+'_seed', case, hphi, launcher, env)
        free_seed = seed_cg(root, model+'_free_seed', case, hphi, launcher, env, pairs=False)
        drives = ['static', 'free_quench', 'onebody'] + (['twobody'] if model == 'HubbardGC' else [])
        for drive in drives:
            result = evolve(root, f'{model}_{drive}', case, seed, 'zvo_eigenvec_0', drive,
                            hphi, launcher, env, interval=2 if drive == 'twobody' else 1)
            if drive == 'onebody':
                te_to_cg(root, case, result, hphi, launcher, env)
        evolve(root, f'{model}_pair_quench', case, free_seed, 'zvo_eigenvec_0', 'pair_quench',
               hphi, launcher, env)
        other = 'replicated' if layout == 'distributed' else 'distributed'
        evolve(root, f'{model}_onebody_{other}', case, seed, 'zvo_eigenvec_0', 'onebody',
               hphi, launcher, layout_env(other))
        fd = seed_fulldiag(root, model+'_fd_seed', case, hphi, env)
        evolve(root, f'{model}_fd_static', case, fd, 'fdseed', 'static', hphi, launcher, env)
    print('maximum errors:', {k: float(v) for k, v in ERRORS.items()})
    print('Fermion GC TPQ/cTPQ/TE passed')


if __name__ == '__main__':
    main()
