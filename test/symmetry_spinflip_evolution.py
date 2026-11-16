#!/usr/bin/env python3
"""Px sectors: TPQ/cTPQ samples and TE trajectories against explicit updates of the same vector."""
from pathlib import Path
import os
import shlex
import shutil
import sys
import numpy as np
import symmetry_spingc_common as c
from symmetry_kondo_common import read_manifest, read_vector_parts
from symmetry_kondo_observables import CORRELATION_FILES, add_correlation_requests
from symmetry_kondo_te import attach_seed, exact_endpoint_step, te_step
from symmetry_kondo_thermal import (EXPLICIT_BETAS, EXPLICIT_FLAGS, NSTEPS, NUM_AVE,
                                    UNIFORM_BETAS, captured_initials, check_correlation_block,
                                    check_data, check_metadata, ctpq_step, finite_data, tpq_step)
from symmetry_spingc_common import checkpoint_digest, mpi_size, run_case
from symmetry_spingc_te import AMPLITUDES, TIMES, checkpoints, schedule
from symmetry_spinflip_common import check_manifest, expert_case, sy_sz_rows
from symmetry_spinflip_core import REFERENCE_TOL, correlation_requests
from symmetry_spinflip_reference import (SZ, Case, make_reference, nearest_bonds, raw_operator,
                                         sector_matrix, spin_operators)

ERRORS = dict(vector=0., endpoint=0., data=0.)
CASES = {'SpinGC': Case('SpinGC', 8, 1, 1), 'Spin': Case('Spin', 8, 1, 2)}
# Px-even drives: the transverse field hx(t), and two-body terms.
DRIVES = {'SpinGC': ('onebody', 'twobody'), 'Spin': ('twobody',)}


def layout_env(layout, **extra):
    return dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout, HPHI_TEST_SYMMETRY_CAPTURE='1',
                **extra)


def drive_operator(case, kind, amplitude=1.):
    """Raw 2^L matrix of a drive, from spin operators."""
    s = spin_operators(case.sites)
    bonds = nearest_bonds(case.sites)
    if kind == 'onebody':
        return -amplitude*.31*sum(s['x'])
    if kind == 'twobody' and case.model == 'SpinGC':
        return amplitude*.17*sum(s['y'][i]@s['z'][j] + s['z'][i]@s['y'][j] for i, j in bonds)
    if kind == 'twobody':
        return amplitude*.23*sum(s['z'][i]@s['z'][j] + (s['plus'][i]@s['minus'][j] +
                                                       s['minus'][i]@s['plus'][j])/2 for i, j in bonds)
    raise ValueError(kind)


def drive_rows(case, kind):
    """TEOneBody (Transfer sign) or TETwoBody (InterAll) rows of drive_operator."""
    n, bonds = case.sites, nearest_bonds(case.sites)
    if kind == 'onebody':
        return [(i, a, i, b, .31/2, 0.) for i in range(n) for a, b in ((1, 0), (0, 1))]
    if kind == 'twobody' and case.model == 'SpinGC':
        return [row for i, j in bonds for row in sy_sz_rows(i, j, .17)]
    if kind == 'twobody':
        rows = [(i, a, i, a, j, b, j, b, .23*(SZ[a, a]*SZ[b, b]).real, 0.)
                for i, j in bonds for a in range(2) for b in range(2)]
        return rows + [row for i, j in bonds for row in c.interall_pair((i, 1, i, 0, j, 0, j, 1), .23/2)]
    if kind == 'hz':
        return [(i, a, i, a, -.21*SZ[a, a].real, 0.) for i in range(n) for a in range(2)]
    raise ValueError(kind)


def drive_matrix(ref, case, kind, amplitude=1.):
    index = ref.words.astype(np.int64)
    raw = drive_operator(case, kind, amplitude)[np.ix_(index, index)]
    return ref.basis.conj().T @ raw @ ref.basis


def check_drive_rows():
    """The written rows are the drive operators (an independent tensor reading)."""
    for case in CASES.values():
        index = make_reference(case).words.astype(np.int64)
        for kind in DRIVES[case.model]:
            family = 'Trans' if kind == 'onebody' else 'InterAll'
            rows = c.families_matrix(case.sites, {family: drive_rows(case, kind)})
            error = np.linalg.norm((rows-drive_operator(case, kind))[np.ix_(index, index)], 2)
            assert error <= 1e-13, (case, kind, error)


def operators_for(case, ref, requests):
    out = {}
    for kind, items in requests.items():
        out[kind] = []
        for row in items:
            data = row[1:] if kind == 'nbody' else row
            out[kind].append(raw_operator(case, ref.words, [tuple(data[f:f+4])
                                                           for f in range(0, len(data), 4)]))
    return out


def expectations(ref, vector, operators):
    psi = ref.basis @ vector
    return {kind: np.array([np.vdot(psi, op @ psi) for op in ops]) for kind, ops in operators.items()}


def moments(case, ref, vector):
    """Flct/SS moments: N = L, no doublons, <Sz> = 0 in a Px sector, <Sz^2> from the state."""
    weights = abs(ref.basis @ vector)**2
    sz, sz2 = float(weights @ ref.moments['Sz']), float(weights @ ref.moments['Sz2'])
    assert abs(sz) < 1e-13, sz
    return dict(N=case.sites, N2=case.sites**2, D=0., D2=0., Sz=sz, Sz2=sz2)


def thermal_history(kind, case, h, initial, prenorm, ref, betas):
    """mTPQ (LargeValue 4) or cTPQ (12th-order exp(-dbeta H/2)) from the captured state."""
    vector, norm = initial.copy(), prenorm
    out = {'SS': [], 'Norm': [], 'Flct': [], 'vectors': []}
    n = case.sites
    for step in range(NSTEPS):
        out['vectors'].append(vector.copy())
        m = moments(case, ref, vector)
        hv = h @ vector
        energy = float(np.vdot(vector, hv).real)
        beta = 2*step/(n*4-energy) if kind == 'tpq' else float(betas[step])
        out['SS'].append([beta, energy, float(np.vdot(hv, hv).real), m['D'], m['N'], step])
        out['Norm'].append([beta, norm, prenorm, step])
        if kind != 'tpq' or step:
            out['Flct'].append([beta, m['N'], m['N2'], m['D'], m['D2'], m['Sz'], m['Sz2'], step])
        if step+1 < NSTEPS:
            if kind == 'tpq':
                vector, norm = tpq_step(h, vector, n)
            else:
                vector, norm = ctpq_step(h, vector, float(betas[step+1]-betas[step]))
    return {key: value if key == 'vectors' else np.asarray(value) for key, value in out.items()}


def check_thermal_correlations(path, aggregate, emitted, history, ref, operators, requests):
    for kind, items in requests.items():
        stem = CORRELATION_FILES[kind]
        expected = [[expectations(ref, v, operators)[kind] for v in item['vectors']]
                    for item in history]
        if aggregate:
            data = finite_data(path/'output'/f'zvo_{stem}_tpq.dat')
            prefixes = [(s, step) for s in range(NUM_AVE) for step in emitted for _ in items]
            np.testing.assert_array_equal(data[:, :2], prefixes)
            offset = 0
            for sample in range(NUM_AVE):
                for step in emitted:
                    check_correlation_block(data[offset:offset+len(items)], items,
                                            expected[sample][step], 2, f'{kind} {sample} {step}',
                                            atol=REFERENCE_TOL)
                    offset += len(items)
        else:
            for sample in range(NUM_AVE):
                for step in range(NSTEPS):
                    file = path/'output'/f'zvo_{stem}_set{sample}step{step}.dat'
                    if step not in emitted:
                        assert not file.exists(), file
                    else:
                        check_correlation_block(finite_data(file), items, expected[sample][step], 0,
                                                f'{kind} {sample} {step}', atol=REFERENCE_TOL)


def run_thermal(root, hphi, launcher, case, kind, explicit, aggregate):
    ref = make_reference(case)
    h = sector_matrix(ref)
    schedule_name = 'explicit' if explicit else 'uniform'
    path = expert_case(root/f'{case.model}_{kind}_{schedule_name}', case,
                       method=1 if kind == 'tpq' else 5, options={
        'CalcMod': {'InitialVecType': 1 if case.momentum == 0 else 0,
                    'OutputGreenFormat': int(aggregate)},
        'ModPara': {'Lanczos_max': NSTEPS, 'initial_iv': 7, 'NumAve': NUM_AVE,
                    'ExpecInterval': 2, 'LargeValue': 4 if kind == 'tpq' else 50,
                    **({'ExpandCoef': 12} if kind == 'ctpq' else {})}})
    requests = correlation_requests(case.model)
    add_correlation_requests(path, requests, aggregate)
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
    check_manifest(read_manifest(path), case, len(h))
    vectors, prenorms = captured_initials(path, len(h), ref.representatives, mpi_size(launcher))
    history = [thermal_history(kind, case, h, v, p, ref, betas) for v, p in zip(vectors, prenorms)]
    check_metadata(path, case, kind, mpi_size(launcher), explicit, len(h))
    check_data(path, aggregate, history)
    check_thermal_correlations(path, aggregate, emitted, history, ref,
                               operators_for(case, ref, requests), requests)
    if case.model == 'SpinGC':
        # HPhi reports <Sz> = 0 exactly; Sz^2 is a state value, not <Sz>^2.
        files = ([path/'output/zvo_Flct_tpq.dat'] if aggregate else
                 [path/'output'/f'zvo_Flct_rand{s}.dat' for s in range(NUM_AVE)])
        flct = np.vstack([finite_data(file)[:, -8:] for file in files])
        assert np.all(flct[:, 5] == 0) and np.min(flct[:, 6]) > .1, flct
    print(f'{case.model} {kind} {schedule_name} dim={len(h)} passed', flush=True)


def seed_cg(root, name, case, hphi, launcher, env):
    path = expert_case(root/name, case, method=3,
                       options={'CalcMod': {'OutputEigenVec': 1}, 'ModPara': {'LanczosEps': 18}})
    run_case(path, hphi, 'cg_seed', launcher, env)
    return path


def seed_fulldiag(root, name, case, hphi, env):
    """LAPACK FullDiag ground state on one rank, read by TE with any rank count."""
    path = expert_case(root/name, case, method=2, options={'CalcMod': {'OutputEigenVec': 1}})
    run_case(path, hphi, 'fd_seed', [], dict(env, HPHI_SYMMETRY_BASIS_LAYOUT='replicated'))
    (path/'output'/'zvo_eigenvec_0_rank_0.dat').rename(path/'output'/'fdseed_rank_0.dat')
    return path


QUENCH = dict(hx=.41)


def te_hamiltonians(case, drive):
    """quench evolves a ground state of the default H with another Px-even H."""
    parameters = QUENCH if drive == 'quench' else None
    ref = make_reference(case, parameters=parameters)
    base = sector_matrix(ref)
    if drive in ('static', 'quench'):
        return ref, [base]*len(TIMES)
    return ref, [base+drive_matrix(ref, case, drive, amplitude) for amplitude in AMPLITUDES]


def evolve(root, name, case, source, label, drive, hphi, launcher, env, interval=1):
    ref, hamiltonians = te_hamiltonians(case, drive)
    path = expert_case(root/name, case, method=4,
                       parameters=QUENCH if drive == 'quench' else None, options={
        'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1},
        'ModPara': {'Lanczos_max': len(TIMES), 'ExpandCoef': 8, 'ExpecInterval': interval,
                    'OutputInterval': interval}})
    requests = correlation_requests(case.model)
    add_correlation_requests(path, requests, False)
    attach_seed(path, source, label)
    family = 'TETwoBody' if drive == 'twobody' else 'TEOneBody'
    schedule(path, family, drive_rows(case, drive) if drive in ('onebody', 'twobody') else [])
    text = run_case(path, hphi, 'evolve', launcher, env)
    ranks = mpi_size(launcher)
    dim = ref.basis.shape[1]
    initial = c.join_rank_vectors(checkpoints(source, label), dim)
    np.testing.assert_array_equal(read_vector_parts(path, 'sector_initial_sample0_step0.rank', dim),
                                  initial)
    meta = read_manifest(path)
    check_manifest(meta, case, dim)
    dynamic = drive in ('onebody', 'twobody')
    assert meta['te_hamiltonian'] == ('time_dependent' if dynamic else 'static')
    if dynamic:
        assert f'all {len(TIMES)} Hamiltonians passed preflight' in text
        assert meta['te_plan_update'] == 'rebuild_with_fixed_basis'
    if drive == 'quench':
        assert 'hamiltonian_changed=yes' in text
    operators = operators_for(case, ref, requests)
    vector, exact = initial.copy(), initial.copy()
    ss, norms, flct, sx = [], [], [], []
    for step, (time, h) in enumerate(zip(TIMES, hamiltonians)):
        np.testing.assert_allclose(h, h.conj().T, atol=1e-12, rtol=0)
        dt = time-(TIMES[step-1] if step else TIMES[0])
        work, norm = te_step(h, vector, dt)
        vector = work/norm
        exact = exact_endpoint_step(h, exact, dt)
        ERRORS['endpoint'] = max(ERRORS['endpoint'], np.linalg.norm(vector-exact))
        assert np.linalg.norm(vector-exact) <= 1e-8 and abs(norm-1) <= 1e-8
        hv = h @ vector
        m = moments(case, ref, vector)
        ss.append([time, np.vdot(vector, hv).real, np.vdot(hv, hv).real, m['D'], m['N'], step])
        norms.append([time, norm, step])
        flct.append([time, m['N'], m['N2'], m['D'], m['D2'], m['Sz'], m['Sz2'], step])
        files = checkpoints(path, f'zvo_eigenvec_{step}')
        if step % interval:
            assert not files
            continue
        assert len(files) == ranks
        error = np.linalg.norm(c.join_rank_vectors(files, dim)-vector)
        ERRORS['vector'] = max(ERRORS['vector'], error)
        assert error <= 3e-11, (name, step, error)
        for file in files:
            header, _ = c.read_checkpoint(file)
            assert tuple(header[22:25]) == (4, 0, step)
        values = expectations(ref, vector, operators)
        for kind, items in requests.items():
            check_correlation_block(finite_data(path/'output'/f'zvo_{CORRELATION_FILES[kind]}_step{step}.dat'),
                                    items, values[kind], 0, f'{name} {kind} step{step}',
                                    atol=REFERENCE_TOL)
        sx.append(values['one'][0])
    for family_name, expected in (('SS', ss), ('Norm', norms), ('Flct', flct)):
        actual = finite_data(path/'output'/f'zvo_{family_name}.dat')
        expected = np.asarray(expected)
        np.testing.assert_array_equal(actual[:, -1], expected[:, -1])
        error = np.max(abs(actual-expected))
        ERRORS['data'] = max(ERRORS['data'], error)
        assert error <= 4e-11, (name, family_name, error)
    if case.model == 'SpinGC' and drive in ('onebody', 'quench'):
        # The field moves <S+> = <Sx>; it stays real since <Sy> = 0.
        sx = np.asarray(sx)
        assert np.ptp(sx.real) > 1e-6 and np.max(abs(sx.imag)) < 1e-12, sx
    finals = checkpoints(path, 'zvo_eigenvec_final')
    parts = [c.read_checkpoint(file) for file in finals]
    assert len(parts) == ranks and all(tuple(h[26:28]) == checkpoint_digest(parts) for h, _ in parts)
    assert np.linalg.norm(c.join_rank_vectors(finals, dim)-vector) <= 3e-11
    print(f'{name} trajectory/outputs passed', flush=True)
    return path


def te_to_cg(root, case, source, hphi, launcher, env):
    """The final TE state is a CG input of the same (k, p) sector; the payload is unchanged."""
    path = expert_case(root/(source.name+'_to_cg'), case, method=3,
                       options={'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1}})
    (path/'output').mkdir()
    for rank, file in enumerate(checkpoints(source, 'zvo_eigenvec_final')):
        shutil.copy2(file, path/'output'/f'zvo_eigenvec_0_rank_{rank}.dat')
    run_case(path, hphi, 'te_to_cg', launcher, env)
    dim = make_reference(case).basis.shape[1]
    np.testing.assert_array_equal(c.join_rank_vectors(checkpoints(path, 'zvo_eigenvec_0'), dim),
                                  c.join_rank_vectors(checkpoints(source, 'zvo_eigenvec_final'), dim))
    for file in checkpoints(path, 'zvo_eigenvec_0'):
        assert c.read_checkpoint(file)[0][22] == 3


def preflight_negatives(root, hphi, env):
    """A drive that breaks Px at one intermediate step stops before any TE output."""
    for model, family, kind in (('SpinGC', 'TEOneBody', 'hz'), ('Spin', 'TETwoBody', 'dm')):
        case = CASES[model]
        seed = seed_cg(root, f'negative_seed_{model}', case, hphi, [], env)
        path = expert_case(root/f'reject_{model}_{kind}', case, method=4, options={
            'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1},
            'ModPara': {'Lanczos_max': len(TIMES), 'ExpandCoef': 8, 'ExpecInterval': 1,
                        'OutputInterval': 1}})
        attach_seed(path, seed, 'zvo_eigenvec_0')
        if kind == 'dm':
            rows = [row for i, j in nearest_bonds(case.sites)
                    for row in c.interall_pair((i, 1, i, 0, j, 0, j, 1), .13j)]
        else:
            rows = drive_rows(case, kind)
        schedule(path, family, rows, amplitudes=[0, 0, 0, 1, 0])
        run_case(path, hphi, 'reject', [], env, 'preflight failed at step 3')
        log = (path/'reject.log').read_text()
        # The step-3 H fails under (e, Px), the first appended operation.
        assert f'invariance failed for normalized terms under op {case.sites}.' in log, log
        assert not (path/'output/zvo_SS.dat').exists()
        assert not list((path/'output').glob('zvo_eigenvec_[0-9]*'))
        print(f'rejected the {model} {kind} step before TE output', flush=True)


def main():
    hphi = Path(sys.argv[1]).resolve()
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = mpi_size(launcher)
    assert ranks is not None
    root = Path('symmetry_spinflip_evolution')
    if root.exists():
        shutil.rmtree(root)
    root.mkdir()
    check_drive_rows()
    for model, case in CASES.items():
        aggregate = model == 'Spin'
        run_thermal(root, hphi, launcher, case, 'tpq', False, aggregate)
        run_thermal(root, hphi, launcher, case, 'ctpq', False, aggregate)
        run_thermal(root, hphi, launcher, case, 'ctpq', True, aggregate)
    layout = os.environ.get('HPHI_SYMMETRY_BASIS_LAYOUT', 'distributed')
    env = layout_env(layout)
    if ranks == 1:
        preflight_negatives(root, hphi, env)
    for model, case in CASES.items():
        seed = seed_cg(root, model+'_seed', case, hphi, launcher, env)
        for drive in ('static', *DRIVES[model]):
            result = evolve(root, f'{model}_{drive}', case, seed, 'zvo_eigenvec_0', drive,
                            hphi, launcher, env, interval=2 if drive == 'twobody' else 1)
            if drive == DRIVES[model][0]:
                te_to_cg(root, case, result, hphi, launcher, env)
        if model == 'SpinGC':
            evolve(root, f'{model}_quench', case, seed, 'zvo_eigenvec_0', 'quench',
                   hphi, launcher, env)
        other = 'replicated' if layout == 'distributed' else 'distributed'
        evolve(root, f'{model}_{DRIVES[model][0]}_{other}', case, seed, 'zvo_eigenvec_0',
               DRIVES[model][0], hphi, launcher, layout_env(other))
        fd = seed_fulldiag(root, model+'_fd_seed', case, hphi, env)
        evolve(root, f'{model}_fd_static', case, fd, 'fdseed', 'static', hphi, launcher, env)
    print('maximum errors:', {k: float(v) for k, v in ERRORS.items()})
    print('Spin-flip TPQ/cTPQ/TE passed')


if __name__ == '__main__':
    main()
