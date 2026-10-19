#!/usr/bin/env python3
"""Fixed-N spin-mixing sectors: TPQ/cTPQ samples and TE trajectories vs explicit updates."""
from pathlib import Path
import os
import shlex
import shutil
import sys
import numpy as np
from symmetry_kondo_common import read_manifest, read_vector_parts
from symmetry_kondo_observables import CORRELATION_FILES, add_correlation_requests
from symmetry_kondo_te import attach_seed, exact_endpoint_step, te_step
from symmetry_kondo_thermal import (EXPLICIT_BETAS, EXPLICIT_FLAGS, NSTEPS, NUM_AVE,
                                    UNIFORM_BETAS, captured_initials, check_correlation_block,
                                    check_data, check_metadata, ctpq_step, finite_data, tpq_step)
from symmetry_spingc_common import (checkpoint_digest, join_rank_vectors, mpi_size,
                                    read_checkpoint, run_case)
from symmetry_spingc_te import AMPLITUDES, TIMES, checkpoints, schedule
from symmetry_nconserved_common import drive_rows, expert_case, hamiltonian_families
from symmetry_nconserved_core import (REFERENCE_TOL, correlation_operators, correlation_requests,
                                      expectations)
from symmetry_nconserved_reference import (Case, DEFAULT, MODELS, drive_matrix, make_reference,
                                           sector_matrix)

MODEL_ID = {'HubbardNConserved': 6, 'tJNConserved': 11}
MOMENTS = ('N', 'N2', 'D', 'D2', 'Sz', 'Sz2')
ERRORS = dict(vector=0., endpoint=0., data=0.)


def layout_env(layout, **extra):
    return dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout, HPHI_TEST_SYMMETRY_CAPTURE='1',
                **extra)


def moments(ref, vector):
    weights = abs(ref.basis @ vector)**2
    return {key: float(weights @ ref.moments[key]) for key in MOMENTS}


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


def check_thermal_correlations(path, aggregate, emitted, history, ref, operators, requests):
    for kind, rows in requests.items():
        stem = CORRELATION_FILES[kind]
        expected = [[expectations(ref, v, operators)[kind] for v in item['vectors']]
                    for item in history]
        if aggregate:
            data = finite_data(path/'output'/f'zvo_{stem}_tpq.dat')
            prefixes = [(s, step) for s in range(NUM_AVE) for step in emitted for _ in rows]
            np.testing.assert_array_equal(data[:, :2], prefixes)
            offset = 0
            for sample in range(NUM_AVE):
                for step in emitted:
                    check_correlation_block(data[offset:offset+len(rows)], rows,
                                            expected[sample][step], 2, f'{kind} {sample} {step}',
                                            atol=REFERENCE_TOL)
                    offset += len(rows)
        else:
            for sample in range(NUM_AVE):
                for step in range(NSTEPS):
                    file = path/'output'/f'zvo_{stem}_set{sample}step{step}.dat'
                    if step not in emitted:
                        assert not file.exists(), file
                    else:
                        check_correlation_block(finite_data(file), rows, expected[sample][step],
                                                0, f'{kind} {sample} {step}', atol=REFERENCE_TOL)


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
    requests = correlation_requests()
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
    vectors, prenorms = captured_initials(path, len(h), ref.representatives, mpi_size(launcher))
    history = [thermal_history(kind, h, v, p, ref, case.sites, betas)
               for v, p in zip(vectors, prenorms)]
    check_metadata(path, case, kind, mpi_size(launcher), explicit, len(h))
    check_data(path, aggregate, history)
    check_thermal_correlations(path, aggregate, emitted, history, ref,
                               correlation_operators(case, ref, requests), requests)
    # Sz fluctuates: Flct reports the state value, not a fixed input quantum number.
    flct = np.vstack([item['Flct'] for item in history])
    assert np.max(flct[:, 6]-flct[:, 5]**2) > 1e-4
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


def te_hamiltonians(case, ref, drive):
    base = sector_matrix(ref)
    if drive == 'static':
        return [base]*len(TIMES)
    if drive == 'quench':
        return [base+drive_matrix(ref, case, 'onebody')]*len(TIMES)
    return [base+drive_matrix(ref, case, drive, amplitude) for amplitude in AMPLITUDES]


def evolve(root, name, case, source, label, drive, hphi, launcher, env, interval=1):
    ref = make_reference(case)
    families = None
    if drive == 'quench':
        families = hamiltonian_families(case.sites)
        families['Trans'] = families['Trans'] + drive_rows(case.sites, 'onebody')
    path = expert_case(root/name, case, method=4, options={
        'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1},
        'ModPara': {'Lanczos_max': len(TIMES), 'ExpandCoef': 8, 'ExpecInterval': interval,
                    'OutputInterval': interval},
        **({'families': families} if families else {})})
    requests = correlation_requests()
    add_correlation_requests(path, requests, False)
    attach_seed(path, source, label)
    family = 'TETwoBody' if drive == 'twobody' else 'TEOneBody'
    rows = drive_rows(case.sites, drive) if drive in ('onebody', 'twobody') else []
    schedule(path, family, rows)
    text = run_case(path, hphi, 'evolve', launcher, env)
    ranks = mpi_size(launcher)
    dim = ref.basis.shape[1]
    initial = join_rank_vectors(checkpoints(source, label), dim)
    np.testing.assert_array_equal(read_vector_parts(path, 'sector_initial_sample0_step0.rank', dim),
                                  initial)
    if len(checkpoints(source, label)) != ranks:
        assert 'cross-rank v1/v2 read: payload checksum verified' in text, text
    meta = read_manifest(path)
    dynamic = drive in ('onebody', 'twobody')
    assert meta['te_hamiltonian'] == ('time_dependent' if dynamic else 'static')
    if dynamic:
        assert f'all {len(TIMES)} Hamiltonians passed preflight' in text
    if drive == 'quench':
        assert 'hamiltonian_changed=yes' in text
    operators = correlation_operators(case, ref, requests)
    vector, exact = initial.copy(), initial.copy()
    ss, norms, flct = [], [], []
    for step, (time, h) in enumerate(zip(TIMES, te_hamiltonians(case, ref, drive))):
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
        error = np.linalg.norm(join_rank_vectors(files, dim)-vector)
        ERRORS['vector'] = max(ERRORS['vector'], error)
        assert error <= 3e-11, (name, step, error)
        for file in files:
            header, _ = read_checkpoint(file)
            assert tuple(header[4:9]) == (MODEL_ID[case.model], case.sites, 0, 0, case.n)
            assert tuple(header[22:25]) == (4, 0, step)
        values = expectations(ref, vector, operators)
        for kind, rows_ in requests.items():
            check_correlation_block(finite_data(path/'output'/f'zvo_{CORRELATION_FILES[kind]}_step{step}.dat'),
                                    rows_, values[kind], 0, f'{name} {kind} step{step}',
                                    atol=REFERENCE_TOL)
    for family_name, expected in (('SS', ss), ('Norm', norms), ('Flct', flct)):
        actual = finite_data(path/'output'/f'zvo_{family_name}.dat')
        expected = np.asarray(expected)
        np.testing.assert_array_equal(actual[:, -1], expected[:, -1])
        error = np.max(abs(actual-expected))
        ERRORS['data'] = max(ERRORS['data'], error)
        assert error <= 4e-11, (name, family_name, error)
    # N is conserved by every drive while Sz is not.
    flct = np.asarray(flct)
    assert np.max(abs(flct[:, 1]-case.n)) < 1e-12
    if drive in ('onebody', 'twobody', 'quench'):
        assert np.ptp(flct[:, 5]) > 1e-6, (name, flct[:, 5])
    finals = checkpoints(path, 'zvo_eigenvec_final')
    parts = [read_checkpoint(file) for file in finals]
    assert len(parts) == ranks and all(tuple(h[26:28]) == checkpoint_digest(parts) for h, _ in parts)
    assert np.linalg.norm(join_rank_vectors(finals, dim)-vector) <= 3e-11
    print(f'{name} trajectory/outputs passed', flush=True)
    return path


def te_to_cg(root, case, source, hphi, launcher, env):
    """The final TE state is a CG input of the same space; the payload is unchanged."""
    path = expert_case(root/(source.name+'_to_cg'), case, method=3,
                       options={'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1}})
    (path/'output').mkdir()
    for rank, file in enumerate(checkpoints(source, 'zvo_eigenvec_final')):
        shutil.copy2(file, path/'output'/f'zvo_eigenvec_0_rank_{rank}.dat')
    run_case(path, hphi, 'te_to_cg', launcher, env)
    dim = make_reference(case).basis.shape[1]
    np.testing.assert_array_equal(join_rank_vectors(checkpoints(path, 'zvo_eigenvec_0'), dim),
                                  join_rank_vectors(checkpoints(source, 'zvo_eigenvec_final'), dim))
    for file in checkpoints(path, 'zvo_eigenvec_0'):
        assert read_checkpoint(file)[0][22] == 3


def preflight_negatives(root, hphi, env):
    """Every TE step is checked before evolution, on the shortest serial path."""
    case = Case('HubbardNConserved', 4, 3, 1)
    seed = seed_cg(root, 'negative_seed', case, hphi, [], env)
    local = [(0, 0, 0, 1, .3, 0.), (0, 1, 0, 0, .3, 0.)]
    canonical = dict(DEFAULT, soc=0, soc2=0, hx=0, w=0)
    for name, rows, message, options, parameters in [
            ('noninvariant', local, 'preflight failed at step 3', {}, None),
            ('fixed_sz_spin_mixing', drive_rows(4, 'onebody'), 'conserve',
             {'ModPara': {'2Sz': 1}}, canonical)]:
        source = seed
        if parameters is not None:
            source = expert_case(root/(name+'_seed'), case, method=3, parameters=parameters,
                                 options={'CalcMod': {'OutputEigenVec': 1}, **options})
            run_case(source, hphi, 'cg_seed', [], env)
        path = expert_case(root/('reject_'+name), case, method=4, parameters=parameters, options={
            'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1},
            'ModPara': {'Lanczos_max': len(TIMES), 'ExpandCoef': 8, 'ExpecInterval': 1,
                        'OutputInterval': 1, **options.get('ModPara', {})}})
        attach_seed(path, source, 'zvo_eigenvec_0')
        amplitudes = [0, 0, 0, 1, 0] if name == 'noninvariant' else [0, 1, 0, 1, 0]
        schedule(path, 'TEOneBody', rows, amplitudes=amplitudes)
        run_case(path, hphi, 'reject', [], env, message)
        assert not (path/'output/zvo_SS.dat').exists()
        assert not list((path/'output').glob('zvo_eigenvec_[0-9]*'))
        print('rejected', name, 'before TE output', flush=True)


def main():
    hphi = Path(sys.argv[1]).resolve()
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = mpi_size(launcher)
    assert ranks is not None
    root = Path('symmetry_nconserved_evolution')
    if root.exists():
        shutil.rmtree(root)
    root.mkdir()
    for model in MODELS:
        case = Case(model, 4, 3, 1)
        aggregate = model == 'tJNConserved'
        run_thermal(root, hphi, launcher, case, 'tpq', False, aggregate)
        run_thermal(root, hphi, launcher, case, 'ctpq', False, aggregate)
        run_thermal(root, hphi, launcher, case, 'ctpq', True, aggregate)
    layout = os.environ.get('HPHI_SYMMETRY_BASIS_LAYOUT', 'distributed')
    env = layout_env(layout)
    if ranks == 1:
        preflight_negatives(root, hphi, env)
    for model in MODELS:
        case = Case(model, 4, 3, 1)
        seed = seed_cg(root, model+'_seed', case, hphi, launcher, env)
        for drive in ('static', 'quench', 'onebody', 'twobody'):
            result = evolve(root, f'{model}_{drive}', case, seed, 'zvo_eigenvec_0', drive,
                            hphi, launcher, env, interval=2 if drive == 'twobody' else 1)
            if drive == 'onebody':
                te_to_cg(root, case, result, hphi, launcher, env)
        other = 'replicated' if layout == 'distributed' else 'distributed'
        evolve(root, f'{model}_onebody_{other}', case, seed, 'zvo_eigenvec_0', 'onebody',
               hphi, launcher, layout_env(other))
        fd = seed_fulldiag(root, model+'_fd_seed', case, hphi, env)
        evolve(root, f'{model}_fd_static', case, fd, 'fdseed', 'static', hphi, launcher, env)
    print('maximum errors:', {k: float(v) for k, v in ERRORS.items()})
    print('NConserved TPQ/cTPQ/TE passed')


if __name__ == '__main__':
    main()
