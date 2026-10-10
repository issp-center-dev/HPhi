#!/usr/bin/env python3
"""Kondo all-times guards and independent right-endpoint TE trajectories."""
import os
from pathlib import Path
import shlex
import shutil
import struct
import sys
import tempfile

import numpy as np

from symmetry_kondo_common import (expert_case, drive_rows, read_kondo_checkpoint,
                                   read_manifest, read_vector_parts)
from symmetry_kondo_reference import (Case, make_reference, sector_hamiltonian,
                                      drive_hamiltonian, _drive_action)
from symmetry_kondo_checkpoint import joined
from symmetry_kondo_observables import (add_correlation_requests,
                                        correlation_requests, CORRELATION_FILES)
from symmetry_kondo_thermal import (finite_data, thermal_correlation_reference,
                                    check_correlation_block)
from symmetry_spingc_common import (definition, checkpoint_digest, mpi_size, run_case)
from symmetry_spingc_te import checkpoints, schedule, TIMES, AMPLITUDES


def preflight_negatives(root, executable, env):
    case = Case('Kondo', 3, 'block', 2, 1, 1)
    seed = expert_case(root/'seed', case, method=3,
                       options={'CalcMod': {'OutputEigenVec': 1}})
    run_case(seed, executable, 'seed', [], env)
    # Removing the composed-view input-form guard must admit this forbidden
    # local hopping, whose physical single-occupancy projection is zero.
    for name, rows, diagnostic in [
        ('local_hopping', [(1, 0, 0, 0, .3, 0), (0, 0, 1, 0, .3, 0)],
         'requires onsite local bilinears'),
        ('nonhermitian', [(3, 0, 4, 0, .3, 0)], 'NonHermite'),
        ('bad_site', [(6, 0, 6, 0, .3, 0)], 'invalid site/spin'),
        ('noninvariant', [(0, 0, 0, 0, .3, 0)], 'preflight failed at step 3'),
        ('nonconserving', [(j, a, j, b, .3, 0) for j in range(3)
                           for a, b in ((0, 1), (1, 0))], 'conserve')]:
        path = expert_case(root/name, case, method=4,
                           options={'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1},
                                    'ModPara': {'Lanczos_max': 5, 'ExpandCoef': 8,
                                                'ExpecInterval': 1, 'OutputInterval': 1}})
        (path/'output').mkdir(exist_ok=True)
        shutil.copy2(seed/'output/zvo_eigenvec_0_rank_0.dat',
                     path/'output/seed_rank_0.dat')
        with (path/'sym.def').open('a') as handle:
            handle.write('SpectrumVec seed\n')
        # Step 3 shares time .025 with step 2: dt=0 must not bypass validation.
        schedule(path, 'TEOneBody', rows, amplitudes=[0, 0, 0, 1, 0])
        run_case(path, executable, 'reject', [], env, diagnostic)
        assert not (path/'output/zvo_Time_TE_Step.dat').exists()
        assert not (path/'output/zvo_SS.dat').exists()
        assert not list((path/'output').glob('zvo_eigenvec*'))
        assert not list(path.glob('sector_initial*'))
        print('rejected late '+name+' before TE state/SS output', flush=True)


ERRORS = dict(vector=0., data=0., endpoint=0., norm=0., commutator=0., residual=0.)


def te_step(h_right, vector, dt):
    term, work = vector.copy(), vector.copy()
    for order in range(1, 9):
        term = (-1j*dt/order)*(h_right@term)
        work += term
    return work, np.linalg.norm(work)


def exact_endpoint_step(h_right, vector, dt):
    energies, rotation = np.linalg.eigh(h_right)
    return rotation@(np.exp(-1j*dt*energies)*(rotation.conj().T@vector))


def attach_seed(path, source, label):
    (path/'output').mkdir(exist_ok=True)
    for rank, file in enumerate(checkpoints(source, label)):
        shutil.copy2(file, path/'output'/f'seed_rank_{rank}.dat')
    with (path/'sym.def').open('a') as handle:
        handle.write('SpectrumVec seed\n')


def laser_definition(path):
    values = [1, .07, .8, 1, 0, 3, 2, 1, 0]
    (path/'laser.def').write_text('====\nNLaser 9\n====\n====\n====\n'+
                                ''.join(f'p{i} {v}\n' for i, v in enumerate(values)))
    with (path/'sym.def').open('a') as handle:
        handle.write('Laser laser.def\n')


def check_identity(source, target):
    # Source solver provenance, H, and payload checksum may change; physical
    # space, group, ordering and rank ownership must agree exactly.
    np.testing.assert_array_equal(source[:21], target[:21])
    np.testing.assert_array_equal(source[28:], target[28:])


def capture_matches_import(path, source, label, dimension, ranks):
    stem = 'sector_initial_sample0_step0.rank'
    files = sorted(path.glob(stem+'*.dat'))
    assert len(files) == ranks, ('missing immediate TE import capture', files)
    actual = read_vector_parts(path, stem, dimension)
    expected = joined(source, label)
    np.testing.assert_array_equal(actual, expected)
    infos = [read_manifest(file.with_suffix('.info')) for file in files]
    assert sum(int(info['local_dim']) for info in infos) == dimension
    assert all(int(info['dim']) == dimension for info in infos)


def evolve(root, name, case, reference, source, source_label, drive,
           executable, launcher, env, interval=1, times=TIMES):
    path = expert_case(root/name, case, method=4,
                       options={'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1},
                                'ModPara': {'Lanczos_max': 5, 'ExpandCoef': 8,
                                            'ExpecInterval': interval,
                                            'OutputInterval': interval,
                                            'Tinit': 0, 'TimeSlice': .01}})
    if drive == 'quench':
        file = path/'Trans.def'
        static_rows = [tuple(map(float, line.split())) for line in file.read_text().splitlines()[5:]]
        static_rows = [tuple(int(x) for x in row[:4])+tuple(row[4:]) for row in static_rows]
        definition(path, 'Trans.def', static_rows+drive_rows(case, 'onebody'))
    requests = correlation_requests() if case.layout == 'block' else None
    if requests:
        add_correlation_requests(path, requests, False)
    attach_seed(path, source, source_label)
    if drive == 'laser':
        times = np.arange(5)*.01
        laser_definition(path)
    elif drive == 'explicit_laser':
        times = np.arange(5)*.01
        # Each row is a hopping increment, so the static twist is retained.
        text = '====\nNTimeSteps 5\n====\n====\n====\n'
        for time in times:
            rows = drive_rows(case, 'laser', .07*np.sin(.8*time))
            text += f'{time} {len(rows)}\n'
            text += ''.join(' '.join(map(str, row))+'\n' for row in rows)
        (path/'TEOneBody.def').write_text(text)
        with (path/'sym.def').open('a') as handle:
            handle.write('TEOneBody TEOneBody.def\n')
    else:
        schedule(path, 'TETwoBody' if drive == 'twobody' else 'TEOneBody',
                 drive_rows(case, drive) if drive in ('onebody', 'twobody') else [], times)
    text = run_case(path, executable, 'evolve', launcher, env)
    ranks = mpi_size(launcher)
    dimension = len(reference.representatives)
    capture_matches_import(path, source, source_label, dimension, ranks)
    metadata = read_manifest(path)
    source_files = checkpoints(source, source_label)
    source_header, _ = read_kondo_checkpoint(source_files[0])
    assert int(metadata['source_method']) == source_header[22]
    assert int(metadata['source_step']) == source_header[24]
    assert float(metadata['source_time']) == struct.unpack('<d', struct.pack('<Q', source_header[25]))[0]
    assert int(metadata['source_hamiltonian_digest'], 16) == source_header[21]
    dynamic = drive in ('onebody', 'twobody', 'laser', 'explicit_laser')
    assert metadata['te_hamiltonian'] == ('time_dependent' if dynamic else 'static')
    if drive == 'quench':
        assert 'hamiltonian_changed=yes' in text
    if dynamic:
        assert 'all 5 Hamiltonians passed preflight' in text
    initial = joined(source, source_label)
    vector, exact = initial.copy(), initial.copy()
    base = sector_hamiltonian(reference)
    ss, norms, flct = [], [], []
    saved_vectors = []
    for step, time in enumerate(times):
        kind = ('laser' if drive == 'explicit_laser' else
                'onebody' if drive == 'quench' else drive)
        amplitude = (1. if drive == 'quench' else .07*np.sin(.8*time)
                     if kind == 'laser' else AMPLITUDES[step])
        h = base if drive == 'static' else drive_hamiltonian(case, kind, amplitude)
        np.testing.assert_allclose(h, h.conj().T, atol=1e-12, rtol=0)
        if kind == 'laser':
            # For Hermitian physical H, zero (I-P)HB is equivalent to [H,P]=0.
            # Test every sector column, not merely its projected matrix B†HB.
            action = _drive_action(case, 'laser', amplitude)
            hb = np.column_stack([action(column) for column in reference.basis.T])
            error = np.linalg.norm(hb-reference.basis@h)
            ERRORS['commutator'] = max(ERRORS['commutator'], error)
            assert error < 1e-12
        dt = time-(times[step-1] if step else times[0])
        assert dt*np.linalg.norm(h, 2) <= .25
        work, norm = te_step(h, vector, dt)
        vector = work/norm  # TEM normalizes after the entire Taylor polynomial.
        exact = exact_endpoint_step(h, exact, dt)
        assert np.isfinite(vector).all() and np.isfinite(norm)
        ERRORS['endpoint'] = max(ERRORS['endpoint'], np.linalg.norm(vector-exact))
        ERRORS['norm'] = max(ERRORS['norm'], abs(norm-1))
        assert np.linalg.norm(vector-exact) <= 1e-8 and abs(norm-1) <= 1e-8
        saved_vectors.append(vector.copy())
        hv = h@vector
        weights = abs(reference.basis@vector)**2
        m = {key: weights@reference.moments[key] for key in ('N', 'N2', 'D', 'D2', 'Sz', 'Sz2')}
        ss.append([time, np.vdot(vector, hv).real, np.vdot(hv, hv).real, m['D'], m['N'], step])
        norms.append([time, norm, step])
        flct.append([time, m['N'], m['N2'], m['D'], m['D2'], m['Sz'], m['Sz2'], step])
        files = checkpoints(path, f'zvo_eigenvec_{step}')
        if step % interval:
            assert not files
        else:
            assert len(files) == ranks
            actual = joined(path, f'zvo_eigenvec_{step}')
            error = np.linalg.norm(actual-vector)
            ERRORS['vector'] = max(ERRORS['vector'], error)
            assert error <= 3e-11, (name, step, error)
            for old_file, file in zip(source_files, files):
                old, _ = read_kondo_checkpoint(old_file)
                header, _ = read_kondo_checkpoint(file)
                check_identity(old, header)
                assert tuple(header[22:25]) == (4, 0, step)
                assert struct.unpack('<d', struct.pack('<Q', header[25]))[0] == time
                digest = (int(metadata[f'te_hamiltonian_{step}'], 16) if dynamic else
                          int(metadata['hamiltonian_digest'].split(':')[1], 16))
                assert header[21] == digest
        if requests:
            values = thermal_correlation_reference(case, reference, vector, requests)
            for kind, rows in requests.items():
                file = path/'output'/f'zvo_{CORRELATION_FILES[kind]}_step{step}.dat'
                if step % interval:
                    assert not file.exists()
                else:
                    check_correlation_block(finite_data(file), rows, values[kind], 0, str(file))
    for family, expected in [('SS', ss), ('Norm', norms), ('Flct', flct)]:
        actual = finite_data(path/'output'/f'zvo_{family}.dat')
        expected = np.asarray(expected)
        np.testing.assert_array_equal(actual[:, -1], expected[:, -1])
        error = np.max(abs(actual-expected))
        ERRORS['data'] = max(ERRORS['data'], error)
        assert error <= 4e-11, (name, family, error)
    finals = checkpoints(path, 'zvo_eigenvec_final')
    assert all(final.read_bytes() == last.read_bytes()
               for final, last in zip(finals, checkpoints(path, 'zvo_eigenvec_4')))
    parts = [read_kondo_checkpoint(file) for file in finals]
    assert len(parts) == ranks and all(tuple(h[26:28]) == checkpoint_digest(parts) for h, _ in parts)
    error = np.linalg.norm(joined(path, 'zvo_eigenvec_final')-vector)
    ERRORS['vector'] = max(ERRORS['vector'], error)
    assert error <= 3e-11
    if drive in ('onebody', 'twobody'):
        assert metadata['te_hamiltonian_0'] == metadata['te_hamiltonian_2']
        assert len({metadata[f'te_hamiltonian_{i}'] for i in range(5)}) == 4
        np.testing.assert_allclose(saved_vectors[2], saved_vectors[3], atol=5e-15, rtol=0)
    print(name+' trajectory/import/outputs passed', flush=True)
    return path


def prepare(root, name, case, reference, executable, launcher, env):
    path = expert_case(root/name, case, method=3,
                       options={'CalcMod': {'OutputEigenVec': 1},
                                'ModPara': {'LanczosEps': 18}})
    run_case(path, executable, 'cg_seed', launcher, env)
    vector = joined(path)
    h = sector_hamiltonian(reference)
    energy = np.vdot(vector, h@vector).real
    residual = np.linalg.norm(h@vector-energy*vector)
    ERRORS['residual'] = max(ERRORS['residual'], residual)
    assert residual < 2e-8, (name, residual)
    return path


def imports(root, case, reference, source, executable, launcher, env):
    # A new initial time and grid demonstrate state import, not solver restart.
    target = evolve(root, source.name+'_reimport', case, reference, source,
                    'zvo_eigenvec_final', 'twobody', executable, launcher, env,
                    times=np.array([.1, .103, .109, .109, .12]))
    meta = read_manifest(target)
    assert meta['source_step'] == '4' and meta['source_method'] == '4'
    assert float(meta['te_time_0']) == .1
    path = expert_case(root/(source.name+'_to_cg'), case, method=3,
                       options={'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1}})
    (path/'output').mkdir(exist_ok=True)
    for rank, file in enumerate(checkpoints(source, 'zvo_eigenvec_final')):
        shutil.copy2(file, path/'output'/f'zvo_eigenvec_0_rank_{rank}.dat')
    run_case(path, executable, 'te_to_cg', launcher, env)
    np.testing.assert_array_equal(joined(path), joined(source, 'zvo_eigenvec_final'))
    for old, new in zip(checkpoints(source, 'zvo_eigenvec_final'), checkpoints(path, 'zvo_eigenvec_0')):
        check_identity(read_kondo_checkpoint(old)[0], read_kondo_checkpoint(new)[0])
    print(path.name+' immediate TE→CG payload passed', flush=True)


def reject_laser_geometry(root, executable, env):
    case = Case('KondoNConserved', 3, 'alternating', 2, None, 1)
    path = expert_case(root/'reject_laser_geometry', case, method=4,
                       options={'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1},
                                'ModPara': {'Lanczos_max': 5, 'ExpandCoef': 8,
                                            'ExpecInterval': 1, 'OutputInterval': 1,
                                            'Tinit': 0, 'TimeSlice': .01}})
    laser_definition(path)
    # Lx=3 yields the same displacement on all alternating bonds. Lx=6
    # instead gives -4,-4,2 after the existing wrap rule: 
    # at t=0 it agrees, but t=.01 breaks translation.
    laser = path/'laser.def'
    laser.write_text(laser.read_text().replace('p5 3\n', 'p5 6\n'))
    run_case(path, executable, 'reject', [], env, 'preflight failed at step 1')
    assert not (path/'output/zvo_Time_TE_Step.dat').exists()
    assert not (path/'output/zvo_SS.dat').exists()
    assert not list((path/'output').glob('zvo_eigenvec*'))
    print('rejected late noninvariant Laser geometry before output', flush=True)


def main():
    executable = Path(sys.argv[1]).resolve()
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = mpi_size(launcher)
    assert ranks in (1, 2, 4, 16)
    root = Path(tempfile.mkdtemp(prefix='symmetry_kondo_te_', dir='.')).resolve()
    print('artifacts:', root, flush=True)
    scope = os.environ.get('HPHI_TEST_KONDO_TE_SCOPE', 'all')
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='distributed', HPHI_TEST_SYMMETRY_CAPTURE='1')
    if ranks == 1 and scope in ('all', 'preflight'):
        preflight_negatives(root, executable, env)
        reject_laser_geometry(root, executable, env)
    if scope == 'preflight':
        return
    for layout in ('distributed', 'replicated'):
        env['HPHI_SYMMETRY_BASIS_LAYOUT'] = layout
        for model in ('Kondo', 'KondoNConserved', 'KondoGC'):
            for momentum in (0, 1):
                case = Case(model, 3, 'block', None if model == 'KondoGC' else 2,
                            1 if model == 'Kondo' else None, momentum)
                reference = make_reference(case)
                prefix = f'{model}_k{momentum}_{layout}'
                source = prepare(root, prefix+'_seed', case, reference, executable, launcher, env)
                for drive in ('static', 'quench', 'onebody', 'twobody', 'laser'):
                    result = evolve(root, prefix+'_'+drive, case, reference, source,
                                    'zvo_eigenvec_0', drive, executable, launcher, env,
                                    interval=momentum+1)
                    if scope == 'capture':
                        return
                    if drive == 'onebody' and momentum == 1 and layout == 'distributed':
                        imports(root, case, reference, result, executable, launcher, env)
    env['HPHI_SYMMETRY_BASIS_LAYOUT'] = 'distributed'
    for model in ('Kondo', 'KondoNConserved', 'KondoGC'):
        case = Case(model, 3, 'alternating', None if model == 'KondoGC' else 2,
                    1 if model == 'Kondo' else None, 1)
        reference = make_reference(case)
        prefix = model+'_alternating'
        source = prepare(root, prefix+'_seed', case, reference, executable, launcher, env)
        evolve(root, prefix+'_explicit', case, reference, source, 'zvo_eigenvec_0',
               'explicit_laser', executable, launcher, env)
    print('maximum errors:', {key: float(value) for key, value in ERRORS.items()}, flush=True)


if __name__ == '__main__':
    main()
