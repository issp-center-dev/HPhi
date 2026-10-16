"""CG no-iteration evaluation uses current convergence, not saved FLAGS alone."""
import shutil
import struct

import numpy as np
import restart_common as rc

ROOT = rc.fresh_dir('restart_cg_noop')
STAN = ('L=8\nmodel="Spin"\nmethod="CG"\nlattice="chain"\nJ=1\n2Sz=0\n'
        'exct=2\nLanczos_max=80\nLanczosEps=10\noutputmode="None"\n')


def clone(source, name):
    path = ROOT / name
    shutil.copytree(source, path)
    return path


def energies(path):
    return np.array([float(line.split()[1]) for line in
                     (path / 'output/zvo_energy.dat').read_text().splitlines()
                     if line.lstrip().startswith('Energy')])


def history(path):
    return (path / 'output/zvo_Lanczos_Step.dat').read_bytes()


def checkpoint(path):
    return {p.name: p.read_bytes() for p in (path / 'output').glob('zvo_restart_*')}


def output_vectors(path, step, sector):
    files = list((path / 'output').glob('zvo_eigenvec_*_rank_*.dat'))
    ranks = rc.manifest(path)[1][0]['ranks']
    assert len(files) == 2*ranks
    for f in files:
        data = f.read_bytes()
        if sector:
            assert struct.unpack_from('<Q', data, 8)[0] in (1, 2)
            assert struct.unpack_from('<Q', data, 24*8)[0] == step
        else:
            assert struct.unpack_from('=i', data)[0] == step
            dim = struct.unpack_from('=Q', data, 4)[0]
            assert len(data) == 12 + 16*(dim+1)
    return {f.name: f.read_bytes() for f in files}


def exercise(template, tag, calc, mod, env=None, sector=False):
    def configure(path, mode, limit, eps):
        rc.set_key(path, calc, 'ReStart', mode)
        rc.set_key(path, mod, 'Lanczos_max', limit)
        rc.set_key(path, mod, 'LanczosEps', eps)

    saved = clone(template, tag+'_saved')
    configure(saved, 1, 80, 10)
    rc.run(saved, 'converged', env=env)
    row = rc.manifest(saved)[1][0]
    step = row['step']
    header, _ = rc.v3_header(rc.restart_files(saved, 'zvo', 0, 0, row['gen'])[0])
    assert header[39] == 1 and 2 < step < 80
    values = np.loadtxt(saved / 'output/zvo_Lanczos_Step.dat', skiprows=1, ndmin=2)
    # Establish that the saved residual really fails the new threshold.
    tight_threshold = 1e-7 * max(1, np.max(np.abs(values[-1, 3:])))
    assert values[-1, 1] > tight_threshold
    full_history = history(saved)
    canonical = b''.join(line for line in full_history.splitlines(True)
                         if line.lstrip().startswith(b'Step') or int(line.split()[0]) < step)
    saved_energy = energies(saved)

    def check_noop(path, saved_step, expected_history, currently_converged, expected_energy=None):
        rc.set_key(path, calc, 'OutputEigenVec', 1)
        saved_checkpoint = checkpoint(path)
        with (path / 'output/zvo_Lanczos_Step.dat').open('ab') as stream:
            stream.write(('{} 0 0 0 0\npartial'.format(saved_step+1)).encode())
        first_energy = first_vectors = None
        for repeat in (0, 1):
            (path / 'output/zvo_energy.dat').write_text('stale sentinel\n')
            text = rc.run(path, 'noop_{}'.format(repeat), env=env)
            messages = [line for line in text.splitlines() if 'LOBPCG: no remaining steps' in line]
            assert len(messages) == 1, text
            assert ('saved step {},'.format(saved_step)) in messages[0]
            assert (', converged' in messages[0]) == currently_converged, messages[0]
            assert ('convergence not confirmed at current tolerance' in messages[0]) != currently_converged
            assert history(path) == expected_history
            assert checkpoint(path) == saved_checkpoint
            result = energies(path)
            assert len(result) == 2 and np.isfinite(result).all()
            if expected_energy is not None:
                np.testing.assert_allclose(result, expected_energy, rtol=0, atol=1e-10)
            vectors = output_vectors(path, saved_step, sector)
            if repeat == 0:
                first_energy, first_vectors = result, vectors
            else:
                np.testing.assert_allclose(result, first_energy, rtol=0, atol=1e-10)
                assert vectors == first_vectors
        # Independent evaluation of exactly the vectors produced by the caller.
        rc.set_key(path, calc, 'ReStart', 0)
        rc.set_key(path, calc, 'OutputEigenVec', 0)
        rc.set_key(path, calc, 'InputEigenVec', 1)
        rc.run(path, 'import_evaluation', env=env)
        np.testing.assert_allclose(energies(path), first_energy, rtol=0, atol=1e-10)
        assert history(path) == expected_history and checkpoint(path) == saved_checkpoint

    for mode in (2, 3):
        # I1: no room remains, but saved FLAGS=1 no longer means convergence.
        path = clone(saved, '{}_tight_low_mode{}'.format(tag, mode))
        configure(path, mode, 1, 14)
        check_noop(path, step, canonical, False, saved_energy)
        # Ordinary converged no-op must still retain the converged STEP row.
        path = clone(saved, '{}_same_eps_mode{}'.format(tag, mode))
        configure(path, mode, 1, 10)
        check_noop(path, step, full_history, True, saved_energy)

    # The existing FLAGS=0 limit-only no-op remains a successful evaluation.
    limited = clone(template, tag+'_limited')
    configure(limited, 1, 2, 10)
    rc.run(limited, 'limited', env=env)
    row = rc.manifest(limited)[1][0]
    header, _ = rc.v3_header(rc.restart_files(limited, 'zvo', 0, 0, row['gen'])[0])
    assert row['step'] == 3 and header[39] == 0
    for mode in (2, 3):
        path = clone(limited, '{}_limit_mode{}'.format(tag, mode))
        configure(path, mode, 1, 14)
        check_noop(path, 3, history(limited), False)

    # Tightening alone must still execute iterations and reprint saved STEP.
    tight_ref = clone(template, tag+'_tight_ref')
    configure(tight_ref, 0, 80, 14)
    rc.run(tight_ref, 'reference', env=env)
    for mode in (2, 3):
        path = clone(saved, '{}_tight_continue_mode{}'.format(tag, mode))
        configure(path, mode, 80, 14)
        before = checkpoint(path)
        text = rc.run(path, 'continued', env=env)
        assert 'no remaining steps' not in text
        resumed = np.loadtxt(path / 'output/zvo_Lanczos_Step.dat', skiprows=1, ndmin=2)
        assert resumed[-1, 0] > step and history(path).startswith(canonical)
        np.testing.assert_array_equal(resumed[:, 0], np.arange(1, resumed[-1, 0]+1))
        assert resumed[step-1, 2] < values[-1, 2]
        np.testing.assert_allclose(energies(path), energies(tight_ref), rtol=0, atol=1e-9)
        if mode == 3:
            assert checkpoint(path) == before
        else:
            assert rc.manifest(path)[1][0]['step'] == int(resumed[-1, 0])
    print('{}: saved STEP {}; modes 2/3 tight+low, converged/upper-limit no-op, '
          'repeat/import and tighter continuation passed'.format(tag, step), flush=True)


raw = ROOT / 'raw_template'
rc.generate(raw, STAN)
exercise(raw, 'raw', 'calcmod.def', 'modpara.def')

import symmetry_general_terms as sg
sg.ROOT = ROOT / 'sector_templates'
sg.ROOT.mkdir()


def sector_case(path, model, length, momentum, states, raw, projector):
    rc.set_key(path, 'mod.def', 'exct', 2)
    (path / 'namelist.def').write_bytes((path / 'sym.def').read_bytes())
    for layout in ('replicated', 'distributed'):
        exercise(path, 'sector_'+layout, 'calc.def', 'mod.def',
                 env={'HPHI_SYMMETRY_BASIS_LAYOUT': layout}, sector=True)


sg.prepare('Spin', 8, sector_test=sector_case, sector_momenta=(0,))
print('CG current-convergence / exhausted-budget no-op contract ok', flush=True)
