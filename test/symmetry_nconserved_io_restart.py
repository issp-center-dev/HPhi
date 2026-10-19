#!/usr/bin/env python3
"""Fixed-N physical identity of state, matrix and restart files, and their reuse."""
from pathlib import Path
import os
import shlex
import shutil
import struct
import sys
import numpy as np
from symmetry_kondo_common import read_manifest
from symmetry_kondo_observables import check_eigenvector, parse_energy
from symmetry_spingc_checkpoint import checkpoint_files, import_copy
from symmetry_spingc_common import checkpoint_digest, mpi_size, read_checkpoint, run_case
from symmetry_spingc_te import checkpoints
from symmetry_nconserved_common import drive_rows, expert_case
from symmetry_nconserved_reference import (Case, DEFAULT, make_reference, raw_dimension,
                                           sector_matrix)

MODEL_ID = {'HubbardNConserved': 6, 'tJNConserved': 11}
IDENTITY_FAILED = 'checkpoint header / sector / layout validation failed'


def v1_parts(path, case, ranks, layout):
    """Words 4-8 carry the model and N only; Nup/Ndown stay unfixed (0)."""
    parts = [read_checkpoint(file) for file in checkpoint_files(path)]
    assert len(parts) == ranks
    offset = 0
    for rank, (header, vector) in enumerate(parts):
        assert tuple(header[4:9]) == (MODEL_ID[case.model], case.sites, 0, 0, case.n), header
        assert header[9] == raw_dimension(case.model, case.sites, case.n)
        assert tuple(header[11:15]) == (ranks, rank, int(layout == 'distributed'), offset)
        offset += len(vector)
    assert offset == parts[0][0][10]
    assert all(tuple(h[26:28]) == checkpoint_digest(parts) for h, _ in parts)
    return np.concatenate([v for _, v in parts])


def save(root, name, case, hphi, launcher, env, options=None):
    calc = dict({'OutputEigenVec': 1}, **(options or {}))
    path = expert_case(root/name, case, method=3,
                       options={'CalcMod': calc, 'ModPara': {'LanczosEps': 18}})
    run_case(path, hphi, 'save', launcher, env)
    return path


def import_into(root, name, source, target_case, mod=None, parameters=None):
    """Another physical space with the source state files as its input vector."""
    path = expert_case(root/name, target_case, method=3, parameters=parameters,
                       options={'CalcMod': {'InputEigenVec': 1}, 'ModPara': mod or {}})
    (path/'output').mkdir()
    for file in checkpoint_files(source):
        shutil.copy2(file, path/'output'/file.name)
    return path


def check_states(root, hphi, launcher, layout, models=tuple(MODEL_ID)):
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout)
    ranks = mpi_size(launcher)
    for model in models:
        case = Case(model, 4, 3, 1)
        h = sector_matrix(make_reference(case))
        source = save(root, f'{model}_{layout}', case, hphi, launcher, env)
        vector = v1_parts(source, case, ranks, layout)
        check_eigenvector(vector, h, np.linalg.eigvalsh(h)[0])
        assert read_manifest(source)['model'] == model
        same = import_copy(root, f'{model}_{layout}_cg', source)
        with (same/'calc.def').open('a') as f:
            f.write('OutputEigenVec 1\n')
        run_case(same, hphi, 'cg_to_cg', launcher, env)
        assert all(a.read_bytes() == b.read_bytes()
                   for a, b in zip(checkpoint_files(source), checkpoint_files(same)))
        # A quench keeps the space, so a changed H is accepted for the same vector.
        quench = import_copy(root, f'{model}_{layout}_quench', source)
        lines = (quench/'Trans.def').read_text().splitlines()
        for i in range(5, len(lines)):
            row = lines[i].split()
            row[-2:] = [repr(float(value)*1.1) for value in row[-2:]]
            lines[i] = ' '.join(row)
        (quench/'Trans.def').write_text('\n'.join(lines)+'\n')
        with (quench/'calc.def').open('a') as f:
            f.write('OutputEigenVec 1\n')
        assert 'hamiltonian_changed=yes' in run_case(quench, hphi, 'quench', launcher, env)
        np.testing.assert_array_equal(v1_parts(quench, case, ranks, layout), vector)
        print(model, layout, 'v1 identity, CG reuse and quench passed', flush=True)


def check_mismatches(root, hphi, launcher):
    """Equal raw/sector dimensions, different physical space or character."""
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated')
    hubbard1 = Case('HubbardNConserved', 4, 1, 1)
    source = save(root, 'source_hubbard_n1', hubbard1, hphi, launcher, env)
    # At N=1 no site is doubly occupied, so Hubbard and tJ words coincide.
    tj1 = Case('tJNConserved', 4, 1, 1)
    assert np.array_equal(make_reference(tj1).representatives, make_reference(hubbard1).representatives)
    for name, target, detail in [('model', tj1, 'model'),
                                 ('particle_number', Case('HubbardNConserved', 4, 7, 1),
                                  'particle number'),
                                 ('momentum', Case('HubbardNConserved', 4, 1, 3),
                                  'group/character digest')]:
        assert make_reference(target).basis.shape[1] == make_reference(hubbard1).basis.shape[1]
        assert raw_dimension(target.model, 4, target.n) == raw_dimension(hubbard1.model, 4, 1)
        path = import_into(root, 'reject_'+name, source, target)
        text = run_case(path, hphi, 'reject', launcher, env, IDENTITY_FAILED)
        assert detail+' differs' in text, text
        print('rejected', name, 'with equal dimensions', flush=True)
    # The fixed-Sz sector of the same N is a different space, never an implicit import.
    path = import_into(root, 'reject_fixed_sz', source, hubbard1, mod={'2Sz': 1},
                       parameters=dict(DEFAULT, soc=0, soc2=0, hx=0, w=0))
    text = run_case(path, hphi, 'reject', launcher, env, IDENTITY_FAILED)
    assert 'model differs' in text, text


def read_matrix(path):
    lines = path.read_text().splitlines()
    meta = dict(line[len('% HPHI '):].split('=', 1) for line in lines if line.startswith('% HPHI '))
    body = [line for line in lines if not line.startswith('%')]
    rows, cols, count = map(int, body[0].split())
    matrix = np.zeros((rows, cols), complex)
    for line in body[1:1+count]:
        i, j, re, im = line.split()
        matrix[int(i)-1, int(j)-1] = complex(float(re), float(im))
    lower = np.tril(matrix, -1)
    return meta, matrix + lower.conj().T


def check_fulldiag_files(root, hphi):
    """LAPACK FullDiag states, matrix export and import on one rank."""
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated')
    for model in MODEL_ID:
        case = Case(model, 4, 3, 1)
        h = sector_matrix(make_reference(case))
        exporter = expert_case(root/f'{model}_ham', case, method=2,
                               options={'CalcMod': {'OutputHam': 1}})
        run_case(exporter, hphi, 'export', [], env)
        meta, matrix = read_matrix(exporter/'output/zvo_ham_sector.mtx')
        assert meta['model'] == model and meta['nup'] == meta['ndown'] == '0'
        assert int(meta['ne']) == case.n and meta['fixed_flags'] == '0'
        assert int(meta['full_dim']) == raw_dimension(model, 4, 3)
        np.testing.assert_allclose(matrix, h, rtol=0, atol=1e-13)
        importer = expert_case(root/f'{model}_ham_in', case, method=2,
                               options={'CalcMod': {'InputHam': 1, 'OutputEigenVec': 1}})
        (importer/'output').mkdir()
        shutil.copy2(exporter/'output/zvo_ham_sector.mtx', importer/'output')
        run_case(importer, hphi, 'import', [], env)
        energies = np.loadtxt(importer/'output/zvo_energy_sector.dat', ndmin=2)[:, 1]
        np.testing.assert_allclose(energies, np.linalg.eigvalsh(h), rtol=0, atol=1e-10)
        # The exported ground state re-enters CG through the v1 reader.
        consumer = expert_case(root/f'{model}_fd_to_cg', case, method=3,
                               options={'CalcMod': {'InputEigenVec': 1}})
        (consumer/'output').mkdir()
        shutil.copy2(importer/'output/zvo_eigenvec_0_rank_0.dat', consumer/'output')
        run_case(consumer, hphi, 'fd_to_cg', [], env)
        assert abs(parse_energy(consumer)[0]-energies[0]) <= 1e-10
        print(model, 'OutputHam/InputHam and FullDiag to CG passed', flush=True)
    # The matrix identity rejects the tJ space even when the dimension agrees.
    source = expert_case(root/'ham_hubbard_n1', Case('HubbardNConserved', 4, 1, 1), method=2,
                         options={'CalcMod': {'OutputHam': 1}})
    run_case(source, hphi, 'export', [], env)
    target = expert_case(root/'ham_into_tj', Case('tJNConserved', 4, 1, 1), method=2,
                         options={'CalcMod': {'InputHam': 1}})
    (target/'output').mkdir()
    shutil.copy2(source/'output/zvo_ham_sector.mtx', target/'output')
    run_case(target, hphi, 'reject', [], env,
             'metadata model does not match this sector (file=HubbardNConserved, expected=tJNConserved)')


def check_restart_identity(root, hphi, launcher):
    """v3 solver files reuse the 30 identity words with FIXED_FLAGS 0."""
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='distributed')
    ranks = mpi_size(launcher)
    for model in MODEL_ID:
        case = Case(model, 4, 3, 1)
        path = expert_case(root/f'{model}_restart', case, method=3,
                           options={'CalcMod': {'ReStart': 1}})
        run_case(path, hphi, 'restart', launcher, env)
        files = sorted((path/'output').glob('zvo_restart_set0_comp*_gen0_rank_*.dat'))
        assert len(files) == 2*ranks
        for file in files:
            header = struct.unpack('<40Q', file.read_bytes()[:320])
            assert header[1] == 3 and tuple(header[4:9]) == (MODEL_ID[model], 4, 0, 0, 3)
            assert header[29] == 0 and header[30] == 1
        assert 'basis sector' in (path/'output/zvo_restart_manifest.dat').read_text()


# Reader rank counts for each writer: 1<->4 and 4<->16 across the CI jobs.
READERS = {1: (1, 4), 4: (4, 1, 16), 16: (4,)}
TE_TIMES = (0., .03, .09, .09, .13, .2, .25, .3)
TE_AMPLITUDES = (0., .5, -.3, -.3, .2, .4, -.1, .3)


def with_ranks(launcher, ranks):
    """MPIRUN with another process count; None if it cannot be expressed."""
    if ranks == 1:
        return []
    command = list(launcher)
    for i, arg in enumerate(command[:-1]):
        if arg in ('-np', '-n'):
            command[i+1] = str(ranks)
            return command
    return None


def readers(launcher):
    found = []
    for ranks in READERS.get(mpi_size(launcher), (mpi_size(launcher),)):
        command = with_ranks(launcher, ranks)
        if command is None:
            print(f'rank change to {ranks} not exercised: MPIRUN has no -np/-n', flush=True)
        else:
            found.append((ranks, command))
    return found


def set_key(path, name, key, value=None):
    """Replace, append (value given) or drop (value None) a `key value` line."""
    lines = [line for line in (path/name).read_text().splitlines()
             if line.split()[:1] != [key]]
    if value is not None:
        lines.append(f'{key} {value}')
    (path/name).write_text('\n'.join(lines)+'\n')


def history(path, method):
    if method == 3:
        energy = [float(line.split()[1]) for line in
                  (path/'output/zvo_energy.dat').read_text().splitlines()
                  if line.strip().startswith('Energy')]
        return {'energy': np.array(energy),
                'steps': np.loadtxt(path/'output/zvo_Lanczos_Step.dat', skiprows=1, ndmin=2)}
    names = (['zvo_SS.dat', 'zvo_Norm.dat', 'zvo_Flct.dat'] if method == 4 else
             [f'zvo_{family}_rand{s}.dat' for family in ('SS', 'Norm', 'Flct') for s in (0, 1)])
    return {name: np.loadtxt(path/'output'/name, ndmin=2) for name in names}


def same_history(actual, expected, method, label):
    assert actual.keys() == expected.keys()
    for key, value in expected.items():
        assert actual[key].shape == value.shape, (label, key)
        if key == 'steps':
            # Residual columns of a continued CG agree to the solver precision.
            np.testing.assert_allclose(actual[key][:, 3:], value[:, 3:], rtol=0, atol=1e-9,
                                       err_msg=label)
            np.testing.assert_allclose(actual[key][:, 1], value[:, 1], rtol=1e-6, atol=1e-12,
                                       err_msg=label)
        else:
            np.testing.assert_allclose(actual[key], value, rtol=0,
                                       atol=5e-12 if method == 4 else 1e-10, err_msg=label)


def te_schedule(path, case):
    rows = drive_rows(case.sites, 'onebody')
    text = f'====\nNTimeSteps {len(TE_TIMES)}\n====\n====\n====\n'
    for time, amplitude in zip(TE_TIMES, TE_AMPLITUDES):
        active = rows if amplitude else []
        text += f'{time} {len(active)}\n' + ''.join(
            ' '.join(map(str, [*row[:-2], amplitude*row[-2], amplitude*row[-1]]))+'\n'
            for row in active)
    (path/'TEOneBody.def').write_text(text)
    with (path/'sym.def').open('a') as f:
        f.write('TEOneBody TEOneBody.def\nSpectrumVec seed\n')


def continuation_case(root, name, case, method, hphi, launcher, env):
    """Inputs for one solver; TE starts from a CG state of the writer rank count."""
    mod = {3: {'LanczosEps': 12},
           1: {'NumAve': 2, 'Lanczos_max': 20, 'LargeValue': 4, 'ExpecInterval': 1000},
           5: {'NumAve': 2, 'Lanczos_max': 20, 'LargeValue': 50, 'ExpandCoef': 12,
               'ExpecInterval': 1000},
           4: {'Lanczos_max': len(TE_TIMES), 'ExpandCoef': 8, 'OutputInterval': 2,
               'ExpecInterval': 1}}[method]
    calc = {'InputEigenVec': 1, 'OutputEigenVec': 1} if method == 4 else {}
    path = expert_case(root/name, case, method=method,
                       options={'CalcMod': calc, 'ModPara': mod})
    if method == 4:
        seed = save(root, name+'_seed', case, hphi, launcher, env)
        (path/'output').mkdir()
        for rank, file in enumerate(checkpoints(seed, 'zvo_eigenvec_0')):
            shutil.copy2(file, path/'output'/f'seed_rank_{rank}.dat')
        te_schedule(path, case)
    return path


def check_continuations(root, hphi, probe, launcher):
    """Interrupted v3 runs continue with other rank counts and layouts."""
    plans = [('HubbardNConserved', 3, '0:2', 1), ('tJNConserved', 3, '0:2', 1),
             ('HubbardNConserved', 1, '1:13', 4), ('tJNConserved', 5, '1:13', 4),
             ('HubbardNConserved', 4, '0:6', None), ('tJNConserved', 4, '0:6', None)]
    writer = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='distributed')
    for index, (model, method, abort, interval) in enumerate(plans):
        case = Case(model, 4, 3, 1)
        label = f'{model}_method{method}'
        ref = continuation_case(root, label+'_ref', case, method, hphi, launcher, writer)
        run_case(ref, hphi, 'ref', launcher, writer)
        expected = history(ref, method)
        saved = continuation_case(root, label+'_saved', case, method, hphi, launcher, writer)
        set_key(saved, 'calc.def', 'ReStart', 1)
        if interval:
            set_key(saved, 'mod.def', 'CheckpointInterval', interval)
        text = run_case(saved, probe, 'save', launcher,
                        dict(writer, HPHI_TEST_RESTART_ABORT_AFTER=abort))
        assert 'Test hook: stopping after checkpoint publication' in text, text
        for header in (struct.unpack('<40Q', file.read_bytes()[:320]) for file in
                       (saved/'output').glob('zvo_restart_set*_rank_*.dat')):
            assert tuple(header[4:9]) == (MODEL_ID[model], 4, 0, 0, 3) and header[29] == 0
        for position, (ranks, command) in enumerate(readers(launcher)):
            # Alternate the reader layout; each method sees both across the CI jobs
            # and the writer layout is always distributed.
            for layout in [('replicated', 'distributed')[(index+position) % 2]]:
                path = root/f'{label}_to{ranks}_{layout}'
                shutil.copytree(saved, path)
                set_key(path, 'calc.def', 'ReStart', 3)
                set_key(path, 'mod.def', 'CheckpointInterval')
                if method == 4:
                    set_key(path, 'calc.def', 'InputEigenVec', 0)
                text = run_case(path, hphi, 'resume', command,
                                dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout))
                # TPQ set 0 finished before the interruption; set 1 must continue.
                if method in (1, 5):
                    assert 'no remaining steps (set 1' not in text, text
                else:
                    assert 'no remaining steps' not in text, text
                same_history(history(path, method), expected, method, path.name)
        print(label, 'continued with', [r for r, _ in readers(launcher)], 'ranks', flush=True)
    # A saved checkpoint is never continued in another physical space or sector.
    source = root/'HubbardNConserved_method3_saved'
    for name, target, detail in [('model', Case('tJNConserved', 4, 3, 1), 'model'),
                                 ('momentum', Case('HubbardNConserved', 4, 3, 3),
                                  'group/character digest')]:
        path = expert_case(root/('resume_reject_'+name), target, method=3,
                           options={'CalcMod': {'ReStart': 3}, 'ModPara': {'LanczosEps': 12}})
        shutil.copytree(source/'output', path/'output')
        run_case(path, hphi, 'reject', launcher, writer, f'restart checkpoint {detail} differs')
        print('resume rejected', name, flush=True)


def main():
    hphi, probe = (Path(p).resolve() for p in sys.argv[1:3])
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    assert mpi_size(launcher) is not None
    root = Path('symmetry_nconserved_io_restart')
    if root.exists():
        shutil.rmtree(root)
    root.mkdir()
    check_states(root, hphi, launcher, 'replicated', ('HubbardNConserved',))
    check_states(root, hphi, launcher, 'distributed')
    check_mismatches(root, hphi, launcher)
    check_fulldiag_files(root, hphi)
    check_restart_identity(root, hphi, launcher)
    check_continuations(root, hphi, probe, launcher)
    print('NConserved state/matrix/restart identity and continuation passed')


if __name__ == '__main__':
    main()
