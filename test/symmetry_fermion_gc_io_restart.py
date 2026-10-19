#!/usr/bin/env python3
"""GC/parity physical identity of state, matrix and restart files, and their reuse."""
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
from symmetry_spingc_common import checkpoint_digest, mpi_size, run_case
from symmetry_spingc_te import checkpoints
from symmetry_fermion_gc_common import (INPUT_MODEL, PARITY_FLAGS, expert_case,
                                        read_gc_checkpoint)
from symmetry_fermion_gc_evolution import drive_rows
from symmetry_fermion_gc_reference import (Case, make_reference, pair_rows, raw_words,
                                           sector_matrix, width)
from symmetry_nconserved_io_restart import history, read_matrix, readers, same_history, set_key

IDENTITY_FAILED = 'checkpoint header / sector / layout validation failed'
TE_TIMES = (0., .03, .09, .09, .13, .2, .25, .3)
TE_AMPLITUDES = (0., .5, -.3, -.3, .2, .4, -.1, .3)
CASES = {'HubbardGC': Case('HubbardGC', 4, 1, 1), 'SpinlessFermionGC': Case('SpinlessFermionGC', 6, -1, 1)}


def gc_case(root, name, case, method, options=None, pairs=True):
    return expert_case(root/name, case, method=method, options=options,
                       pairs=pair_rows(case.model, case.sites) if pairs else None)


def v2_parts(path, case, ranks, layout):
    """Version 2, phase 1, zero mask, no fixed N, and the parity in FIXED_FLAGS."""
    parts = [read_gc_checkpoint(file, case) for file in checkpoint_files(path)]
    assert len(parts) == ranks
    offset = 0
    for rank, (header, vector) in enumerate(parts):
        assert header[9] == len(raw_words(case))
        assert tuple(header[11:15]) == (ranks, rank, int(layout == 'distributed'), offset)
        offset += len(vector)
    assert offset == parts[0][0][10]
    assert all(tuple(h[26:28]) == checkpoint_digest(parts) for h, _ in parts)
    return np.concatenate([v for _, v in parts])


def save(root, name, case, hphi, launcher, env, options=None, pairs=True):
    calc = dict({'OutputEigenVec': 1}, **(options or {}))
    path = gc_case(root, name, case, 3, options={'CalcMod': calc, 'ModPara': {'LanczosEps': 18}},
                   pairs=pairs)
    run_case(path, hphi, 'save', launcher, env)
    return path


def import_into(root, name, source, target_case, mod=None, calc=None, pairs=True):
    """Another physical space with the source state files as its input vector."""
    path = gc_case(root, name, target_case, 3, pairs=pairs,
                   options={'CalcMod': dict({'InputEigenVec': 1}, **(calc or {})),
                            'ModPara': mod or {}})
    (path/'output').mkdir()
    for file in checkpoint_files(source):
        shutil.copy2(file, path/'output'/file.name)
    return path


def scale_pairs(path, factor):
    """Change only the pair coefficients; each Hermitian partner changes with its row."""
    lines = (path/'AnomalousTerm.def').read_text().splitlines()
    for i in range(5, len(lines)):
        row = lines[i].split()
        row[-2:] = [repr(float(value)*factor) for value in row[-2:]]
        lines[i] = ' '.join(row)
    (path/'AnomalousTerm.def').write_text('\n'.join(lines)+'\n')


def check_states(root, hphi, launcher, layout, models=tuple(CASES)):
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout)
    ranks = mpi_size(launcher)
    for model in models:
        case = CASES[model]
        h = sector_matrix(make_reference(case, pairs=pair_rows(model, case.sites)))
        source = save(root, f'{model}_{layout}', case, hphi, launcher, env)
        vector = v2_parts(source, case, ranks, layout)
        check_eigenvector(vector, h, np.linalg.eigvalsh(h)[0])
        assert read_manifest(source)['model'] == model
        same = import_copy(root, f'{model}_{layout}_cg', source)
        with (same/'calc.def').open('a') as f:
            f.write('OutputEigenVec 1\n')
        run_case(same, hphi, 'cg_to_cg', launcher, env)
        assert all(a.read_bytes() == b.read_bytes()
                   for a, b in zip(checkpoint_files(source), checkpoint_files(same)))
        # Changing only the pairs keeps the space: an accepted quench of the same vector.
        quench = import_copy(root, f'{model}_{layout}_quench', source)
        scale_pairs(quench, 1.1)
        with (quench/'calc.def').open('a') as f:
            f.write('OutputEigenVec 1\n')
        assert 'hamiltonian_changed=yes' in run_case(quench, hphi, 'quench', launcher, env)
        np.testing.assert_array_equal(v2_parts(quench, case, ranks, layout), vector)
        print(model, layout, 'v2 identity, CG reuse and pair quench passed', flush=True)


def check_mismatches(root, hphi, launcher):
    """Equal dimensions in another space, parity or character, and corrupted flags."""
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated')
    even = Case('HubbardGC', 4, 1, 1)
    source = save(root, 'source_even', even, hphi, launcher, env)
    for name, target, detail in [
            ('parity', Case('HubbardGC', 4, -1, 1), None),
            ('full', Case('HubbardGC', 3, 0, 1), None),
            ('momentum', Case('HubbardGC', 4, 1, 3), 'group/character digest'),
            ('model', Case('SpinlessFermionGC', 8, 1, 1), None)]:
        if name in ('parity', 'momentum'):
            assert make_reference(target).basis.shape[1] == make_reference(even).basis.shape[1]
        path = import_into(root, 'reject_'+name, source, target)
        text = run_case(path, hphi, 'reject', launcher, env, IDENTITY_FAILED)
        assert (detail or '') + ' differs' in text, text
        print('rejected', name, flush=True)
    # The flags word alone separates the parities and rejects invalid values.
    for name, value in (('odd', PARITY_FLAGS[-1]), ('full', PARITY_FLAGS[0]), ('bit3', 8)):
        path = import_into(root, 'reject_flags_'+name, source, even)
        for file in checkpoint_files(path):
            data = bytearray(file.read_bytes())
            data[29*8:30*8] = struct.pack('<Q', value)
            file.write_bytes(bytes(data))
        text = run_case(path, hphi, 'reject', launcher, env, IDENTITY_FAILED)
        assert 'fixed quantity flags differs' in text, text
    # A fixed-N sector of the same sites is never an implicit import.
    path = gc_case(root, 'reject_fixed_n', even, 3, pairs=False,
                   options={'CalcMod': {'CalcModel': 0, 'InputEigenVec': 1},
                            'ModPara': {'Ncond': 4, 'ParticleParity': None}})
    (path/'output').mkdir()
    for file in checkpoint_files(source):
        shutil.copy2(file, path/'output'/file.name)
    assert 'differs' in run_case(path, hphi, 'reject', launcher, env, IDENTITY_FAILED)
    print('rejected flags and fixed-N imports', flush=True)


def check_fulldiag_files(root, hphi):
    """LAPACK FullDiag states, matrix export and import on one rank."""
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated')
    for model, case in CASES.items():
        h = sector_matrix(make_reference(case, pairs=pair_rows(model, case.sites)))
        exporter = gc_case(root, f'{model}_ham', case, 2, options={'CalcMod': {'OutputHam': 1}})
        run_case(exporter, hphi, 'export', [], env)
        meta, matrix = read_matrix(exporter/'output/zvo_ham_sector.mtx')
        assert meta['model'] == model and meta['nup'] == meta['ndown'] == meta['ne'] == '0'
        assert meta['fixed_flags'] == str(PARITY_FLAGS[case.parity]) and meta['phase'] == '1'
        assert int(meta['full_dim']) == len(raw_words(case))
        np.testing.assert_allclose(matrix, h, rtol=0, atol=1e-13)
        importer = gc_case(root, f'{model}_ham_in', case, 2,
                           options={'CalcMod': {'InputHam': 1, 'OutputEigenVec': 1}})
        (importer/'output').mkdir()
        shutil.copy2(exporter/'output/zvo_ham_sector.mtx', importer/'output')
        run_case(importer, hphi, 'import', [], env)
        energies = np.loadtxt(importer/'output/zvo_energy_sector.dat', ndmin=2)[:, 1]
        np.testing.assert_allclose(energies, np.linalg.eigvalsh(h), rtol=0, atol=1e-10)
        # The exported ground state re-enters CG through the v2 reader.
        consumer = gc_case(root, f'{model}_fd_to_cg', case, 3,
                           options={'CalcMod': {'InputEigenVec': 1}})
        (consumer/'output').mkdir()
        shutil.copy2(importer/'output/zvo_eigenvec_0_rank_0.dat', consumer/'output')
        run_case(consumer, hphi, 'fd_to_cg', [], env)
        assert abs(parse_energy(consumer)[0]-energies[0]) <= 1e-10
        print(model, 'OutputHam/InputHam and FullDiag to CG passed', flush=True)
    # The matrix identity rejects the other parity of equal dimension.
    source = gc_case(root, 'ham_even', Case('HubbardGC', 4, 1, 1), 2,
                     options={'CalcMod': {'OutputHam': 1}})
    run_case(source, hphi, 'export', [], env)
    target = gc_case(root, 'ham_into_odd', Case('HubbardGC', 4, -1, 1), 2,
                     options={'CalcMod': {'InputHam': 1}})
    (target/'output').mkdir()
    shutil.copy2(source/'output/zvo_ham_sector.mtx', target/'output')
    run_case(target, hphi, 'reject', [], env,
             'metadata fixed_flags does not match this sector (file=4, expected=12)')


def check_restart_identity(root, hphi, launcher):
    """v3 solver files reuse the 30 identity words with the parity flags."""
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='distributed')
    ranks = mpi_size(launcher)
    # The even/odd fixtures and the full GC space (flags 0) of each model.
    cases = [*CASES.values(), Case('HubbardGC', 3, 0, 1), Case('SpinlessFermionGC', 6, 0, 1)]
    for case in cases:
        model = case.model
        path = gc_case(root, f'{model}_parity{case.parity}_restart', case, 3,
                       options={'CalcMod': {'ReStart': 1}})
        run_case(path, hphi, 'restart', launcher, env)
        files = sorted((path/'output').glob('zvo_restart_set0_comp*_gen0_rank_*.dat'))
        assert len(files) == 2*ranks
        for file in files:
            header = struct.unpack('<40Q', file.read_bytes()[:320])
            assert header[1] == 3 and header[2] == 1
            assert tuple(header[4:9]) == (INPUT_MODEL[model], case.sites, 0, 0, 0)
            assert header[28] == 0 and header[29] == PARITY_FLAGS[case.parity] and header[30] == 1
        assert 'basis sector' in (path/'output/zvo_restart_manifest.dat').read_text()


def te_schedule(path, case):
    rows = drive_rows(case, 'onebody')
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
    path = gc_case(root, name, case, method, options={'CalcMod': calc, 'ModPara': mod})
    if method == 4:
        seed = save(root, name+'_seed', case, hphi, launcher, env)
        (path/'output').mkdir()
        for rank, file in enumerate(checkpoints(seed, 'zvo_eigenvec_0')):
            shutil.copy2(file, path/'output'/f'seed_rank_{rank}.dat')
        te_schedule(path, case)
    return path


def check_continuations(root, hphi, probe, launcher):
    """Interrupted v3 runs continue with other rank counts and layouts."""
    plans = [('HubbardGC', 3, '0:2', 1), ('SpinlessFermionGC', 3, '0:2', 1),
             ('HubbardGC', 1, '1:13', 4), ('SpinlessFermionGC', 5, '1:13', 4),
             ('HubbardGC', 4, '0:6', None), ('SpinlessFermionGC', 4, '0:6', None)]
    writer = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='distributed')
    for index, (model, method, abort, interval) in enumerate(plans):
        case = CASES[model]
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
            assert tuple(header[4:9]) == (INPUT_MODEL[model], case.sites, 0, 0, 0)
            assert header[29] == PARITY_FLAGS[case.parity]
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
                # A TE checkpoint records the driven H of its step, not the base H.
                if method != 4:
                    assert 'hamiltonian_changed=no' in text, text
                # TPQ set 0 finished before the interruption; set 1 must continue.
                if method in (1, 5):
                    assert 'no remaining steps (set 1' not in text, text
                else:
                    assert 'no remaining steps' not in text, text
                same_history(history(path, method), expected, method, path.name)
        print(label, 'continued with', [r for r, _ in readers(launcher)], 'ranks', flush=True)
    source = root/'HubbardGC_method3_saved'
    # The restart design reports, not rejects, a changed H; a pair-only change is detected.
    path = root/'resume_changed_pairs'
    shutil.copytree(source, path)
    set_key(path, 'calc.def', 'ReStart', 3)
    set_key(path, 'mod.def', 'CheckpointInterval')
    scale_pairs(path, 1.1)
    assert 'hamiltonian_changed=yes' in run_case(path, hphi, 'resume', launcher, writer)
    # A saved checkpoint is never continued in another parity or sector.
    for name, target, detail in [('parity', Case('HubbardGC', 4, -1, 1), 'fixed quantity flags'),
                                 ('momentum', Case('HubbardGC', 4, 1, 3), 'group/character digest'),
                                 ('model', Case('SpinlessFermionGC', 8, 1, 1), 'model')]:
        path = gc_case(root, 'resume_reject_'+name, target, 3,
                       options={'CalcMod': {'ReStart': 3}, 'ModPara': {'LanczosEps': 12}})
        shutil.copytree(source/'output', path/'output')
        text = run_case(path, hphi, 'reject', launcher, writer, 'restart checkpoint')
        assert 'differs' in text, text
        print('resume rejected', name, flush=True)


def main():
    hphi, probe = (Path(p).resolve() for p in sys.argv[1:3])
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    assert mpi_size(launcher) is not None
    root = Path('symmetry_fermion_gc_io_restart')
    if root.exists():
        shutil.rmtree(root)
    root.mkdir()
    check_states(root, hphi, launcher, 'replicated', ('HubbardGC',))
    check_states(root, hphi, launcher, 'distributed')
    check_mismatches(root, hphi, launcher)
    check_fulldiag_files(root, hphi)
    check_restart_identity(root, hphi, launcher)
    check_continuations(root, hphi, probe, launcher)
    print('Fermion GC state/matrix/restart identity and continuation passed')


if __name__ == '__main__':
    main()
