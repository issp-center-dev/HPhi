#!/usr/bin/env python3
"""Px identity of state, matrix and restart files, their reuse, and equal-dimension rejections."""
from pathlib import Path
import os
import shlex
import shutil
import struct
import sys
import numpy as np
from symmetry_kondo_common import read_manifest
from symmetry_kondo_observables import parse_energy
from symmetry_fermion_gc_core import check_ground_state
from symmetry_spingc_checkpoint import checkpoint_files, import_copy
from symmetry_spingc_common import checkpoint_digest, mpi_size, read_checkpoint, run_case
from symmetry_spingc_te import checkpoints
from symmetry_nconserved_io_restart import history, read_matrix, readers, same_history, set_key
from symmetry_spinflip_common import INPUT_MODEL, check_manifest, expert_case
from symmetry_spinflip_evolution import drive_rows
from symmetry_spinflip_reference import Case, make_reference, raw_words, sector_matrix

IDENTITY_FAILED = 'checkpoint header / sector / layout validation failed'
TE_TIMES = (0., .03, .09, .09, .13, .2, .25, .3)
TE_AMPLITUDES = (0., .5, -.3, -.3, .2, .4, -.1, .3)
CASES = {'SpinGC': Case('SpinGC', 8, -1, 1), 'Spin': Case('Spin', 8, -1, 1)}
# (name, source, target) of equal sector dimension: another p at the same k,
# another k at the same p, and Px off/on for the same group and k (Spin L=4
# k=1 has D(k,+)=0, so D(k,-)=D(k)).
EQUAL_DIMENSION = (('parity', Case('SpinGC', 5, 1, 1), Case('SpinGC', 5, -1, 1)),
                   ('momentum', Case('SpinGC', 5, 1, 1), Case('SpinGC', 5, 1, 2)),
                   ('no_px', Case('Spin', 4, -1, 1), Case('Spin', 4, 0, 1)),
                   ('px_from_off', Case('Spin', 4, 0, 1), Case('Spin', 4, -1, 1)))


def equal_dimension_pairs():
    for name, source, target in EQUAL_DIMENSION:
        assert make_reference(source).basis.shape[1] == make_reference(target).basis.shape[1], name
        yield name, source, target


def manifest_group(path):
    return int(read_manifest(path)['group_digest'].split(':')[1], 16)


def fixed_numbers(case):
    return (case.sites//2,)*3 if case.model == 'Spin' else (0, 0, 0)


def v1_parts(path, case, ranks, layout):
    """Version 1, phase 1: the Px sector keeps the v1 words; GROUP carries Px and p."""
    parts = [read_checkpoint(file) for file in checkpoint_files(path)]
    assert len(parts) == ranks
    group = manifest_group(path)
    offset = 0
    for rank, (header, vector) in enumerate(parts):
        assert tuple(header[4:9]) == (INPUT_MODEL[case.model], case.sites, *fixed_numbers(case)), header
        assert header[9] == len(raw_words(case)) and header[16] == group
        assert tuple(header[11:15]) == (ranks, rank, int(layout == 'distributed'), offset)
        offset += len(vector)
    assert offset == parts[0][0][10]
    assert all(tuple(h[26:28]) == checkpoint_digest(parts) for h, _ in parts)
    return np.concatenate([v for _, v in parts])


def save(root, name, case, hphi, launcher, env, options=None, parameters=None):
    calc = dict({'OutputEigenVec': 1}, **(options or {}))
    path = expert_case(root/name, case, method=3, parameters=parameters,
                       options={'CalcMod': calc, 'ModPara': {'LanczosEps': 22}})
    run_case(path, hphi, 'save', launcher, env)
    return path


def import_into(root, name, source, target_case, mod=None):
    """Another sector with the source state files as its input vector."""
    path = expert_case(root/name, target_case, method=3,
                       options={'CalcMod': {'InputEigenVec': 1}, 'ModPara': mod or {}})
    (path/'output').mkdir()
    for file in checkpoint_files(source):
        shutil.copy2(file, path/'output'/file.name)
    return path


def check_states(root, hphi, launcher, layout, models=tuple(CASES)):
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout)
    ranks = mpi_size(launcher)
    for model in models:
        case = CASES[model]
        h = sector_matrix(make_reference(case))
        source = save(root, f'{model}_{layout}', case, hphi, launcher, env)
        vector = v1_parts(source, case, ranks, layout)
        check_ground_state(vector, h, np.linalg.eigvalsh(h)[0], source.name)
        check_manifest(read_manifest(source), case, len(h))
        same = import_copy(root, f'{model}_{layout}_cg', source)
        with (same/'calc.def').open('a') as f:
            f.write('OutputEigenVec 1\n')
        run_case(same, hphi, 'cg_to_cg', launcher, env)
        assert all(a.read_bytes() == b.read_bytes()
                   for a, b in zip(checkpoint_files(source), checkpoint_files(same)))
        # Another Px-even H of the same sector: an accepted quench of the same vector.
        quench = import_copy(root, f'{model}_{layout}_quench', source)
        families = expert_case(root/f'{model}_{layout}_quench_h', case, method=3,
                               parameters=dict(jz=.5) if model == 'SpinGC' else dict(j2=.4))
        for file in families.glob('*.def'):
            if file.name not in ('calc.def', 'mod.def'):
                shutil.copy2(file, quench/file.name)
        with (quench/'calc.def').open('a') as f:
            f.write('OutputEigenVec 1\n')
        assert 'hamiltonian_changed=yes' in run_case(quench, hphi, 'quench', launcher, env)
        np.testing.assert_array_equal(v1_parts(quench, case, ranks, layout), vector)
        print(model, layout, 'v1 identity, CG reuse and Px-even quench passed', flush=True)


def check_mismatches(root, hphi, launcher):
    """Equal dimensions in another p, k or without Px are never imported."""
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated')
    for name, source_case, target in equal_dimension_pairs():
        source = save(root, 'source_'+name, source_case, hphi, launcher, env)
        path = import_into(root, 'reject_'+name, source, target)
        text = run_case(path, hphi, 'reject', launcher, env, IDENTITY_FAILED)
        assert 'group/character digest differs' in text, text
        print('rejected an equal-dimension state of', name, flush=True)
    # Spin and SpinGC sectors never share files.
    path = import_into(root, 'reject_model', root/'source_parity', Case('Spin', 6, 1, 0))
    assert 'differs' in run_case(path, hphi, 'reject', launcher, env, IDENTITY_FAILED)


def check_fulldiag_files(root, hphi):
    """LAPACK FullDiag states, matrix export and import on one rank."""
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated')
    for case in (Case('SpinGC', 6, -1, 1), Case('Spin', 8, 1, 0)):
        h = sector_matrix(make_reference(case))
        exporter = expert_case(root/f'{case.model}_ham', case, method=2,
                               options={'CalcMod': {'OutputHam': 1}})
        run_case(exporter, hphi, 'export', [], env)
        text = (exporter/'output/zvo_ham_sector.mtx').read_text()
        assert f'% spin_flip_operator=product_sigma_x spin_flip_parity={case.parity:+d}\n' in text
        meta, matrix = read_matrix(exporter/'output/zvo_ham_sector.mtx')
        assert meta['model'] == case.model and meta['phase'] == '1' and meta['fixed_flags'] == '0'
        assert (int(meta['nup']), int(meta['ndown']), int(meta['ne'])) == fixed_numbers(case), meta
        assert int(meta['full_dim']) == len(raw_words(case))
        assert int(meta['group_digest'], 16) == manifest_group(exporter), meta
        np.testing.assert_allclose(matrix, h, rtol=0, atol=1e-13)
        importer = expert_case(root/f'{case.model}_ham_in', case, method=2,
                               options={'CalcMod': {'InputHam': 1, 'OutputEigenVec': 1}})
        (importer/'output').mkdir()
        shutil.copy2(exporter/'output/zvo_ham_sector.mtx', importer/'output')
        run_case(importer, hphi, 'import', [], env)
        energies = np.loadtxt(importer/'output/zvo_energy_sector.dat', ndmin=2)[:, 1]
        np.testing.assert_allclose(energies, np.linalg.eigvalsh(h), rtol=0, atol=1e-10)
        # The FullDiag ground state re-enters CG through the v1 reader.
        consumer = expert_case(root/f'{case.model}_fd_to_cg', case, method=3,
                               options={'CalcMod': {'InputEigenVec': 1}})
        (consumer/'output').mkdir()
        shutil.copy2(importer/'output/zvo_eigenvec_0_rank_0.dat', consumer/'output')
        run_case(consumer, hphi, 'fd_to_cg', [], env)
        assert abs(parse_energy(consumer)[0]-energies[0]) <= 1e-10
        print(case.model, 'OutputHam/InputHam and FullDiag to CG passed', flush=True)
    # The matrix identity rejects another p, k or Px setting of equal dimension.
    for name, source_case, target in equal_dimension_pairs():
        source = expert_case(root/f'ham_{name}', source_case, method=2,
                             options={'CalcMod': {'OutputHam': 1}})
        run_case(source, hphi, 'export', [], env)
        path = expert_case(root/f'ham_into_{name}', target, method=2,
                           options={'CalcMod': {'InputHam': 1}})
        (path/'output').mkdir()
        shutil.copy2(source/'output/zvo_ham_sector.mtx', path/'output')
        run_case(path, hphi, 'reject', [], env, 'metadata group_digest does not match this sector')
    print('matrix files of another p, k or without Px rejected', flush=True)


def check_restart_identity(root, hphi, launcher):
    """v3 solver files reuse the v1 identity words, the group digest included."""
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='distributed')
    ranks = mpi_size(launcher)
    for case in (*CASES.values(), Case('Spin', 4, -1, 1)):
        path = expert_case(root/f'{case.model}_L{case.sites}_restart', case, method=3,
                           options={'CalcMod': {'ReStart': 1}})
        run_case(path, hphi, 'restart', launcher, env)
        files = sorted((path/'output').glob('zvo_restart_set0_comp*_gen0_rank_*.dat'))
        assert len(files) == 2*ranks
        for file in files:
            header = struct.unpack('<40Q', file.read_bytes()[:320])
            assert header[1] == 3 and header[2] == 1
            assert tuple(header[4:9]) == (INPUT_MODEL[case.model], case.sites, *fixed_numbers(case))
            assert header[16] == manifest_group(path)
            assert header[28] == header[29] == 0 and header[30] == 1
        assert 'basis sector' in (path/'output/zvo_restart_manifest.dat').read_text()


def te_schedule(path, case):
    rows = drive_rows(case, 'onebody' if case.model == 'SpinGC' else 'twobody')
    family = 'TEOneBody' if case.model == 'SpinGC' else 'TETwoBody'
    text = f'====\nNTimeSteps {len(TE_TIMES)}\n====\n====\n====\n'
    for time, amplitude in zip(TE_TIMES, TE_AMPLITUDES):
        active = rows if amplitude else []
        text += f'{time} {len(active)}\n' + ''.join(
            ' '.join(map(str, [*row[:-2], amplitude*row[-2], amplitude*row[-1]]))+'\n'
            for row in active)
    (path/f'{family}.def').write_text(text)
    with (path/'sym.def').open('a') as f:
        f.write(f'{family} {family}.def\nSpectrumVec seed\n')


def continuation_case(root, name, case, method, hphi, launcher, env):
    """Inputs for one solver; TE starts from a CG state of the writer rank count."""
    mod = {3: {'LanczosEps': 12},
           1: {'NumAve': 2, 'Lanczos_max': 20, 'LargeValue': 4, 'ExpecInterval': 1000},
           5: {'NumAve': 2, 'Lanczos_max': 20, 'LargeValue': 50, 'ExpandCoef': 12,
               'ExpecInterval': 1000},
           4: {'Lanczos_max': len(TE_TIMES), 'ExpandCoef': 8, 'OutputInterval': 2,
               'ExpecInterval': 1}}[method]
    calc = {'InputEigenVec': 1, 'OutputEigenVec': 1} if method == 4 else {}
    path = expert_case(root/name, case, method=method, options={'CalcMod': calc, 'ModPara': mod})
    if method == 4:
        seed = save(root, name+'_seed', case, hphi, launcher, env)
        (path/'output').mkdir()
        for rank, file in enumerate(checkpoints(seed, 'zvo_eigenvec_0')):
            shutil.copy2(file, path/'output'/f'seed_rank_{rank}.dat')
        te_schedule(path, case)
    return path


def check_continuations(root, hphi, probe, launcher):
    """Interrupted v3 runs continue with other rank counts and layouts."""
    plans = [('SpinGC', 3, '0:2', 1), ('Spin', 3, '0:2', 1),
             ('SpinGC', 1, '1:13', 4), ('Spin', 5, '1:13', 4),
             ('SpinGC', 4, '0:6', None), ('Spin', 4, '0:6', None)]
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
        group = manifest_group(saved)
        for header in (struct.unpack('<40Q', file.read_bytes()[:320]) for file in
                       (saved/'output').glob('zvo_restart_set*_rank_*.dat')):
            assert tuple(header[4:9]) == (INPUT_MODEL[model], case.sites, *fixed_numbers(case))
            assert header[16] == group
        for position, (ranks, command) in enumerate(readers(launcher)):
            # Alternate the reader layout; each method sees both across the CI jobs
            # and the writer layout is always distributed.
            layout = ('replicated', 'distributed')[(index+position) % 2]
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
    # A saved checkpoint is never continued in another p, k or Px setting,
    # even with an equal sector dimension.
    for name, source_case, target in equal_dimension_pairs():
        source = continuation_case(root, f'resume_source_{name}', source_case, 3, hphi,
                                   launcher, writer)
        set_key(source, 'calc.def', 'ReStart', 1)
        run_case(source, hphi, 'source', launcher, writer)
        path = expert_case(root/f'resume_reject_{name}', target, method=3,
                           options={'CalcMod': {'ReStart': 3}, 'ModPara': {'LanczosEps': 12}})
        shutil.copytree(source/'output', path/'output')
        text = run_case(path, hphi, 'reject', launcher, writer, 'restart checkpoint')
        assert 'group/character digest differs' in text, text
        print('resume rejected', name, flush=True)


def main():
    hphi, probe = (Path(p).resolve() for p in sys.argv[1:3])
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    assert mpi_size(launcher) is not None
    root = Path('symmetry_spinflip_io_restart')
    if root.exists():
        shutil.rmtree(root)
    root.mkdir()
    check_states(root, hphi, launcher, 'replicated', ('SpinGC',))
    check_states(root, hphi, launcher, 'distributed')
    check_mismatches(root, hphi, launcher)
    check_fulldiag_files(root, hphi)
    check_restart_identity(root, hphi, launcher)
    check_continuations(root, hphi, probe, launcher)
    print('Spin-flip state/matrix/restart identity and continuation passed')


if __name__ == '__main__':
    main()
