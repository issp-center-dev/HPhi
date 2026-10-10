#!/usr/bin/env python3
"""Kondo checkpoint physical identity, collective rejection and state imports."""
import os
from pathlib import Path
import shlex
import shutil
import struct
import sys
import tempfile
import numpy as np
from symmetry_kondo_common import expert_case, read_kondo_checkpoint, read_manifest
from symmetry_kondo_reference import Case
from symmetry_spingc_common import run_case, mpi_size, checkpoint_digest
from symmetry_spingc_checkpoint import checkpoint_files, import_copy


def header_check(path, case, ranks, layout):
    parts = [read_kondo_checkpoint(file) for file in checkpoint_files(path)]
    assert len(parts) == ranks
    local_sites = range(case.cells) if case.layout == 'block' else range(0, 2*case.cells, 2)
    model, fixed, quantities = {
        'Kondo': (2, 3, (3, 2, 5)),
        'KondoNConserved': (12, 1, (0, 0, 5)),
        'KondoGC': (5, 0, (0, 0, 0)),
    }[case.model]
    offset = 0
    for rank, (header, payload) in enumerate(parts):
        assert len(header) == 30 and tuple(header[1:4]) == (2, 2, 128)
        assert tuple(header[4:9]) == (model, 6, *quantities)
        assert tuple(header[11:15]) == (ranks, rank, int(layout == 'distributed'), offset)
        assert header[28] == sum(1 << i for i in local_sites)
        assert header[29] == fixed
        assert len(payload) == header[15] and np.isfinite(payload).all()
        offset += len(payload)
    assert offset == parts[0][0][10]
    assert all(tuple(h[26:28]) == checkpoint_digest(parts) for h, _ in parts)
    assert abs(sum(np.vdot(v, v).real for _, v in parts)-1) <= 1e-8
    manifest = read_manifest(path)
    keys = list(manifest)
    extra = ['n_local_spin', 'n_conduction_sites', 'local_site_mask',
             'translation_convention', 'particle_number_convention', 'basis_space_digest',
             'fixed_quantities']
    fixed_keys = {'Kondo': ['fixed_ncond', 'fixed_2sz', 'fixed_nup', 'fixed_ndown', 'fixed_ne'],
                  'KondoNConserved': ['fixed_ncond', 'fixed_ne'], 'KondoGC': []}[case.model]
    start = keys.index('n_local_spin')
    assert keys[start:start+len(extra+fixed_keys)+1] == extra+fixed_keys+['momentum_index']
    assert [key for key in keys if key.startswith('fixed_')] == ['fixed_quantities']+fixed_keys
    return parts



def store_parts(files, parts, refresh=False):
    digest = checkpoint_digest(parts) if refresh else None
    for file, (header, vector) in zip(files, parts):
        header = list(header)
        if digest is not None:
            header[26:28] = digest
        file.write_bytes(struct.pack('<30Q', *header)+np.asarray(vector, dtype='<c16').tobytes())


def joined(path, label='zvo_eigenvec_0'):
    from symmetry_spingc_te import checkpoints
    parts = [read_kondo_checkpoint(file) for file in checkpoints(path, label)]
    return np.concatenate([v for _, v in parts])


def imports(root, source, case, hphi, launcher, env):
    cg = import_copy(root, source.name+'_cg', source)
    with (cg/'calc.def').open('a') as f:
        f.write('OutputEigenVec 1\n')
    run_case(cg, hphi, 'cg_to_cg', launcher, env)
    assert all(a.read_bytes() == b.read_bytes() for a, b in
               zip(checkpoint_files(source), checkpoint_files(cg)))
    quench = import_copy(root, source.name+'_quench', source)
    lines = (quench/'Trans.def').read_text().splitlines()
    for i in range(5, len(lines)):
        row = lines[i].split()
        row[-2:] = [str(float(value)*1.1) for value in row[-2:]]
        lines[i] = ' '.join(row)
    (quench/'Trans.def').write_text('\n'.join(lines)+'\n')
    with (quench/'calc.def').open('a') as f:
        f.write('OutputEigenVec 1\n')
    text = run_case(quench, hphi, 'quench', launcher, env)
    assert 'hamiltonian_changed=yes' in text
    np.testing.assert_array_equal(joined(quench), joined(source))
    for a, b in zip(checkpoint_files(source), checkpoint_files(quench)):
        old, _ = read_kondo_checkpoint(a)
        new, _ = read_kondo_checkpoint(b)
        assert old[21] != new[21]

    from symmetry_spingc_te import checkpoints, schedule
    def evolve(name, seed, label, times):
        path = expert_case(root/name, case, method=4,
                          options={'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1},
                                   'ModPara': {'Lanczos_max': 5, 'ExpandCoef': 8,
                                               'ExpecInterval': 1, 'OutputInterval': 1}})
        (path/'output').mkdir()
        for rank, file in enumerate(checkpoints(seed, label)):
            shutil.copy2(file, path/'output'/f'seed_rank_{rank}.dat')
        with (path/'sym.def').open('a') as f:
            f.write('SpectrumVec seed\n')
        schedule(path, 'TEOneBody', [], times=times)
        run_case(path, hphi, name, launcher, env)
        np.testing.assert_allclose(joined(path), joined(seed, label), atol=3e-11, rtol=0)
        header, _ = read_kondo_checkpoint(checkpoints(path, 'zvo_eigenvec_0')[0])
        assert header[22] == 4 and header[24] == 0
        assert struct.unpack('<d', struct.pack('<Q', header[25]))[0] == times[0]
        return path
    te = evolve(source.name+'_cg_te', source, 'zvo_eigenvec_0', [0., .01, .02, .03, .04])
    te2 = evolve(source.name+'_te_te', te, 'zvo_eigenvec_final', [0., .002, .004, .006, .008])
    meta = read_manifest(te2)
    assert meta['source_method'] == '4' and meta['source_step'] == '4'
    assert float(meta['source_time']) == .04
    te_cg = import_copy(root, source.name+'_te_cg', source)
    for rank, file in enumerate(checkpoints(te, 'zvo_eigenvec_final')):
        shutil.copy2(file, te_cg/'output'/f'zvo_eigenvec_0_rank_{rank}.dat')
    with (te_cg/'calc.def').open('a') as f:
        f.write('OutputEigenVec 1\n')
    text = run_case(te_cg, hphi, 'te_to_cg', launcher, env)
    assert 'step=4 time=0.040000000000000001' in text
    np.testing.assert_allclose(joined(te_cg), joined(te, 'zvo_eigenvec_final'), atol=3e-11, rtol=0)
    if mpi_size(launcher) > 1:
        mixed = import_copy(root, source.name+'_mixed_actual', source)
        shutil.copy2(checkpoints(te, 'zvo_eigenvec_final')[-1], checkpoint_files(mixed)[-1])
        run_case(mixed, hphi, 'mixed_actual', launcher, env,
                 'checkpoint cross-rank metadata consistency failed')
    print(source.name, 'CG→CG, CG→TE, TE→CG, TE→TE and quench passed', flush=True)


def corruption_checks(root, source, hphi, launcher, env):
    def reject(name, mutate, diagnostic, detail=None):
        path = import_copy(root, name, source)
        files = checkpoint_files(path)
        parts = [(list(h), v.copy()) for h, v in map(read_kondo_checkpoint, files)]
        mutate(files, parts)
        text = run_case(path, hphi, name, launcher, env, 'checkpoint '+diagnostic+' failed')
        if detail:
            assert detail+' differs' in text, text
        assert not (path/'output/zvo_energy.dat').exists()
        print('rejected', name, diagnostic, flush=True)
        return text
    identity = 'header / sector / layout validation'
    for word, name, detail in [(0, 'magic', None), (1, 'unknown_version', None),
            (2, 'phase', 'phase convention'), (3, 'scalar', 'scalar encoding'),
            (4, 'unknown_model', 'model'), (5, 'nsite', 'nsite'),
            (6, 'nup', 'nup'), (7, 'ndown', 'ndown'), (8, 'ne', 'particle number'),
            (9, 'full_dim', 'raw dimension'), (10, 'dimension', 'sector dimension'),
            (11, 'rank_count', 'MPI rank count'), (12, 'rank', 'MPI rank'),
            (13, 'layout', 'basis layout'), (14, 'offset', 'local offset'),
            (15, 'local_dim', 'local dimension'), (16, 'group', 'group/character digest'),
            (17, 'sector_count', 'sector count'), (18, 'sector_xor', 'sector XOR'),
            (19, 'sector_sum', 'sector sum'), (20, 'order', 'basis order'),
            (28, 'mask', 'local site mask'), (29, 'flags', 'fixed quantity flags')]:
        def mutate(files, parts, word=word):
            parts[-1][0][word] = 99 if word in (1, 4) else int(parts[-1][0][word]) ^ 1
            store_parts(files, parts)
        reject(name, mutate, identity, detail)
    for name, change, stage in [
            ('truncated_prefix', lambda b: b[:8], identity),
            ('truncated_header', lambda b: b[:29*8], identity),
            ('truncated_payload', lambda b: b[:-1], 'payload length'),
            ('trailing_payload', lambda b: b+b'x', 'payload length')]:
        reject(name, lambda files, parts, change=change: files[0].write_bytes(change(files[0].read_bytes())), stage)
    reject('missing_rank', lambda files, parts: files[-1].unlink(), identity)
    def version_one(files, parts, mixed=False):
        for file, (header, vector) in list(zip(files, parts))[-1:] if mixed else zip(files, parts):
            header[1] = 1
            file.write_bytes(struct.pack('<28Q', *header[:28])+vector.astype('<c16').tobytes())
    reject('v1_kondo', version_one, identity, 'version')
    if mpi_size(launcher) > 1:
        reject('mixed_versions', lambda f, p: version_one(f, p, True), 'version agreement')
        for word, name in [(21, 'H'), (22, 'method'), (23, 'state'), (24, 'step'), (25, 'time')]:
            def mixed(files, parts, word=word):
                parts[-1][0][word] = 4 if word == 22 else int(parts[-1][0][word]) ^ 1
                store_parts(files, parts)
            reject('mixed_'+name, mixed, 'cross-rank metadata consistency')
    def phase(files, parts):
        v = parts[0][1]
        v[np.argmax(abs(v))] *= -1
        assert abs(sum(np.vdot(v, v).real for _, v in parts)-1) <= 1e-8
        store_parts(files, parts)
    reject('checksum', phase, 'payload checksum')
    def nonfinite(files, parts):
        parts[0][1][0] = np.nan
        store_parts(files, parts, True)
    reject('nan', nonfinite, 'finite vector validation')
    for delta, name in [(5e-9, 'norm_inside'), (2e-8, 'norm_outside')]:
        path = import_copy(root, name, source)
        files = checkpoint_files(path)
        parts = [(list(h), v.copy()) for h, v in map(read_kondo_checkpoint, files)]
        total = sum(np.vdot(v, v).real for _, v in parts)
        local = np.vdot(parts[0][1], parts[0][1]).real
        parts[0][1][:] *= np.sqrt((local+1+delta-total)/local)
        store_parts(files, parts, True)
        run_case(path, hphi, name, launcher, env,
                 None if name == 'norm_inside' else 'checkpoint global norm validation failed')
        print(name, 'independent checksum refreshed', flush=True)


def io_failures(root, source, hphi, launcher, env):
    for stage in ('write', 'publish'):
        path = import_copy(root, 'io_'+stage, source)
        for file in checkpoint_files(path):
            file.unlink()
        (path/'calc.def').write_text((path/'calc.def').read_text().replace('InputEigenVec 1', 'OutputEigenVec 1'))
        victim = path/'output'/f'zvo_eigenvec_0_rank_{mpi_size(launcher)-1}.dat'
        temporary = Path(str(victim)+'.part')
        if stage == 'write':
            (path/'output'/'regular').write_text('not a directory')
            temporary.symlink_to('regular/payload')
        else:
            victim.mkdir()
            (victim/'sentinel').write_text('retain')
        run_case(path, hphi, stage, launcher, env, 'checkpoint '+stage+' failed')
        if stage == 'write':
            assert temporary.is_symlink()
            temporary.unlink()
        else:
            assert (victim/'sentinel').read_text() == 'retain' and not temporary.exists()
            (victim/'sentinel').unlink()
            victim.rmdir()
        assert not list((path/'output').glob('*.part'))
        run_case(path, hphi, 'retry', launcher, env)
        parts = [read_kondo_checkpoint(file) for file in checkpoint_files(path)]
        assert all(tuple(h[26:28]) == checkpoint_digest(parts) for h, _ in parts)
        print('rank-local', stage, 'failure and retry passed', flush=True)


def mask_only(root, hphi, launcher, env):
    # Fully polarized L1/C1: either local-site assignment has raw/sector dim 1
    # and representative word 0101; only the explicit mask distinguishes them.
    case = Case('Kondo', 1, 'block', 1, 2, 0)
    source = expert_case(root/'mask_source', case, method=3, empty=True,
                         options={'CalcMod': {'OutputEigenVec': 1}})
    run_case(source, hphi, 'save', launcher, env)
    target = import_copy(root, 'mask_other_physical_space', source)
    from symmetry_spingc_common import definition
    definition(target, 'loc.def', [(0, 0), (1, 1)], 1)
    # First save the other real space, proving every legacy identity word is equal.
    other = expert_case(root/'mask_other_save', case, method=3, empty=True,
                        options={'CalcMod': {'OutputEigenVec': 1}})
    definition(other, 'loc.def', [(0, 0), (1, 1)], 1)
    run_case(other, hphi, 'save', launcher, env)
    for a, b in zip(checkpoint_files(source), checkpoint_files(other)):
        ha, _ = read_kondo_checkpoint(a)
        hb, _ = read_kondo_checkpoint(b)
        assert np.array_equal(ha[:21], hb[:21]) and ha[28] != hb[28]
    text = run_case(target, hphi, 'mask_reject', launcher, env,
                    'checkpoint header / sector / layout validation failed')
    assert 'local site mask differs' in text
    print('same dimension/representative/order mask-only rejection passed', flush=True)


def reject_non_kondo_v2(root, hphi, launcher, env):
    from symmetry_spingc_common import write_case, read_checkpoint
    from symmetry_spingc_observables import fixture
    nsite, permutations, characters, families, _, _, _ = fixture('D')
    source = root/'non_kondo_v1'
    write_case(source, nsite, permutations, characters, families, 3,
               {'CalcMod': {'OutputEigenVec': 1}})
    run_case(source, hphi, 'save_v1', launcher, env)
    target = import_copy(root, 'non_kondo_v2', source)
    for file in checkpoint_files(target):
        header, vector = read_checkpoint(file)
        header = list(header)+[0, 0]
        header[1] = 2
        file.write_bytes(struct.pack('<30Q', *header)+vector.astype('<c16').tobytes())
    text = run_case(target, hphi, 'reject_v2', launcher, env,
                    'checkpoint header / sector / layout validation failed')
    assert 'version differs' in text
    print('v2 non-Kondo rejected', flush=True)


def main():
    hphi = Path(sys.argv[1]).resolve()
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = mpi_size(launcher)
    assert ranks in (1, 2, 4, 16)
    root = Path(tempfile.mkdtemp(prefix='symmetry_kondo_checkpoint_', dir='.')).resolve()
    print('artifacts:', root, flush=True)
    for layout in ('distributed', 'replicated'):
        env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout)
        for model in ('Kondo', 'KondoNConserved', 'KondoGC'):
            case = Case(model, 3, 'block', None if model == 'KondoGC' else 2,
                        1 if model == 'Kondo' else None, 1)
            path = expert_case(root/(model+'_'+layout), case, method=3,
                               options={'CalcMod': {'OutputEigenVec': 1}})
            run_case(path, hphi, 'cg_save', launcher, env)
            header_check(path, case, ranks, layout)
            print(model, layout, 'v2 header passed', flush=True)
            imports(root, path, case, hphi, launcher, env)
            if model == 'Kondo':
                negative = root/('negative_'+layout)
                negative.mkdir()
                corruption_checks(negative, path, hphi, launcher, env)
                io_failures(negative, path, hphi, launcher, env)
        mask_root = root/('mask_'+layout)
        mask_root.mkdir()
        mask_only(mask_root, hphi, launcher, env)
        reject_non_kondo_v2(mask_root, hphi, launcher, env)


if __name__ == '__main__':
    main()
