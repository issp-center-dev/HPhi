"""SpinGC version-1 checkpoint identity, corruption guards, and CG imports.

The embedded canonical payloads below are verbatim output of HPhi commit
279ede8fdd9aacde001d706f6545369ad5c84e08, before SpinGC sector support.
Generation: expert Spin N=6, 2Sz=0, cyclic k=0, Ising=.37, Exchange=.23
on every oriented periodic bond; CG exct=1, initial_iv=-1, Lanczos_max=400,
LanczosEps=12, LanczosTarget=1, LargeValue=100, PreCG=0, OutputEigenVec=1.
One MPI rank, replicated/distributed layouts. No header or payload rewriting.
"""
import base64
import hashlib
import os
from pathlib import Path
import shlex
import shutil
import struct
import sys
import tempfile

import numpy as np
import symmetry_spingc_common as c
from symmetry_spingc_observables import fixture, parse_energy, rank_slice

# BASELINE_PAYLOADS are actual baseline output, not current-writer fixtures.
BASELINE_PAYLOADS = {
    'replicated': (
        'e6ffa61837f7a9aa2597c8ad1d6cde8131c909dbdf75e3d8951a36bce7953809',
        'SFBISVNWMQoBAAAAAAAAAAEAAAAAAAAAgAAAAAAAAAABAAAAAAAAAAYAAAAAAAAAAwAAAAAAAAADAAAAAAAAAAMA'
        'AAAAAAAAFAAAAAAAAAAEAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAAAAAAAAAA9ng41'
        'H172rQQAAAAAAAAAkqKlryTbRHCyn6JUFEHMLvM/3x7thSrwB7buDV7o2xIDAAAAAAAAAAAAAAAAAAAACQAAAAAA'
        'AAAAAAAAAAAAAOBsbzEObsHM4GxvMQ5uwcyafSGeksbBPzGqqWPekqg/yvWcBnU41L+RKw58GfS7v/HdcAJ0ONS/'
        'M9I4Nh/0u78q1uSZRkfqP7PxPJgCKtI/'),
    'distributed': (
        '99b4e96ff5c4ebc12c2e0614534353ab73cd7d485fe533bfcd81e189c28e8e30',
        'SFBISVNWMQoBAAAAAAAAAAEAAAAAAAAAgAAAAAAAAAABAAAAAAAAAAYAAAAAAAAAAwAAAAAAAAADAAAAAAAAAAMA'
        'AAAAAAAAFAAAAAAAAAAEAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAABAAAAAAAAAAAAAAAAAAAABAAAAAAAAAA9ng41'
        'H172rQQAAAAAAAAAkqKlryTbRHCyn6JUFEHMLvM/3x7thSrwB7buDV7o2xIDAAAAAAAAAAAAAAAAAAAACQAAAAAA'
        'AAAAAAAAAAAAAOBsbzEObsHM4GxvMQ5uwcyafSGeksbBPzGqqWPekqg/yvWcBnU41L+RKw58GfS7v/HdcAJ0ONS/'
        'M9I4Nh/0u78q1uSZRkfqP7PxPJgCKtI/'),
}


def checkpoint_files(path):
    return sorted((path/'output').glob('zvo_eigenvec_0_rank_*.dat'),
                  key=lambda p: int(p.stem.rsplit('_', 1)[1]))


def store_parts(paths, parts, refresh_digest=False):
    digest = c.checkpoint_digest(parts) if refresh_digest else None
    for path, (header, vector) in zip(paths, parts):
        header = list(header)
        if digest is not None:
            header[26:28] = digest
        path.write_bytes(struct.pack('<28Q', *header)
                         + np.asarray(vector, dtype='<c16').tobytes())


def assert_digest(parts):
    digest = c.checkpoint_digest(parts)
    assert all(tuple(header[26:28]) == digest for header, _ in parts)


def import_copy(root, name, source):
    path = root/name
    path.mkdir()
    for item in source.glob('*.def'):
        shutil.copy2(item, path/item.name)
    (path/'output').mkdir()
    for item in checkpoint_files(source):
        shutil.copy2(item, path/'output'/item.name)
    calc = (path/'calc.def').read_text().replace('OutputEigenVec 1\n', '')
    (path/'calc.def').write_text(calc+'InputEigenVec 1\n')
    return path


def check_vector(path, nsite, basis, hamiltonian, ranks, layout):
    files = checkpoint_files(path)
    assert len(files) == ranks
    parts = [c.read_checkpoint(p) for p in files]
    dimension = basis.shape[1]
    for rank, (header, vector) in enumerate(parts):
        offset, size = rank_slice(dimension, ranks, rank)
        assert header[4:10] == (4, nsite, 0, 0, 0, 1 << nsite)
        assert header[10:16] == (dimension, ranks, rank,
                                int(layout == 'distributed'), offset, size)
        assert header[22:24] == (3, 0) and header[25] == 0
        assert len(vector) == size
    assert_digest(parts)
    vector = c.join_rank_vectors(files, dimension)
    assert abs(np.vdot(vector, vector).real-1) <= 1e-8
    energy, _ = parse_energy(path)
    assert np.isfinite(energy)
    np.testing.assert_allclose(energy, np.linalg.eigvalsh(hamiltonian)[0],
                               atol=3e-8, rtol=0)
    residual = np.linalg.norm(hamiltonian@vector-energy*vector)
    assert residual/max(1, np.linalg.norm(hamiltonian, 2)) <= 1e-8
    # Reconstruct the raw state with the independently constructed basis.
    raw = basis@vector
    assert np.isfinite(raw).all()
    np.testing.assert_allclose(basis.conj().T@raw, vector, atol=1e-12, rtol=0)
    return parts, vector, energy


def corruption_checks(root, source, hphi, launcher, env):
    def reject(name, mutate, diagnostic):
        path = import_copy(root, name, source)
        files = checkpoint_files(path)
        parts = [(list(h), v.copy()) for h, v in map(c.read_checkpoint, files)]
        mutate(files, parts)
        text = c.run_case(path, hphi, name, launcher, env,
                          'checkpoint '+diagnostic+' failed')
        assert not (path/'output/zvo_energy.dat').exists()
        print('rejected {}: {}'.format(name, diagnostic), flush=True)
        return text

    for word, name in [(4, 'model'), (6, 'nup'), (7, 'ndown'), (8, 'ne'),
                       (9, 'full_dim'), (11, 'rank_count'), (12, 'rank'),
                       (13, 'layout'), (14, 'offset'), (15, 'local_dim'),
                       (16, 'sector'), (20, 'order')]:
        def mutate(files, parts, word=word):
            parts[-1][0][word] ^= 1
            store_parts(files, parts)
        reject(name, mutate, 'header / sector / layout validation')

    for name, change in [('truncated_header', lambda b: b[:32]),
                         ('truncated_payload', lambda b: b[:-1]),
                         ('trailing_payload', lambda b: b+b'x')]:
        # Rank zero always owns payload, including the D np16 empty-rank case.
        reject(name, lambda files, parts, change=change:
               files[0].write_bytes(change(files[0].read_bytes())),
               'header / sector / layout validation' if name == 'truncated_header'
               else 'payload length')
    reject('missing_rank', lambda files, parts: files[-1].unlink(),
           'header / sector / layout validation')

    # D k=1 and k=5 have the same dimension: reject the actual character
    # change, independently of size and manually damaged header tests.
    path = import_copy(root, 'other_sector', source)
    lines = (path/'group.def').read_text().splitlines()
    for index in range(5, 11):
        group, real, imag = lines[index].split()
        lines[index] = '{} {} {}'.format(group, real, -float(imag))
    (path/'group.def').write_text('\n'.join(lines)+'\n')
    text = c.run_case(path, hphi, 'other_sector', launcher, env,
                      'checkpoint header / sector / layout validation failed')
    assert 'group/character digest differs' in text
    path = import_copy(root, 'other_layout', source)
    other_layout = ('replicated' if env['HPHI_SYMMETRY_BASIS_LAYOUT'] == 'distributed'
                    else 'distributed')
    text = c.run_case(path, hphi, 'other_layout', launcher,
                      dict(env, HPHI_SYMMETRY_BASIS_LAYOUT=other_layout),
                      'checkpoint header / sector / layout validation failed')
    assert 'basis layout differs' in text

    def phase(files, parts):
        vector = parts[0][1]
        index = int(np.argmax(abs(vector)))
        assert abs(vector[index]) > 1e-6
        old_digest = c.checkpoint_digest(parts)
        vector[index] *= -1  # exact phase change preserves norm, changes bits
        norm = sum(np.vdot(v, v).real for _, v in parts)
        assert np.isfinite(norm) and abs(norm-1) <= 1e-8
        assert c.checkpoint_digest(parts) != old_digest
        store_parts(files, parts)
    text = reject('phase_checksum', phase, 'payload checksum')
    assert 'norm validation failed' not in text and 'finite vector validation failed' not in text

    for value, name in [(float('nan'), 'nan'), (float('inf'), 'inf')]:
        def nonfinite(files, parts, value=value):
            parts[0][1][0] = complex(value, 0)
            store_parts(files, parts, refresh_digest=True)
        text = reject(name, nonfinite, 'finite vector validation')
        assert 'payload checksum failed' not in text

    for delta, name in [(5e-9, 'norm_inside'), (2e-8, 'norm_outside')]:
        path = import_copy(root, name, source)
        files = checkpoint_files(path)
        parts = [(list(h), v.copy()) for h, v in map(c.read_checkpoint, files)]
        total = sum(np.vdot(v, v).real for _, v in parts)
        local = np.vdot(parts[0][1], parts[0][1]).real
        assert local > 0
        parts[0][1][:] *= np.sqrt((local+1+delta-total)/local)
        norm = sum(np.vdot(v, v).real for _, v in parts)
        assert (abs(norm-1) <= 1e-8) == (name == 'norm_inside')
        store_parts(files, parts, refresh_digest=True)
        assert_digest([c.read_checkpoint(p) for p in files])
        expected = None if name == 'norm_inside' else 'checkpoint global norm validation failed'
        text = c.run_case(path, hphi, name, launcher, env, expected)
        assert 'payload checksum failed' not in text
        print('{} norm2={:.17g}'.format(name, norm), flush=True)

    if c.mpi_size(launcher) > 1:
        for word, name in [(21, 'source_hamiltonian'), (23, 'source_state'),
                           (24, 'source_step'), (25, 'source_time')]:
            def mixed(files, parts, word=word):
                parts[-1][0][word] ^= 1
                store_parts(files, parts)
            reject('mixed_'+name, mixed, 'cross-rank metadata consistency')


def io_failures(root, source, hphi, launcher, env):
    original = {p: p.read_bytes() for p in checkpoint_files(source)}
    for stage in ('write', 'publish'):
        path = import_copy(root, 'io_'+stage, source)
        # Clear only these fixture copies; the successful source set is intact.
        for p in checkpoint_files(path):
            p.unlink()
        calc = (path/'calc.def').read_text().replace('InputEigenVec 1', 'OutputEigenVec 1')
        (path/'calc.def').write_text(calc)
        victim = path/'output'/('zvo_eigenvec_0_rank_{}.dat'.format(c.mpi_size(launcher)-1))
        temporary = Path(str(victim)+'.part')
        if stage == 'write':
            blocked = path/'output'/'regular_file_parent'
            blocked.write_text('not a directory\n')
            # fopen follows this link and fails with ENOTDIR on just one rank.
            temporary.symlink_to('regular_file_parent/payload')
        else:
            victim.mkdir()
            (victim/'sentinel').write_text('nonempty directory\n')
        c.run_case(path, hphi, 'fail_'+stage, launcher, env,
                   'checkpoint '+stage+' failed')
        if stage == 'write':
            assert temporary.is_symlink() and blocked.is_file()
            temporary.unlink()
        else:
            assert (victim/'sentinel').read_text() == 'nonempty directory\n'
            assert not temporary.exists()
            (victim/'sentinel').unlink()
            victim.rmdir()
        assert not list((path/'output').glob('*.part'))
        c.run_case(path, hphi, 'retry_'+stage, launcher, env)
        assert len(checkpoint_files(path)) == c.mpi_size(launcher)
        assert_digest([c.read_checkpoint(p) for p in checkpoint_files(path)])
        assert all(p.read_bytes() == data for p, data in original.items())
        print('rank-local {} failure and retry passed'.format(stage), flush=True)


def exercise_fixture(root, label, layout, hphi, probe, launcher):
    nsite, permutations, characters, families, basis, _, hamiltonian = fixture(label)
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout)
    path = root/(label+'_'+layout)
    c.write_case(path, nsite, permutations, characters, families, 3,
                 {'CalcMod': {'OutputEigenVec': 1}, 'ModPara': {'LanczosEps': 18}})
    text = c.run_case(path, probe, 'poison_write', launcher,
                      dict(env, HPHI_TEST_SPINGC_POISON_FIXED='1'))
    assert 'SpinGCProbe poisoned fixed fields after basis setup.' in text
    parts, vector, energy = check_vector(path, nsite, basis, hamiltonian,
                                         c.mpi_size(launcher), layout)
    if label == 'D':
        assert basis.shape[1] == 9
        assert sum(len(v) == 0 for _, v in parts) == max(0, c.mpi_size(launcher)-9)
    reload = import_copy(root, label+'_'+layout+'_reload', path)
    with (reload/'calc.def').open('a') as handle:
        handle.write('OutputEigenVec 1\n')
    c.run_case(reload, hphi, 'cg_import', launcher, env)
    for old, new in zip(checkpoint_files(path), checkpoint_files(reload)):
        assert old.read_bytes() == new.read_bytes(), 'CG import changed checkpoint bytes'
    np.testing.assert_allclose(parse_energy(reload)[0], energy, atol=2e-12, rtol=0)

    quench = import_copy(root, label+'_'+layout+'_quench', path)
    # Trans enters H with a minus sign; add +.21 Sy to the physical H.
    rows = [(i, a, j, b, re-(.21*c.Sy[a, b]).real,
             im-(.21*c.Sy[a, b]).imag) for i, a, j, b, re, im in families['Trans']]
    c.definition(quench, 'Trans.def', rows)
    text = c.run_case(quench, hphi, 'quench_import', launcher, env)
    source_digest = parts[0][0][21]
    assert 'source_hamiltonian={:016x} hamiltonian_changed=yes'.format(source_digest) in text
    manifest = dict(line.split('=', 1) for line in
                    (quench/'output/zvo_symmetry_sector.dat').read_text().splitlines())
    assert manifest['hamiltonian_digest'].startswith('hphi-parsed-hamiltonian-fnv1a64-v3:')
    assert int(manifest['hamiltonian_digest'].split(':')[1], 16) != source_digest
    changed_h = hamiltonian+.21*(basis.conj().T@sum(c.spin_operators(nsite)['y'])@basis)
    expected = np.vdot(vector, changed_h@vector)
    assert np.isfinite(expected) and abs(expected.imag) <= 1e-12
    np.testing.assert_allclose(parse_energy(quench)[0], expected.real, atol=1e-8, rtol=0)
    for old, new in zip(checkpoint_files(path), checkpoint_files(quench)):
        assert old.read_bytes() == new.read_bytes(), 'import must preserve source payload'
    if label == 'D':
        negative = root/('negative_'+layout)
        negative.mkdir()
        if c.mpi_size(launcher) > 1:
            # Two individually valid saved sets, identical basis but different H.
            # Only one rank is replaced: cross-rank provenance must reject it.
            with (quench/'calc.def').open('a') as handle:
                handle.write('OutputEigenVec 1\n')
            c.run_case(quench, hphi, 'save_quenched_state', launcher, env)
            changed = [c.read_checkpoint(p) for p in checkpoint_files(quench)]
            assert_digest(changed)
            assert all(h[21] != source_digest for h, _ in changed)
            mixed = import_copy(negative, 'mixed_actual_saved_sets', path)
            victim = checkpoint_files(mixed)[-1]
            shutil.copy2(checkpoint_files(quench)[-1], victim)
            c.run_case(mixed, hphi, 'mixed_actual_saved_sets', launcher, env,
                       'checkpoint cross-rank metadata consistency failed')
        corruption_checks(negative, path, hphi, launcher, env)
        io_failures(negative, path, hphi, launcher, env)
    print('{} {} CG checkpoint/import/quench passed'.format(label, layout), flush=True)


def old_payload_compatibility(root, hphi):
    nsite = 6
    permutations = [[(i+g) % nsite for i in range(nsite)] for g in range(nsite)]
    families = {key: [(i, (i+1) % nsite, value) for i in range(nsite)]
                for key, value in [('Ising', .37), ('Exchange', .23)]}
    for layout, (checksum, encoded) in BASELINE_PAYLOADS.items():
        payload = base64.b64decode(encoded)
        assert hashlib.sha256(payload).hexdigest() == checksum
        path = root/('old_canonical_'+layout)
        c.write_case(path, nsite, permutations, [1]*nsite, families, 3,
                     {'CalcMod': {'CalcModel': 1, 'InputEigenVec': 1},
                      'ModPara': {'2Sz': 0}})
        (path/'output').mkdir()
        file = path/'output/zvo_eigenvec_0_rank_0.dat'
        file.write_bytes(payload)
        header, _ = c.read_checkpoint(file)
        assert header[4:9] == (1, 6, 3, 3, 3)
        text = c.run_case(path, hphi, 'baseline_import', [],
                          dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout))
        assert 'hamiltonian_changed=no' in text
        manifest = (path/'output/zvo_symmetry_sector.dat').read_text()
        assert 'hamiltonian_digest=hphi-parsed-hamiltonian-fnv1a64-v2:12dbe85e0deeb607' in manifest
        np.testing.assert_allclose(parse_energy(path)[0], -.8615361366332042, atol=2e-12, rtol=0)
        assert file.read_bytes() == payload
        print('actual baseline canonical {} payload accepted'.format(layout), flush=True)


def main():
    hphi = Path(sys.argv[1]).resolve()
    probe = hphi.with_name('unittest_symmetry_spingc_probe')
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = c.mpi_size(launcher)
    assert ranks in (1, 4, 16), 'checkpoint suite requires explicit np1, np4, or np16'
    root = Path(tempfile.mkdtemp(prefix='symmetry_spingc_checkpoint_', dir='.')).resolve()
    print('artifacts: {}'.format(root), flush=True)
    for layout in ('distributed', 'replicated') if ranks <= 4 else ('distributed',):
        for label in ('B0', 'B1', 'D'):
            exercise_fixture(root, label, layout, hphi, probe, launcher)
    if ranks == 1:
        old_payload_compatibility(root, hphi)


if __name__ == '__main__':
    main()
