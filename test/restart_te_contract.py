"""Same-rank TE continuation, saved clocks, recovery and collective rejection."""
import shutil
import struct

import numpy as np
import restart_common as rc

TIMES = [0, .03, .09, .09, .13, .2, .25, .3]


def write_times(path, times):
    (path / 'times.def').write_text('====\nNTimeSteps {}\n====\n====\n====\n'.format(len(times)) +
                                  ''.join('{} 0\n'.format(t) for t in times))


def build_time_table_case(root, times):
    """Reusable raw fixture (Task 21): CG seed followed by a quenched TE grid."""
    rc.generate(root, 'L=8\nmodel="Spin"\nmethod="CG"\nlattice="chain"\nJ=1\n2Sz=0\nexct=1\noutputmode="None"\n')
    rc.set_key(root, 'calcmod.def', 'OutputEigenVec', 1)
    rc.run(root, 'cg')
    for p in (root / 'output').glob('zvo_eigenvec_0_rank_*.dat'):
        shutil.copyfile(p, p.with_name(p.name.replace('zvo_', 'seed_')))
    p = root / 'exchange.def'
    lines = p.read_text().splitlines()
    for i in range(5, len(lines)):
        row = lines[i].split()
        if row:
            row[-1] = str(float(row[-1]) * 1.3)
            lines[i] = ' '.join(row)
    p.write_text('\n'.join(lines) + '\n')
    configure(root, times)
    return root


def configure(path, times):
    for key, value in [('CalcType', 4), ('InputEigenVec', 1), ('OutputEigenVec', 1),
                       ('ReStart', 0), ('OutputGreenFormat', 1)]:
        rc.set_key(path, 'calcmod.def', key, value)
    for key, value in [('Lanczos_max', len(times)), ('ExpandCoef', 8),
                       ('OutputInterval', 2), ('ExpecInterval', 1)]:
        rc.set_key(path, 'modpara.def', key, value)
    names = path / 'namelist.def'
    names.write_text(''.join(line + '\n' for line in names.read_text().splitlines() if not line.lstrip().startswith(('SpectrumVec', 'OneBodyG'))) + 'SpectrumVec seed_eigenvec_0\nTEOneBody times.def\nOneBodyG greenone.def\n')
    (path / 'greenone.def').write_text('====\nNCisAjs 2\n====\n====\n====\n0 0 0 0\n1 0 1 0\n')
    write_times(path, times)


def laser(path, length):
    # DC Laser fixture from symmetry_te_dynamic.py (mode 2, field .4).
    p = path / 'namelist.def'
    p.write_text(p.read_text().replace('TEOneBody times.def', 'Laser laser.def'))
    (path / 'laser.def').write_text('====\nNLaser 9\n====\n====\n====\n' +
        ''.join('p{} {}\n'.format(i, v) for i, v in enumerate([2, .4, 1, 1, 0, length, 1, 1, 0])))
    rc.set_key(path, 'modpara.def', 'Tinit', .01)
    rc.set_key(path, 'modpara.def', 'TimeSlice', .035)


def build_raw_laser(root):
    rc.generate(root, 'W=2\nL=2\nmodel="Hubbard"\nmethod="CG"\nlattice="square"\nt=1\nU=4\nnelec=4\n2Sz=0\nexct=1\noutputmode="None"\n')
    rc.set_key(root, 'calcmod.def', 'OutputEigenVec', 1)
    rc.run(root, 'cg')
    for p in (root / 'output').glob('zvo_eigenvec_0_rank_*.dat'):
        shutil.copyfile(p, p.with_name(p.name.replace('zvo_', 'seed_')))
    configure(root, TIMES)
    laser(root, 2)
    # Raw TE diagonal kernels require allocated empty per-step term tables.
    p = root / "namelist.def"
    p.write_text(p.read_text() + "TEOneBody times.def\n")
    return root


def series(path):
    return {p.name: p.read_bytes() for p in (path / 'output').glob('*.dat')
            if p.name in ('SS.dat', 'Norm.dat', 'Flct.dat', 'zvo_SS.dat', 'zvo_Norm.dat', 'zvo_Flct.dat')
            or p.name.endswith('_te.dat')}


def checkpoints(path):
    return {p.name: p.read_bytes() for p in (path / 'output').glob('zvo_restart*')}


def raw_final_vector(path):
    """Read the unchanged 12-byte raw OutputEigenVec format as one vector."""
    files = sorted((path / 'output').glob('zvo_eigenvec_0_rank_*.dat'),
                   key=lambda p: int(p.stem.split('rank_')[1]))
    assert files, path
    pieces = []
    next_steps = set()
    for source in files:
        data = source.read_bytes()
        next_step, dim = struct.unpack_from('=iQ', data)
        assert len(data) == 12 + 16 * (dim + 1), (source, dim, len(data))
        next_steps.add(next_step)
        pieces.append(np.frombuffer(data, dtype='=c16', offset=12)[1:].copy())
    assert len(next_steps) == 1, next_steps
    return next_steps.pop(), np.concatenate(pieces)


def compare_raw_final_vectors(reference, actual):
    expected_step, expected = raw_final_vector(reference)
    actual_step, found = raw_final_vector(actual)
    assert expected.size > 0 and found.size > 0
    assert actual_step == expected_step and found.shape == expected.shape
    overlap = np.vdot(expected, found)
    assert abs(overlap) > 1e-12
    phase = overlap / abs(overlap)
    np.testing.assert_allclose(found / phase, expected, rtol=0, atol=1e-11)


def append_tail(path):
    for name in series(path):
        p = path / 'output' / name
        rows = [line.split() for line in p.read_text().splitlines() if line and not line.startswith('#')]
        assert rows, name
        row = rows[-1]
        row[0 if name.endswith('_te.dat') else -1] = '6'
        p.write_text(p.read_text() + ' '.join(row) + '\n')


def compare(reference, path):
    assert set(series(reference)) == set(series(path))
    for name in series(reference):
        x, y = rc.rows(reference, name), rc.rows(path, name)
        if name.endswith('_te.dat'):
            x = rc.sorted_by_key(x, tuple(range(x.shape[1] - 2)))
            y = rc.sorted_by_key(y, tuple(range(y.shape[1] - 2)))
            assert np.any(abs(x[:, -2:]) > 1e-8)
        np.testing.assert_allclose(y, x, atol=5e-12, rtol=0, err_msg=name)


def clone(path, name):
    target = path.with_name(name)
    shutil.copytree(path, target)
    return target


def exercise(seed, sector=False, layout=None, is_laser=False):
    env = {'HPHI_SYMMETRY_BASIS_LAYOUT': layout} if layout else {}
    run = lambda p, label, **kw: rc.run(p, label, env=env, **kw)
    ref = clone(seed, seed.name + '_ref')
    run(ref, 'reference')
    if sector:
        greens = rc.rows(ref, 'zvo_cisajs_te.dat')
        assert np.max(abs(greens[-2:, -2:] - greens[:2, -2:])) > 1e-6
    cont = clone(seed, seed.name + '_cont')
    rc.set_key(cont, 'calcmod.def', 'ReStart', 1)
    rc.run(cont, 'interrupted', binary=rc.probe(), env=dict(env, HPHI_TEST_RESTART_ABORT_AFTER='0:6'))
    published = rc.manifest(cont)
    assert published is not None, 'TE ReStart=1 did not publish a manifest'
    assert published[1][0]['step'] == 6
    files = rc.restart_files(cont, 'zvo', 0, 0, rc.manifest(cont)[1][0]['gen'])
    expected_time = .01 + (5 if sector else 6) * .035 if is_laser else TIMES[5]
    for p in files:
        h, _ = rc.v3_header(p)
        np.testing.assert_allclose(struct.unpack('<d', struct.pack('<Q', h[25]))[0], expected_time, atol=1e-15, rtol=0)
    saved_payloads = [rc.v3_header(p)[1].tobytes() for p in files]
    saved_times = [rc.v3_header(p)[0][25] for p in files]
    saved = checkpoints(cont)
    if not sector and not is_laser:
        import restart_observable_failures as output_failures
        for aggregate in (False, True):
            output_failures.exercise(cont, "one", aggregate, 6, te=True)
    prefix = series(cont)
    rc.set_key(cont, 'calcmod.def', 'InputEigenVec', 0)
    rc.set_key(cont, 'calcmod.def', 'ReStart', 3)
    append_tail(cont)
    stopped = clone(cont, seed.name + '_stopped')
    run(cont, 'continued')
    compare(ref, cont)
    if not sector:
        compare_raw_final_vectors(ref, cont)
    assert checkpoints(cont) == saved
    if sector:
        record = (cont / 'output/zvo_symmetry_sector.dat').read_text()
        assert 'continuation=yes' in record and 'source_step=6' in record
    # Equal and lowered limits; both read modes preserve saved checkpoint bytes.
    for restart in (2, 3):
        noop = clone(stopped, seed.name + '_noop' + str(restart))
        rc.set_key(noop, 'calcmod.def', 'ReStart', restart)
        for attempt, limit in enumerate((6, 4, 4)):
            rc.set_key(noop, 'modpara.def', 'Lanczos_max', limit)
            if is_laser and not sector and limit == 4:
                write_times(noop, TIMES[:4])  # allocation table does not constrain saved Laser STEP
            assert 'no remaining steps' in run(noop, 'noop{}_{}_{}'.format(restart, limit, attempt))
            assert series(noop) == prefix and checkpoints(noop) == saved
            if sector:
                assert len(list((noop / 'output').glob('zvo_eigenvec_final_rank_*.dat'))) == len(files)
                for rank, p in enumerate(sorted((noop / 'output').glob('zvo_eigenvec_final_rank_*.dat'),
                                               key=lambda p: int(p.stem.split('rank_')[1]))):
                    data = p.read_bytes()
                    assert data[224:] == saved_payloads[rank]
                    assert struct.unpack_from('<Q', data, 25*8)[0] == saved_times[rank]
                    assert struct.unpack_from('<Q', data, 24*8)[0] == 6
                    np.testing.assert_allclose(struct.unpack_from('<d', data, 25*8)[0], expected_time, atol=1e-15, rtol=0)
            else:
                for rank, p in enumerate(sorted((noop / 'output').glob('zvo_eigenvec_6_rank_*.dat'),
                                               key=lambda p: int(p.stem.split('rank_')[1]))):
                    data = p.read_bytes()
                    assert struct.unpack_from('=i', data)[0] == 6
                    assert data[28:] == saved_payloads[rank]
                assert len(list((noop / 'output').glob('zvo_eigenvec_6_rank_*.dat'))) == len(files)
    # STEP 7 is not periodic (interval 2): this tests the terminal writer itself.
    terminal = clone(seed, seed.name + '_terminal')
    rc.set_key(terminal, 'modpara.def', 'Lanczos_max', 7)
    rc.set_key(terminal, 'calcmod.def', 'ReStart', 1)
    run(terminal, 'terminal')
    assert rc.manifest(terminal)[1][0]['step'] == 7
    rc.set_key(terminal, 'modpara.def', 'Lanczos_max', 8)
    rc.set_key(terminal, 'calcmod.def', 'ReStart', 3)
    rc.set_key(terminal, 'calcmod.def', 'InputEigenVec', 0)
    run(terminal, 'extension')
    compare(ref, terminal)
    # A missing v3 checkpoint retains the input-vector requirement, also in sector mode.
    fallback = clone(seed, seed.name + '_fallback')
    rc.set_key(fallback, 'calcmod.def', 'ReStart', 2)
    rc.set_key(fallback, 'calcmod.def', 'InputEigenVec', 0)
    run(fallback, 'fallback', fail='A file of Inputvector is not inputted')
    if sector:
        old_import = clone(seed, seed.name + '_old_import')
        rc.set_key(old_import, 'calcmod.def', 'ReStart', 3)
        run(old_import, 'old_import')
        compare(ref, old_import)
        record = (old_import / 'output/zvo_symmetry_sector.dat').read_text()
        assert 'source_method=3' in record and 'continuation=yes' not in record
    if not is_laser:
        for label, times, limit, failure in [
            ('grid', TIMES[:5] + [.201] + TIMES[6:], 8, 'time grid changed'),
            ('short', TIMES[:4], 4, 'time table shorter than the saved checkpoint')]:
            bad = clone(stopped, seed.name + '_' + label)
            write_times(bad, times)
            rc.set_key(bad, 'modpara.def', 'Lanczos_max', limit)
            before = series(bad)
            run(bad, label, fail=failure)
            assert series(bad) == before and checkpoints(bad) == saved
    print(seed.name + ': continuation/Green/tail/no-op/terminal/fallback PASS', flush=True)
    return stopped


def main():
    root = rc.fresh_dir('restart_te_contract')
    raw = build_time_table_case(root / 'raw_table', TIMES)
    stopped = exercise(raw)
    exercise(build_raw_laser(root / 'raw_laser'), is_laser=True)
    # Import the sector preparation helper without executing its numerical suite.
    import symmetry_te
    import symmetry_general_terms as fixture
    fixture.ROOT = root / 'sector_fixture'
    fixture.ROOT.mkdir()
    def sector_case(path, model, length, momentum, states, raw_h, projector):
        for layout in ('replicated', 'distributed'):
            for drive in ('table', 'laser'):
                target = root / ('sector_' + layout + '_' + drive)
                symmetry_te.build_restart_seed(path, target, layout)
                configure(target, TIMES)
                # Nonlocal hopping observes the actual quenched state evolution.
                (target / 'greenone.def').write_text('====\nNCisAjs 2\n====\n====\n====\n0 0 1 0\n1 0 0 0\n')
                if drive == 'laser':
                    laser(target, length)
                exercise(target, sector=True, layout=layout, is_laser=drive == 'laser')
    fixture.prepare('SpinlessFermion', 6, sector_test=sector_case, sector_momenta=(0,))
    # Mutate all rank headers and the manifest coherently, isolating STEP range validation.
    for value in (0, 2**31-1, 2**31):
        bad = clone(stopped, 'bad_step_' + str(value))
        m = bad / 'output/zvo_restart_manifest.dat'
        m.write_text(m.read_text().replace('step 6 ', 'step {} '.format(value)))
        row = rc.manifest(bad)[1][0]
        for p in rc.restart_files(bad, 'zvo', 0, 0, row['gen']):
            data = bytearray(p.read_bytes()); struct.pack_into('<Q', data, 24*8, value); p.write_bytes(data)
        before = series(bad)
        rc.run(bad, 'bad_step', fail='step')
        assert series(bad) == before
    # Positive publication interval is required even without OutputEigenVec.
    missing = clone(raw, 'missing_interval')
    rc.set_key(missing, 'calcmod.def', 'ReStart', 1)
    rc.set_key(missing, 'calcmod.def', 'OutputEigenVec', 0)
    rc.drop_key(missing, 'modpara.def', 'OutputInterval')
    rc.run(missing, 'missing_interval', mpi=False, fail='TE requires positive ExpandCoef, ExpecInterval and vector OutputInterval')
    blocked = clone(raw, 'blocked_output')
    (blocked / 'output/SS.dat').mkdir()
    rc.run(blocked, 'blocked_output', fail='failed to open TE output', timeout=60)
    # Small TIME isolates the absolute 1e-15 floor (old illustrative floor ~1e-12 is wrong).
    small_times = [0, 1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6]
    small = clone(raw, 'small_time')
    write_times(small, small_times)
    rc.set_key(small, 'calcmod.def', 'ReStart', 1)
    rc.run(small, 'interrupted', binary=rc.probe(), env={'HPHI_TEST_RESTART_ABORT_AFTER': '0:6'})
    rc.set_key(small, 'calcmod.def', 'ReStart', 3)
    for label, delta, fail in [('inside', .5e-15, None), ('outside', 2e-15, 'time grid changed')]:
        bad = clone(small, 'small_' + label)
        changed = small_times.copy(); changed[5] += delta
        write_times(bad, changed)
        rc.run(bad, label, fail=fail)
    print('TE STEP bounds, interval guard and small-TIME tolerance PASS', flush=True)


if __name__ == '__main__':
    main()
