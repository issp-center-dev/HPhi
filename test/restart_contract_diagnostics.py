"""Observable diagnostics and v3/legacy compatibility, without MPI redistribution."""
import shutil
import numpy as np
import restart_common as rc
import restart_import as ri
import restart_te_contract as te


def clone(source, target):
    shutil.copytree(source, target)
    return target


def fallback(seed, root, sector=False):
    for mode in (2, 3):
        path = clone(seed, root / ('import_' + str(mode)))
        assert rc.manifest(path) is None
        rc.run(path, 'old_history')
        for name in te.series(path):
            p = path / 'output' / name
            p.write_bytes(b'# OLD_HISTORY\n' + p.read_bytes())
        inputs = {p.name: p.read_bytes() for p in (path / 'output').glob('seed_eigenvec_0_rank_*.dat')}
        rc.set_key(path, 'calcmod.def', 'ReStart', mode)
        rc.set_key(path, 'modpara.def', 'Lanczos_max', 1)
        log = rc.run(path, 'fallback')
        assert 'SpectrumVec is imported at step 0 and existing TE time-series files are replaced' in log
        assert all(b'OLD_HISTORY' not in data for data in te.series(path).values())
        rows = rc.rows(path, ('zvo_' if sector else '') + 'SS.dat')
        assert rows.shape[0] == 1 and rows[0, -1] == 0
        if mode == 3:
            assert rc.manifest(path) is None
        else:
            assert rc.manifest(path)[1][0]['step'] == 1
        assert inputs == {p.name: p.read_bytes() for p in (path / 'output').glob('seed_eigenvec_0_rank_*.dat')}


def main():
    root = rc.fresh_dir('restart_contract_diagnostics').resolve()
    cg = root / 'cg'
    rc.generate(cg, ri.CG)
    rc.set_key(cg, 'calcmod.def', 'ReStart', 1)
    rc.set_key(cg, 'calcmod.def', 'OutputEigenVec', 1)
    rc.run(cg, 'cg')
    legacy = root / 'legacy_cg'
    rc.generate(legacy, ri.CG)
    for p in (cg / 'output').glob('zvo_eigenvec_0_rank_*.dat'):
        shutil.copyfile(p, legacy / 'output' / p.name.replace('zvo_eigenvec_0', 'tmpvec_set0'))
    rc.set_key(legacy, 'calcmod.def', 'ReStart', 2)
    log = rc.run(legacy, 'legacy')
    assert 'legacy CG continuation loads only the states, resets search directions' in log
    assert int(np.loadtxt(legacy / 'output/zvo_Lanczos_Step.dat', skiprows=1, ndmin=2)[0, 0]) == 1

    raw = te.build_time_table_case(root / 'raw', te.TIMES)
    # Ordinary raw CG->TE import remains legacy and must not become a rejection.
    rc.run(raw, 'legacy_import')
    legacy_te = clone(raw, root / 'legacy_te')
    # A shorter TE export contains a valid next step for continuation.
    rc.set_key(legacy_te, 'modpara.def', 'Lanczos_max', 6)
    rc.run(legacy_te, 'short_export')
    for p in (legacy_te / 'output').glob('zvo_eigenvec_0_rank_*.dat'):
        shutil.copyfile(p, p.with_name(p.name.replace('zvo_', 'seed_')))
    rc.set_key(legacy_te, 'modpara.def', 'Lanczos_max', 8)
    rc.set_key(legacy_te, 'calcmod.def', 'ReStart', 3)
    log = rc.run(legacy_te, 'legacy_continuation')
    assert 'legacy TE continuation restores only the step' in log
    assert 'existing TE time-series files are replaced' not in log

    ri.copy_state(cg, raw, 'seed_eigenvec_0')
    raw_cases = root / 'raw_imports'; raw_cases.mkdir()
    fallback(raw, raw_cases)
    stopped = clone(raw, root / 'stopped')
    rc.set_key(stopped, 'calcmod.def', 'ReStart', 1)
    rc.run(stopped, 'stop', binary=rc.probe(), env={'HPHI_TEST_RESTART_ABORT_AFTER': '0:6'})
    assert rc.manifest(stopped)[1][0]['step'] == 6
    bad = clone(stopped, root / 'bad_boundary')
    rc.set_key(bad, 'calcmod.def', 'ReStart', 2)
    times = list(te.TIMES); times[5] = .201
    te.write_times(bad, times)
    before = te.series(bad), te.checkpoints(bad)
    rc.run(bad, 'reject', fail='time grid changed')
    assert before == (te.series(bad), te.checkpoints(bad))
    for change_h in (False, True):
        case = clone(stopped, root / ('write_' + str(change_h)))
        rc.set_key(case, 'calcmod.def', 'ReStart', 2)
        if change_h:
            p = case / 'exchange.def'
            lines = p.read_text().splitlines()
            for i in range(5, len(lines)):
                row = lines[i].split()
                if row:
                    row[-1] = str(float(row[-1]) * 1.4)
                    lines[i] = ' '.join(row)
            p.write_text('\n'.join(lines) + '\n')
        before = te.checkpoints(case)
        log = rc.run(case, 'continue')
        assert 'v3 continuation does not verify the Hamiltonian or propagation history' in log
        assert rc.manifest(case)[1][0]['step'] == 8 and te.checkpoints(case) != before

    import symmetry_general_terms as fixture
    import symmetry_te
    fixture.ROOT = root / 'fixture'; fixture.ROOT.mkdir()

    def sector_case(path, model, length, momentum, states, h, projector):
        seed = root / 'sector'
        symmetry_te.build_restart_seed(path, seed, 'distributed')
        te.configure(seed, te.TIMES)
        cases = root / 'sector_imports'; cases.mkdir()
        fallback(seed, cases, sector=True)
        names = seed / 'namelist.def'
        names.write_text(names.read_text().replace('TEOneBody times.def', 'TETwoBody drive.def'))
        table = '====\nNTimeSteps 8\n====\n====\n====\n'
        for time, amplitude in zip(te.TIMES, [0, .21, -.17, 0, .31, .44, .2, .1]):
            table += '{} {}\n'.format(time, length if amplitude else 0)
            if amplitude:
                table += ''.join('{} 0 {} 0 {} 0 {} 0 {} 0\n'.format(
                    i, i, (i+1) % length, (i+1) % length, amplitude) for i in range(length))
        (seed / 'drive.def').write_text(table)
        reference = clone(seed, root / 'dynamic_reference')
        log = rc.run(reference, 'reference')
        assert 'hamiltonian_changed=unknown' in log and 'hamiltonian_changed=yes' not in log
        case = clone(seed, root / 'dynamic_continued')
        rc.set_key(case, 'calcmod.def', 'ReStart', 1)
        rc.run(case, 'stop', binary=rc.probe(), env={'HPHI_TEST_RESTART_ABORT_AFTER': '0:6'})
        rc.set_key(case, 'calcmod.def', 'ReStart', 2)
        rc.set_key(case, 'calcmod.def', 'InputEigenVec', 0)
        log = rc.run(case, 'resume')
        assert 'hamiltonian_changed=unknown' in log and 'hamiltonian_changed=yes' not in log
        assert 'SpectrumVec is imported at step 0' not in log
        te.compare(reference, case)
        np.testing.assert_allclose(symmetry_te.vector(reference, 'zvo_eigenvec_final', 7, .3),
                                   symmetry_te.vector(case, 'zvo_eigenvec_final', 7, .3), atol=1e-12, rtol=0)

    fixture.prepare('Spin', 8, sector_test=sector_case, sector_momenta=(0,))
    print('Restart diagnostics: legacy/v3 import, mode-2 rejection/publication and dynamic sector TE PASS')


if __name__ == '__main__':
    main()
