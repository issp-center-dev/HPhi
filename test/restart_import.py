"""Import uses checkpoint identity, not the source solver clock."""
import shutil
import struct
import numpy as np
import restart_common as rc
import restart_te_contract as te

CG = 'L=8\nmodel="Spin"\nmethod="CG"\nlattice="chain"\nJ=1\n2Sz=0\nexct=1\nLanczos_max=100\noutputmode="None"\n'


def copy_state(source, target, base):
    row = rc.manifest(source)[1][0]
    (target / 'output').mkdir(exist_ok=True)
    files = []
    for p in rc.restart_files(source, 'zvo', 0, 0, row['gen']):
        q = target / 'output' / ('{}_rank_{}.dat'.format(base, p.stem.split('rank_')[1]))
        shutil.copyfile(p, q)
        files.append(q)
    return files


def clone(source, target):
    shutil.copytree(source, target)
    return target


def mutate(files, words):
    for p in files:
        data = bytearray(p.read_bytes())
        for word, value in words.items():
            struct.pack_into('<d' if isinstance(value, float) else '<Q', data, word*8, value)
        p.write_bytes(data)


def snapshot(path):
    return {p.name: p.read_bytes() for p in (path / 'output').glob('*.dat')
            if '_rank_' in p.name or '_manifest' in p.name or p.name.endswith('_te.dat') or
            p.name in ('SS.dat', 'Norm.dat', 'Flct.dat', 'zvo_SS.dat', 'zvo_Norm.dat', 'zvo_Flct.dat', 'zvo_Lanczos_Step.dat')}


def run(path, label, layout=None, **kwargs):
    return rc.run(path, label, env={'HPHI_SYMMETRY_BASIS_LAYOUT': layout} if layout else None,
                  timeout=60, **kwargs)


def check_clock(path, times, sector=False):
    ss = rc.rows(path, ('zvo_' if sector else '') + 'SS.dat')
    np.testing.assert_array_equal(ss[:, -1], np.arange(len(times)))
    np.testing.assert_allclose(ss[:, 0], times, rtol=0, atol=1e-12)
    green = rc.rows(path, 'zvo_cisajs_te.dat')
    np.testing.assert_array_equal(green[:, 0], np.repeat(np.arange(len(times)), 2))


def check_initial_vector(path, sector):
    if not sector:
        # Raw final OutputEigenVec intentionally overwrites label0. Isolate
        # the initial dt0 row so the final file still contains the import.
        path = clone(path, path.with_name(path.name + '_first_step'))
        rc.set_key(path, 'modpara.def', 'Lanczos_max', 1)
        run(path, 'first_step')
    for seed in (path / 'output').glob('seed_eigenvec_0_rank_*.dat'):
        rank = seed.stem.split('rank_')[1]
        _, expected = rc.v3_header(seed)
        data = (path / 'output' / ('zvo_eigenvec_0_rank_{}.dat'.format(rank))).read_bytes()
        actual = np.frombuffer(data[224 if sector else 28:], dtype='<c16')
        np.testing.assert_allclose(actual, expected, atol=3e-13, rtol=0)


def energy(path):
    text = (path / 'output/zvo_energy.dat').read_text()
    return float(next(line.split()[1] for line in text.splitlines() if line.strip().startswith('Energy ')))


def reject_mutation(source, root, label, words, diagnostic, layout=None, one_rank=False):
    bad = clone(source, root / label)
    files = sorted((bad / 'output').glob('seed_eigenvec_0_rank_*.dat'))
    mutate(files[-1:] if one_rank else files, words)
    before = snapshot(bad)
    run(bad, label, layout, fail=diagnostic)
    assert snapshot(bad) == before, label


def exercise(seed, root, layout=None):
    sector = layout is not None
    times = [0, .025, .05, .075, .1]
    # Save a real CG-v3 checkpoint in the destination basis/layout.
    cg = clone(seed, root / 'cg')
    rc.set_key(cg, 'calcmod.def', 'CalcType', 3)
    rc.set_key(cg, 'calcmod.def', 'InputEigenVec', 0)
    rc.set_key(cg, 'calcmod.def', 'ReStart', 1)
    run(cg, 'cg', layout)
    target = clone(seed, root / 'cg_to_te')
    te.configure(target, times)
    files = copy_state(cg, target, 'seed_eigenvec_0')
    # Every provenance word differs from the new solver. Oversized STEP is not
    # cast to int on import; defined FLAGS=1 is valid for a TE destination.
    provenance = {22: 3, 23: 6, 24: 2**63+5, 25: 12.5, 35: 9,
                  36: -1.25, 37: 2.5, 38: .125, 39: 1}
    mutate(files, provenance)
    # No manifest in seed: ReStart=3 still imports v3 rather than raw-classifying it.
    rc.set_key(target, 'calcmod.def', 'ReStart', 3)
    text = run(target, 'import', layout)
    expected = ['source_method=3', 'source_state=6', 'source_step=9223372036854775813',
                'source_time=12.5', 'source_nset=9', 'source_flags=1',
                'source_scalar_a=-1.25', 'source_scalar_b=2.5', 'source_scalar_c=0.125']
    assert 'Restart checkpoint imported' in text
    for item in expected:
        assert item in text, (item, text)
    if sector:
        metadata = (target / 'output/zvo_symmetry_sector.dat').read_text()
        for item in expected:
            assert item in metadata, (item, metadata)
    check_clock(target, times, sector)
    check_initial_vector(target, sector)
    np.testing.assert_allclose(rc.rows(target, ('zvo_' if sector else '')+'SS.dat')[0, 1],
                               energy(cg), atol=2e-8, rtol=0)
    # STEP=0 is likewise a valid provenance value; only continuation has bounds.
    zero = clone(target, root / 'zero_source_step')
    mutate(list((zero / 'output').glob('seed_eigenvec_0_rank_*.dat')), {24: 0})
    run(zero, 'zero_step', layout)
    check_clock(zero, times, sector)
    # Genuine TE-v3 -> TE import uses a different time origin and fresh indices.
    evolved = clone(target, root / 'te_source')
    rc.set_key(evolved, 'calcmod.def', 'ReStart', 1)
    run(evolved, 'save_te', layout)
    restarted = clone(seed, root / 'te_new_origin')
    new_times = [2, 2.03, 2.08, 2.12, 2.2]
    te.configure(restarted, new_times)
    copy_state(evolved, restarted, 'seed_eigenvec_0')
    text = run(restarted, 'new_origin', layout)
    assert 'source_method=4' in text and 'source_step=5' in text and 'source_time=0.1' in text
    check_clock(restarted, new_times, sector)
    check_initial_vector(restarted, sector)
    # Importing TE into CG evaluates the state without a solver iteration.
    evaluate = clone(seed, root / 'cg_evaluate')
    rc.set_key(evaluate, 'calcmod.def', 'InputEigenVec', 1)
    rc.set_key(evaluate, 'calcmod.def', 'OutputEigenVec', 1)
    copy_state(evolved, evaluate, 'zvo_eigenvec_0')
    text = run(evaluate, 'evaluate', layout)
    assert 'source_method=4' in text and rc.manifest(evaluate) is None
    expected_energy = rc.rows(evolved, ('zvo_' if sector else '')+'SS.dat')[-1, 1]
    np.testing.assert_allclose(energy(evaluate), expected_energy, atol=2e-8, rtol=0)
    for p in (evaluate / 'output').glob('zvo_eigenvec_0_rank_*.dat'):
        data = p.read_bytes()
        if sector:
            h = struct.unpack('<28Q', data[:224])
            assert h[1] == 1 and h[24] == 0 and len(data) == 224 + 16*h[15]
        else:
            step, dim = struct.unpack('=iQ', data[:12])
            assert step == 0 and len(data) == 12+16*(dim+1)
    for label, words, diagnostic in [
        ('phase', {2: 7}, 'phase convention'),
        ('model', {4: 999}, 'model differs'),
        ('component', {31: 1}, 'identity / unit validation'),
        ('flags', {39: 2}, 'flags differs'),
        ('nonfinite', {36: float('nan')}, 'identity / unit validation'),
        ('version', {1: 99}, 'import format / version'),
        ('basis_identity', {16 if sector else 33: 123}, 'digest differs')]:
        reject_mutation(target, root, label, words, diagnostic, layout)
    if len(files) > 1:
        reject_mutation(target, root, 'rank_provenance', {25: 17.0}, 'header scan', layout, True)
    mismatch = clone(target, root / 'source_count')
    rc.set_key(mismatch, 'modpara.def', 'RestartSourceRanks', len(files)+1)
    before = snapshot(mismatch)
    run(mismatch, 'ranks', layout, fail='RestartSourceRanks does not match')
    assert snapshot(mismatch) == before
    continuing = clone(evolved, root / 'continuation_count')
    rc.set_key(continuing, 'calcmod.def', 'ReStart', 3)
    rc.set_key(continuing, 'modpara.def', 'RestartSourceRanks', len(files)+1)
    series_name = ('zvo_' if sector else '') + 'SS.dat'
    expected_series = rc.rows(continuing, series_name).copy()
    checkpoint_before = {
        p.name: p.read_bytes()
        for p in (continuing / 'output').glob('zvo_restart*')
    }
    text = run(continuing, 'ranks', layout)
    assert 'no remaining steps' in text
    np.testing.assert_allclose(rc.rows(continuing, series_name), expected_series,
                               rtol=0, atol=0)
    assert {p.name: p.read_bytes()
            for p in (continuing / 'output').glob('zvo_restart*')} == checkpoint_before
    # A genuinely different destination process count, not a forged header.
    # Both named-state import and manifest continuation redistribute v3 files.
    if not sector and len(files) > 1:
        changed = clone(target, root / 'changed_np_import')
        text = run(changed, 'changed_np_import', mpi=False)
        assert 'Restart read: header files {}'.format(len(files)) in text
        check_clock(changed, times)
        continued = clone(evolved, root / 'changed_np_continuation')
        rc.set_key(continued, 'calcmod.def', 'ReStart', 3)
        expected_series = rc.rows(continued, 'SS.dat').copy()
        text = run(continued, 'changed_np_continuation', mpi=False)
        assert ('no remaining steps' in text and
                'Restart read: header files {}'.format(len(files)) in text)
        np.testing.assert_allclose(rc.rows(continued, 'SS.dat'), expected_series,
                                   rtol=0, atol=0)
    print(('sector '+layout if sector else 'raw') + ': CG/TE import, CG evaluation, provenance, identity and ranks PASS', flush=True)


def main():
    root = rc.fresh_dir('restart_import')
    seed = root / 'raw_seed'
    rc.generate(seed, CG)
    exercise(seed, root)
    import symmetry_te
    import symmetry_general_terms as fixture
    fixture.ROOT = root / 'sector_fixture'
    fixture.ROOT.mkdir()
    def sector_case(path, model, length, momentum, states, raw, projector):
        for layout in ('replicated', 'distributed'):
            sub = root / layout
            sub.mkdir()
            seed = sub / 'seed'
            symmetry_te.build_restart_seed(path, seed, layout)
            exercise(seed, sub, layout)
    fixture.prepare('SpinlessFermion', 6, sector_test=sector_case, sector_momenta=(0,))
    print('ok')


if __name__ == '__main__':
    main()
