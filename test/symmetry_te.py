"""Static sector quench versus independently projected dense propagation."""
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess

import numpy as np
import symmetry_general_terms as fixture

fixture.ROOT = Path('symmetry_te')
fixture.ROOT.mkdir(exist_ok=True)


def vector(path, label, expected_step, expected_time):
    parts = []
    files = sorted((path / 'output').glob(label + '_rank_*.dat'),
                   key=lambda p: int(p.stem.split('rank_')[1]))
    assert files
    for p in files:
        data = p.read_bytes()
        h = struct.unpack('<28Q', data[:224])
        assert h[22] == 4 and h[24] == expected_step
        assert struct.unpack('<d', data[200:208])[0] == expected_time
        assert h[14] == len(parts)
        parts.extend(np.frombuffer(data[224:], dtype='<c16'))
    return np.array(parts)


def check(path, model, length, momentum, states, raw, projector):
    columns, covered = [], set()
    for i in range(len(states)):
        col = projector[:, i]
        if i in covered or np.linalg.norm(col) < 1e-12:
            continue
        covered.update(np.flatnonzero(abs(col) > 1e-12))
        columns.append(col / np.linalg.norm(col))
    basis = np.column_stack(columns)
    hi = basis.conj().T @ raw @ basis
    dim = len(hi)
    calc, mod, names = [(path / f).read_text() for f in ['calc.def', 'mod.def', 'sym.def']]
    inter = (path / 'InterAll.def').read_text()
    # Change a translation-invariant density interaction, not just an energy offset.
    rows = [line.split() for line in inter.splitlines()[5:] if line.strip()]
    q = np.zeros(len(states))
    for i in range(length):
        j = (i+1) % length
        rows.append([i, 0, i, 0, j, 0, j, 0, .27, 0])
        for r, state in enumerate(states):
            scale = 2 if model in ('Hubbard', 'tJ') else 1
            ni, nj = (state >> (scale*i)) & 1, (state >> (scale*j)) & 1
            if model == 'Spin':
                ni, nj = 1-ni, 1-nj
            q[r] += .27*ni*nj
    hf = hi + basis.conj().T @ (q[:, None]*basis)
    eig, rot = np.linalg.eigh(hf)
    times = [0., .03, .09, .09, .13]
    n = length if model == 'Spin' else length//2 if model == 'SpinlessFermion' else 3
    sz = .5 if model in ('Hubbard', 'tJ') else 0
    d_raw = np.array([sum(((s >> (2*i)) & 3) == 3 for i in range(length))
                      if model in ('Hubbard', 'tJ') else 0 for s in states])
    d1 = basis.conj().T @ (d_raw[:, None]*basis)
    d2 = basis.conj().T @ ((d_raw**2)[:, None]*basis)
    for layout in ['distributed', 'replicated']:
        (path / 'InterAll.def').write_text(inter)
        (path / 'calc.def').write_text(calc + 'OutputEigenVec 1\n')
        (path / 'mod.def').write_text(mod)
        (path / 'sym.def').write_text(names)
        if (path / 'output').exists():
            shutil.rmtree(path / 'output')
        fixture.run(path, 'seed_k{}_{}'.format(momentum, layout), layout=layout)
        files = sorted((path / 'output').glob('zvo_eigenvec_0_rank_*.dat'),
                       key=lambda p: int(p.stem.split('rank_')[1]))
        initial = []
        for p in files:
            initial.extend(np.frombuffer(p.read_bytes()[224:], dtype='<c16'))
            shutil.copyfile(p, p.with_name(p.name.replace('zvo_', 'seed_')))
        initial = np.array(initial)
        np.testing.assert_allclose(hi @ initial, np.linalg.eigvalsh(hi)[0]*initial, atol=3e-7)
        fixture.definition(path, 'InterAll.def', rows)
        te_calc = calc.replace('CalcType 3', 'CalcType 4') + 'InputEigenVec 1\nOutputEigenVec 1\n'
        (path / 'calc.def').write_text(te_calc)
        (path / 'mod.def').write_text(mod.replace('Lanczos_max 400', 'Lanczos_max 5') + 'ExpandCoef 8\nOutputInterval 1\n')
        (path / 'sym.def').write_text(names + 'SpectrumVec seed_eigenvec_0\nTEOneBody times.def\n')
        (path / 'times.def').write_text('====\nNTimeSteps 5\n====\n====\n====\n' + ''.join('{} 0\n'.format(t) for t in times))
        fixture.run(path, 'evolve_k{}_{}'.format(momentum, layout), layout=layout)
        expected_ss, expected_norm, expected_flct = [], [], []
        state = initial.copy()
        for step, t in enumerate(times):
            dt = t-times[step-1] if step else 0
            x = -1j*dt*eig
            poly = sum(x**k / math.factorial(k) for k in range(9))
            state = rot @ (poly*(rot.conj().T @ state))
            norm = np.linalg.norm(state)
            state /= norm
            hv = hf @ state
            doublon, doublon2 = np.vdot(state, d1 @ state).real, np.vdot(state, d2 @ state).real
            expected_ss.append([t, np.vdot(state, hv).real, np.vdot(hv, hv).real, doublon, n, step])
            expected_norm.append([t, norm, step])
            expected_flct.append([t, n, n*n, doublon, doublon2, sz, sz*sz, step])
            np.testing.assert_allclose(vector(path, 'zvo_eigenvec_{}'.format(step), step, t), state, atol=2e-11, rtol=0)
        np.testing.assert_allclose(vector(path, 'zvo_eigenvec_final', 4, times[-1]), state, atol=2e-11, rtol=0)
        for family, expected in [('SS', expected_ss), ('Norm', expected_norm), ('Flct', expected_flct)]:
            np.testing.assert_allclose(np.loadtxt(path / 'output' / ('zvo_'+family+'.dat')), expected, atol=3e-11, rtol=0)
        exact = rot @ (np.exp(-1j*times[-1]*eig)*(rot.conj().T @ initial))
        np.testing.assert_allclose(state, exact, atol=2e-7, rtol=0)
        record = dict(line.split('=', 1) for line in (path / 'output/zvo_symmetry_sector.dat').read_text().splitlines())
        assert record['calc_type'] == 'TimeEvolution' and record['te_hamiltonian'] == 'static'
        assert record['source_method'] == '3' and record['expand_coef'] == '8'
        for step, t in enumerate(times):
            assert float(record['te_time_'+str(step)]) == t
        shutil.copytree(path / 'output', path / 'k{}_{}'.format(momentum, layout))
        if model == 'Spin' and momentum == 0:
            # A TE checkpoint imports a state into a NEW clock, with the default
            # distributed layout as well as the explicit rollback layout.
            (path / 'sym.def').write_text(names + 'SpectrumVec zvo_eigenvec_final\nTEOneBody times.def\n')
            env = dict(os.environ)
            if layout == 'distributed':
                env.pop('HPHI_SYMMETRY_BASIS_LAYOUT', None)
            else:
                env['HPHI_SYMMETRY_BASIS_LAYOUT'] = layout
            with (path / ('import_te_'+layout+'.log')).open('w') as fp:
                result = subprocess.run(fixture.MPI + [fixture.HPHI, '-e', 'sym.def'], cwd=path,
                                        env=env, stdout=fp, stderr=subprocess.STDOUT, timeout=120)
            assert result.returncode == 0, (path / ('import_te_'+layout+'.log')).read_text()
            np.testing.assert_allclose(vector(path, 'zvo_eigenvec_0', 0, 0), state, atol=2e-11, rtol=0)
            second = state.copy()
            for dt in np.diff(times):
                poly = sum((-1j*dt*eig)**k / math.factorial(k) for k in range(9))
                second = rot @ (poly*(rot.conj().T @ second))
                second /= np.linalg.norm(second)
            np.testing.assert_allclose(vector(path, 'zvo_eigenvec_final', 4, times[-1]), second, atol=3e-11, rtol=0)
            record = dict(line.split('=', 1) for line in (path / 'output/zvo_symmetry_sector.dat').read_text().splitlines())
            assert record['source_method'] == '4' and float(record['source_time']) == times[-1]
            (path / 'sym.def').write_text(names + 'SpectrumVec seed_eigenvec_0\nTEOneBody times.def\n')
            # These fail before propagation; do not depend on a particular MPI partition.
            (path / 'calc.def').write_text(te_calc.replace('InputEigenVec 1', 'InputEigenVec 2'))
            fixture.run(path, 'reject_input_'+layout, layout=layout, fail='requires sector checkpoint InputEigenVec=1')
            (path / 'calc.def').write_text(te_calc + 'ReStart 1\n')
            fixture.run(path, 'reject_restart_'+layout, layout=layout, fail='ReStart')
            (path / 'calc.def').write_text(te_calc)
            (path / 'times.def').write_text('====\nNTimeSteps 5\n====\n====\n====\n' + ''.join('{} 1\n0 0 0 0 .1 0\n'.format(t) for t in times))
            fixture.run(path, 'reject_dynamic_'+layout, layout=layout, fail='requires a static Hamiltonian')
    (path / 'InterAll.def').write_text(inter)
    (path / 'calc.def').write_text(calc)
    (path / 'mod.def').write_text(mod)
    (path / 'sym.def').write_text(names)
    print('{} k={}: static quench, every vector and observable PASS'.format(model, momentum))


for model, length in [('Spin', 6), ('SpinlessFermion', 6), ('Hubbard', 4), ('tJ', 4)]:
    fixture.prepare(model, length, sector_test=check)
