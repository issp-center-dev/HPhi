"""Independent finite-sample TPQ/cTPQ validation for SpinGC sectors."""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile

import numpy as np

import symmetry_spingc_common as c
import symmetry_spingc_observables as observable


ATOL_DATA = 2e-11
ATOL_CORRELATION = 1e-8
NUM_AVE = 2
NSTEPS = 5


def translation(nsite, momentum):
    permutations = [[(site+shift) % nsite for site in range(nsite)]
                    for shift in range(nsite)]
    characters = np.exp(-2j*np.pi*momentum*np.arange(nsite)/nsite)
    return permutations, characters


def fixture(label, momentum):
    nsite = 6 if label == 'D' else 8
    permutations, characters = translation(nsite, momentum)
    families = c.mixed_families(nsite)
    basis, representatives, _ = c.sector_basis(
        nsite, permutations, characters)
    hamiltonian = basis.conj().T@c.mixed_hamiltonian(nsite)@basis
    np.testing.assert_allclose(hamiltonian, hamiltonian.conj().T,
                               atol=1e-12, rtol=0)
    return (nsite, permutations, characters, families, basis,
            representatives, hamiltonian)


def invoke(path, command, label, env, expected_error=None):
    result = subprocess.run(command, cwd=path, env=env, capture_output=True,
                            text=True, timeout=120)
    text = result.stdout+result.stderr
    (path/(label+'.log')).write_text(text)
    if expected_error is None:
        assert result.returncode == 0, text
    else:
        assert result.returncode != 0 and expected_error in text, text
    return text


def initial_vectors(path, helper, launcher, env, dimension, dtype):
    invoke(path, launcher+[str(helper), '--dump-initial', str(dimension),
                          str(dtype), '7', str(NUM_AVE)], 'initial', env)
    ranks = c.mpi_size(launcher)
    vectors, norms = [], []
    for sample in range(NUM_AVE):
        vector = []
        first_norm = None
        for rank in range(ranks):
            file = path/'initial_{}_rank{}.dat'.format(sample, rank)
            lines = file.read_text().splitlines()
            rank_norm = float(lines[0])
            if first_norm is None:
                first_norm = rank_norm
            else:
                np.testing.assert_allclose(rank_norm, first_norm,
                                           atol=0, rtol=0)
            vector.extend(complex(*map(float, line.split()))
                          for line in lines[1:])
        vector = np.asarray(vector)
        assert len(vector) == dimension and np.isfinite(vector).all()
        np.testing.assert_allclose(np.vdot(vector, vector), 1,
                                   atol=1e-12, rtol=0)
        vectors.append(vector)
        norms.append(first_norm)
    return vectors, norms


def correlation_operators(nsite, basis):
    return {
        kind: [basis.conj().T@observable.correlation_operator(nsite, kind, row)@basis
               for row in rows]
        for kind, rows in observable.CORRELATION_REQUESTS.items()
    }


def tpq_step(h, v, nsite):
    w = 4*v-h@v/nsite
    norm = np.linalg.norm(w)
    assert np.isfinite(norm) and norm > 0
    return w/norm, norm


def expected_history(kind, hamiltonian, initial, first_norm, nsite,
                     representatives, betas=None):
    vector = initial.copy()
    norm = first_norm
    ss, norms, flct, vectors = [], [], [], []
    magnetization = np.array([int(rep).bit_count()-nsite/2
                              for rep in representatives])
    for step in range(NSTEPS):
        vectors.append(vector.copy())
        hv = hamiltonian@vector
        energy = np.vdot(vector, hv).real
        energy2 = np.vdot(hv, hv).real
        sz = np.vdot(vector, magnetization*vector).real
        sz2 = np.vdot(vector, magnetization*magnetization*vector).real
        beta = (2*step/(4*nsite-energy) if kind == 'tpq'
                else betas[step])
        ss.append([beta, energy, energy2, 0, nsite, step])
        norms.append([beta, norm, first_norm, step])
        if kind != 'tpq' or step:
            flct.append([beta, nsite, nsite**2, 0, 0, sz, sz2, step])
        if step+1 < NSTEPS:
            if kind == 'tpq':
                vector, norm = tpq_step(hamiltonian, vector, nsite)
            else:
                vector, norm = c.taylor_step(
                    hamiltonian, vector, -(betas[step+1]-betas[step])/2, 12)
    return dict(SS=np.asarray(ss), Norm=np.asarray(norms),
                Flct=np.asarray(flct), vectors=vectors)


def load_finite(path):
    data = np.loadtxt(path, ndmin=2)
    assert np.isfinite(data).all(), path
    return data


def check_data_files(path, aggregate, expected):
    for family in ('SS', 'Norm', 'Flct'):
        if aggregate:
            file = path/'output'/('zvo_{}_tpq.dat'.format(family))
            data = load_finite(file)
            sequence = [(sample, int(row[-1]))
                        for sample in range(NUM_AVE)
                        for row in expected[sample][family]]
            np.testing.assert_array_equal(data[:, :2], sequence)
            joined = np.vstack([value[family][:, :-1]
                                for value in expected])
            np.testing.assert_allclose(data[:, 2:], joined,
                                       atol=ATOL_DATA, rtol=ATOL_DATA,
                                       err_msg=family)
            assert not list((path/'output').glob('zvo_{}_rand*.dat'.format(family)))
        else:
            assert not (path/'output'/('zvo_{}_tpq.dat'.format(family))).exists()
            for sample in range(NUM_AVE):
                data = load_finite(
                    path/'output'/('zvo_{}_rand{}.dat'.format(family, sample)))
                np.testing.assert_array_equal(
                    data[:, -1], expected[sample][family][:, -1])
                np.testing.assert_allclose(data[:, :-1],
                                           expected[sample][family][:, :-1],
                                           atol=ATOL_DATA, rtol=ATOL_DATA,
                                           err_msg='{} sample{}'.format(family, sample))


def check_correlation_block(data, rows, operators, vector, prefix, label):
    expected_shape = (len(rows), prefix+len(rows[0])+2)
    assert data.shape == expected_shape, (label, data.shape, expected_shape)
    np.testing.assert_array_equal(data[:, prefix:-2], np.asarray(rows))
    values = data[:, -2]+1j*data[:, -1]
    expected = np.asarray([np.vdot(vector, operator@vector)
                           for operator in operators])
    assert np.isfinite(values).all() and np.isfinite(expected).all()
    np.testing.assert_allclose(values, expected, atol=ATOL_CORRELATION,
                               rtol=0, err_msg=label)
    return expected


def check_correlations(path, aggregate, emitted_steps, expected, operators):
    all_values = {kind: [] for kind in observable.CORRELATION_REQUESTS}
    for kind, rows in observable.CORRELATION_REQUESTS.items():
        stem = observable.CORRELATION_FILES[kind]
        if aggregate:
            file = path/'output'/('zvo_{}_tpq.dat'.format(stem))
            data = load_finite(file)
            prefixes = [(sample, step) for sample in range(NUM_AVE)
                        for step in emitted_steps for _ in rows]
            np.testing.assert_array_equal(data[:, :2], prefixes)
            assert data.shape[0] == len(prefixes)
            offset = 0
            for sample in range(NUM_AVE):
                for step in emitted_steps:
                    block = data[offset:offset+len(rows)]
                    all_values[kind].extend(check_correlation_block(
                        block, rows, operators[kind],
                        expected[sample]['vectors'][step], 2,
                        '{} sample{} step{}'.format(kind, sample, step)))
                    offset += len(rows)
            assert not list((path/'output').glob(
                'zvo_{}_set*step*.dat'.format(stem)))
        else:
            assert not (path/'output'/('zvo_{}_tpq.dat'.format(stem))).exists()
            for sample in range(NUM_AVE):
                for step in range(NSTEPS):
                    file = path/'output'/('zvo_{}_set{}step{}.dat'.format(
                        stem, sample, step))
                    if step not in emitted_steps:
                        assert not file.exists(), file
                        continue
                    all_values[kind].extend(check_correlation_block(
                        load_finite(file), rows, operators[kind],
                        expected[sample]['vectors'][step], 0,
                        '{} sample{} step{}'.format(kind, sample, step)))
    for kind, values in all_values.items():
        values = np.asarray(values)
        assert np.max(np.abs(values)) > 1e-6, (kind, values)
        assert np.max(np.abs(values.imag)) > 1e-6, (kind, values)


def metadata(path, kind, dimension, dtype, ranks, explicit):
    values = dict(line.split('=', 1) for line in
                  (path/'output/zvo_symmetry_sector.dat').read_text().splitlines())
    assert values['ensemble'] == 'single_symmetry_sector'
    assert values['calc_type'] == ('TPQ' if kind == 'tpq' else 'cTPQ')
    assert int(values['sector_dim']) == dimension
    assert int(values['num_ave']) == NUM_AVE
    assert int(values['initial_vec_type']) == dtype
    assert int(values['mpi_ranks']) == ranks
    if kind == 'tpq':
        assert float(values['large_value']) == 4
    else:
        assert int(values['canonical_tpq_steps']) == NSTEPS
        assert values['beta_schedule'] == ('explicit' if explicit else 'uniform')
        if explicit:
            for step, beta in enumerate([0, .02, .05, .05, .08]):
                row = list(map(float, values['invtemp_row_{}'.format(step)].split()))
                np.testing.assert_allclose(row, [beta, 12], atol=0, rtol=0)
        else:
            assert float(values['large_value']) == 50
            assert float(values['beta_step']) == .02
            assert int(values['expand_coef']) == 12


def run_thermal(root, label, momentum, kind, explicit, aggregate, interval,
                layout, hphi, helper, launcher, base_env):
    (nsite, permutations, characters, families, basis, representatives,
     hamiltonian) = fixture(label, momentum)
    dimension = len(representatives)
    if label == 'D':
        assert dimension == 9 and c.mpi_size(launcher) == 16
    dtype = 1 if momentum == 0 else 0
    schedule = 'explicit' if explicit else 'uniform'
    path = root/'{}_{}_k{}_{}_{}_{}_i{}'.format(
        label, layout, momentum, kind, schedule,
        'aggregate' if aggregate else 'legacy', interval)
    calc_type = 1 if kind == 'tpq' else 5
    c.write_case(path, nsite, permutations, characters, families, calc_type,
                 {'CalcMod': {'InitialVecType': dtype,
                              'OutputGreenFormat': int(aggregate)},
                  'ModPara': {'Lanczos_max': NSTEPS, 'initial_iv': 7,
                              'NumAve': NUM_AVE, 'ExpecInterval': interval,
                              'LargeValue': 4 if kind == 'tpq' else 50,
                              **({'ExpandCoef': 12} if kind == 'ctpq' else {})}})
    observable.add_correlation_requests(path)
    if explicit:
        betas = np.array([0, .02, .05, .05, .08])
        flags = [0, 1, 0, 1, 0]
        (path/'beta.def').write_text(''.join(
            '{} 12 {} 0\n'.format(beta, flag)
            for beta, flag in zip(betas, flags)))
        (path/'sym.def').write_text((path/'sym.def').read_text()+'InvTemp beta.def\n')
        emitted_steps = [0, 1, 3]
    elif kind == 'ctpq':
        betas = np.arange(NSTEPS)*.02
        emitted_steps = [0]+[step for step in range(1, NSTEPS)
                             if step % interval == 0]
    else:
        betas = None
        emitted_steps = [0]+[step for step in range(2, NSTEPS)
                             if step % interval == 0]
    env = dict(base_env, HPHI_SYMMETRY_BASIS_LAYOUT=layout)
    initials, first_norms = initial_vectors(
        path, helper, launcher, env, dimension, dtype)
    expected = [expected_history(kind, hamiltonian, vector, norm, nsite,
                                 representatives, betas)
                for vector, norm in zip(initials, first_norms)]
    operators = correlation_operators(nsite, basis)
    text = invoke(path, launcher+[str(hphi), '-e', 'sym.def'], 'run', env)
    assert 'Symmetry basis layout: {}'.format(layout) in text
    if 'OMP_NUM_THREADS' in env:
        assert 'OpenMP threads : {}'.format(env['OMP_NUM_THREADS']) in text
    metadata(path, kind, dimension, dtype, c.mpi_size(launcher), explicit)
    check_data_files(path, aggregate, expected)
    check_correlations(path, aggregate, emitted_steps, expected, operators)
    print('{} k={} {} {} {} interval{} dim={} passed'.format(
        label, momentum, kind, schedule, layout, interval, dimension), flush=True)
    return path


def reject_invtemp_eigen(root, source, hphi, base_env):
    for target in (0, 2, 4):
        path = root/'reject_eigen_{}'.format(target)
        shutil.copytree(source, path)
        rows = (path/'beta.def').read_text().splitlines()
        fields = rows[target].split()
        fields[-1] = '1'
        rows[target] = ' '.join(fields)
        (path/'beta.def').write_text('\n'.join(rows)+'\n')
        output = path/'output'
        shutil.rmtree(output)
        invoke(path, [str(hphi), '-e', 'sym.def'],
               'reject_eigen_{}'.format(target), base_env,
               'does not support InvTemp eigenvector output')
        assert not list(output.glob('*eigen*'))


def main():
    hphi, helper = [Path(arg).resolve() for arg in sys.argv[1:]]
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = c.mpi_size(launcher)
    assert ranks in (1, 4, 16), 'thermal suite requires np1, np4, or np16'
    root = Path(tempfile.mkdtemp(prefix='symmetry_spingc_thermal_', dir='.'))
    print('artifacts: {}'.format(root.resolve()), flush=True)
    base_env = dict(os.environ)
    layouts = ['distributed']+(['replicated'] if ranks in (1, 4) else [])
    explicit_rejection_source = None
    for layout in layouts:
        for momentum, aggregate, interval in ((0, False, 1), (1, True, 2)):
            run_thermal(root, 'B', momentum, 'tpq', False, aggregate,
                        interval, layout, hphi, helper, launcher, base_env)
            run_thermal(root, 'B', momentum, 'ctpq', False, aggregate,
                        interval, layout, hphi, helper, launcher, base_env)
            explicit_path = run_thermal(
                root, 'B', momentum, 'ctpq', True, aggregate,
                interval, layout, hphi, helper, launcher, base_env)
            if ranks == 1 and layout == 'distributed' and momentum == 0:
                explicit_rejection_source = explicit_path
    if ranks == 1:
        run_thermal(root, 'B', 1, 'tpq', False, False, 1,
                    'distributed', hphi, helper, launcher, base_env)
        run_thermal(root, 'B', 1, 'ctpq', False, False, 1,
                    'distributed', hphi, helper, launcher, base_env)
        assert explicit_rejection_source is not None
        reject_invtemp_eigen(root, explicit_rejection_source, hphi, base_env)
    if ranks == 16:
        run_thermal(root, 'D', 1, 'tpq', False, True, 2,
                    'distributed', hphi, helper, launcher, base_env)
        run_thermal(root, 'D', 1, 'ctpq', False, True, 2,
                    'distributed', hphi, helper, launcher, base_env)


if __name__ == '__main__':
    main()
