"""Independent right-endpoint Taylor and state-import tests for SpinGC TE."""
import os
from pathlib import Path
import shlex
import shutil
import struct
import sys
import tempfile

import numpy as np
import symmetry_spingc_common as c
import symmetry_spingc_observables as observable
import symmetry_spingc_thermal as thermal

TIMES = np.array([0., .01, .025, .025, .04])
AMPLITUDES = [0, .7, 0, -.4, .2]


def checkpoints(path, label):
    return sorted((path/'output').glob(label+'_rank_*.dat'),
                  key=lambda p: int(p.stem.rsplit('_', 1)[1]))


def schedule(path, family, rows, times=TIMES, amplitudes=AMPLITUDES):
    text = '====\nNTimeSteps 5\n====\n====\n====\n'
    for time, amplitude in zip(times, amplitudes):
        active = rows if amplitude else []
        text += '{} {}\n'.format(time, len(active))
        for row in active:
            text += ' '.join(map(str, list(row[:-2])+[
                amplitude*row[-2], amplitude*row[-1]]))+'\n'
    (path/(family+'.def')).write_text(text)
    with (path/'sym.def').open('a') as handle:
        handle.write('{} {}.def\n'.format(family, family))


def prepare(root, name, fixture, hphi, launcher, env):
    n, permutations, characters, families, basis, reps, h = fixture
    path = root/name
    c.write_case(path, n, permutations, characters, families, 3,
                 {'CalcMod': {'OutputEigenVec': 1},
                  'ModPara': {'LanczosEps': 18}})
    c.run_case(path, hphi, 'cg_seed', launcher, env)
    return path, c.join_rank_vectors(checkpoints(path, 'zvo_eigenvec_0'), len(reps))


def evolve(root, name, fixture, source, source_label, drive, interval,
           output_interval, hphi, launcher, env, times=TIMES):
    n, permutations, characters, original, basis, reps, h = fixture
    families = {key: list(rows) for key, rows in original.items()}
    ops = c.spin_operators(n)
    sy = basis.conj().T@sum(ops['y'])@basis
    pair = basis.conj().T@sum(ops['plus'][i]@ops['plus'][(i+1)%n]
                             +ops['minus'][i]@ops['minus'][(i+1)%n]
                             for i in range(n))@basis
    one = [(i,a,i,b,-(.21*c.Sy[a,b]).real,-(.21*c.Sy[a,b]).imag)
           for i in range(n) for a,b in ((0,1),(1,0))]
    two = [row for i in range(n) for row in
           c.interall_pair((i,1,i,0,(i+1)%n,1,(i+1)%n,0), .19)]
    if drive == 'quench':
        families['Trans'] += one
        h = h+.21*sy
    hs = [h+amplitude*(.21*sy if drive == 'one' else .19*pair)
          if drive in ('one','two') else h for amplitude in AMPLITUDES]
    path = root/name
    c.write_case(path, n, permutations, characters, families, 4,
                 {'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1},
                  'ModPara': {'Lanczos_max': 5, 'ExpandCoef': 8,
                              'ExpecInterval': interval, 'OutputInterval': output_interval}})
    observable.add_correlation_requests(path)
    (path/'output').mkdir(exist_ok=True)
    inputs = checkpoints(source, source_label)
    for rank, file in enumerate(inputs):
        shutil.copy2(file, path/'output'/('seed_rank_{}.dat'.format(rank)))
    with (path/'sym.def').open('a') as handle:
        handle.write('SpectrumVec seed\n')
    schedule(path, 'TETwoBody' if drive == 'two' else 'TEOneBody',
             two if drive == 'two' else one if drive == 'one' else [], times)
    text = c.run_case(path, hphi, 'evolve', launcher, env)
    assert 'Symmetry basis layout: '+env['HPHI_SYMMETRY_BASIS_LAYOUT'] in text
    if 'OMP_NUM_THREADS' in env:
        assert 'OpenMP threads : '+env['OMP_NUM_THREADS'] in text
    initial = c.join_rank_vectors(inputs, len(reps))
    vector, exact = initial.copy(), initial.copy()
    sz = np.array([int(rep).bit_count()-n/2 for rep in reps])
    operators = thermal.correlation_operators(n, basis)
    ss, norms, flct, vectors = [], [], [], []
    metadata = dict(line.split('=',1) for line in
                    (path/'output/zvo_symmetry_sector.dat').read_text().splitlines())
    base_digest = int(metadata['hamiltonian_digest'].split(':')[1],16)
    source_header, _ = c.read_checkpoint(inputs[0])
    assert int(metadata['source_method']) == source_header[22]
    assert int(metadata['source_step']) == source_header[24]
    assert float(metadata['source_time']) == struct.unpack('<d',struct.pack('<Q',source_header[25]))[0]
    assert int(metadata['source_hamiltonian_digest'],16) == source_header[21]
    assert metadata['te_hamiltonian'] == ('time_dependent' if drive in ('one','two') else 'static')
    correlation_values = {kind: [] for kind in operators}
    for step, (time, hstep) in enumerate(zip(times, hs)):
        dt = time-(times[step-1] if step else times[0])
        assert dt*np.linalg.norm(hstep,2) <= .25
        vector, norm = c.taylor_step(hstep, vector, -1j*dt, 8)
        eig, rot = np.linalg.eigh(hstep)
        exact = rot@(np.exp(-1j*dt*eig)*(rot.conj().T@exact))
        assert abs(norm-1) <= 1e-8 and np.linalg.norm(vector-exact) <= 1e-8
        vectors.append(vector.copy())
        hv = hstep@vector
        ss.append([time,np.vdot(vector,hv).real,np.vdot(hv,hv).real,0,n,step])
        norms.append([time,norm,step])
        flct.append([time,n,n*n,0,0,np.vdot(vector,sz*vector).real,
                     np.vdot(vector,sz*sz*vector).real,step])
        files = checkpoints(path, 'zvo_eigenvec_'+str(step))
        if step % output_interval == 0:
            assert len(files) == c.mpi_size(launcher)
            actual = c.join_rank_vectors(files,len(reps))
            assert np.linalg.norm(actual-vector) <= 3e-11
            digest = int(metadata['te_hamiltonian_'+str(step)],16) if drive in ('one','two') else base_digest
            for file in files:
                header, _ = c.read_checkpoint(file)
                assert header[21] == digest and header[22:25] == (4,0,step)
                assert struct.unpack('<d',struct.pack('<Q',header[25]))[0] == time
        else:
            assert not files
        for kind, rows in observable.CORRELATION_REQUESTS.items():
            file = path/'output'/('zvo_{}_step{}.dat'.format(observable.CORRELATION_FILES[kind],step))
            if step % interval:
                assert not file.exists(), file
            else:
                data = thermal.load_finite(file)
                expected = thermal.check_correlation_block(data,rows,
                                                operators[kind],vector,0,str(file))
                correlation_values[kind].extend(expected)
    for family, expected in [('SS',ss),('Norm',norms),('Flct',flct)]:
        actual = thermal.load_finite(path/'output'/('zvo_'+family+'.dat'))
        np.testing.assert_array_equal(actual[:,-1],np.asarray(expected)[:,-1])
        np.testing.assert_allclose(actual[:,:-1],np.asarray(expected)[:,:-1],atol=4e-11,rtol=0)
    for kind, values in correlation_values.items():
        assert np.max(np.abs(values)) > 1e-6, kind
        assert np.max(np.abs(np.imag(values))) > 1e-6, kind
    final_files = checkpoints(path,'zvo_eigenvec_final')
    assert len(final_files) == c.mpi_size(launcher)
    parts = [c.read_checkpoint(file) for file in final_files]
    if n == 6 and c.mpi_size(launcher) == 16:
        assert len(reps) == 9 and sum(len(v) == 0 for _, v in parts) == 7
    assert all(header[26:28] == c.checkpoint_digest(parts) for header, _ in parts)
    for header, _ in parts:
        assert header[22:25] == (4,0,4)
        assert struct.unpack('<d',struct.pack('<Q',header[25]))[0] == times[-1]
    final = c.join_rank_vectors(checkpoints(path,'zvo_eigenvec_final'),len(reps))
    assert np.linalg.norm(final-vector) <= 3e-11
    if drive in ('one','two'):
        assert metadata['te_hamiltonian'] == 'time_dependent'
        assert metadata['te_plan_update'] == 'rebuild_with_fixed_basis'
        assert metadata['te_hamiltonian_0'] == metadata['te_hamiltonian_2']
        assert len({metadata['te_hamiltonian_'+str(i)] for i in range(5)}) == 4
    if drive == 'quench':
        assert 'hamiltonian_changed=yes' in text
    print(name+' passed',flush=True)
    return path, vectors



def preflight_negatives(root, source, hphi, env):
    for name, diagnostic in [('late_noninvariant','preflight failed at step 4'),
                             ('nan_time','finite and nondecreasing'),
                             ('inf_time','finite and nondecreasing'),
                             ('decreasing_time','finite and nondecreasing'),
                             ('nonhermitian','NonHermite'),
                             ('bad_site','preflight failed'),
                             ('mixed_drives','invalid or mixed symmetry TE'),
                             ('laser','SpinGC TransSym does not support Laser')]:
        path = root/('reject_'+name)
        path.mkdir()
        for file in source.glob('*.def'):
            shutil.copy2(file,path/file.name)
        (path/'output').mkdir()
        for file in (source/'output').glob('seed_rank_*.dat'):
            shutil.copy2(file,path/'output'/file.name)
        names = (path/'sym.def').read_text()
        (path/'sym.def').write_text(''.join(line+'\n' for line in names.splitlines()
                                            if not line.startswith('TEOneBody ')))
        if name in ('late_noninvariant','nonhermitian','bad_site'):
            rows = [(0,0,0,0,.3,0)]
            if name == 'nonhermitian':
                rows = [(0,1,0,0,.3,0)]
            if name == 'bad_site':
                rows = [(8,0,8,0,.3,0)]
            schedule(path,'TEOneBody',rows,amplitudes=[0,0,0,0,1])
        elif name.endswith('time'):
            times = TIMES.copy()
            times[4] = {'nan_time':float('nan'),'inf_time':float('inf'),
                        'decreasing_time':.02}[name]
            schedule(path,'TEOneBody',[],times)
        elif name == 'mixed_drives':
            schedule(path,'TEOneBody',[(i,0,i,0,.1,0) for i in range(8)])
            schedule(path,'TETwoBody',[(i,0,i,0,(i+1)%8,0,(i+1)%8,0,.1,0)
                                       for i in range(8)])
        else:
            (path/'laser.def').write_text('====\nNLaser 9\n====\n====\n====\n'+
                ''.join('p{} {}\n'.format(i,v) for i,v in
                        enumerate([2,.4,1,1,0,8,1,1,0])))
            with (path/'sym.def').open('a') as handle:
                handle.write('Laser laser.def\n')
        c.run_case(path,hphi,'reject',[],env,diagnostic)
        assert not (path/'output/zvo_SS.dat').exists()
        assert not list((path/'output').glob('zvo_eigenvec*'))
        print('rejected '+name,flush=True)


def main():
    hphi = Path(sys.argv[1]).resolve()
    launcher = shlex.split(os.environ.get('MPIRUN',''))
    ranks = c.mpi_size(launcher)
    assert ranks in (1,4,16)
    root = Path(tempfile.mkdtemp(prefix='symmetry_spingc_te_',dir='.')).resolve()
    print('artifacts: '+str(root),flush=True)
    for layout in (('distributed','replicated') if ranks <= 4 else ('distributed',)):
        for label, momentum in [('B',0),('B',1)]+([('D',1)] if ranks == 16 else []):
            fixture = thermal.fixture(label,momentum)
            env = dict(os.environ,HPHI_SYMMETRY_BASIS_LAYOUT=layout)
            prefix = '{}_k{}_{}'.format(label,momentum,layout)
            seed, _ = prepare(root,prefix+'_seed',fixture,hphi,launcher,env)
            for drive in ('static','quench','one','two'):
                for interval, output_interval in ((1,1),(2,2)):
                    result, _ = evolve(root,prefix+'_{}_i{}'.format(drive,interval),fixture,
                           seed,'zvo_eigenvec_0',drive,interval,output_interval,
                           hphi,launcher,env)
                    if ranks == 1 and layout == 'distributed' and momentum == 0 and drive == 'one' and interval == 1:
                        preflight_negatives(root,result,hphi,env)


if __name__ == '__main__':
    main()
