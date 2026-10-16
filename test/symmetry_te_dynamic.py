"""All time slices against independent dense operators, including changing sparsity."""
import math
from pathlib import Path
import shutil
import struct

import numpy as np
import symmetry_correlation as correlation
import symmetry_general_terms as fixture

fixture.ROOT = Path('symmetry_te_dynamic')
fixture.ROOT.mkdir(exist_ok=True)


def checkpoint(path, label):
    parts, headers = [], []
    for p in sorted((path / 'output').glob(label+'_rank_*.dat'), key=lambda p: int(p.stem.split('rank_')[1])):
        data = p.read_bytes()
        headers.append(struct.unpack('<28Q', data[:224]))
        parts.extend(np.frombuffer(data[224:], dtype='<c16'))
    assert headers
    return np.array(parts), headers


def check(path, model, length, momentum, states, raw, projector):
    cols, seen = [], set()
    for i in range(len(states)):
        col = projector[:, i]
        if i in seen or np.linalg.norm(col) < 1e-12:
            continue
        seen.update(np.flatnonzero(abs(col) > 1e-12))
        cols.append(col / np.linalg.norm(col))
    basis = np.column_stack(cols)
    dim = len(cols)
    scale = 2 if model in ('Hubbard', 'tJ') else 1
    width = length*scale
    annihilate, cache = {}, {}
    if model != 'Spin':
        for i in range(width):
            annihilate[i] = fixture.tensor_operator(width, i, np.array([[0,1],[0,0]]), True)

    def op(indices):
        key = tuple(indices)
        if key not in cache:
            i,s,j,t = key
            if model == 'Spin':
                assert i == j
                mat = np.zeros((2,2)); mat[s,t] = 1
                cache[key] = fixture.tensor_operator(width, i, mat)
            else:
                cache[key] = annihilate[scale*i+s].T @ annihilate[scale*j+t]
        return cache[key]

    def projected(matrix):
        return basis.conj().T @ matrix[np.ix_(states,states)] @ basis

    hi = basis.conj().T @ raw @ basis
    calc, mod, names = [(path/f).read_text() for f in ('calc.def','mod.def','sym.def')]
    requests = correlation.step_requests(model, length)
    operators = correlation.sector_operators(model, length, states, basis, requests)
    correlation.add_requests(path, 'sym.def', requests)
    correlation_names = (path/'sym.def').read_text()
    (path/'sym.def').write_text(names)
    amplitudes = [0, .21, -.17, 0, .31]
    families = ['TEOneBody', 'TETwoBody'] + ([] if model == 'Spin' else ['Laser'])
    for family in families:
        times = np.arange(5)*.035 if family == 'Laser' else np.array([0., .03, .09, .12, .17])
        rows = []
        delta = np.zeros_like(raw)
        for i in range(length):
            j, k = (i+1) % length, (i+2) % length
            if family == 'TEOneBody':
                if model == 'Spin':
                    rows.append([i,0,i,0,.23,0])
                    delta -= .23*op([i,0,i,0])[np.ix_(states,states)]
                else:
                    # Complex hopping and a diagonal potential, with both spin flavours.
                    for spin in range(scale):
                        for a,b,v in [(i,j,.31+.12j),(j,i,.31-.12j),(i,i,.13+0j)]:
                            row=[a,spin,b,spin,v.real,v.imag];rows.append(row)
                            delta -= v*op(row[:4])[np.ix_(states,states)]
            elif family == 'TETwoBody':
                ix=[i,0,i,0,j,0,j,0]
                rows.append(ix+[.37,0]);delta += .37*(op(ix[:4])@op(ix[4:]))[np.ix_(states,states)]
                # Same-orbital contraction n_i n_i, exercising TEChemi bookkeeping.
                ix=[i,0,i,0,i,0,i,0]
                rows.append(ix+[.11,0]);delta += .11*(op(ix[:4])@op(ix[4:]))[np.ix_(states,states)]
                ix=([i,0,i,1,j,1,j,0] if model == 'Spin' else [i,0,j,0,k,0,k,0])
                a,s,b,t,c,u,d,v=ix
                dagger=[d,v,c,u,b,t,a,s]
                for indices,value in [(ix,.19+.07j),(dagger,.19-.07j)]:
                    rows.append(indices+[value.real,value.imag])
                    delta += value*(op(indices[:4])@op(indices[4:]))[np.ix_(states,states)]
        hs = []
        for step,time in enumerate(times):
            h = hi + amplitudes[step]*(basis.conj().T @ delta @ basis)
            if family == 'Laser':
                h=hi.copy()
                for row in (path/'Trans.def').read_text().splitlines()[5:]:
                    f=row.split()
                    if not f:continue
                    ix=list(map(int,f[:4]));value=complex(*map(float,f[4:]))
                    distance=(ix[0]-ix[2]+length//2) % length-length//2
                    phase=np.exp(-1j*.4*time*distance)
                    h -= value*(phase-1)*projected(op(ix))
            np.testing.assert_allclose(h,h.conj().T,atol=1e-12)
            hs.append(h)
        schedule='====\nNTimeSteps 5\n====\n====\n====\n'
        for step,time in enumerate(times):
            # An actually empty term list tests plan graph removal, not only zero values.
            active=rows if amplitudes[step] else []
            schedule+='{} {}\n'.format(time,len(active))
            for row in active:
                schedule+=' '.join(map(str,row[:-2]+[amplitudes[step]*row[-2],amplitudes[step]*row[-1]]))+'\n'
        for layout in ['distributed','replicated']:
            (path/'calc.def').write_text(calc+'OutputEigenVec 1\n')
            (path/'mod.def').write_text(mod)
            (path/'sym.def').write_text(names)
            if (path/'output').exists():shutil.rmtree(path/'output')
            fixture.run(path,'seed_{}_k{}_{}'.format(family,momentum,layout),layout=layout)
            state,_=checkpoint(path,'zvo_eigenvec_0')
            for p in (path/'output').glob('zvo_eigenvec_0_rank_*.dat'):
                shutil.copyfile(p,p.with_name(p.name.replace('zvo_','seed_')))
            (path/'calc.def').write_text(calc.replace('CalcType 3','CalcType 4')+'InputEigenVec 1\nOutputEigenVec 1\n')
            (path/'mod.def').write_text(mod.replace('Lanczos_max 400','Lanczos_max 5')+'ExpandCoef 8\nOutputInterval 1\nTimeSlice .035\nTinit 0\n')
            (path/'sym.def').write_text(correlation_names+'SpectrumVec seed_eigenvec_0\n'+family+' drive.def\n')
            (path/'drive.def').write_text(schedule if family != 'Laser' else '====\nNLaser 9\n====\n====\n====\n'+''.join('p{} {}\n'.format(i,v) for i,v in enumerate([2,.4,1,1,0,length,1,1,0])))
            fixture.run(path,'evolve_{}_k{}_{}'.format(family,momentum,layout),layout=layout)
            number = length if model == 'Spin' else length//2 if model == 'SpinlessFermion' else 3
            sz = .5 if model in ('Hubbard','tJ') else 0
            d_raw = np.array([sum(((s >> (2*i)) & 3) == 3 for i in range(length))
                              if model in ('Hubbard','tJ') else 0 for s in states])
            d1 = basis.conj().T @ (d_raw[:,None]*basis)
            d2 = basis.conj().T @ ((d_raw**2)[:,None]*basis)
            ss,norm,flct,states_by_step=[],[],[],[]
            for step,(time,h) in enumerate(zip(times,hs)):
                dt=time-times[step-1] if step else 0
                eig,rot=np.linalg.eigh(h)
                poly=sum((-1j*dt*eig)**n/math.factorial(n) for n in range(9))
                state=rot @ (poly*(rot.conj().T @ state))
                before=np.linalg.norm(state);state/=before
                states_by_step.append(state.copy())
                actual,headers=checkpoint(path,'zvo_eigenvec_'+str(step))
                np.testing.assert_allclose(actual,state,atol=3e-11,rtol=0)
                assert all(hd[24]==step and struct.unpack('<d',struct.pack('<Q',hd[25]))[0]==time for hd in headers)
                hv=h@state
                d,dd = np.vdot(state,d1@state).real,np.vdot(state,d2@state).real
                ss.append([time,np.vdot(state,hv).real,np.vdot(hv,hv).real,d,number,step])
                flct.append([time,number,number**2,d,dd,sz,sz**2,step])
                norm.append([time,before,step])
            np.testing.assert_allclose(np.loadtxt(path/'output/zvo_SS.dat'),ss,atol=4e-11,rtol=0)
            np.testing.assert_allclose(np.loadtxt(path/'output/zvo_Flct.dat'),flct,atol=4e-11,rtol=0)
            np.testing.assert_allclose(np.loadtxt(path/'output/zvo_Norm.dat'),norm,atol=4e-11,rtol=0)
            correlation.assert_nonzero_operators(operators, states_by_step)
            for step in range(5):
                correlation.check_step_file(
                    path, 'zvo_{{}}_step{}.dat'.format(step), requests,
                    operators, states_by_step[step],
                    label='{} {} k{} {} step{}'.format(
                        model, family, momentum, layout, step
                    ),
                )
            np.testing.assert_allclose(checkpoint(path,'zvo_eigenvec_final')[0],state,atol=3e-11,rtol=0)
            record=dict(line.split('=',1) for line in (path/'output/zvo_symmetry_sector.dat').read_text().splitlines())
            assert record['te_hamiltonian']=='time_dependent'
            assert record['te_plan_update']=='rebuild_with_fixed_basis'
            assert len({record['te_hamiltonian_'+str(i)] for i in range(5)}) > 1
            for step in range(5):
                assert checkpoint(path,'zvo_eigenvec_'+str(step))[1][0][21]==int(record['te_hamiltonian_'+str(step)],16)
            shutil.copytree(path/'output',path/'{}_k{}_{}'.format(family,momentum,layout))
            if model=='Spin' and momentum==0 and family=='TEOneBody':
                # A symmetry break at a LATE row must be rejected before any time output.
                shutil.rmtree(path/'output')
                (path/'output').mkdir()
                bad='====\nNTimeSteps 5\n====\n====\n====\n'
                for step,time in enumerate(times):
                    bad+='{} {}\n'.format(time,int(step==4))
                    if step==4:bad+='0 0 0 0 .3 0\n'
                (path/'drive.def').write_text(bad)
                fixture.run(path,'reject_late_'+layout,layout=layout,fail='preflight failed at step 4')
                assert not (path/'output/zvo_SS.dat').exists()
        print('{} k={} {}: every slice/vector, sparsity changes and digest PASS'.format(model,momentum,family))
    (path/'calc.def').write_text(calc);(path/'mod.def').write_text(mod);(path/'sym.def').write_text(names)


for model,length in [('Spin',6),('SpinlessFermion',6),('Hubbard',4),('tJ',4)]:
    fixture.prepare(model,length,sector_test=check)
