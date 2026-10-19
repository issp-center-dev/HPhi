"""Portable expert input and v2 vectors of fermion GC sectors; the oracle is in reference.py."""
from pathlib import Path
import numpy as np
from symmetry_spingc_common import definition
from symmetry_nconserved_common import hamiltonian_families as hubbard_families
from symmetry_fermion_gc_reference import Case, SPINLESS_DEFAULT

INPUT_MODEL = {'HubbardGC': 3, 'SpinlessFermionGC': 8}
# FIXED_FLAGS word: bit2 marks a ParticleParity, bit3 the odd one.
PARITY_FLAGS = {0: 0, 1: 4, -1: 12}


def spinless_families(sites, parameters=None):
    """Expert rows for reference.spinless_terms; Transfer carries HPhi's minus sign."""
    p = parameters or SPINLESS_DEFAULT
    t = complex(p['t'])
    transfer, coulomb = [], []
    for j in range(sites):
        nxt = (j+1) % sites
        transfer += [(nxt, 0, j, 0, t.real, t.imag), (j, 0, nxt, 0, t.real, -t.imag),
                     (j, 0, j, 0, float(p['mu']), 0.)]
        coulomb.append((j, nxt, float(p['v'])))
    families = {'Trans': transfer}
    if p['v']:
        families['CoulombInter'] = coulomb
    return families


def hamiltonian_families(case: Case, parameters=None):
    if case.model == 'HubbardGC':
        return hubbard_families(case.sites, parameters)
    return spinless_families(case.sites, parameters)


def anomalous_rows(rows):
    """AnomalousTerm lines from reference.pair_rows (type i si j sj Re Im)."""
    return [(kind, i, si, j, sj, complex(v).real, complex(v).imag)
            for kind, i, si, j, sj, v in rows]


def expert_case(path, case: Case, *, method: int, options: dict | None = None,
                empty: bool = False, parameters=None, pairs=None) -> Path:
    """TransSym translation sector of a GC fermion model, ParticleParity if nonzero."""
    options = options or {}
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    calc = dict(CalcType=method, CalcModel=INPUT_MODEL[case.model], OutputMode=0, OutputDataHead=1)
    calc.update(options.get('CalcMod', {}))
    (path/'calc.def').write_text(''.join(f'{k} {v}\n' for k, v in calc.items()))
    mod = dict(Nsite=case.sites, Lanczos_max=400, initial_iv=-1, exct=1,
               LanczosEps=12, LanczosTarget=1, LargeValue=100, PreCG=0)
    if case.parity:
        mod['ParticleParity'] = case.parity
    mod.update(options.get('ModPara', {}))
    mod = {k: v for k, v in mod.items() if v is not None}
    (path/'mod.def').write_text('====\nModel_Parameters 0\n====\n====\n====\n'
        'CDataFileHead zvo\nCParaFileHead zqp\n====\n'
        + ''.join(f'{k} {v}\n' for k, v in mod.items()))
    p = case.sites
    definition(path, 'loc.def', [(i, 0) for i in range(p)], 0)
    rows = [(g, np.cos(-2*np.pi*case.momentum*g/p), np.sin(-2*np.pi*case.momentum*g/p))
            for g in range(p)]
    rows += [(g, j, (j+g) % p, 1) for g in range(p) for j in range(p)]
    definition(path, 'group.def', rows, p, 'NQPTrans')
    with (path/'group.def').open('a') as f:
        f.write(f'# MomentumIndex {case.momentum}\n')
    families = {} if empty else hamiltonian_families(case, parameters)
    if pairs:
        families['AnomalousTerm'] = anomalous_rows(pairs)
    families.update(options.get('families', {}))
    names = ['CalcMod calc.def', 'ModPara mod.def', 'LocSpin loc.def', 'TransSym group.def']
    for key, rows in families.items():
        definition(path, key+'.def', rows)
        names.append(f'{key} {key}.def')
    names += options.get('names', [])
    (path/'sym.def').write_text('\n'.join(names)+'\n')
    return path


def read_gc_checkpoint(path, case: Case | None = None):
    """Version 2, phase 1: zero mask, no fixed N/Nup/Ndown, parity in the flags."""
    data = Path(path).read_bytes()
    assert len(data) >= 30*8, 'truncated GC checkpoint header'
    header = np.frombuffer(data[:30*8], dtype='<u8').copy()
    assert tuple(header[:4]) == (0x0a31565349485048, 2, 1, 128), tuple(header[:4])
    assert tuple(header[6:9]) == (0, 0, 0) and header[28] == 0, header
    if case is not None:
        assert header[4] == INPUT_MODEL[case.model] and header[5] == case.sites, header[4:6]
        assert header[29] == PARITY_FLAGS[case.parity], header[29]
    assert len(data) == 30*8+int(header[15])*16
    payload = np.frombuffer(data[30*8:], dtype='<c16').copy()
    assert np.isfinite(payload).all()
    return header, payload


def join_gc_vectors(paths, dimension, case: Case | None = None):
    parts = sorted((read_gc_checkpoint(p, case) for p in paths),
                   key=lambda part: (part[0][14], part[0][12]))
    result, offset = np.empty(dimension, complex), 0
    for header, vector in parts:
        assert header[10] == dimension and header[11] == len(parts)
        assert header[14] == offset
        result[offset:offset+len(vector)] = vector
        offset += len(vector)
    assert offset == dimension and np.isfinite(result).all()
    return result
