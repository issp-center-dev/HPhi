"""Portable expert input for fixed-N sectors; the physics oracle is in reference.py."""
from pathlib import Path
import numpy as np
from symmetry_spingc_common import definition
from symmetry_kondo_common import merged_rows
from symmetry_nconserved_reference import Case, DEFAULT, DRIVE

INPUT_MODEL = {'HubbardNConserved': 0, 'tJNConserved': 9}


def _adjoint(key):
    pairs = list(zip(key[::2], key[1::2]))
    return tuple(x for pair in reversed(pairs) for x in pair)


def hamiltonian_families(sites, parameters=None):
    """Expert rows for reference.physical_terms; Transfer carries HPhi's minus sign."""
    p = parameters or DEFAULT
    transfer, inter, intra, coulomb = [], [], [], []

    def pair(rows, key, value):
        rows.extend([(key, complex(value)), (_adjoint(key), complex(value).conjugate())])
    for j in range(sites):
        nxt = (j+1) % sites
        for s in (0, 1):
            pair(transfer, (nxt, s, j, s), p['t'])
            transfer.append(((j, s, j, s), complex(p['mu'] - (-.5 if s else .5)*p['hz'])))
        pair(transfer, (nxt, 0, j, 1), p['soc'])
        pair(transfer, (nxt, 1, j, 0), p['soc2'])
        pair(transfer, (j, 0, j, 1), p['hx'])
        pair(inter, (j, 0, j, 1, nxt, 0, nxt, 1), p['w'])
        intra.append((j, p['u']))
        coulomb.append((j, nxt, p['v']))
    families = {}
    for key, rows in (('Trans', transfer), ('InterAll', inter)):
        families[key] = [(*indices, z.real, z.imag) for indices, z in merged_rows(rows) if z]
    if any(abs(u) for _, u in intra):
        families['CoulombIntra'] = [(j, float(u)) for j, u in intra]
    if any(abs(v) for _, _, v in coulomb):
        families['CoulombInter'] = [(i, j, float(v)) for i, j, v in coulomb]
    return {key: rows for key, rows in families.items() if rows}


def drive_rows(sites, kind, amplitude=1.):
    """TEOneBody (Transfer sign) or TETwoBody (InterAll) rows for reference.drive_terms."""
    rows = []

    def pair(key, value):
        rows.extend([(key, complex(value)), (_adjoint(key), complex(value).conjugate())])
    for j in range(sites):
        nxt = (j+1) % sites
        if kind == 'onebody':
            pair((j, 0, j, 1), amplitude*DRIVE['field'])
        elif kind == 'twobody':
            pair((j, 0, j, 1, nxt, 0, nxt, 1), amplitude*DRIVE['pair'])
        else:
            raise ValueError(kind)
    return [(*indices, z.real, z.imag) for indices, z in merged_rows(rows)]


def expert_case(path, case: Case, *, method: int, options: dict | None = None,
                empty: bool = False, parameters=None) -> Path:
    """TransSym translation sector with Ncond only (no 2Sz/Nup/Ndown) unless overridden."""
    options = options or {}
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    calc = dict(CalcType=method, CalcModel=INPUT_MODEL[case.model], OutputMode=0, OutputDataHead=1)
    calc.update(options.get('CalcMod', {}))
    (path/'calc.def').write_text(''.join(f'{k} {v}\n' for k, v in calc.items()))
    mod = dict(Nsite=case.sites, Ncond=case.n, Lanczos_max=400, initial_iv=-1, exct=1,
               LanczosEps=12, LanczosTarget=1, LargeValue=100, PreCG=0)
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
    families = {} if empty else hamiltonian_families(p, parameters)
    families.update(options.get('families', {}))
    names = ['CalcMod calc.def', 'ModPara mod.def', 'LocSpin loc.def', 'TransSym group.def']
    for key, rows in families.items():
        definition(path, key+'.def', rows)
        names.append(f'{key} {key}.def')
    names += options.get('names', [])
    (path/'sym.def').write_text('\n'.join(names)+'\n')
    return path
