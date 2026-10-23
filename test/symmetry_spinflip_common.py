"""Portable expert input of Px sectors; the physics oracle is in symmetry_spinflip_reference.py."""
from pathlib import Path
import symmetry_spingc_common as c
from symmetry_spinflip_reference import (SPIN_DEFAULT, SPINGC_DEFAULT, SY, SZ, Case,
                                         input_group, nearest_bonds, second_bonds)

INPUT_MODEL = {'Spin': 1, 'SpinGC': 4}
GROUP_DOMAIN = {False: 'hphi-group-fnv1a64-v1', True: 'hphi-group-spinflip-fnv1a64-v1'}


def sy_sz_rows(i, j, value):
    """value*(Sy_i Sz_j + Sz_i Sy_j) as InterAll rows with their adjoints."""
    rows = []
    for cc in range(2):
        rows += c.interall_pair((i, 1, i, 0, j, cc, j, cc), value*SY[1, 0]*SZ[cc, cc])
        rows += c.interall_pair((i, cc, i, cc, j, 1, j, 0), value*SZ[cc, cc]*SY[1, 0])
    return rows


def hamiltonian_families(case: Case, parameters=None, extra=()):
    """Expert rows of reference.physical_hamiltonian, including the Px-breaking extras."""
    n = case.sites
    families = {}
    if case.model == 'SpinGC':
        q = dict(SPINGC_DEFAULT, **(parameters or {}))
        bonds = nearest_bonds(n)
        families['Ising'] = [(i, j, q['jz']) for i, j in bonds]
        families['Exchange'] = [(i, j, q['jxy']/2) for i, j in bonds]
        families['PairLift'] = [(i, j, q['pair_lift']) for i, j in bonds]
        # Transfer carries a minus sign: -hx Sx = -(hx/2)(E10 + E01).
        families['Trans'] = [(i, a, i, b, q['hx']/2, 0.) for i in range(n) for a, b in ((1, 0), (0, 1))]
        families['InterAll'] = [row for i, j in bonds for row in sy_sz_rows(i, j, q['yz'])]
    else:
        q = dict(SPIN_DEFAULT, **(parameters or {}))
        couplings = [(i, j, q['jz'], q['jxy']) for i, j in nearest_bonds(n)]
        couplings += [(i, j, q['j2'], q['j2']) for i, j in second_bonds(n)]
        families['Ising'] = [(i, j, jz) for i, j, jz, _ in couplings]
        families['Exchange'] = [(i, j, jxy/2) for i, j, _, jxy in couplings]
    for name in extra:
        if name == 'hz':
            families.setdefault('Trans', [])
            families['Trans'] += [(i, a, i, a, -.21*SZ[a, a].real, 0.) for i in range(n) for a in range(2)]
        elif name == 'sy':
            families.setdefault('Trans', [])
            families['Trans'] += [(i, a, i, b, -.21*SY[a, b].real, -.21*SY[a, b].imag)
                                  for i in range(n) for a, b in ((1, 0), (0, 1))]
        elif name in ('dm', 'complex_pair'):
            spins = (1, 0, 0, 1) if name == 'dm' else (1, 0, 1, 0)
            families.setdefault('InterAll', [])
            for i, j in nearest_bonds(n):
                families['InterAll'] += c.interall_pair(
                    (i, spins[0], i, spins[1], j, spins[2], j, spins[3]), .13j)
        else:
            raise ValueError(name)
    return {key: rows for key, rows in families.items() if rows}


def group_rows(case: Case):
    rows = input_group(case)
    out = [(g, complex(chi).real, complex(chi).imag) for g, (_, chi) in enumerate(rows)]
    out += [(g, i, target, 1) for g, (perm, _) in enumerate(rows) for i, target in enumerate(perm)]
    return out, len(rows)


def expert_case(path, case: Case, *, method: int, options: dict | None = None,
                parameters=None, extra=(), empty=False) -> Path:
    """TransSym sector of Spin(2Sz=0)/SpinGC; SpinFlipParity when the case has a parity."""
    options = options or {}
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    calc = dict(CalcType=method, CalcModel=INPUT_MODEL[case.model], OutputMode=0, OutputDataHead=1)
    calc.update(options.get('CalcMod', {}))
    (path/'calc.def').write_text(''.join(f'{k} {v}\n' for k, v in calc.items()))
    mod = dict(Nsite=case.sites, Lanczos_max=400, initial_iv=-1, exct=1,
               LanczosEps=12, LanczosTarget=1, LargeValue=100, PreCG=0)
    if case.model == 'Spin':
        mod['2Sz'] = 0
    if case.parity:
        mod['SpinFlipParity'] = case.parity
    mod.update(options.get('ModPara', {}))
    mod = {k: v for k, v in mod.items() if v is not None}
    (path/'mod.def').write_text('====\nModel_Parameters 0\n====\n====\n====\n'
        'CDataFileHead zvo\nCParaFileHead zqp\n====\n'
        + ''.join(f'{k} {v}\n' for k, v in mod.items()))
    c.definition(path, 'loc.def', [(i, 1) for i in range(case.sites)], case.sites)
    rows, order = group_rows(case)
    c.definition(path, 'group.def', rows, order, 'NQPTrans')
    if case.group.startswith('translation'):
        with (path/'group.def').open('a') as f:
            f.write(f'# MomentumIndex {case.momentum}\n')
    families = {} if empty else hamiltonian_families(case, parameters, extra)
    families.update(options.get('families', {}))
    names = ['CalcMod calc.def', 'ModPara mod.def', 'LocSpin loc.def', 'TransSym group.def']
    for key, values in families.items():
        if key == 'InterAll':
            c.validate_interall(values)
        c.definition(path, key+'.def', values)
        names.append(f'{key} {key}.def')
    names += options.get('names', [])
    (path/'sym.def').write_text('\n'.join(names)+'\n')
    return path


def check_manifest(data, case: Case, dim: int):
    """Model, fixed quantities, Px keys and the group digest domain of a sector manifest."""
    assert data['model'] == case.model, data
    assert int(data['nsite']) == case.sites and int(data['sector_dim']) == dim, data
    order = len(input_group(case))
    if case.model == 'SpinGC':
        assert data['fixed_quantities'] == 'none', data
    else:
        assert data['fixed_2sz'] == '0', data
    if case.parity:
        assert data['spin_flip_operator'] == 'product_sigma_x', data
        assert data['spin_flip_parity'] == f'{case.parity:+d}', data
        assert int(data['input_group_order']) == order and int(data['group_order']) == 2*order, data
    else:
        for key in ('spin_flip_operator', 'spin_flip_parity', 'input_group_order'):
            assert key not in data, (key, data)
        assert int(data['group_order']) == order, data
    assert data['group_digest'].startswith(GROUP_DOMAIN[bool(case.parity)]+':'), data
    if case.group.startswith('translation'):
        assert int(data['momentum_index']) == case.momentum, data
