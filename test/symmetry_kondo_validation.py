#!/usr/bin/env python3
"""Exercise each malformed input through a serial expert entry, even under np16 CI."""
from pathlib import Path
import os
from itertools import permutations
import sys
from symmetry_kondo_reference import Case
from symmetry_kondo_common import expert_case, definition, run_case


def main():
    hphi, probe = map(lambda p: Path(p).resolve(), sys.argv[1:])
    base = Case('Kondo', 2, 'block', 2, 0, 0)
    gc = Case('KondoGC', 2, 'block', None, None, 0)
    root = Path('symmetry_kondo_validation')
    count = 0
    def negative(name, diagnostic, case=base, options=None, mutate=None, method=3):
        nonlocal count
        path = expert_case(root/name, case, method=method, options=options or {}, empty=True)
        if mutate:
            mutate(path)
        text = run_case(path, hphi, 'validation', [], dict(os.environ), expected_error=diagnostic)
        results = [f for f in (path/'output').glob('*') if any(k in f.name for k in ('energy', 'eigen', 'SS', 'Flct', 'Norm', 'symmetry_sector'))]
        assert not results, results
        messages = [line for line in text.splitlines() if diagnostic in line]
        print(f'{name}: {messages[0]}')
        count += 1
    negative('ncond_pair_conflict', 'Ncond conflicts', options={'ModPara': {'Nup': 3, 'Ndown': 3}})
    negative('nup_only', 'Nup and Ndown', options={'ModPara': {'Nup': 2}})
    for field in ('Ncond', '2Sz', 'Nup', 'Ndown'):
        negative('gc_explicit_'+field, 'Nup and Ndown' if field in ('Nup', 'Ndown') else 'does not accept explicit', gc, {'ModPara': {field: 0}})
    negative('parity', 'incompatible parity', options={'ModPara': {'2Sz': 1}})
    negative('unphysical', 'physical', options={'ModPara': {'Ncond': 8}})
    negative('general_spin', 'Spin-1/2', gc, mutate=lambda path: definition(path, 'loc.def', [(0, 2), (1, 1), (2, 0), (3, 0)], 2))
    negative('local_count', 'LocSpn', gc, mutate=lambda path: definition(path, 'loc.def', [(0, 1), (1, 0), (2, 0), (3, 0)], 2))
    def bad_group(path):
        definition(path, 'group.def', [(0, 1, 0), (1, 1, 0)]+[(g, i, i if g == 0 else i ^ 2, 1) for g in range(2) for i in range(4)], 2, 'NQPTrans')
    def too_wide(path):
        import struct
        width = struct.calcsize('L')*8//2
        nsite = width+1
        (path/'mod.def').write_text((path/'mod.def').read_text().replace('Nsite 4', f'Nsite {nsite}'))
        definition(path, 'loc.def', [(i, int(i < 2)) for i in range(nsite)], 2)
        definition(path, 'group.def', [(0, 1, 0)]+[(0, i, i, 1) for i in range(nsite)], 1, 'NQPTrans')
    negative('word_width', 'site count exceeds the state word', mutate=too_wide)
    negative('placement' , 'site type', gc, mutate=bad_group)
    for name, rows in [('local_hopping', [(0, 0, 1, 0, .2, 0), (1, 0, 0, 0, .2, 0)]),
                       ('hybridization', [(0, 0, 2, 0, .2, 0), (2, 0, 0, 0, .2, 0)])]:
        negative(name, 'onsite', gc, {'families': {'Trans': rows}})
    negative('local_pairhop', 'PairHop requires conduction sites', gc, {'families': {'PairHop': [(0, 2, .1)]}})
    # InterAll is reader scratch: mixed diagonal/offdiagonal rows must be
    # rejected before Hermitian-pair packing can overwrite the crossed row.
    crossed_rows = [(0, 0, 1, 0, 1, 0, 0, 0, .2, 0),
                    (0, 0, 0, 1, 1, 0, 1, 0, .1, 0),
                    (1, 0, 1, 0, 0, 1, 0, 0, .1, 0)]
    for order, rows in enumerate(permutations(crossed_rows)):
        negative('crossed_interall_mixed_'+str(order), 'Site component of (i, j, k, l)',
                 Case('KondoGC', 1, 'block', None, None, 0),
                 {'families': {'InterAll': list(rows)}}, method=2)
    negative('crossed_interall', 'Site component of (i, j, k, l)', gc, {'families': {'InterAll': [(0, 0, 2, 0, 2, 0, 0, 0, .2, 0)]}})
    negative('local_correlation', 'onsite', gc, {'families': {'OneBodyG': [(0, 0, 2, 0)]}})
    negative('pairlift', 'PairLift is active only in SpinGC', gc, {'families': {'PairLift': [(0, 1, .1)]}})
    negative('nbody_interall', 'unsupported term families', gc,
             {'families': {'NBodyInterAll': [(3, 0, 0, 0, 0, 1, 0, 1, 0, 2, 0, 2, 0, .2, 0)]}})
    negative('anomalous_term', 'AnomalousTerm is currently supported only for HubbardGC', gc,
             {'families': {'AnomalousTerm': [(0, 2, 0, 2, 1, .3, 0), (1, 2, 1, 2, 0, .3, 0)]}})
    negative('anomalous_g', 'AnomalousG is currently supported only for HubbardGC', gc,
             {'families': {'AnomalousG': [(0, 2, 0, 2, 1)]}})
    negative('nonhermitian' , 'Hermite', gc, {'families': {'Trans': [(2, 0, 3, 0, .2, .1)]}})
    negative('nontranslation', 'invariance', gc, {'families': {'Trans': [(0, 0, 0, 0, .2, 0)]}})
    for field in ('OutputHam', 'InputHam', 'ReStart', 'CalcSpec'):
        negative('option_'+field, 'does not support '+field if field != 'CalcSpec' else 'spectrum calculations', gc, {'CalcMod': {field: 1}}, method=2 if field == 'OutputHam' else 3)
    for field in ('OutputEigenVec', 'InputEigenVec'):
        negative('fulldiag_'+field, field, gc, {'CalcMod': {field: 1}}, method=2)
    # The new early guard must not restrict raw KondoGC InterAll input.
    raw = expert_case(root/'raw_gc_crossed_interall',
                      Case('KondoGC', 1, 'block', None, None, 0), method=2,
                      options={'families': {'InterAll': crossed_rows}}, empty=True)
    (raw/'sym.def').write_text((raw/'sym.def').read_text().replace('TransSym group.def\n', ''))
    run_case(raw, hphi, 'raw_gc', [], dict(os.environ))
    assert (raw/'output/Eigenvalue.dat').exists()
    # The TransSym-only vacuum relaxation must survive the full expert reader.
    for name, mod in [('ncond', {'Ncond': 0, '2Sz': None}),
                      ('explicit_pair', {'Ncond': None, '2Sz': None, 'Nup': 0, 'Ndown': 0})]:
        path = expert_case(root/('vacuum_'+name), Case('Kondo', 1, 'block', 0, 0, 0),
                           method=2, options={'ModPara': mod}, empty=True)
        definition(path, 'loc.def', [(0, 0), (1, 0)], 0)
        run_case(path, hphi, 'vacuum', [], dict(os.environ))
        values = (path/'output/zvo_energy_sector.dat').read_text().splitlines()
        assert len([line for line in values if line.strip() and not line.startswith('#')]) == 1
    print(f'Kondo validation passed; {count} guards, launcher=[]')


if __name__ == '__main__':
    main()
