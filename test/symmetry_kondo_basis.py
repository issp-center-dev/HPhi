#!/usr/bin/env python3
"""Catch wrong Kondo gates, physical phases, matvec coefficients and space identity."""
from pathlib import Path
import os
import struct
import sys
import os
import numpy as np
from symmetry_kondo_reference import Case, make_reference, sector_hamiltonian, _translate, apply_operator
from symmetry_kondo_common import expert_case, run_probe


def manifest_check(data, case, ref):
    p = case.cells
    local = range(p) if case.layout == 'block' else range(0, 2*p, 2)
    mask = sum(1 << i for i in local)
    flags = {'Kondo': 3, 'KondoNConserved': 1, 'KondoGC': 0}[case.model]
    ne = p+case.ncond if flags else 0
    up = (ne+case.sz2)//2 if flags == 3 else 0
    down = (ne-case.sz2)//2 if flags == 3 else 0
    values = [2 if flags == 3 else 12 if flags == 1 else 5, 2*p, mask, flags, up, down, ne, 2]
    digest = 14695981039346656037
    for byte in b'hphi-kondo-space-fnv1a64-v1\0'+struct.pack('<8Q', *values):
        digest = ((digest ^ byte)*1099511628211) & ((1 << 64)-1)
    assert data['model'] == case.model
    assert int(data['n_local_spin']) == p and int(data['n_conduction_sites']) == p
    assert int(data['local_site_mask'], 16) == mask
    assert data['translation_convention'] == 'kondo-physical-v1'
    assert data['particle_number_convention'] == 'total_fermions_including_local_spins'
    assert data['basis_space_digest'] == f'hphi-kondo-space-fnv1a64-v1:{digest:016x}'
    assert int(data['full_dim']) == len(ref.words)
    assert int(data['sector_dim']) == ref.basis.shape[1]
    assert data['fixed_quantities'] == {3: 'ncond,2sz', 1: 'ncond', 0: 'none'}[flags]
    expected = {} if not flags else dict(fixed_ncond=case.ncond, fixed_ne=ne)
    if flags == 3:
        expected.update(fixed_2sz=case.sz2, fixed_nup=up, fixed_ndown=down)
    for key in ('fixed_ncond', 'fixed_ne', 'fixed_2sz', 'fixed_nup', 'fixed_ndown'):
        assert (int(data[key]) == expected[key]) if key in expected else key not in data
    ordered = ['n_local_spin', 'n_conduction_sites', 'local_site_mask', 'translation_convention',
               'particle_number_convention', 'basis_space_digest', 'fixed_quantities']
    ordered += ['fixed_ncond', 'fixed_2sz', 'fixed_nup', 'fixed_ndown', 'fixed_ne'] if flags == 3 else ['fixed_ncond', 'fixed_ne'] if flags else []
    ordered += ['momentum_index']
    keys = list(data)
    positions = [keys.index(k) for k in ordered]
    assert positions == list(range(positions[0], positions[0]+len(positions)))


def family_check(root, probe, model):
    """Catch runtime rejection of supported families and crossed-Exchange projection."""
    case = Case(model, 2, 'block', None if model == 'KondoGC' else 2,
                0 if model == 'Kondo' else None, 0)
    ref = make_reference(case, empty=True)
    families = {
        'CoulombIntra': [(i, .41) for i in range(4)],
        'CoulombInter': [(j, j+2, .13) for j in range(2)],
        'Hund': [(j, j+2, .17) for j in range(2)],
        'Exchange': [(j, j+2, .19) for j in range(2)]+[(j, j, .37) for j in range(2)],
        'Ising': [(0, 1, .23)],
        'PairHop': [(2, 3, .29)],
    }
    # Use onsite local matrix units and conduction CAR, never crossed local strings.
    terms = []
    for j in range(2):
        c = j+2
        terms.append((.41, ((c, 0, c, 0), (c, 1, c, 1))))
        for a in (0, 1):
            for b in (0, 1):
                terms.append((.13, ((j, a, j, a), (c, b, c, b))))
            terms.append((-.17, ((j, a, j, a), (c, a, c, a))))
            terms.append((-.19, ((j, a, j, 1-a), (c, 1-a, c, a))))
    for a in (0, 1):
        for b in (0, 1):
            terms.append((.23*(.5-a)*(.5-b), ((0, a, 0, a), (1, b, 1, b))))
    for i, j in ((2, 3), (3, 2)):
        terms.append((.29, ((i, 0, j, 0), (i, 1, j, 1))))
    expected = np.zeros((ref.basis.shape[1],)*2, complex)
    for col, v in enumerate(ref.basis.T):
        out = sum(value*apply_operator(case, factors, v) for value, factors in terms)
        expected[:, col] = ref.basis.conj().T @ out
    path = expert_case(root/('families_'+model), case, method=3,
                       options={'families': families}, empty=True)
    actual = run_probe(path, case, probe, action='matvec', layout=os.environ.get('HPHI_SYMMETRY_BASIS_LAYOUT', 'distributed'))['matrix']
    error = np.linalg.norm(actual-expected, 2)
    assert error <= 1e-11*max(1, np.linalg.norm(expected, 2)), (model, error)
    return error


def main():
    from symmetry_kondo_common import selected_cases, selected_layouts
    from symmetry_kondo_empty import run_empty_rank_suite
    selection = selected_cases()
    if selection in ('empty-ranks', 'all'):
        args = [Path(arg).resolve() for arg in sys.argv[1:]]
        run_empty_rank_suite('basis', args[0], args[1] if len(args) > 1 else None)
    if selection == 'empty-ranks':
        return

    hphi, probe = map(lambda p: Path(p).resolve(), sys.argv[1:])
    root = Path('symmetry_kondo_basis')
    max_error = 0.
    max_character = 0.
    for model in ('Kondo', 'KondoNConserved', 'KondoGC'):
        for p, sites in ((3, 'block'), (4, 'block'), (4, 'alternating'), (2, 'block'), (1, 'block')):
            total = 0
            empty = p == 1
            for k in range(p):
                case = Case(model, p, sites, None if model == 'KondoGC' else 1 if empty else 2,
                            (0 if p != 3 else 1) if model == 'Kondo' else None, k)
                ref = make_reference(case, empty=empty)
                dim = ref.basis.shape[1]
                total += dim
                # Orbit supports are disjoint; normalize each column independently.
                assert np.max(np.count_nonzero(ref.basis, axis=1)) <= 1
                assert np.linalg.norm(np.sum(abs(ref.basis)**2, axis=0)-1) <= 1e-12
                character_error = np.linalg.norm(_translate(case, ref.basis, 1)-
                    np.exp(-2j*np.pi*k/p)*ref.basis)
                assert character_error <= 1e-12
                max_character = max(max_character, character_error)
                # P4 only needs a single vector; all columns for P3 and empty H.
                action = 'apply' if p == 4 else 'matvec'
                vector = np.arange(1, dim+1)+1j*np.arange(dim, 0, -1)
                vector /= np.linalg.norm(vector)
                expected = ref.basis.conj().T @ ref.apply_h(ref.basis @ vector)
                h = sector_hamiltonian(ref) if p != 4 or k == 1 else None
                for layout in selected_layouts():
                    path = root/f'{model}_P{p}_{sites}_k{k}_{layout}'
                    expert_case(path, case, method=3, options={}, empty=empty)
                    result = run_probe(path, case, probe, action=action, layout=layout,
                                       vector=vector if action == 'apply' else None)
                    manifest_check(result['manifest'], case, ref)
                    np.testing.assert_array_equal(result['representatives'], ref.representatives)
                    if action == 'matvec':
                        error = np.linalg.norm(result['matrix']-h, 2)
                        scale = max(1, np.linalg.norm(h, 2))
                        if empty:
                            assert np.count_nonzero(result['matrix']) == 0
                    else:
                        error = np.linalg.norm(result['vector']-expected)
                        scale = max(1, np.linalg.norm(h, 2)) if h is not None else 1
                    assert np.isfinite(error) and error <= 1e-11*scale, (case, error)
                    max_error = max(max_error, error)
            assert total == len(ref.words)
            if empty:
                assert total == {'Kondo': 2, 'KondoNConserved': 4, 'KondoGC': 8}[model]
    for model in ('Kondo', 'KondoNConserved', 'KondoGC'):
        max_error = max(max_error, family_check(root, probe, model))
    print(f'Kondo basis/matvec/manifest passed; maximum error {max_error:.3e}, character error {max_character:.3e}')


if __name__ == '__main__':
    main()
