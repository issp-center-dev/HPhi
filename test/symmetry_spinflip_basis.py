#!/usr/bin/env python3
"""Catch wrong G x Z2 orbits, Px phases, matvec, manifest keys and SpinFlipParity input gates."""
from pathlib import Path
import os
import shlex
import shutil
import sys
import numpy as np
import symmetry_spingc_common as c
from symmetry_kondo_common import run_probe, selected_layouts
from symmetry_spingc_common import mpi_size, run_case
from symmetry_spinflip_common import check_manifest, expert_case, hamiltonian_families
from symmetry_spinflip_reference import (Case, PARITIES, make_reference, physical_hamiltonian,
                                         raw_words, sector_matrix, self_test)
from symmetry_fermion_gc_common import expert_case as fermion_expert_case
from symmetry_fermion_gc_reference import Case as FermionCase

INVARIANCE_FAILED = 'Error: TransSym Hamiltonian invariance failed for normalized terms under op'
# Px-off fixtures keep the digests of HPhi before SpinFlipParity (develop
# 148308d75): k=1 translation sectors of L=4 with the default H.
LEGACY_GROUP_DIGEST = 'hphi-group-fnv1a64-v1:c58cc92ecdaea6b9'
LEGACY_SECTOR_DIGESTS = {
    'SpinGC': 'hphi-sector-multiset-v1:3:3cf952b5e7ca3bc3:f2081d519835d24d',
    'Spin': 'hphi-sector-multiset-v1:1:5aa3464dd465db45:5aa3464dd465db45'}


def families_match_oracle():
    """The expert rows written by the suites are the oracle H, including the Px-breaking extras."""
    for model, sites, extras in (('SpinGC', 4, ((), ('hz',), ('sy',), ('complex_pair',))),
                                 ('Spin', 6, ((), ('dm',), ('hz',)))):
        case = Case(model, sites, 0)
        index = raw_words(case).astype(np.int64)
        for extra in extras:
            rows = c.families_matrix(sites, hamiltonian_families(case, extra=extra))
            expected = physical_hamiltonian(case, extra=extra)
            error = np.linalg.norm((rows-expected)[np.ix_(index, index)], 2)
            assert error <= 1e-13, (model, extra, error)
    print('expert rows reproduce the oracle Hamiltonians', flush=True)


def probe_sector(root, name, case, probe, layout, ref):
    h = sector_matrix(ref)
    path = expert_case(root/name, case, method=3)
    result = run_probe(path, case, probe, action='matvec', layout=layout)
    check_manifest(result['manifest'], case, ref.basis.shape[1])
    np.testing.assert_array_equal(result['representatives'], ref.representatives)
    error = np.linalg.norm(result['matrix']-h, 2)
    assert np.isfinite(error) and error <= 1e-11*max(1, np.linalg.norm(h, 2)), (case, layout, error)
    return result, error


def sector_suite(root, hphi, probe, launcher):
    """Every momentum and parity: representatives, Q^dagger H Q, D(k,+)+D(k,-)=D(k)."""
    max_error, empty = 0., 0
    for model, sites in (('SpinGC', 3), ('SpinGC', 4), ('SpinGC', 5), ('Spin', 4), ('Spin', 6)):
        spectrum, total = [], 0
        for k in range(sites):
            control = make_reference(Case(model, sites, 0, k))
            block = []
            for parity in PARITIES:
                case = Case(model, sites, parity, k)
                ref = make_reference(case)
                dim = ref.basis.shape[1]
                total += dim
                if dim == 0:
                    # The empty G x Z2 sector is diagnosed, not run with no basis.
                    path = expert_case(root/f'{model}_L{sites}_k{k}_p{parity}_empty', case, method=3)
                    run_case(path, hphi, 'empty', launcher, dict(os.environ),
                             expected_error='Error: TransSym sector has zero basis dimension.')
                    empty += 1
                    continue
                # Both layouts at one complex momentum per size bound the CI time.
                layouts = selected_layouts() if k == 1 else selected_layouts()[:1]
                for layout in layouts:
                    _, error = probe_sector(root, f'{model}_L{sites}_k{k}_p{parity}_{layout}',
                                            case, probe, layout, ref)
                    max_error = max(max_error, error)
                block += list(np.linalg.eigvalsh(sector_matrix(ref)))
            assert len(block) == control.basis.shape[1], (model, sites, k)
            np.testing.assert_allclose(np.sort(block), np.linalg.eigvalsh(sector_matrix(control)),
                                       rtol=0, atol=1e-10)
            spectrum += block
        assert total == len(raw_words(Case(model, sites, 0))), (model, sites, total)
        print(f'{model} L={sites}: all k and p sectors, dimension and spectrum sums passed', flush=True)
    return max_error, empty


def group_suite(root, probe):
    """Px alone, with a reflection, and translation x reflection; Neel stabilizers."""
    max_error = 0.
    layout = selected_layouts()[0]
    for control in (Case('SpinGC', 3, 0, 0, 'identity'), Case('Spin', 4, 0, 0, 'identity'),
                    Case('Spin', 4, 0, 0, 'reflection', -1),
                    Case('SpinGC', 4, 0, 0, 'translation_reflection', 1),
                    Case('SpinGC', 4, 0, 2, 'translation_reflection', 1)):
        block = []
        for parity in PARITIES:
            case = Case(control.model, control.sites, parity, control.momentum,
                        control.group, control.reflection)
            ref = make_reference(case)
            name = f'{case.model}_{case.group}_k{case.momentum}_r{case.reflection}_p{parity}'
            result, error = probe_sector(root, name, case, probe, layout, ref)
            max_error = max(max_error, error)
            if case.group == 'identity':
                # Px alone halves the space exactly: |w> +- |~w>.
                assert ref.basis.shape[1]*2 == len(ref.words), case
            block += list(np.linalg.eigvalsh(sector_matrix(ref)))
        np.testing.assert_allclose(np.sort(block), np.linalg.eigvalsh(sector_matrix(
            make_reference(control))), rtol=0, atol=1e-10)
    # The Neel word is fixed by T Px; its sector survives only where the
    # stabilizer characters add up: (k=0, p=+1) and (k=pi, p=-1).
    for k, parity, present in ((0, 1, True), (0, -1, False), (2, 1, False), (2, -1, True)):
        case = Case('SpinGC', 4, parity, k)
        ref = make_reference(case)
        result, error = probe_sector(root, f'neel_k{k}_p{parity}', case, probe, layout, ref)
        assert (0b0101 in result['representatives']) == present, (k, parity)
        max_error = max(max_error, error)
    print('Px alone, reflections and Neel stabilizers passed', flush=True)
    return max_error


def legacy_suite(root, probe):
    """Without SpinFlipParity the group digest and representatives are those of HPhi 3.7."""
    layout = selected_layouts()[0]
    for model in ('SpinGC', 'Spin'):
        case = Case(model, 4, 0, 1)
        result, _ = probe_sector(root, f'legacy_{model}', case, probe, layout, make_reference(case))
        data = result['manifest']
        assert data['group_digest'] == LEGACY_GROUP_DIGEST, data
        assert data['sector_digest'] == LEGACY_SECTOR_DIGESTS[model], data
        # An explicit 0 is the same input.
        path = expert_case(root/f'legacy_{model}_zero', case, method=3,
                           options={'ModPara': {'SpinFlipParity': 0}})
        zero = run_probe(path, case, probe, action='matvec', layout=layout)
        assert zero['manifest']['group_digest'] == LEGACY_GROUP_DIGEST
        assert zero['manifest']['sector_digest'] == LEGACY_SECTOR_DIGESTS[model]
        np.testing.assert_array_equal(zero['representatives'], result['representatives'])
        np.testing.assert_array_equal(zero['matrix'], result['matrix'])
    print('Px-off group digests and representatives are unchanged', flush=True)


def rejection_suite(root, hphi, launcher, env):
    """Each gate is reached by the shortest input path."""
    def reject(name, case, message, **kwargs):
        path = expert_case(root/('reject_'+name), case, method=3, **kwargs)
        return run_case(path, hphi, 'reject', launcher, env, expected_error=message)

    gc4, spin4 = Case('SpinGC', 4, 1, 0), Case('Spin', 4, -1, 0)
    for value in (2, .5, -2):
        reject(f'value_{value}', gc4,
               'SpinFlipParity must be 0 (no Px sector), 1 (Px=+1) or -1 (Px=-1).',
               options={'ModPara': {'SpinFlipParity': value}})
    for name, mod in (('2Sz', {'2Sz': 2}), ('NupNdown', {'2Sz': None, 'Nup': 3, 'Ndown': 1})):
        reject(f'spin_{name}', spin4, 'SpinFlipParity with Spin requires 2Sz=0.',
               options={'ModPara': mod})
    # No TransSym: Px needs the sector machinery, also with an identity group.
    path = expert_case(root/'reject_raw', gc4, method=3)
    names = (path/'sym.def').read_text().replace('TransSym group.def\n', '')
    (path/'sym.def').write_text(names)
    run_case(path, hphi, 'reject', launcher, env,
             expected_error='SpinFlipParity requires TransSym with Spin or SpinGC.')
    for model in ('HubbardGC', 'SpinlessFermionGC'):
        path = fermion_expert_case(root/f'reject_{model}', FermionCase(model, 3, 0, 0), method=3,
                                   options={'ModPara': {'SpinFlipParity': 1}})
        run_case(path, hphi, 'reject', launcher, env,
                 expected_error='SpinFlipParity requires TransSym with Spin or SpinGC.')
    path = fermion_expert_case(root/'reject_Hubbard', FermionCase('HubbardGC', 2, 0, 0), method=3,
                               options={'CalcMod': {'CalcModel': 0},
                                        'ModPara': {'SpinFlipParity': -1, 'Nup': 1, 'Ndown': 1}})
    run_case(path, hphi, 'reject', launcher, env,
             expected_error='SpinFlipParity requires TransSym with Spin or SpinGC.')
    # General spin: S=1 local spins are rejected by the TransSym sector itself.
    path = expert_case(root/'reject_general_spin', Case('SpinGC', 3, 1, 0), method=3, empty=True)
    c.definition(path, 'loc.def', [(i, 2) for i in range(3)], 3)
    run_case(path, hphi, 'reject', launcher, env,
             expected_error='TransSym symmetry basis v1 supports only Spin-1/2.')
    # Hamiltonians that keep the translations but break Px. The first failing
    # operation is (e, Px), the first one of the appended half.
    for model, extra in (('SpinGC', 'hz'), ('SpinGC', 'sy'), ('SpinGC', 'complex_pair'),
                         ('Spin', 'dm'), ('Spin', 'hz')):
        case = Case(model, 4, 1, 0)
        text = reject(f'{model}_{extra}', case, INVARIANCE_FAILED, extra=(extra,))
        assert INVARIANCE_FAILED + ' 4.' in text, text
        assert 'op 4 applies the spin flip Px of SpinFlipParity' in text, text
        # The same H is a valid Px-off sector (serial LAPACK FullDiag).
        path = expert_case(root/f'accept_{model}_{extra}', Case(model, 4, 0, 0), method=2,
                           extra=(extra,))
        run_case(path, hphi, 'accept', [], dict(env, HPHI_SYMMETRY_BASIS_LAYOUT='replicated'))
    # hx keeps Px but changes Sz: never a fixed-Sz Spin sector.
    reject('spin_hx', spin4, 'Error: TransSym Hamiltonian does not conserve fixed Sz/Nup/Ndown.',
           options={'families': {'Trans': [(i, a, i, b, .3, 0.) for i in range(4)
                                           for a, b in ((1, 0), (0, 1))]}})
    print('SpinFlipParity input gates passed', flush=True)


def main():
    hphi, probe = (Path(p).resolve() for p in sys.argv[1:3])
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    env = dict(os.environ)
    # Collective rejection is covered up to np4; larger jobs reject serially.
    reject_launcher = launcher if (mpi_size(launcher) or 1) <= 4 else []
    root = Path('symmetry_spinflip_basis')
    if root.exists():
        shutil.rmtree(root)
    self_test()
    families_match_oracle()
    rejection_suite(root, hphi, reject_launcher, env)
    legacy_suite(root, probe)
    error = group_suite(root, probe)
    sector_error, empty = sector_suite(root, hphi, probe, reject_launcher)
    error = max(error, sector_error)
    assert empty > 0, 'no empty G x Z2 sector was diagnosed'
    print(f'Spin-flip basis/matvec/manifest/gates passed; maximum error {error:.3e}')


if __name__ == '__main__':
    main()
