#!/usr/bin/env python3
"""Catch wrong GC/parity enumeration, group signs, matvec, manifest and input gates."""
from pathlib import Path
import os
import shlex
import shutil
import sys
import numpy as np
from symmetry_kondo_common import run_probe, selected_layouts
from symmetry_spingc_common import mpi_size, run_case
from symmetry_fermion_gc_common import expert_case, hamiltonian_families
from symmetry_fermion_gc_reference import (Case, PARITIES, equivalent_pair_rows,
                                           make_reference, pair_rows, raw_words,
                                           sector_matrix, self_test, width)

PARITY_NAME = {0: 'full', 1: 'even', -1: 'odd'}


def manifest_check(data, case, ref, pairs):
    assert data['model'] == case.model, data
    assert data['fixed_quantities'] == 'none', data
    assert data['particle_parity'] == PARITY_NAME[case.parity], data
    for key in ('fixed_ne', 'fixed_nup', 'fixed_ndown', 'fixed_2sz', 'fixed_ncond'):
        assert key not in data, (key, data)
    assert int(data['nsite']) == case.sites
    assert int(data['full_dim']) == len(ref.words) == (1 << (width(case.model)*case.sites)) >> (case.parity != 0)
    assert int(data['sector_dim']) == ref.basis.shape[1]
    assert int(data['momentum_index']) == case.momentum
    assert data['hamiltonian_digest'].startswith('hphi-parsed-hamiltonian-fnv1a64-v4:'), data
    assert data['term_scope'].endswith(f' anomalous:{len(pairs or [])}'), data['term_scope']


def sector_suite(root, probe):
    """All momenta and parities: representatives, Q^dagger H Q and dimension sums."""
    max_error = 0.
    for model, sites in (('HubbardGC', 3), ('HubbardGC', 4), ('SpinlessFermionGC', 4),
                         ('SpinlessFermionGC', 5)):
        pairs = pair_rows(model, sites)
        spectra = {}
        for parity in PARITIES:
            total, spectrum = 0, []
            for k in range(sites):
                case = Case(model, sites, parity, k)
                ref = make_reference(case, pairs=pairs)
                dim = ref.basis.shape[1]
                total += dim
                h = sector_matrix(ref)
                scale = max(1, np.linalg.norm(h, 2))
                spectrum += list(np.linalg.eigvalsh(h))
                # Both layouts at a complex momentum of the larger fixtures;
                # one elsewhere bounds CI time.
                layouts = (selected_layouts() if (sites, k) in ((4, 1), (5, 2))
                           else selected_layouts()[:1])
                for layout in layouts:
                    name = f'{model}_L{sites}_{PARITY_NAME[parity]}_k{k}_{layout}'
                    path = expert_case(root/name, case, method=3, pairs=pairs)
                    result = run_probe(path, case, probe, action='matvec', layout=layout)
                    manifest_check(result['manifest'], case, ref, pairs)
                    np.testing.assert_array_equal(result['representatives'], ref.representatives)
                    error = np.linalg.norm(result['matrix']-h, 2)
                    assert np.isfinite(error) and error <= 1e-11*scale, (case, layout, error)
                    max_error = max(max_error, error)
            assert total == len(raw_words(Case(model, sites, parity, 0))), (model, sites, parity)
            spectra[parity] = np.sort(spectrum)
        # Each momentum block splits into its even and odd parts.
        np.testing.assert_allclose(np.sort(np.r_[spectra[1], spectra[-1]]), spectra[0],
                                   rtol=0, atol=1e-10)
        print(f'{model} L={sites} full/even/odd sectors passed', flush=True)
    return max_error


def equivalent_input_suite(root, probe):
    """A rewritten pair input of the same H passes the invariance check and gives Q^dagger H Q."""
    max_error = 0.
    for case in (Case('HubbardGC', 4, 1, 1), Case('SpinlessFermionGC', 5, -1, 2)):
        pairs = pair_rows(case.model, case.sites)
        rewritten = equivalent_pair_rows(pairs)
        ref = make_reference(case, pairs=pairs)
        h = sector_matrix(ref)
        path = expert_case(root/f'equivalent_{case.model}', case, method=3, pairs=rewritten)
        result = run_probe(path, case, probe, action='matvec', layout=selected_layouts()[0])
        manifest_check(result['manifest'], case, ref, rewritten)
        np.testing.assert_array_equal(result['representatives'], ref.representatives)
        error = np.linalg.norm(result['matrix']-h, 2)
        assert np.isfinite(error) and error <= 1e-11*max(1, np.linalg.norm(h, 2)), (case, error)
        max_error = max(max_error, error)
    print('rewritten pair inputs passed', flush=True)
    return max_error


def rejection_suite(root, hphi, launcher, env):
    """Each gate is reached by the shortest input path."""
    def reject(name, case, message, **kwargs):
        path = expert_case(root/('reject_'+name), case, method=3, **kwargs)
        run_case(path, hphi, 'reject', launcher, env, expected_error=message)

    def drop_transsym(path):
        names = (path/'sym.def').read_text().splitlines()
        (path/'sym.def').write_text('\n'.join(n for n in names if not n.startswith('TransSym'))+'\n')

    hubbard, spinless = Case('HubbardGC', 3, 1, 0), Case('SpinlessFermionGC', 4, -1, 0)
    # The reader already requires Nup and Ndown together.
    for name, mod in (('Ncond', {'Ncond': 2}), ('2Sz', {'2Sz': 0}),
                      ('NupNdown', {'Nup': 1, 'Ndown': 1})):
        reject(f'hubbard_{name}', hubbard,
               'HubbardGC TransSym does not accept explicit 2Sz/Nup/Ndown/Ncond.',
               options={'ModPara': mod})
    reject('spinless_ncond', spinless,
           'SpinlessFermionGC TransSym does not accept explicit 2Sz/Nup/Ndown/Ncond.',
           options={'ModPara': {'Ncond': 2}})
    for value in (2, .5, -2):
        reject(f'parity_{value}', hubbard,
               'ParticleParity must be 0 (all N), 1 (even N) or -1 (odd N).',
               options={'ModPara': {'ParticleParity': value}})
    reject('parity_canonical', hubbard,
           'ParticleParity requires TransSym with HubbardGC or SpinlessFermionGC.',
           options={'CalcMod': {'CalcModel': 0},
                    'ModPara': {'Nup': 2, 'Ndown': 1, 'Ncond': 3}})
    # Raw GC inputs keep their space; a parity there is an error, 0 is not.
    path = expert_case(root/'reject_parity_raw', hubbard, method=3)
    drop_transsym(path)
    run_case(path, hphi, 'reject', launcher, env,
             expected_error='ParticleParity requires TransSym with HubbardGC or SpinlessFermionGC.')
    path = expert_case(root/'reject_raw_spinless_fulldiag', Case('SpinlessFermionGC', 4, 0, 0), method=2)
    drop_transsym(path)
    run_case(path, hphi, 'reject', [], env,
             expected_error='SpinlessFermion / SpinlessFermionGC are not currently supported for FullDiag')
    reject('spinless_hund', spinless, 'Hund is not active in SpinlessFermion/SpinlessFermionGC.',
           options={'families': {'Hund': [(0, 1, .1)]}})
    broken = hamiltonian_families(Case('SpinlessFermionGC', 4, 0, 0))
    # One bond and its adjoint change, so the input stays Hermitian.
    broken['Trans'] = [row if row[:4] not in ((1, 0, 0, 0), (0, 0, 1, 0))
                       else (*row[:4], 1.5*row[4], 1.5*row[5]) for row in broken['Trans']]
    reject('nonuniform_hopping', spinless, 'Error: TransSym Hamiltonian invariance failed',
           options={'families': broken})

    # Pair terms: the reader's row checks, then the sector term invariance.
    def pairs(rows):
        return {'families': {'AnomalousTerm': rows}}
    reject('pair_spinless_spin', spinless, 'Error: Spin index of AnomalousTerm is incorrect.',
           options=pairs([(1, 0, 1, 1, 0, .1, 0.), (0, 1, 0, 0, 1, .1, 0.)]))
    reject('pair_same_orbital', hubbard,
           'Error: AnomalousTerm cannot use the same fermion operator twice in one pair.',
           options=pairs([(1, 0, 1, 0, 1, .1, 0.), (0, 0, 1, 0, 1, .1, 0.)]))
    reject('pair_type', hubbard, 'Error: AnomalousTerm type must be 0 or 1.',
           options=pairs([(2, 0, 0, 0, 1, .1, 0.), (0, 0, 1, 0, 0, .1, 0.)]))
    reject('pair_missing_adjoint', hubbard,
           'Error: AnomalousTerm terms must appear as adjacent Hermite pairs.',
           options=pairs([(1, 0, 0, 0, 1, .1, 0.)]))
    reject('pair_wrong_adjoint', hubbard,
           'Error: AnomalousTerm Hermite pair has inconsistent coefficients.',
           options=pairs([(1, 0, 0, 0, 1, .1, .2), (0, 0, 1, 0, 0, .1, .2)]))
    # Rows come in Hermitian pairs; change one bond with its adjoint.
    uniform = pair_rows('SpinlessFermionGC', 4)
    nonuniform = [row if n//2 != 2 else (*row[:5], 1.5*row[5]) for n, row in enumerate(uniform)]
    reject('pair_nonuniform', spinless, 'Error: TransSym Hamiltonian invariance failed',
           pairs=nonuniform)
    reject('pair_nonuniform_rewritten', spinless, 'Error: TransSym Hamiltonian invariance failed',
           pairs=equivalent_pair_rows(nonuniform))


def main():
    hphi, probe = (Path(p).resolve() for p in sys.argv[1:3])
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    env = dict(os.environ)
    # Collective rejection is covered up to np4; larger jobs reject serially.
    reject_launcher = launcher if (mpi_size(launcher) or 1) <= 4 else []
    root = Path('symmetry_fermion_gc_basis')
    if root.exists():
        shutil.rmtree(root)
    # The oracle's own CAR, character, parity and spectrum-sum checks run first.
    self_test()
    rejection_suite(root, hphi, reject_launcher, env)
    error = max(sector_suite(root, probe), equivalent_input_suite(root, probe))
    print(f'Fermion GC basis/matvec/manifest/gates passed; maximum error {error:.3e}')


if __name__ == '__main__':
    main()
