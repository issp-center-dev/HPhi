#!/usr/bin/env python3
"""Catch wrong fixed-N enumeration, group signs, matvec, manifest and input gates."""
from pathlib import Path
import os
import shlex
import shutil
import sys
import numpy as np
from symmetry_kondo_common import run_probe, selected_layouts
from symmetry_spingc_common import mpi_size, run_case
from symmetry_nconserved_common import expert_case, hamiltonian_families
from symmetry_nconserved_reference import (Case, DEFAULT, make_reference, raw_dimension,
                                           sector_matrix, self_test)


def manifest_check(data, case, ref):
    assert data['model'] == case.model, data
    assert data['fixed_quantities'] == 'ncond'
    assert int(data['fixed_ne']) == case.n
    assert 'fixed_nup' not in data and 'fixed_ndown' not in data and 'fixed_2sz' not in data
    assert int(data['nsite']) == case.sites
    assert int(data['full_dim']) == len(ref.words) == raw_dimension(case.model, case.sites, case.n)
    assert int(data['sector_dim']) == ref.basis.shape[1]
    assert int(data['momentum_index']) == case.momentum


def sector_suite(root, probe):
    """All momenta: representatives, Q^dagger H Q and dimension sums."""
    max_error = 0.
    groups = [('HubbardNConserved', 4, (1, 3, 4)), ('tJNConserved', 4, (1, 3, 4)),
              ('HubbardNConserved', 2, (4,)), ('tJNConserved', 3, (3,))]
    for model, sites, numbers in groups:
        for n in numbers:
            total = 0
            for k in range(sites):
                case = Case(model, sites, n, k)
                ref = make_reference(case)
                dim = ref.basis.shape[1]
                total += dim
                if dim == 0:
                    continue
                h = sector_matrix(ref)
                scale = max(1, np.linalg.norm(h, 2))
                # Both layouts on the largest fixture at a complex momentum; one
                # elsewhere bounds CI time.
                layouts = (selected_layouts() if (sites, n, k) == (4, 3, 1)
                           else selected_layouts()[:1])
                for layout in layouts:
                    path = expert_case(root/f'{model}_L{sites}_N{n}_k{k}_{layout}', case, method=3)
                    result = run_probe(path, case, probe, action='matvec', layout=layout)
                    manifest_check(result['manifest'], case, ref)
                    np.testing.assert_array_equal(result['representatives'], ref.representatives)
                    error = np.linalg.norm(result['matrix']-h, 2)
                    assert np.isfinite(error) and error <= 1e-11*scale, (case, layout, error)
                    max_error = max(max_error, error)
            assert total == raw_dimension(model, sites, n), (model, sites, n, total)
    return max_error


def rejection_suite(root, hphi, launcher, env):
    """Each gate is reached by the shortest input path."""
    def reject(name, case, message, **kwargs):
        path = expert_case(root/('reject_'+name), case, method=3, **kwargs)
        run_case(path, hphi, 'reject', launcher, env, expected_error=message)

    hubbard = Case('HubbardNConserved', 2, 2, 0)
    reject('hubbard_vacuum', Case('HubbardNConserved', 2, 0, 0),
           'TransSym HubbardNConserved requires 1 <= Ncond <= 4.')
    reject('hubbard_overfill', hubbard, 'TransSym HubbardNConserved requires 1 <= Ncond <= 4.',
           options={'ModPara': {'Ncond': 5}})
    reject('tj_overfill', Case('tJNConserved', 2, 2, 0),
           'TransSym tJNConserved requires 1 <= Ncond <= 2.', options={'ModPara': {'Ncond': 3}})
    reject('half_spin_numbers', hubbard, 'Nup and Ndown must be specified together.',
           options={'ModPara': {'Nup': 1}})
    reject('inconsistent_spin_numbers', hubbard,
           'TransSym Hubbard Ncond=2 conflicts with Nup+Ndown=3.',
           options={'ModPara': {'Nup': 2, 'Ndown': 1}})
    # Consistent Nup/Ndown keep the fixed-Sz model; spin mixing is not silently accepted,
    # by the reader for InterAll and by the sector term check for Transfer.
    reject('fixed_sz_with_spin_flip_interall', hubbard, 'This operator breaks Sz Symmetry',
           options={'ModPara': {'Nup': 1, 'Ndown': 1}})
    reject('fixed_sz_with_spin_mixing_transfer', hubbard,
           'Error: TransSym Hamiltonian does not conserve fixed Sz/Nup/Ndown.',
           options={'ModPara': {'Nup': 1, 'Ndown': 1}}, parameters=dict(DEFAULT, w=0))
    broken = hamiltonian_families(3)
    broken['Trans'] = [row if (row[0], row[1], row[2], row[3]) not in ((1, 0, 0, 1), (0, 1, 1, 0))
                       else (*row[:4], 1.5*row[4], 1.5*row[5]) for row in broken['Trans']]
    reject('nonuniform_spin_mixing', Case('HubbardNConserved', 3, 2, 0),
           'Error: TransSym Hamiltonian invariance failed', options={'families': broken})
    reject('anomalous_term', hubbard,
           'Error: TransSym AnomalousTerm is supported only for HubbardGC and SpinlessFermionGC.',
           options={'families': {'AnomalousTerm': [(1, 0, 0, 1, 1, .1, 0.), (0, 1, 1, 0, 0, .1, 0.)]}})
    reject('empty_sector', Case('HubbardNConserved', 2, 4, 1),
           'Error: TransSym sector has zero basis dimension.')


def canonical_suite(root, probe):
    """Ncond with 2Sz or Nup/Ndown still selects the fixed-Sz Hubbard/tJ sector."""
    parameters = dict(DEFAULT, soc=0, soc2=0, hx=0, w=0)
    for model, mod in (('HubbardNConserved', {'2Sz': 1}), ('tJNConserved', {'Nup': 2, 'Ndown': 1})):
        case = Case(model, 4, 3, 1)
        ref = make_reference(case, parameters)
        sz = ref.moments['Sz']
        rows = np.flatnonzero(abs(sz-.5) < 1e-12)
        basis = ref.basis[rows]
        # Translation keeps Sz; the fixed-Sz sector is the Sz=1/2 block of the same k.
        keep = np.flatnonzero(np.linalg.norm(basis, axis=0) > 1e-12)
        h = basis[:, keep].conj().T @ ref.hamiltonian[np.ix_(rows, rows)] @ basis[:, keep]
        path = expert_case(root/f'canonical_{model}', case, method=3, parameters=parameters,
                           options={'ModPara': mod})
        result = run_probe(path, case, probe, action='matvec', layout=selected_layouts()[0])
        data = result['manifest']
        assert data['model'] == model.replace('NConserved', ''), data
        assert 'fixed_quantities' not in data
        assert (int(data['fixed_nup']), int(data['fixed_ndown'])) == (2, 1)
        np.testing.assert_allclose(np.linalg.eigvalsh(result['matrix']), np.linalg.eigvalsh(h),
                                   rtol=0, atol=1e-11)


def main():
    hphi, probe = (Path(p).resolve() for p in sys.argv[1:3])
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    env = dict(os.environ)
    # Collective rejection is covered up to np4; larger jobs reject serially.
    reject_launcher = launcher if (mpi_size(launcher) or 1) <= 4 else []
    root = Path('symmetry_nconserved_basis')
    if root.exists():
        shutil.rmtree(root)
    # The oracle's own CAR, character and spectrum-sum checks run first.
    self_test()
    rejection_suite(root, hphi, reject_launcher, env)
    canonical_suite(root, probe)
    error = sector_suite(root, probe)
    print(f'NConserved basis/matvec/manifest/gates passed; maximum error {error:.3e}')


if __name__ == '__main__':
    main()
