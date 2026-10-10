#!/usr/bin/env python3
"""Every P3 momentum of every Kondo model on each available dense backend."""
from pathlib import Path
import os
import shlex
import sys
import numpy as np
from symmetry_kondo_reference import (Case, make_reference, sector_hamiltonian,
                                      make_empty_rank_reference)
from symmetry_kondo_common import (expert_case, expert_empty_rank_case, run_case,
                                   read_manifest, selected_cases)
from symmetry_kondo_basis import manifest_check
from symmetry_spingc_common import mpi_size


def main():
    hphi, probe = map(lambda p: Path(p).resolve(), sys.argv[1:])
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = mpi_size(launcher)
    selection = selected_cases()
    backends = [(0, 'LAPACK')]
    for solver, flag, label in [(1, 'HPHI_HAS_SCALAPACK', 'ScaLAPACK'), (3, 'HPHI_HAS_ELPA', 'ELPA')]:
        if os.environ.get(flag) == '1' and ranks > 1:
            backends.append((solver, label))
        else:
            print(f'{label} not exercised: backend unavailable or launcher has one rank', flush=True)
    for solver, label in backends:
        error = 0.
        # LAPACK's established contract is one rank. Distributed backends
        # must use the requested launcher and prove the actual branch ran.
        backend_launcher = [] if solver == 0 else launcher
        backend_ranks = mpi_size(backend_launcher)
        env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated',
                   HPHI_TEST_SYMMETRY_ACTION='matvec', HPHI_TEST_SYMMETRY_CAPTURE='1')
        for model in ('Kondo', 'KondoNConserved', 'KondoGC'):
            total = 0
            cases = [(k, False) for k in range(3)] if selection != 'empty-ranks' else []
            if selection in ('empty-ranks', 'all'):
                cases += [(1, True)]
            for k, tiny in cases:
                case = Case(model, 3, 'block', None if model == 'KondoGC' else 2,
                            1 if model == 'Kondo' else None, k)
                ref = make_empty_rank_reference(model) if tiny else make_reference(case)
                h = sector_hamiltonian(ref)
                path = Path('symmetry_kondo_fulldiag')/f'{model}_k{k}_solver{solver}_tiny{int(tiny)}'
                options = {'CalcMod': {'Solver': solver, 'NGPU': 0}}
                if tiny:
                    expert_empty_rank_case(path, model, method=2, options=options)
                else:
                    expert_case(path, case, method=2, options=options)
                text = run_case(path, hphi, 'fulldiag', backend_launcher, env)
                actual = np.loadtxt(path/'output/zvo_energy_sector.dat', ndmin=2)
                np.testing.assert_array_equal(actual[:, 0], np.arange(len(h)))
                assert np.isfinite(actual).all()
                difference = np.max(np.abs(actual[:, 1]-np.linalg.eigvalsh(h)))
                assert difference <= 3e-8, (model, k, label, difference)
                error = max(error, difference)
                data = read_manifest(path)
                assert int(data['solver_id']) == solver and data['basis_layout'] == 'replicated'
                assert data['output_scope'] == 'eigenvalues'
                assert f'Solver {solver}' in (path/'calc.def').read_text()
                # The runtime diagnostic distinguishes actual selected backend
                # from the input/manifest selector and catches silent fallback.
                assert f'FullDiag solver: {label}' in text, text
                if not tiny:
                    total += len(h)
                    manifest_check(data, case, ref)
                assert not list(path.glob('sector_*rank*')), 'production must ignore probe env'
                assert not list((path/'output').glob('*eigenvec*'))
                print(f'{model} k={k} {label} Solver{solver} np={backend_ranks} dim={len(h)} error={difference:.3e} passed', flush=True)
            if selection != 'empty-ranks':
                assert total == {'Kondo': 39, 'KondoNConserved': 120, 'KondoGC': 512}[model]
                # Tiny zero Hamiltonians need LAPACK: ELPA rejects an invalid
                # process grid instead of pretending to run it.
                if solver == 0:
                    case = Case(model, 1, 'block', None if model == 'KondoGC' else 1,
                                0 if model == 'Kondo' else None, 0)
                    ref = make_reference(case, empty=True)
                    path = Path('symmetry_kondo_fulldiag')/(model+'_empty')
                    expert_case(path, case, method=2, options=options, empty=True)
                    run_case(path, hphi, 'empty', backend_launcher, env)
                    actual = np.loadtxt(path/'output/zvo_energy_sector.dat', ndmin=2)
                    assert actual.shape == (len(ref.representatives), 2) and np.isfinite(actual).all()
                    np.testing.assert_array_equal(actual[:, 1], np.zeros(len(actual)))
                    manifest_check(read_manifest(path), case, ref)
        print(f'{label} maximum eigenvalue error {error:.3e}', flush=True)


if __name__ == '__main__':
    main()
