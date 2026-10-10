#!/usr/bin/env python3
"""Compare every P3 sector eigenvalue with the independent tensor reference."""
from pathlib import Path
import os
import sys
import numpy as np
from symmetry_kondo_reference import Case, make_reference, sector_hamiltonian
from symmetry_kondo_common import expert_case, run_case, read_manifest
from symmetry_kondo_basis import manifest_check


def main():
    hphi, probe = map(lambda p: Path(p).resolve(), sys.argv[1:])
    error = 0.
    for model in ('Kondo', 'KondoNConserved', 'KondoGC'):
        total = 0
        for k in range(3):
            case = Case(model, 3, 'block', None if model == 'KondoGC' else 2,
                        1 if model == 'Kondo' else None, k)
            ref = make_reference(case)
            h = sector_hamiltonian(ref)
            path = Path('symmetry_kondo_fulldiag')/f'{model}_k{k}'
            expert_case(path, case, method=2, options={'CalcMod': {'Solver': 0}})
            env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated',
                       HPHI_TEST_SYMMETRY_ACTION='matvec')
            run_case(path, hphi, 'fulldiag', [], env)
            actual = np.loadtxt(path/'output/zvo_energy_sector.dat', ndmin=2)
            np.testing.assert_array_equal(actual[:, 0], np.arange(len(h)))
            assert np.isfinite(actual).all()
            difference = np.max(np.abs(actual[:, 1]-np.linalg.eigvalsh(h)))
            assert difference <= 3e-8
            error = max(error, difference)
            total += len(h)
            manifest_check(read_manifest(path), case, ref)
            assert not list(path.glob('*probe*')), 'production must ignore test probe env'
        assert total == {'Kondo': 39, 'KondoNConserved': 120, 'KondoGC': 512}[model]
    print(f'Kondo FullDiag passed; maximum eigenvalue error {error:.3e}')


if __name__ == '__main__':
    main()
