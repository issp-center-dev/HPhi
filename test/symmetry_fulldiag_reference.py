"""Independent canonical U and sector H for all eight models, plus eigenvalue-only FullDiag."""
from pathlib import Path
import sys

import numpy as np

import symmetry_fulldiag_common as common
from symmetry_reference import canonical_basis, check_basis


def rejects(basis, projector):
    try:
        check_basis(projector, basis)
    except AssertionError:
        return True
    return False


def control_basis():
    """Orthonormality alone misses a column swap or a representative phase."""
    # Two orbits of a two-element group: {0, 1} and {2, 3}, trivial character.
    projector = np.zeros((4, 4), complex)
    projector[np.ix_([0, 1], [0, 1])] = .5
    projector[np.ix_([2, 3], [2, 3])] = .5
    basis = canonical_basis(projector)
    check_basis(projector, basis)
    swapped = basis[:, ::-1]
    phased = basis*np.array([1, 1j])
    for label, wrong in (('column swap', swapped), ('phase change', phased)):
        np.testing.assert_allclose(wrong.conj().T @ wrong, np.eye(2), rtol=0, atol=1e-12)
        np.testing.assert_allclose(wrong @ wrong.conj().T, projector, rtol=0, atol=1e-12)
        assert rejects(wrong, projector), label
    print('canonical U controls: column swap and phase change rejected', flush=True)


def main():
    hphi = Path(sys.argv[1]).resolve()
    mode = sys.argv[2] if len(sys.argv) > 2 else 'all'
    assert mode in ('all', 'boundary'), mode
    control_basis()
    root = common.suite_root('reference')
    totals, raw, momenta = {}, {}, {}
    for case in common.cases(root, mode):
        group = (case.model, case.length)
        dim = len(case.hq)
        totals[group] = totals.get(group, 0) + dim
        raw[group] = len(case.states)
        momenta.setdefault(group, set()).add(case.momentum)
        if case.name == 'Spin_L8_k0':
            assert raw[group] == 70 and dim == 10 and case.basis.shape[1] == 10
        common.run_fd(case, hphi, 0, {}, [], 'eigenvalues')
        output = case.path/'output'
        actual = np.loadtxt(output/'zvo_energy_sector.dat', ndmin=2)
        np.testing.assert_array_equal(actual[:, 0], np.arange(dim))
        error = np.max(np.abs(actual[:, 1] - np.linalg.eigvalsh(case.hq)))
        assert error <= 3e-8, (case.name, error)
        assert not list(output.glob('*eigenvec*')) and not list(output.glob('*cisajs*'))
        print(f'{case.name}: D={dim} raw={len(case.states)} LAPACK error={error:.3e}', flush=True)
    for group, seen in momenta.items():
        if len(seen) == group[1]:
            assert totals[group] == raw[group], (group, totals[group], raw[group])
            print(f'{group[0]} L={group[1]}: all {len(seen)} momenta cover {raw[group]} states',
                  flush=True)
    models = {group[0] for group in momenta}
    if mode == 'all':
        assert models == {'Spin', 'SpinlessFermion', 'Hubbard', 'tJ', 'SpinGC',
                          'Kondo', 'KondoNConserved', 'KondoGC'}, models
    print(f'sector FullDiag reference ({mode}): {len(models)} models PASS', flush=True)


if __name__ == '__main__':
    main()
