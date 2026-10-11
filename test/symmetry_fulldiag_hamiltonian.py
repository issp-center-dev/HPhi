"""Sector Hamiltonian identity, canonical storage and digest through the test probe."""
from pathlib import Path
import sys

import numpy as np

import symmetry_fulldiag_common as common

MODEL_IDS = {'Hubbard': 0, 'Spin': 1, 'Kondo': 2, 'SpinGC': 4, 'KondoGC': 5,
             'SpinlessFermion': 7, 'tJ': 9, 'KondoNConserved': 12}
IDENTITY = ('model', 'nsite', 'nup', 'ndown', 'ne', 'local_site_mask', 'fixed_flags',
            'phase', 'full_dim', 'sector_dim', 'group_digest', 'sector_count',
            'sector_xor', 'sector_sum', 'basis_order')


def read_contract(case, ranks):
    contract = {key: int(value) for key, value in
                common.read_key_values(case.path/'matrix-contract.dat').items()}
    assert set(contract) == set(IDENTITY) | {'matrix_digest'}, sorted(contract)
    dim = contract['sector_dim']
    lower = np.zeros((dim, dim), complex)
    seen = np.zeros((dim, dim), bool)
    for rank in range(ranks):
        # Dump columns follow the ELPA panel split NC = ceil(D/P) on every storage.
        chunk = -(-dim // ranks)
        first, last = rank*chunk + 1, min((rank+1)*chunk, dim)
        for line in (case.path/f'matrix-contract.rank{rank}.dat').read_text().splitlines():
            row, col, real, imag = line.split()
            row, col = int(row), int(col)
            assert first <= col <= last and col <= row <= dim, (rank, row, col)
            assert not seen[row-1, col-1], (row, col)
            seen[row-1, col-1] = True
            lower[row-1, col-1] = complex(float(real), float(imag))
    assert seen[np.tril_indices(dim)].all() and not seen[np.triu_indices(dim, 1)].any()
    assert np.isfinite(lower).all()
    assert np.all(np.diag(lower).imag == 0)
    full = lower + np.tril(lower, -1).conj().T
    return contract, lower, full


def check_identity(case, contract, manifest):
    assert contract['model'] == MODEL_IDS[case.model], (case.name, contract['model'])
    assert contract['sector_dim'] == len(case.hq) == contract['sector_count']
    assert contract['full_dim'] == len(case.states) == int(manifest['full_dim'])
    assert contract['nsite'] == int(manifest['nsite'])
    kondo = case.model.startswith('Kondo')
    assert contract['phase'] == (2 if kondo else 1)
    if not kondo:
        assert contract['local_site_mask'] == contract['fixed_flags'] == 0
    else:
        assert contract['local_site_mask'] == int(manifest['local_site_mask'], 16)
    if case.model == 'SpinGC':
        assert contract['nup'] == contract['ndown'] == contract['ne'] == 0
    for key in ('fixed_nup', 'fixed_ndown', 'fixed_ne'):
        if key in manifest:
            assert contract[key[6:]] == int(manifest[key]), (case.name, key)
    assert contract['group_digest'] == int(manifest['group_digest'].rsplit(':', 1)[1], 16)
    count, xor_hash, sum_hash = manifest['sector_digest'].split(':')[1:]
    assert (contract['sector_count'], contract['sector_xor'], contract['sector_sum']) == \
        (int(count), int(xor_hash, 16), int(sum_hash, 16))
    group_order = int(manifest['group_order'])
    assert contract['basis_order'] == common.basis_order(case, group_order), case.name


def known_matrix(dim):
    h = np.diag(np.arange(1, dim+1)/8).astype(complex)
    for i in range(dim-1):
        h[i+1, i] = .25+.125j
        h[i, i+1] = .25-.125j
    return h


def contract_mode(hphi, probe_exe, launcher, ranks):
    root = common.suite_root('hamiltonian-contract')
    backends = common.solver_backends(launcher, ranks)
    negative_done = set()
    for case in common.cases(root, 'all'):
        reference = None
        for solver, label, backend_launcher in backends:
            backend_ranks = common.spingc.mpi_size(backend_launcher)
            common.run_probe(case, probe_exe, solver, backend_launcher, 'fulldiag-matrix', f'matrix_{label}')
            manifest = common.read_key_values(case.path/'output/zvo_symmetry_sector.dat')
            contract, lower, full = read_contract(case, backend_ranks)
            check_identity(case, contract, manifest)
            error = np.max(np.abs(full - case.hq))
            assert error <= 1e-10, (case.name, label, error)
            # The dump is the 17-digit binary64 matrix that was hashed.
            assert contract['matrix_digest'] == common.matrix_digest(lower), (case.name, label)
            if reference is None:
                reference = contract
            else:
                # Generated values, identity and digest are storage independent.
                assert contract == reference, (case.name, label)
            common.run_probe(case, probe_exe, solver, backend_launcher, 'fulldiag-matrix-known',
                  f'known_{label}')
            known, known_lower, _ = read_contract(case, backend_ranks)
            expected = known_matrix(len(case.hq))
            np.testing.assert_array_equal(known_lower, np.tril(expected))
            assert known['matrix_digest'] == common.matrix_digest(expected), (case.name, label)
            print(f'{case.name} {label} np={backend_ranks}: identity, |H-Hq|={error:.2e}, '
                  f'digest {contract["matrix_digest"]:016x}', flush=True)
            if case.name != 'Spin_L6_k0' or label in negative_done:
                continue
            negative_done.add(label)
            for action, message in (
                    ('fulldiag-matrix-nan', 'sector Hamiltonian has a non-finite element'),
                    ('fulldiag-matrix-nonhermitian', 'sector Hamiltonian is not Hermitian'),
                    ('fulldiag-matrix-diagonal', 'sector Hamiltonian is not Hermitian')):
                common.run_probe(case, probe_exe, solver, backend_launcher, action, f'{action}_{label}', message)
                assert not (case.path/'matrix-contract.dat').exists()
            common.run_probe(case, probe_exe, solver, backend_launcher, 'fulldiag-matrix-diagonal-small',
                  f'diagonal_small_{label}')
            small, small_lower, _ = read_contract(case, backend_ranks)
            # Within tolerance: the diagonal imaginary part is canonicalized to +0.
            np.testing.assert_array_equal(small_lower, lower)
            assert small['matrix_digest'] == contract['matrix_digest']
            if solver == 1 and backend_ranks > 1:
                common.run_probe(case, probe_exe, solver, backend_launcher, 'fulldiag-matrix-rank-mismatch',
                      f'rank_mismatch_{label}', 'sector Hamiltonian digest differs across MPI ranks')
            print(f'{case.name} {label}: non-finite / non-Hermitian / tolerance / rank checks PASS',
                  flush=True)
    assert negative_done == {label for _, label, _ in backends}, negative_done


def main():
    hphi = Path(sys.argv[1]).resolve()
    mode = sys.argv[2] if len(sys.argv) > 2 else 'contract'
    launcher, ranks = common.launcher_from_env()
    if mode == 'contract':
        contract_mode(hphi, common.probe_executable(hphi), launcher, ranks)
    else:
        raise ValueError(mode)
    print(f'sector FullDiag hamiltonian ({mode}) PASS', flush=True)


if __name__ == '__main__':
    main()
