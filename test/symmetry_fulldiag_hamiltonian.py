"""Sector Hamiltonian identity, canonical storage and digest through the test
probe (contract), and the Matrix Market OutputHam/InputHam files (io)."""
from pathlib import Path
import shutil
import sys

import numpy as np

import symmetry_fulldiag_common as common
from symmetry_fulldiag_common import W

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


def check_hermitian_boundaries(case, probe_exe, solver, label, launcher):
    # Exercise PrepareSymmetryHamiltonian without asking a dense solver to
    # diagonalize values near DBL_MAX. The expected canonical entries are
    # fixed independently of the production tolerance calculation.
    failures = []
    for suffix, accepted in (
            ('diagonal-below', True), ('diagonal-above', False),
            ('pair-below', True), ('pair-above', False),
            ('large-hermitian', True), ('large-near', True),
            ('large-outside', False), ('large-nonhermitian', False),
            ('large-diagonal', False), ('large-diagonal-below', True),
            ('large-diagonal-above', False)):
        action = f'fulldiag-matrix-{suffix}'
        try:
            common.run_probe(case, probe_exe, solver, launcher, action, f'{suffix}_{label}',
                             None if accepted else 'sector Hamiltonian is not Hermitian')
            if not accepted:
                assert not (case.path/'matrix-contract.dat').exists()
                continue
            contract, lower, _ = read_contract(case, common.spingc.mpi_size(launcher))
            expected = known_matrix(len(case.hq))
            if suffix.startswith('diagonal-'):
                expected[0, 0] = 1.0
            elif suffix.startswith('pair-'):
                expected[1, 0] = 1.07
            else:
                expected[0, 0] = 1.4e308
                expected[1, 0] = complex(1.4e308, 1.4e308)
            np.testing.assert_array_equal(lower, np.tril(expected))
            assert contract['matrix_digest'] == common.matrix_digest(expected)
        except AssertionError as error:
            failures.append(f'{suffix}: {error}')
    assert not failures, '\n'.join(failures)
    print(f'{case.name} {label}: 11 Hermitian numeric boundary cases PASS', flush=True)


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
            check_hermitian_boundaries(case, probe_exe, solver, label, backend_launcher)
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


MATRIX = 'zvo_ham_sector.mtx'
GREEN_GLOBS = ('*cisajs*', '*ThreeBody*', '*FourBody*', '*SixBody*', '*NBodyG*')
ATOL = 1e-8


def metadata_contract(metadata):
    """The 15 identity fields of a matrix file as the probe contract integers."""
    return {key: MODEL_IDS[metadata[key]] if key == 'model'
            else int(metadata[key], 16 if key in common.MATRIX_HEX_KEYS else 10)
            for key in IDENTITY}


def fresh_output(case):
    output = case.path/'output'
    if output.exists():
        shutil.rmtree(output)
    output.mkdir()
    return output


def digest_text(digest):
    return f'hphi-sector-matrix-fnv1a64-v1:{digest:016x}'


def export(case, hphi, solver, launcher, label):
    """OutputHam writes the canonical matrix and stops before diagonalizing."""
    output = fresh_output(case)
    common.run_fd(case, hphi, solver, {'OutputHam': 1}, launcher, f'export_{label}')
    metadata, exported = common.read_sector_matrix(output/MATRIX)
    assert not list(output.glob('*.part'))
    assert not (output/'zvo_energy_sector.dat').exists()
    assert not list(output.glob('*eigenvec*'))
    for pattern in GREEN_GLOBS:
        assert not list(output.glob(pattern)), pattern
    manifest = common.read_manifest(output)
    assert manifest['output_scope'] == 'requested' and manifest['requested_outputs'] == 'hamiltonian'
    assert manifest['matrix_source'] == 'generated' and manifest['fulldiag_status'] == 'complete'
    assert manifest['completed_outputs'] == 'hamiltonian'
    digest = common.matrix_digest(exported)
    assert manifest['matrix_digest'] == digest_text(digest), (manifest['matrix_digest'], digest)
    assert 'eigenvalue_file' not in manifest, manifest
    # Signed zeros are written as +0.
    records = [line.split() for line in (output/MATRIX).read_text().splitlines()
               if not line.startswith('%')][1:]
    assert not [r for r in records if '-0' in r[2:]], records
    assert metadata['format'] == 'sector-hamiltonian-v1'
    assert metadata['scalar'] == 'complex-binary64'
    check_identity(case, metadata_contract(metadata), manifest)
    error = np.max(np.abs(exported - case.hq))
    assert error <= 1e-10, (case.name, label, error)
    return metadata, exported, digest, (output/MATRIX).read_bytes()


def import_matrix(case, hphi, solver, launcher, label, h, writer):
    """InputHam: the file written by writer(path) is the H that is diagonalized."""
    output = fresh_output(case)
    writer(output/MATRIX)
    text = common.run_fd(case, hphi, solver, {'InputHam': 1, 'OutputEigenVec': 1},
                         launcher, f'import_{label}')
    assert f'FullDiag solver: {label.split("_")[0]}' in text, text
    ranks = common.spingc.mpi_size(launcher)
    dim = len(h)
    energies = np.loadtxt(output/'zvo_energy_sector.dat', ndmin=2)[:, 1]
    vectors, headers = common.read_states(output, dim, ranks)
    common.check_eigensystem(h, energies, vectors)
    digest = common.matrix_digest(h)
    for parts in headers:
        for header in parts:
            assert int(header[W.HAMILTONIAN]) == digest
    manifest = common.read_manifest(output)
    requested = 'eigenvalues,correlations,eigenvectors' if case.operators else 'eigenvalues,eigenvectors'
    assert manifest['output_scope'] == 'requested' and manifest['requested_outputs'] == requested
    assert manifest['matrix_source'] == 'input' and manifest['fulldiag_status'] == 'complete'
    assert manifest['completed_outputs'] == requested
    assert manifest['matrix_digest'] == digest_text(digest), (manifest['matrix_digest'], digest)
    for family, requests in case.operators.items():
        actual = common.read_green(output, family, dim, False, [row for row, _ in requests])
        expected = np.array([[np.vdot(v, oq @ v) for _, oq in requests] for v in vectors.T])
        np.testing.assert_allclose(actual.real, expected.real, rtol=0, atol=ATOL)
        np.testing.assert_allclose(actual.imag, expected.imag, rtol=0, atol=ATOL)
    return energies, digest


def edits(case, hphi, solver, launcher, label, metadata, exported, digest):
    """External edits keep the identity comments and change the provenance."""
    shifted = exported + .125*np.eye(len(exported))
    energies, shifted_digest = import_matrix(
        case, hphi, solver, launcher, f'{label}_shift', shifted,
        lambda path: common.write_sector_matrix(path, metadata, shifted))
    np.testing.assert_allclose(energies, np.linalg.eigvalsh(exported) + .125, rtol=0, atol=3e-8)
    assert shifted_digest != digest
    offdiagonal = exported.copy()
    offdiagonal[1, 0] += .07+.03j
    offdiagonal[0, 1] = np.conj(offdiagonal[1, 0])
    _, offdiagonal_digest = import_matrix(
        case, hphi, solver, launcher, f'{label}_offdiagonal', offdiagonal,
        lambda path: common.write_sector_matrix(path, metadata, offdiagonal))
    assert offdiagonal_digest not in (digest, shifted_digest)
    print(f'{case.name} {label}: delta-I and off-diagonal edits PASS', flush=True)


def precision(case, hphi, solver, launcher, label, metadata, exported):
    """17-digit values, explicit zeros, any record order, CRLF and nnz=0."""
    h = exported.copy()
    h[1, 0] = h[0, 1] = 0
    h[2, 0], h[0, 2] = 1e-15+3e-16j, 1e-15-3e-16j
    h[2, 2] = 1/3
    records = common.matrix_records(h, explicit_zeros={(2, 1)})
    assert (2, 1, 0j) in records

    def writer(path):
        common.write_matrix_records(path, metadata, len(h), records[::-1])
        lines = path.read_text().splitlines()
        lines[1:1] = ['% general comment kept by any Matrix Market reader', '']
        path.write_bytes(('\r\n'.join(lines) + '\r\n% trailing comment\r\n\r\n').encode())
    import_matrix(case, hphi, solver, launcher, f'{label}_precision', h, writer)

    def unterminated(path):
        # The longest accepted line (511 bytes) ends the file without a line end.
        common.write_sector_matrix(path, metadata, h)
        lines = path.read_text().splitlines()
        lines[-1] = lines[-1].ljust(511)
        path.write_text('\n'.join(lines))
    import_matrix(case, hphi, solver, launcher, f'{label}_unterminated', h, unterminated)
    zero = np.zeros_like(h)
    import_matrix(case, hphi, solver, launcher, f'{label}_zero', zero,
                  lambda path: common.write_matrix_records(path, metadata, len(h), []))
    print(f'{case.name} {label}: 17-digit/explicit-zero/reversed/CRLF/unterminated/nnz=0 inputs PASS',
          flush=True)


def chunks(case, hphi, solver, launcher, label, metadata, exported):
    """Every lower element written explicitly, in a fixed shuffled order:
    more than three record chunks of 4096 reach every rank."""
    dim = len(exported)
    records = common.matrix_records(
        exported, explicit_zeros={(row, col) for col in range(1, dim+1) for row in range(col, dim+1)})
    assert len(records) > 3*4096, len(records)
    order = np.random.default_rng(20261011).permutation(len(records))
    records = [records[i] for i in order]
    import_matrix(case, hphi, solver, launcher, f'{label}_chunks', exported,
                  lambda path: common.write_matrix_records(path, metadata, dim, records))
    print(f'{case.name} {label}: {len(records)} records in {-(-len(records)//4096)} chunks PASS',
          flush=True)


def io_mode(hphi, launcher, ranks):
    selection = common.selected_cases()
    root = common.suite_root(f'hamiltonian-io-{selection}')
    backends = common.solver_backends(launcher, ranks)
    labels = [label for _, label, _ in backends]
    exercised, edited, precise, chunked = {}, set(), set(), set()
    for case in common.cases(root, selection):
        if not common.attach_correlations(case, False):
            case.operators = {}
        dim = len(case.hq)
        exports = {label: export(case, hphi, solver, backend_launcher, label)
                   for solver, label, backend_launcher in backends}
        # Generated values, identity and digest are rank and storage independent.
        assert len({data[3] for data in exports.values()}) == 1, case.name
        for index, (solver, label, backend_launcher) in enumerate(backends):
            n = common.spingc.mpi_size(backend_launcher)
            # Each backend reads the file of another rank count/backend.
            source = labels[(index+1) % len(labels)]
            metadata, exported, digest, text = exports[source]
            if solver == 3 and dim < common.grid_max(n):
                output = fresh_output(case)
                (output/MATRIX).write_bytes(text)
                common.run_fd(case, hphi, solver, {'InputHam': 1, 'OutputEigenVec': 1},
                              backend_launcher, f'grid_reject_{label}', 'smaller than the process grid')
                assert not list(output.glob('*eigenvec*'))
                continue
            _, imported = import_matrix(case, hphi, solver, backend_launcher, label, exported,
                                        lambda path: path.write_bytes(text))
            # The 17-digit file restores the exact binary64 matrix.
            assert imported == digest
            exercised.setdefault(label, []).append((case.name, dim, n))
            print(f'{case.name} {label} np={n}: export, import of the {source} file, '
                  f'digest {digest:016x} PASS', flush=True)
            if selection != 'all':
                continue
            if dim >= 2 and (case.model, label) not in edited:
                edits(case, hphi, solver, backend_launcher, label, metadata, exported, digest)
                edited.add((case.model, label))
            if dim >= 3 and label not in precise:
                precision(case, hphi, solver, backend_launcher, label, metadata, exported)
                precise.add(label)
            if dim*(dim+1)//2 > 3*4096 and label not in chunked:
                chunks(case, hphi, solver, backend_launcher, label, metadata, exported)
                chunked.add(label)
    for label, runs in exercised.items():
        print(f'{label} exercised: ' + ', '.join(f'{name}(D={dim},np={n})' for name, dim, n in runs),
              flush=True)
    if selection == 'all':
        models = {model for model, _ in edited}
        assert len(models) == 8 and precise == chunked == set(exercised), \
            (sorted(models), precise, chunked)
    common.require_backends(exercised)


def main():
    hphi = Path(sys.argv[1]).resolve()
    mode = sys.argv[2] if len(sys.argv) > 2 else 'contract'
    launcher, ranks = common.launcher_from_env()
    if mode == 'contract':
        contract_mode(hphi, common.probe_executable(hphi), launcher, ranks)
    elif mode == 'io':
        io_mode(hphi, launcher, ranks)
    else:
        raise ValueError(mode)
    print(f'sector FullDiag hamiltonian ({mode}) PASS', flush=True)


if __name__ == '__main__':
    main()
