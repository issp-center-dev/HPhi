"""Sector FullDiag I/O failures must stop every rank without a hang."""
import dataclasses
import os
from pathlib import Path
import shutil
import sys

import symmetry_fulldiag_common as common

TIMEOUT = 30


def blocked_run(case, hphi, solver, launcher, label, blocked, message, flags):
    """Pre-create an output name as a directory so that only its opener fails."""
    output = case.path/'output'
    if output.exists():
        shutil.rmtree(output)
    output.mkdir()
    (output/blocked).mkdir()
    text = common.run_fd(case, hphi, solver, flags, launcher, label, message, timeout=TIMEOUT)
    manifest = (output/'zvo_symmetry_sector.dat').read_text()
    assert 'fulldiag_status=complete' not in manifest, manifest
    return text


def green_mode(hphi, launcher, ranks):
    root = common.suite_root('failures-green')
    cases = {case.name: case for case in common.cases(root, 'all')}
    # Hubbard requests every Green family, so later TwoBodyG files exist.
    case = cases['Hubbard_L4_k0']
    assert common.attach_correlations(case, False)
    families = set(case.operators)
    assert {'one', 'two', 'three', 'four', 'six', 'nbody'} <= families, families
    scenarios = [
        ('onebody', 'zvo_cisajs_eigen0.dat', 0),
        ('twobody_later_file', 'zvo_ThreeBody_eigen0.dat', 0),
        ('twobody_last_file', 'zvo_SixBody_eigen0.dat', 0),
        ('nbody', 'zvo_NBodyG_eigen0.dat', 0),
        ('second_state', 'zvo_cisajs_eigen1.dat', 0),
        ('aggregate', 'zvo_cisajs_eigen.dat', 1),
        ('aggregate_twobody', 'zvo_cisajscktalt_eigen.dat', 1),
    ]
    backends = common.solver_backends(launcher, ranks)
    if ranks > 1 and len(backends) == 1:
        print('green failures not exercised: a multi-rank FullDiag needs ScaLAPACK or ELPA, '
              'which this build lacks', flush=True)
        common.require_backends({})
    for solver, label, backend_launcher in backends:
        if ranks > 1 and solver == 0:
            continue  # LAPACK is serial; the collective cases need several ranks
        for name, blocked, aggregate in scenarios:
            blocked_run(case, hphi, solver, backend_launcher, f'{name}_{label}', blocked,
                        'sector correlation output open failed',
                        {'OutputEigenVec': 0, 'OutputGreenFormat': aggregate})
            print(f'{label} np={common.spingc.mpi_size(backend_launcher)} {name}: '
                  f'{blocked} rejected on every rank', flush=True)


# (stage, rank, backends, diagnostic, flags). Allocation and checkpoint
# stages fail on rank 1, root-only output stages on rank 0.
IO_STAGES = (
    ('vector_alloc', 1, ('ScaLAPACK', 'ELPA'),
     'sector FullDiag state buffer allocation failed', {'OutputEigenVec': 1}),
    ('matrix_alloc', 1, ('ELPA',),
     'sector Hamiltonian column buffer allocation failed', {'OutputEigenVec': 1}),
    ('checkpoint_write', 1, ('ScaLAPACK', 'ELPA'),
     'symmetry checkpoint write failed', {'OutputEigenVec': 1}),
    ('checkpoint_close', 1, ('ScaLAPACK', 'ELPA'),
     'symmetry checkpoint write failed', {'OutputEigenVec': 1}),
    ('checkpoint_rename', 1, ('ScaLAPACK', 'ELPA'),
     'symmetry checkpoint publish failed', {'OutputEigenVec': 1}),
    ('green_close', 0, ('ScaLAPACK', 'ELPA'),
     'sector correlation output close failed', {'OutputEigenVec': 0}),
    ('matrix_write', 0, ('ScaLAPACK', 'ELPA'),
     'sector Hamiltonian output write failed', {'OutputHam': 1}),
    ('matrix_close', 0, ('ScaLAPACK', 'ELPA'),
     'sector Hamiltonian output write failed', {'OutputHam': 1}),
    ('matrix_rename', 0, ('ScaLAPACK', 'ELPA'),
     'sector Hamiltonian output publish failed', {'OutputHam': 1}),
)
# The matrix writer has a replicated (ScaLAPACK) and a column-panel (ELPA)
# path; every available one is exercised.
MATRIX_OUTPUT_STAGES = ('matrix_write', 'matrix_close', 'matrix_rename')


def io_mode(hphi, launcher, ranks):
    """Injected rank-local I/O and allocation failures end every rank.

    matrix_alloc is the ELPA column-panel buffer; it needs an ELPA build."""
    if ranks < 2:
        print('io failures not exercised: injection on rank 1 needs an MPI launcher with '
              'at least two ranks', flush=True)
        common.require_backends({})  # a job that requires backends must not skip
        return
    backends = {label: solver for solver, label, _ in common.solver_backends(launcher, ranks)}
    if set(backends) == {'LAPACK'}:
        print('io failures not exercised: a multi-rank FullDiag needs ScaLAPACK or ELPA, '
              'which this build lacks', flush=True)
        common.require_backends({})
        return
    probe = common.probe_executable(hphi)
    root = common.suite_root('failures-io')
    cases = {case.name: case for case in common.cases(root, 'all')}
    case = cases['Hubbard_L4_k0']
    assert common.attach_correlations(case, False)
    exercised = []
    for stage, rank, labels, message, flags in IO_STAGES:
        usable = [label for label in labels if label in backends]
        if not usable:
            print(f'{stage} not exercised: needs {"/".join(labels)} in this build', flush=True)
            continue
        for label in usable if stage in MATRIX_OUTPUT_STAGES else usable[-1:]:
            io_stage(case, hphi, probe, backends, launcher, ranks, stage, rank, label,
                     message, flags)
            exercised.append((stage, label))
    assert exercised, 'no I/O failure stage was exercised'
    if os.environ.get('HPHI_HAS_ELPA') == '1':
        assert ('matrix_alloc', 'ELPA') in exercised
    assert {stage for stage, _ in exercised} >= set(MATRIX_OUTPUT_STAGES), exercised


def io_stage(case, hphi, probe, backends, launcher, ranks, stage, rank, label, message, flags):
    """One injected failure on one rank, then the production control."""
    env = {'HPHI_TEST_FD_FAIL': stage, 'HPHI_TEST_FD_FAIL_RANK': str(rank)}
    output = case.path/'output'
    if output.exists():
        shutil.rmtree(output)
    sentinel = b'%% a matrix from an earlier run\n'
    if stage in MATRIX_OUTPUT_STAGES:
        output.mkdir()
        (output/MATRIX).write_bytes(sentinel)
    common.run_fd(case, probe, backends[label], {**flags, 'OutputGreenFormat': 0},
                  launcher, f'{stage}_{label}', message, timeout=TIMEOUT, env_extra=env)
    manifest = (output/'zvo_symmetry_sector.dat').read_text()
    assert 'fulldiag_status=complete' not in manifest, (stage, manifest)
    if stage.startswith('checkpoint'):
        # The failed state is never published, and no temporary remains.
        assert not (output/f'zvo_eigenvec_0_rank_{rank}.dat').exists(), stage
        assert not list(output.glob('*.part')), stage
    if stage in MATRIX_OUTPUT_STAGES:
        # The earlier file under the final name is untouched, no temporary is left.
        assert (output/MATRIX).read_bytes() == sentinel, stage
        assert not list(output.glob('*.part')), stage
    print(f'{stage} on rank {rank} with {label} np={ranks}: every rank stopped', flush=True)
    # The same environment does not reach production HPhi.
    if output.exists():
        shutil.rmtree(output)
    common.run_fd(case, hphi, backends[label], {**flags, 'OutputGreenFormat': 0},
                  launcher, f'{stage}_{label}_production', timeout=TIMEOUT, env_extra=env)
    assert common.read_manifest(output)['fulldiag_status'] == 'complete'


MATRIX = 'zvo_ham_sector.mtx'
HEADER_LINES = 1 + len(common.MATRIX_KEYS)  # first line and the HPHI keys


def _key_line(lines, key):
    return next(i for i, line in enumerate(lines) if line.startswith(f'% HPHI {key}='))


def _replace_key(key, value):
    def mutate(lines):
        lines[_key_line(lines, key)] = f'% HPHI {key}={value}'
    return mutate


def _changed_value(key, value):
    if key == 'format':
        return 'sector-hamiltonian-v2'
    if key == 'scalar':
        return 'complex-binary32'
    if key == 'model':
        return 'SpinGC' if value != 'SpinGC' else 'Spin'
    if key in common.MATRIX_HEX_KEYS:
        return f'{int(value, 16) ^ 1:016x}'
    return str(int(value) + 1)


def matrix_negatives(metadata, dim, nnz):
    """(label, mutation of the file lines, expected diagnostic)."""
    dims = HEADER_LINES  # index of the dimension line
    first = dims + 1     # first record
    out = []
    for key in common.MATRIX_KEYS:
        message = {'format': 'unsupported format', 'scalar': 'unsupported scalar'}.get(
            key, f'metadata {key} does not match this sector')
        out.append((f'key_{key}', _replace_key(key, _changed_value(key, metadata[key])), message))

    def edit(index, text):
        """text None: a copy of the last line (a duplicate record)."""
        def mutate(lines):
            lines[index] = lines[-1] if text is None else text
        return mutate

    def insert(index, text):
        def mutate(lines):
            lines.insert(index, text)
        return mutate

    def delete(index):
        def mutate(lines):
            del lines[index]
        return mutate

    def record(row, col, value='0.5 0'):
        return edit(first, f'{row} {col} {value}')

    out += [
        ('header_real', edit(0, '%%MatrixMarket matrix coordinate real symmetric'),
         'first line must be'),
        ('header_missing', delete(0), 'first line must be'),
        ('key_missing', lambda lines: lines.pop(_key_line(lines, 'basis_order')),
         'missing metadata key basis_order'),
        ('key_duplicate', lambda lines: lines.insert(2, lines[_key_line(lines, 'nsite')]),
         'duplicate metadata key nsite'),
        ('key_unknown', insert(2, '% HPHI solver=LAPACK'), 'unknown metadata key solver'),
        ('key_malformed', insert(2, '% HPHI nsite'), 'malformed metadata line'),
        ('model_unknown', _replace_key('model', 'Heisenberg'), 'invalid metadata value for model'),
        ('decimal_sign', _replace_key('nsite', '+' + metadata['nsite']),
         'invalid metadata value for nsite'),
        ('hex_width', _replace_key('group_digest', metadata['group_digest'][1:]),
         'invalid metadata value for group_digest'),
        ('decimal_overflow', _replace_key('full_dim', '18446744073709551616'),
         'invalid metadata value for full_dim'),
        ('dimension_rows', edit(dims, f'{dim+1} {dim+1} {nnz}'), 'does not match sector dimension'),
        ('dimension_square', edit(dims, f'{dim} {dim+1} {nnz}'), 'does not match sector dimension'),
        ('dimension_tokens', edit(dims, f'{dim} {dim}'), 'dimension line must hold three'),
        ('dimension_extra', edit(dims, f'{dim} {dim} {nnz} 1'), 'dimension line must hold three'),
        ('dimension_negative', edit(dims, f'{dim} {dim} -1'), 'dimension line must hold three'),
        ('nnz_triangle', edit(dims, f'{dim} {dim} {dim*(dim+1)//2+1}'),
         'more records than the lower triangle'),
        ('row_zero', record(0, 1), 'index out of range'),
        ('row_large', record(dim+1, 1), 'index out of range'),
        ('col_zero', record(1, 0), 'index out of range'),
        ('upper', record(1, 2), 'upper-triangle record'),
        ('duplicate', lambda lines: lines.insert(first+1, lines[first]), 'duplicated record'),
        # Column D has the highest owner rank; nnz still matches the records.
        ('duplicate_replaced', edit(-2, None), 'duplicated record'),
        ('nan', record(1, 1, 'nan 0'), 'non-finite coefficient'),
        ('inf', record(2, 1, '0 inf'), 'non-finite coefficient'),
        ('diagonal_imag', record(1, 1, '0.5 1e-300'), 'nonzero diagonal imaginary part'),
        ('missing_record', lambda lines: lines.pop(), 'fewer records than declared'),
        ('extra_record', lambda lines: lines.append(lines[first]), 'unexpected content after'),
        ('extra_token', record(1, 1, '0.5 0 7'), 'malformed record'),
        # fgets/strlen can hide the suffix at EOF instead of rejecting it.
        ('nul_last_record', lambda lines: lines.__setitem__(-1, lines[-1] + '\0 77'),
         'NUL'),
        ('nul_trailer', lambda lines: lines.append('\0 99 99 invalid extra record'), 'NUL'),
        ('bad_number', record(1, 1, '0.5x 0'), 'malformed record'),
        ('comment_in_records', insert(first, '% a comment between records'), 'malformed record'),
        ('long_line', record(1, 1, '0.5 ' + '0'*600), 'too long'),
        ('line_limit_lf', lambda lines: lines.__setitem__(-1, lines[-1].ljust(511)), 'too long'),
        ('line_limit_crlf', lambda lines: lines.__setitem__(-1, lines[-1].ljust(510) + '\r'),
         'too long'),
    ]
    return out


def matrix_case(cases):
    """Spin L6 k1: complex couplings and D >= 2 records on both triangles."""
    case = cases['Spin_L6_k1']
    assert len(case.hq) >= 3, len(case.hq)
    return case


def matrix_run(case, hphi, solver, launcher, label, text, message, flags=None):
    output = case.path/'output'
    if output.exists():
        shutil.rmtree(output)
    output.mkdir()
    if text is not None:
        (output/MATRIX).write_text(text)
    common.run_fd(case, hphi, solver, flags or {'InputHam': 1, 'OutputEigenVec': 1}, launcher,
                  label, message, timeout=TIMEOUT)
    assert not (output/'zvo_energy_sector.dat').exists(), label
    assert not list(output.glob('*eigenvec*')), label
    manifest = (output/'zvo_symmetry_sector.dat').read_text()
    assert 'fulldiag_status=complete' not in manifest, (label, manifest)


def raw_elpa_guard(case, hphi, launcher, label):
    """Without TransSym, multi-rank ELPA keeps rejecting the raw Hamiltonian
    I/O, which needs a replicated matrix; the guard now runs after TransSym."""
    raw = dataclasses.replace(case, path=case.path.parent/f'{case.name}_raw')
    if raw.path.exists():
        shutil.rmtree(raw.path)
    shutil.copytree(case.path, raw.path, ignore=shutil.ignore_patterns('output', '*.log'))
    names = [line for line in (raw.path/'sym.def').read_text().splitlines()
             if not line.split() or line.split()[0] != 'TransSym']
    (raw.path/'sym.def').write_text('\n'.join(names) + '\n')
    for flag in ('OutputHam', 'InputHam'):
        (raw.path/'output').mkdir(exist_ok=True)
        common.run_fd(raw, hphi, 3, {flag: 1}, launcher, f'raw_{flag}_{label}',
                      'incompatible with OutputHam/InputHam', timeout=TIMEOUT)
    print(f'raw FullDiag {label}: multi-rank ELPA OutputHam/InputHam still rejected', flush=True)


def cross_chunk_duplicate(case, hphi, solver, launcher, label):
    """More than one record chunk (4096): every lower element is written
    explicitly, and the (D,D) record of the last column owner is repeated in
    the first chunk in place of another record."""
    output = case.path/'output'
    if output.exists():
        shutil.rmtree(output)
    common.run_fd(case, hphi, solver, {'OutputHam': 1}, launcher, f'export_{label}')
    metadata, h = common.read_sector_matrix(output/MATRIX)
    dim = len(h)
    records = common.matrix_records(
        h, explicit_zeros={(row, col) for col in range(1, dim+1) for row in range(col, dim+1)})
    assert len(records) > 3*4096 and records[-1][:2] == (dim, dim), (dim, len(records))
    records = [records[-1]] + records[1:]
    text_file = case.path/'cross_chunk.mtx'
    common.write_matrix_records(text_file, metadata, dim, records)
    matrix_run(case, hphi, solver, launcher, f'cross_chunk_{label}', text_file.read_text(),
               'duplicated record')


def matrix_mode(hphi, launcher, ranks):
    """Every malformed sector matrix stops all ranks before diagonalization."""
    root = common.suite_root('failures-matrix')
    cases = {case.name: case for case in common.cases(root, 'all')}
    case = matrix_case(cases)
    # A sector of the same model and dimension: only its identity differs.
    other = next(c for c in cases.values() if c.model == case.model and c is not case
                 and len(c.hq) == len(case.hq))
    exercised = []
    for solver, label, backend_launcher in common.solver_backends(launcher, ranks):
        n = common.spingc.mpi_size(backend_launcher)
        output = case.path/'output'
        if output.exists():
            shutil.rmtree(output)
        common.run_fd(case, hphi, solver, {'OutputHam': 1}, backend_launcher, f'export_{label}')
        valid = (output/MATRIX).read_text()
        metadata, h = common.read_sector_matrix(output/MATRIX)
        dim, nnz = len(h), len(common.matrix_records(h))
        assert nnz >= 3 and h[1, 0] != 0, (dim, nnz)
        common.run_fd(other, hphi, solver, {'OutputHam': 1}, backend_launcher, f'export_{label}')
        matrix_run(case, hphi, solver, backend_launcher, f'other_sector_{label}',
                   (other.path/'output'/MATRIX).read_text(), 'does not match this sector')
        for name, mutate, message in matrix_negatives(metadata, dim, nnz):
            lines = valid.splitlines()
            mutate(lines)
            matrix_run(case, hphi, solver, backend_launcher, f'{name}_{label}',
                       '\n'.join(lines) + '\n', message)
        matrix_run(case, hphi, solver, backend_launcher, f'input_missing_{label}', None,
                   'could not open sector Hamiltonian input')
        # The temporary name is taken by a directory: the root cannot open it.
        output = case.path/'output'
        shutil.rmtree(output)
        output.mkdir()
        (output/(MATRIX + '.part')).mkdir()
        common.run_fd(case, hphi, solver, {'OutputHam': 1}, backend_launcher,
                      f'output_open_{label}', 'could not open sector Hamiltonian output',
                      timeout=TIMEOUT)
        assert not (output/MATRIX).exists()
        assert 'fulldiag_status=complete' not in (output/'zvo_symmetry_sector.dat').read_text()
        exercised.append((label, n))
        print(f'{label} np={n}: {len(matrix_negatives(metadata, dim, nnz)) + 4} malformed or '
              f'unavailable matrix files rejected on every rank', flush=True)
        cross_chunk_duplicate(cases['KondoGC_P3_k0'], hphi, solver, backend_launcher, label)
        if label == 'ELPA' and n > 1:
            raw_elpa_guard(case, hphi, backend_launcher, label)
    assert exercised


def main():
    hphi = Path(sys.argv[1]).resolve()
    mode = sys.argv[2] if len(sys.argv) > 2 else 'green'
    launcher, ranks = common.launcher_from_env()
    if mode == 'green':
        green_mode(hphi, launcher, ranks)
    elif mode == 'io':
        io_mode(hphi, launcher, ranks)
    elif mode == 'matrix':
        matrix_mode(hphi, launcher, ranks)
    else:
        raise ValueError(mode)
    print(f'sector FullDiag failures ({mode}) PASS', flush=True)


if __name__ == '__main__':
    main()
