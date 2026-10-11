"""Independent sector FullDiag fixtures, checkpoint decoding and run helpers.

Every suite is started as ``script.py HPhi mode``: symmetry_general_terms
reads sys.argv[1] as the HPhi executable when it is imported.
"""
from dataclasses import dataclass, field
from enum import IntEnum
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

import numpy as np

import symmetry_general_terms as general
import symmetry_spingc_common as spingc
from symmetry_kondo_common import expert_case, expert_empty_rank_case
from symmetry_kondo_reference import (Case, make_empty_rank_reference,
                                      make_reference, sector_hamiltonian)
from symmetry_reference import canonical_basis, check_basis


@dataclass
class FDCase:
    name: str
    model: str
    path: Path
    states: np.ndarray
    basis: np.ndarray
    hq: np.ndarray
    # family -> [(input integer row, sector operator matrix), ...]
    operators: dict[str, list[tuple[tuple[int, ...], np.ndarray]]]
    length: int = 0
    momentum: int = 0
    # Kondo reference case; None for the empty-rank KondoGC geometry.
    kondo: Case | None = None


class CheckpointWord(IntEnum):
    MAGIC = 0; VERSION = 1; PHASE = 2; SCALAR = 3
    MODEL = 4; NSITE = 5; NUP = 6; NDOWN = 7; NE = 8
    FULL_DIM = 9; DIM = 10; RANKS = 11; RANK = 12; LAYOUT = 13
    OFFSET = 14; LOCAL_DIM = 15; GROUP = 16
    SECTOR_COUNT = 17; SECTOR_XOR = 18; SECTOR_SUM = 19; ORDER = 20
    HAMILTONIAN = 21; METHOD = 22; STATE = 23; STEP = 24; TIME = 25
    DATA_XOR = 26; DATA_SUM = 27; LOCAL_SITE_MASK = 28; FIXED_FLAGS = 29


W = CheckpointWord
CHECKPOINT_MAGIC = 0x0a31565349485048
KONDO_MODELS = ('Kondo', 'KondoNConserved', 'KondoGC')


def suite_root(label: str) -> Path:
    """A fresh directory below the CTest working directory; never reused."""
    root = Path(tempfile.mkdtemp(prefix=f'sector-fd-{label}-', dir=Path.cwd()))
    print(f'sector FullDiag {label} root: {root}', flush=True)
    return root


def launcher_from_env() -> tuple[list[str], int]:
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = spingc.mpi_size(launcher)
    assert ranks is not None and ranks >= 1, ('MPIRUN names no process count', launcher)
    return launcher, ranks


def block_range(dim: int, ranks: int, rank: int) -> tuple[int, int]:
    """SymmetryBlockRange: the first remainder ranks own one extra row."""
    base, remainder = divmod(dim, ranks)
    start = base*rank + min(rank, remainder)
    return start, base + (rank < remainder)


def _canonical(root, model, length, momenta, result, nup=2, ndown=1):
    general.ROOT = root/f'generator_{model}_{length}'
    general.ROOT.mkdir()

    def capture(path, model, length, momentum, states, raw, projector):
        # prepare() rewrites the same directory for the next momentum.
        name = f'{model}_L{length}_k{momentum}'
        target = root/name
        shutil.copytree(path, target)
        basis = canonical_basis(projector)
        check_basis(projector, basis)
        result.append(FDCase(name, model, target, np.asarray(states), basis,
                             basis.conj().T @ raw @ basis, {}, length, momentum))

    general.prepare(model, length, nup, ndown, sector_test=capture,
                    sector_momenta=momenta)


def _spingc(root, nsite, result):
    families = spingc.mixed_families(nsite)
    h = spingc.mixed_hamiltonian(nsite)
    np.testing.assert_allclose(spingc.families_matrix(nsite, families), h,
                               rtol=0, atol=1e-12)
    perms = [[(i+g) % nsite for i in range(nsite)] for g in range(nsite)]
    for k in range(nsite):
        chars = np.exp(-2j*np.pi*k*np.arange(nsite)/nsite)
        basis, _, projector = spingc.sector_basis(nsite, perms, chars)
        check_basis(projector, basis)
        name = f'SpinGC_L{nsite}_k{k}'
        spingc.write_case(root/name, nsite, perms, chars, families, 2, {})
        result.append(FDCase(name, 'SpinGC', root/name, np.arange(1 << nsite),
                             basis, basis.conj().T @ h @ basis, {}, nsite, k))


def _kondo_basis(reference):
    check_basis(None, reference.basis)
    rows = np.searchsorted(reference.words, reference.representatives)
    first = [np.flatnonzero(abs(column) > 1e-12)[0] for column in reference.basis.T]
    np.testing.assert_array_equal(rows, first)
    return reference.basis


def _kondo(root, result):
    for model in KONDO_MODELS:
        for k in range(3):
            case = Case(model, 3, 'block', None if model == 'KondoGC' else 2,
                        1 if model == 'Kondo' else None, k)
            reference = make_reference(case)
            name = f'{model}_P3_k{k}'
            expert_case(root/name, case, method=2, options={})
            result.append(FDCase(name, model, root/name, reference.words,
                                 _kondo_basis(reference), sector_hamiltonian(reference),
                                 {}, 3, k, case))


def _kondo_empty_rank(root, result):
    for model in KONDO_MODELS:
        reference = make_empty_rank_reference(model)
        case = None if model == 'KondoGC' else Case(
            model, 2, 'block', 2, 0 if model == 'Kondo' else None, 1)
        name = f'{model}_empty_rank'
        expert_empty_rank_case(root/name, model, method=2, options={})
        result.append(FDCase(name, model, root/name, reference.words,
                             _kondo_basis(reference), sector_hamiltonian(reference),
                             {}, 2, 1, case))


def build_cases(root: Path, selection: str) -> list[FDCase]:
    assert selection in ('all', 'boundary'), selection
    result: list[FDCase] = []
    if selection == 'all':
        for model, length in (('Spin', 6), ('SpinlessFermion', 6),
                              ('Hubbard', 4), ('tJ', 4)):
            _canonical(root, model, length, range(length), result)
        _spingc(root, 4, result)
        _kondo(root, result)
    else:
        _canonical(root, 'Spin', 4, range(4), result)
        _canonical(root, 'Spin', 8, [0], result)
        _canonical(root, 'SpinlessFermion', 8, [1], result)
        _kondo_empty_rank(root, result)
    for case in result:
        np.testing.assert_allclose(case.hq, case.hq.conj().T, rtol=0, atol=1e-12)
        assert case.basis.shape == (len(case.states), len(case.hq)), case.name
    return result


def cases(root: Path, selection: str = 'all') -> list[FDCase]:
    # selection is all or boundary; the generators are fixed above.
    return build_cases(root, selection)


def _decoded(output):
    if output is None:
        return ''
    return output if isinstance(output, str) else output.decode(errors='replace')


def execute_case(case, executable, solver, flags, launcher, label, expected_error,
                 timeout=120, env_extra=None):
    values = dict(line.split() for line in (case.path/'calc.def').read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith('#'))
    values.update(CalcType=2, Solver=solver, NGPU=0, OutputEigenVec=0,
                  InputEigenVec=0, InputHam=0, OutputHam=0)
    values.update(flags)
    (case.path/'calc.def').write_text(''.join(f'{key} {value}\n' for key, value in values.items()))
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated', **(env_extra or {}))
    launcher = list(launcher)
    # Same established rule as symmetry_spingc_common.run_case: a single-rank
    # launcher can split an MPI_Abort diagnostic between the two streams.
    if expected_error is not None and spingc.mpi_size(launcher) == 1:
        launcher = []
    command = [*launcher, str(Path(executable).resolve()), '-e', 'sym.def']
    log = case.path/f'{label}.log'
    try:
        result = subprocess.run(command, cwd=case.path, env=env, capture_output=True,
                                text=True, timeout=timeout)
    except subprocess.TimeoutExpired as error:
        log.write_text(_decoded(error.stdout) + _decoded(error.stderr))
        print(log.read_text(), flush=True)
        raise
    text = result.stdout + result.stderr
    log.write_text(text)
    if expected_error is None:
        assert result.returncode == 0, f'{command}\n{text}'
    else:
        assert result.returncode != 0 and expected_error in text, f'{command}\n{text}'
    return text


def run_fd(case: FDCase, executable: Path, solver: int,
           flags: dict[str, int], launcher: list[str], label: str,
           expected_error: str | None = None, timeout: int = 120,
           env_extra: dict[str, str] | None = None) -> str:
    return execute_case(case, executable, solver, flags, launcher, label,
                        expected_error, timeout, env_extra)


def probe_executable(hphi: Path) -> Path:
    return Path(hphi).resolve().with_name('unittest_symmetry_sector_probe')


def read_key_values(path: Path) -> dict[str, str]:
    return dict(line.split('=', 1) for line in Path(path).read_text().splitlines() if line)


def run_probe(case: FDCase, executable: Path, solver: int, launcher: list[str],
              action: str, label: str, expected_error: str | None = None,
              clean: bool = True, flags: dict[str, int] | None = None) -> str:
    """Run a test-probe action; clean drops the previous output and dumps."""
    if clean:
        for stale in case.path.glob('matrix-contract*.dat'):
            stale.unlink()
        for stale in case.path.glob('output/contract_rank_*'):
            stale.unlink()
        if (case.path/'output').exists():
            shutil.rmtree(case.path/'output')
    return run_fd(case, executable, solver, flags or {}, launcher, label, expected_error,
                  env_extra={'HPHI_TEST_SYMMETRY_ACTION': action})


def solver_backends(launcher: list[str], ranks: int) -> list[tuple[int, str, list[str]]]:
    """LAPACK runs serially (HPhi's one-process contract); distributed
    backends need the requested multi-rank launcher and a configured build."""
    backends = [(0, 'LAPACK', [])]
    if ranks > 1 and os.environ.get('HPHI_HAS_SCALAPACK') == '1':
        backends.append((1, 'ScaLAPACK', launcher))
    if ranks > 1 and os.environ.get('HPHI_HAS_ELPA') == '1':
        backends.append((3, 'ELPA', launcher))
    return backends


def matrix_digest(h):
    """hphi-sector-matrix-fnv1a64-v1 of the lower triangle, column-major."""
    import struct
    value = 14695981039346656037

    def feed(word):
        nonlocal value
        for byte in struct.pack('<Q', word):
            value = ((value ^ byte) * 1099511628211) & ((1 << 64)-1)
    feed(1); feed(len(h))
    for col in range(len(h)):
        for row in range(col, len(h)):
            for part in (h[row, col].real, h[row, col].imag):
                bits = struct.unpack('<Q', struct.pack('<d', float(part) if part else 0.))[0]
                feed(bits)
    return value


def basis_order(case: FDCase, group_order: int) -> int:
    """Independent FNV-1a stream of index, representative, orbit, stabilizer."""
    import struct
    value = 14695981039346656037
    for index, column in enumerate(case.basis.T, start=1):
        support = np.flatnonzero(abs(column) > 1e-12)
        orbit = len(support)
        assert group_order % orbit == 0, (case.name, orbit, group_order)
        for word in (index, int(case.states[support[0]]), orbit, group_order // orbit):
            for byte in struct.pack('<Q', word):
                value = ((value ^ byte) * 1099511628211) & ((1 << 64)-1)
    return value


def checkpoint_part(file):
    data = Path(file).read_bytes()
    prefix = np.frombuffer(data[:16], dtype='<u8')
    assert len(prefix) == 2 and int(prefix[W.VERSION]) in (1, 2), file
    assert int(prefix[W.MAGIC]) == CHECKPOINT_MAGIC, file
    words = 28 if int(prefix[W.VERSION]) == 1 else 30
    header = np.frombuffer(data[:8*words], dtype='<u8').copy()
    assert len(header) == words, file
    size = int(header[W.LOCAL_DIM])
    assert len(data) == 8*words + 16*size, file
    vector = np.frombuffer(data[8*words:], dtype='<c16').copy()
    assert np.isfinite(vector).all(), file
    return header, vector


def read_fd_state(path: Path, state: int, ranks: int,
                  expect_state: int | None = -1) -> tuple[np.ndarray, list[np.ndarray]]:
    """Join the rank files <path>/zvo_eigenvec_<state>_rank_<r>.dat.

    expect_state=-1 requires header STATE == state; None skips the check
    (TE names files by step but records the source state)."""
    files = {int(file.name.rsplit('_rank_', 1)[1][:-4]): file
             for file in Path(path).glob(f'zvo_eigenvec_{state}_rank_*.dat')}
    assert sorted(files) == list(range(ranks)), (path, state, sorted(files))
    parts = [checkpoint_part(files[rank]) for rank in range(ranks)]
    headers = [header for header, _ in parts]
    dim = int(headers[0][W.DIM])
    vector = np.empty(dim, complex)
    offset = 0
    for rank, (header, part) in enumerate(parts):
        assert int(header[W.RANKS]) == ranks and int(header[W.RANK]) == rank
        assert int(header[W.DIM]) == dim and int(header[W.LAYOUT]) == 0
        start, size = block_range(dim, ranks, rank)
        assert int(header[W.OFFSET]) == offset == start, (rank, header[W.OFFSET], offset)
        assert int(header[W.LOCAL_DIM]) == len(part) == size
        for word in (W.METHOD, W.STATE, W.STEP, W.TIME, W.HAMILTONIAN,
                     W.DATA_XOR, W.DATA_SUM, W.GROUP, W.SECTOR_COUNT,
                     W.SECTOR_XOR, W.SECTOR_SUM, W.VERSION, W.PHASE):
            assert int(header[word]) == int(headers[0][word]), (rank, word.name)
        if expect_state == -1:
            assert int(header[W.STATE]) == state, (rank, header[W.STATE], state)
        elif expect_state is not None:
            assert int(header[W.STATE]) == expect_state, (rank, header[W.STATE])
        vector[offset:offset+size] = part
        offset += size
    assert offset == dim and np.isfinite(vector).all()
    return vector, headers


def check_eigensystem(hq, energies, vectors):
    """Orthonormal eigenvectors of the independent Hq; degenerate clusters
    are compared as subspaces, isolated states up to a global phase."""
    d = len(hq)
    np.testing.assert_allclose(vectors.conj().T @ vectors, np.eye(d),
                               rtol=0, atol=1e-10)
    residuals = np.linalg.norm(hq @ vectors - vectors * energies, axis=0)
    scale = max(1., np.linalg.norm(hq, 'fro'))
    assert np.max(residuals / (scale*np.linalg.norm(vectors, axis=0))) <= 1e-10
    expected, reference = np.linalg.eigh(hq)
    np.testing.assert_allclose(energies, expected, rtol=0, atol=3e-8)
    for start, end in cluster_edges(hq):
        actual = vectors[:, start:end]
        target = reference[:, start:end]
        np.testing.assert_allclose(actual @ actual.conj().T,
                                   target @ target.conj().T, rtol=0, atol=2e-8)
        if end-start == 1:
            v, u = actual[:, 0], target[:, 0]
            overlap = np.vdot(u, v)
            assert abs(overlap) > 0
            assert np.linalg.norm(v*np.conj(overlap)/abs(overlap)-u) <= 2e-8


def cluster_edges(hq):
    """Degenerate clusters: adjacent independent eigenvalue gaps at most
    1e-10*max(1, ||Hq||_2)."""
    expected = np.linalg.eigvalsh(hq)
    gaps = np.diff(expected)
    edges = [0, *(np.flatnonzero(gaps > 1e-10*max(1., np.linalg.norm(hq, 2)))+1), len(hq)]
    return list(zip(edges, edges[1:]))


def read_states(output: Path, dim: int, ranks: int):
    columns, headers = [], []
    for state in range(dim):
        vector, parts = read_fd_state(output, state, ranks)
        columns.append(vector)
        headers.append(parts)
    return np.column_stack(columns), headers


def read_manifest(output: Path) -> dict[str, str]:
    """Last occurrence wins: the FullDiag status is appended on completion."""
    return read_key_values(Path(output)/'zvo_symmetry_sector.dat')


def grid_max(ranks: int) -> int:
    divisors = [d for d in range(1, int(ranks**0.5)+1) if ranks % d == 0]
    return ranks // max(divisors)


def selected_cases() -> str:
    selection = os.environ.get('HPHI_FD_CASES', 'all')
    assert selection in ('all', 'boundary'), selection
    return selection


GREEN_STEMS = {'one': 'cisajs', 'two': 'cisajscktalt', 'three': 'ThreeBody',
               'four': 'FourBody', 'six': 'SixBody', 'nbody': 'NBodyG'}
GREEN_KEYWORDS = {'one': 'OneBodyG', 'two': 'TwoBodyG', 'three': 'ThreeBodyG',
                  'four': 'FourBodyG', 'six': 'SixBodyG', 'nbody': 'NBodyG'}
GREEN_FACTORS = {'one': 1, 'two': 2, 'three': 3, 'four': 4, 'six': 6}


def row_factors(kind, row):
    data, count = (row[1:], row[0]) if kind == 'nbody' else (row, GREEN_FACTORS[kind])
    assert len(data) == 4*count, (kind, row)
    return tuple(tuple(data[4*f:4*f+4]) for f in range(count))


def _with_duplicate(requests):
    """Every family also carries a repeated row."""
    out = {}
    for kind, rows in requests.items():
        rows = [list(row) for row in rows]
        if not any(rows[i] == rows[j] for i in range(len(rows)) for j in range(i)):
            rows.append(list(rows[0]))
        out[kind] = rows
    return out


def _canonical_requests(case):
    saved = general.ROOT
    import symmetry_correlation as correlation  # sets its own fixture ROOT on import
    general.ROOT = saved
    requests = _with_duplicate(correlation.green_requests(case.model, case.length))
    states = [int(s) for s in case.states]
    operators = {kind: [(tuple(row), case.basis.conj().T @ correlation.matrix(
                    case.model, case.length, states, row_factors(kind, row)) @ case.basis)
                        for row in rows] for kind, rows in requests.items()}
    return requests, operators


def _spingc_requests(case):
    import symmetry_spingc_observables as observables
    nsite = case.length
    requests = {}
    for kind, rows in observables.CORRELATION_REQUESTS.items():
        kept = [row for row in rows
                if all(max(f[0], f[2]) < nsite for f in row_factors(kind, row))]
        assert kept, kind
        requests[kind] = kept
    requests = _with_duplicate(requests)
    operators = {kind: [(tuple(row), case.basis.conj().T @
                         observables.correlation_operator(nsite, kind, row) @ case.basis)
                        for row in rows] for kind, rows in requests.items()}
    return requests, operators


def _kondo_requests(case):
    from symmetry_kondo_reference import apply_operator
    if case.kondo is None:
        return None, None  # explicit-mask KondoGC geometry: no Case for the oracle
    if case.kondo.cells == 3:
        import symmetry_kondo_observables as observables
        requests = observables.correlation_requests()
    else:
        # Two-cell block geometry: local 0, 1 and conduction 2, 3.
        hop, flip, up = [2, 0, 3, 0], [0, 0, 0, 1], [1, 0, 1, 0]
        requests = {'one': [hop, flip, hop],
                    'two': [hop + up, flip + [3, 1, 3, 0], [0, 0, 0, 1, 0, 1, 0, 0], hop + up],
                    'nbody': [[3] + hop + up + [2, 1, 2, 1], [3] + hop + up + [2, 1, 2, 1]]}
    requests = _with_duplicate(requests)
    operators = {}
    for kind, rows in requests.items():
        operators[kind] = []
        for row in rows:
            factors = row_factors(kind, row)
            applied = np.column_stack([apply_operator(case.kondo, factors, column)
                                       for column in case.basis.T])
            operators[kind].append((tuple(row), case.basis.conj().T @ applied))
    return requests, operators


def attach_correlations(case: FDCase, aggregate: bool) -> bool:
    """Write independent correlation requests and their sector operators.

    Returns False when the fixture has no independent operator oracle."""
    if case.model in KONDO_MODELS:
        requests, operators = _kondo_requests(case)
    elif case.model == 'SpinGC':
        requests, operators = _spingc_requests(case)
    else:
        requests, operators = _canonical_requests(case)
    if requests is None:
        return False
    namelist = (case.path/'sym.def').read_text()
    for kind, rows in requests.items():
        spingc.definition(case.path, GREEN_STEMS[kind]+'.def', rows)
        line = f'{GREEN_KEYWORDS[kind]} {GREEN_STEMS[kind]}.def\n'
        if line not in namelist:
            namelist += line
    (case.path/'sym.def').write_text(namelist)
    set_calc_keys(case, {'OutputGreenFormat': int(aggregate)})
    case.operators = operators
    return True


def set_calc_keys(case: FDCase, updates: dict):
    values = dict(line.split() for line in (case.path/'calc.def').read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith('#'))
    values.update({key: str(value) for key, value in updates.items()})
    (case.path/'calc.def').write_text(''.join(f'{k} {v}\n' for k, v in values.items()))


def read_green(output: Path, family: str, states: int, aggregate: bool,
               rows: list[tuple[int, ...]]) -> np.ndarray:
    """(states, rows) complex values; index columns must equal the requests."""
    stem = GREEN_STEMS[family]
    expected = np.asarray(rows, dtype=float)
    values = np.empty((states, len(rows)), complex)
    if aggregate:
        data = np.loadtxt(Path(output)/f'zvo_{stem}_eigen.dat', ndmin=2)
        assert data.shape == (states*len(rows), expected.shape[1]+3), (stem, data.shape)
        np.testing.assert_array_equal(data[:, 0], np.repeat(np.arange(states), len(rows)))
        blocks = data[:, 1:].reshape(states, len(rows), -1)
    else:
        blocks = []
        for state in range(states):
            data = np.loadtxt(Path(output)/f'zvo_{stem}_eigen{state}.dat', ndmin=2)
            assert data.shape == (len(rows), expected.shape[1]+2), (stem, state, data.shape)
            blocks.append(data)
        blocks = np.asarray(blocks)
    for state in range(states):
        np.testing.assert_array_equal(blocks[state][:, :-2], expected)
        values[state] = blocks[state][:, -2] + 1j*blocks[state][:, -1]
    return values


def zero_hamiltonian_case(root: Path) -> FDCase:
    """SpinGC L4 k0 without terms: every sector state is degenerate."""
    nsite = 4
    perms = [[(i+g) % nsite for i in range(nsite)] for g in range(nsite)]
    chars = np.ones(nsite)
    basis, _, projector = spingc.sector_basis(nsite, perms, chars)
    check_basis(projector, basis)
    name = 'SpinGC_zero_L4_k0'
    spingc.write_case(root/name, nsite, perms, chars, {}, 2, {})
    dim = basis.shape[1]
    return FDCase(name, 'SpinGC', root/name, np.arange(1 << nsite), basis,
                  np.zeros((dim, dim), complex), {}, nsite, 0)


def set_calc_file(path: Path, updates: dict):
    """Update CalcMod keys; TE/CG consumers always use Solver 0 and NGPU 0."""
    values = dict(line.split() for line in (Path(path)/'calc.def').read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith('#'))
    if int(updates.get('CalcType', values.get('CalcType', 2))) != 2:
        values.update(Solver='0', NGPU='0')
    values.update({key: str(value) for key, value in updates.items()})
    (Path(path)/'calc.def').write_text(''.join(f'{k} {v}\n' for k, v in values.items()))


def set_mod_file(path: Path, updates: dict):
    lines = (Path(path)/'mod.def').read_text().splitlines()
    head, body = lines[:8], lines[8:]
    values = dict(line.split() for line in body if line.strip())
    values.update({key: str(value) for key, value in updates.items()})
    (Path(path)/'mod.def').write_text('\n'.join(head) + '\n'
                                     + ''.join(f'{k} {v}\n' for k, v in values.items()))


def read_green_file(file: Path, rows: list[tuple[int, ...]]) -> np.ndarray:
    data = np.loadtxt(file, ndmin=2)
    expected = np.asarray(rows, dtype=float)
    assert data.shape == (len(rows), expected.shape[1]+2), (file, data.shape)
    np.testing.assert_array_equal(data[:, :-2], expected)
    return data[:, -2] + 1j*data[:, -1]


def uniform_quench(case: FDCase, coefficient: float):
    """Translation-invariant diagonal InterAll term and its independent Oq.

    Spin/SpinGC/tJ: sum Sz_i Sz_i+1; SpinlessFermion: sum n_i n_i+1;
    Hubbard: sum n_i n_i+1 with n = n_up + n_down; Kondo: sum over cells
    Sloc_z Scond_z (spin index 0 is Sz=+1/2). For three tJ electrons on the
    four-site ring, sum n_i n_i+1 is the constant 2, a mere energy shift."""
    rows, terms = [], []
    if case.model in KONDO_MODELS:
        cells = case.kondo.cells
        for j in range(cells):
            l, c = j, cells + j  # block layout of the fixtures
            for s, sz in ((0, .5), (1, -.5)):
                for t, tz in ((0, .5), (1, -.5)):
                    rows.append((l, s, l, s, c, t, c, t, coefficient*sz*tz, 0))
    else:
        L = case.length
        spins = (0,) if case.model == 'SpinlessFermion' else (0, 1)
        for i in range(L):
            j = (i+1) % L
            for s in spins:
                for t in spins:
                    if case.model in ('Spin', 'SpinGC', 'tJ'):
                        value = coefficient*.25*(1 if s == t else -1)
                    else:
                        value = coefficient
                    rows.append((i, s, i, s, j, t, j, t, value, 0))
    if case.model == 'SpinGC':
        raw = spingc.families_matrix(case.length, {'InterAll': rows})
        delta = case.basis.conj().T @ raw @ case.basis
    elif case.model in KONDO_MODELS:
        from symmetry_kondo_reference import apply_operator
        applied = np.zeros(case.basis.shape, complex)
        for row in rows:
            factors = (tuple(row[0:4]), tuple(row[4:8]))
            applied += row[8]*np.column_stack([apply_operator(case.kondo, factors, column)
                                                for column in case.basis.T])
        delta = case.basis.conj().T @ applied
    else:
        saved = general.ROOT
        import symmetry_correlation as correlation
        general.ROOT = saved
        states = [int(s) for s in case.states]
        raw = sum(correlation.matrix(case.model, case.length, states,
                                     (tuple(row[0:4]), tuple(row[4:8])), row[8]) for row in rows)
        delta = case.basis.conj().T @ raw @ case.basis
    np.testing.assert_allclose(delta, delta.conj().T, rtol=0, atol=1e-12)
    return rows, delta
