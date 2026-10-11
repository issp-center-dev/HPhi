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
                 timeout=120):
    values = dict(line.split() for line in (case.path/'calc.def').read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith('#'))
    values.update(CalcType=2, Solver=solver, NGPU=0, OutputEigenVec=0,
                  InputEigenVec=0, InputHam=0, OutputHam=0)
    values.update(flags)
    (case.path/'calc.def').write_text(''.join(f'{key} {value}\n' for key, value in values.items()))
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT='replicated')
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
           expected_error: str | None = None, timeout: int = 120) -> str:
    return execute_case(case, executable, solver, flags, launcher, label,
                        expected_error, timeout)


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
