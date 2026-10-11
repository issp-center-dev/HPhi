"""Sector FullDiag eigenstate checkpoints: provenance contract and all-state output."""
from pathlib import Path
import shutil
import struct
import sys

import numpy as np

import symmetry_fulldiag_common as common
from symmetry_fulldiag_common import W

FULLDIAG = 2
KIND_LINE = 'source_hamiltonian_kind=sector_matrix_v1 hamiltonian_changed=unknown'
INTACT = 'FullDiagProbe sentinel intact after rejected checkpoint'
HEADER_ERROR = 'checkpoint header / sector / layout validation failed'


def state_files(output: Path, state: int = 0) -> list[Path]:
    return sorted(output.glob(f'zvo_eigenvec_{state}_rank_*.dat'),
                  key=lambda file: int(file.stem.rsplit('_', 1)[1]))


def write_checkpoint(case, probe, solver, launcher, label, action='fulldiag-checkpoint'):
    text = common.run_probe(case, probe, solver, launcher, action, label)
    assert 'FullDiagProbe checkpoint round-trip ok' in text, text
    return {key: int(value) for key, value in
            common.read_key_values(case.path/'matrix-contract.dat').items()}


def check_headers(case, contract, ranks):
    output = case.path/'output'
    assert len(list(output.glob('contract_rank_*.dat'))) == ranks
    vector, headers = common.read_fd_state(output, 0, ranks)
    expected = np.zeros(len(case.hq), complex)
    expected[0] = 1
    np.testing.assert_array_equal(vector, expected)
    assert all(int(h[W.METHOD]) == FULLDIAG and int(h[W.STATE]) == 0 for h in headers)
    assert all(int(h[W.STEP]) == 0 and int(h[W.TIME]) == 0 for h in headers)
    assert all(int(h[W.DIM]) == len(case.hq) and int(h[W.RANKS]) == ranks for h in headers)
    assert all(int(h[W.HAMILTONIAN]) == int(headers[0][W.HAMILTONIAN]) for h in headers)
    assert int(headers[0][W.HAMILTONIAN]) == contract['matrix_digest']
    kondo = case.model.startswith('Kondo')
    header = headers[0]
    assert int(header[W.VERSION]) == (2 if kondo else 1) == int(header[W.PHASE])
    # The matrix identity is the checkpoint normalization word for word.
    for word, key in ((W.MODEL, 'model'), (W.NSITE, 'nsite'), (W.NUP, 'nup'),
                      (W.NDOWN, 'ndown'), (W.NE, 'ne'), (W.FULL_DIM, 'full_dim'),
                      (W.DIM, 'sector_dim'), (W.GROUP, 'group_digest'),
                      (W.SECTOR_COUNT, 'sector_count'), (W.SECTOR_XOR, 'sector_xor'),
                      (W.SECTOR_SUM, 'sector_sum'), (W.PHASE, 'phase')):
        assert int(header[word]) == contract[key], (case.name, word.name)
    if kondo:
        assert int(header[W.LOCAL_SITE_MASK]) == contract['local_site_mask']
        assert int(header[W.FIXED_FLAGS]) == contract['fixed_flags']
    if ranks == 1:
        # One rank: the local ORDER stream is the global basis_order stream.
        assert int(header[W.ORDER]) == contract['basis_order']
    return headers


def flip_word(file: Path, word: int, value: int | None = None):
    data = bytearray(file.read_bytes())
    old = struct.unpack_from('<Q', data, 8*word)[0]
    struct.pack_into('<Q', data, 8*word, old ^ 1 if value is None else value)
    file.write_bytes(bytes(data))


def reject(case, probe, solver, launcher, label, message):
    text = common.run_probe(case, probe, solver, launcher, 'fulldiag-checkpoint-read',
                            label, message, clean=False)
    assert INTACT in text, text


def restore(files, saved):
    for file in files:
        file.write_bytes(saved[file])


def negative_cases(by_name, probe, solver, launcher, label, ranks):
    """Rejected FullDiag checkpoints leave the destination vector untouched."""
    source, other, kondo = by_name['Spin_L6_k1'], by_name['Spin_L6_k2'], by_name['Kondo_P3_k0']
    assert len(source.hq) == len(other.hq), 'same-dimension sectors are required'
    write_checkpoint(source, probe, solver, launcher, f'negative_source_{label}')
    text = common.run_probe(source, probe, solver, launcher, 'fulldiag-checkpoint-read',
                            f'read_ok_{label}', clean=False)
    assert 'FullDiagProbe checkpoint read ok state=0 kind=1' in text and KIND_LINE in text, text
    files = state_files(source.path/'output')
    saved = {file: file.read_bytes() for file in files}
    # Same dimension, different momentum sector.
    write_checkpoint(other, probe, solver, launcher, f'negative_other_{label}')
    for file in state_files(other.path/'output'):
        file.unlink()
    for file in files:
        shutil.copy2(file, other.path/'output'/file.name)
    reject(other, probe, solver, launcher, f'reject_other_sector_{label}', HEADER_ERROR)
    victim = files[-1]
    for word in (W.PHASE, W.ORDER, W.LAYOUT, W.RANKS, W.RANK, W.OFFSET):
        flip_word(victim, word)
        reject(source, probe, solver, launcher, f'reject_{word.name.lower()}_{label}', HEADER_ERROR)
        restore(files, saved)
    # FullDiag provenance: step/time must be zero and the state inside D.
    for word, value in ((W.STEP, 1), (W.TIME, struct.unpack('<Q', struct.pack('<d', 1.))[0]),
                        (W.STATE, len(source.hq))):
        for file in files:
            flip_word(file, word, value)
        reject(source, probe, solver, launcher, f'reject_fulldiag_{word.name.lower()}_{label}',
               HEADER_ERROR)
        restore(files, saved)
    owned = next(file for file in files if common.checkpoint_part(file)[0][W.LOCAL_DIM] > 0)
    for name, payload, message in (
            ('truncated', saved[owned][:-1], 'checkpoint payload length failed'),
            ('trailing', saved[owned] + b'x', 'checkpoint payload length failed'),
            ('nan', None, 'checkpoint finite vector validation failed')):
        if payload is None:
            data = bytearray(saved[owned])
            words = 28 if common.checkpoint_part(owned)[0][W.VERSION] == 1 else 30
            struct.pack_into('<d', data, 8*words, float('nan'))
            payload = bytes(data)
        owned.write_bytes(payload)
        reject(source, probe, solver, launcher, f'reject_{name}_{label}', message)
        restore(files, saved)
    if ranks > 1:
        # Rank files of another state or another H must not be combined.
        for word in (W.STATE, W.HAMILTONIAN):
            flip_word(victim, word)
            reject(source, probe, solver, launcher, f'reject_mixed_{word.name.lower()}_{label}',
                   'checkpoint cross-rank metadata consistency failed')
            restore(files, saved)
        write_checkpoint(source, probe, solver, launcher, f'negative_known_{label}',
                         'fulldiag-checkpoint-known')
        known = {file: file.read_bytes() for file in state_files(source.path/'output')}
        restore(files, saved)
        victim.write_bytes(known[victim])
        reject(source, probe, solver, launcher, f'reject_mixed_matrix_{label}',
               'checkpoint cross-rank metadata consistency failed')
        restore(files, saved)
    write_checkpoint(kondo, probe, solver, launcher, f'negative_kondo_{label}')
    kondo_files = state_files(kondo.path/'output')
    kondo_saved = {file: file.read_bytes() for file in kondo_files}
    for word in (W.LOCAL_SITE_MASK, W.FIXED_FLAGS):
        flip_word(kondo_files[-1], word)
        reject(kondo, probe, solver, launcher, f'reject_{word.name.lower()}_{label}', HEADER_ERROR)
        restore(kondo_files, kondo_saved)
    print(f'{label}: FullDiag checkpoint negatives PASS (np={ranks})', flush=True)
    return files, saved


def checkpoint_mode(hphi, probe, launcher, ranks):
    root = common.suite_root('states-checkpoint')
    backends = common.solver_backends(launcher, ranks)
    cases = common.cases(root, 'all')
    by_name = {case.name: case for case in cases}
    for case in cases:
        for solver, label, backend_launcher in backends:
            backend_ranks = common.spingc.mpi_size(backend_launcher)
            contract = write_checkpoint(case, probe, solver, backend_launcher, f'checkpoint_{label}')
            check_headers(case, contract, backend_ranks)
            text = common.run_probe(case, probe, solver, backend_launcher,
                                    'fulldiag-checkpoint-read', f'read_{label}', clean=False)
            assert KIND_LINE in text, text
            print(f'{case.name} {label} np={backend_ranks}: FullDiag checkpoint provenance PASS',
                  flush=True)
    for solver, label, backend_launcher in backends:
        backend_ranks = common.spingc.mpi_size(backend_launcher)
        files, saved = negative_cases(by_name, probe, solver, backend_launcher, label,
                                      backend_ranks)
        if backend_ranks > 1:
            # Files written by more ranks are rejected by a one-rank reader.
            source = by_name['Spin_L6_k1']
            restore(files, saved)
            reject(source, probe, 0, [], f'reject_rank_count_{label}', HEADER_ERROR)


def main():
    hphi = Path(sys.argv[1]).resolve()
    mode = sys.argv[2] if len(sys.argv) > 2 else 'states'
    launcher, ranks = common.launcher_from_env()
    if mode == 'checkpoint':
        checkpoint_mode(hphi, common.probe_executable(hphi), launcher, ranks)
    else:
        raise ValueError(mode)
    print(f'sector FullDiag states ({mode}) PASS', flush=True)


if __name__ == '__main__':
    main()
