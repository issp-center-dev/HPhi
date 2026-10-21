"""Independent small Kondo Hamiltonians for CalcHS=0/1 and raw site MPI."""
import argparse
from itertools import product
import os
from pathlib import Path
import shlex
import struct
import subprocess
import tempfile

import numpy as np


def basis(nsite, local, ncond, sz2):
    words = []
    for digits in product(range(4), repeat=nsite):
        if any(digits[i] not in (1, 2) for i in local):
            continue
        if sum(d.bit_count() for i, d in enumerate(digits) if i not in local) != ncond:
            continue
        if sum((d & 1) - ((d >> 1) & 1) for d in digits) != sz2:
            continue
        words.append(sum(d << (2*i) for i, d in enumerate(digits)))
    return sorted(words)


def terms(nsite, local):
    """Itinerant hopping, diagonal fields and local/conduction spin exchange."""
    transfer, interall = [], []
    conduction = [i for i in range(nsite) if i not in local]
    for i in range(nsite):
        for spin in (0, 1):
            potential = 0.037*(i+1) + (1-2*spin)*0.023*(i+1)
            transfer.append((i, spin, i, spin, -potential, 0.0))
    for i, j in zip(conduction, conduction[1:]):
        for spin in (0, 1):
            transfer += [(i, spin, j, spin, 0.7, 0.13),
                         (j, spin, i, spin, 0.7, -0.13)]
    for k, i in enumerate(local):
        j = conduction[k % len(conduction)]
        coupling = 0.81 + 0.17*k
        for s, t in product((0, 1), repeat=2):
            interall.append((i, s, i, s, j, t, j, t,
                             coupling*(1-2*s)*(1-2*t)/4, 0.0))
        interall += [(i, 0, i, 1, j, 1, j, 0, coupling/2, 0.0),
                     (j, 0, j, 1, i, 1, i, 0, coupling/2, 0.0)]
    return transfer, interall


def hamiltonian(words, transfer, interall):
    # Apply the literal input operator strings before projecting. Do not use
    # the parser's canonicalization or any production basis/counting helper.
    lookup = {w: i for i, w in enumerate(words)}
    h = np.zeros((len(words), len(words)), dtype=complex)
    for rows, sign in ((transfer, -1), (interall, 1)):
        for row in rows:
            indices = row[:-2]
            operators = [(j % 2 == 0, 2*indices[2*j] + indices[2*j+1])
                         for j in range(len(indices)//2)]
            for col, state in enumerate(words):
                value = sign*complex(*row[-2:])
                for create, orbital in reversed(operators):
                    if bool(state & (1 << orbital)) == create:
                        break
                    value *= (-1)**((state & ((1 << orbital)-1)).bit_count())
                    state ^= 1 << orbital
                else:
                    if state in lookup:
                        h[lookup[state], col] += value
    np.testing.assert_allclose(h, h.conj().T, atol=1e-13, rtol=0)
    return h


def definition(path, label, rows, count=None):
    path.write_text(f"====\n{label} {len(rows) if count is None else count}\n====\n====\n====\n"
                    + ''.join(' '.join(map(str, row))+'\n' for row in rows))


def prepare(path, nsite, local, ncond, sz2, hs):
    path.mkdir()
    transfer, interall = terms(nsite, local)
    definition(path/'loc.def', 'NlocalSpin',
               [(i, int(i in local)) for i in range(nsite)], len(local))
    definition(path/'trans.def', 'NTransfer', transfer)
    definition(path/'interall.def', 'NInterAll', interall)
    (path/'namelist.def').write_text('CalcMod calc.def\nModPara mod.def\n'
                                   'LocSpin loc.def\nTrans trans.def\nInterAll interall.def\n')
    (path/'mod.def').write_text(
        '====\nModel_Parameters 0\n====\n====\n====\n'
        'CDataFileHead zvo\nCParaFileHead zqp\n====\n'
        f'Nsite {nsite}\nNcond {ncond}\n2Sz {sz2}\n'
        + ('' if hs is None else f'CalcHS {hs}\n')
        + 'initial_iv -1\nexct 1\nLanczosEps 22\nLanczos_max 500\n'
          'LanczosTarget 1\nLargeValue 100\nPreCG 0\n')
    return transfer, interall


def run(path, hphi, method, launcher, env, timeout=180):
    (path/'calc.def').write_text(f'CalcModel 2\nCalcType {method}\nOutputMode 0\n'
                               + ('OutputEigenVec 1\n' if method == 3 else 'Solver 0\n'))
    logpath = path/('cg.log' if method == 3 else 'fulldiag.log')
    try:
        with logpath.open('w') as log:
            result = subprocess.run(launcher + [hphi, '-e', 'namelist.def'],
                                    cwd=path, env=env, stdout=log,
                                    stderr=subprocess.STDOUT, timeout=timeout)
        if result.returncode:
            raise RuntimeError(f'{path.name}: method={method} exit={result.returncode}')
    except (RuntimeError, subprocess.TimeoutExpired):
        print(logpath.read_text(), flush=True)
        raise


def vector(path, words, nsite, ranks):
    files = sorted((path/'output').glob('zvo_eigenvec_0_rank_*.dat'),
                   key=lambda p: int(p.stem.rsplit('_', 1)[1]))
    assert [int(p.stem.rsplit('_', 1)[1]) for p in files] == list(range(ranks))
    inter_sites, remaining = 0, ranks
    while remaining > 1:
        assert remaining % 4 == 0
        inter_sites += 1
        remaining //= 4
    shift = 2*(nsite-inter_sites)
    psi = np.zeros(len(words), dtype=complex)
    for rank, file in enumerate(files):
        data = file.read_bytes()
        _, size = struct.unpack_from('=iQ', data)
        part = np.frombuffer(data[12:], dtype=np.complex128)
        indices = [i for i, w in enumerate(words) if w >> shift == rank]
        assert size == len(indices) and len(part) == size+1, (file, size, len(indices))
        psi[indices] = part[1:]  # HPhi raw vectors use a dummy zeroth element.
    assert np.isfinite(psi).all()
    np.testing.assert_allclose(np.linalg.norm(psi), 1, atol=1e-12, rtol=0)
    return psi


def check_result(path, hphi, case, hs, launcher, ranks, env, dense, timeout=180):
    nsite, local, ncond, sz2 = case
    words = basis(*case)
    transfer, interall = prepare(path, *case, hs)
    h = hamiltonian(words, transfer, interall)
    expected = np.linalg.eigvalsh(h)
    tolerance = 1e-10*max(1, np.max(abs(expected)))
    if dense:
        run(path, hphi, 2, [], env, timeout)
        actual = np.loadtxt(path/'output/Eigenvalue.dat', ndmin=2)[:, 1]
        assert actual.shape == expected.shape, (path, actual.shape, expected.shape)
        np.testing.assert_allclose(actual, expected, atol=tolerance, rtol=0)
    run(path, hphi, 3, launcher, env, timeout)
    psi = vector(path, words, nsite, ranks)
    energy = np.vdot(psi, h @ psi).real
    residual = np.linalg.norm(h @ psi-energy*psi)/max(1, np.max(abs(expected)))
    assert abs(energy-expected[0]) <= tolerance and residual <= 1e-10, (path, energy, expected[0], residual)
    reported = [float(line.split()[1]) for line in (path/'output/zvo_energy.dat').read_text().splitlines()
                if line.split() and line.split()[0] == 'Energy']
    assert len(reported) == 1 and abs(reported[0]-energy) <= tolerance, (path, reported)
    print(f'{path.name}: D={len(words)} ranks={ranks} energy={energy:.14g} residual={residual:.3e} passed', flush=True)


def check(path, *args, **kwargs):
    try:
        check_result(path, *args, **kwargs)
    except Exception:
        print(f'Failed case: {path}', flush=True)
        for name in ('calc.def', 'mod.def', 'loc.def', 'trans.def', 'interall.def',
                     'fulldiag.log', 'cg.log'):
            file = path/name
            if file.exists():
                print(f'--- {name} ---\n{file.read_text()}', flush=True)
        raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('hphi', type=lambda p: str(Path(p).resolve()))
    parser.add_argument('mode', choices=('serial', 'mpi', 'reproduce'))
    args = parser.parse_args()
    root = Path(tempfile.mkdtemp(prefix=f'kondo_basis_{args.mode}_', dir='.')).resolve()
    print(f'Kondo basis artifacts: {root}', flush=True)
    env = dict(os.environ)
    env.setdefault('OMP_NUM_THREADS', '1')
    env['OPENBLAS_NUM_THREADS'] = env['VECLIB_MAXIMUM_THREADS'] = '1'
    if args.mode == 'reproduce':
        path = root/'low_filling'
        prepare(path, 4, (0,), 2, 1, None)
        run(path, args.hphi, 2, [], env)
        return
    if args.mode == 'serial':
        cases = [(4, (0,), 2, 1), (4, (3,), 2, 1), (4, (3,), 2, -1),
                 (5, (2,), 2, 1), (5, (0,), 2, 1)]
        for i, case in enumerate(cases):
            for hs in (0, None):  # Omitted CalcHS exercises the default fast path.
                name = f'case{i}_hs{0 if hs == 0 else "default"}'
                check(root/name, args.hphi, case, hs, [], 1, env, dense=True)
    else:
        launcher = shlex.split(os.environ.get('MPIRUN', ''))
        ranks = 1
        if launcher:
            positions = [i for i, value in enumerate(launcher) if value in ('-n', '-np')]
            assert len(positions) == 1, launcher
            ranks = int(launcher[positions[0]+1])
        case = (6, (0, 1), 2, 0)
        for hs in (0, 1):
            check(root/f'mpi_hs{hs}', args.hphi, case, hs, launcher, ranks, env, dense=False)


if __name__ == '__main__':
    main()
