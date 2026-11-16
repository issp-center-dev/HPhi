"""Check canonical Kondo empty MPI ranks against exhaustive occupations."""
import argparse
import os
from pathlib import Path
import re
import shlex
import tempfile

import kondo_basis_bounds as kondo


def check_dimensions(path, case, ranks):
    words = kondo.basis(*case)
    inter_sites, remaining = 0, ranks
    while remaining > 1:
        assert remaining % 4 == 0, ranks
        inter_sites += 1
        remaining //= 4
    shift = 2*(case[0]-inter_sites)
    expected = [sum(w >> shift == rank for w in words) for rank in range(ranks)]
    log = (path/'cg.log').read_text()
    table = log.split('Process       Dimension', 1)[1].split('Total dimension', 1)[0]
    rows = re.findall(r'^\s*(\d+)\s+(\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+',
                      table, re.MULTILINE)
    assert [int(row[0]) for row in rows] == list(range(ranks)), rows
    actual = [int(row[1]) for row in rows]
    assert actual == expected, (path, actual, expected)
    assert sum(actual) == len(words), (path, actual, len(words))
    print(f'{path.name}: rank dimensions={actual}, empty={actual.count(0)} passed', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('hphi', type=lambda p: str(Path(p).resolve()))
    args = parser.parse_args()
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = 1
    if launcher:
        positions = [i for i, value in enumerate(launcher) if value in ('-n', '-np')]
        assert len(positions) == 1, launcher
        ranks = int(launcher[positions[0]+1])
    root = Path(tempfile.mkdtemp(prefix='kondo_mpi_dimension_', dir='.')).resolve()
    print(f'Kondo dimension artifacts: {root}', flush=True)
    env = dict(os.environ)
    env.setdefault('OMP_NUM_THREADS', '1')
    env['OPENBLAS_NUM_THREADS'] = env['VECLIB_MAXIMUM_THREADS'] = '1'
    # With np16, low_balanced has wrapped Nup on rank 5 and wrapped Ndown
    # on rank 10. The polarized pair also covers each spin at np4.
    # high_balanced tests Nup/Ndown > local Nsite without unsigned wrap.
    cases = [('low_balanced', (6, (0,), 1, 0), 10),
             ('low_up', (6, (0,), 1, 2), 5),
             ('low_down', (6, (0,), 1, -2), 5),
             ('high_balanced', (6, (0,), 9, 0), 10)]
    for name, case, dimension in cases:
        assert len(kondo.basis(*case)) == dimension
        for hs in (0, 1):
            path = root/f'{name}_hs{hs}'
            kondo.check(path, args.hphi, case, hs, launcher, ranks, env,
                        dense=False, timeout=30)
            check_dimensions(path, case, ranks)


if __name__ == '__main__':
    main()
