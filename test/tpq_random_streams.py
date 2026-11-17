"""Check production initial states for reused streams across samples/workers.

No stochastic acceptance criterion: the old bug repeats whole component blocks.
The helper calls MakeIniVec directly; sector solver tests separately compare
subsequent propagation against independently constructed dense Hamiltonians.
"""
import itertools
import os
from pathlib import Path
import shlex
import struct
import subprocess
import sys
import tempfile

import numpy as np

helper = str(Path(sys.argv[1]).resolve())
mpi = shlex.split(os.environ.get('MPIRUN', '')) if '--mpi' in sys.argv else []
ranks = 1
for i, arg in enumerate(mpi[:-1]):
    if arg in ('-np', '-n'):
        ranks = int(mpi[i+1])
if '--mpi' in sys.argv and ranks < 2:
    print('Skipping: a counted MPIRUN with at least two ranks is required')
    sys.exit(77)
root = Path(tempfile.mkdtemp(prefix='tpq-random-streams-', dir='.')).resolve()
long_bits = 8*struct.calcsize('l')


def initial(kind, seed, threads, samples=3, suffix=''):
    path = root / 'type{}_seed{}_threads{}_samples{}{}'.format(kind, seed, threads, samples, suffix)
    path.mkdir()
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='FALSE')
    with (path / 'run.log').open('w') as out:
        result = subprocess.run(mpi + [helper, '--dump-initial', str(32*ranks), str(kind), str(seed), str(samples)],
                                cwd=path, env=env, stdout=out, stderr=subprocess.STDOUT, timeout=90)
    assert result.returncode == 0, (path / 'run.log').read_text()
    states, chunks = [], []
    for sample in range(samples):
        parts = []
        for rank in range(ranks):
            lines = (path / 'initial_{}_rank{}.dat'.format(sample, rank)).read_text().splitlines()
            norm = float(lines[0])
            values = np.array([complex(*map(float, line.split())) for line in lines[1:]])
            assert len(values) == 32 and np.isfinite(values).all() and np.isfinite(norm) and norm > 0
            parts.append(values)
            # Static OpenMP ownership; all workers have equal-sized nonempty blocks.
            chunks.extend((sample, rank, thread, block*norm)
                          for thread, block in enumerate(np.split(values, threads)))
        state = np.concatenate(parts)
        np.testing.assert_allclose(np.vdot(state, state), 1, rtol=0, atol=1e-12)
        states.append(state)
    for left, right in itertools.combinations(chunks, 2):
        assert not np.allclose(left[3], right[3], rtol=0, atol=1e-12), (
            'reused random stream', kind, seed, threads, left[:3], right[:3])
    return np.array(states)


for kind in (0, 1, -1):
    for threads in (2, 1):
        positive = initial(kind, 1, threads)
        # MPI coverage isolates rank/thread separation for each distribution.
        # Seed limits/signs and NumAve-prefix cases below need no extra ranks.
        if mpi:
            continue
        negative = initial(kind, -1, threads)
        np.testing.assert_allclose(positive, negative, rtol=0, atol=1e-14)
        repeated = initial(kind, 1, threads, samples=4)
        np.testing.assert_allclose(positive, repeated[:3], rtol=0, atol=1e-14)
        initial(kind, 0, threads)
        initial(kind, -(1 << (long_bits-1)), threads)
        if long_bits > 32:
            high = initial(kind, (1 << 32)+1, threads)
            assert not np.allclose(positive, high, rtol=0, atol=1e-12), 'seed high bits were discarded'
print('TPQ RNG: no repeated sample/rank/thread blocks; deterministic prefixes, seed signs, zero and limits PASS ({} ranks)'.format(ranks))
