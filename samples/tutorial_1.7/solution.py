import subprocess

import numpy as np
from numpy.linalg import norm

import read_gs


def create_input():
    with open("std.in", "w") as f:
        print('lattice = "chain"', file=f)
        print('L = 8', file=f)
        print('model = "spin"', file=f)
        print('2S = 1', file=f)
        print('2Sz = 0', file=f)
        print('J = 1', file=f)
        print('method = "cg"', file=f)
        print('exct = 2', file=f)
        print('eigenvecio = "out"', file=f)


if __name__ == "__main__":
    import sys
    HPhi = sys.argv[1] if len(sys.argv) > 1 else "HPhi"
    create_input()
    cmd = [HPhi, "-s", "std.in"]
    subprocess.call(cmd)
    print()

    gs = read_gs.read_gs()
    exc = read_gs.read_gs(index=1)

    print("norm of gs    = ", norm(gs))
    print("norm of exc   = ", norm(exc))
    print("inner product = ", np.abs(np.dot(gs.conj(), exc)))
