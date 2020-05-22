import subprocess

import numpy as np
from numpy.linalg import norm


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


def read_gs(*, exct=0, rank=0):
    from os.path import join
    from struct import unpack

    filename = join("output", "zvo_eigenvec_{}_rank_{}.dat".format(exct, rank))
    with open(filename, "rb") as f:
        f.read(4)
        nelems = unpack("L", f.read(8))[0]
        ret = np.zeros(nelems, dtype=np.complex)
        f.read(16)
        for i in range(nelems):
            re = unpack("d", f.read(8))[0]
            im = unpack("d", f.read(8))[0]
            ret[i] = np.complex(re, im)
        return ret


if __name__ == "__main__":
    create_input()
    cmd = ["HPhi", "-s", "std.in"]
    subprocess.call(cmd)
    print()

    gs = read_gs()
    exc = read_gs(exct=1)

    print("norm of gs    = ", norm(gs))
    print("norm of exc   = ", norm(exc))
    print("inner product = ", np.abs(np.dot(gs.conj(), exc)))
