def read_gs(*, index=0, rank=0):
    import numpy as np
    from os.path import join
    from struct import unpack

    filename = join("output",
                    "zvo_eigenvec_{}_rank_{}.dat".format(index,
                                                         rank))
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

vec = read_gs()
#print(vec)
