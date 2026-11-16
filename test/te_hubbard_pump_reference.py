"""Reproduce the AC/DC/pulse Hubbard-square reference tables independently.

After running the corresponding serial CTest, run this script with its case
directory as argument. It constructs fermion matrix elements in NumPy, checks
SS/Norm/Flct against HPhi, and writes reference-current-h.dat (Flct columns).
No HPhi matvec or output observables are used to generate the reference.
"""
from pathlib import Path
import struct
import sys

import numpy as np


def rows(path):
    return [line.split() for line in path.read_text().splitlines()[5:] if line.strip()]


def main(path):
    names = dict(line.split() for line in (path / "namelist.def").read_text().splitlines())
    # Fail rather than silently omit any unimplemented Hamiltonian term.
    allowed = {"ModPara", "LocSpin", "Trans", "CoulombIntra", "OneBodyG",
               "TwoBodyG", "CalcMod", "PairExcitation", "TEOneBody", "SpectrumVec"}
    assert set(names) <= allowed, set(names)-allowed
    params = dict(line.split()[:2] for line in (path / names["ModPara"]).read_text().splitlines()
                  if len(line.split()) >= 2)
    nsite, number, sz, order = [int(params[k]) for k in ("Nsite", "Ncond", "2Sz", "ExpandCoef")]
    assert int(params["Lanczos_max"]) == 40
    assert all(int(r[1]) == 0 for r in rows(path / names["LocSpin"]))
    all_words = np.arange(1 << (2*nsite))
    up = sum((all_words >> (2*i)) & 1 for i in range(nsite))
    down = sum((all_words >> (2*i+1)) & 1 for i in range(nsite))
    words = all_words[(up+down == number) & (up-down == sz)]
    dim = len(words)
    lookup = np.full(len(all_words), -1)
    lookup[words] = np.arange(dim)
    doubles = sum(((words >> (2*i)) & 3) == 3 for i in range(nsite))
    diagonal = sum(float(r[1])*(((words >> (2*int(r[0]))) & 3) == 3)
                   for r in rows(path / names["CoulombIntra"]))
    seed = path / "output" / (names["SpectrumVec"] + "_rank_0.dat")
    data = seed.read_bytes()
    assert struct.unpack_from("=iQ", data)[1] == dim, "Run the fixture in serial first"
    psi = np.frombuffer(data, dtype="=c16", offset=12)[1:].copy()
    assert abs(np.linalg.norm(psi)-1) < 1e-10
    static = rows(path / names["Trans"])
    schedule = iter(rows(path / names["TEOneBody"]))
    cache = {}

    def edges(indices):
        key = tuple(indices)
        if key not in cache:
            a, sa, b, sb = map(int, indices)
            target = words.copy()
            valid = np.ones(dim, bool)
            signs = np.ones(dim)
            for create, orbital in ((False, 2*b+sb), (True, 2*a+sa)):
                valid &= ((target >> orbital) & 1) != int(create)
                parity = sum((target >> bit) & 1 for bit in range(orbital))
                signs *= 1-2*(parity % 2)
                target ^= 1 << orbital
            cols = np.flatnonzero(valid)
            dest = lookup[target[valid]]
            assert np.all(dest >= 0)
            cache[key] = dest, cols, signs[valid]
        return cache[key]

    flct, energies, norms = [], [], []
    previous = None
    for step in range(40):
        time, count = next(schedule)
        time, count = float(time), int(count)
        transfer = static + [next(schedule) for _ in range(count)]
        dests, cols, coefficients = [], [], []
        for row in transfer:
            dest, col, sign = edges(row[:4])
            dests.append(dest)
            cols.append(col)
            coefficients.append(-complex(float(row[4]), float(row[5]))*sign)
        dests, cols, coefficients = map(np.concatenate, (dests, cols, coefficients))

        def apply(state):
            values = coefficients*state[cols]
            return diagonal*state + np.bincount(dests, values.real, minlength=dim) + \
                1j*np.bincount(dests, values.imag, minlength=dim)

        dt = 0 if previous is None else time-previous
        term, evolved = psi.copy(), psi.copy()
        for k in range(1, order+1):
            term = (-1j*dt/k)*apply(term)
            evolved += term
        norm = np.linalg.norm(evolved)
        psi = evolved/norm
        weight = abs(psi)**2
        doublon, doublon2 = weight @ doubles, weight @ doubles**2
        hpsi = apply(psi)
        flct.append([time, number, number**2, doublon, doublon2, sz/2, (sz/2)**2, step])
        energies.append([np.vdot(psi, hpsi).real, np.vdot(hpsi, hpsi).real])
        norms.append(norm)
        previous = time
    flct = np.array(flct)
    actual = np.loadtxt(path / "output/Flct.dat")
    energy = np.loadtxt(path / "output/SS.dat")[:, 1:3]
    norm = np.loadtxt(path / "output/Norm.dat")[:, 1]
    print(f"{path.name}: dim={dim}, Flct error={np.max(abs(actual-flct)):.3e}, "
          f"SS error={np.max(abs(energy-energies)):.3e}, "
          f"Norm error={np.max(abs(norm-norms)):.3e}", flush=True)
    np.testing.assert_allclose(actual, flct, atol=2e-10, rtol=0)
    np.testing.assert_allclose(energy, energies, atol=2e-10, rtol=0)
    np.testing.assert_allclose(norm, norms, atol=2e-10, rtol=0)
    np.savetxt(path / "reference-current-h.dat", flct, fmt=["%.16f"]*7+["%d"])


if __name__ == "__main__":
    main(Path(sys.argv[1]))
