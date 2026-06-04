#!/usr/bin/env python3
# HPhi  -  Quantum Lattice Model Simulator
# Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Independent exact-diagonalization check of the off-diagonal dynamical Green
function references used by the HPhi regression tests

    test/spectrum_hubbard_offdiag_single.sh        (single excitation)
    test/spectrum_hubbard_offdiag_pair_density.sh  (pair / density excitation)
    test/spectrum_hubbard_offdiag_single_mpi.sh    (MPI consistency)

Model: 1D Hubbard 4-site ring (t=1, U=4, PBC, Nup=Ndown=2), E0 = -2.1027484835.

Convention (matches HPhi): G_BA(omega) = <gs| B^dagger (z - H)^-1 A |gs>
with z = omega + OmegaOrg + i*OmegaIm, OmegaOrg = E0.  numpy.vdot conjugates its
first argument, mirroring HPhi's VecProdMPI(<B phi|, .) projection.

Run (requires numpy; see test/tools/README.md):
    pyenv activate venv312
    python test/tools/exact_diag_hubbard_offdiag.py

The printed values agree with HPhi (BiCG) to ~1e-8 (single) / ~1e-6 (pair).
"""
import numpy as np
import itertools

L = 4
t = 1.0
U = 4.0


def popcount(x):
    return bin(x).count("1")


def orb(site, spin):
    # spin-orbital index: o = 2*site + spin (spin 0 = up, 1 = down)
    return 2 * site + spin


def fermion_sign(state, o):
    """Sign from anticommuting c_o^(dag) past lower-index occupied orbitals."""
    return -1 if (popcount(state & ((1 << o) - 1)) & 1) else 1


def annihilate(state, o):
    if not (state >> o) & 1:
        return None, 0
    return state & ~(1 << o), fermion_sign(state, o)


def create(state, o):
    if (state >> o) & 1:
        return None, 0
    return state | (1 << o), fermion_sign(state, o)


def build_basis(nup, ndown):
    up_orbs = [orb(s, 0) for s in range(L)]
    dn_orbs = [orb(s, 1) for s in range(L)]
    states = []
    for ups in itertools.combinations(up_orbs, nup):
        for dns in itertools.combinations(dn_orbs, ndown):
            st = 0
            for o in ups + dns:
                st |= (1 << o)
            states.append(st)
    states.sort()
    index = {st: i for i, st in enumerate(states)}
    return states, index


def build_H(states, index):
    n = len(states)
    H = np.zeros((n, n), dtype=complex)
    bonds = [(s, (s + 1) % L) for s in range(L)]  # PBC chain
    for i, st in enumerate(states):
        # U n_up n_down
        diag = 0.0
        for s in range(L):
            if ((st >> orb(s, 0)) & 1) and ((st >> orb(s, 1)) & 1):
                diag += U
        H[i, i] += diag
        # -t hopping c^dag_{i s} c_{j s} + h.c.
        for (a, b) in bonds:
            for spin in (0, 1):
                oa, ob = orb(a, spin), orb(b, spin)
                for (o1, o2) in ((oa, ob), (ob, oa)):  # both directions
                    s1, sgn1 = annihilate(st, o2)
                    if s1 is None:
                        continue
                    s2, sgn2 = create(s1, o1)
                    if s2 is None:
                        continue
                    j = index.get(s2)
                    if j is not None:
                        H[j, i] += -t * sgn1 * sgn2
    return H


def apply_single(vec, states_from, idx_to, states_to, site, spin, itype):
    """itype==1: creation c^dag; else annihilation c. Returns vector in target basis."""
    out = np.zeros(len(states_to), dtype=complex)
    o = orb(site, spin)
    for i, st in enumerate(states_from):
        if vec[i] == 0:
            continue
        if itype == 1:
            st2, sgn = create(st, o)
        else:
            st2, sgn = annihilate(st, o)
        if st2 is None:
            continue
        j = idx_to.get(st2)
        if j is not None:
            out[j] += sgn * vec[i]
    return out


def apply_number(vec, states, site, spin):
    """n_{site,spin} (diagonal, preserves sector)."""
    o = orb(site, spin)
    out = np.zeros(len(states), dtype=complex)
    for i, st in enumerate(states):
        if (st >> o) & 1:
            out[i] = vec[i]
    return out


def main():
    # Ground state in (Nup=2, Ndown=2)
    states0, idx0 = build_basis(2, 2)
    H0 = build_H(states0, idx0)
    evals, evecs = np.linalg.eigh(H0)
    E0 = evals[0]
    gs = evecs[:, 0]
    print(f"E0 = {E0:.10f}  (dim {len(states0)})")

    omega_org = E0
    eta = 0.1
    omegas = [0.0, 0.1, 0.2, 0.3, 0.4]

    def green(a_ket, b_bra, states_mid, idx_mid):
        """G(omega) = <gs| B^dag (z-H_mid)^-1 A |gs>, z = omega + OmegaOrg + i eta."""
        h_mid = build_H(states_mid, idx_mid)
        out = []
        for w in omegas:
            z = w + omega_org + 1j * eta
            x = np.linalg.solve(z * np.eye(len(states_mid)) - h_mid, a_ket)
            out.append(np.vdot(b_bra, x))  # <B bra| x>
        return out

    # single: A = c_{1up} (annihilate), B = c_{0up}; mid sector (Nup=1, Ndown=2)
    states_m, idx_m = build_basis(1, 2)
    a_single = apply_single(gs, states0, idx_m, states_m, site=1, spin=0, itype=0)
    b_single = apply_single(gs, states0, idx_m, states_m, site=0, spin=0, itype=0)
    print("\n# single  G_BA = <gs| c0up^dag (z-H)^-1 c1up |gs>")
    for w, g in zip(omegas, green(a_single, b_single, states_m, idx_m)):
        print(f"{w:.4f}  {g.real:.10f}  {g.imag:.10f}")

    # pair-density: A = n_{0up}, B = n_{1up}; same (Nup=2, Ndown=2) sector
    a_pair = apply_number(gs, states0, site=0, spin=0)
    b_pair = apply_number(gs, states0, site=1, spin=0)
    print("\n# pair-density  G_BA = <gs| n1up (z-H)^-1 n0up |gs>")
    for w, g in zip(omegas, green(a_pair, b_pair, states0, idx0)):
        print(f"{w:.4f}  {g.real:.10f}  {g.imag:.10f}")


if __name__ == "__main__":
    main()
