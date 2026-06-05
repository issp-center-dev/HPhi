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
"""Independent exact-diagonalization check of the off-diagonal Sz-Sz dynamical
Green-function references used by the HPhi regression tests

    test/spectrum_spin_offdiag_szsz.sh     (canonical Spin: S=1/2 + S=1)
    test/spectrum_spingc_offdiag_szsz.sh   (SpinGC:        S=1/2 + S=1)

Model: Heisenberg ring (chain, PBC), H = J sum_<ij> S_i . S_j, J=1, L=4.
Ground-state energies: E0 = -2 (S=1/2), E0 = -6 (S=1).

G_BA(omega) = <gs| Sz_1 (z - H)^-1 Sz_0 |gs>, z = omega + E0 + i*eta.
Sz conserves total Sz and the ground state is in Sz=0, so the canonical Spin and
SpinGC results coincide (the two HPhi tests exercise different leaves but expect
the same values).  Works for any spin via twoS = 2S.

Run (requires numpy; see test/tools/README.md):
    python3 test/tools/exact_diag_spin_offdiag.py

The printed values agree with HPhi (BiCG) to ~1e-8.
"""
import numpy as np


def run(L, twoS, bonds, omegas, eta=0.1, J=1.0):
    d = twoS + 1                       # local dimension (2S+1)
    S = twoS / 2.0
    m = np.array([S - k for k in range(d)])   # Sz eigenvalues: S, S-1, ..., -S
    Sz = np.diag(m)
    # S^+ raises Sz by 1, S^- lowers it: <m+-1|S^+-|m> = sqrt(S(S+1) - m(m+-1))
    Sp = np.zeros((d, d))
    Sm = np.zeros((d, d))
    for k in range(d):
        mm = m[k]
        if k > 0:               # raise m -> m+1 (index k-1)
            Sp[k - 1, k] = np.sqrt(S * (S + 1) - mm * (mm + 1))
        if k < d - 1:           # lower m -> m-1 (index k+1)
            Sm[k + 1, k] = np.sqrt(S * (S + 1) - mm * (mm - 1))
    I = np.eye(d)

    def op_at(site, M):
        out = np.array([[1.0]])
        for s in range(L):
            out = np.kron(out, M if s == site else I)
        return out

    dim = d ** L
    H = np.zeros((dim, dim))
    for (a, b) in bonds:
        Sza, Szb = op_at(a, Sz), op_at(b, Sz)
        Spa, Sma = op_at(a, Sp), op_at(a, Sm)
        Spb, Smb = op_at(b, Sp), op_at(b, Sm)
        H += J * (Sza @ Szb + 0.5 * (Spa @ Smb + Sma @ Spb))

    evals, evecs = np.linalg.eigh(H)
    E0 = evals[0]
    gs = evecs[:, 0]
    A = op_at(0, Sz) @ gs       # Sz_0 |gs>
    B = op_at(1, Sz) @ gs       # Sz_1 |gs>
    res = []
    for w in omegas:
        z = w + E0 + 1j * eta
        x = np.linalg.solve(z * np.eye(dim) - H, A)
        res.append(np.vdot(B, x))   # <Sz_1 gs| x>
    return E0, res


def main():
    L = 4
    bonds = [(s, (s + 1) % L) for s in range(L)]  # PBC ring
    omegas = [0.0, 0.1, 0.2, 0.3, 0.4]
    for twoS, label in [(1, "S=1/2 (Half leaf)"), (2, "S=1 (General leaf)")]:
        E0, g = run(L, twoS, bonds, omegas)
        print(f"\n# Heisenberg ring L=4, {label};  E0 = {E0:.10f}")
        for w, val in zip(omegas, g):
            print(f"{w:.4f}  {val.real:.10f}  {val.imag:.10f}")


if __name__ == "__main__":
    main()
