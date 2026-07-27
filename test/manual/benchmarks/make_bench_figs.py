#!/usr/bin/env python3
"""Build the manual-appendix benchmark figures and RST table rows
from the timers collected on clavius (bench_times.csv:
L,tag,lapackdiag,calcphys)."""
import csv
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv_path, figdir = sys.argv[1], sys.argv[2]

t = defaultdict(dict)  # t[L][tag] = (diag, phys)
with open(csv_path) as f:
    for row in csv.reader(f):
        L, tag, diag, phys = row
        t[int(L)][tag] = (float(diag), float(phys))

Ls = sorted(t)
Ns = [2**L for L in Ls]

diag_series = [
    ("s0", "Solver 0 (LAPACK, 1proc x 16threads)", "o-"),
    ("s1", "Solver 1 (ScaLAPACK, 16procs)", "s-"),
    ("s3cpu", "Solver 3 (ELPA CPU, 16procs)", "^-"),
    ("g1", "Solver 3 (ELPA GPU, 1x A100)", "D-"),
    ("g2", "Solver 3 (ELPA GPU, 2x A100)", "v-"),
    ("g4", "Solver 3 (ELPA GPU, 4x A100 / 2 nodes)", "*-"),
]
expec_series = [
    ("m0", "ExpecMode 0", "o-"),
    ("m1", "ExpecMode 1", "s-"),
    ("m2", "ExpecMode 2", "^-"),
]


def plot(series, idx, ylabel, outfile):
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for tag, label, style in series:
        ys = [t[L][tag][idx] for L in Ls]
        ax.plot(Ns, ys, style, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Matrix dimension $N$")
    ax.set_ylabel(ylabel)
    ax.set_xticks(Ns)
    ax.set_xticklabels([str(n) for n in Ns])
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{figdir}/{outfile}", dpi=150)


plot(diag_series, 0, "Diagonalization time (s)", "fulldiag_solver_bench.png")
plot(expec_series, 1, "Observable-evaluation time (s)", "fulldiag_expecmode_bench.png")


def fmt(x):
    return f"{x:.2f}" if x < 100 else f"{x:.1f}"


print("DIAG TABLE ROWS:")
for L in Ls:
    cells = [str(2**L)] + [fmt(t[L][tag][0]) for tag, _, _ in diag_series]
    print("   " + ", ".join(f'"{c}"' for c in cells))
print("EXPEC TABLE ROWS:")
for L in Ls:
    cells = [str(2**L)] + [fmt(t[L][tag][1]) for tag, _, _ in expec_series]
    print("   " + ", ".join(f'"{c}"' for c in cells))
