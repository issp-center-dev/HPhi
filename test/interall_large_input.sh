#!/bin/sh
# Reading a large InterAll.def must scale better than O(N^2).
#
# Before, every row was compared with all earlier rows (to sum repeated rows)
# and with all rows (to find its Hermitian partner), so ~5e5 rows took
# minutes. Both lookups are now indexed. This test writes a Hermitian
# two-body term on 16 sites, 262,144 rows, twice:
#
#  (1) in natural order;
#  (2) shuffled, with every element split over two repeated rows (524,288 rows).
#
# One electron (16 states) keeps the diagonalization cheap. Both must give
# the same FullDiag spectrum, and the test has a timeout
# (CMakeLists.txt) that an O(N^2) reader exceeds.

set -e

HPHI="$(cd "$(pwd)/.." && pwd)/src/HPhi"
mkdir -p interall_large_input
cd interall_large_input

python3 - <<'PYEOF'
import os, random

N = 16
rng = random.Random(12345)

def orbit(a, b, c, d):
    return [((a, b, c, d), False), ((b, a, d, c), False),
            ((c, d, a, b), True), ((d, c, b, a), True)]

V = {}
for a in range(N):
    for b in range(N):
        for c in range(N):
            for d in range(N):
                if (a, b, c, d) in V:
                    continue
                v = complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) * 0.01
                keys = orbit(a, b, c, d)
                if any(k == (a, b, c, d) for k, conj in keys if conj):
                    v = complex(v.real, 0.0)
                for k, conj in keys:
                    V[k] = v.conjugate() if conj else v

def rows():
    for (a, b, c, d), v in V.items():
        for s in (0, 1):
            for t in (0, 1):
                yield (a, s, c, s, b, t, d, t), 0.5 * v

def write_case(name, lines):
    os.makedirs(name, exist_ok=True)
    with open(os.path.join(name, "interall.def"), "w") as f:
        f.write("====\nNInterAll %d\n====\n====\n====\n" % len(lines))
        f.writelines(lines)
    trans = [(i, s, (i + 1) % N, s) for i in range(N) for s in (0, 1)]
    trans += [(j, s, i, s) for i, s, j, _ in trans]
    with open(os.path.join(name, "trans.def"), "w") as f:
        f.write("====\nNTransfer %d\n====\n====\n====\n" % len(trans))
        f.writelines("%d %d %d %d 1.0 0.0\n" % r for r in trans)
    with open(os.path.join(name, "locspn.def"), "w") as f:
        f.write("====\nNlocalSpin 0\n====\n====\n====\n")
        f.writelines("%d 0\n" % i for i in range(N))
    with open(os.path.join(name, "modpara.def"), "w") as f:
        f.write("--------------------\nModel_Parameters   0\n--------------------\n"
                "HPhi_Cal_Parameters\n--------------------\nCDataFileHead  zvo\n"
                "CParaFileHead  zqp\n--------------------\nNsite          %d\n"
                "Ncond          1\n2Sz            1\nLanczos_max    2000\n"
                "initial_iv     -1\nexct           1\nLanczosEps     14\n" % N)
    with open(os.path.join(name, "calcmod.def"), "w") as f:
        f.write("CalcType 2\nCalcModel 0\nReStart 0\nCalcSpec 0\nCalcEigenVec 0\n"
                "InitialVecType 0\nInputEigenVec 0\nOutputEigenVec 0\n")
    with open(os.path.join(name, "namelist.def"), "w") as f:
        f.write("ModPara modpara.def\nCalcMod calcmod.def\nLocSpin locspn.def\n"
                "Trans trans.def\nInterAll interall.def\n")

def line(r, v):
    return " ".join(map(str, r)) + " %.17e %.17e\n" % (v.real, v.imag)

all_rows = list(rows())
write_case("ordered", [line(r, v) for r, v in all_rows])
split = []
for r, v in all_rows:
    split.append(line(r, 0.25 * v))
    split.append(line(r, 0.75 * v))
rng.shuffle(split)
write_case("shuffled_split", split)
PYEOF

for d in ordered shuffled_split; do
  ( cd "$d" && "${HPHI}" -e namelist.def > log.txt 2>&1 ) || {
    echo "[interall_large_input] HPhi failed in $d"; tail -20 "$d/log.txt"; exit 1; }
done

python3 - <<'PYEOF'
import sys
def read(p):
    return sorted(float(l.split()[1]) for l in open(p) if l.strip())
a = read("ordered/output/Eigenvalue.dat")
b = read("shuffled_split/output/Eigenvalue.dat")
if len(a) != 16 or len(a) != len(b):
    sys.exit("[interall_large_input] %d and %d eigenvalues, expected 16" % (len(a), len(b)))
d = max(abs(x - y) for x, y in zip(a, b))
print("[interall_large_input] 16 eigenvalues, max |diff| = %.3e" % d)
if d > 1e-9:
    sys.exit("[interall_large_input] spectra differ")
PYEOF
echo "interall_large_input: OK"
