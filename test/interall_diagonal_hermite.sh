#!/bin/sh
# Diagonal InterAll rows (n_a n_b and c+_a c_b c+_b c_a) must have real total
# coefficients. Their imaginary parts used to be dropped without a check, so
# a non-Hermitian input ran as if it were real. Rows whose imaginary parts
# cancel, such as (a a b b) and (b b a a) with conjugate values, are Hermitian
# and must still be accepted.
#
# Two sites, one up and one down electron (4 states). Each case is run, and:
#  - "reject" cases must stop with a NonHermite message;
#  - "accept" cases must give the spectrum of the matching real input.

set -e

HPHI="$(cd "$(pwd)/.." && pwd)/src/HPhi"
mkdir -p interall_diagonal_hermite
cd interall_diagonal_hermite

python3 - <<'PYEOF'
import os

CASES = {
    # reference: real rows only
    "real":        ["0 0 0 0 1 1 1 1 2.0 0.0", "0 0 1 1 1 1 0 0 0.5 0.0"],
    # n_0up n_1dn with a complex value
    "diag_complex": ["0 0 0 0 1 1 1 1 2.0 5.0", "0 0 1 1 1 1 0 0 0.5 0.0"],
    # c+_0up c_1dn c+_1dn c_0up (= n_0up - n_0up n_1dn) with a complex value
    "hund_complex": ["0 0 0 0 1 1 1 1 2.0 0.0", "0 0 1 1 1 1 0 0 0.5 3.0"],
    # n_0up n_0up (= n_0up) with a complex value
    "self_complex": ["0 0 0 0 1 1 1 1 2.0 0.0", "0 0 1 1 1 1 0 0 0.5 0.0",
                     "0 0 0 0 0 0 0 0 0.0 1.0"],
    # the n_0up n_1dn parts cancel, but n_0up keeps the imaginary part 3
    "hund_partial": ["0 0 0 0 1 1 1 1 2.0 3.0", "0 0 1 1 1 1 0 0 0.5 3.0"],
    # (a a b b) and (b b a a) with conjugate values: Hermitian, same as "real"
    "pair_cancel":  ["0 0 0 0 1 1 1 1 1.0 5.0", "1 1 1 1 0 0 0 0 1.0 -5.0",
                     "0 0 1 1 1 1 0 0 0.5 0.0"],
    # the hund row and a direct row cancel each other's n_0up n_1dn part
    "mixed_cancel": ["0 0 0 0 1 1 1 1 2.0 -3.0", "0 0 1 1 1 1 0 0 0.5 -3.0",
                     "0 0 0 0 0 0 0 0 0.0 3.0"],
}

def write(name, rows):
    os.makedirs(name, exist_ok=True)
    def w(f, s):
        with open(os.path.join(name, f), "w") as fh:
            fh.write(s)
    w("interall.def", "====\nNInterAll %d\n====\n====\n====\n" % len(rows)
      + "".join(r + "\n" for r in rows))
    trans = [(0, s, 1, s) for s in (0, 1)] + [(1, s, 0, s) for s in (0, 1)]
    w("trans.def", "====\nNTransfer %d\n====\n====\n====\n" % len(trans)
      + "".join("%d %d %d %d 1.0 0.0\n" % t for t in trans))
    w("locspn.def", "====\nNlocalSpin 0\n====\n====\n====\n0 0\n1 0\n")
    w("modpara.def", "--------------------\nModel_Parameters   0\n--------------------\n"
      "HPhi_Cal_Parameters\n--------------------\nCDataFileHead  zvo\n"
      "CParaFileHead  zqp\n--------------------\nNsite          2\n"
      "Ncond          2\n2Sz            0\nLanczos_max    2000\n"
      "initial_iv     -1\nexct           1\nLanczosEps     14\n")
    w("calcmod.def", "CalcType 2\nCalcModel 0\nReStart 0\nCalcSpec 0\nCalcEigenVec 0\n"
      "InitialVecType 0\nInputEigenVec 0\nOutputEigenVec 0\n")
    w("namelist.def", "ModPara modpara.def\nCalcMod calcmod.def\nLocSpin locspn.def\n"
      "Trans trans.def\nInterAll interall.def\n")

for name, rows in CASES.items():
    write(name, rows)
PYEOF

fail=0
for d in real pair_cancel mixed_cancel; do
  if ! ( cd "$d" && "${HPHI}" -e namelist.def > log.txt 2>&1 ); then
    echo "[interall_diagonal_hermite] $d: HPhi failed on a Hermitian input"
    tail -5 "$d/log.txt"; fail=1
  fi
done
for d in diag_complex hund_complex self_complex hund_partial; do
  if ( cd "$d" && "${HPHI}" -e namelist.def > log.txt 2>&1 ); then
    echo "[interall_diagonal_hermite] $d: a non-Hermitian input was accepted"; fail=1
  elif ! grep -q "NonHermite" "$d/log.txt"; then
    echo "[interall_diagonal_hermite] $d: failed without a NonHermite message"
    tail -5 "$d/log.txt"; fail=1
  fi
done
[ $fail -eq 0 ] || exit 1

python3 - <<'PYEOF'
import math, sys
def read(p):
    return sorted(float(l.split()[1]) for l in open(p) if l.strip())
ref = read("real/output/Eigenvalue.dat")
if len(ref) != 4 or not all(math.isfinite(x) for x in ref):
    sys.exit("[interall_diagonal_hermite] expected 4 finite eigenvalues, got %r" % ref)
for name in ("pair_cancel", "mixed_cancel"):
    e = read(name + "/output/Eigenvalue.dat")
    d = max(abs(x - y) for x, y in zip(ref, e)) if len(e) == 4 else float("inf")
    print("[interall_diagonal_hermite] %s: max |diff| = %.3e" % (name, d))
    if d > 1e-10:
        sys.exit("[interall_diagonal_hermite] %s: spectrum differs from the real input" % name)
PYEOF
echo "interall_diagonal_hermite: OK"
