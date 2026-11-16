#!/bin/sh
# FullDiag for SpinlessFermion (CalcModel 7) and SpinlessFermionGC (CalcModel 8).
#
# makeHam had no spinless branch (issue #293), so the hopping was dropped from
# the FullDiag Hamiltonian and every eigenvalue came out wrong. This test
# compares the full spectrum against two independent references:
#
#  (1) canonical: the same Hamiltonian written as a Hubbard model with every
#      particle in the up spin (CalcModel 0, 2Sz = Ncond), whose FullDiag
#      branch is long established;
#  (2) grand canonical: the union of the canonical spectra over Ncond = 1..L
#      plus the vacuum energy 0.
#
# The Hamiltonian has complex open-chain hopping, a chemical potential on one
# site (a diagonal Trans row), nearest-neighbour CoulombInter, and a complex
# off-diagonal NBodyInterAll pair term, so the hopping, its Hermitian partner,
# the diagonal part and the NBodyInterAll branch are all exercised.
# The particle number written to the phys file is checked as well.

set -e

L=5
NCOND=2

mkdir -p fulldiag_spinless
cd fulldiag_spinless

write_common() {  # $1 = dir
  cat > "$1/locspn.def" <<EOF
================================
NlocalSpin     0
================================
========i_1LocSpn_0IteElc ======
================================
    0      0
    1      0
    2      0
    3      0
    4      0
EOF
  cat > "$1/trans.def" <<EOF
========================
NTransfer      9
========================
========i_j_s_tijs======
========================
0 0 1 0 1.0 0.3
1 0 0 0 1.0 -0.3
1 0 2 0 0.8 0.0
2 0 1 0 0.8 0.0
2 0 3 0 1.1 -0.2
3 0 2 0 1.1 0.2
3 0 4 0 0.9 0.1
4 0 3 0 0.9 -0.1
2 0 2 0 0.35 0.0
EOF
  cat > "$1/coulombinter.def" <<EOF
========================
NCoulombInter 4
========================
========CoulombInter====
========================
0 1 0.7
1 2 0.4
2 3 0.9
3 4 0.6
EOF
  cat > "$1/nbodyinterall.def" <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
2 0 0 2 0 1 0 3 0 0.45 0.15
2 3 0 1 0 2 0 0 0 0.45 -0.15
EOF
  cat > "$1/namelist.def" <<EOF
         ModPara  modpara.def
         LocSpin  locspn.def
           Trans  trans.def
    CoulombInter  coulombinter.def
   NBodyInterAll  nbodyinterall.def
         CalcMod  calcmod.def
EOF
}

write_calcmod() {  # $1 = dir, $2 = CalcModel
  cat > "$1/calcmod.def" <<EOF
CalcType        2
CalcModel       $2
ReStart         0
CalcSpec        0
CalcEigenVec    0
InitialVecType  0
InputEigenVec   0
OutputEigenVec  0
EOF
}

write_modpara() {  # $1 = dir, $2 = extra lines (Ncond / 2Sz), may be empty
  cat > "$1/modpara.def" <<EOF
--------------------
Model_Parameters   0
--------------------
HPhi_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
Nsite          ${L}
$2
Lanczos_max    2000
initial_iv     -1
exct           1
LanczosEps     14
LanczosTarget  2
LargeValue     12.0
NumAve         5
ExpecInterval  20
EOF
}

run_case() {  # $1 = dir
  ( cd "$1" && ${MPIRUNFC} ../../../src/HPhi -e namelist.def > log.txt 2>&1 ) || {
    echo "[fulldiag_spinless] HPhi failed in $1"; cat "$1/log.txt"; exit 1; }
}

compare_spectra() {  # $1 = label, $2 = file A, $3 = file B
  python3 - "$1" "$2" "$3" <<'EOF'
import sys
label, fa, fb = sys.argv[1:4]
a = sorted(float(l.split()[1]) for l in open(fa) if l.strip())
b = sorted(float(l.split()[0]) for l in open(fb) if l.strip())
if len(a) != len(b):
    sys.exit("[fulldiag_spinless] %s: %d eigenvalues, expected %d" % (label, len(a), len(b)))
d = max(abs(x - y) for x, y in zip(a, b))
print("[fulldiag_spinless] %s: %d eigenvalues, max |diff| = %.3e" % (label, len(a), d))
if d > 1e-8:
    sys.exit("[fulldiag_spinless] %s: spectra differ" % label)
EOF
}

# (1) canonical spinless vs. all-up Hubbard
for d in spinless hubbard_allup; do
  mkdir -p "$d"
  write_common "$d"
done
write_calcmod spinless 7
write_modpara spinless "Ncond          ${NCOND}"
write_calcmod hubbard_allup 0
write_modpara hubbard_allup "Ncond          ${NCOND}
2Sz            ${NCOND}"
run_case spinless
run_case hubbard_allup
awk '{print $2}' hubbard_allup/output/Eigenvalue.dat > ref_canonical.dat
compare_spectra "CalcModel 7 vs all-up Hubbard" spinless/output/Eigenvalue.dat ref_canonical.dat

# The phys file lists N = Ncond for every state.
awk -v n="${NCOND}" 'NR > 1 { d = $2 - n; if (d < 0) d = -d; if (d > 1e-8) bad = 1 }
  END { exit bad }' "spinless/output/zvo_phys_Nup${NCOND}_Ndown0.dat" || {
  echo "[fulldiag_spinless] canonical phys file does not report N = ${NCOND}"
  cat "spinless/output/zvo_phys_Nup${NCOND}_Ndown0.dat"; exit 1; }

# (2) grand canonical spinless vs. the union of canonical sectors
mkdir -p spinless_gc
write_common spinless_gc
write_calcmod spinless_gc 8
write_modpara spinless_gc ""
run_case spinless_gc
echo "0.0" > ref_gc.dat
n=1
while [ "$n" -le "$L" ]; do
  d="sector_N$n"
  mkdir -p "$d"
  write_common "$d"
  write_calcmod "$d" 7
  write_modpara "$d" "Ncond          $n"
  run_case "$d"
  awk '{print $2}' "$d/output/Eigenvalue.dat" >> ref_gc.dat
  n=$((n + 1))
done
compare_spectra "CalcModel 8 vs union of CalcModel 7 sectors" spinless_gc/output/Eigenvalue.dat ref_gc.dat

# Trace of N over the full Fock space: sum_N N * C(L, N) = L * 2^(L-1).
awk -v want=$((L * (1 << (L - 1)))) 'NR > 1 { s += $2 }
  END { d = s - want; if (d < 0) d = -d; printf "[fulldiag_spinless] GC sum of N = %.8f (expected %d)\n", s, want; exit (d > 1e-6) }' \
  spinless_gc/output/zvo_phys.dat || {
  echo "[fulldiag_spinless] GC phys file has the wrong particle numbers"; exit 1; }

echo "fulldiag_spinless: OK"
