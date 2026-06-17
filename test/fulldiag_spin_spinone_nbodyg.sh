#!/bin/sh
set -e

testname="fulldiag_spin_spinone_nbodyg"
hphi="../../src/HPhi"
tol="0.00000001"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

cat > stan.in <<EOF
model = "Spin"
method = "Lanczos"
lattice = "chain"
L = 4
J = 1.0
2S = 2
2Sz = 0
Lanczos_max = 120
initial_iv = 1
EOF

rm -rf output
run_hphi log_sdry.txt "${hphi}" -sdry stan.in
printf '    NBodyG  nbodyg.def\n' >> namelist.def

cat > nbodyg.def <<EOF
========================
NNBodyG 5
========================
========NBodyG==========
========================
1 0 2 0 2
1 1 1 1 1
2 0 2 0 0 1 0 1 2
2 0 0 0 2 1 2 1 0
2 0 2 0 1 0 2 0 1
EOF

run_hphi log_lanczos.txt "${hphi}" -e namelist.def
test -f output/zvo_NBodyG.dat || { echo "zvo_NBodyG.dat was not generated"; exit 1; }
cp output/zvo_NBodyG.dat nbodyg_lanczos.dat

sed -e 's/^CalcType.*/CalcType   2/' -e 's/^OutputHam.*/OutputHam   0/' calcmod.def > calcmod.fulldiag
mv calcmod.def calcmod.lanczos
mv calcmod.fulldiag calcmod.def
rm -rf output
run_hphi log_fulldiag.txt "${hphi}" -e namelist.def
test -f output/zvo_NBodyG_eigen0.dat || {
  echo "zvo_NBodyG_eigen0.dat was not generated"
  ls output
  exit 1
}
cp output/zvo_NBodyG_eigen0.dat nbodyg_fulldiag.dat

awk -v t="${tol}" '
  function prefix(line, out) {
    out = line;
    sub(/[[:space:]]+[-+0-9.eE]+[[:space:]]+[-+0-9.eE]+$/, "", out);
    return out;
  }
  NR == FNR {
    s_prefix[NR] = prefix($0);
    s_re[NR] = $(NF - 1);
    s_im[NR] = $NF;
    n = NR;
    next;
  }
  {
    m = FNR;
    if (prefix($0) != s_prefix[m]) bad = 1;
    dr = $(NF - 1) - s_re[m];
    di = $NF - s_im[m];
    if (dr < 0) dr = -dr;
    if (di < 0) di = -di;
    if (dr >= t || di >= t) bad = 1;
  }
  END { exit (n == m && bad != 1) ? 0 : 1; }
' nbodyg_lanczos.dat nbodyg_fulldiag.dat || {
  echo "canonical Spin spin-one NBodyG FullDiag output mismatch"
  paste nbodyg_lanczos.dat nbodyg_fulldiag.dat
  exit 1
}

echo "canonical Spin spin-one NBodyG FullDiag output matches Lanczos output."
