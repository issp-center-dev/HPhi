#!/bin/sh
set -e

testname="lanczos_spingc_spinone_nbodyg"
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
model = "SpinGC"
method = "Lanczos"
lattice = "chain"
L = 4
J = 1.0
2S = 2
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
1 1 1 1 0
2 0 2 0 0 1 0 1 2
2 0 0 0 2 1 2 1 0
2 0 2 0 1 0 2 0 1
EOF

run_hphi log_lanczos.txt "${hphi}" -e namelist.def
test -f output/zvo_NBodyG.dat || { echo "zvo_NBodyG.dat was not generated"; exit 1; }

awk 'NF >= 7 {count++} END{exit count == 5 ? 0 : 1}' output/zvo_NBodyG.dat || {
  echo "Unexpected NBodyG output line count"
  cat output/zvo_NBodyG.dat
  exit 1
}

awk -v t="${tol}" '
  $1 == 1 && $2 == 0 && $3 == 2 && $4 == 0 && $5 == 2 {
    re = $6; if (re < 0) re = -re;
    found_diag = (re > t);
  }
  $1 == 2 && $2 == 0 && $3 == 2 && $4 == 0 && $5 == 1 && $6 == 0 && $7 == 2 && $8 == 0 && $9 == 1 {
    re = $(NF - 1); im = $NF;
    if (re < 0) re = -re;
    if (im < 0) im = -im;
    found_zero = (re < t && im < t);
  }
  END { exit (found_diag && found_zero) ? 0 : 1; }
' output/zvo_NBodyG.dat || {
  echo "Expected nonzero diagonal and zero same-site NBodyG outputs were not found"
  cat output/zvo_NBodyG.dat
  exit 1
}

echo "SpinGC spin-one NBodyG serial output checks passed."
