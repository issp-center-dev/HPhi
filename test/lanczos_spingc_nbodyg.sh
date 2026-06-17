#!/bin/sh
set -e

testname="lanczos_spingc_nbodyg"
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
2S = 1
Lanczos_max = 200
initial_iv = 1
EOF

rm -rf output
run_hphi log_sdry.txt "${hphi}" -sdry stan.in
printf '    NBodyG  nbodyg.def\n' >> namelist.def

cat > greenone.def <<EOF
========================
NCisAjs 2
========================
========GreenOne========
========================
0 1 0 1
1 0 1 0
EOF

cat > greentwo.def <<EOF
========================
NCisAjsCktAlt 1
========================
========GreenTwo========
========================
0 1 0 1 1 0 1 0
EOF

cat > nbodyg.def <<EOF
========================
NNBodyG 4
========================
========NBodyG==========
========================
1 0 1 0 1
1 1 0 1 0
2 0 1 0 1 1 0 1 0
2 0 1 0 0 0 1 0 0
EOF

run_hphi log_lanczos.txt "${hphi}" -e namelist.def

test -f output/zvo_NBodyG.dat || { echo "zvo_NBodyG.dat was not generated"; exit 1; }

awk -v t="${tol}" '
  NR == FNR {
    if ($1 == 0 && $2 == 1 && $3 == 0 && $4 == 1) {
      ref_re = $5;
      ref_im = $6;
      found_ref = 1;
    }
    next;
  }
  $1 == 1 && $2 == 0 && $3 == 1 && $4 == 0 && $5 == 1 {
    dr = $6 - ref_re;
    di = $7 - ref_im;
    if (dr < 0) dr = -dr;
    if (di < 0) di = -di;
    found_nbody = 1;
    ok = (found_ref && dr < t && di < t);
  }
  END { exit (found_ref && found_nbody && ok) ? 0 : 1; }
' output/zvo_cisajs.dat output/zvo_NBodyG.dat || {
  echo "NBodyG one-body value does not match cisajs"
  cat output/zvo_cisajs.dat
  cat output/zvo_NBodyG.dat
  exit 1
}

awk -v t="${tol}" '
  NR == FNR {
    if ($1 == 0 && $2 == 1 && $3 == 0 && $4 == 1 &&
        $5 == 1 && $6 == 0 && $7 == 1 && $8 == 0) {
      ref_re = $9;
      ref_im = $10;
      found_ref = 1;
    }
    next;
  }
  $1 == 2 && $2 == 0 && $3 == 1 && $4 == 0 && $5 == 1 &&
      $6 == 1 && $7 == 0 && $8 == 1 && $9 == 0 {
    dr = $10 - ref_re;
    di = $11 - ref_im;
    if (dr < 0) dr = -dr;
    if (di < 0) di = -di;
    found_nbody = 1;
    ok = (found_ref && dr < t && di < t);
  }
  END { exit (found_ref && found_nbody && ok) ? 0 : 1; }
' output/zvo_cisajscktalt.dat output/zvo_NBodyG.dat || {
  echo "NBodyG two-body value does not match cisajscktalt"
  cat output/zvo_cisajscktalt.dat
  cat output/zvo_NBodyG.dat
  exit 1
}

awk -v t="${tol}" '
  $1 == 2 && $2 == 0 && $3 == 1 && $4 == 0 && $5 == 0 &&
      $6 == 0 && $7 == 1 && $8 == 0 && $9 == 0 {
    re = $10;
    im = $11;
    if (re < 0) re = -re;
    if (im < 0) im = -im;
    found = (re < t && im < t);
  }
  END { exit found ? 0 : 1; }
' output/zvo_NBodyG.dat || {
  echo "Zero same-site NBodyG product was not written as zero"
  cat output/zvo_NBodyG.dat
  exit 1
}

echo "SpinGC NBodyG one-body, two-body, and zero-product checks passed."
