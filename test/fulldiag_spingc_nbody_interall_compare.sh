#!/bin/sh
set -e

testname="fulldiag_spingc_nbody_interall_compare"
hphi="../../../src/HPhi"
tol="0.00000001"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

make_base() {
  dir="$1"
  mkdir -p "${dir}"
  cd "${dir}"
  cat > stan.in <<EOF
model = "SpinGC"
method = "FullDiag"
lattice = "chain"
L = 4
J = 0.0
2S = 1
Lanczos_max = 100
initial_iv = 1
EOF
  rm -rf output
  run_hphi log_sdry.txt "${hphi}" -sdry stan.in
  cd ..
}

rm -rf nbody interall
make_base nbody
make_base interall

cd nbody
printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 3
========================
========NBodyInterAll===
========================
2 0 1 0 1 1 0 1 0 0.2100000000000000 0.0000000000000000
2 0 1 0 0 1 0 1 1 0.3700000000000000 0.1100000000000000
2 0 0 0 1 1 1 1 0 0.3700000000000000 -0.1100000000000000
EOF
run_hphi log_nbody.txt "${hphi}" -e namelist.def
cd ..

cd interall
printf '    InterAll  interall.def\n' >> namelist.def
cat > interall.def <<EOF
========================
NInterAll 3
========================
========zInterAll=======
========================
0 1 0 1 1 0 1 0 0.2100000000000000 0.0000000000000000
0 1 0 0 1 0 1 1 0.3700000000000000 0.1100000000000000
1 1 1 0 0 0 0 1 0.3700000000000000 -0.1100000000000000
EOF
run_hphi log_interall.txt "${hphi}" -e namelist.def
cd ..

diff=$(paste nbody/output/zvo_phys.dat interall/output/zvo_phys.dat \
  | awk 'NR>1{d=$1-$6; if(d<0)d=-d; if(d>m)m=d} END{printf "%.12g", m+0}')
awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d < t) ? 0 : 1}' || {
  echo "NBodyInterAll N=2 spectrum differs from legacy InterAll: max diff ${diff}"
  paste nbody/output/zvo_phys.dat interall/output/zvo_phys.dat
  exit 1
}

echo "SpinGC NBodyInterAll N=2 matches legacy InterAll in FullDiag."
