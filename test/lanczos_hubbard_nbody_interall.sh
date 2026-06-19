#!/bin/sh
set -e

testname="lanczos_hubbard_nbody_interall"
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
model = "Hubbard"
method = "Lanczos"
lattice = "chain"
L = 4
t = 1.0
U = 0.0
nelec = 4
2Sz = 0
Lanczos_max = 120
initial_iv = 1
EOF
  run_hphi log_sdry.txt "${hphi}" -sdry stan.in
  cd ..
}

rm -rf nbody legacy
make_base nbody
make_base legacy

cd nbody
printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
2 3 0 0 0 1 1 1 1 0.3700000000000000 0.1100000000000000
2 1 1 1 1 0 0 3 0 0.3700000000000000 -0.1100000000000000
EOF
run_hphi log_nbody.txt "${hphi}" -e namelist.def
cp output/zvo_energy.dat ../energy_nbody.dat
cd ..

cd legacy
printf '    InterAll  interall.def\n' >> namelist.def
cat > interall.def <<EOF
========================
NInterAll 2
========================
========zInterAll=======
========================
3 0 0 0 1 1 1 1 0.3700000000000000 0.1100000000000000
1 1 1 1 0 0 3 0 0.3700000000000000 -0.1100000000000000
EOF
run_hphi log_legacy.txt "${hphi}" -e namelist.def
cp output/zvo_energy.dat ../energy_legacy.dat
cd ..

diff=$(paste energy_nbody.dat energy_legacy.dat \
  | awk '$1 == "Energy" && $3 == "Energy" {d=$2-$4; if(d<0)d=-d; if(d>m)m=d} END{printf "%.12g", m+0}')
awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d < t) ? 0 : 1}' || {
  echo "Hubbard NBodyInterAll N=2 energy differs from legacy InterAll: max diff ${diff}"
  paste energy_nbody.dat energy_legacy.dat
  exit 1
}

echo "Hubbard NBodyInterAll N=2 matches legacy InterAll in Lanczos."
