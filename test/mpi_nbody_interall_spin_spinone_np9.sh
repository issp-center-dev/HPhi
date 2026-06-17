#!/bin/sh
set -e

testname="mpi_nbody_interall_spin_spinone_np9"
hphi="../../../src/HPhi"
tol="0.00000001"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

rm -rf serial mpi
mkdir -p serial mpi

cd serial
cat > stan.in <<EOF
model = "Spin"
method = "Lanczos"
lattice = "chain"
L = 3
J = 0.0
2S = 2
2Sz = 0
Lanczos_max = 100
initial_iv = 1
EOF
run_hphi log_sdry.txt "${hphi}" -sdry stan.in
printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
2 1 2 1 0 2 0 2 2 0.2000000000000000 0.0000000000000000
2 1 0 1 2 2 2 2 0 0.2000000000000000 0.0000000000000000
EOF
run_hphi log_serial.txt "${hphi}" -e namelist.def
cp output/zvo_energy.dat ../energy_serial.dat
cd ..

cp serial/*.def mpi/
cd mpi
run_hphi log_mpi.txt ${MPIRUN} "${hphi}" -e namelist.def
cp output/zvo_energy.dat ../energy_mpi.dat
cd ..

diff=$(paste energy_serial.dat energy_mpi.dat \
  | awk '$1 == "Energy" && $3 == "Energy" {count++; d=$2-$4; if(d<0)d=-d; if(d>m)m=d} END{if(count==0) print "missing"; else printf "%.12g", m+0}')
awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d != "missing" && d + 0 < t) ? 0 : 1}' || {
  echo "canonical Spin spin-one np=9 NBodyInterAll energy mismatch: max diff ${diff}"
  paste energy_serial.dat energy_mpi.dat
  exit 1
}

echo "canonical Spin spin-one NBodyInterAll MPI np=9 matches serial energy."
