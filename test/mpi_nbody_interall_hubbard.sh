#!/bin/sh
set -e

testname="mpi_nbody_interall_hubbard"
hphi="../../../src/HPhi"
tol="0.00000001"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

compare_energy() {
  label="$1"
  left="$2"
  right="$3"
  diff=$(paste "${left}" "${right}" \
    | awk '$1 == "Energy" && $3 == "Energy" {d=$2-$4; if(d<0)d=-d; if(d>m)m=d} END{printf "%.12g", m+0}')
  awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d < t) ? 0 : 1}' || {
    echo "${label}: max diff ${diff}"
    paste "${left}" "${right}"
    exit 1
  }
}

rm -rf serial mpi ncond_serial ncond_mpi
mkdir -p serial mpi ncond_serial ncond_mpi

cd serial
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
printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 3
========================
========NBodyInterAll===
========================
1 3 0 3 0 0.1250000000000000 0.0000000000000000
2 3 0 0 0 1 1 1 1 0.3700000000000000 0.1100000000000000
2 1 1 1 1 0 0 3 0 0.3700000000000000 -0.1100000000000000
EOF
run_hphi log_serial.txt "${hphi}" -e namelist.def
cp output/zvo_energy.dat ../energy_serial.dat
cd ..

cp serial/*.def mpi/
cd mpi
run_hphi log_mpi.txt ${MPIRUN} "${hphi}" -e namelist.def
cp output/zvo_energy.dat ../energy_mpi.dat
cd ..

compare_energy "Serial/MPI Hubbard NBodyInterAll energy mismatch" energy_serial.dat energy_mpi.dat

grep -q "INTER process site" mpi/log_mpi.txt || {
  echo "MPI run did not print an inter-process site summary"
  cat mpi/log_mpi.txt
  exit 1
}

cd ncond_serial
cat > stan.in <<EOF
model = "Hubbard"
method = "Lanczos"
lattice = "chain"
L = 4
t = 1.0
U = 0.0
ncond = 4
Lanczos_max = 120
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
2 3 0 3 1 0 1 0 0 0.2300000000000000 0.0700000000000000
2 0 0 0 1 3 1 3 0 0.2300000000000000 -0.0700000000000000
EOF
run_hphi log_ncond_serial.txt "${hphi}" -e namelist.def
cp output/zvo_energy.dat ../energy_ncond_serial.dat
cd ..

cp ncond_serial/*.def ncond_mpi/
cd ncond_mpi
run_hphi log_ncond_mpi.txt ${MPIRUN} "${hphi}" -e namelist.def
cp output/zvo_energy.dat ../energy_ncond_mpi.dat
cd ..

compare_energy "Serial/MPI HubbardNConserved NBodyInterAll energy mismatch" \
  energy_ncond_serial.dat energy_ncond_mpi.dat

grep -q "INTER process site" ncond_mpi/log_ncond_mpi.txt || {
  echo "HubbardNConserved MPI run did not print an inter-process site summary"
  cat ncond_mpi/log_ncond_mpi.txt
  exit 1
}

echo "Hubbard and HubbardNConserved NBodyInterAll MPI np=4 match serial energies."
