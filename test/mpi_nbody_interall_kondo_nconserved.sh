#!/bin/sh
set -e

testname="mpi_nbody_interall_kondo_nconserved"
hphi="../../../src/HPhi"
tol="0.00000001"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

write_input() {
  cat > stan.in <<EOF
model = "Kondo"
method = "Lanczos"
lattice = "chain"
L = 4
t = 0.0
J = 0.0
ncond = 2
Lanczos_max = 1000
initial_iv = 1
EOF
}

write_nbody() {
  cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
1 0 0 0 1 0.1300000000000000 0.0200000000000000
1 0 1 0 0 0.1300000000000000 -0.0200000000000000
EOF
}

assert_total_dimension() {
  label="$1"
  log="$2"
  expected="$3"
  total=$(sed -n 's/.*Total dimension[[:space:]]*:[[:space:]]*\([0-9][0-9]*\).*/\1/p' "${log}" | tail -1)
  [ "x${total}" = "x${expected}" ] || {
    echo "${label}: total dimension ${total:-<missing>} != expected ${expected}"
    cat "${log}"
    exit 1
  }
}

rm -rf serial mpi
mkdir -p serial mpi

cd serial
write_input
run_hphi log_sdry.txt "${hphi}" -sdry stan.in
printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
write_nbody
run_hphi log_serial.txt "${hphi}" -e namelist.def
assert_total_dimension serial log_serial.txt 448
cp output/zvo_energy.dat ../energy_serial.dat
cd ..

cp serial/*.def mpi/
cd mpi
run_hphi log_mpi.txt ${MPIRUN} "${hphi}" -e namelist.def
assert_total_dimension mpi log_mpi.txt 448
cp output/zvo_energy.dat ../energy_mpi.dat
cd ..

diff=$(paste energy_serial.dat energy_mpi.dat \
  | awk '$1 == "Energy" && $3 == "Energy" {d=$2-$4; if(d<0)d=-d; if(d>m)m=d} END{printf "%.12g", m+0}')
awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d < t) ? 0 : 1}' || {
  echo "Serial/MPI KondoNConserved NBodyInterAll energy mismatch: max diff ${diff}"
  paste energy_serial.dat energy_mpi.dat
  exit 1
}

echo "KondoNConserved NBodyInterAll MPI energy matches serial energy."
