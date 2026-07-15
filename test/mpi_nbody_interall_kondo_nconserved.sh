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
t = 1.0
J = 0.5
ncond = 2
Lanczos_max = 2000
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
1 7 0 7 1 0.3000000000000000 0.1000000000000000
1 7 1 7 0 0.3000000000000000 -0.1000000000000000
EOF
}

append_transfer_reference() {
  awk '
    $1 == "NTransfer" {
      printf "%s      %d\n", $1, $2 + 2
      next
    }
    { print }
    END {
      printf "7 0 7 1 0.3000000000000000 0.1000000000000000\n"
      printf "7 1 7 0 0.3000000000000000 -0.1000000000000000\n"
    }
  ' trans.def > trans.def.tmp
  mv trans.def.tmp trans.def
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

compare_energy() {
  label="$1"
  lhs="$2"
  rhs="$3"
  diff=$(paste "${lhs}" "${rhs}" \
    | awk '$1 == "Energy" && $3 == "Energy" {d=$2-$4; if(d<0)d=-d; if(d>m)m=d} END{printf "%.12g", m+0}')
  awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d < t) ? 0 : 1}' || {
    echo "${label}: energy mismatch: max diff ${diff}"
    paste "${lhs}" "${rhs}"
    exit 1
  }
}

rm -rf serial mpi legacy
mkdir -p serial mpi legacy

cd serial
write_input
run_hphi log_sdry.txt "${hphi}" -sdry stan.in
printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
write_nbody
run_hphi log_serial.txt "${hphi}" -e namelist.def
assert_total_dimension serial log_serial.txt 448
cp output/zvo_energy.dat ../energy_serial.dat
cd ..

cd legacy
write_input
run_hphi log_sdry.txt "${hphi}" -sdry stan.in
append_transfer_reference
run_hphi log_legacy.txt "${hphi}" -e namelist.def
assert_total_dimension legacy log_legacy.txt 448
cp output/zvo_energy.dat ../energy_legacy.dat
cd ..

cp serial/*.def mpi/
cd mpi
run_hphi log_mpi.txt ${MPIRUN} "${hphi}" -e namelist.def
assert_total_dimension mpi log_mpi.txt 448
cp output/zvo_energy.dat ../energy_mpi.dat
cd ..

grep -q "INTER process site" mpi/log_mpi.txt || {
  echo "MPI run did not print an inter-process site summary"
  cat mpi/log_mpi.txt
  exit 1
}

compare_energy "KondoNConserved NBodyInterAll serial vs legacy Trans" energy_serial.dat energy_legacy.dat
compare_energy "KondoNConserved NBodyInterAll serial vs MPI" energy_serial.dat energy_mpi.dat

echo "KondoNConserved inter-process NBodyInterAll MPI energy matches serial and legacy Trans."
