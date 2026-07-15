#!/bin/sh -e

# MPI consistency test for Kondo + ncond without 2Sz in CG mode.
# The input is internally promoted to KondoNConserved.  Compare a serial run
# against a 4-rank MPI run so the particle-number-conserved, Sz-unconserved
# basis is exercised in the normal iterative-solver path, not only FullDiag.

testname="mpi_consistency_kondo_nconserved_cg"
tolerance="0.000001"

if [ -z "${MPIRUN}" ]; then
  echo "Error: MPIRUN is not set. Please set MPIRUN to run MPI tests."
  exit 1
fi

mkdir -p "${testname}"
cd "${testname}"

fail() {
  echo "FAILED (${testname}): $1" >&2
  exit 1
}

write_input() {
  cat > stan.in <<EOF
model = "Kondo"
method = "CG"
lattice = "chain"
L = 4
t = 1.0
J = 4.0
ncond = 2
exct = 1
EOF
}

extract_energy() {
  awk '$1 == "Energy" { print $2; exit }' "$1"
}

assert_total_dimension() {
  label=$1
  log=$2
  expected=$3

  total=$(sed -n 's/.*Total dimension[[:space:]]*:[[:space:]]*\([0-9][0-9]*\).*/\1/p' "${log}" | tail -1)
  [ "x${total}" = "x${expected}" ] || {
    cat "${log}"
    fail "${label}: total dimension ${total:-<missing>} != expected ${expected}"
  }
}

write_input

rm -rf output
../../src/HPhi -s stan.in > serial.log 2>&1 || { cat serial.log; fail "serial HPhi failed"; }
assert_total_dimension serial serial.log 448
serial_energy=$(extract_energy output/zvo_energy.dat)
[ -n "${serial_energy}" ] || fail "serial energy was not written"

rm -rf output
${MPIRUN} ../../src/HPhi -s stan.in > mpi.log 2>&1 || { cat mpi.log; fail "MPI HPhi failed"; }
assert_total_dimension mpi mpi.log 448
mpi_energy=$(extract_energy output/zvo_energy.dat)
[ -n "${mpi_energy}" ] || fail "MPI energy was not written"

diff=$(awk -v a="${serial_energy}" -v b="${mpi_energy}" 'BEGIN {
  d = a - b
  if (d < 0) d = -d
  printf "%.12f", d
}')

echo "KondoNConserved CG serial energy: ${serial_energy}"
echo "KondoNConserved CG MPI energy:    ${mpi_energy}"
echo "Energy difference: ${diff}"

awk -v d="${diff}" -v t="${tolerance}" 'BEGIN { exit (d < t) ? 0 : 1 }' || {
  fail "energy mismatch: |serial - MPI| = ${diff} >= ${tolerance}"
}

echo "KondoNConserved CG MPI consistency check passed."
