#!/bin/sh -e

# Regression test for Kondo + ncond without 2Sz.
#
# This standard input is internally promoted to KondoNConserved.  The test is
# intentionally serial-only even when the surrounding CI job defines MPIRUN:
# it verifies basis construction and CalcHS dispatch, not MPI decomposition.

testname="kondo_nconserved_dimension"
hphi="../../../src/HPhi"

mkdir -p "${testname}"
cd "${testname}"

fail() {
  echo "FAILED (${testname}): $1" >&2
  exit 1
}

write_standard_input() {
  ncond=$1
  cat > stan.in <<EOF
L = 4
model = "Kondo"
method = "FullDiag"
lattice = "chain"
t = 0.0
J = 0.0
ncond = ${ncond}
exct = 1
EOF
}

assert_dimension() {
  label=$1
  expected=$2
  log=$3

  if grep -q "Error: in sz" "${log}"; then
    cat "${log}"
    fail "${label}: sz() enumeration is inconsistent with idim_max"
  fi

  idim=$(sed -n 's/.*idim_max=[[:space:]]*\([0-9][0-9]*\).*/\1/p' "${log}" | tail -1)
  total=$(sed -n 's/.*Total dimension[[:space:]]*:[[:space:]]*\([0-9][0-9]*\).*/\1/p' "${log}" | tail -1)

  [ "x${idim}" = "x${expected}" ] || {
    cat "${log}"
    fail "${label}: idim_max ${idim:-<missing>} != expected ${expected}"
  }
  [ "x${total}" = "x${expected}" ] || {
    cat "${log}"
    fail "${label}: total dimension ${total:-<missing>} != expected ${expected}"
  }
}

run_standard_case() {
  label=$1
  calchs=$2

  rm -rf "${label}"
  mkdir -p "${label}"
  (
    cd "${label}"
    write_standard_input 2
    "${hphi}" -sdry stan.in > gen.log 2>&1 || { cat gen.log; fail "${label}: standard input generation failed"; }
    echo "CalcHS         ${calchs}" >> modpara.def

    rm -rf output
    mkdir -p output
    "${hphi}" -e namelist.def > run.log 2>&1 || { cat run.log; fail "${label}: HPhi failed"; }
    assert_dimension "${label}" 448 run.log
  )
}

run_sparse_local_spin_case() {
  label="sparse_local_spins"

  rm -rf "${label}"
  mkdir -p "${label}"
  (
    cd "${label}"
    write_standard_input 3
    "${hphi}" -sdry stan.in > gen.log 2>&1 || { cat gen.log; fail "${label}: standard input generation failed"; }

    # Replace the all-local-spin standard Kondo layout by a dilute Kondo
    # expert layout: Nsite=8, NlocalSpin=2, NsCond=6, ncond=3.
    # Expected dimension: 2^2 * sum_k C(6,k) C(6,3-k) = 880.
    cat > locspn.def <<EOF
================================
NlocalSpin     2
================================
========i_1LocSpn_0IteElc ======
================================
    0      1
    1      1
    2      0
    3      0
    4      0
    5      0
    6      0
    7      0
EOF
    echo "CalcHS         1" >> modpara.def

    rm -rf output
    mkdir -p output
    "${hphi}" -e namelist.def > run.log 2>&1 || { cat run.log; fail "${label}: HPhi failed"; }
    assert_dimension "${label}" 880 run.log
  )
}

run_standard_case standard_calchs0 0
run_standard_case standard_calchs1 1
run_sparse_local_spin_case

echo "KondoNConserved dimension and CalcHS=0/1 checks passed."
