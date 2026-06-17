#!/bin/sh -e

if [ -z "${MPIRUN}" ]; then
  MPIRUN=""
fi

testname="ctpq_invtemp_validation"

mkdir -p "${testname}"
cd "${testname}"
hphi="$(pwd)/../../src/HPhi"

prepare_case() {
  name="$1"
  restart_value="$2"
  invtemp_file="$3"

  rm -rf "${name}"
  mkdir -p "${name}"
  (
    cd "${name}"
    cat > stan.in <<EOF
L = 8
model = "SpinGC"
method = "TPQ"
lattice = "chain"
J = 1.0
NumAve = 1
Lanczos_max = 5
outputmode = "None"
EOF
    "${hphi}" -sdry stan.in > gen.log 2>&1

    cat > calcmod.def <<EOF
CalcType   5
CalcModel   4
ReStart   ${restart_value}
CalcSpec   0
CalcEigenVec   0
InitialVecType   0
InputEigenVec   0
OutputEigenVec   0
InputHam   0
OutputHam   0
OutputExVec   0
EOF
    printf "InvTemp         %s\n" "${invtemp_file}" >> namelist.def
  )
}

expect_reject() {
  name="$1"
  expected_msg="$2"

  (
    cd "${name}"
    set +e
    ${MPIRUN} "${hphi}" -e namelist.def > run.log 2>&1
    rc=$?
    set -e

    if [ "${rc}" = "0" ]; then
      echo "[${name}] ERROR: run succeeded but was expected to fail" >&2
      exit 1
    fi
    if ! grep -q "${expected_msg}" run.log; then
      echo "[${name}] ERROR: expected error message '${expected_msg}' not found" >&2
      tail -40 run.log >&2
      exit 1
    fi
  )
}

expect_accept() {
  name="$1"

  (
    cd "${name}"
    ${MPIRUN} "${hphi}" -e namelist.def > run.log 2>&1 || {
      echo "[${name}] ERROR: run failed but was expected to succeed" >&2
      tail -40 run.log >&2
      exit 1
    }
  )
}

prepare_case empty_file 0 list_inv_temp_empty.def
: > empty_file/list_inv_temp_empty.def
expect_reject empty_file "must contain at least one complete row"

prepare_case malformed_rows 0 list_inv_temp_bad.def
cat > malformed_rows/list_inv_temp_bad.def <<EOF
0.0 4 1 0
1.0 4
EOF
expect_reject malformed_rows "invalid row"

prepare_case extra_tuple 0 list_inv_temp_extra.def
cat > extra_tuple/list_inv_temp_extra.def <<EOF
0.0 4 1 0 1.0 4 1 0
EOF
expect_reject extra_tuple "invalid row"

prepare_case restart_in 3 list_inv_temp.def
cat > restart_in/list_inv_temp.def <<EOF
0.0 4 1 0
1.0 4 1 0
EOF
expect_reject restart_in "cannot be combined with Restart=2 or Restart=3"

prepare_case trailing_blank 0 list_inv_temp_trailing_blank.def
cat > trailing_blank/list_inv_temp_trailing_blank.def <<EOF
0.0 4 1 0
1.0 4 1 0

EOF
expect_accept trailing_blank

echo "cTPQ InvTemp validation rejects empty, malformed, extra-field, and restart-combined inputs, and accepts trailing blanks: OK"
