#!/bin/sh -e

testname="threebody_calcmodel_validation"

mkdir -p "${testname}"
cd "${testname}"
hphi="$(pwd)/../../src/HPhi"

expect_reject() {
  name="$1"
  expected_msg="$2"
  shift 2

  rm -rf "${name}"
  mkdir -p "${name}"
  (
    cd "${name}"
    "$@"

    set +e
    "${hphi}" -e namelist.def > run.log 2>&1
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

prepare_spingc_bad_threebody() {
  cat > stan.in <<EOF
L = 4
model = "SpinGC"
method = "Lanczos"
lattice = "chain"
J = 1.0
h = 0.1
2S = 1
outputmode = "None"
EOF
  "${hphi}" -sdry stan.in > gen.log 2>&1
  cat > green3.def <<EOF
===================
num       1
===================
===================
===================
0 0 0 1 1 1 1 0 2 0 3 0
EOF
  printf "ThreeBodyG  green3.def\n" >> namelist.def
}

prepare_spin_threebody() {
  cat > stan.in <<EOF
L = 4
model = "Spin"
method = "Lanczos"
lattice = "chain"
J = 1.0
2Sz = 0
outputmode = "None"
EOF
  "${hphi}" -sdry stan.in > gen.log 2>&1
  cat > green3.def <<EOF
===================
num       1
===================
===================
===================
0 0 0 1 1 1 1 0 2 0 2 0
EOF
  printf "ThreeBodyG  green3.def\n" >> namelist.def
}

prepare_internal_calcmodel() {
  cat > namelist.def <<EOF
CalcMod   calcmod.def
ModPara   modpara.def
LocSpin   locspn.def
Trans     trans.def
InterAll  interall.def
OneBodyG  greenone.def
TwoBodyG  greentwo.def
EOF
  cat > calcmod.def <<EOF
CalcType        0
CalcModel       11
ReStart         0
CalcSpec        0
CalcEigenVec    0
InitialVecType  0
InputEigenVec   0
OutputEigenVec  0
InputHam        0
OutputHam       0
EOF
  cat > modpara.def <<EOF
--------------------
Model_Parameters   0
--------------------
HPhi_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
Nsite          4
Ncond          2
Lanczos_max    5
initial_iv     1
exct           1
LanczosEps     10
LanczosTarget  2
LargeValue     4.5
NumAve         1
ExpecInterval  20
EOF
  cat > locspn.def <<EOF
================================
NlocalSpin     0
================================
========i_1LocSpn_0IteElc ======
================================
0 0
1 0
2 0
3 0
EOF
  printf "========================\nNTransfer      0\n========================\n========i_j_s_tijs======\n========================\n" > trans.def
  printf "======================\nNInterAll      0\n======================\n========zInterAll=====\n======================\n" > interall.def
  printf "===========\nNCisAjs          0\n===========\n===========\n===========\n" > greenone.def
  printf "===========\nNCisAjsCktAlt          0\n===========\n===========\n===========\n" > greentwo.def
}

expect_reject spingc_threebody_site_mismatch "requires the 5th and 6th operators to be on the same site" prepare_spingc_bad_threebody
expect_reject spin_threebody_unsupported "canonical Spin is not supported" prepare_spin_threebody
expect_reject internal_tjn_calcmodel "reserved for internal use" prepare_internal_calcmodel

echo "ThreeBodyG and internal CalcModel validation rejects invalid inputs: OK"
