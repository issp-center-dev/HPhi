#!/bin/sh
set -e

testname="green_output_format"
srcdir="$1"
hphi="../../../src/HPhi"

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || {
    echo "Command failed: $*" >&2
    tail -80 "${log}" >&2
    exit 1
  }
}

assert_file() {
  file="$1"
  if [ ! -s "${file}" ]; then
    fail "expected non-empty file ${file}"
  fi
}

assert_no_match() {
  pattern="$1"
  set -- ${pattern}
  if [ -e "$1" ]; then
    ls -l ${pattern} >&2
    fail "unexpected file matching ${pattern}"
  fi
}

append_aggregate_mode() {
  file="$1"
  printf "OutputGreenFormat 1\n" >> "${file}"
}

check_tpq_rows() {
  file="$1"
  awk '
    NF != 8 { exit 1 }
    $1 !~ /^[0-9]+$/ { exit 1 }
    $2 !~ /^[0-9]+$/ { exit 1 }
    { seen[$1 ":" $2] = 1 }
    END {
      for (key in seen) count++
      if (count < 2) exit 1
    }
  ' "${file}" || fail "TPQ aggregate rows must have set and step columns: ${file}"
}

check_tpq_nbody_rows() {
  file="$1"
  awk '
    NF != 9 { exit 1 }
    $1 !~ /^[0-9]+$/ { exit 1 }
    $2 !~ /^[0-9]+$/ { exit 1 }
    $3 != 1 { exit 1 }
    { seen[$1 ":" $2] = 1 }
    END {
      for (key in seen) count++
      if (count < 2) exit 1
    }
  ' "${file}" || fail "TPQ NBodyG aggregate rows must have set, step, and n columns: ${file}"
}

check_single_index_rows() {
  file="$1"
  label="$2"
  awk '
    NF != 7 { exit 1 }
    $1 !~ /^[0-9]+$/ { exit 1 }
  ' "${file}" || fail "${label} aggregate rows must have one leading index column: ${file}"
}

rm -rf "${testname}"
mkdir -p "${testname}"
cd "${testname}"

mkdir tpq
(
  cd tpq
  cat > stan.in <<EOF
L = 8
model = "Spin"
method = "TPQ"
lattice = "chain"
J = 1.0
2Sz = 0
Lanczos_max = 4
LargeValue = 5
NumAve = 1
ExpecInterval = 1
outputmode = "correlation"
EOF
  run_hphi sdry.log "${hphi}" -sdry stan.in
  append_aggregate_mode calcmod.def
  printf '    NBodyG  nbodyg.def\n' >> namelist.def
  cat > nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 0 0 0
EOF
  run_hphi run.log "${hphi}" -e namelist.def
  assert_file output/zvo_cisajs_tpq.dat
  assert_file output/zvo_cisajscktalt_tpq.dat
  assert_file output/zvo_NBodyG_tpq.dat
  assert_no_match "output/zvo_cisajs_set*step*.dat"
  assert_no_match "output/zvo_cisajscktalt_set*step*.dat"
  assert_no_match "output/zvo_NBodyG_set*step*.dat"
  check_tpq_rows output/zvo_cisajs_tpq.dat
  check_tpq_nbody_rows output/zvo_NBodyG_tpq.dat
)

mkdir te
(
  cd te
  python3 "${srcdir}/test/testTECalc.py" -p "${hphi}" -m "Spin" > gen.log 2>&1
  append_aggregate_mode calcmod2.def
  sed -e 's/^Lanczos_max.*/Lanczos_max    3/' modpara2.def > modpara2.tmp
  mv modpara2.tmp modpara2.def
  rm -f output/zvo_cisajs_step*.dat output/zvo_cisajscktalt_step*.dat
  rm -f output/zvo_cisajs_te.dat output/zvo_cisajscktalt_te.dat
  run_hphi run_te.log "${hphi}" -e namelist2.def
  assert_file output/zvo_cisajs_te.dat
  assert_file output/zvo_cisajscktalt_te.dat
  assert_no_match "output/zvo_cisajs_step*.dat"
  assert_no_match "output/zvo_cisajscktalt_step*.dat"
  check_single_index_rows output/zvo_cisajs_te.dat "TimeEvolution"
)

mkdir fulldiag
(
  cd fulldiag
  cat > stan.in <<EOF
L = 2
model = "Hubbard"
method = "FullDiag"
lattice = "chain"
t = 1.0
U = 1.0
nelec = 2
2Sz = 0
outputmode = "correlation"
EOF
  run_hphi sdry.log "${hphi}" -sdry stan.in
  append_aggregate_mode calcmod.def
  run_hphi run.log "${hphi}" -e namelist.def
  assert_file output/zvo_cisajs_eigen.dat
  assert_file output/zvo_cisajscktalt_eigen.dat
  assert_no_match "output/zvo_cisajs_eigen[0-9]*.dat"
  assert_no_match "output/zvo_cisajscktalt_eigen[0-9]*.dat"
  check_single_index_rows output/zvo_cisajs_eigen.dat "FullDiag"
)

mkdir invalid
(
  cd invalid
  cat > stan.in <<EOF
L = 4
model = "Spin"
method = "Lanczos"
lattice = "chain"
J = 1.0
2Sz = 0
outputmode = "None"
EOF
  run_hphi sdry.log "${hphi}" -sdry stan.in
  printf "OutputGreenFormat 2\n" >> calcmod.def
  set +e
  "${hphi}" -e namelist.def > invalid.log 2>&1
  rc=$?
  set -e
  if [ "${rc}" -eq 0 ]; then
    fail "OutputGreenFormat=2 should be rejected"
  fi
  grep -q "OutputGreenFormat" invalid.log || {
    cat invalid.log >&2
    fail "OutputGreenFormat validation message was not found"
  }
)

echo "OutputGreenFormat aggregate output creates indexed aggregate files and rejects invalid mode: OK"
