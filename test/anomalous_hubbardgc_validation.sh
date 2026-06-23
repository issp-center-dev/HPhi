#!/bin/sh
set -e

testname="anomalous_hubbardgc_validation"
hphi="../../../src/HPhi"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

expect_fail() {
  dir="$1"
  pattern="$2"
  cd "${dir}"
  set +e
  "${hphi}" -e namelist.def > log.txt 2>&1
  rc=$?
  set -e
  if [ "${rc}" -eq 0 ]; then
    echo "Expected ${dir} to fail, but it succeeded"
    exit 1
  fi
  grep -q "${pattern}" log.txt || {
    echo "Expected pattern '${pattern}' was not found in ${dir}/log.txt"
    cat log.txt
    exit 1
  }
  cd ..
}

make_base() {
  dir="$1"
  model="$2"
  mkdir -p "${dir}"
  cd "${dir}"
  if [ "${model}" = "SpinGC" ]; then
    cat > stan.in <<EOF
model = "${model}"
method = "Lanczos"
lattice = "chain"
L = 4
J = 0.0
2S = 1
Lanczos_max = 50
initial_iv = 1
EOF
  else
    cat > stan.in <<EOF
model = "${model}"
method = "Lanczos"
lattice = "chain"
L = 1
t = 0.0
U = 0.0
Lanczos_max = 50
initial_iv = 1
EOF
  fi
  run_hphi log_sdry.txt "${hphi}" -sdry stan.in
  cd ..
}

write_good_term() {
  dir="$1"
  printf '    AnomalousTerm  anomalousterm.def\n' >> "${dir}/namelist.def"
  cat > "${dir}/anomalousterm.def" <<EOF
========================
NAnomalousTerm 2
========================
========AnomalousTerm===
========================
0 0 0 0 1 0.3000000000000000 0.0000000000000000
1 0 1 0 0 0.3000000000000000 0.0000000000000000
EOF
}

write_good_g() {
  dir="$1"
  printf '    AnomalousG     anomalousg.def\n' >> "${dir}/namelist.def"
  cat > "${dir}/anomalousg.def" <<EOF
========================
NAnomalousG 1
========================
========AnomalousG=======
========================
0 0 0 0 1
EOF
}

rm -rf bad_model_g bad_same_term bad_hermite bad_calcspec_g bad_te_term

make_base bad_model_g "SpinGC"
write_good_g bad_model_g

make_base bad_same_term "HubbardGC"
printf '    AnomalousTerm  anomalousterm.def\n' >> bad_same_term/namelist.def
cat > bad_same_term/anomalousterm.def <<EOF
========================
NAnomalousTerm 1
========================
========AnomalousTerm===
========================
0 0 0 0 0 0.3000000000000000 0.0000000000000000
EOF

make_base bad_hermite "HubbardGC"
printf '    AnomalousTerm  anomalousterm.def\n' >> bad_hermite/namelist.def
cat > bad_hermite/anomalousterm.def <<EOF
========================
NAnomalousTerm 2
========================
========AnomalousTerm===
========================
0 0 0 0 1 0.3000000000000000 0.0000000000000000
1 0 0 0 1 0.3000000000000000 0.0000000000000000
EOF

make_base bad_calcspec_g "HubbardGC"
write_good_g bad_calcspec_g
sed -e 's/^CalcSpec.*/CalcSpec   1/' bad_calcspec_g/calcmod.def > bad_calcspec_g/calcmod.tmp
mv bad_calcspec_g/calcmod.tmp bad_calcspec_g/calcmod.def

make_base bad_te_term "HubbardGC"
write_good_term bad_te_term
sed -e 's/^CalcType.*/CalcType   4/' bad_te_term/calcmod.def > bad_te_term/calcmod.tmp
mv bad_te_term/calcmod.tmp bad_te_term/calcmod.def

expect_fail bad_model_g "AnomalousG is currently supported only for HubbardGC"
expect_fail bad_same_term "cannot use the same fermion operator twice"
expect_fail bad_hermite "Hermite pair"
expect_fail bad_calcspec_g "AnomalousG does not support CalcSpec"
expect_fail bad_te_term "AnomalousTerm is not supported in TimeEvolution"

echo "HubbardGC anomalous validation rejects unsupported or malformed inputs."
