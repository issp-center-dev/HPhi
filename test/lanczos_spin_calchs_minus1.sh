#!/bin/sh -e

testname="lanczos_spin_calchs_minus1"

mkdir -p "${testname}"
cd "${testname}"
hphi="$(pwd)/../../src/HPhi"

fail() {
  echo "FAILED (${testname}): $1"
  exit 1
}

run_case() {
  dir="$1"
  calc_hs_line="$2"

  rm -rf "${dir}"
  mkdir -p "${dir}"
  cd "${dir}"

  cat > stan.in <<EOF
L = 4
model = "Spin"
method = "Lanczos"
lattice = "chain"
J = 1.0
2Sz = 0
initial_iv = 1
Lanczos_max = 100
LanczosEps = 12
outputmode = "None"
EOF

  "${hphi}" -sdry stan.in > stdface.log 2>&1 || { cat stdface.log; fail "${dir}: input generation failed"; }
  if [ -n "${calc_hs_line}" ]; then
    printf "%s\n" "${calc_hs_line}" >> modpara.def
  fi

  "${hphi}" -e namelist.def > run.log 2>&1 || { cat run.log; fail "${dir}: HPhi failed"; }
  grep -q "Error: in sz" run.log && { cat run.log; fail "${dir}: Error in sz"; }
  [ -s output/zvo_energy.dat ] || { cat run.log; fail "${dir}: missing zvo_energy.dat"; }
  awk 'NR==1{print $2}' output/zvo_energy.dat > energy.txt

  cd ..
}

run_case calchs1 ""
run_case calchs_minus1 "CalcHS         -1"

e1=$(cat calchs1/energy.txt)
em1=$(cat calchs_minus1/energy.txt)

echo "E(CalcHS=1)=${e1}  E(CalcHS=-1)=${em1}"
awk -v a="${e1}" -v b="${em1}" 'BEGIN{
  if (a == "" || b == "") { print "missing energy"; exit 1 }
  d = a - b; if (d < 0) d = -d
  printf "max |E(CalcHS=-1) - E(CalcHS=1)| = %.3e\n", d
  exit (d < 1e-10) ? 0 : 1
}' || { echo "MISMATCH"; exit 1; }

echo "Spin CalcHS=-1 == CalcHS=1: OK"
