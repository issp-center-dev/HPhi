#!/bin/sh
set -e

testname="lanczos_tj_kondo_nbody_interall"
tol="0.00000001"

mkdir -p "${testname}"
cd "${testname}"
hphi="$(pwd)/../../src/HPhi"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

make_tj_base() {
  dir="$1"
  calcmodel="$2"
  write_2sz="${3:-yes}"

  mkdir -p "${dir}"
  cd "${dir}"
  cat > namelist.def <<EOF
         ModPara  modpara.def
         LocSpin  locspn.def
           Trans  trans.def
        InterAll  interall.def
        OneBodyG  greenone.def
        TwoBodyG  greentwo.def
         CalcMod  calcmod.def
EOF
  cat > calcmod.def <<EOF
CalcType        0
CalcModel       ${calcmodel}
ReStart         0
CalcSpec        0
CalcEigenVec    0
InitialVecType  0
InputEigenVec   0
OutputEigenVec  0
InputHam        0
OutputHam       0
OutputExVec     0
EOF
  {
    printf -- "--------------------\nModel_Parameters   0\n--------------------\n"
    printf "HPhi_Cal_Parameters\n--------------------\n"
    printf "CDataFileHead  zvo\nCParaFileHead  zqp\n--------------------\n"
    printf "Nsite             4\n"
    if [ "${calcmodel}" = "9" ]; then
      printf "Ncond             2\n"
      if [ "${write_2sz}" = "yes" ]; then
        printf "2Sz               0\n"
      fi
    fi
    printf "Lanczos_max       2000\ninitial_iv        1\nexct              1\n"
    printf "LanczosEps        14\nLanczosTarget     2\nLargeValue        20.0\n"
    printf "NumAve            1\nExpecInterval     20\n"
  } > modpara.def
  cat > locspn.def <<EOF
================================
NlocalSpin     0
================================
========i_1LocSpn_0IteElc ======
================================
    0      0
    1      0
    2      0
    3      0
EOF
  {
    printf "========================\nNTransfer      12\n========================\n"
    printf "========i_j_s_tijs======\n========================\n"
    for b in 0 1 2; do
      j=$((b + 1))
      for s in 0 1; do
        printf "    %d     %d     %d     %d    1.000000000000000    0.0\n" "${b}" "${s}" "${j}" "${s}"
        printf "    %d     %d     %d     %d    1.000000000000000    0.0\n" "${j}" "${s}" "${b}" "${s}"
      done
    done
  } > trans.def
  cat > interall.def <<EOF
======================
NInterAll      0
======================
========zInterAll=====
======================
EOF
  printf "===========\nNCisAjs          0\n===========\n===========\n===========\n" > greenone.def
  printf "===========\nNCisAjsCktAlt          0\n===========\n===========\n===========\n" > greentwo.def
  cd ..
}

make_kondo_base() {
  dir="$1"
  model="$2"

  mkdir -p "${dir}"
  cd "${dir}"
  cat > stan.in <<EOF
model = "${model}"
method = "Lanczos"
lattice = "chain"
L = 2
t = 1.0
J = 0.0
Lanczos_max = 1000
initial_iv = 1
EOF
  if [ "${model}" = "Kondo" ]; then
    printf "nelec = 2\n2Sz = 0\n" >> stan.in
  fi
  run_hphi log_sdry.txt "${hphi}" -sdry stan.in
  cd ..
}

make_kondon_base() {
  dir="$1"

  mkdir -p "${dir}"
  cd "${dir}"
  cat > stan.in <<EOF
model = "Kondo"
method = "Lanczos"
lattice = "chain"
L = 2
t = 0.0
J = 0.0
ncond = 2
Lanczos_max = 1000
initial_iv = 1
EOF
  run_hphi log_sdry.txt "${hphi}" -sdry stan.in
  cd ..
}

write_tj_nbody() {
  cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
2 0 0 0 1 1 1 1 0 0.1700000000000000 0.0300000000000000
2 1 0 1 1 0 1 0 0 0.1700000000000000 -0.0300000000000000
EOF
}

write_tj_interall() {
  cat > interall.def <<EOF
======================
NInterAll      2
======================
========zInterAll=====
======================
0 0 0 1 1 1 1 0 0.1700000000000000 0.0300000000000000
1 0 1 1 0 1 0 0 0.1700000000000000 -0.0300000000000000
EOF
}

write_kondo_nbody() {
  cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
2 2 0 2 1 0 1 0 0 0.1700000000000000 0.0300000000000000
2 0 0 0 1 2 1 2 0 0.1700000000000000 -0.0300000000000000
EOF
}

write_tjn_nbody() {
  cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
1 2 0 2 1 0.1300000000000000 0.0200000000000000
1 2 1 2 0 0.1300000000000000 -0.0200000000000000
EOF
}

append_tjn_transfer() {
  awk '
    $1 == "NTransfer" {
      printf "%s      %d\n", $1, $2 + 2
      next
    }
    { print }
    END {
      printf "2 0 2 1 0.1300000000000000 0.0200000000000000\n"
      printf "2 1 2 0 0.1300000000000000 -0.0200000000000000\n"
    }
  ' trans.def > trans.def.tmp
  mv trans.def.tmp trans.def
}

write_kondo_interall() {
  cat > interall.def <<EOF
======================
NInterAll      2
======================
========zInterAll=====
======================
2 0 2 1 0 1 0 0 0.1700000000000000 0.0300000000000000
0 0 0 1 2 1 2 0 0.1700000000000000 -0.0300000000000000
EOF
}

write_kondon_nbody() {
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

write_kondon_transfer() {
  cat > trans.def <<EOF
========================
NTransfer      2
========================
========i_j_s_tijs======
========================
0 0 0 1 0.1300000000000000 0.0200000000000000
0 1 0 0 0.1300000000000000 -0.0200000000000000
EOF
}

compare_case() {
  label="$1"
  diff=$(paste "${label}_nbody_energy.dat" "${label}_legacy_energy.dat" \
    | awk '$1 == "Energy" && $3 == "Energy" {d=$2-$4; if(d<0)d=-d; if(d>m)m=d} END{printf "%.12g", m+0}')
  awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d < t) ? 0 : 1}' || {
    echo "${label}: NBodyInterAll Lanczos energy differs from legacy InterAll: max diff ${diff}"
    paste "${label}_nbody_energy.dat" "${label}_legacy_energy.dat"
    exit 1
  }
}

run_tj_case() {
  label="$1"
  calcmodel="$2"

  rm -rf "${label}"
  mkdir -p "${label}"
  cd "${label}"
  make_tj_base nbody "${calcmodel}"
  make_tj_base legacy "${calcmodel}"

  cd nbody
  printf '   NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  write_tj_nbody
  run_hphi log_nbody.txt ${MPIRUN} "${hphi}" -e namelist.def
  cp output/zvo_energy.dat "../${label}_nbody_energy.dat"
  cd ..

  cd legacy
  write_tj_interall
  run_hphi log_legacy.txt ${MPIRUN} "${hphi}" -e namelist.def
  cp output/zvo_energy.dat "../${label}_legacy_energy.dat"
  cd ..

  compare_case "${label}"
  cd ..
}

run_kondo_case() {
  label="$1"
  model="$2"

  rm -rf "${label}"
  mkdir -p "${label}"
  cd "${label}"
  make_kondo_base nbody "${model}"
  make_kondo_base legacy "${model}"

  cd nbody
  printf '   NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  write_kondo_nbody
  run_hphi log_nbody.txt ${MPIRUN} "${hphi}" -e namelist.def
  cp output/zvo_energy.dat "../${label}_nbody_energy.dat"
  cd ..

  cd legacy
  printf '        InterAll  interall.def\n' >> namelist.def
  write_kondo_interall
  run_hphi log_legacy.txt ${MPIRUN} "${hphi}" -e namelist.def
  cp output/zvo_energy.dat "../${label}_legacy_energy.dat"
  cd ..

  compare_case "${label}"
  cd ..
}

run_tjn_case() {
  label="$1"

  rm -rf "${label}"
  mkdir -p "${label}"
  cd "${label}"
  make_tj_base nbody 9 no
  make_tj_base legacy 9 no

  cd nbody
  printf '   NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  write_tjn_nbody
  run_hphi log_nbody.txt ${MPIRUN} "${hphi}" -e namelist.def
  cp output/zvo_energy.dat "../${label}_nbody_energy.dat"
  cd ..

  cd legacy
  append_tjn_transfer
  run_hphi log_legacy.txt ${MPIRUN} "${hphi}" -e namelist.def
  cp output/zvo_energy.dat "../${label}_legacy_energy.dat"
  cd ..

  compare_case "${label}"
  cd ..
}

run_kondon_case() {
  label="$1"

  rm -rf "${label}"
  mkdir -p "${label}"
  cd "${label}"
  make_kondon_base nbody
  make_kondon_base legacy

  cd nbody
  printf '   NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  write_kondon_nbody
  run_hphi log_nbody.txt ${MPIRUN} "${hphi}" -e namelist.def
  cp output/zvo_energy.dat "../${label}_nbody_energy.dat"
  cd ..

  cd legacy
  write_kondon_transfer
  run_hphi log_legacy.txt ${MPIRUN} "${hphi}" -e namelist.def
  cp output/zvo_energy.dat "../${label}_legacy_energy.dat"
  cd ..

  compare_case "${label}"
  cd ..
}

run_tj_case tj 9
run_tjn_case tjn
run_tj_case tjgc 10
run_kondo_case kondo Kondo
run_kondo_case kondogc KondoGC
run_kondon_case kondon

echo "tJ/tJNConserved/tJGC/Kondo/KondoGC/KondoNConserved NBodyInterAll mltply terms match legacy operators in Lanczos."
