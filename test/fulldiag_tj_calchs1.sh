#!/bin/sh -e

testname="fulldiag_tj_calchs1"

mkdir -p ${testname}
cd ${testname}
hphi="$(pwd)/../../src/HPhi"

fail() {
  echo "FAILED (${testname}): $1"
  exit 1
}

gen_input() {
  dir=$1
  calc_hs=$2
  write_2sz=$3
  nsite=$4
  ncond=$5
  twosz=$6

  rm -rf "${dir}"
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
CalcType        2
CalcModel       9
ReStart         0
CalcSpec        0
CalcEigenVec    0
InitialVecType  0
InputEigenVec   0
OutputEigenVec  0
InputHam        0
OutputHam       0
EOF

  {
    printf -- "--------------------\n"
    printf "Model_Parameters   0\n"
    printf -- "--------------------\n"
    printf "HPhi_Cal_Parameters\n"
    printf -- "--------------------\n"
    printf "CDataFileHead  zvo\n"
    printf "CParaFileHead  zqp\n"
    printf -- "--------------------\n"
    printf "Nsite             %d\n" "${nsite}"
    printf "Ncond             %d\n" "${ncond}"
    if [ "${write_2sz}" = "yes" ]; then
      printf "2Sz               %d\n" "${twosz}"
    fi
    printf "CalcHS            %d\n" "${calc_hs}"
    printf "Lanczos_max       20\n"
    printf "initial_iv        1\n"
    printf "exct              1\n"
    printf "LanczosEps        10\n"
    printf "LanczosTarget     2\n"
    printf "LargeValue        4.5\n"
    printf "NumAve            1\n"
    printf "ExpecInterval     20\n"
  } > modpara.def

  cat > locspn.def <<EOF
================================
NlocalSpin     0
================================
========i_1LocSpn_0IteElc ======
================================
EOF
  isite=0
  while [ "${isite}" -lt "${nsite}" ]; do
    printf "%5d      0\n" "${isite}" >> locspn.def
    isite=$((isite + 1))
  done

  cat > trans.def <<EOF
========================
NTransfer      0
========================
========i_j_s_tijs======
========================
EOF

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

run_case() {
  dir=$1
  expected_dim=$2

  cd "${dir}"
  rm -rf output
  mkdir -p output
  "${hphi}" -e namelist.def > run.log 2>&1 || { cat run.log; fail "${dir}: HPhi failed"; }
  grep -q "Error: in sz" run.log && { cat run.log; fail "${dir}: Error in sz"; }

  dim=$(grep "Total dimension :" run.log | tail -1 | awk '{print $NF}')
  [ "x${dim}" = "x${expected_dim}" ] || { cat run.log; fail "${dir}: dimension ${dim} != ${expected_dim}"; }
  cd ..
}

gen_input tj_even_hs0 0 yes 4 2 0
gen_input tj_even_hs1 1 yes 4 2 0
gen_input tjn_even_hs0 0 no 4 2 0
gen_input tjn_even_hs1 1 no 4 2 0
gen_input tj_odd_hs0 0 yes 5 3 1
gen_input tj_odd_hs1 1 yes 5 3 1
gen_input tjn_odd_hs0 0 no 5 3 0
gen_input tjn_odd_hs1 1 no 5 3 0

run_case tj_even_hs0 12
run_case tj_even_hs1 12
run_case tjn_even_hs0 24
run_case tjn_even_hs1 24
run_case tj_odd_hs0 30
run_case tj_odd_hs1 30
run_case tjn_odd_hs0 80
run_case tjn_odd_hs1 80

diff -u tj_even_hs0/output/Eigenvalue.dat tj_even_hs1/output/Eigenvalue.dat
diff -u tj_even_hs0/output/zvo_phys_Nup1_Ndown1.dat tj_even_hs1/output/zvo_phys_Nup1_Ndown1.dat
diff -u tjn_even_hs0/output/Eigenvalue.dat tjn_even_hs1/output/Eigenvalue.dat
diff -u tjn_even_hs0/output/zvo_phys_Nup0_Ndown0.dat tjn_even_hs1/output/zvo_phys_Nup0_Ndown0.dat
diff -u tj_odd_hs0/output/Eigenvalue.dat tj_odd_hs1/output/Eigenvalue.dat
diff -u tj_odd_hs0/output/zvo_phys_Nup2_Ndown1.dat tj_odd_hs1/output/zvo_phys_Nup2_Ndown1.dat
diff -u tjn_odd_hs0/output/Eigenvalue.dat tjn_odd_hs1/output/Eigenvalue.dat
diff -u tjn_odd_hs0/output/zvo_phys_Nup0_Ndown0.dat tjn_odd_hs1/output/zvo_phys_Nup0_Ndown0.dat

echo "PASSED (${testname})"
exit 0
