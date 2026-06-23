#!/bin/sh
set -e

testname="nbody_tj_kondo_validation"
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

check_nbodyg_output() {
  dir="$1"
  if ! ls "${dir}"/output/zvo_NBodyG*.dat >/dev/null 2>&1; then
    echo "Expected ${dir} to write NBodyG output"
    ls -R "${dir}"/output
    exit 1
  fi
}

make_tj_base() {
  dir="$1"
  calcmodel="$2"
  write_2sz="$3"

  mkdir -p "${dir}"
  cd "${dir}"
  cat > namelist.def <<EOF
         ModPara  modpara.def
         LocSpin  locspn.def
           Trans  trans.def
        InterAll  interall.def
        OneBodyG  greenone.def
        TwoBodyG  greentwo.def
   NBodyInterAll  nbodyinterall.def
          NBodyG  nbodyg.def
         CalcMod  calcmod.def
EOF

  cat > calcmod.def <<EOF
CalcType        2
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
    printf -- "--------------------\n"
    printf "Model_Parameters   0\n"
    printf -- "--------------------\n"
    printf "HPhi_Cal_Parameters\n"
    printf -- "--------------------\n"
    printf "CDataFileHead  zvo\n"
    printf "CParaFileHead  zqp\n"
    printf -- "--------------------\n"
    printf "Nsite             4\n"
    if [ "${calcmodel}" = "9" ]; then
      printf "Ncond             2\n"
      if [ "${write_2sz}" = "yes" ]; then
        printf "2Sz               0\n"
      fi
    fi
    printf "Lanczos_max       50\n"
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
    0      0
    1      0
    2      0
    3      0
EOF

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

make_kondo_base() {
  dir="$1"
  model="$2"
  write_sector="$3"

  mkdir -p "${dir}"
  cd "${dir}"
  cat > stan.in <<EOF
model = "${model}"
method = "FullDiag"
lattice = "chain"
L = 2
t = 0.0
J = 0.0
Lanczos_max = 50
initial_iv = 1
EOF
  if [ "${write_sector}" = "yes" ]; then
    {
      printf "nelec = 2\n"
      printf "2Sz = 0\n"
    } >> stan.in
  elif [ "${write_sector}" = "ncond" ]; then
    printf "nelec = 2\n" >> stan.in
  fi
  run_hphi log_sdry.txt "${hphi}" -sdry stan.in
  printf '   NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  printf '          NBodyG  nbodyg.def\n' >> namelist.def
  cd ..
}

rm -rf accept_tj accept_tjgc accept_kondo accept_kondogc \
  bad_tj_particle bad_kondo_nbodyg_particle bad_tjn bad_kondon bad_kondo_local_cross

make_tj_base accept_tj 9 yes
cat > accept_tj/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 0 0 0 0.1000000000000000 0.0000000000000000
EOF
cat > accept_tj/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 0 0 0
EOF

make_tj_base accept_tjgc 10 no
cat > accept_tjgc/nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
1 1 0 0 0 0.0500000000000000 0.0000000000000000
1 0 0 1 0 0.0500000000000000 0.0000000000000000
EOF
cat > accept_tjgc/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 1 0 0 1
EOF

make_kondo_base accept_kondo Kondo yes
cat > accept_kondo/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 0 0 0 0.1000000000000000 0.0000000000000000
EOF
cat > accept_kondo/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 0 0 0
EOF

make_kondo_base accept_kondogc KondoGC no
cat > accept_kondogc/nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
2 2 0 2 1 0 1 0 0 0.1000000000000000 0.0000000000000000
2 0 0 0 1 2 1 2 0 0.1000000000000000 0.0000000000000000
EOF
cat > accept_kondogc/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 1 0 0
EOF

make_tj_base bad_tj_particle 9 yes
cat > bad_tj_particle/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 1 0 1 1 0.1000000000000000 0.0000000000000000
EOF
cat > bad_tj_particle/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_kondo_base bad_kondo_nbodyg_particle Kondo yes
cat > bad_kondo_nbodyg_particle/nbodyinterall.def <<EOF
========================
NNBodyInterAll 0
========================
========NBodyInterAll===
========================
EOF
cat > bad_kondo_nbodyg_particle/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 1 0 0
EOF

make_tj_base bad_tjn 9 no
cat > bad_tjn/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 0 0 0 0.1000000000000000 0.0000000000000000
EOF
cat > bad_tjn/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_kondo_base bad_kondon Kondo ncond
cat > bad_kondon/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 0 0 0 0.1000000000000000 0.0000000000000000
EOF
cat > bad_kondon/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_kondo_base bad_kondo_local_cross KondoGC no
cat > bad_kondo_local_cross/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 0 1 0 0.1000000000000000 0.0000000000000000
EOF
cat > bad_kondo_local_cross/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

cd accept_tj
run_hphi log_accept.txt "${hphi}" -e namelist.def
cd ..
check_nbodyg_output accept_tj
cd accept_tjgc
run_hphi log_accept.txt "${hphi}" -e namelist.def
cd ..
check_nbodyg_output accept_tjgc
cd accept_kondo
run_hphi log_accept.txt "${hphi}" -e namelist.def
cd ..
check_nbodyg_output accept_kondo
cd accept_kondogc
run_hphi log_accept.txt "${hphi}" -e namelist.def
cd ..
check_nbodyg_output accept_kondogc

expect_fail bad_tj_particle "does not conserve particle numbers"
expect_fail bad_kondo_nbodyg_particle "does not conserve particle numbers"
expect_fail bad_tjn "NBodyInterAll does not support tJNConserved"
expect_fail bad_kondon "NBodyInterAll does not support tJNConserved or KondoNConserved"
expect_fail bad_kondo_local_cross "Kondo local-spin NBodyInterAll factors require site_out == site_in"

echo "tJ/Kondo NBody validation accepts supported sectors and rejects unsupported N-conserved/local-spin cases."
