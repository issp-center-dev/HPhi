#!/bin/sh
set -e

testname="nbody_hubbard_validation"
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

set_calcspec() {
  dir="$1"
  value="$2"
  awk -v value="${value}" '
    $1 == "CalcSpec" {
      printf "CalcSpec        %s\n", value
      next
    }
    { print }
  ' "${dir}/calcmod.def" > "${dir}/calcmod.def.tmp"
  mv "${dir}/calcmod.def.tmp" "${dir}/calcmod.def"
}

make_hubbard_base() {
  dir="$1"
  sector="${2:-fixed}"
  mkdir -p "${dir}"
  cd "${dir}"
  cat > stan.in <<EOF
model = "Hubbard"
method = "FullDiag"
lattice = "chain"
L = 4
t = 1.0
U = 0.0
Lanczos_max = 50
initial_iv = 1
EOF
  if [ "${sector}" = "ncond" ]; then
    printf "ncond = 4\n" >> stan.in
  else
    printf "nelec = 4\n2Sz = 0\n" >> stan.in
  fi
  run_hphi log_sdry.txt "${hphi}" -sdry stan.in
  printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  printf '    NBodyG  nbodyg.def\n' >> namelist.def
  cd ..
}

rm -rf accept accept_ncond bad_interall_particle bad_nbodyg_particle \
  bad_spin_interall bad_pair diag_im bad_ncond_calcspec_interall \
  bad_ncond_calcspec_nbodyg

make_hubbard_base accept
cat > accept/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 0 0 0 0.1000000000000000 0.0000000000000000
EOF
cat > accept/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 1 0 0 0
EOF

make_hubbard_base accept_ncond ncond
cat > accept_ncond/nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
1 0 0 0 1 0.0500000000000000 0.0000000000000000
1 0 1 0 0 0.0500000000000000 0.0000000000000000
EOF
cat > accept_ncond/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 0 0 1
EOF

make_hubbard_base bad_interall_particle
cat > bad_interall_particle/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 1 0 0 1 0.1000000000000000 0.0000000000000000
EOF
cat > bad_interall_particle/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_hubbard_base bad_nbodyg_particle
cat > bad_nbodyg_particle/nbodyinterall.def <<EOF
========================
NNBodyInterAll 0
========================
========NBodyInterAll===
========================
EOF
cat > bad_nbodyg_particle/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 1 0 0 1
EOF

make_hubbard_base bad_spin_interall
cat > bad_spin_interall/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 2 0 0 0.1000000000000000 0.0000000000000000
EOF
cat > bad_spin_interall/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_hubbard_base bad_pair
cat > bad_pair/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 1 0 0 0 0.1000000000000000 0.0000000000000000
EOF
cat > bad_pair/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_hubbard_base diag_im
cat > diag_im/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 0 0 0 0.1000000000000000 0.1000000000000000
EOF
cat > diag_im/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_hubbard_base bad_ncond_calcspec_interall ncond
set_calcspec bad_ncond_calcspec_interall 1
cat > bad_ncond_calcspec_interall/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 0 0 0 0.1000000000000000 0.0000000000000000
EOF
cat > bad_ncond_calcspec_interall/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_hubbard_base bad_ncond_calcspec_nbodyg ncond
set_calcspec bad_ncond_calcspec_nbodyg 1
cat > bad_ncond_calcspec_nbodyg/nbodyinterall.def <<EOF
========================
NNBodyInterAll 0
========================
========NBodyInterAll===
========================
EOF
cat > bad_ncond_calcspec_nbodyg/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 0 0 0
EOF

cd accept
run_hphi log_accept.txt "${hphi}" -e namelist.def
cd ..
cd accept_ncond
run_hphi log_accept_ncond.txt "${hphi}" -e namelist.def
cd ..
check_nbodyg_output accept_ncond

expect_fail bad_interall_particle "does not conserve particle numbers"
expect_fail bad_nbodyg_particle "does not conserve particle numbers"
expect_fail bad_spin_interall "Spin index of NBodyInterAll is incorrect"
expect_fail bad_pair "Off-diagonal NBodyInterAll terms must appear as adjacent Hermite pairs"
expect_fail diag_im "Diagonal NBodyInterAll term has a finite imaginary part"
expect_fail bad_ncond_calcspec_interall "NBodyInterAll does not support HubbardNConserved with CalcSpec"
expect_fail bad_ncond_calcspec_nbodyg "NBodyG does not support HubbardNConserved with CalcSpec"

echo "canonical and NConserved Hubbard NBody validation rejects malformed and unsupported inputs."
