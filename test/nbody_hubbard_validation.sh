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

make_hubbard_base() {
  dir="$1"
  mkdir -p "${dir}"
  cd "${dir}"
  cat > stan.in <<EOF
model = "Hubbard"
method = "FullDiag"
lattice = "chain"
L = 4
t = 1.0
U = 0.0
nelec = 4
2Sz = 0
Lanczos_max = 50
initial_iv = 1
EOF
  run_hphi log_sdry.txt "${hphi}" -sdry stan.in
  printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  printf '    NBodyG  nbodyg.def\n' >> namelist.def
  cd ..
}

rm -rf accept bad_interall_particle bad_nbodyg_particle bad_spin_interall bad_pair diag_im

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

cd accept
run_hphi log_accept.txt "${hphi}" -e namelist.def
cd ..

expect_fail bad_interall_particle "does not conserve particle numbers"
expect_fail bad_nbodyg_particle "does not conserve particle numbers"
expect_fail bad_spin_interall "Spin index of NBodyInterAll is incorrect"
expect_fail bad_pair "Off-diagonal NBodyInterAll terms must appear as adjacent Hermite pairs"
expect_fail diag_im "Diagonal NBodyInterAll term has a finite imaginary part"

echo "canonical Hubbard NBody validation rejects malformed and nonconserving inputs."
