#!/bin/sh
set -e

testname="nbody_spingc_spinone_validation"
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

make_spingc_spinone_base() {
  dir="$1"
  mkdir -p "${dir}"
  cd "${dir}"
  cat > stan.in <<EOF
model = "SpinGC"
method = "Lanczos"
lattice = "chain"
L = 4
J = 0.0
2S = 2
Lanczos_max = 50
initial_iv = 1
EOF
  run_hphi log_sdry.txt "${hphi}" -sdry stan.in
  printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  printf '    NBodyG  nbodyg.def\n' >> namelist.def
  cd ..
}

make_spin_canonical_general_base() {
  dir="$1"
  mkdir -p "${dir}"
  cd "${dir}"
  cat > stan.in <<EOF
model = "Spin"
method = "Lanczos"
lattice = "chain"
L = 4
J = 0.0
2S = 2
2Sz = 0
Lanczos_max = 50
initial_iv = 1
EOF
  run_hphi log_sdry.txt "${hphi}" -sdry stan.in
  printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  printf '    NBodyG  nbodyg.def\n' >> namelist.def
  cd ..
}

rm -rf bad_spin zero_interall bad_pair spin_general_accept \
  spin_general_bad_interall_sz spin_general_bad_nbodyg_sz nbodyg_bad_spin

make_spingc_spinone_base bad_spin
cat > bad_spin/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 3 0 3 0.1000000000000000 0.0000000000000000
EOF
cat > bad_spin/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_spingc_spinone_base zero_interall
cat > zero_interall/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
2 0 2 0 1 0 0 0 2 0.1000000000000000 0.0000000000000000
EOF
cat > zero_interall/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_spingc_spinone_base bad_pair
cat > bad_pair/nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
1 0 2 0 0 0.1000000000000000 0.0000000000000000
1 1 0 1 2 0.1000000000000000 0.0000000000000000
EOF
cat > bad_pair/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_spin_canonical_general_base spin_general_accept
cat > spin_general_accept/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 2 0 2 0.1000000000000000 0.0000000000000000
EOF
cat > spin_general_accept/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 2 0 2
EOF

make_spin_canonical_general_base spin_general_bad_interall_sz
cat > spin_general_bad_interall_sz/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 2 0 0 0.1000000000000000 0.0000000000000000
EOF
cat > spin_general_bad_interall_sz/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_spin_canonical_general_base spin_general_bad_nbodyg_sz
cat > spin_general_bad_nbodyg_sz/nbodyinterall.def <<EOF
========================
NNBodyInterAll 0
========================
========NBodyInterAll===
========================
EOF
cat > spin_general_bad_nbodyg_sz/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 2 0 0
EOF

make_spingc_spinone_base nbodyg_bad_spin
cat > nbodyg_bad_spin/nbodyinterall.def <<EOF
========================
NNBodyInterAll 0
========================
========NBodyInterAll===
========================
EOF
cat > nbodyg_bad_spin/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 3 0 3
EOF

expect_fail bad_spin "Spin index of NBodyInterAll is incorrect"
expect_fail zero_interall "zero same-site operator product"
expect_fail bad_pair "Hermite pair has inconsistent factors"
cd spin_general_accept
run_hphi log_accept.txt "${hphi}" -e namelist.def
cd ..
expect_fail spin_general_bad_interall_sz "does not conserve total Sz"
expect_fail spin_general_bad_nbodyg_sz "does not conserve total Sz"
expect_fail nbodyg_bad_spin "Spin index of NBodyG is incorrect"

echo "SpinGC/canonical Spin spin-one NBody validation rejects malformed and Sz-nonconserving inputs."
