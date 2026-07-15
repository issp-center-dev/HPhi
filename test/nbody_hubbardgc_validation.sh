#!/bin/sh
set -e

testname="nbody_hubbardgc_validation"
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

make_hubbardgc_base() {
  dir="$1"
  mkdir -p "${dir}"
  cd "${dir}"
  cat > stan.in <<EOF
model = "HubbardGC"
method = "Lanczos"
lattice = "chain"
L = 4
t = 0.0
U = 0.0
Lanczos_max = 50
initial_iv = 1
EOF
  run_hphi log_sdry.txt "${hphi}" -sdry stan.in
  printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  printf '    NBodyG  nbodyg.def\n' >> namelist.def
  cd ..
}

make_spingc_base() {
  dir="$1"
  mkdir -p "${dir}"
  cd "${dir}"
  cat > stan.in <<EOF
model = "SpinGC"
method = "Lanczos"
lattice = "chain"
L = 4
J = 0.0
2S = 1
Lanczos_max = 50
initial_iv = 1
EOF
  run_hphi log_sdry.txt "${hphi}" -sdry stan.in
  printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  printf '    NBodyG  nbodyg.def\n' >> namelist.def
  cd ..
}

rm -rf accept bad_spin_interall bad_spin_nbodyg bad_pair diag_im spin_cross_site

make_hubbardgc_base accept
cat > accept/nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
2 1 0 0 0 2 1 2 1 0.1000000000000000 0.0000000000000000
2 2 1 2 1 0 0 1 0 0.1000000000000000 0.0000000000000000
EOF
cat > accept/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 1 0 0 0
EOF

make_hubbardgc_base bad_spin_interall
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

make_hubbardgc_base bad_spin_nbodyg
cat > bad_spin_nbodyg/nbodyinterall.def <<EOF
========================
NNBodyInterAll 0
========================
========NBodyInterAll===
========================
EOF
cat > bad_spin_nbodyg/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 2 0 0
EOF

make_hubbardgc_base bad_pair
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

make_hubbardgc_base diag_im
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

make_spingc_base spin_cross_site
cat > spin_cross_site/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 1 1 0 1 0.1000000000000000 0.0000000000000000
EOF
cat > spin_cross_site/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 1 1 0 1
EOF

cd accept
run_hphi log_accept.txt "${hphi}" -e namelist.def
cd ..
expect_fail bad_spin_interall "Spin index of NBodyInterAll is incorrect"
expect_fail bad_spin_nbodyg "Spin index of NBodyG is incorrect"
expect_fail bad_pair "adjacent Hermite"
expect_fail diag_im "finite imaginary"
expect_fail spin_cross_site "requires site_out == site_in"

echo "HubbardGC NBody validation accepts cross-site fermion factors and rejects malformed inputs."
