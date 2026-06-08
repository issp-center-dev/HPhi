#!/bin/sh
set -e

testname="nbodyg_validation"
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
  printf '    NBodyG  nbodyg.def\n' >> namelist.def
  cd ..
}

rm -rf bad_spin bad_site cross_site too_few unsupported
make_spingc_base bad_spin
make_spingc_base bad_site
make_spingc_base cross_site
make_spingc_base too_few

# Spin index outside {0,1}.
cat > bad_spin/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 2 0 2
EOF

# Site index outside [0, Nsite).
cat > bad_site/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 9 1 9 1
EOF

# site_out != site_in is not allowed for a single factor.
cat > cross_site/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 1 1 1
EOF

# A factor must carry 4 integer fields (site_out spin_out site_in spin_in).
cat > too_few/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 1 0
EOF

# NBodyG is currently restricted to spin-1/2 SpinGC; canonical Spin is rejected.
mkdir -p unsupported
cd unsupported
cat > stan.in <<EOF
model = "Spin"
method = "Lanczos"
lattice = "chain"
L = 4
J = 0.0
2S = 1
2Sz = 0
Lanczos_max = 50
initial_iv = 1
EOF
run_hphi log_sdry.txt "${hphi}" -sdry stan.in
printf '    NBodyG  nbodyg.def\n' >> namelist.def
cat > nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 1 0 1
EOF
cd ..

expect_fail bad_spin "Spin index of NBodyG is incorrect"
expect_fail bad_site "Site index of NBodyG is incorrect"
expect_fail cross_site "requires site_out == site_in"
expect_fail too_few "too few integer fields"
expect_fail unsupported "supported only for SpinGC"

echo "NBodyG validation rejects unsupported and malformed inputs."
