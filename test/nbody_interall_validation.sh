#!/bin/sh
set -e

testname="nbody_interall_validation"
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
  printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  cd ..
}

rm -rf diag_im unpaired zero_product unsupported
make_spingc_base diag_im
make_spingc_base unpaired
make_spingc_base zero_product

cat > diag_im/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 1 0 1 0.0000000000000000 0.1000000000000000
EOF

cat > unpaired/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 1 0 0 1.0000000000000000 0.0000000000000000
EOF

cat > zero_product/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
2 0 1 0 0 0 1 0 0 1.0000000000000000 0.0000000000000000
EOF

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
printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 1 0 1 0.1000000000000000 0.0000000000000000
EOF
cd ..

expect_fail diag_im "finite imaginary"
expect_fail unpaired "adjacent Hermite"
expect_fail zero_product "zero same-site operator product"
expect_fail unsupported "supported only for SpinGC"

echo "NBodyInterAll validation rejects unsupported and ambiguous inputs."
