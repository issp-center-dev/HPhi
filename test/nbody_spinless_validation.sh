#!/bin/sh
set -e

testname="nbody_spinless_validation"
hphi="../../../src/HPhi"
SRCDIR="$1"

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

make_spinless_base() {
  dir="$1"
  mkdir -p "${dir}"
  cd "${dir}"
  python3 "${SRCDIR}/test/testSpinlessCalc.py" -p /bin/true \
    -m SpinlessFermion -s 8 -n 3 -V 0.5 > log_generate.txt 2>&1
  printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  printf '    NBodyG  nbodyg.def\n' >> namelist.def
  cd ..
}

rm -rf accept bad_spin_interall bad_spin_nbodyg bad_pair diag_im bad_fulldiag

make_spinless_base accept
cat > accept/nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
1 0 0 1 0 0.1000000000000000 0.0000000000000000
1 1 0 0 0 0.1000000000000000 0.0000000000000000
EOF
cat > accept/nbodyg.def <<EOF
========================
NNBodyG 1
========================
========NBodyG==========
========================
1 0 0 1 0
EOF

make_spinless_base bad_spin_interall
cat > bad_spin_interall/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 1 0 0 0.1000000000000000 0.0000000000000000
EOF
cat > bad_spin_interall/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_spinless_base bad_spin_nbodyg
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
1 0 1 0 0
EOF

make_spinless_base bad_pair
cat > bad_pair/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 0 1 0 0.1000000000000000 0.0000000000000000
EOF
cat > bad_pair/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

make_spinless_base diag_im
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

make_spinless_base bad_fulldiag
awk '{if ($1 == "CalcType") print "CalcType   2"; else print}' \
  bad_fulldiag/calcmod.def > bad_fulldiag/calcmod.tmp
mv bad_fulldiag/calcmod.tmp bad_fulldiag/calcmod.def
cat > bad_fulldiag/nbodyinterall.def <<EOF
========================
NNBodyInterAll 1
========================
========NBodyInterAll===
========================
1 0 0 0 0 0.1000000000000000 0.0000000000000000
EOF
cat > bad_fulldiag/nbodyg.def <<EOF
========================
NNBodyG 0
========================
========NBodyG==========
========================
EOF

cd accept
run_hphi log_accept.txt "${hphi}" -e namelist.def
cd ..

expect_fail bad_spin_interall "Spin index of NBodyInterAll is incorrect"
expect_fail bad_spin_nbodyg "Spin index of NBodyG is incorrect"
expect_fail bad_pair "Off-diagonal NBodyInterAll terms must appear as adjacent Hermite pairs"
expect_fail diag_im "Diagonal NBodyInterAll term has a finite imaginary part"
expect_fail bad_fulldiag "NBodyInterAll is not yet supported in FullDiag for SpinlessFermion"

echo "Spinless NBody validation accepts cross-site fermion factors and rejects malformed inputs."
