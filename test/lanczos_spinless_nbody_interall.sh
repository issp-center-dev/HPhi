#!/bin/sh
set -e

testname="lanczos_spinless_nbody_interall"
hphi="../../../src/HPhi"
SRCDIR="$1"
tol="0.00000001"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

make_base() {
  dir="$1"
  model="$2"
  mkdir -p "${dir}"
  cd "${dir}"
  if [ "${model}" = "SpinlessFermion" ]; then
    python3 "${SRCDIR}/test/testSpinlessCalc.py" -p /bin/true \
      -m "${model}" -s 6 -n 3 --hopping 0.0 > log_generate.txt 2>&1
  else
    python3 "${SRCDIR}/test/testSpinlessCalc.py" -p /bin/true \
      -m "${model}" -s 6 --hopping 0.0 > log_generate.txt 2>&1
  fi
  cp namelist.def namelist.base
  cd ..
}

run_case() {
  tag="$1"
  model="$2"

  rm -rf "${tag}_nbody" "${tag}_legacy"
  make_base "${tag}_nbody" "${model}"
  make_base "${tag}_legacy" "${model}"

  cd "${tag}_nbody"
  cp namelist.base namelist.def
  printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 3
========================
========NBodyInterAll===
========================
1 0 0 1 0 0.3700000000000000 0.1100000000000000
1 1 0 0 0 0.3700000000000000 -0.1100000000000000
2 2 0 2 0 4 0 4 0 0.2300000000000000 0.0000000000000000
EOF
  run_hphi log_nbody.txt "${hphi}" -e namelist.def
  cp output/zvo_energy.dat "../${tag}_energy_nbody.dat"
  cd ..

  cd "${tag}_legacy"
  cp namelist.base namelist.def
  cat > trans.def <<EOF
========================
NTransfer      2
========================
========i_j_s_tijs======
========================
0 0 1 0 0.3700000000000000 0.1100000000000000
1 0 0 0 0.3700000000000000 -0.1100000000000000
EOF
  printf '    CoulombInter  coulombinter_eq.def\n' >> namelist.def
  cat > coulombinter_eq.def <<EOF
=============================================
NCoulombInter          1
=============================================
================CoulombInter=================
=============================================
2 4 0.2300000000000000
EOF
  run_hphi log_legacy.txt "${hphi}" -e namelist.def
  cp output/zvo_energy.dat "../${tag}_energy_legacy.dat"
  cd ..

  diff=$(paste "${tag}_energy_nbody.dat" "${tag}_energy_legacy.dat" \
    | awk '$1 == "Energy" && $3 == "Energy" {d=$2-$4; if(d<0)d=-d; if(d>m)m=d} END{printf "%.12g", m+0}')
  awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d < t) ? 0 : 1}' || {
    echo "${model} NBodyInterAll energy differs from legacy Transfer/CoulombInter: max diff ${diff}"
    paste "${tag}_energy_nbody.dat" "${tag}_energy_legacy.dat"
    exit 1
  }
}

run_case spinless SpinlessFermion
run_case spinlessgc SpinlessFermionGC

echo "Spinless NBodyInterAll N=1 hopping and N=2 density terms match legacy operators."
