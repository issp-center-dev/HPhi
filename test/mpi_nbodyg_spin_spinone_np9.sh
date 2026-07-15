#!/bin/sh
set -e

testname="mpi_nbodyg_spin_spinone_np9"
hphi="../../../src/HPhi"
tol="0.00000001"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

rm -rf serial mpi
mkdir -p serial mpi

cd serial
cat > stan.in <<EOF
model = "Spin"
method = "Lanczos"
lattice = "chain"
L = 3
J = 1.0
2S = 2
2Sz = 0
Lanczos_max = 100
initial_iv = 1
EOF
run_hphi log_sdry.txt "${hphi}" -sdry stan.in
printf '    NBodyG  nbodyg.def\n' >> namelist.def
cat > nbodyg.def <<EOF
========================
NNBodyG 4
========================
========NBodyG==========
========================
1 2 2 2 2
2 1 2 1 0 2 0 2 2
2 1 0 1 2 2 2 2 0
3 0 2 0 2 1 2 1 0 2 0 2 2
EOF
run_hphi log_serial.txt "${hphi}" -e namelist.def
cp output/zvo_NBodyG.dat ../nbodyg_serial.dat
cd ..

cp serial/*.def mpi/
cd mpi
run_hphi log_mpi.txt ${MPIRUN} "${hphi}" -e namelist.def
cp output/zvo_NBodyG.dat ../nbodyg_mpi.dat
cd ..

awk -v t="${tol}" '
  function prefix(line, out) {
    out = line;
    sub(/[[:space:]]+[-+0-9.eE]+[[:space:]]+[-+0-9.eE]+$/, "", out);
    return out;
  }
  NR == FNR {
    s_prefix[NR] = prefix($0);
    s_re[NR] = $(NF - 1);
    s_im[NR] = $NF;
    n = NR;
    next;
  }
  {
    m = FNR;
    if (prefix($0) != s_prefix[m]) bad = 1;
    dr = $(NF - 1) - s_re[m];
    di = $NF - s_im[m];
    if (dr < 0) dr = -dr;
    if (di < 0) di = -di;
    if (dr >= t || di >= t) bad = 1;
  }
  END { exit (n == m && bad != 1) ? 0 : 1; }
' nbodyg_serial.dat nbodyg_mpi.dat || {
  echo "canonical Spin spin-one np=9 NBodyG output mismatch"
  paste nbodyg_serial.dat nbodyg_mpi.dat
  exit 1
}

echo "canonical Spin spin-one NBodyG MPI np=9 output matches serial output."
