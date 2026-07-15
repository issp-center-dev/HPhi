#!/bin/sh
set -e

testname="mpi_nbodyg_spin"
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
L = 8
J = 1.0
2S = 1
2Sz = 0
Lanczos_max = 200
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
1 7 1 7 1
1 6 1 6 1
2 6 1 6 0 7 0 7 1
2 6 0 6 1 7 1 7 0
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
  echo "Serial/MPI canonical Spin NBodyG output mismatch"
  paste nbodyg_serial.dat nbodyg_mpi.dat
  exit 1
}

grep -q "INTER process site" mpi/log_mpi.txt || {
  echo "MPI run did not print an inter-process site summary"
  cat mpi/log_mpi.txt
  exit 1
}

awk -v t="${tol}" '
  $1 == 1 && $2 == 7 && $3 == 1 && $4 == 7 && $5 == 1 {
    re = $6; if (re < 0) re = -re;
    found = (re > t);
  }
  END { exit found ? 0 : 1; }
' nbodyg_mpi.dat || {
  echo "Inter-process diagonal canonical Spin NBodyG operator <n_{7 up}> was zero or missing"
  cat nbodyg_mpi.dat
  exit 1
}

echo "canonical Spin NBodyG MPI np=4 output matches serial output."
