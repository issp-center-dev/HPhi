#!/bin/sh
set -e

testname="mpi_anomalous_hubbardgc"
hphi="../../../src/HPhi"
tol="0.00000001"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

expect_mpi_fail() {
  dir="$1"
  pattern="$2"
  cd "${dir}"
  set +e
  ${MPIRUN} "${hphi}" -e namelist.def > log_mpi.txt 2>&1
  rc=$?
  set -e
  if [ "${rc}" -eq 0 ]; then
    echo "Expected ${dir} to fail, but it succeeded"
    exit 1
  fi
  grep -q "${pattern}" log_mpi.txt || {
    echo "Expected pattern '${pattern}' was not found in ${dir}/log_mpi.txt"
    cat log_mpi.txt
    exit 1
  }
  cd ..
}

compare_anomalousg() {
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
  ' "$@"
}

write_anomalous_term() {
  dir="$1"
  printf '    AnomalousTerm  anomalousterm.def\n' >> "${dir}/namelist.def"
  cat > "${dir}/anomalousterm.def" <<EOF
========================
NAnomalousTerm 4
========================
========AnomalousTerm===
========================
1 3 0 0 1  0.2100000000000000  0.1300000000000000
0 0 1 3 0  0.2100000000000000 -0.1300000000000000
1 3 0 3 1  0.0700000000000000 -0.0400000000000000
0 3 1 3 0  0.0700000000000000  0.0400000000000000
EOF
}

write_anomalousg() {
  dir="$1"
  printf '    AnomalousG     anomalousg.def\n' >> "${dir}/namelist.def"
  cat > "${dir}/anomalousg.def" <<EOF
========================
NAnomalousG 4
========================
========AnomalousG=======
========================
1 3 0 0 1
0 0 1 3 0
1 3 0 3 1
0 3 1 3 0
EOF
}

rm -rf serial mpi fulldiag_term fulldiag_g
mkdir -p serial mpi fulldiag_term fulldiag_g

cd serial
cat > stan.in <<EOF
model = "HubbardGC"
method = "Lanczos"
lattice = "chain"
L = 4
t = 1.0
U = 2.0
Lanczos_max = 120
initial_iv = 1
EOF
run_hphi log_sdry.txt "${hphi}" -sdry stan.in
cd ..
write_anomalous_term serial
write_anomalousg serial

cd serial
run_hphi log_serial.txt "${hphi}" -e namelist.def
cp output/zvo_energy.dat ../energy_serial.dat
cp output/zvo_AnomalousG.dat ../anomalousg_serial.dat
cd ..

cp serial/*.def mpi/
cd mpi
run_hphi log_mpi.txt ${MPIRUN} "${hphi}" -e namelist.def
cp output/zvo_energy.dat ../energy_mpi.dat
cp output/zvo_AnomalousG.dat ../anomalousg_mpi.dat
cd ..

diff=$(paste energy_serial.dat energy_mpi.dat \
  | awk '$1 == "Energy" && $3 == "Energy" {d=$2-$4; if(d<0)d=-d; if(d>m)m=d} END{printf "%.12g", m+0}')
awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d < t) ? 0 : 1}' || {
  echo "Serial/MPI HubbardGC AnomalousTerm energy mismatch: max diff ${diff}"
  paste energy_serial.dat energy_mpi.dat
  exit 1
}

compare_anomalousg anomalousg_serial.dat anomalousg_mpi.dat || {
  echo "Serial/MPI HubbardGC AnomalousG output mismatch"
  paste anomalousg_serial.dat anomalousg_mpi.dat
  exit 1
}

grep -q "INTER process site" mpi/log_mpi.txt || {
  echo "MPI run did not print an inter-process site summary"
  cat mpi/log_mpi.txt
  exit 1
}

awk -v t="${tol}" '
  $1 == 1 && $2 == 3 && $3 == 0 && $4 == 3 && $5 == 1 {
    re = $6; im = $7;
    if (re < 0) re = -re;
    if (im < 0) im = -im;
    found = (re > t || im > t);
  }
  END { exit found ? 0 : 1; }
' anomalousg_mpi.dat || {
  echo "Two-rank-bit HubbardGC AnomalousG operator was zero or missing"
  cat anomalousg_mpi.dat
  exit 1
}

cp serial/*.def fulldiag_term/
sed -e 's/^CalcType.*/CalcType   2/' fulldiag_term/calcmod.def > fulldiag_term/calcmod.tmp
mv fulldiag_term/calcmod.tmp fulldiag_term/calcmod.def

cp serial/*.def fulldiag_g/
rm -f fulldiag_g/anomalousterm.def
grep -v 'AnomalousTerm' fulldiag_g/namelist.def > fulldiag_g/namelist.tmp
mv fulldiag_g/namelist.tmp fulldiag_g/namelist.def
sed -e 's/^CalcType.*/CalcType   2/' fulldiag_g/calcmod.def > fulldiag_g/calcmod.tmp
mv fulldiag_g/calcmod.tmp fulldiag_g/calcmod.def

expect_mpi_fail fulldiag_term "AnomalousTerm does not support MPI FullDiag"
expect_mpi_fail fulldiag_g "AnomalousG does not support MPI FullDiag"

echo "HubbardGC AnomalousTerm/AnomalousG MPI np=4 matches serial and rejects MPI FullDiag."
