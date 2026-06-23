#!/bin/sh
set -e

testname="lanczos_hubbardgc_anomalous_pair"
hphi="../../src/HPhi"
tol="0.00000001"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

check_energy_file() {
  file="$1"
  awk -v t="${tol}" '
    /^Energy/ {
      d = $2 + 0.3;
      if (d < 0) d = -d;
      ok = (d < t);
    }
    END { exit ok ? 0 : 1; }
  ' "${file}" || {
    echo "Unexpected Lanczos anomalous-pair ground-state energy"
    cat "${file}"
    exit 1
  }
}

check_phys_file() {
  file="$1"
  awk -v t="${tol}" '
    NR == 2 {
      d = $1 + 0.3;
      if (d < 0) d = -d;
      ok = (d < t);
    }
    END { exit ok ? 0 : 1; }
  ' "${file}" || {
    echo "Unexpected FullDiag anomalous-pair ground-state energy"
    cat "${file}"
    exit 1
  }
}

check_anomalousg_file() {
  file="$1"
  awk -v t="${tol}" '
    function abs(x) { return x < 0 ? -x : x }
    NR == 1 {
      ok1 = ($1 == 0 && $2 == 0 && $3 == 0 && $4 == 0 && $5 == 1 &&
             abs($6 + 0.5) < t && abs($7) < t)
    }
    NR == 2 {
      ok2 = ($1 == 1 && $2 == 0 && $3 == 1 && $4 == 0 && $5 == 0 &&
             abs($6 + 0.5) < t && abs($7) < t)
    }
    END { exit (NR == 2 && ok1 && ok2) ? 0 : 1; }
  ' "${file}" || {
    echo "Unexpected AnomalousG output"
    cat "${file}"
    exit 1
  }
}

cat > stan.in <<EOF
model = "HubbardGC"
method = "Lanczos"
lattice = "chain"
L = 1
t = 0.0
U = 0.0
Lanczos_max = 50
initial_iv = 1
EOF

rm -rf output
run_hphi log_sdry.txt "${hphi}" -sdry stan.in
printf '    AnomalousTerm  anomalousterm.def\n' >> namelist.def
printf '    AnomalousG     anomalousg.def\n' >> namelist.def
sed -e 's/^OutputEigenVec.*/OutputEigenVec   1/' calcmod.def > calcmod.outvec
mv calcmod.outvec calcmod.def

cat > anomalousterm.def <<EOF
========================
NAnomalousTerm 2
========================
========AnomalousTerm===
========================
0 0 0 0 1 0.3000000000000000 0.0000000000000000
1 0 1 0 0 0.3000000000000000 0.0000000000000000
EOF

cat > anomalousg.def <<EOF
========================
NAnomalousG 2
========================
========AnomalousG=======
========================
0 0 0 0 1
1 0 1 0 0
EOF

# ED reference for H = 0.3 * (c_up c_down + c_down^dagger c_up^dagger):
# <0|H|up down> = -0.3, so E0 = -0.3 and both anomalous expectations are -0.5.
run_hphi log_lanczos.txt "${hphi}" -e namelist.def
test -f output/zvo_AnomalousG.dat || { echo "zvo_AnomalousG.dat was not generated"; exit 1; }
check_energy_file output/zvo_energy.dat
check_anomalousg_file output/zvo_AnomalousG.dat

cat > teonebody.def <<EOF
========================
NTimeSteps    2
========================
=========  OneBody Time Evolution  ==========
========================
0.00  0
0.01  0
EOF
grep -v 'AnomalousTerm' namelist.def > namelist_te.def
printf '       TEOneBody  teonebody.def\n' >> namelist_te.def
cp calcmod.def calcmod.lanczos
sed -e 's/^CalcType.*/CalcType   4/' \
    -e 's/^InputEigenVec.*/InputEigenVec   1/' \
    -e 's/^OutputEigenVec.*/OutputEigenVec   0/' \
    calcmod.lanczos > calcmod.def
cp modpara.def modpara.lanczos
sed -e 's/^Lanczos_max.*/Lanczos_max    2/' modpara.lanczos > modpara.def
printf 'ExpandCoef     10\n' >> modpara.def
run_hphi log_te.txt "${hphi}" -e namelist_te.def
test -f output/zvo_AnomalousG_step0.dat || { echo "zvo_AnomalousG_step0.dat was not generated"; exit 1; }
check_anomalousg_file output/zvo_AnomalousG_step0.dat
mv modpara.lanczos modpara.def

sed -e 's/^CalcType.*/CalcType   2/' \
    -e 's/^OutputEigenVec.*/OutputEigenVec   0/' \
    -e 's/^OutputHam.*/OutputHam   0/' \
    calcmod.lanczos > calcmod.fulldiag
mv calcmod.fulldiag calcmod.def
rm -rf output
run_hphi log_fulldiag.txt "${hphi}" -e namelist.def
test -f output/zvo_AnomalousG_eigen0.dat || { echo "zvo_AnomalousG_eigen0.dat was not generated"; exit 1; }
check_phys_file output/zvo_phys.dat
check_anomalousg_file output/zvo_AnomalousG_eigen0.dat

sed -e 's/^OutputHam.*/OutputHam   1/' calcmod.def > calcmod.outputham
mv calcmod.outputham calcmod.def
rm -rf output
run_hphi log_outputham.txt "${hphi}" -e namelist.def
test -f output/zvo_Ham.dat || { echo "zvo_Ham.dat was not generated"; exit 1; }
awk -v t="${tol}" '
  NR > 2 && $1 != $2 {
    d = $3 + 0.3;
    if (d < 0) d = -d;
    ok = (d < t && $4 < t && $4 > -t);
  }
  END { exit ok ? 0 : 1; }
' output/zvo_Ham.dat || {
  echo "FullDiag Hamiltonian output has no anomalous off-diagonal entry"
  cat output/zvo_Ham.dat
  exit 1
}

echo "HubbardGC AnomalousTerm/AnomalousG Lanczos and FullDiag checks passed."
