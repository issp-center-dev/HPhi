#!/bin/sh
set -e

testname="lanczos_spin_nbody_interall"
hphi="../../src/HPhi"
tol="0.000001"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

cat > stan.in <<EOF
model = "Spin"
method = "Lanczos"
lattice = "chain"
L = 8
J = 1.0
2S = 1
2Sz = 0
Lanczos_max = 200
exct = 1
LanczosTarget = 0
initial_iv = 1
EOF

rm -rf output
run_hphi log_sdry.txt "${hphi}" -sdry stan.in
printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
printf '    NBodyG  nbodyg.def\n' >> namelist.def

cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 3
========================
========NBodyInterAll===
========================
1 2 1 2 1 0.2000000000000000 0.0000000000000000
2 0 1 0 0 1 0 1 1 0.3000000000000000 0.0000000000000000
2 0 0 0 1 1 1 1 0 0.3000000000000000 0.0000000000000000
EOF

cat > nbodyg.def <<EOF
========================
NNBodyG 2
========================
========NBodyG==========
========================
1 2 1 2 1
2 0 1 0 0 1 0 1 1
EOF

run_hphi log_lanczos.txt "${hphi}" -e namelist.def
e_lanczos=$(awk '/^Energy/{print $2; exit}' output/zvo_energy.dat)
test -f output/zvo_NBodyG.dat || { echo "zvo_NBodyG.dat was not generated"; exit 1; }
awk 'NF >= 7 {count++} END{exit count == 2 ? 0 : 1}' output/zvo_NBodyG.dat || {
  echo "Unexpected NBodyG output"
  cat output/zvo_NBodyG.dat
  exit 1
}

sed -e 's/^CalcType.*/CalcType   2/' -e 's/^OutputHam.*/OutputHam   0/' calcmod.def > calcmod.fulldiag
mv calcmod.def calcmod.lanczos
mv calcmod.fulldiag calcmod.def
rm -rf output
run_hphi log_fulldiag.txt "${hphi}" -e namelist.def
if [ -f output/zvo_phys.dat ]; then
  e_fulldiag=$(awk 'NR==2{print $1; exit}' output/zvo_phys.dat)
else
  e_fulldiag=$(awk 'NR==1{print $2; exit}' output/Eigenvalue.dat)
fi

diff=$(awk -v a="${e_lanczos}" -v b="${e_fulldiag}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%.12g", d}')
awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d < t) ? 0 : 1}' || {
  echo "canonical Spin NBodyInterAll Lanczos/FullDiag energy mismatch: ${e_lanczos} vs ${e_fulldiag}"
  exit 1
}

sed -e 's/^OutputHam.*/OutputHam   1/' calcmod.def > calcmod.outputham
mv calcmod.outputham calcmod.def
rm -rf output
run_hphi log_outputham.txt "${hphi}" -e namelist.def

test -f output/zvo_Ham.dat || { echo "zvo_Ham.dat was not generated"; exit 1; }
awk 'NR>2 && $1 != $2 {found=1} END{exit found ? 0 : 1}' output/zvo_Ham.dat || {
  echo "canonical Spin NBodyInterAll produced no off-diagonal Hamiltonian entry"
  cat output/zvo_Ham.dat
  exit 1
}
awk 'NR>2 && $1 == $2 && $3 == "0.200000" {found=1} END{exit found ? 0 : 1}' output/zvo_Ham.dat || {
  echo "canonical Spin diagonal NBodyInterAll term was not folded into the diagonal"
  cat output/zvo_Ham.dat
  exit 1
}

echo "canonical Spin NBodyInterAll/NBodyG Lanczos, FullDiag, and OutputHam checks passed."
