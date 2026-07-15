#!/bin/sh
set -e

testname="lanczos_spingc_nbody_interall"
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
model = "SpinGC"
method = "Lanczos"
lattice = "chain"
L = 4
J = 0.0
2S = 1
Lanczos_max = 200
initial_iv = 1
EOF

rm -rf output
run_hphi log_sdry.txt "${hphi}" -sdry stan.in
printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def

cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 4
========================
========NBodyInterAll===
========================
3 0 1 0 1 1 0 1 0 2 1 2 1 0.3000000000000000 0.0000000000000000
2 0 1 0 0 0 0 0 1 0.2000000000000000 0.0000000000000000
3 0 1 0 0 1 0 1 1 2 1 2 1 0.4000000000000000 0.2000000000000000
3 0 0 0 1 1 1 1 0 2 1 2 1 0.4000000000000000 -0.2000000000000000
EOF

run_hphi log_lanczos.txt "${hphi}" -e namelist.def
e_lanczos=$(awk '/^Energy/{print $2; exit}' output/zvo_energy.dat)

sed -e 's/^CalcType.*/CalcType   2/' -e 's/^OutputHam.*/OutputHam   0/' calcmod.def > calcmod.fulldiag
mv calcmod.def calcmod.lanczos
mv calcmod.fulldiag calcmod.def
rm -rf output
run_hphi log_fulldiag.txt "${hphi}" -e namelist.def
e_fulldiag=$(awk 'NR==2{print $1; exit}' output/zvo_phys.dat)

diff=$(awk -v a="${e_lanczos}" -v b="${e_fulldiag}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%.12g", d}')
awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d < t) ? 0 : 1}' || {
  echo "NBodyInterAll Lanczos/FullDiag energy mismatch: ${e_lanczos} vs ${e_fulldiag}"
  exit 1
}

sed -e 's/^OutputHam.*/OutputHam   1/' calcmod.def > calcmod.outputham
mv calcmod.outputham calcmod.def
rm -rf output
run_hphi log_outputham.txt "${hphi}" -e namelist.def

test -f output/zvo_Ham.dat || { echo "zvo_Ham.dat was not generated"; exit 1; }
awk 'NR>2 && $1 != $2 {found=1} END{exit found ? 0 : 1}' output/zvo_Ham.dat || {
  echo "FullDiag Hamiltonian output has no off-diagonal NBodyInterAll entry"
  cat output/zvo_Ham.dat
  exit 1
}
awk 'NR>2 && $1 == $2 && $3 == "0.200000" {found=1} END{exit found ? 0 : 1}' output/zvo_Ham.dat || {
  echo "Net-diagonal same-site NBodyInterAll term was not folded into the diagonal"
  cat output/zvo_Ham.dat
  exit 1
}

echo "SpinGC NBodyInterAll Lanczos, FullDiag, and OutputHam checks passed."
