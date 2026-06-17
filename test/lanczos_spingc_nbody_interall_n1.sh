#!/bin/sh
set -e

# Exercise the N=1 NBodyInterAll path end to end (a single-factor term).
# The Hamiltonian is, on a J=0 SpinGC L=4 chain,
#   H = 0.2 n_{3 up} + 0.3 (S^+_0 + S^-_0) = 0.2 n_{3 up} + 0.3 sigma^x_0 .
# Sites 1 and 2 carry no term. Site 3 contributes 0 (down) or 0.2 (up); site 0
# is a free spin in a transverse field with eigenvalues +/-0.3. Hence the exact
# ground-state energy is -0.3, independent of sites 1, 2, 3 (their minimum is 0).

testname="lanczos_spingc_nbody_interall_n1"
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

# Net-diagonal N=1 field on site 3, plus an N=1 Hermitian off-diagonal pair
# (S^+ then its conjugate S^- on the immediately following line) on site 0.
cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 3
========================
========NBodyInterAll===
========================
1 3 1 3 1 0.2000000000000000 0.0000000000000000
1 0 1 0 0 0.3000000000000000 0.0000000000000000
1 0 0 0 1 0.3000000000000000 0.0000000000000000
EOF

run_hphi log_lanczos.txt "${hphi}" -e namelist.def
e_lanczos=$(awk '/^Energy/{print $2; exit}' output/zvo_energy.dat)

# Exact ground-state energy for this single-factor Hamiltonian is -0.3.
awk -v a="${e_lanczos}" -v t="${tol}" \
  'BEGIN{d=a-(-0.3); if(d<0)d=-d; exit (d < t) ? 0 : 1}' || {
  echo "N=1 NBodyInterAll ground-state energy is not -0.3: ${e_lanczos}"
  exit 1
}

sed -e 's/^CalcType.*/CalcType   2/' -e 's/^OutputHam.*/OutputHam   0/' calcmod.def > calcmod.fulldiag
mv calcmod.def calcmod.lanczos
mv calcmod.fulldiag calcmod.def
rm -rf output
run_hphi log_fulldiag.txt "${hphi}" -e namelist.def
e_fulldiag=$(awk 'NR==2{print $1; exit}' output/zvo_phys.dat)

diff=$(awk -v a="${e_lanczos}" -v b="${e_fulldiag}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%.12g", d}')
awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d < t) ? 0 : 1}' || {
  echo "N=1 NBodyInterAll Lanczos/FullDiag energy mismatch: ${e_lanczos} vs ${e_fulldiag}"
  exit 1
}

sed -e 's/^OutputHam.*/OutputHam   1/' calcmod.def > calcmod.outputham
mv calcmod.outputham calcmod.def
rm -rf output
run_hphi log_outputham.txt "${hphi}" -e namelist.def

test -f output/zvo_Ham.dat || { echo "zvo_Ham.dat was not generated"; exit 1; }
# The N=1 S^+/S^- pair must produce off-diagonal matrix elements.
awk 'NR>2 && $1 != $2 {found=1} END{exit found ? 0 : 1}' output/zvo_Ham.dat || {
  echo "N=1 off-diagonal NBodyInterAll term produced no off-diagonal entry"
  cat output/zvo_Ham.dat
  exit 1
}
# The N=1 net-diagonal field must be folded into the diagonal as 0.2.
awk 'NR>2 && $1 == $2 && $3 == "0.200000" {found=1} END{exit found ? 0 : 1}' output/zvo_Ham.dat || {
  echo "N=1 net-diagonal NBodyInterAll field was not folded into the diagonal"
  cat output/zvo_Ham.dat
  exit 1
}

echo "SpinGC N=1 NBodyInterAll energy, Lanczos/FullDiag, and OutputHam checks passed."
