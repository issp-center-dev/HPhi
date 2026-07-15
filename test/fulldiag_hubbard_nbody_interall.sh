#!/bin/sh
set -e

testname="fulldiag_hubbard_nbody_interall"
hphi="../../src/HPhi"
tol="0.000001"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

compare_scalar() {
  label="$1"
  left="$2"
  right="$3"
  diff=$(awk -v a="${left}" -v b="${right}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%.12g", d}')
  awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d < t) ? 0 : 1}' || {
    echo "${label}: ${left} vs ${right}"
    exit 1
  }
}

cat > stan.in <<EOF
model = "Hubbard"
method = "Lanczos"
lattice = "chain"
L = 4
t = 1.0
U = 0.0
nelec = 4
2Sz = 0
Lanczos_max = 120
initial_iv = 1
EOF

rm -rf output
run_hphi log_sdry.txt "${hphi}" -sdry stan.in
printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def

cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 3
========================
========NBodyInterAll===
========================
1 0 0 0 0 0.1250000000000000 0.0000000000000000
2 3 0 0 0 1 1 1 1 0.2500000000000000 0.0000000000000000
2 1 1 1 1 0 0 3 0 0.2500000000000000 0.0000000000000000
EOF

run_hphi log_lanczos.txt "${hphi}" -e namelist.def
e_lanczos=$(awk '/^Energy/{print $2; exit}' output/zvo_energy.dat)

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

compare_scalar "Hubbard NBodyInterAll Lanczos/FullDiag energy mismatch" "${e_lanczos}" "${e_fulldiag}"

sed -e 's/^OutputHam.*/OutputHam   1/' calcmod.def > calcmod.outputham
mv calcmod.outputham calcmod.def
rm -rf output
run_hphi log_outputham.txt "${hphi}" -e namelist.def

test -f output/zvo_Ham.dat || { echo "zvo_Ham.dat was not generated"; exit 1; }
awk 'NR>2 && $1 != $2 {found=1} END{exit found ? 0 : 1}' output/zvo_Ham.dat || {
  echo "Hubbard NBodyInterAll produced no off-diagonal Hamiltonian entry"
  cat output/zvo_Ham.dat
  exit 1
}
awk 'NR>2 && $1 == $2 {
  d = $3 - 0.125000;
  if (d < 0) d = -d;
  if (d < 0.000001) found=1;
} END{exit found ? 0 : 1}' output/zvo_Ham.dat || {
  echo "Hubbard diagonal NBodyInterAll term was not folded into the diagonal"
  cat output/zvo_Ham.dat
  exit 1
}

cat > stan_ncond.in <<EOF
model = "Hubbard"
method = "Lanczos"
lattice = "chain"
L = 4
t = 1.0
U = 0.0
ncond = 4
Lanczos_max = 120
initial_iv = 1
EOF

rm -rf output
run_hphi log_ncond_sdry.txt "${hphi}" -sdry stan_ncond.in
printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
2 0 0 0 1 1 1 1 0 0.1700000000000000 0.0300000000000000
2 1 0 1 1 0 1 0 0 0.1700000000000000 -0.0300000000000000
EOF

run_hphi log_ncond_lanczos.txt "${hphi}" -e namelist.def
e_ncond_lanczos=$(awk '/^Energy/{print $2; exit}' output/zvo_energy.dat)

sed -e 's/^CalcType.*/CalcType   2/' -e 's/^OutputHam.*/OutputHam   0/' calcmod.def > calcmod.ncond.fulldiag
mv calcmod.def calcmod.ncond.lanczos
mv calcmod.ncond.fulldiag calcmod.def
rm -rf output
run_hphi log_ncond_fulldiag.txt "${hphi}" -e namelist.def
if [ -f output/zvo_phys.dat ]; then
  e_ncond_fulldiag=$(awk 'NR==2{print $1; exit}' output/zvo_phys.dat)
else
  e_ncond_fulldiag=$(awk 'NR==1{print $2; exit}' output/Eigenvalue.dat)
fi

compare_scalar "HubbardNConserved NBodyInterAll Lanczos/FullDiag energy mismatch" \
  "${e_ncond_lanczos}" "${e_ncond_fulldiag}"

echo "Hubbard and HubbardNConserved NBodyInterAll Lanczos, FullDiag, and OutputHam checks passed."
