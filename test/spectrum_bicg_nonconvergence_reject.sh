#!/bin/sh -e

# Regression test for BiCG spectrum non-convergence handling.
# A deliberately tiny Lanczos_max must fail clearly instead of writing a
# near-zero DynamicalGreen file and returning success.

HPHI=../../src/HPhi

rm -rf spectrum_bicg_nonconvergence_reject/
mkdir -p spectrum_bicg_nonconvergence_reject/
cd spectrum_bicg_nonconvergence_reject

cat > stan_gs.in <<EOF
model = "Hubbard"
method = "CG"
lattice = "chain"
L = 4
t = 1.0
U = 4.0
nelec = 4
2Sz = 0
EigenVecIO = "out"
EOF

${HPHI} -s stan_gs.in > gs.log 2>&1

cat > calcmod_cg.def <<EOF
CalcType   3
CalcModel   0
ReStart   0
CalcSpec   1
CalcEigenVec   0
InitialVecType   0
InputEigenVec   1
OutputEigenVec   0
InputHam   0
OutputHam   0
OutputExVec   0
EOF

cat > SpectrumModpara <<EOF
--------------------
Model_Parameters   0
--------------------
HPhi_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
Nsite          4
2Sz            0
Ncond          4
Lanczos_max    1
initial_iv     -1
exct           1
LanczosEps     14
LanczosTarget  2
LargeValue     8.0
NumAve         5
ExpecInterval  20
NOmega         5
OmegaOrg       -2.1027484834620633
OmegaMin       0.0
OmegaMax       0.5
OmegaIm        0.1
EOF

cat > singleexcitation.def <<EOF
=============================================
NSingleExcitation 1
=============================================
======== Single Excitation ==================
=============================================
1 0 0         1.000000000000000         0.000000000000000
EOF

cat > singleexcitationbra.def <<EOF
=============================================
NSingleExcitationBra 1
=============================================
======== Single Excitation Bra ==============
=============================================
0 0 0         1.000000000000000         0.000000000000000
EOF

cat > namelist_cg.def <<EOF
         ModPara  SpectrumModpara
         LocSpin  locspn.def
           Trans  trans.def
    CoulombIntra  coulombintra.def
        OneBodyG  greenone.def
        TwoBodyG  greentwo.def
         CalcMod  calcmod_cg.def
  SingleExcitation  singleexcitation.def
  SingleExcitationBra  singleexcitationbra.def
     SpectrumVec  zvo_eigenvec_0
EOF

set +e
${HPHI} -e namelist_cg.def > run.log 2>&1
rc=$?
set -e

if [ "${rc}" = "0" ]; then
    echo "ERROR: BiCG spectrum run succeeded but was expected to fail" >&2
    tail -40 run.log >&2
    exit 1
fi

grep -q "BiCG spectrum did not finish successfully" run.log
grep -q "BiCG iteration-1 diagnostic:" run.log
grep -q "v2_weighted_sum=" run.log
grep -q "Hv2_weighted_sum=" run.log

if [ -e output/zvo_DynamicalGreen.dat ]; then
    echo "ERROR: non-converged BiCG spectrum wrote DynamicalGreen output" >&2
    cat output/zvo_DynamicalGreen.dat >&2
    exit 1
fi

rm -f output/zvo_TMComponents.dat output/zvo_DynamicalGreen.dat

sed -e 's/^CalcSpec.*/CalcSpec   3/' calcmod_cg.def > calcmod_restart_out.def
sed -e 's/CalcMod  calcmod_cg.def/CalcMod  calcmod_restart_out.def/' \
    -e '/SingleExcitationBra/d' \
    namelist_cg.def > namelist_restart_out.def

set +e
${HPHI} -e namelist_restart_out.def > restart_out.log 2>&1
rc=$?
set -e

if [ "${rc}" = "0" ]; then
    echo "ERROR: BiCG spectrum restart_out run succeeded but was expected to fail" >&2
    tail -40 restart_out.log >&2
    exit 1
fi

grep -q "BiCG spectrum did not finish successfully" restart_out.log

if [ ! -s output/zvo_TMComponents.dat ]; then
    echo "ERROR: non-converged CalcSpec=restart_out did not save TMComponents" >&2
    tail -40 restart_out.log >&2
    exit 1
fi

if [ -e output/zvo_DynamicalGreen.dat ]; then
    echo "ERROR: non-converged BiCG restart_out wrote DynamicalGreen output" >&2
    cat output/zvo_DynamicalGreen.dat >&2
    exit 1
fi

echo "BiCG spectrum non-convergence is rejected: OK"
