#!/bin/sh -e

# Regression coverage for SpectrumLoopExct + SpectrumNumOp additional
# single_ex_<op>.def readers. The additional files must accept the normal
# SingleExcitation header style and reject invalid operator rows before the
# spectrum calculation starts.

HPHI=../../src/HPhi

mkdir -p spectrum_loop_multiop_reader_validation/
cd spectrum_loop_multiop_reader_validation

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

${HPHI} -s stan_gs.in

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
Lanczos_max    100
initial_iv     -1
exct           1
LanczosEps     14
LanczosTarget  2
LargeValue     8.0
NumAve         5
ExpecInterval  20
NOmega         2
OmegaOrg       -2.1027484834620633
OmegaMin       0.0
OmegaMax       0.2
OmegaIm        0.1
SpectrumLoopExct 1
SpectrumNumOp 2
EOF

cat > singleexcitation.def <<EOF
=============================================
NSingleExcitation 1
=============================================
======== Single Excitation ==================
=============================================
1 0 0         1.000000000000000         0.000000000000000
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
     SpectrumVec  zvo_eigenvec
EOF

write_single_ex_1() {
  site="$1"
  spin="$2"
  type="$3"
  cat > single_ex_1.def <<EOF
=============================================
NSingleExcitation 1
=============================================
======== Single Excitation ==================
=============================================
${site} ${spin} ${type}         1.000000000000000         0.000000000000000
EOF
}

expect_reject() {
  name="$1"
  pat="$2"
  if ${HPHI} -e namelist_cg.def > "reject_${name}.log" 2>&1; then
    echo "FAIL: invalid case '${name}' was ACCEPTED"
    cat "reject_${name}.log"
    exit 1
  fi
  if ! grep -q "${pat}" "reject_${name}.log"; then
    echo "FAIL: '${name}' was not rejected by the expected guard '${pat}'"
    cat "reject_${name}.log"
    exit 1
  fi
  echo "OK: '${name}' rejected by the expected guard"
}

# The additional file intentionally uses the standard SingleExcitation keyword,
# not the DCore-specific "NSingle" spelling.
write_single_ex_1 1 0 0
${HPHI} -e namelist_cg.def > accept_nsingleexcitation.log 2>&1
test -s output/zvo_DynamicalGreen_0_0.dat
test -s output/zvo_DynamicalGreen_0_1.dat

write_single_ex_1 99 0 0
expect_reject "invalid_site" "invalid site index"

write_single_ex_1 1 0 2
expect_reject "invalid_type" "invalid single-excitation type"

exit 0
