#!/bin/sh -e

# Regression for HubbardNConserved single-excitation spectrum with CalcHS=2.
# The vacuum excited sector (Ncond=1, annihilation -> Ne=0) used to call
# snoob(0) in sz_hacker_for_large_systems(). The full sector is checked too.

HPHI=../../../src/HPhi

mkdir -p spectrum_hubbard_nconserved_calchs2_boundary/
cd spectrum_hubbard_nconserved_calchs2_boundary

run_case() {
  name="$1"
  nelec="$2"
  op_type="$3"

  rm -rf "${name}"
  mkdir "${name}"
  cd "${name}"

  cat > stan_gs.in <<EOF
model = "Hubbard"
method = "CG"
lattice = "chain"
L = 4
t = 1.0
U = 4.0
nelec = ${nelec}
EigenVecIO = "out"
EOF

  ${HPHI} -s stan_gs.in > log_gs.txt 2>&1

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
Ncond          ${nelec}
CalcHS         2
Lanczos_max    100
initial_iv     -1
exct           1
LanczosEps     14
LanczosTarget  2
LargeValue     8.0
NumAve         5
ExpecInterval  20
NOmega         2
OmegaOrg       -2.0
OmegaMin       0.0
OmegaMax       0.2
OmegaIm        0.1
EOF

  cat > singleexcitation.def <<EOF
=============================================
NSingleExcitation 2
=============================================
======== Single Excitation ==================
=============================================
0 0 ${op_type}         1.000000000000000         0.000000000000000
0 1 ${op_type}         1.000000000000000         0.000000000000000
EOF

  cat > singleexcitationbra.def <<EOF
=============================================
NSingleExcitationBra 2
=============================================
======== Single Excitation Bra ==============
=============================================
0 0 ${op_type}         1.000000000000000         0.000000000000000
0 1 ${op_type}         1.000000000000000         0.000000000000000
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

  ${HPHI} -e namelist_cg.def > log_spectrum.txt 2>&1
  test -s output/zvo_DynamicalGreen.dat
  cd ..
}

run_case annihilate_to_vacuum 1 0
run_case create_to_full 7 1

exit 0
