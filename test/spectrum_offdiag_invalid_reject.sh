#!/bin/sh -e

# Negative tests for the off-diagonal (bra) dynamical-Green guards.
# A 4-site Hubbard chain ground state is computed once; then several invalid
# bra-excitation inputs are fed to the spectrum run, each of which MUST be
# rejected (HPhi exits non-zero). The test fails if any invalid case is accepted.

HPHI=../../src/HPhi

mkdir -p spectrum_offdiag_invalid_reject/
cd spectrum_offdiag_invalid_reject

#
# Ground state + base lattice files
#
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
EOF

# Excitation operator building blocks
single_ket() {  # c_{1up} annihilation
  printf '%s\n' '=====' 'NSingleExcitation 1' '=====' '=====' '=====' \
    '1 0 0  1.0 0.0' > singleexcitation.def
}
pair_ket() {    # n_{0up}
  printf '%s\n' '=====' 'NPairExcitation 1' '=====' '=====' '=====' \
    '0 0 0 0 1  1.0 0.0' > pairexcitation.def
}
calcmod() {     # $1 = CalcType, $2 = CalcSpec
  cat > calcmod_cg.def <<EOF
CalcType   $1
CalcModel   0
ReStart   0
CalcSpec   $2
CalcEigenVec   0
InitialVecType   0
InputEigenVec   1
OutputEigenVec   0
InputHam   0
OutputHam   0
OutputExVec   0
EOF
}
nl() {  # build namelist with the given excitation keyword lines (passed as args)
  {
    echo '         ModPara  SpectrumModpara'
    echo '         LocSpin  locspn.def'
    echo '           Trans  trans.def'
    echo '    CoulombIntra  coulombintra.def'
    echo '        OneBodyG  greenone.def'
    echo '        TwoBodyG  greentwo.def'
    echo '         CalcMod  calcmod_cg.def'
    for line in "$@"; do echo "  $line"; done
    echo '     SpectrumVec  zvo_eigenvec_0'
  } > namelist_cg.def
}

fail=0
expect_reject() {  # $1 = case name, $2 = expected guard-message substring
  name="$1"
  pat="$2"
  if ${HPHI} -e namelist_cg.def > "reject_${name}.log" 2>&1; then
    echo "FAIL: invalid case '${name}' was ACCEPTED (should be rejected)"
    cat "reject_${name}.log"
    fail=1
  elif ! grep -q "${pat}" "reject_${name}.log"; then
    echo "FAIL: '${name}' rejected, but not by the expected guard ('${pat}')"
    cat "reject_${name}.log"
    fail=1
  else
    echo "OK: '${name}' rejected by the expected guard"
  fi
}

#
# 1. bra requires method=CG  (CalcType = Lanczos)
#
single_ket
printf '%s\n' '=====' 'NSingleExcitationBra 1' '=====' '=====' '=====' \
  '0 0 0  1.0 0.0' > singleexcitationbra.def
calcmod 0 1
nl 'SingleExcitation  singleexcitation.def' 'SingleExcitationBra  singleexcitationbra.def'
expect_reject "calctype_not_cg" "requires method"

#
# 2. bra requires CalcSpec = Normal  (CalcSpec = Restart)
#
calcmod 3 4
nl 'SingleExcitation  singleexcitation.def' 'SingleExcitationBra  singleexcitationbra.def'
expect_reject "calcspec_not_normal" "requires CalcSpec"

#
# 3. ket/bra excitation type mismatch  (ket single, bra pair)
#
calcmod 3 1
pair_ket  # not used by ket here, but present
printf '%s\n' '=====' 'NPairExcitationBra 1' '=====' '=====' '=====' \
  '1 0 1 0 1  1.0 0.0' > pairexcitationbra.def
nl 'SingleExcitation  singleexcitation.def' 'PairExcitationBra  pairexcitationbra.def'
expect_reject "type_mismatch" "must be the same type"

#
# 4. both SingleExcitationBra and PairExcitationBra specified
#
nl 'SingleExcitation  singleexcitation.def' \
   'SingleExcitationBra  singleexcitationbra.def' \
   'PairExcitationBra  pairexcitationbra.def'
expect_reject "both_bra_types" "cannot be used together"

#
# 5. ket/bra map to different Hilbert sectors
#    (ket c_{1up} annihilation, bra c_{0up} creation -> dNe = -1 vs +1)
#
printf '%s\n' '=====' 'NSingleExcitationBra 1' '=====' '=====' '=====' \
  '0 0 1  1.0 0.0' > singleexcitationbra.def
nl 'SingleExcitation  singleexcitation.def' 'SingleExcitationBra  singleexcitationbra.def'
expect_reject "sector_mismatch" "different Hilbert sectors"

#
# 6. bra operator set mixes inconsistent sector shifts
#    (annihilation + creation in one bra file)
#
printf '%s\n' '=====' 'NSingleExcitationBra 2' '=====' '=====' '=====' \
  '0 0 0  1.0 0.0' '0 0 1  1.0 0.0' > singleexcitationbra.def
nl 'SingleExcitation  singleexcitation.def' 'SingleExcitationBra  singleexcitationbra.def'
expect_reject "set_inconsistent" "mixes operators with different"

test "${fail}" = "0"
exit $?
