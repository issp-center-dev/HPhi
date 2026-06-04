#!/bin/sh -e

# Off-diagonal dynamical spin structure factor (Sz-Sz cross term) in the
# grand-canonical spin ensemble:  G(z) = <gs| Sz_1 (1/(z-H)) Sz_0 |gs>
# Heisenberg ring (chain, PBC), H = J sum S_i.S_j, J=1, method=CG.
# Covers the SpinGC leaves:
#   - S=1/2 : GetPairExcitedStateHalfSpinGC
#   - S=1   : GetPairExcitedStateGeneralSpinGC
# References were cross-checked against exact diagonalization (~1e-8).

HPHI=../../src/HPhi

mkdir -p spectrum_spingc_offdiag_szsz/
cd spectrum_spingc_offdiag_szsz

#
# --- S = 1/2 (HalfSpinGC leaf) ---
#
cat > stan_gs.in <<EOF
model = "SpinGC"
method = "CG"
lattice = "chain"
L = 4
J = 1.0
EigenVecIO = "out"
EOF
${HPHI} -s stan_gs.in

cat > calcmod_cg.def <<EOF
CalcType   3
CalcModel   4
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
Lanczos_max    2000
initial_iv     -1
exct           1
LanczosEps     14
LanczosTarget  2
LargeValue     8.0
NumAve         5
ExpecInterval  20
NOmega         5
OmegaOrg       -2.0
OmegaMin       0.0
OmegaMax       0.5
OmegaIm        0.1
EOF
cat > pairexcitation.def <<EOF
=============================================
NPairExcitation 2
=============================================
=============== Pair Excitation =============
=============================================
0 0 0 0 1        -0.500000000000000         0.000000000000000
0 1 0 1 1         0.500000000000000         0.000000000000000
EOF
cat > pairexcitationbra.def <<EOF
=============================================
NPairExcitationBra 2
=============================================
=============== Pair Excitation Bra =========
=============================================
1 0 1 0 1        -0.500000000000000         0.000000000000000
1 1 1 1 1         0.500000000000000         0.000000000000000
EOF
cat > namelist_cg.def <<EOF
         ModPara  SpectrumModpara
         LocSpin  locspn.def
           Trans  trans.def
    CoulombInter  coulombinter.def
            Hund  hund.def
        Exchange  exchange.def
        OneBodyG  greenone.def
        TwoBodyG  greentwo.def
         CalcMod  calcmod_cg.def
  PairExcitation  pairexcitation.def
  PairExcitationBra  pairexcitationbra.def
     SpectrumVec  zvo_eigenvec_0
EOF
${HPHI} -e namelist_cg.def
cat > reference.dat <<EOF
0.0000000000 0.0000000000 0.1650164998 0.0165016500
0.1000000000 0.0000000000 0.1829268272 0.0203252030
0.2000000000 0.0000000000 0.2051282028 0.0256410254
0.3000000000 0.0000000000 0.2333333308 0.0333333330
0.4000000000 0.0000000000 0.2702702673 0.0450450446
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste1.dat`

#
# --- S = 1 (GeneralSpinGC leaf) ---
#
cat > stan_gs.in <<EOF
model = "SpinGC"
method = "CG"
lattice = "chain"
L = 4
2S = 2
J = 1.0
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
Lanczos_max    2000
initial_iv     -1
exct           1
LanczosEps     14
LanczosTarget  2
LargeValue     8.0
NumAve         5
ExpecInterval  20
NOmega         5
OmegaOrg       -6.0
OmegaMin       0.0
OmegaMax       0.5
OmegaIm        0.1
EOF
cat > pairexcitation.def <<EOF
=============================================
NPairExcitation 3
=============================================
=============== Pair Excitation =============
=============================================
0 0 0 0 1        -1.000000000000000         0.000000000000000
0 1 0 1 1         0.000000000000000         0.000000000000000
0 2 0 2 1         1.000000000000000         0.000000000000000
EOF
cat > pairexcitationbra.def <<EOF
=============================================
NPairExcitationBra 3
=============================================
=============== Pair Excitation Bra =========
=============================================
1 0 1 0 1        -1.000000000000000         0.000000000000000
1 1 1 1 1         0.000000000000000         0.000000000000000
1 2 1 2 1         1.000000000000000         0.000000000000000
EOF
cat > namelist_cg.def <<EOF
         ModPara  SpectrumModpara
         LocSpin  locspn.def
           Trans  trans.def
        InterAll  interall.def
        OneBodyG  greenone.def
        TwoBodyG  greentwo.def
         CalcMod  calcmod_cg.def
  PairExcitation  pairexcitation.def
  PairExcitationBra  pairexcitationbra.def
     SpectrumVec  zvo_eigenvec_0
EOF
${HPHI} -e namelist_cg.def
cat > reference.dat <<EOF
0.0000000000 0.0000000000 0.4950495098 0.0495049516
0.1000000000 0.0000000000 0.5487804932 0.0609756111
0.2000000000 0.0000000000 0.6153846214 0.0769230786
0.3000000000 0.0000000000 0.7000000069 0.1000000020
0.4000000000 0.0000000000 0.8108108188 0.1351351377
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste2.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste2.dat`

test "${diff}" = "0.000000"

exit $?
