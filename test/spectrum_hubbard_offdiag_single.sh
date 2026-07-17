#!/bin/sh -e

# Off-diagonal dynamical Green function
#   G_BA(z) = <gs| c_{0,up}^dagger (1/(z-H)) c_{1,up} |gs>
# on a 1D Hubbard 4-site cluster (t=1, U=4, half filling, 2Sz=0), method=CG (BiCG).
# ket A = c_{1,up} (SingleExcitation), bra B = c_{0,up} (SingleExcitationBra).
# The reference values reproduce exact diagonalization (new_spectrum validation).
# Serial test (MPI consistency is covered separately).

HPHI=../../src/HPhi

mkdir -p spectrum_hubbard_offdiag_single/
cd spectrum_hubbard_offdiag_single

#
# Ground state (standard mode; StdFace writes the lattice .def files and the eigenvector)
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

#
# Off-diagonal spectrum (expert mode: reuses the StdFace lattice files,
# adds the ket/bra single-excitation operators)
#
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
Lanczos_max    2000
initial_iv     -1
exct           1
LanczosEps     14
LanczosTarget  2
# The ignored third field must not leak into later one-value Omega lines.
LargeValue     8.0  3.25
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

${HPHI} -e namelist_cg.def

cat > reference.dat <<EOF
0.0000000000 0.0000000000 -0.3354182045 -0.0537129306
0.1000000000 0.0000000000 -0.3952267119 -0.0749316242
0.2000000000 0.0000000000 -0.4790334010 -0.1113242043
0.3000000000 0.0000000000 -0.6027761001 -0.1810475321
0.4000000000 0.0000000000 -0.7942785902 -0.3376755873
EOF

paste output/zvo_DynamicalGreen.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste1.dat`

echo "spectrum_hubbard_offdiag_single: accumulated L1 diff vs reference = ${diff}"
# Reference values were generated with a different BLAS (Accelerate on macOS);
# allow a small cross-platform tolerance on the accumulated L1 difference.
# A genuine regression shifts the spectrum by O(0.1) or more, far above this.
awk -v d="${diff}" 'BEGIN { exit (d < 1.0e-3) ? 0 : 1 }'

exit $?
