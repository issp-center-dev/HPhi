#!/bin/sh -e

# Off-diagonal density-density dynamical Green function
#   G_BA(z) = <gs| n_{1,up} (1/(z-H)) n_{0,up} |gs>
# on a 1D Hubbard 4-site cluster (t=1, U=4, half filling, 2Sz=0), method=CG (BiCG).
# ket A = n_{0,up} (PairExcitation), bra B = n_{1,up} (PairExcitationBra).
# The reference values were cross-checked against exact diagonalization (~1e-6).
# Serial test (MPI consistency is covered separately).

HPHI=../../src/HPhi

mkdir -p spectrum_hubbard_offdiag_pair_density/
cd spectrum_hubbard_offdiag_pair_density

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
# Off-diagonal density-density spectrum (expert mode)
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
LargeValue     8.0
NumAve         5
ExpecInterval  20
NOmega         5
OmegaOrg       -2.1027484834620633
OmegaMin       0.0
OmegaMax       0.5
OmegaIm        0.1
EOF

cat > pairexcitation.def <<EOF
=============================================
NPairExcitation 1
=============================================
=============== Pair Excitation =============
=============================================
0 0 0 0 1         1.000000000000000         0.000000000000000
EOF

cat > pairexcitationbra.def <<EOF
=============================================
NPairExcitationBra 1
=============================================
=============== Pair Excitation Bra =========
=============================================
1 0 1 0 1         1.000000000000000         0.000000000000000
EOF

cat > namelist_cg.def <<EOF
         ModPara  SpectrumModpara
         LocSpin  locspn.def
           Trans  trans.def
    CoulombIntra  coulombintra.def
        OneBodyG  greenone.def
        TwoBodyG  greentwo.def
         CalcMod  calcmod_cg.def
  PairExcitation  pairexcitation.def
  PairExcitationBra  pairexcitationbra.def
     SpectrumVec  zvo_eigenvec_0
EOF

${HPHI} -e namelist_cg.def

cat > reference.dat <<EOF
0.0000000000 0.0000000000 0.3775198728 -2.3744943609
0.1000000000 0.0000000000 1.7520780671 -0.9972366790
0.2000000000 0.0000000000 1.6189779535 0.1362815906
0.3000000000 0.0000000000 0.7113192440 0.9748842355
0.4000000000 0.0000000000 -0.0181002508 0.4441586632
EOF

paste output/zvo_DynamicalGreen.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste1.dat`

echo "spectrum_hubbard_offdiag_pair_density: accumulated L1 diff vs reference = ${diff}"
# Reference values were generated with a different BLAS (Accelerate on macOS);
# allow a small cross-platform tolerance on the accumulated L1 difference.
# A genuine regression shifts the spectrum by O(0.1) or more, far above this.
awk -v d="${diff}" 'BEGIN { exit (d < 1.0e-3) ? 0 : 1 }'

exit $?
