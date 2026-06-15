#!/bin/sh -e

# MPI consistency for the off-diagonal dynamical Green function:
# the same G_BA = <gs| c_{0up}^dag (z-H)^-1 c_{1up} |gs> on a 4-site Hubbard
# chain (t=1, U=4) must reproduce the serial/exact result when run under MPI.
# A 4-site Hubbard requires the MPI process count to be a power of 4, so this
# test is registered to run only at np = 4 (skipped otherwise via the precheck).
# MPI introduces ~1e-7 BiCG-convergence differences, so a small tolerance is used.

if [ -z "${MPIRUN}" ]; then
    MPIRUN=""
fi

mkdir -p spectrum_hubbard_offdiag_single_mpi/
cd spectrum_hubbard_offdiag_single_mpi

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
${MPIRUN} ../../src/HPhi -s stan_gs.in

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
${MPIRUN} ../../src/HPhi -e namelist_cg.def

# Serial / exact reference (validated against exact diagonalization).
cat > reference.dat <<EOF
0.0000000000 0.0000000000 -0.3354182045 -0.0537129306
0.1000000000 0.0000000000 -0.3952267119 -0.0749316242
0.2000000000 0.0000000000 -0.4790334010 -0.1113242043
0.3000000000 0.0000000000 -0.6027761001 -0.1810475321
0.4000000000 0.0000000000 -0.7942785902 -0.3376755873
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%e", diff}' paste1.dat`
echo "MPI-vs-serial L2 diff = ${diff}"

ok=`awk "BEGIN{print (${diff} < 1.0e-5) ? 1 : 0}"`
test "${ok}" = "1"

exit $?
