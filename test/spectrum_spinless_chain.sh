#!/bin/sh -e
# Spinless-fermion pair-excitation dynamical Green's function (spectrum).
#
# Regression test for the spinless pair-excitation spectrum, which was both
# unreachable (no dispatch in GetPairExcitedState) and, once reached, produced
# an all-zero DynamicalGreen because:
#   (A) CisAjt_Hermite (the local spinless hopping kernel) wrote the output
#       vector tmp_v0 only in M_MLTPLY mode, never in M_CALCSPEC, so the excited
#       state c^+_i c_j |GS> (i != j) was never built during a spectrum run.
#   (B) GetPairExcitedState_SpinlessFermion had no isite1 == isite2 branch, so
#       the number operator n_i = c^+_i c_i was routed through the hopping kernel,
#       which returns 0 for creation and annihilation on the same site.
#
# Runs a small spinless chain in expert mode (single process; this exercises the
# local kernels), computes the ground state, then a density-density (diagonal
# n_i) dynamical Green's function. Before the fix every value is exactly 0;
# after the fix the spectrum is non-trivial.
#
# $1 = CMAKE_SOURCE_DIR (passed by ctest via add_hphi_test_with_srcdir)

HPHI="../../src/HPhi"

mkdir -p spectrum_spinless_chain/
cd spectrum_spinless_chain

python3 "$1/test/testSpinlessCalc.py" -p "${HPHI}" -mpi "" -m "SpinlessFermion" -s 6 -n 3 -V 1.0 > /dev/null 2>&1

cat > modpara.def <<'EOF'
--------------------
Model_Parameters   0
--------------------
HPhi_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
Nsite          6
Ncond          3
Lanczos_max    2000
initial_iv     1
exct           1
LanczosEps     14
LanczosTarget  2
LargeValue     12.0
NumAve         5
ExpecInterval  20
NOmega         5
OmegaMax       8.0    0.2
OmegaMin       -8.0   0.2
OmegaOrg       0.0    0.0
EOF

# density-density excitation operators: n_i = c^+_i c_i for each site
cat > pair.def <<'EOF'
=============================================
NPair 6
=============================================
=============== Pair Excitation =============
=============================================
0 0 0 0 1         1.000000000000000         0.000000000000000
1 0 1 0 1         1.000000000000000         0.000000000000000
2 0 2 0 1         1.000000000000000         0.000000000000000
3 0 3 0 1         1.000000000000000         0.000000000000000
4 0 4 0 1         1.000000000000000         0.000000000000000
5 0 5 0 1         1.000000000000000         0.000000000000000
EOF

# Step 1: ground state, output the eigenvector
cat > calcmod.def <<'EOF'
CalcType   0
CalcModel   7
ReStart   0
CalcSpec   0
CalcEigenVec   0
InitialVecType   0
InputEigenVec   0
OutputEigenVec   1
EOF
cat > namelist.def <<'EOF'
         ModPara  modpara.def
         CalcMod  calcmod.def
         LocSpin  locspn.def
           Trans  trans.def
    CoulombInter  coulombinter.def
        OneBodyG  greenone.def
        TwoBodyG  greentwo.def
EOF
rm -rf output
${HPHI} -e namelist.def

# Step 2: density-density dynamical Green's function (CG/BiCG spectrum)
cat > calcmod.def <<'EOF'
CalcType   3
CalcModel   7
ReStart   0
CalcSpec   1
CalcEigenVec   0
InitialVecType   0
InputEigenVec   0
OutputEigenVec   0
EOF
cat > namelist.def <<'EOF'
         ModPara  modpara.def
         CalcMod  calcmod.def
         LocSpin  locspn.def
           Trans  trans.def
    CoulombInter  coulombinter.def
        OneBodyG  greenone.def
        TwoBodyG  greentwo.def
  PairExcitation  pair.def
     SpectrumVec  zvo_eigenvec_0
EOF
${HPHI} -e namelist.def

# --- check: the dynamical Green's function must be non-trivial ---
# Column layout: Re(omega) Im(omega) Re(G) Im(G)
weight=$(awk 'BEGIN{s=0.0} {s += ($3<0?-$3:$3) + ($4<0?-$4:$4)} END{printf "%.10f", s}' output/zvo_DynamicalGreen.dat)
echo "Sum |Re(G)|+|Im(G)| over omega = ${weight}"

result=$(awk -v w="${weight}" 'BEGIN{ print (w > 0.000001) ? "PASS" : "FAIL" }')
if [ "${result}" = "PASS" ]; then
    echo "Spinless pair-excitation spectrum is non-trivial: PASS"
    exit 0
else
    echo "Spinless pair-excitation spectrum is all-zero (excited state not built): FAIL"
    cat output/zvo_DynamicalGreen.dat
    exit 1
fi
