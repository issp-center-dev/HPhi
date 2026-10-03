#!/bin/sh
set -eu
mkdir -p tpq_zero_norm
cd tpq_zero_norm
cat > calc.def <<'EOF'
CalcType 1
CalcModel 1
OutputMode 0
EOF
cat > mod.def <<'EOF'
====
Model_Parameters 0
====
====
====
CDataFileHead zvo
CParaFileHead zqp
====
Nsite 4
2Sz 0
Lanczos_max 5
LanczosEps 12
exct 1
initial_iv 1
NumAve 1
LargeValue 1
EOF
cat > loc.def <<'EOF'
====
NlocalSpin 4
====
====
====
0 1
1 1
2 1
3 1
EOF
cat > coulomb.def <<'EOF'
====
NCoulombInter 1
====
====
====
0 1 4
EOF
cat > namelist.def <<'EOF'
CalcMod calc.def
ModPara mod.def
LocSpin loc.def
CoulombInter coulomb.def
EOF
# For Spin, density is the identity: H=4I, hence (l-H/Nsite)=0.
# This tests caller error propagation in serial; MPI reductions are covered
# separately by unittest_tpq_failure_mpi with zero-row ranks.
if ../../src/HPhi -e namelist.def > failure.log 2>&1; then
    cat failure.log
    echo 'mTPQ accepted an annihilated initial vector' >&2
    exit 1
fi
if ! grep -q 'mTPQ first step has zero or non-finite global norm' failure.log; then
    cat failure.log
    exit 1
fi
