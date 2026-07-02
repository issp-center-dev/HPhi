#!/bin/sh -e

mkdir -p symmetry_spin_chain_validation
cd symmetry_spin_chain_validation

write_common_defs() {
cat > calcmod.def <<EOF
CalcType 0
CalcModel 1
OutputMode 0
CalcEigenVec 0
InitialVecType 0
OutputEigenVec 0
InputEigenVec 0
OutputHam 0
InputHam 0
ReStart 0
CalcSpec 0
EOF
cat > modpara.def <<EOF
--------------------
Model_Parameters 0
--------------------
--------------------
--------------------
CDataFileHead zvo
CParaFileHead zqp
--------------------
Nsite 4
2Sz 0
Lanczos_max 8
initial_iv 1
exct 1
LanczosEps 12
LanczosTarget 2
LargeValue 50
EOF
cat > locspn.def <<EOF
================
NlocalSpin 4
================
========i_1LocSpn ======
================
0 1
1 1
2 1
3 1
EOF
cat > exchange.def <<EOF
================
NExchange 4
================
========i_j_J ======
================
0 1 1.0
1 2 1.0
2 3 1.0
3 0 1.0
EOF
}

write_valid_transsym() {
cat > qptransidx.def <<EOF
=============================================
NQPTrans          4
=============================================
======== TrIdx_TrWeight_and_TrIdx_i_xi ======
=============================================
0 1.0
1 1.0
2 1.0
3 1.0
0 0 0 1
0 1 1 1
0 2 2 1
0 3 3 1
1 0 1 1
1 1 2 1
1 2 3 1
1 3 0 1
2 0 2 1
2 1 3 1
2 2 0 1
2 3 1 1
3 0 3 1
3 1 0 1
3 2 1 1
3 3 2 1
EOF
}

write_common_defs
write_valid_transsym
cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
TransSym qptransidx.def
EOF

cat > qptransidx.def <<EOF
=============================================
NQPTrans          2
=============================================
======== TrIdx_TrWeight_and_TrIdx_i_xi ======
=============================================
0 1.0
1 1.0
0 0 0 1
0 1 1 1
0 2 2 1
0 3 3 1
1 0 1 1
1 1 2 1
1 2 3 1
1 3 0 1
EOF

if ../../src/HPhi -e namelist.def > nonclosed.log 2>&1; then
    cat nonclosed.log
    exit 1
fi
grep -q "not closed" nonclosed.log

write_valid_transsym
perl -0pi -e 's/1 0 1 1/1 0 1 -1/' qptransidx.def
if ../../src/HPhi -e namelist.def > anti.log 2>&1; then
    cat anti.log
    exit 1
fi
grep -q "Anti must be 1" anti.log

write_valid_transsym
perl -0pi -e 's/1 1 2 1/1 0 2 1/' qptransidx.def
if ../../src/HPhi -e namelist.def > duplicate.log 2>&1; then
    cat duplicate.log
    exit 1
fi
grep -q "duplicate TransSym permutation entry" duplicate.log

exit 0
