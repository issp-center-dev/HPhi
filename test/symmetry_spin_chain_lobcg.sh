#!/bin/sh -e

mkdir -p symmetry_spin_chain_lobcg
cd symmetry_spin_chain_lobcg

cat > calcmod.def <<EOF
CalcType 3
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
Lanczos_max 20
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

cat > namelist_ref.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
EOF

../../src/HPhi -e namelist_ref.def > reference.log 2>&1
ref_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
rm -rf output

cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
TransSym qptransidx.def
EOF

../../src/HPhi -e namelist.def > symmetry.log 2>&1
sym_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
test -n "${ref_energy}"
test -n "${sym_energy}"
diff=`awk -v a="${sym_energy}" -v b="${ref_energy}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
test "${diff}" = "0.000000"

grep -q "Symmetry basis: raw_dim=6 sector_dim=2 group_order=4" symmetry.log

exit $?
