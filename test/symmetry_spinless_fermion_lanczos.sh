#!/bin/sh -e

mkdir -p symmetry_spinless_fermion_lanczos
cd symmetry_spinless_fermion_lanczos

write_calcmod() {
cat > calcmod.def <<EOF
CalcType 0
CalcModel 7
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
}

write_modpara() {
    ncond="$1"
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
Ncond ${ncond}
Lanczos_max 20
initial_iv -1
exct 1
LanczosEps 12
LanczosTarget 1
LargeValue 50
EOF
}

write_locspn() {
cat > locspn.def <<EOF
================================
NlocalSpin     0
================================
========i_0LocSpn_Sr=Sr_i=======
================================
EOF
}

write_transfer_ring() {
cat > transfer.def <<EOF
================
NTransfer 8
================
========i s j t t_ij======
================
0 0 1 0 1.0 0.0
1 0 0 0 1.0 0.0
1 0 2 0 1.0 0.0
2 0 1 0 1.0 0.0
2 0 3 0 1.0 0.0
3 0 2 0 1.0 0.0
3 0 0 0 1.0 0.0
0 0 3 0 1.0 0.0
EOF
}

write_coulombinter() {
cat > coulombinter.def <<EOF
================
NCoulombInter 1
================
========i_j_V ======
================
0 1 0.25
EOF
}

write_k0_transsym() {
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

write_kpi2_transsym() {
cat > qptransidx.def <<EOF
=============================================
NQPTrans          4
=============================================
======== TrIdx_TrWeight_and_TrIdx_i_xi ======
=============================================
0 1.0 0.0
1 0.0 -1.0
2 -1.0 0.0
3 0.0 1.0
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

write_namelist() {
cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Trans transfer.def
TransSym qptransidx.def
EOF
}

assert_energy() {
    expected="$1"
    log="$2"
    energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
    test -n "${energy}"
    diff=`awk -v a="${energy}" -v b="${expected}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
    if [ "${diff}" != "0.000000" ]; then
        cat "${log}"
        echo "Energy mismatch: got ${energy}, expected ${expected}"
        exit 1
    fi
}

expect_failure() {
    pattern="$1"
    log="$2"
    shift 2
    if "$@" > "${log}" 2>&1; then
        cat "${log}"
        exit 1
    fi
    if ! grep -q "${pattern}" "${log}"; then
        cat "${log}"
        echo "Expected pattern not found: ${pattern}"
        exit 1
    fi
}

write_calcmod
write_locspn
write_modpara 1
write_transfer_ring
write_k0_transsym
write_namelist

../../src/HPhi -e namelist.def > spinless_k0.log 2>&1
assert_energy "-2.0" spinless_k0.log
grep -q "Symmetry basis: raw_dim=4 sector_dim=1 group_order=4" spinless_k0.log

rm -rf output
write_modpara 2
write_kpi2_transsym

../../src/HPhi -e namelist.def > spinless_kpi2.log 2>&1
assert_energy "-2.0" spinless_kpi2.log
grep -q "Symmetry basis: raw_dim=6 sector_dim=2 group_order=4" spinless_kpi2.log

rm -rf output
write_k0_transsym
perl -0pi -e 's/2 0 3 0 1.0 0.0/2 0 3 0 0.5 0.0/' transfer.def
perl -0pi -e 's/3 0 2 0 1.0 0.0/3 0 2 0 0.5 0.0/' transfer.def
expect_failure "SpinlessFermion Transfer invariance failed" \
    noninvariant_transfer.log ../../src/HPhi -e namelist.def

write_transfer_ring
perl -0pi -e 's/CalcType 0/CalcType 1/' calcmod.def
expect_failure "supports only Lanczos and CG" \
    unsupported_method.log ../../src/HPhi -e namelist.def

write_calcmod
write_coulombinter
cat >> namelist.def <<EOF
CoulombInter coulombinter.def
EOF
expect_failure "SpinlessFermion symmetry basis supports Transfer terms only" \
    unsupported_term.log ../../src/HPhi -e namelist.def
