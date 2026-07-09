#!/bin/sh -e

mkdir -p symmetry_hubbard_lanczos
cd symmetry_hubbard_lanczos

write_calcmod() {
cat > calcmod.def <<EOF
CalcType 0
CalcModel 0
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
Nup 1
Ndown 1
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
NTransfer 16
================
========i s j t t_ij======
================
0 0 1 0 1.0 0.0
1 0 0 0 1.0 0.0
0 1 1 1 1.0 0.0
1 1 0 1 1.0 0.0
1 0 2 0 1.0 0.0
2 0 1 0 1.0 0.0
1 1 2 1 1.0 0.0
2 1 1 1 1.0 0.0
2 0 3 0 1.0 0.0
3 0 2 0 1.0 0.0
2 1 3 1 1.0 0.0
3 1 2 1 1.0 0.0
3 0 0 0 1.0 0.0
0 0 3 0 1.0 0.0
3 1 0 1 1.0 0.0
0 1 3 1 1.0 0.0
EOF
}

write_coulombintra() {
cat > coulombintra.def <<EOF
================
NCoulombIntra 4
================
========i U_i ======
================
0 0.5
1 0.5
2 0.5
3 0.5
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

write_ref_namelist() {
cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Trans transfer.def
CoulombIntra coulombintra.def
EOF
}

write_sym_namelist() {
    with_coulomb="$1"
cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Trans transfer.def
TransSym qptransidx.def
EOF
    if [ "${with_coulomb}" = "yes" ]; then
cat >> namelist.def <<EOF
CoulombIntra coulombintra.def
EOF
    fi
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

assert_energy_matches_reference() {
    expected="$1"
    log="$2"
    energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
    test -n "${energy}"
    diff=`awk -v a="${energy}" -v b="${expected}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
    if [ "${diff}" != "0.000000" ]; then
        cat "${log}"
        echo "Energy mismatch: got ${energy}, reference ${expected}"
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
write_modpara
write_locspn
write_transfer_ring
write_coulombintra
write_ref_namelist

../../src/HPhi -e namelist.def > hubbard_ref.log 2>&1
ref_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
test -n "${ref_energy}"
rm -rf output

write_k0_transsym
write_sym_namelist yes
../../src/HPhi -e namelist.def > hubbard_k0.log 2>&1
assert_energy_matches_reference "${ref_energy}" hubbard_k0.log
grep -q "Symmetry basis: raw_dim=16 sector_dim=4 group_order=4" hubbard_k0.log

rm -rf output
write_kpi2_transsym
write_sym_namelist no
../../src/HPhi -e namelist.def > hubbard_kpi2.log 2>&1
assert_energy "-2.0" hubbard_kpi2.log
grep -q "Symmetry basis: raw_dim=16 sector_dim=4 group_order=4" hubbard_kpi2.log

rm -rf output
write_k0_transsym
write_sym_namelist yes
perl -0pi -e 's/2 0 3 0 1.0 0.0/2 0 3 0 0.5 0.0/' transfer.def
perl -0pi -e 's/3 0 2 0 1.0 0.0/3 0 2 0 0.5 0.0/' transfer.def
expect_failure "Hubbard Transfer invariance failed" \
    noninvariant_transfer.log ../../src/HPhi -e namelist.def

write_transfer_ring
perl -0pi -e 's/2 0.5/2 0.25/' coulombintra.def
expect_failure "Hubbard CoulombIntra invariance failed" \
    noninvariant_coulombintra.log ../../src/HPhi -e namelist.def

write_coulombintra
write_coulombinter
cat >> namelist.def <<EOF
CoulombInter coulombinter.def
EOF
expect_failure "Hubbard symmetry basis supports Transfer and CoulombIntra terms only" \
    unsupported_term.log ../../src/HPhi -e namelist.def

write_sym_namelist yes
if [ -n "${MPIRUN}" ]; then
    MPI_NP=`printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}'`
    if printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$" && [ "${MPI_NP}" -gt 1 ]; then
        expect_failure "Hubbard symmetry basis is serial-only" \
            hubbard_mpi_reject.log ${MPIRUN} ../../src/HPhi -e namelist.def
    fi
fi
