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

write_coulombinter_onsite() {
cat > coulombinter.def <<EOF
================
NCoulombInter 1
================
========i_j_V ======
================
0 0 0.50
EOF
}

write_coulombinter_ring() {
cat > coulombinter.def <<EOF
================
NCoulombInter 4
================
========i_j_V ======
================
0 1 0.25
1 2 0.25
2 3 0.25
3 0 0.25
EOF
}

write_coulombinter_nonuniform_ring() {
cat > coulombinter.def <<EOF
================
NCoulombInter 4
================
========i_j_V ======
================
0 1 0.25
1 2 0.50
2 3 0.25
3 0 0.25
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

run_mpi_symmetry_case() {
    label="$1"
    expected_energy="$2"
    expected_dim="$3"
    log_file="spinless_${label}_mpi.log"
    rm -rf output
    if ! ${MPIRUN} ../../src/HPhi -e namelist.def > "${log_file}" 2>&1; then
        cat "${log_file}"
        exit 1
    fi
    mpi_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
    if [ -z "${mpi_energy}" ]; then
        cat "${log_file}"
        echo "MPI energy was not written to output/zvo_energy.dat"
        exit 1
    fi
    mpi_diff=`awk -v a="${mpi_energy}" -v b="${expected_energy}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
    if [ "${mpi_diff}" != "0.000000" ]; then
        cat "${log_file}"
        echo "MPI energy mismatch: got ${mpi_energy}, expected ${expected_energy}"
        exit 1
    fi
    if ! grep -q "Symmetry basis: raw_dim=.* sector_dim=${expected_dim} group_order=4" "${log_file}"; then
        cat "${log_file}"
        echo "Expected SpinlessFermion symmetry sector_dim=${expected_dim} was not found"
        exit 1
    fi
    grep -q "vector_exchange=halo" "${log_file}"
    grep -q "columns=local/ghost-slots" "${log_file}"
    if grep -q "MPI site separation summary" "${log_file}"; then
        echo "TransSym SpinlessFermion MPI path unexpectedly used site decomposition."
        exit 1
    fi

    log_file="spinless_${label}_allgather_mpi.log"
    rm -rf output
    if ! env HPHI_SYMMETRY_VECTOR_EXCHANGE=allgather \
        ${MPIRUN} ../../src/HPhi -e namelist.def > "${log_file}" 2>&1; then
        cat "${log_file}"
        exit 1
    fi
    mpi_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
    test -n "${mpi_energy}"
    mpi_diff=`awk -v a="${mpi_energy}" -v b="${expected_energy}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
    if [ "${mpi_diff}" != "0.000000" ]; then
        cat "${log_file}"
        echo "MPI allgather energy mismatch: got ${mpi_energy}, expected ${expected_energy}"
        exit 1
    fi
    grep -q "Symmetry basis: raw_dim=.* sector_dim=${expected_dim} group_order=4" "${log_file}"
    grep -q "vector_exchange=allgather" "${log_file}"
    grep -q "columns=global" "${log_file}"
}

run_mpi_if_available() {
    label="$1"
    expected_energy="$2"
    expected_dim="$3"
    if [ -n "${MPIRUN}" ]; then
        MPI_NP=`printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}'`
        if printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$" && [ "${MPI_NP}" -gt 1 ]; then
            run_mpi_symmetry_case "$label" "$expected_energy" "$expected_dim"
        fi
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
grep -q "vector_exchange=halo" spinless_k0.log
grep -q "columns=local/ghost-slots" spinless_k0.log
run_mpi_if_available k0 "-2.0" 1

rm -rf output
write_modpara 2
write_kpi2_transsym

../../src/HPhi -e namelist.def > spinless_kpi2.log 2>&1
assert_energy "-2.0" spinless_kpi2.log
grep -q "Symmetry basis: raw_dim=6 sector_dim=2 group_order=4" spinless_kpi2.log
grep -q "vector_exchange=halo" spinless_kpi2.log
grep -q "columns=local/ghost-slots" spinless_kpi2.log
rm -rf output
env HPHI_SYMMETRY_VECTOR_EXCHANGE=allgather \
    ../../src/HPhi -e namelist.def > spinless_kpi2_allgather.log 2>&1
assert_energy "-2.0" spinless_kpi2_allgather.log
grep -q "vector_exchange=allgather" spinless_kpi2_allgather.log
grep -q "columns=global" spinless_kpi2_allgather.log
run_mpi_if_available kpi2 "-2.0" 2

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
write_coulombinter_ring
cat >> namelist.def <<EOF
CoulombInter coulombinter.def
EOF
rm -rf output
../../src/HPhi -e namelist.def > spinless_coulombinter.log 2>&1
assert_energy "0.25" spinless_coulombinter.log
grep -q "Symmetry basis: raw_dim=6 sector_dim=1 group_order=4" spinless_coulombinter.log
grep -q "vector_exchange=halo" spinless_coulombinter.log
test -s output/zvo_energy.dat
run_mpi_if_available coulombinter "0.25" 1
write_coulombinter
expect_failure "SpinlessFermion CoulombInter term 0 maps to a missing pair" \
    missing_coulombinter.log ../../src/HPhi -e namelist.def
write_coulombinter_onsite
expect_failure "SpinlessFermion CoulombInter term 0 is on-site" \
    onsite_coulombinter.log ../../src/HPhi -e namelist.def
write_coulombinter_nonuniform_ring
expect_failure "SpinlessFermion CoulombInter invariance failed" \
    noninvariant_coulombinter.log ../../src/HPhi -e namelist.def
if [ -n "${MPIRUN}" ]; then
    MPI_NP=`printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}'`
    if printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$" && [ "${MPI_NP}" -gt 1 ]; then
        write_coulombinter
        expect_failure "SpinlessFermion CoulombInter term 0 maps to a missing pair" \
            missing_coulombinter_mpi.log ${MPIRUN} ../../src/HPhi -e namelist.def
        write_coulombinter_onsite
        expect_failure "SpinlessFermion CoulombInter term 0 is on-site" \
            onsite_coulombinter_mpi.log ${MPIRUN} ../../src/HPhi -e namelist.def
        write_coulombinter_nonuniform_ring
        expect_failure "SpinlessFermion CoulombInter invariance failed" \
            noninvariant_coulombinter_mpi.log ${MPIRUN} ../../src/HPhi -e namelist.def
    fi
fi
