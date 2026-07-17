#!/bin/sh -e

mkdir -p symmetry_spin_chain_lanczos
cd symmetry_spin_chain_lanczos

run_hphi() {
    log="$1"
    shift
    "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

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
Nsite 6
2Sz 0
Lanczos_max 20
initial_iv -1
exct 1
LanczosEps 12
LanczosTarget 1
LargeValue 50
EOF

cat > locspn.def <<EOF
================
NlocalSpin 6
================
========i_1LocSpn ======
================
0 1
1 1
2 1
3 1
4 1
5 1
EOF

cat > exchange.def <<EOF
================
NExchange 6
================
========i_j_J ======
================
0 1 1.0
1 2 1.0
2 3 1.0
3 4 1.0
4 5 1.0
5 0 1.0
EOF

cat > ising.def <<EOF
================
NIsing 6
================
========i_j_J ======
================
0 1 1.0
1 2 1.0
2 3 1.0
3 4 1.0
4 5 1.0
5 0 1.0
EOF

write_kpi_transsym_l6() {
{
cat <<EOF
=============================================
NQPTrans          6
=============================================
======== TrIdx_TrWeight_and_TrIdx_i_xi ======
=============================================
EOF
op=0
while [ "${op}" -lt 6 ]; do
    if [ $((op % 2)) -eq 0 ]; then
        weight="1.0"
    else
        weight="-1.0"
    fi
    printf "%d %s\n" "${op}" "${weight}"
    op=$((op + 1))
done
op=0
while [ "${op}" -lt 6 ]; do
    site=0
    while [ "${site}" -lt 6 ]; do
        printf "%d %d %d 1\n" "${op}" "${site}" "$(((site + op) % 6))"
        site=$((site + 1))
    done
    op=$((op + 1))
done
} > qptransidx.def
}

write_kpi_over_3_transsym_l6() {
{
cat <<EOF
=============================================
NQPTrans          6
=============================================
======== TrIdx_TrWeight_and_TrIdx_i_xi ======
=============================================
0 1.0 0.0
1 0.5 -0.8660254037844386
2 -0.5 -0.8660254037844386
3 -1.0 0.0
4 -0.5 0.8660254037844386
5 0.5 0.8660254037844386
EOF
op=0
while [ "${op}" -lt 6 ]; do
    site=0
    while [ "${site}" -lt 6 ]; do
        printf "%d %d %d 1\n" "${op}" "${site}" "$(((site + op) % 6))"
        site=$((site + 1))
    done
    op=$((op + 1))
done
} > qptransidx.def
}

write_kpi_transsym_l6

cat > namelist_ref.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
EOF

run_hphi reference.log ../../src/HPhi -e namelist_ref.def
ref_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
rm -rf output

cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
TransSym qptransidx.def
EOF

run_hphi symmetry.log ../../src/HPhi -e namelist.def
sym_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
test -n "${ref_energy}"
test -n "${sym_energy}"
diff=`awk -v a="${sym_energy}" -v b="${ref_energy}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
test "${diff}" = "0.000000"

grep -q "Symmetry basis: raw_dim=20 sector_dim=4 group_order=6" symmetry.log

cat > namelist_heisenberg_ref.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
Ising ising.def
EOF

rm -rf output
run_hphi reference_heisenberg.log ../../src/HPhi -e namelist_heisenberg_ref.def
ref_heisenberg_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
rm -rf output

write_kpi_transsym_l6
cat > namelist_heisenberg.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
Ising ising.def
TransSym qptransidx.def
EOF

run_hphi symmetry_heisenberg.log ../../src/HPhi -e namelist_heisenberg.def
sym_heisenberg_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
test -n "${ref_heisenberg_energy}"
test -n "${sym_heisenberg_energy}"
heisenberg_diff=`awk -v a="${sym_heisenberg_energy}" -v b="${ref_heisenberg_energy}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
test "${heisenberg_diff}" = "0.000000"
grep -q "Symmetry basis: raw_dim=20 sector_dim=4 group_order=6" symmetry_heisenberg.log

rm -rf output
write_kpi_over_3_transsym_l6
run_hphi symmetry_complex.log ../../src/HPhi -e namelist.def
complex_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
test -n "${complex_energy}"
complex_diff=`awk -v a="${complex_energy}" 'BEGIN{d=a+1.0; if(d<0)d=-d; printf "%8.6f", d}'`
test "${complex_diff}" = "0.000000"
grep -q "Symmetry basis: raw_dim=20 sector_dim=3 group_order=6" symmetry_complex.log

run_mpi_symmetry_case() {
    label="$1"
    namelist="$2"
    expected_energy="$3"
    expected_dim="$4"
    log_file="symmetry_${label}_mpi.log"
    rm -rf output
    ${MPIRUN} ../../src/HPhi -e "${namelist}" > "${log_file}" 2>&1
    mpi_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
    test -n "${mpi_energy}"
    mpi_diff=`awk -v a="${mpi_energy}" -v b="${expected_energy}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
    test "${mpi_diff}" = "0.000000"
    grep -q "Symmetry basis: raw_dim=20 sector_dim=${expected_dim} group_order=6" "${log_file}"
    if grep -q "MPI site separation summary" "${log_file}"; then
        echo "TransSym MPI path unexpectedly used site decomposition."
        exit 1
    fi
}

if [ -n "${MPIRUN}" ]; then
    MPI_NP=`printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}'`
    if printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$" && [ "${MPI_NP}" -gt 1 ]; then
        write_kpi_transsym_l6
        run_mpi_symmetry_case heisenberg namelist_heisenberg.def "${sym_heisenberg_energy}" 4
        write_kpi_over_3_transsym_l6
        run_mpi_symmetry_case complex namelist.def "${complex_energy}" 3
    fi
fi

exit $?
