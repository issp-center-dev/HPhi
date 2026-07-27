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
PreCG 0
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
grep -q "Symmetry allocation: raw_dim=6 global_dim=2 local_dim=2 raw_basis_list_elements=0 raw_diagonal_elements=0 initial_vector_elements=9" symmetry.log
grep -q "Symmetry LOBPCG allocation: local_dim=2 exct=1 workspace_vector_elements=18" symmetry.log
grep -q "Symmetry distributed matvec: global_rows=2" symmetry.log
grep -q "columns=local/ghost-slots" symmetry.log
grep -q \
    "Symmetry basis layout: distributed (default for TransSym CG)." \
    symmetry.log

rm -rf output
env HPHI_SYMMETRY_BASIS_LAYOUT=replicated \
    ../../src/HPhi -e namelist.def > symmetry_replicated.log 2>&1
replicated_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
test -n "${replicated_energy}"
replicated_diff=`awk -v a="${replicated_energy}" -v b="${sym_energy}" \
    'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
test "${replicated_diff}" = "0.000000"
grep -q "Symmetry matvec: mode=plan vector_exchange=halo" \
    symmetry_replicated.log
grep -q \
    "Symmetry basis layout: replicated (explicit rollback for TransSym CG)." \
    symmetry_replicated.log

rm -rf output
env HPHI_SYMMETRY_BASIS_LAYOUT=distributed \
    ../../src/HPhi -e namelist.def > symmetry_distributed.log 2>&1
distributed_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
test -n "${distributed_energy}"
distributed_diff=`awk -v a="${distributed_energy}" -v b="${sym_energy}" \
    'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
test "${distributed_diff}" = "0.000000"
grep -q "Symmetry distributed matvec:" symmetry_distributed.log
grep -q \
    "Symmetry basis layout: distributed (explicit environment)." \
    symmetry_distributed.log

if [ -n "${MPIRUN}" ]; then
    MPI_NP=`printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}'`
    if printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$" && [ "${MPI_NP}" -gt 1 ]; then
        rm -rf output
        ${MPIRUN} ../../src/HPhi -e namelist.def > symmetry_mpi.log 2>&1
        mpi_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
        test -n "${mpi_energy}"
        mpi_diff=`awk -v a="${mpi_energy}" -v b="${sym_energy}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
        test "${mpi_diff}" = "0.000000"
        grep -q "Symmetry basis: raw_dim=6 sector_dim=2 group_order=4" symmetry_mpi.log
        grep -q "raw_basis_list_elements=0 raw_diagonal_elements=0" symmetry_mpi.log
        grep -Eq "Symmetry LOBPCG allocation: local_dim=[01] exct=1 workspace_vector_elements=(6|12)" symmetry_mpi.log
        grep -q "Symmetry distributed matvec: global_rows=2" symmetry_mpi.log
        grep -q "columns=local/ghost-slots" symmetry_mpi.log
        grep -q \
            "Symmetry basis layout: distributed (default for TransSym CG)." \
            symmetry_mpi.log
        if grep -q "MPI site separation summary" symmetry_mpi.log; then
            echo "TransSym MPI path unexpectedly used site decomposition."
            exit 1
        fi
    fi
fi

exit $?
