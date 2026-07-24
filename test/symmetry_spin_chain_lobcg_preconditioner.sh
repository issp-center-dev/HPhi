#!/bin/sh -eu

mkdir -p symmetry_spin_chain_lobcg_preconditioner
cd symmetry_spin_chain_lobcg_preconditioner

HPHI=../../src/HPhi
RUNNER=${MPIRUN:-}

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

cat > locspn.def <<EOF
================
NlocalSpin 8
================
========i_1LocSpn ======
================
0 1
1 1
2 1
3 1
4 1
5 1
6 1
7 1
EOF

cat > exchange.def <<EOF
================
NExchange 8
================
========i_j_J ======
================
0 1 1.0
1 2 1.0
2 3 1.0
3 4 1.0
4 5 1.0
5 6 1.0
6 7 1.0
7 0 1.0
EOF

{
    echo "============================================="
    echo "NQPTrans          8"
    echo "============================================="
    echo "======== TrIdx_TrWeight_and_TrIdx_i_xi ======"
    echo "============================================="
    translation=0
    while [ "${translation}" -lt 8 ]; do
        echo "${translation} 1.0"
        translation=$((translation + 1))
    done
    translation=0
    while [ "${translation}" -lt 8 ]; do
        site=0
        while [ "${site}" -lt 8 ]; do
            target=$(((site + translation) % 8))
            echo "${translation} ${site} ${target} 1"
            site=$((site + 1))
        done
        translation=$((translation + 1))
    done
} > qptransidx.def

cat > namelist_normal.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
EOF

cat > namelist_symmetry.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
TransSym qptransidx.def
EOF

run_case()
{
    label=$1
    precondition=$2
    namelist=$3

    cat > modpara.def <<EOF
--------------------
Model_Parameters 0
--------------------
--------------------
--------------------
CDataFileHead zvo
CParaFileHead zqp
--------------------
Nsite 8
2Sz 0
Lanczos_max 100
initial_iv 1
exct 1
LanczosEps 12
LanczosTarget 2
LargeValue 50
PreCG ${precondition}
EOF

    rm -rf output
    ${RUNNER} "${HPHI}" -e "${namelist}" > "${label}.log" 2>&1
    test -s output/zvo_energy.dat
    test -s output/zvo_Lanczos_Step.dat
    cp output/zvo_energy.dat "${label}_energy.dat"
    cp output/zvo_Lanczos_Step.dat "${label}_steps.dat"
}

check_convergence()
{
    step_file=$1
    awk '
        $1 ~ /^[0-9]+$/ {
            count++
            residual = $2
            threshold = $3
        }
        END {
            exit !(count >= 2 && residual < threshold)
        }
    ' "${step_file}"
}

extract_energy()
{
    awk '$1 == "Energy" {print $2; exit}' "$1"
}

check_energy_close()
{
    left=$1
    right=$2
    awk -v left="${left}" -v right="${right}" '
        BEGIN {
            difference = left - right
            if (difference < 0) difference = -difference
            exit !(difference <= 1.0e-10)
        }
    '
}

run_case normal_precg0 0 namelist_normal.def
run_case symmetry_precg0 0 namelist_symmetry.def
run_case symmetry_precg1 1 namelist_symmetry.def

check_convergence normal_precg0_steps.dat
check_convergence symmetry_precg0_steps.dat
check_convergence symmetry_precg1_steps.dat

normal_energy=$(extract_energy normal_precg0_energy.dat)
symmetry_precg0_energy=$(extract_energy symmetry_precg0_energy.dat)
symmetry_precg1_energy=$(extract_energy symmetry_precg1_energy.dat)
test -n "${normal_energy}"
test -n "${symmetry_precg0_energy}"
test -n "${symmetry_precg1_energy}"

check_energy_close "${symmetry_precg0_energy}" "${symmetry_precg1_energy}"
check_energy_close "${normal_energy}" "${symmetry_precg0_energy}"
check_energy_close "${normal_energy}" "${symmetry_precg1_energy}"

for log in symmetry_precg0.log symmetry_precg1.log; do
    grep -q \
        "Symmetry basis: raw_dim=70 sector_dim=10 group_order=8" "${log}"
    grep -q "Symmetry matvec: mode=plan vector_exchange=halo" "${log}"
    grep -q "columns=local/ghost-slots" "${log}"
    if grep -q "MPI site separation summary" "${log}"; then
        echo "TransSym MPI path unexpectedly used site decomposition."
        exit 1
    fi
done
