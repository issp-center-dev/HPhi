#!/bin/sh -e

testname="lanczos_small_dimension_regression"
hphi="../../src/HPhi"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
    log="$1"
    shift
    "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

expect_fail() {
    log="$1"
    shift
    if "$@" > "${log}" 2>&1; then
        cat "${log}"
        exit 1
    fi
}

write_calcmod() {
    restart="$1"
    calc_eigenvec="${2:--1}"
    cat > calcmod.def <<EOF
CalcType 0
CalcModel 1
OutputMode 0
CalcEigenVec ${calc_eigenvec}
InitialVecType 0
OutputEigenVec 0
InputEigenVec 0
OutputHam 0
InputHam 0
ReStart ${restart}
CalcSpec 0
EOF
}

write_modpara() {
    nsite="$1"
    twosz="$2"
    initial_iv="$3"
    exct="$4"
    target="$5"
    nvec="${6:-}"
    cat > modpara.def <<EOF
--------------------
Model_Parameters 0
--------------------
--------------------
--------------------
CDataFileHead zvo
CParaFileHead zqp
--------------------
Nsite ${nsite}
2Sz ${twosz}
Lanczos_max 4
initial_iv ${initial_iv}
exct ${exct}
LanczosEps 12
LanczosTarget ${target}
LargeValue 50
EOF
    if [ -n "${nvec}" ]; then
        printf "nvec %s\n" "${nvec}" >> modpara.def
    fi
}

write_locspn() {
    nsite="$1"
    {
        cat <<EOF
================
NlocalSpin ${nsite}
================
========i_1LocSpn ======
================
EOF
        site=0
        while [ "${site}" -lt "${nsite}" ]; do
            printf "%d 1\n" "${site}"
            site=$((site + 1))
        done
    } > locspn.def
}

write_exchange_l2() {
    cat > exchange.def <<EOF
================
NExchange 1
================
========i_j_J ======
================
0 1 1.0
EOF
    cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
EOF
}

write_ising_l2() {
    cat > ising.def <<EOF
================
NIsing 1
================
========i_j_J ======
================
0 1 1.0
EOF
    cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Ising ising.def
EOF
}

write_calcmod 0 0
write_modpara 2 0 1 1 1
write_locspn 2
write_ising_l2
rm -rf output
run_hphi beta_zero.log "${hphi}" -e namelist.def
grep -q "Lanczos Krylov space exhausted" beta_zero.log
if grep -qi "nan" beta_zero.log; then
    cat beta_zero.log
    exit 1
fi
if grep -qi "nan" output/zvo_energy.dat; then
    cat output/zvo_energy.dat
    exit 1
fi

write_calcmod 0
write_modpara 2 2 1 1 2
write_locspn 2
write_exchange_l2
rm -rf output
run_hphi dim1.log "${hphi}" -e namelist.def
grep -q "Lanczos Hilbert space exhausted" dim1.log
grep -q "LanczosTarget=2 is outside Hilbert dimension 1; use 0" dim1.log

write_calcmod 0
write_modpara 2 0 -1 3 3
write_locspn 2
write_exchange_l2
rm -rf output
expect_fail exct_too_large.log "${hphi}" -e namelist.def
grep -q "Hilbert dimension 2 is smaller than exct=3" exct_too_large.log

write_calcmod 0
write_modpara 2 0 -1 1 -1
write_locspn 2
write_exchange_l2
rm -rf output
expect_fail negative_target.log "${hphi}" -e namelist.def
grep -q "LanczosTarget=-1 must be non-negative" negative_target.log

write_calcmod 1
write_modpara 2 0 -1 1 1
write_locspn 2
write_exchange_l2
rm -rf output
run_hphi restart_out_full_krylov.log "${hphi}" -e namelist.def
grep -q "Lanczos Krylov space exhausted" restart_out_full_krylov.log
test ! -f output/zvo_TMComponents.dat
test ! -f output/zvo_recalcvec_rank_0.dat

write_calcmod 0
write_modpara 2 0 -1 1 1 1
write_locspn 2
write_exchange_l2
rm -rf output
run_hphi nvec_target_energy.log "${hphi}" -e namelist.def
grep -q "Lanczos EigenValue = 1.0000000000" nvec_target_energy.log

write_calcmod 0
write_modpara 2 0 1 1 1
write_locspn 2
write_ising_l2
rm -rf output
run_hphi unreachable_target_warning.log "${hphi}" -e namelist.def
grep -q "LanczosTarget=1 is unreachable from this initial vector" unreachable_target_warning.log

exit 0
