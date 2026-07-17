#!/bin/sh -e

mkdir -p symmetry_hubbard_standard_generation
cd symmetry_hubbard_standard_generation

run_hphi() {
    log="$1"
    shift
    "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

write_hubbard_stan() {
    momentum_line="$1"
    u_value="$2"
    method_value="$3"
    if [ -z "${method_value}" ]; then
        method_value="Lanczos"
    fi
    cat > stan.in <<EOF
L = 4
model = FermionHubbard
method = ${method_value}
lattice = chain
outputmode = none
t = 1.0
U = ${u_value}
nelec = 2
2Sz = 0
Lanczos_max = 20
initial_iv = -1
exct = 1
LanczosEps = 12
LanczosTarget = 1
${momentum_line}
EOF
}

write_hubbard_stan_no_2sz() {
    cat > stan.in <<EOF
L = 4
model = FermionHubbard
method = Lanczos
lattice = chain
outputmode = none
t = 1.0
U = 0.5
nelec = 2
MomentumIndex = 0
EOF
}

write_square_hubbard_stan() {
    cat > stan.in <<EOF
L = 2
W = 2
model = FermionHubbard
method = Lanczos
lattice = square
outputmode = none
t = 1.0
U = 0.5
nelec = 2
2Sz = 0
MomentumIndex = 0
EOF
}

write_k0_transsym() {
cat > qptransidx.def <<EOF
=============================================
NQPTrans          4
=============================================
======== TrIdx_TrWeight_and_TrIdx_i_xi ======
=============================================
0 1.0 0.0
1 1.0 0.0
2 1.0 0.0
3 1.0 0.0
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

read_energy() {
    awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat
}

assert_energy_matches_reference() {
    expected="$1"
    log="$2"
    energy=`read_energy`
    test -n "${energy}"
    diff=`awk -v a="${energy}" -v b="${expected}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
    if [ "${diff}" != "0.000000" ]; then
        cat "${log}"
        echo "Energy mismatch: got ${energy}, reference ${expected}"
        exit 1
    fi
}

assert_symmetry_log() {
    log="$1"
    if ! grep -q "Symmetry basis: raw_dim=16 sector_dim=4 group_order=4" "${log}"; then
        cat "${log}"
        echo "Expected Hubbard symmetry sector_dim=4 was not found"
        exit 1
    fi
    if grep -q "MPI site separation summary" "${log}"; then
        cat "${log}"
        echo "TransSym Hubbard path unexpectedly used site decomposition."
        exit 1
    fi
}

assert_generated_files() {
    grep -q "qptransidx.def is written for MomentumIndex" "$1"
    grep -Eq "NQPTrans[[:space:]]+4" qptransidx.def
    grep -q "TransSym  qptransidx.def" namelist.def
    grep -q "Ncond" modpara.def
    grep -q "2Sz" modpara.def
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

run_mpi_generated_if_available() {
    label="$1"
    expected_energy="$2"
    if [ -n "${MPIRUN}" ]; then
        MPI_NP=`printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}'`
        if printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$" && [ "${MPI_NP}" -gt 1 ]; then
            log_file="auto_${label}_mpi.log"
            rm -rf output
            if ! ${MPIRUN} ../../src/HPhi -e namelist.def > "${log_file}" 2>&1; then
                cat "${log_file}"
                exit 1
            fi
            assert_energy_matches_reference "${expected_energy}" "${log_file}"
            assert_symmetry_log "${log_file}"
        fi
    fi
}

write_hubbard_stan "" "0.5" "Lanczos"
rm -rf output
run_hphi ref.log ../../src/HPhi -s stan.in
ref_energy=`read_energy`
test -n "${ref_energy}"

write_hubbard_stan "MomentumIndex = 0" "0.5" "Lanczos"
rm -rf output
run_hphi auto_k0.log ../../src/HPhi -s stan.in
auto_k0_energy=`read_energy`
test -n "${auto_k0_energy}"
assert_energy_matches_reference "${ref_energy}" auto_k0.log
assert_symmetry_log auto_k0.log
assert_generated_files auto_k0.log
run_mpi_generated_if_available k0 "${auto_k0_energy}"

rm -rf expert_k0
mkdir expert_k0
cd expert_k0
write_hubbard_stan "" "0.5" "Lanczos"
run_hphi expert_sdry.log ../../../src/HPhi -sdry stan.in
write_k0_transsym
printf "        TransSym  qptransidx.def\n" >> namelist.def
rm -rf output
run_hphi expert_k0.log ../../../src/HPhi -e namelist.def
assert_energy_matches_reference "${auto_k0_energy}" expert_k0.log
assert_symmetry_log expert_k0.log
cd ..

write_hubbard_stan "MomentumIndex = 1" "0.0" "Lanczos"
rm -rf output
run_hphi auto_k1.log ../../src/HPhi -s stan.in
auto_k1_energy=`read_energy`
test -n "${auto_k1_energy}"
assert_energy_matches_reference "-2.0" auto_k1.log
assert_symmetry_log auto_k1.log
assert_generated_files auto_k1.log
awk 'NF == 3 && $1 == 1 {found=1; ok=($2 < 0.000001 && $2 > -0.000001 && $3 + 1.0 < 0.000001 && $3 + 1.0 > -0.000001)} END{exit found && ok ? 0 : 1}' qptransidx.def
run_mpi_generated_if_available k1 "${auto_k1_energy}"

rm -rf expert_k1
mkdir expert_k1
cd expert_k1
write_hubbard_stan "" "0.0" "Lanczos"
run_hphi expert_sdry.log ../../../src/HPhi -sdry stan.in
write_kpi2_transsym
printf "        TransSym  qptransidx.def\n" >> namelist.def
rm -rf output
run_hphi expert_k1.log ../../../src/HPhi -e namelist.def
assert_energy_matches_reference "${auto_k1_energy}" expert_k1.log
assert_symmetry_log expert_k1.log
cd ..

write_hubbard_stan "MomentumIndex = 4" "0.5" "Lanczos"
expect_failure "MomentumIndex must satisfy" invalid_index.log ../../src/HPhi -sdry stan.in

write_hubbard_stan "MomentumIndex = 0" "0.5" "Lanczos"
cat >> stan.in <<EOF
phase0 = 180.0
EOF
expect_failure "MomentumIndex does not support boundary phase" boundary_phase.log ../../src/HPhi -sdry stan.in

write_hubbard_stan_no_2sz
expect_failure "MomentumIndex for Hubbard requires nelec and 2Sz" missing_2sz.log ../../src/HPhi -sdry stan.in

write_hubbard_stan "MomentumIndex = 0" "0.5" "Lanczos"
cat >> stan.in <<EOF
V = 0.25
EOF
expect_failure "MomentumIndex for Hubbard does not support V/CoulombInter" coulombinter_reject.log ../../src/HPhi -sdry stan.in

write_square_hubbard_stan
expect_failure "MomentumIndex currently supports only chain lattice" non_chain.log ../../src/HPhi -sdry stan.in

write_hubbard_stan "MomentumIndex = 0" "0.5" "fulldiag"
expect_failure "MomentumIndex currently supports only Lanczos and CG" fulldiag_reject.log ../../src/HPhi -sdry stan.in

exit 0
