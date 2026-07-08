#!/bin/sh -e

mkdir -p symmetry_spinless_fermion_standard_generation
cd symmetry_spinless_fermion_standard_generation

run_hphi() {
    log="$1"
    shift
    "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

assert_close() {
    expected="$1"
    actual="$2"
    awk -v a="${actual}" -v b="${expected}" 'BEGIN{d=a-b; if(d<0)d=-d; exit d < 0.000001 ? 0 : 1}'
}

write_spinless_stan() {
    ncond="$1"
    momentum_line="$2"
    cat > stan.in <<EOF
L = 4
model = SpinlessFermion
method = Lanczos
lattice = chain
outputmode = none
t = 1.0
ncond = ${ncond}
Lanczos_max = 20
initial_iv = -1
exct = 1
LanczosEps = 12
LanczosTarget = 1
${momentum_line}
EOF
}

write_spinless_transsym_l4() {
    momentum_index="$1"
    {
    cat <<EOF
=============================================
NQPTrans          4
=============================================
======== TrIdx_TrWeight_and_TrIdx_i_xi ======
=============================================
EOF
    case "${momentum_index}" in
        0)
            printf "0 1.0 0.0\n"
            printf "1 1.0 0.0\n"
            printf "2 1.0 0.0\n"
            printf "3 1.0 0.0\n"
            ;;
        1)
            printf "0 1.0 0.0\n"
            printf "1 0.0 -1.0\n"
            printf "2 -1.0 0.0\n"
            printf "3 0.0 1.0\n"
            ;;
        *)
            echo "unsupported momentum index ${momentum_index}" >&2
            exit 1
            ;;
    esac
    op=0
    while [ "${op}" -lt 4 ]; do
        site=0
        while [ "${site}" -lt 4 ]; do
            printf "%d %d %d 1\n" "${op}" "${site}" "$(((site + op) % 4))"
            site=$((site + 1))
        done
        op=$((op + 1))
    done
    } > qptransidx.def
}

run_case() {
    label="$1"
    ncond="$2"
    momentum_index="$3"
    expected_raw_dim="$4"
    expected_sector_dim="$5"
    expected_energy="$6"
    char_op1_re="$7"
    char_op1_im="$8"

    rm -rf "${label}"
    mkdir "${label}"
    cd "${label}"

    write_spinless_stan "${ncond}" "MomentumIndex = ${momentum_index}"
    run_hphi auto.log ../../../src/HPhi -s stan.in
    auto_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
    test -n "${auto_energy}" || { cat auto.log; exit 1; }
    assert_close "${expected_energy}" "${auto_energy}"
    grep -q "qptransidx.def is written for MomentumIndex = ${momentum_index}" auto.log
    grep -q "TransSym  qptransidx.def" namelist.def
    grep -q "CalcModel   7" calcmod.def
    grep -q "Ncond          ${ncond}" modpara.def
    grep -q "Symmetry basis: raw_dim=${expected_raw_dim} sector_dim=${expected_sector_dim} group_order=4" auto.log
    awk -v re0="${char_op1_re}" -v im0="${char_op1_im}" '
      NF == 3 && $1 == 1 {
        found = 1;
        dre = $2 - re0; if (dre < 0) dre = -dre;
        dim = $3 - im0; if (dim < 0) dim = -dim;
        ok = (dre < 0.000001 && dim < 0.000001);
      }
      END { exit found && ok ? 0 : 1 }
    ' qptransidx.def

    mkdir expert
    cd expert
    write_spinless_stan "${ncond}" ""
    run_hphi expert_sdry.log ../../../../src/HPhi -sdry stan.in
    write_spinless_transsym_l4 "${momentum_index}"
    printf "        TransSym  qptransidx.def\n" >> namelist.def
    run_hphi expert.log ../../../../src/HPhi -e namelist.def
    expert_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
    test -n "${expert_energy}" || { cat expert.log; exit 1; }
    assert_close "${auto_energy}" "${expert_energy}"
    grep -q "Symmetry basis: raw_dim=${expected_raw_dim} sector_dim=${expected_sector_dim} group_order=4" expert.log
    cd ../..
}

rm -rf no_momentum
mkdir no_momentum
cd no_momentum
write_spinless_stan 1 ""
run_hphi no_momentum.log ../../../src/HPhi -s stan.in
no_momentum_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
test -n "${no_momentum_energy}" || { cat no_momentum.log; exit 1; }
assert_close "-2.0" "${no_momentum_energy}"
grep -q "CalcModel   7" calcmod.def
grep -q "Ncond          1" modpara.def
if grep -q "TransSym" namelist.def; then
    cat namelist.def
    exit 1
fi
cd ..

run_case k0_n1 1 0 4 1 -2.0 1.0 0.0
run_case kpi2_n2 2 1 6 2 -2.0 0.0 -1.0

rm -rf invalid_index
mkdir invalid_index
cd invalid_index
write_spinless_stan 1 "MomentumIndex = 4"
if ../../../src/HPhi -sdry stan.in > invalid_index.log 2>&1; then
    cat invalid_index.log
    exit 1
fi
grep -q "MomentumIndex must satisfy" invalid_index.log
cd ..

rm -rf boundary_phase
mkdir boundary_phase
cd boundary_phase
cat > stan.in <<EOF
L = 4
model = SpinlessFermion
method = Lanczos
lattice = chain
outputmode = none
t = 1.0
ncond = 1
phase0 = 180.0
MomentumIndex = 1
EOF
if ../../../src/HPhi -sdry stan.in > boundary_phase.log 2>&1; then
    cat boundary_phase.log
    exit 1
fi
grep -q "MomentumIndex does not support boundary phase" boundary_phase.log
cd ..

rm -rf density_reject
mkdir density_reject
cd density_reject
cat > stan.in <<EOF
L = 4
model = SpinlessFermion
method = Lanczos
lattice = chain
outputmode = none
t = 1.0
V = 0.5
ncond = 1
MomentumIndex = 1
EOF
if ../../../src/HPhi -sdry stan.in > density_reject.log 2>&1; then
    cat density_reject.log
    exit 1
fi
grep -q "V is SPECIFIED but will NOT be USED" density_reject.log
cd ..

rm -rf non_chain_reject
mkdir non_chain_reject
cd non_chain_reject
cat > stan.in <<EOF
L = 2
W = 2
model = SpinlessFermion
method = Lanczos
lattice = square
outputmode = none
t = 1.0
ncond = 1
EOF
if ../../../src/HPhi -sdry stan.in > non_chain_reject.log 2>&1; then
    cat non_chain_reject.log
    exit 1
fi
grep -q "SpinlessFermion Standard mode currently supports only chain lattice" non_chain_reject.log
cd ..

exit 0
