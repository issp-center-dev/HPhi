#!/bin/sh -e

mkdir -p symmetry_spin_chain_standard_generation
cd symmetry_spin_chain_standard_generation

run_hphi() {
    log="$1"
    shift
    "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

write_exchange_only_stan() {
    momentum_line="$1"
    cat > stan.in <<EOF
L = 6
model = Spin
method = Lanczos
lattice = chain
outputmode = none
Jx = 1.0
Jy = 1.0
Jz = 0.0
2Sz = 0
Lanczos_max = 20
initial_iv = -1
exct = 1
LanczosEps = 12
LanczosTarget = 1
${momentum_line}
EOF
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

write_exchange_only_stan "MomentumIndex = 1"
rm -rf output
run_hphi auto.log ../../src/HPhi -s stan.in
auto_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
test -n "${auto_energy}"
grep -q "qptransidx.def is written for MomentumIndex = 1" auto.log
grep -q "TransSym  qptransidx.def" namelist.def
grep -q "Symmetry basis: raw_dim=20 sector_dim=3 group_order=6" auto.log
awk 'NF == 3 && $1 == 1 {found=1; re=$2; im=$3; if(re<0) re=-re; if(im+0.8660254037844386<0) d=-(im+0.8660254037844386); else d=(im+0.8660254037844386); ok=(re-0.5 < 0.000001 && re-0.5 > -0.000001 && d < 0.000001)} END{exit found && ok ? 0 : 1}' qptransidx.def

rm -rf expert output
mkdir expert
cd expert
write_exchange_only_stan ""
run_hphi expert_sdry.log ../../../src/HPhi -sdry stan.in
write_kpi_over_3_transsym_l6
printf "        TransSym  qptransidx.def\n" >> namelist.def
run_hphi expert.log ../../../src/HPhi -e namelist.def
expert_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
test -n "${expert_energy}"
diff=`awk -v a="${auto_energy}" -v b="${expert_energy}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
test "${diff}" = "0.000000"
grep -q "Symmetry basis: raw_dim=20 sector_dim=3 group_order=6" expert.log
cd ..

write_exchange_only_stan "MomentumIndex = 6"
if ../../src/HPhi -sdry stan.in > invalid_index.log 2>&1; then
    cat invalid_index.log
    exit 1
fi
grep -q "MomentumIndex must satisfy" invalid_index.log

cat > stan.in <<EOF
L = 6
model = Spin
method = Lanczos
lattice = chain
outputmode = none
Jx = 1.0
Jy = 1.0
Jz = 0.0
phase0 = 180.0
2Sz = 0
MomentumIndex = 1
EOF
if ../../src/HPhi -sdry stan.in > boundary_phase.log 2>&1; then
    cat boundary_phase.log
    exit 1
fi
grep -q "MomentumIndex does not support boundary phase" boundary_phase.log

exit 0
