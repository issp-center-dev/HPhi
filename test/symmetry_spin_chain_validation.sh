#!/bin/sh -e

mkdir -p symmetry_spin_chain_validation
cd symmetry_spin_chain_validation

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

write_common_defs() {
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
Nsite 4
2Sz 0
Lanczos_max 8
initial_iv -1
exct 1
LanczosEps 12
LanczosTarget 1
LargeValue 50
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
}

write_valid_transsym() {
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

write_c4_kpi2_transsym() {
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

write_base_namelist() {
cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
TransSym qptransidx.def
EOF
}

write_common_defs
write_valid_transsym
write_base_namelist

cat > qptransidx.def <<EOF
=============================================
NQPTrans          2
=============================================
======== TrIdx_TrWeight_and_TrIdx_i_xi ======
=============================================
0 1.0
1 1.0
0 0 0 1
0 1 1 1
0 2 2 1
0 3 3 1
1 0 1 1
1 1 2 1
1 2 3 1
1 3 0 1
EOF

if ../../src/HPhi -e namelist.def > nonclosed.log 2>&1; then
    cat nonclosed.log
    exit 1
fi
grep -q "not closed" nonclosed.log

write_valid_transsym
perl -0pi -e 's/1 0 1 1/1 0 1 -1/' qptransidx.def
if ../../src/HPhi -e namelist.def > anti.log 2>&1; then
    cat anti.log
    exit 1
fi
grep -q "Anti must be 1" anti.log

write_valid_transsym
perl -0pi -e 's/1 1 2 1/1 0 2 1/' qptransidx.def
if ../../src/HPhi -e namelist.def > duplicate.log 2>&1; then
    cat duplicate.log
    exit 1
fi
grep -q "duplicate TransSym permutation entry" duplicate.log

write_common_defs
write_valid_transsym
perl -0pi -e 's/CalcType 0/CalcType 2/' calcmod.def
if ../../src/HPhi -e namelist.def > fulldiag.log 2>&1; then
    cat fulldiag.log
    exit 1
fi
grep -q "does not support FullDiag" fulldiag.log

for calc_type in 1 4 5; do
    write_common_defs
    write_valid_transsym
    perl -0pi -e "s/CalcType 0/CalcType ${calc_type}/" calcmod.def
    log="unsupported_calctype_${calc_type}.log"
    if ../../src/HPhi -e namelist.def > "${log}" 2>&1; then
        cat "${log}"
        exit 1
    fi
    grep -q "supports only Lanczos and CG" "${log}"
done

write_common_defs
write_valid_transsym
perl -0pi -e 's/2 3 1.0/2 3 0.75/' exchange.def
if ../../src/HPhi -e namelist.def > noninvariant.log 2>&1; then
    cat noninvariant.log
    exit 1
fi
grep -q "Hamiltonian invariance failed" noninvariant.log

write_common_defs
write_valid_transsym
perl -ni -e 'print unless /^2Sz[[:space:]]/' modpara.def
if ../../src/HPhi -e namelist.def > nosz.log 2>&1; then
    cat nosz.log
    exit 1
fi
grep -Eq "2Sz is not defined|requires fixed 2Sz" nosz.log

write_common_defs
write_valid_transsym
cat > interall.def <<EOF
================
NInterAll 1
================
========zInterAll ======
================
0 0 0 0 1 0 1 0 1.0 0.0
EOF
cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
InterAll interall.def
TransSym qptransidx.def
EOF
if ../../src/HPhi -e namelist.def > interall.log 2>&1; then
    cat interall.log
    exit 1
fi
grep -q "Exchange terms only" interall.log

write_common_defs
write_valid_transsym
cat > ising.def <<EOF
================
NIsing 4
================
========i_j_J ======
================
0 1 1.0
1 2 1.0
2 3 1.0
3 0 1.0
EOF
cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
Ising ising.def
TransSym qptransidx.def
EOF
if ../../src/HPhi -e namelist.def > ising.log 2>&1; then
    cat ising.log
    exit 1
fi
grep -q "Exchange terms only" ising.log

write_common_defs
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
3 0 1 1
3 1 2 1
3 2 3 1
3 3 0 1
EOF
write_base_namelist
if ../../src/HPhi -e namelist.def > duplicate_operation.log 2>&1; then
    cat duplicate_operation.log
    exit 1
fi
grep -q "duplicate TransSym operation" duplicate_operation.log

write_common_defs
write_valid_transsym
cat > exchange.def <<EOF
================
NExchange 5
================
========i_j_J ======
================
0 1 1.0
0 1 1.0
1 2 1.0
2 3 1.0
3 0 1.0
EOF
write_base_namelist
if ../../src/HPhi -e namelist.def > duplicate_bond.log 2>&1; then
    cat duplicate_bond.log
    exit 1
fi
grep -q "Hamiltonian invariance failed" duplicate_bond.log

for flag in OutputEigenVec InputEigenVec ReStart; do
    write_common_defs
    write_valid_transsym
    write_base_namelist
    perl -0pi -e "s/${flag} 0/${flag} 1/" calcmod.def
    log="unsupported_${flag}.log"
    if ../../src/HPhi -e namelist.def > "${log}" 2>&1; then
        cat "${log}"
        exit 1
    fi
    grep -q "EigenVec/ReStart" "${log}"
done

write_common_defs
write_valid_transsym
write_base_namelist
perl -0pi -e 's/^1 1\.0$/1 nan/m' qptransidx.def
if ../../src/HPhi -e namelist.def > nan_character.log 2>&1; then
    cat nan_character.log
    exit 1
fi
grep -q "character must be finite" nan_character.log

write_common_defs
write_valid_transsym
write_base_namelist
perl -0pi -e 's/^2Sz 0$/Nup 2\nNdown 2/m' modpara.def
rm -rf output
run_hphi nup_ndown.log ../../src/HPhi -e namelist.def
grep -q "Symmetry basis: raw_dim=6 sector_dim=2 group_order=4" nup_ndown.log
awk '$1 == "Sz" {found=1; d=$2; if(d<0)d=-d; ok=(d < 0.000001)} END{exit found && ok ? 0 : 1}' output/zvo_energy.dat

write_common_defs
write_valid_transsym
write_base_namelist
perl -0pi -e 's/^2Sz 0$/Ndown 2\nNup 2/m' modpara.def
rm -rf output
run_hphi ndown_nup.log ../../src/HPhi -e namelist.def
grep -q "Symmetry basis: raw_dim=6 sector_dim=2 group_order=4" ndown_nup.log
awk '$1 == "Sz" {found=1; d=$2; if(d<0)d=-d; ok=(d < 0.000001)} END{exit found && ok ? 0 : 1}' output/zvo_energy.dat

write_common_defs
write_valid_transsym
write_base_namelist
perl -0pi -e 's/^2Sz 0$/Nup 3\nNdown 1\n2Sz 0/m' modpara.def
rm -rf output
expect_fail nup_ndown_2sz_conflict.log ../../src/HPhi -e namelist.def
grep -q "conflicts with Nup-Ndown" nup_ndown_2sz_conflict.log

write_common_defs
write_valid_transsym
write_base_namelist
perl -0pi -e 's/^2Sz 0$/Ndown 4/m' modpara.def
rm -rf output
expect_fail ndown_without_nup.log ../../src/HPhi -e namelist.def
grep -q "Nup and Ndown must be specified together" ndown_without_nup.log

write_common_defs
write_c4_kpi2_transsym
write_base_namelist
perl -0pi -e 's/CalcType 0/CalcType 3/' calcmod.def
perl -0pi -e 's/^exct 1$/exct 2/m' modpara.def
if ../../src/HPhi -e namelist.def > exct_too_large.log 2>&1; then
    cat exct_too_large.log
    exit 1
fi
grep -q "smaller than exct" exct_too_large.log

write_common_defs
write_valid_transsym
write_base_namelist
perl -0pi -e 's/^LanczosTarget 1$/LanczosTarget 2/m' modpara.def
rm -rf output
run_hphi target_clamped.log ../../src/HPhi -e namelist.def
grep -q "LanczosTarget=2 is outside Hilbert dimension 2; use 1" target_clamped.log
awk '$1 == "Sz" {found=1; d=$2; if(d<0)d=-d; ok=(d < 0.000001)} END{exit found && ok ? 0 : 1}' output/zvo_energy.dat

exit 0
