#!/bin/sh -e
set -e

# "# MomentumIndex m" metadata in the TransSym file (qptransidx.def).
# Comment lines are skipped by the reader, so the metadata never changes the
# calculation; it is scanned separately and reported in the log as
# "TransSym metadata: MomentumIndex=m". Malformed metadata is rejected.

mkdir -p symmetry_transsym_metadata
cd symmetry_transsym_metadata

run_hphi() {
    log="$1"
    shift
    "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

expect_failure() {
    pattern="$1"
    log="$2"
    shift 2
    if "$@" > "${log}" 2>&1; then
        cat "${log}"
        echo "Expected failure did not happen: ${pattern}"
        exit 1
    fi
    if ! grep -q "${pattern}" "${log}"; then
        cat "${log}"
        echo "Expected pattern not found: ${pattern}"
        exit 1
    fi
}

energy() {
    awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat
}

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

cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
TransSym qptransidx.def
EOF

# Translations of the 6-site ring in the k = 0 sector (all characters 1),
# without any metadata.
write_k0_body() {
    awk 'BEGIN {
        L = 6
        print "============================================="
        printf "NQPTrans %10d\n", L
        print "============================================="
        print "======== TrIdx_TrWeight_and_TrIdx_i_xi ======"
        print "============================================="
        for (op = 0; op < L; op++) printf "%d 1.0 0.0\n", op
        for (op = 0; op < L; op++)
            for (site = 0; site < L; site++)
                printf "%d %d %d 1\n", op, site, (site + op) % L
    }'
}

# write_with_metadata <first line> [<last line>]: the body with a comment
# line before it, a second comment and an empty line, and optionally a
# trailing line.
write_with_metadata() {
    {
        printf "%s\n" "$1"
        printf "# generated for the metadata test\n\n"
        write_k0_body
        if [ -n "$2" ]; then printf "%s\n" "$2"; fi
    } > qptransidx.def
}

write_k0_body > qptransidx.def
rm -rf output
run_hphi no_metadata.log ../../src/HPhi -e namelist.def
ref_energy=`energy`
test -n "${ref_energy}"
grep -q "Symmetry basis: raw_dim=20 sector_dim=4 group_order=6" no_metadata.log
if grep -q "TransSym metadata" no_metadata.log; then
    cat no_metadata.log
    echo "Metadata was reported although the file has none."
    exit 1
fi

# Accepted spellings. The metadata must not change the result.
check_accepted() {
    label="$1"
    expected="$2"
    shift 2
    rm -rf output
    run_hphi "${label}.log" "$@" ../../src/HPhi -e namelist.def
    grep -q "TransSym metadata: MomentumIndex=${expected}$" "${label}.log"
    grep -q "Symmetry basis: raw_dim=20 sector_dim=4 group_order=6" "${label}.log"
    this_energy=`energy`
    test -n "${this_energy}"
    diff=`awk -v a="${this_energy}" -v b="${ref_energy}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
    test "${diff}" = "0.000000"
}

write_with_metadata "# MomentumIndex 2"
check_accepted first_line 2
write_with_metadata "#MomentumIndex 2"
check_accepted no_space 2
write_with_metadata "#   momentumindex 2"
check_accepted lower_case 2
write_with_metadata "# a comment that is not metadata" "# MomentumIndex 2"
check_accepted last_line 2
write_with_metadata "# MomentumIndex 2" "# MomentumIndex 2"
check_accepted repeated_same_value 2
write_with_metadata "# MomentumIndex +2  "
check_accepted explicit_plus 2

# Check both int bounds and values that overflow even a 64-bit long.
# The metadata is descriptive: INT_MAX need not be a valid momentum sector.
check_integer_boundaries() {
    prefix="$1"
    shift
    for value in 0 2147483647; do
        write_with_metadata "# MomentumIndex ${value}"
        check_accepted "${prefix}_boundary_${value}" "${value}" "$@"
    done
    for value in 2147483648 4294967296 -4294967296 -4294967294 \
                 9223372036854775808 -9223372036854775809; do
        write_with_metadata "# MomentumIndex ${value}"
        expect_failure "TransSym metadata must be" "${prefix}_out_of_range_${value}.log" \
            "$@" ../../src/HPhi -e namelist.def
    done
}

check_integer_boundaries serial

if [ -n "${MPIRUN}" ]; then
    MPI_NP=`printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}'`
    if printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$" && [ "${MPI_NP}" -gt 1 ]; then
        write_with_metadata "# MomentumIndex 2"
        check_accepted mpi 2 ${MPIRUN}
        check_integer_boundaries mpi ${MPIRUN}
        write_with_metadata "# MomentumIndex two"
        rm -rf output
        expect_failure "TransSym metadata must be" mpi_malformed.log ${MPIRUN} ../../src/HPhi -e namelist.def
    fi
fi

# Rejected metadata.
write_with_metadata "# MomentumIndex two"
expect_failure "TransSym metadata must be" malformed_value.log ../../src/HPhi -e namelist.def
write_with_metadata "# MomentumIndex"
expect_failure "TransSym metadata must be" missing_value.log ../../src/HPhi -e namelist.def
write_with_metadata "# MomentumIndex 2.5"
expect_failure "TransSym metadata must be" fractional_value.log ../../src/HPhi -e namelist.def
write_with_metadata "# MomentumIndex 2x"
expect_failure "TransSym metadata must be" suffixed_value.log ../../src/HPhi -e namelist.def
write_with_metadata "# MomentumIndex=2"
expect_failure "TransSym metadata must be" malformed_key.log ../../src/HPhi -e namelist.def
write_with_metadata "# MomentumIndex -1"
expect_failure "TransSym metadata must be" negative.log ../../src/HPhi -e namelist.def
write_with_metadata "# MomentumIndex 2 3"
expect_failure "TransSym metadata must be" trailing.log ../../src/HPhi -e namelist.def
write_with_metadata "# MomentumIndex 2" "# MomentumIndex 3"
expect_failure "given twice with different values" conflicting.log ../../src/HPhi -e namelist.def

exit 0
