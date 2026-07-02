#!/bin/sh -e
# Regression test for missing restart-vector fallback in mTPQ and cTPQ.
#
# Restart_in without tmpvec_set*_rank_*.dat should fall back to a fresh random
# initial vector.  This used to dereference a NULL FILE pointer in the TPQ
# restart path.

if [ -z "${MPIRUN}" ]; then
    MPIRUN=""
fi

testname="tpq_restart_missing_vector_fallback"
rm -rf "${testname}"
mkdir -p "${testname}"
cd "${testname}"
hphi="$(pwd)/../../src/HPhi"

run_mtpq_missing_restart() {
    mkdir -p mtpq
    cd mtpq

    cat > stan.in <<EOF
L = 8
model = "Spin"
method = "TPQ"
lattice = "chain"
J = 1.0
2Sz = 0
Lanczos_max = 3
LargeValue = 5
NumAve = 1
Restart = "Restart_in"
outputmode = "None"
EOF

    if ! ${MPIRUN} "${hphi}" -s stan.in > run.log 2>&1; then
        echo "ERROR: mTPQ missing restart-vector fallback failed" >&2
        tail -40 run.log >&2
        exit 1
    fi
    grep -q "Start to calculate in normal procedure." run.log
    test "$(wc -l < output/Norm_rand0.dat)" -gt 1

    cd ..
}

run_ctpq_missing_restart() {
    mkdir -p ctpq
    cd ctpq

    cat > stan.in <<EOF
L = 8
model = "SpinGC"
method = "TPQ"
lattice = "chain"
J = 1.0
NumAve = 1
Lanczos_max = 3
LargeValue = 5
outputmode = "None"
EOF

    "${hphi}" -sdry stan.in > gen.log 2>&1

    cat > calcmod.def <<EOF
CalcType   5
CalcModel   4
ReStart   3
CalcSpec   0
CalcEigenVec   0
InitialVecType   0
InputEigenVec   0
OutputEigenVec   0
InputHam   0
OutputHam   0
OutputExVec   0
EOF

    if ! ${MPIRUN} "${hphi}" -e namelist.def > run.log 2>&1; then
        echo "ERROR: cTPQ missing restart-vector fallback failed" >&2
        tail -40 run.log >&2
        exit 1
    fi
    grep -q "Start to calculate in normal procedure." run.log
    test "$(wc -l < output/Norm_rand0.dat)" -gt 1

    cd ..
}

run_mtpq_missing_restart
run_ctpq_missing_restart

echo "mTPQ/cTPQ missing restart-vector fallback completes: OK"
