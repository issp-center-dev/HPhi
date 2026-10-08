#!/bin/sh
set -e

# Regression test for exitMPI(): when every rank aborts on the same input
# error, the diagnostic that only rank 0 prints must survive even if rank 0
# is the slowest rank to reach exitMPI(). The driver delays rank 0 before it
# prints the marker (see unittest_exitmpi_abort_grace.c); without the grace
# period in exitMPI() the launcher kills rank 0 first and the marker is lost.

if [ -z "${MPIRUN}" ]; then
    echo "Error: MPIRUN is not set."
    exit 2
fi

if [ ! -x "./unittest_exitmpi_abort_grace" ]; then
    echo "Error: unittest_exitmpi_abort_grace is not available in $(pwd)."
    exit 2
fi

testname="exitmpi_abort_grace_mpi"
bin="$(pwd)/unittest_exitmpi_abort_grace"
marker="exitMPI abort grace marker"

rm -rf "${testname}"
mkdir -p "${testname}"
cd "${testname}"

for trial in 1 2 3; do
    set +e
    ${MPIRUN} "${bin}" > "run_${trial}.log" 2>&1
    rc=$?
    set -e
    if [ "${rc}" = "0" ]; then
        echo "[trial ${trial}] ERROR: run succeeded but every rank was expected to abort" >&2
        tail -40 "run_${trial}.log" >&2
        exit 1
    fi
    if ! grep -q "${marker}" "run_${trial}.log"; then
        echo "[trial ${trial}] ERROR: the rank-0 diagnostic '${marker}' was lost" >&2
        tail -40 "run_${trial}.log" >&2
        exit 1
    fi
    echo "[trial ${trial}] rank-0 diagnostic present (launcher exit code ${rc})"
done
