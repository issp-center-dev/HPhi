#!/bin/sh
set -e

if [ -z "${MPIRUN}" ]; then
    echo "Error: MPIRUN is not set."
    exit 2
fi

if [ ! -x "./unittest_symmetry_distribution_warning" ]; then
    echo "Error: unittest_symmetry_distribution_warning is not available in $(pwd)."
    exit 2
fi

if ! output=$(${MPIRUN} ./unittest_symmetry_distribution_warning --mpi 2>&1); then
    printf "%s\n" "${output}"
    exit 1
fi
printf "%s\n" "${output}"
printf "%s\n" "${output}" |
    grep -F "exceeds warning threshold=2; continuing with one sample per nonempty MPI rank." >/dev/null
