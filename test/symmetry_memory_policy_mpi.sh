#!/bin/sh
set -e

if [ -z "${MPIRUN}" ]; then
    echo "Error: MPIRUN is not set."
    exit 2
fi

if [ ! -x "./unittest_symmetry_memory_policy" ]; then
    echo "Error: unittest_symmetry_memory_policy is not available in $(pwd)."
    exit 2
fi

${MPIRUN} ./unittest_symmetry_memory_policy --mpi
