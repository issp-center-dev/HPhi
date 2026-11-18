#!/bin/sh
set -eu
# These diagnostics concern selection/recovery, not MPI redistribution.
# Use one process even in the MPI baseline; collective cases have dedicated tests.
export MPIRUN=''
exec python3 "$(dirname "$0")/restart_contract_diagnostics.py" ../src/HPhi
