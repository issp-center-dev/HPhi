#!/bin/sh
set -eu
exec python3 "$(dirname "$0")/kondo_mpi_dimension.py" ../src/HPhi
