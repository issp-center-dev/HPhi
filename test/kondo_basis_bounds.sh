#!/bin/sh
set -eu
exec python3 "$(dirname "$0")/kondo_basis_bounds.py" ../src/HPhi serial
