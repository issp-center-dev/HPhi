#!/bin/sh
set -eu
exec python3 "$(dirname "$0")/observable_output.py" ../src/HPhi tj_sz ../src/unittest_symmetry_sector_probe
