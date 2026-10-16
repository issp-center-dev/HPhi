#!/bin/sh
set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec python3 "$script_dir/symmetry_kondo_basis.py" ../src/HPhi ../src/unittest_symmetry_sector_probe
