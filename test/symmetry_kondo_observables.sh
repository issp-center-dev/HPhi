#!/bin/sh
set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
python3 "$script_dir/symmetry_kondo_observables.py" \
  ../src/HPhi ../src/unittest_symmetry_sector_probe
