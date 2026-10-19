#!/bin/sh
set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec python3 "$script_dir/symmetry_fermion_gc_io_restart.py" ../src/HPhi ../src/unittest_symmetry_sector_probe
