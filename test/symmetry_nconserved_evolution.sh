#!/bin/sh
set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec python3 "$script_dir/symmetry_nconserved_evolution.py" ../src/unittest_symmetry_sector_probe
