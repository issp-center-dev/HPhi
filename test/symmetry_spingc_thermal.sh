#!/bin/sh
set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
python3 "$script_dir/symmetry_spingc_thermal.py" \
  ../src/HPhi ./unittest_tpq_failure
