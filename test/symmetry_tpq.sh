#!/bin/sh
set -eu
python3 "$(dirname "$0")/symmetry_tpq.py" ../src/HPhi ./unittest_tpq_failure
