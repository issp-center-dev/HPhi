#!/bin/sh
set -eu
python3 "$(dirname "$0")/symmetry_ctpq.py" ../src/HPhi ./unittest_tpq_failure
