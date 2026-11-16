#!/bin/sh
set -eu
python3 "$(dirname "$0")/tpq_legacy_restart.py" ../src/HPhi
