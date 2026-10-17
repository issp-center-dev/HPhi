#!/bin/sh
set -eu
PYTHONPATH="$(dirname "$0")" python3 "$(dirname "$0")/restart_tpq_contract.py" ../src/HPhi "${1:-tpq}"
