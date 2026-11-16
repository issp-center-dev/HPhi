#!/bin/sh
set -eu
PYTHONPATH="$(dirname "$0")${PYTHONPATH:+:$PYTHONPATH}" python3 "$(dirname "$0")/restart_cg_contract.py" ../src/HPhi
