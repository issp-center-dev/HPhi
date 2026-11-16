#!/bin/sh
set -eu
python3 "$(dirname "$0")/restart_te_contract.py" ../src/HPhi
