#!/bin/sh
set -eu
PYTHONPATH="$(dirname "$0")${PYTHONPATH:+:$PYTHONPATH}" python3 "$(dirname "$0")/restart_cg_noop.py" ../src/HPhi
