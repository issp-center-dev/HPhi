#!/bin/sh
set -eu
python3 "$(dirname "$0")/lanczos_restart_tm_sequence.py" ../src/HPhi
