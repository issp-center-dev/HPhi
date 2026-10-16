#!/bin/sh
set -eu
PYTHONPATH="$(dirname "$0")" python3 "$(dirname "$0")/restart_input_validation.py" ../src/HPhi
