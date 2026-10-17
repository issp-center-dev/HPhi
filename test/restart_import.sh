#!/bin/sh
set -eu
python3 "$(dirname "$0")/restart_import.py" ../src/HPhi
