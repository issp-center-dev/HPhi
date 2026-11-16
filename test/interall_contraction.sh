#!/bin/sh
set -eu
exec python3 "$(dirname "$0")/interall_contraction.py" ../src/HPhi
