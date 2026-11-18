#!/bin/sh
set -eu
exec python3 "$(dirname "$0")/observable_output.py" ../src/HPhi spinless_number
