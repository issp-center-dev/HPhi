#!/bin/sh
set -eu
python3 "$(dirname "$0")/namelist_parse_validation.py" ../src/HPhi
