#!/bin/sh
set -eu
PYTHONPATH="$(dirname "$0")${PYTHONPATH:+:$PYTHONPATH}" \
    python3 "$(dirname "$0")/restart_publish_failures.py" ../src/HPhi
