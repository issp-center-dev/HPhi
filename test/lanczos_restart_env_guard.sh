#!/bin/sh
set -eu
python3 "$(dirname "$0")/lanczos_restart_env_guard.py" ../src/HPhi
