#!/bin/sh
set -eu
exec python3 "$(dirname "$0")/tpq_random_streams.py" ./unittest_tpq_failure
