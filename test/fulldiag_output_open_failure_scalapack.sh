#!/bin/sh -e

# Reuse the backend-neutral failure/synchronization scenarios with Solver 1.
exec "$(dirname "$0")/fulldiag_output_open_failure.sh"
