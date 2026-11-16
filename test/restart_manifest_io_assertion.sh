#!/bin/sh
set -u

output=$("$@" 2>&1)
status=$?
printf '%s\n' "$output"

if [ "$status" -eq 0 ]; then
    echo "intentional MPI assertion unexpectedly succeeded" >&2
    exit 1
fi
case "$output" in
    *"intentional rank-local assertion"*) ;;
    *)
        echo "intentional MPI assertion marker was not reported" >&2
        exit 1
        ;;
esac
