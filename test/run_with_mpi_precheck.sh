#!/bin/sh
set -e

SKIP_RC=77

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <min:N|exact:N|pow:N> <test_script> [args...]"
    echo "  min:N   : run if MPI ranks >= N, skip otherwise (skip if MPIRUN unset)"
    echo "  exact:N : run if MPI ranks == N, skip otherwise (skip if MPIRUN unset)"
    echo "  pow:N   : run if MPI ranks is a power of N (N^k, k>=0), skip otherwise"
    echo "            (run in serial if MPIRUN unset)"
    exit 2
fi

mode="$1"
shift
test_script="$1"
shift

# Strip leading zeros so that the value is never read as octal by $((...)).
to_decimal() {
    d=$(printf "%s" "$1" | sed 's/^0*//')
    [ -z "${d}" ] && d=0
    printf "%s" "${d}"
}

# Validate the mode first, independently of MPIRUN.
kind="${mode%%:*}"
req="${mode#*:}"
case "${kind}" in
    min|exact|pow) ;;
    *)
        echo "Internal error: unknown mode '${mode}'."
        exit 2
        ;;
esac
if [ "${req}" = "${mode}" ] || ! printf "%s\n" "${req}" | grep -Eq "^[0-9]+$"; then
    echo "Internal error: invalid mode '${mode}'."
    exit 2
fi
req=$(to_decimal "${req}")
if [ "${kind}" = "pow" ] && [ "${req}" -lt 2 ]; then
    echo "Internal error: invalid mode '${mode}' (base must be >= 2)."
    exit 2
fi

if [ -z "${MPIRUN}" ]; then
    if [ "${kind}" = "pow" ]; then
        # Serial run is always allowed for pow:N tests.
        exec "${test_script}" "$@"
    fi
    echo "Skipping (MPI-only test): MPIRUN is not set (required by $(basename "${test_script}"))."
    exit "${SKIP_RC}"
fi

MPI_NP=$(printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}')
if ! printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$"; then
    echo "Skipping (MPI-only test): Could not parse -np/-n from MPIRUN='${MPIRUN}'."
    exit "${SKIP_RC}"
fi
MPI_NP=$(to_decimal "${MPI_NP}")

case "${kind}" in
    min)
        if [ "${MPI_NP}" -lt "${req}" ]; then
            echo "Skipping (MPI-only test): requires MPI ranks >= ${req} (current: ${MPI_NP})."
            exit "${SKIP_RC}"
        fi
        ;;
    exact)
        if [ "${MPI_NP}" -ne "${req}" ]; then
            echo "Skipping (MPI-only test): requires MPI ranks = ${req} (current: ${MPI_NP})."
            exit "${SKIP_RC}"
        fi
        ;;
    pow)
        n="${MPI_NP}"
        while [ "${n}" -gt 1 ] && [ $((n % req)) -eq 0 ]; do
            n=$((n / req))
        done
        if [ "${n}" -ne 1 ]; then
            echo "Skipping (MPI rank constraint): requires MPI ranks = ${req}^k (current: ${MPI_NP})."
            exit "${SKIP_RC}"
        fi
        ;;
esac

exec "${test_script}" "$@"
