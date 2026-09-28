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

if [ -z "${MPIRUN}" ]; then
    case "${mode}" in
        pow:*)
            # Serial run is always allowed for pow:N tests.
            exec "${test_script}" "$@"
            ;;
    esac
    echo "Skipping (MPI-only test): MPIRUN is not set (required by $(basename "${test_script}"))."
    exit "${SKIP_RC}"
fi

MPI_NP=$(printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}')
if ! printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$"; then
    echo "Skipping (MPI-only test): Could not parse -np/-n from MPIRUN='${MPIRUN}'."
    exit "${SKIP_RC}"
fi

case "${mode}" in
    min:*)
        req="${mode#min:}"
        if ! printf "%s\n" "${req}" | grep -Eq "^[0-9]+$"; then
            echo "Internal error: invalid mode '${mode}'."
            exit 2
        fi
        if [ "${MPI_NP}" -lt "${req}" ]; then
            echo "Skipping (MPI-only test): requires MPI ranks >= ${req} (current: ${MPI_NP})."
            exit "${SKIP_RC}"
        fi
        ;;
    exact:*)
        req="${mode#exact:}"
        if ! printf "%s\n" "${req}" | grep -Eq "^[0-9]+$"; then
            echo "Internal error: invalid mode '${mode}'."
            exit 2
        fi
        if [ "${MPI_NP}" -ne "${req}" ]; then
            echo "Skipping (MPI-only test): requires MPI ranks = ${req} (current: ${MPI_NP})."
            exit "${SKIP_RC}"
        fi
        ;;
    pow:*)
        req="${mode#pow:}"
        if ! printf "%s\n" "${req}" | grep -Eq "^[2-9][0-9]*$"; then
            echo "Internal error: invalid mode '${mode}'."
            exit 2
        fi
        n="${MPI_NP}"
        while [ "${n}" -gt 1 ] && [ $((n % req)) -eq 0 ]; do
            n=$((n / req))
        done
        if [ "${n}" -ne 1 ]; then
            echo "Skipping (MPI rank constraint): requires MPI ranks = ${req}^k (current: ${MPI_NP})."
            exit "${SKIP_RC}"
        fi
        ;;
    *)
        echo "Internal error: unknown mode '${mode}'."
        exit 2
        ;;
esac

exec "${test_script}" "$@"
