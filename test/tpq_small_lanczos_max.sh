#!/bin/sh -e
# Regression test for the modulo-by-zero (% 0) UB in CalcByTPQ.c.
#
# The mTPQ progress printout used `step_i % ((Lanczos_max - step_iO) / 10)`.
# When (Lanczos_max - step_iO) < 10 the divisor is 0, so the modulo is
# undefined behavior: on x86 it raises SIGFPE and aborts the run. (On arm64
# integer % 0 happens not to trap, so this crash guard is effective on the
# x86 CI jobs.) The fix clamps the print interval to >= 1. This runs mTPQ with
# lanczos_max=5 and asserts the run completes instead of crashing.

# If MPIRUN is unset, run in serial mode.
if [ -z "${MPIRUN}" ]; then
    MPIRUN=""
fi

mkdir -p tpq_small_lanczos_max
cd tpq_small_lanczos_max

# L = 8 (not 4): the baseline CI runs this at mpisize=16, and the MPI site
# partition for the Spin model needs enough sites for nproc=2^k=16. L=4 aborts
# with "The number of PROCESS should be 2-exponent ...".
cat > stan.in <<EOF
L = 8
model = "Spin"
method = "TPQ"
lattice = "chain"
J = 1.0
2Sz = 0
lanczos_max = 5
LargeValue = 5
NumAve = 1
EOF

if ! ${MPIRUN} ../../src/HPhi -s stan.in > tpq_small.log 2>&1; then
    echo "ERROR: mTPQ with lanczos_max=5 failed (modulo-by-zero regression?)" >&2
    tail -40 tpq_small.log >&2
    exit 1
fi

echo "mTPQ lanczos_max=5 (<10) completes without modulo-by-zero: OK"
