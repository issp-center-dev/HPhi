#!/bin/sh -e
# Regression test for the modulo-by-zero (% 0) UB in CalcByTEM.c.
#
# The TE progress printout used `step_i % (Lanczos_max / 10)`. For
# Lanczos_max < 10 the divisor is 0, so the modulo is undefined behavior:
# on x86 it raises SIGFPE and aborts the run. (On arm64 integer % 0 happens
# not to trap, so this crash guard is effective on the x86 CI jobs.) The fix
# clamps the print interval to >= 1. This runs a TE with Lanczos_max=5 and
# asserts the run completes instead of crashing.
# $1 = CMAKE_SOURCE_DIR (for test/testTECalc.py).

# If MPIRUN is unset, run in serial mode.
if [ -z "${MPIRUN}" ]; then
    MPIRUN=""
fi

SRCDIR="$1"
HPHI=../../src/HPhi

mkdir -p te_small_lanczos_max
cd te_small_lanczos_max

# Generate a working Hubbard TE setup (also produces the input eigenvector).
python3 "${SRCDIR}/test/testTECalc.py" -p "${HPHI}" -mpi "${MPIRUN}" > gen.log 2>&1

# Re-run the time-evolution step with a small Lanczos_max (< 10) that makes
# the buggy `Lanczos_max / 10` divisor zero. Edit the ModPara file that
# namelist2.def actually references (testTECalc.py leaves it as modpara.def).
MODPARA=$(awk '$1=="ModPara"{print $2}' namelist2.def)
sed -e 's/^Lanczos_max.*/Lanczos_max    5/' "${MODPARA}" > _t && mv _t "${MODPARA}"

if ! ${MPIRUN} "${HPHI}" -e namelist2.def > te_small.log 2>&1; then
    echo "ERROR: TE with Lanczos_max=5 failed (modulo-by-zero regression?)" >&2
    tail -40 te_small.log >&2
    exit 1
fi

echo "TE Lanczos_max=5 (<10) completes without modulo-by-zero: OK"
