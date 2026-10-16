#!/bin/sh
# Regression test for run_with_mpi_precheck.sh (min:N / exact:N / pow:N).
# $1 = CMAKE_SOURCE_DIR.  No HPhi or MPI process is launched: the precheck
# wraps a dummy script and only its exit code is checked.
set -u
srcdir="$1"
precheck="${srcdir}/test/run_with_mpi_precheck.sh"

mkdir -p precheck_modes
cd precheck_modes
printf '#!/bin/sh\necho "dummy ran: $*"\n' > dummy.sh
chmod +x dummy.sh

fail=0
# check <expected rc> <mode> <MPIRUN value or "unset">
check() {
  expected="$1"; mode="$2"; np="$3"
  if [ "${np}" = "unset" ]; then
    env -u MPIRUN sh "${precheck}" "${mode}" ./dummy.sh > out.txt 2>&1
  else
    MPIRUN="mpirun -np ${np}" sh "${precheck}" "${mode}" ./dummy.sh > out.txt 2>&1
  fi
  rc=$?
  if [ "${rc}" -ne "${expected}" ]; then
    echo "FAIL: mode=${mode} MPIRUN=${np}: expected rc=${expected}, got rc=${rc}"
    cat out.txt
    fail=1
  else
    echo "ok  : mode=${mode} MPIRUN=${np} -> rc=${rc}"
  fi
}

# min:N
check 77 min:4 unset
check 77 min:4 2
check 0  min:4 4
check 0  min:4 16
# exact:N
check 77 exact:4 unset
check 77 exact:4 2
check 0  exact:4 4
check 77 exact:4 16
# pow:N
check 0  pow:3 unset
check 0  pow:3 1
check 77 pow:3 2
check 0  pow:3 3
check 77 pow:3 4
check 0  pow:3 9
check 77 pow:3 16
check 0  pow:3 27
check 0  pow:10 100
check 0  pow:16 256
check 77 pow:16 64
# leading zeros must not be read as octal
check 0  pow:3 09
check 0  pow:3 027
check 0  exact:4 04
# invalid modes are rejected even when MPIRUN is unset
check 2  pow:1 unset
check 2  pow:1 4
check 2  pow:abc unset
check 2  min:abc 4
check 2  bogus:4 4
check 2  pow 4
# unparsable MPIRUN is a skip
env MPIRUN="srun" sh "${precheck}" exact:4 ./dummy.sh > out.txt 2>&1
rc=$?
if [ "${rc}" -ne 77 ]; then echo "FAIL: unparsable MPIRUN: expected rc=77, got ${rc}"; cat out.txt; fail=1; else echo "ok  : unparsable MPIRUN -> rc=77"; fi
# srun-style "-n"
MPIRUN="srun -n 9" sh "${precheck}" pow:3 ./dummy.sh > out.txt 2>&1
rc=$?
if [ "${rc}" -ne 0 ]; then echo "FAIL: srun -n 9 with pow:3: expected rc=0, got ${rc}"; cat out.txt; fail=1; else echo "ok  : srun -n 9 pow:3 -> rc=0"; fi

exit ${fail}
