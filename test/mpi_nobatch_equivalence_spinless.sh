#!/bin/sh -e
# Assert HPHI_MPI_NOBATCH=1 (per-term MPI) reproduces the batched result for the
# SpinlessFermion / SpinlessFermionGC transfer MPIsingle batched path. The
# existing mpi_nobatch_equivalence test excludes the Expert-mode spinless models;
# they are only covered serial-vs-MPI, which cannot isolate a rank-uniform
# batched bug. This compares batched vs HPHI_MPI_NOBATCH=1 under the SAME MPIRUN
# for both the energy and the one-body Green function.
# $1 = CMAKE_SOURCE_DIR (for test/testSpinlessCalc.py).
if [ -z "${MPIRUN}" ]; then echo "MPIRUN not set. Skipping."; exit 0; fi

SRCDIR="$1"
mkdir -p mpi_nobatch_equivalence_spinless
cd mpi_nobatch_equivalence_spinless

spinless_case() {  # $1 = tag, $2... = testSpinlessCalc.py args
  tag="$1"; shift
  python3 "${SRCDIR}/test/testSpinlessCalc.py" -p "../../src/HPhi" "$@" --onebody-offdiag > "log_${tag}_gen.txt" 2>&1
  rm -rf output
  ${MPIRUN} ../../src/HPhi -e namelist.def > "log_${tag}_batched.txt" 2>&1
  # Confirm the batched Spinless MPIsingle path actually fired.
  if ! grep -q "Spinless:.*MPIsingle transfers -> .* groups" "log_${tag}_batched.txt"; then
    echo "[${tag}] ERROR: batched Spinless MPIsingle path did not fire"
    grep -i "MPI Batching\|MPI batching" "log_${tag}_batched.txt" || true; exit 1
  fi
  cp output/zvo_energy.dat "energy_${tag}_batched.dat"
  cp output/zvo_cisajs.dat "cisajs_${tag}_batched.dat"
  rm -rf output
  HPHI_MPI_NOBATCH=1 ${MPIRUN} ../../src/HPhi -e namelist.def > "log_${tag}_nobatch.txt" 2>&1
  cp output/zvo_energy.dat "energy_${tag}_nobatch.dat"
  cp output/zvo_cisajs.dat "cisajs_${tag}_nobatch.dat"
  de=$(paste "energy_${tag}_batched.dat" "energy_${tag}_nobatch.dat" \
       | awk '$2 ~ /^[-+0-9.eE]+$/ && $4 ~ /^[-+0-9.eE]+$/ { x=$2-$4; if(x<0)x=-x; if(x>m)m=x } END{ printf "%.3e", m+0 }')
  dg=$(paste "cisajs_${tag}_batched.dat" "cisajs_${tag}_nobatch.dat" \
       | awk '{ dre=$5-$11; dim=$6-$12; d=sqrt(dre*dre+dim*dim); if(d>m)m=d } END{ printf "%.3e", m+0 }')
  echo "[${tag}] max |batched - nobatch| : energy = ${de} , onebody-Green = ${dg}"
  awk -v de="${de}" -v dg="${dg}" 'BEGIN{ exit (de < 1e-10 && dg < 1e-10) ? 0 : 1 }' \
    || { echo "[${tag}] MISMATCH"; exit 1; }
}

spinless_case spinless    -m SpinlessFermion   -s 8
spinless_case spinless_gc -m SpinlessFermionGC -s 8 -V 0.5

echo "SpinlessFermion(GC): batched == no-batch within tolerance (energy + one-body Green)."
