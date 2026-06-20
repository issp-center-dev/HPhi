#!/bin/sh
set -e

testname="mpi_nbody_interall_spinless"
hphi="../../../src/HPhi"
SRCDIR="$1"
tol="0.00000001"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

make_case() {
  dir="$1"
  model="$2"
  mkdir -p "${dir}"
  cd "${dir}"
  if [ "${model}" = "SpinlessFermion" ]; then
    python3 "${SRCDIR}/test/testSpinlessCalc.py" -p /bin/true \
      -m "${model}" -s 8 -n 3 --hopping 0.0 > log_generate.txt 2>&1
  else
    python3 "${SRCDIR}/test/testSpinlessCalc.py" -p /bin/true \
      -m "${model}" -s 8 --hopping 0.0 > log_generate.txt 2>&1
  fi
  printf '    NBodyInterAll  nbodyinterall.def\n' >> namelist.def
  cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 3
========================
========NBodyInterAll===
========================
1 6 0 6 0 0.1250000000000000 0.0000000000000000
1 7 0 0 0 0.3700000000000000 0.1100000000000000
1 0 0 7 0 0.3700000000000000 -0.1100000000000000
EOF
  cd ..
}

run_case() {
  tag="$1"
  model="$2"

  rm -rf "serial_${tag}" "mpi_${tag}"
  make_case "serial_${tag}" "${model}"

  cd "serial_${tag}"
  run_hphi log_serial.txt "${hphi}" -e namelist.def
  cp output/zvo_energy.dat "../energy_serial_${tag}.dat"
  cd ..

  mkdir -p "mpi_${tag}"
  cp "serial_${tag}"/*.def "mpi_${tag}/"
  cd "mpi_${tag}"
  run_hphi log_mpi.txt ${MPIRUN} "${hphi}" -e namelist.def
  cp output/zvo_energy.dat "../energy_mpi_${tag}.dat"
  cd ..

  diff=$(paste "energy_serial_${tag}.dat" "energy_mpi_${tag}.dat" \
    | awk '$1 == "Energy" && $3 == "Energy" {d=$2-$4; if(d<0)d=-d; if(d>m)m=d} END{printf "%.12g", m+0}')
  awk -v d="${diff}" -v t="${tol}" 'BEGIN{exit (d < t) ? 0 : 1}' || {
    echo "${model} serial/MPI NBodyInterAll energy mismatch: max diff ${diff}"
    paste "energy_serial_${tag}.dat" "energy_mpi_${tag}.dat"
    exit 1
  }

  grep -q "INTER process site" "mpi_${tag}/log_mpi.txt" || {
    echo "${model} MPI run did not print an inter-process site summary"
    cat "mpi_${tag}/log_mpi.txt"
    exit 1
  }
}

run_case spinless SpinlessFermion
run_case spinlessgc SpinlessFermionGC

echo "Spinless NBodyInterAll serial/MPI energies match for inter-process factors."
