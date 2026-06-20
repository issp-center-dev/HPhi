#!/bin/sh
set -e

testname="mpi_nbodyg_spinless"
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
      -m "${model}" -s 8 -n 3 -V 0.5 --onebody-offdiag --offdiag > log_generate.txt 2>&1
  else
    python3 "${SRCDIR}/test/testSpinlessCalc.py" -p /bin/true \
      -m "${model}" -s 8 -V 0.5 --onebody-offdiag --offdiag > log_generate.txt 2>&1
  fi
  printf '    NBodyG  nbodyg.def\n' >> namelist.def
  cat > nbodyg.def <<EOF
========================
NNBodyG 3
========================
========NBodyG==========
========================
1 6 0 6 0
1 7 0 0 0
2 7 0 0 0 6 0 6 0
EOF
  cd ..
}

compare_outputs() {
  serial_file="$1"
  mpi_file="$2"
  label="$3"
  awk -v t="${tol}" '
    function prefix(line, out) {
      out = line;
      sub(/[[:space:]]+[-+0-9.eE]+[[:space:]]+[-+0-9.eE]+$/, "", out);
      return out;
    }
    NR == FNR {
      s_prefix[NR] = prefix($0);
      s_re[NR] = $(NF - 1);
      s_im[NR] = $NF;
      n = NR;
      next;
    }
    {
      m = FNR;
      if (prefix($0) != s_prefix[m]) bad = 1;
      dr = $(NF - 1) - s_re[m];
      di = $NF - s_im[m];
      if (dr < 0) dr = -dr;
      if (di < 0) di = -di;
      if (dr >= t || di >= t) bad = 1;
    }
    END { exit (n == m && bad != 1) ? 0 : 1; }
  ' "${serial_file}" "${mpi_file}" || {
    echo "${label} serial/MPI NBodyG output mismatch"
    paste "${serial_file}" "${mpi_file}"
    exit 1
  }
}

run_case() {
  tag="$1"
  model="$2"

  rm -rf "serial_${tag}" "mpi_${tag}"
  make_case "serial_${tag}" "${model}"

  cd "serial_${tag}"
  run_hphi log_serial.txt "${hphi}" -e namelist.def
  cp output/zvo_NBodyG.dat "../nbodyg_serial_${tag}.dat"
  cd ..

  mkdir -p "mpi_${tag}"
  cp "serial_${tag}"/*.def "mpi_${tag}/"
  cd "mpi_${tag}"
  run_hphi log_mpi.txt ${MPIRUN} "${hphi}" -e namelist.def
  cp output/zvo_NBodyG.dat "../nbodyg_mpi_${tag}.dat"
  cd ..

  compare_outputs "nbodyg_serial_${tag}.dat" "nbodyg_mpi_${tag}.dat" "${model}"

  grep -q "INTER process site" "mpi_${tag}/log_mpi.txt" || {
    echo "${model} MPI run did not print an inter-process site summary"
    cat "mpi_${tag}/log_mpi.txt"
    exit 1
  }

  awk -v t="${tol}" '
    $1 == 1 && $2 == 6 && $3 == 0 && $4 == 6 && $5 == 0 {
      re = $6; if (re < 0) re = -re;
      found = (re > t);
    }
    END { exit found ? 0 : 1; }
  ' "nbodyg_mpi_${tag}.dat" || {
    echo "${model} inter-process diagonal NBodyG operator was zero or missing"
    cat "nbodyg_mpi_${tag}.dat"
    exit 1
  }
}

run_case spinless SpinlessFermion
run_case spinlessgc SpinlessFermionGC

echo "Spinless NBodyG serial/MPI outputs match for inter-process factors."
