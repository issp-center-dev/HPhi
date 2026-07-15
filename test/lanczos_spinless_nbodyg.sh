#!/bin/sh
set -e

testname="lanczos_spinless_nbodyg"
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

compare_one() {
  so="$1"
  si="$2"
  label="$3"
  awk -v t="${tol}" -v so="${so}" -v si="${si}" '
    NR == FNR {
      if ($1 == so && $2 == 0 && $3 == si && $4 == 0) {
        ref_re = $5;
        ref_im = $6;
        found_ref = 1;
      }
      next;
    }
    $1 == 1 && $2 == so && $3 == 0 && $4 == si && $5 == 0 {
      dr = $6 - ref_re;
      di = $7 - ref_im;
      if (dr < 0) dr = -dr;
      if (di < 0) di = -di;
      found_nbody = 1;
      ok = (found_ref && dr < t && di < t);
    }
    END { exit (found_ref && found_nbody && ok) ? 0 : 1; }
  ' output/zvo_cisajs.dat output/zvo_NBodyG.dat || {
    echo "${label} NBodyG one-body value does not match cisajs"
    cat output/zvo_cisajs.dat
    cat output/zvo_NBodyG.dat
    exit 1
  }
}

compare_two() {
  label="$1"
  awk -v t="${tol}" '
    NR == FNR {
      if ($1 == 0 && $2 == 0 && $3 == 1 && $4 == 0 &&
          $5 == 2 && $6 == 0 && $7 == 2 && $8 == 0) {
        ref_re = $9;
        ref_im = $10;
        found_ref = 1;
      }
      next;
    }
    $1 == 2 && $2 == 0 && $3 == 0 && $4 == 1 && $5 == 0 &&
        $6 == 2 && $7 == 0 && $8 == 2 && $9 == 0 {
      dr = $10 - ref_re;
      di = $11 - ref_im;
      if (dr < 0) dr = -dr;
      if (di < 0) di = -di;
      found_nbody = 1;
      ok = (found_ref && dr < t && di < t);
    }
    END { exit (found_ref && found_nbody && ok) ? 0 : 1; }
  ' output/zvo_cisajscktalt.dat output/zvo_NBodyG.dat || {
    echo "${label} NBodyG two-body value does not match cisajscktalt"
    cat output/zvo_cisajscktalt.dat
    cat output/zvo_NBodyG.dat
    exit 1
  }
}

run_case() {
  tag="$1"
  model="$2"

  rm -rf "${tag}"
  mkdir -p "${tag}"
  cd "${tag}"
  if [ "${model}" = "SpinlessFermion" ]; then
    python3 "${SRCDIR}/test/testSpinlessCalc.py" -p /bin/true \
      -m "${model}" -s 6 -n 3 -V 0.5 --onebody-offdiag --offdiag > log_generate.txt 2>&1
  else
    python3 "${SRCDIR}/test/testSpinlessCalc.py" -p /bin/true \
      -m "${model}" -s 6 -V 0.5 --onebody-offdiag --offdiag > log_generate.txt 2>&1
  fi
  printf '    NBodyG  nbodyg.def\n' >> namelist.def
  cat > nbodyg.def <<EOF
========================
NNBodyG 3
========================
========NBodyG==========
========================
1 0 0 0 0
1 0 0 1 0
2 0 0 1 0 2 0 2 0
EOF
  run_hphi log_lanczos.txt "${hphi}" -e namelist.def
  test -f output/zvo_NBodyG.dat || { echo "zvo_NBodyG.dat was not generated"; exit 1; }

  compare_one 0 0 "${model}"
  compare_one 0 1 "${model}"
  compare_two "${model}"
  cd ..
}

run_case spinless SpinlessFermion
run_case spinlessgc SpinlessFermionGC

echo "Spinless NBodyG N=1 and N=2 values match legacy Green's functions."
