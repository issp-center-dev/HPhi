#!/bin/sh -e
# SpinlessFermion TimeEvolution does not implement diagonal TE one-body/two-body
# handlers. A diagonal TEOneBody term must be rejected with a clear message,
# while an off-diagonal-only TEOneBody run must still be accepted.
# $1 = CMAKE_SOURCE_DIR (for test/testSpinlessCalc.py).

SRCDIR="$1"
HPHI=../../../src/HPhi
TE_STEPS=10
EXPECTED_MSG="time evolution with diagonal one-body / two-body terms is not supported for SpinlessFermion / SpinlessFermionGC"

mkdir -p te_spinless_diagonal_reject
cd te_spinless_diagonal_reject

prepare_spinless_case() {
  case_dir="$1"
  model="$2"

  rm -rf "${case_dir}"
  mkdir -p "${case_dir}"
  (
    cd "${case_dir}"
    python3 "${SRCDIR}/test/testSpinlessCalc.py" -p /bin/true -m "${model}" -s 8 > gen.log 2>&1
    printf '     SpectrumVec  zvo_eigenvec_0\n' >> namelist.def

    sed -e 's/^OutputEigenVec.*/OutputEigenVec   1/' calcmod.def > _t && mv _t calcmod.def
    "${HPHI}" -e namelist.def > gs.log 2>&1

    sed -e 's/^CalcType.*/CalcType   4/' \
        -e 's/^InputEigenVec.*/InputEigenVec   1/' \
        -e 's/^OutputEigenVec.*/OutputEigenVec   0/' \
        calcmod.def > _t && mv _t calcmod.def
    sed -e "s/^Lanczos_max.*/Lanczos_max    ${TE_STEPS}/" modpara.def > _t && mv _t modpara.def
    printf 'ExpandCoef     10\n' >> modpara.def
  )
}

write_diagonal_teonebody() {
  case_dir="$1"
  (
    cd "${case_dir}"
    cat > teonebody.def <<EOF
========================
NTimeSteps    ${TE_STEPS}
========================
=========  OneBody Time Evolution  ==========
========================
EOF
    i=0
    while [ "${i}" -lt "${TE_STEPS}" ]; do
      printf '0.%02d  1\n' "${i}" >> teonebody.def
      printf '%s\n' '0  0  0  0  0.2  0.0' >> teonebody.def
      i=$((i + 1))
    done
    printf '       TEOneBody  teonebody.def\n' >> namelist.def
  )
}

write_offdiag_teonebody() {
  case_dir="$1"
  (
    cd "${case_dir}"
    cat > teonebody.def <<EOF
========================
NTimeSteps    ${TE_STEPS}
========================
=========  OneBody Time Evolution  ==========
========================
EOF
    i=0
    while [ "${i}" -lt "${TE_STEPS}" ]; do
      printf '0.%02d  2\n' "${i}" >> teonebody.def
      printf '%s\n' '0  0  1  0  0.2  0.0' >> teonebody.def
      printf '%s\n' '1  0  0  0  0.2  0.0' >> teonebody.def
      i=$((i + 1))
    done
    printf '       TEOneBody  teonebody.def\n' >> namelist.def
  )
}

run_expect_reject() {
  name="$1"
  model="$2"

  prepare_spinless_case "${name}" "${model}"
  write_diagonal_teonebody "${name}"
  (
    cd "${name}"
    set +e
    "${HPHI}" -e namelist.def > te.log 2>&1
    rc=$?
    set -e

    if [ "${rc}" = "0" ]; then
      echo "[${name}] ERROR: diagonal Spinless TE run succeeded but was expected to fail" >&2
      exit 1
    fi
    if ! grep -q "${EXPECTED_MSG}" te.log; then
      echo "[${name}] ERROR: expected Spinless TE diagonal reject message not found" >&2
      tail -40 te.log >&2
      exit 1
    fi
  )
}

run_expect_accept() {
  name="$1"
  model="$2"

  prepare_spinless_case "${name}" "${model}"
  write_offdiag_teonebody "${name}"
  (
    cd "${name}"
    if ! "${HPHI}" -e namelist.def > te.log 2>&1; then
      echo "[${name}] ERROR: off-diagonal-only Spinless TE run failed" >&2
      tail -40 te.log >&2
      exit 1
    fi
  )
}

run_expect_reject spinless_diagonal SpinlessFermion
run_expect_reject spinless_gc_diagonal SpinlessFermionGC
run_expect_accept spinless_offdiag SpinlessFermion

echo "SpinlessFermion(GC) diagonal TEOneBody reject and off-diagonal TEOneBody accept: OK"
