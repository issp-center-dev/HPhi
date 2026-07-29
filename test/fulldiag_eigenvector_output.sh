#!/bin/sh -eu

SRCDIR="$(cd "$1" && pwd)"
TEST_WORK_ROOT="$(pwd)"
HPHI="$(cd "${TEST_WORK_ROOT}/.." && pwd)/src/HPhi"
VERIFY="${SRCDIR}/test/verify_fulldiag_eigenvectors.py"
SOLVER="${HPHI_FULLDIAG_SOLVER:-0}"

make_input() {
  cat > stan.in <<EOF
L = 4
model = "Spin"
method = "FullDiag"
lattice = "chain"
J = 1.0
2Sz = 0
EOF
  "${HPHI}" -sdry stan.in > generate.log 2>&1
  sed 's/^OutputEigenVec.*/OutputEigenVec  1/' calcmod.def > calcmod.tmp
  mv calcmod.tmp calcmod.def
}

mkdir -p fulldiag_eigenvector_output/
cd fulldiag_eigenvector_output

# Solver 0 is the format and eigenspace reference in every build.
mkdir -p lapack/
cd lapack
make_input
printf "Solver  0\nExpecMode  0\n" >> calcmod.def
"${HPHI}" -e namelist.def > run.log 2>&1
python3 "${VERIFY}" output zvo
REFERENCE_OUTPUT="$(pwd)/output"
cd ..

if [ "${SOLVER}" = "0" ]; then
  echo "fulldiag_eigenvector_output: serial LAPACK OK"
  exit 0
fi

case "${SOLVER}" in
  1) solver_name="scalapack" ;;
  3) solver_name="elpa" ;;
  *)
    echo "ERROR: unsupported HPHI_FULLDIAG_SOLVER=${SOLVER}" >&2
    exit 1
    ;;
esac

# The same global-file contract must hold for the gathered path (Mode 0)
# and the state-panel writers (Modes 1 and 2).
for mode in 0 1 2; do
  casedir="${solver_name}_mode${mode}"
  mkdir -p "${casedir}/"
  cd "${casedir}"
  make_input
  printf "Solver  %s\nExpecMode  %s\n" "${SOLVER}" "${mode}" >> calcmod.def
  if [ "${SOLVER}" = "3" ]; then
    printf "NGPU  0\n" >> calcmod.def
  fi
  ${MPIRUN} "${HPHI}" -e namelist.def > run.log 2>&1
  python3 "${VERIFY}" output zvo "${REFERENCE_OUTPUT}"
  cd ..
done

echo "fulldiag_eigenvector_output: ${solver_name} Modes 0/1/2 OK"
