#!/bin/sh -e

SOLVER="${HPHI_FULLDIAG_SOLVER:-3}"
if [ "${SOLVER}" != "1" ] && [ "${SOLVER}" != "3" ]; then
  echo "ERROR: HPHI_FULLDIAG_SOLVER must be 1 or 3" >&2
  exit 1
fi

rootdir="fulldiag_output_open_failure_solver${SOLVER}"
mkdir -p "${rootdir}/"
cd "${rootdir}"

cat > stan.in <<EOF
L = 4
model = "FermionHubbard"
method = "FullDiag"
lattice = "chain"
t = 1.0
U = 4.0
nelec = 4
2Sz = 0
EOF

../../src/HPhi -sdry stan.in
printf "Solver  %s\n" "${SOLVER}" >> calcmod.def
if [ "${SOLVER}" = "3" ]; then
  printf "NGPU  0\n" >> calcmod.def
fi

# Force the rank-0-only Eigenvalue.dat open to fail after diagonalization.
# Every rank must receive the same failure verdict and terminate, rather than
# non-root ranks continuing into collective observable evaluation and hanging.
if [ -d output/Eigenvalue.dat ]; then
  rmdir output/Eigenvalue.dat
else
  rm -f output/Eigenvalue.dat
fi
mkdir output/Eigenvalue.dat
if ${MPIRUN} ../../src/HPhi -e namelist.def > open_failure.log 2>&1; then
  echo "ERROR: an unwritable Eigenvalue.dat path should fail"
  cat open_failure.log
  exit 1
fi
grep -q "FileOpenError" open_failure.log || {
  echo "ERROR: output-open failure diagnostic was not emitted"
  cat open_failure.log
  exit 1
}

# Force the first FullDiag eigenvector file to fail in both distributed
# observable paths. Mode 0 gathers each state to rank 0; Mode 1 redistributes
# states and lets the owning rank write the global file. In either case every
# rank must receive the same I/O verdict and terminate without hanging.
for mode in 0 1; do
  cd ..
  casedir="fulldiag_eigenvector_open_failure_solver${SOLVER}_mode${mode}"
  mkdir -p "${casedir}/"
  cd "${casedir}"
  cp "../${rootdir}/stan.in" .
  ../../src/HPhi -sdry stan.in
  sed 's/^OutputEigenVec.*/OutputEigenVec  1/' calcmod.def > calcmod.tmp
  mv calcmod.tmp calcmod.def
  printf "Solver  %s\nExpecMode  %s\n" "${SOLVER}" "${mode}" >> calcmod.def
  if [ "${SOLVER}" = "3" ]; then
    printf "NGPU  0\n" >> calcmod.def
  fi
  if [ -d output/zvo_eigenvec_0_rank_0.dat ]; then
    rmdir output/zvo_eigenvec_0_rank_0.dat
  else
    rm -f output/zvo_eigenvec_0_rank_0.dat
  fi
  mkdir output/zvo_eigenvec_0_rank_0.dat
  if ${MPIRUN} ../../src/HPhi -e namelist.def > open_failure.log 2>&1; then
    echo "ERROR: an unwritable FullDiag eigenvector path should fail (Mode ${mode})"
    cat open_failure.log
    exit 1
  fi
  grep -q "FileOpenError" open_failure.log || {
    echo "ERROR: eigenvector open failure diagnostic was not emitted (Mode ${mode})"
    cat open_failure.log
    exit 1
  }
done

echo "fulldiag_output_open_failure: OK"
