#!/bin/sh -e

mkdir -p fulldiag_output_open_failure/
cd fulldiag_output_open_failure

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
printf "Solver  3\nNGPU  0\n" >> calcmod.def

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

echo "fulldiag_output_open_failure: OK"
