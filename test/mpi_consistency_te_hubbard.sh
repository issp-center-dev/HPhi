#!/bin/sh -e
# Test MPI consistency for Hubbard time evolution (AC Laser)
# This tests the fix for stale coefficients in batched MPI processing

TOLERANCE="0.000001"

# Check that MPIRUN is set and non-empty
if [ -z "${MPIRUN}" ]; then
    echo "Error: MPIRUN is not set. Please set MPIRUN to run MPI tests."
    echo "Example: MPIRUN=\"mpirun -np 2\" make test"
    exit 1
fi


mkdir -p mpi_consistency_te_hubbard/
cd mpi_consistency_te_hubbard

# First run Lanczos to get ground state
cat > stan.in <<EOF
model = "Hubbard"
method = "lanczos"
lattice = "square"
a0w = 3
a0l = 0
a1w = 0
a1l = 3
t = 1.0
U = 10.0
nelec = 8
2Sz = 0
lanczos_max = 1000
initial_iv = 1
EigenvecIO = "out"
EOF

echo "Running Lanczos without MPI..."
../../src/HPhi -s stan.in

# Time evolution without MPI
cat > stan_te.in <<EOF
model = "Hubbard"
method = "Time-Evolution"
lattice = "square"
a0w = 3
a0l = 0
a1w = 0
a1l = 3
t = 1.0
U = 10.0
nelec = 8
2Sz = 0
lanczos_max = 20
initial_iv = 1
EigenvecIO = "in"
dt = 0.01
tshift = 0.1
freq = 10.0
PumpType = "AC Laser"
VecPotW = 0.5
VecPotL = 0.5
EOF

echo "Running Time Evolution without MPI..."
../../src/HPhi -s stan_te.in
cp output/Flct.dat flct_nompi.dat

# Clean and run with MPI
rm -rf output

echo "Running Lanczos with MPI..."
${MPIRUN} ../../src/HPhi -s stan.in

echo "Running Time Evolution with MPI..."
${MPIRUN} ../../src/HPhi -s stan_te.in
cp output/Flct.dat flct_mpi.dat

# Compare Flct values (columns 2-5 contain physical quantities)
paste flct_nompi.dat flct_mpi.dat > compare.dat
diff=$(awk -v tol=${TOLERANCE} '
BEGIN { maxdiff = 0.0 }
NR > 1 {
    # Compare N (col 2)
    d = sqrt(($2 - $10) * ($2 - $10))
    if (d > maxdiff) maxdiff = d
    # Compare N^2 (col 3)
    d = sqrt(($3 - $11) * ($3 - $11))
    if (d > maxdiff) maxdiff = d
    # Compare D (col 4)
    d = sqrt(($4 - $12) * ($4 - $12))
    if (d > maxdiff) maxdiff = d
    # Compare D^2 (col 5)
    d = sqrt(($5 - $13) * ($5 - $13))
    if (d > maxdiff) maxdiff = d
}
END { printf "%.10f", maxdiff }
' compare.dat)

echo "Max Flct difference: ${diff}"

result=$(awk -v diff=${diff} -v tol=${TOLERANCE} 'BEGIN { print (diff < tol) ? "PASS" : "FAIL" }')

if [ "${result}" = "PASS" ]; then
    echo "MPI consistency test PASSED for Hubbard time evolution"
    exit 0
else
    echo "MPI consistency test FAILED for Hubbard time evolution"
    echo "Expected difference < ${TOLERANCE}, got ${diff}"
    echo "First few lines of comparison:"
    head -5 compare.dat
    exit 1
fi
