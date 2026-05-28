#!/bin/sh -e
# Test MPI consistency for Spin time evolution
# This tests that MPI and non-MPI runs produce identical results

TOLERANCE="0.000001"

# Check that MPIRUN is set and non-empty
if [ -z "${MPIRUN}" ]; then
    echo "Error: MPIRUN is not set. Please set MPIRUN to run MPI tests."
    echo "Example: MPIRUN=\"mpirun -np 2\" make test"
    exit 1
fi

# Require multiple MPI ranks (-np/ -n > 1)
MPI_NP=$(printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}')
if ! printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$"; then
    echo "Error: MPIRUN must include -np or -n with an integer > 1."
    echo "Current MPIRUN: ${MPIRUN}"
    exit 1
fi
if [ "${MPI_NP}" -le 1 ]; then
    echo "Error: MPI consistency tests require more than one rank."
    echo "Current MPIRUN: ${MPIRUN}"
    exit 1
fi


mkdir -p mpi_consistency_te_spin/
cd mpi_consistency_te_spin

# First run Lanczos to get ground state
cat > stan.in <<EOF
model = "Spin"
method = "lanczos"
lattice = "Chain"
L = 8
J = 1.0
2Sz = 0
lanczos_max = 1000
initial_iv = 1
EigenvecIO = "out"
EOF

echo "Running Lanczos without MPI..."
../../src/HPhi -s stan.in

# Time evolution without MPI
cat > stan_te.in <<EOF
model = "Spin"
method = "Time-Evolution"
lattice = "Chain"
L = 8
J = 1.0
2Sz = 0
lanczos_max = 20
initial_iv = 1
EigenvecIO = "in"
dt = 0.01
tshift = 0.1
freq = 5.0
PumpType = "Quench"
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

# Compare Flct values (columns 2-3: Sz and Sz^2)
paste flct_nompi.dat flct_mpi.dat > compare.dat
diff=$(awk -v tol=${TOLERANCE} '
BEGIN { maxdiff = 0.0 }
NR > 1 {
    # Compare Sz (col 2)
    d = sqrt(($2 - $10) * ($2 - $10))
    if (d > maxdiff) maxdiff = d
    # Compare Sz^2 (col 3)
    d = sqrt(($3 - $11) * ($3 - $11))
    if (d > maxdiff) maxdiff = d
}
END { printf "%.10f", maxdiff }
' compare.dat)

echo "Max Flct difference: ${diff}"

result=$(awk -v diff=${diff} -v tol=${TOLERANCE} 'BEGIN { print (diff < tol) ? "PASS" : "FAIL" }')

if [ "${result}" = "PASS" ]; then
    echo "MPI consistency test PASSED for Spin time evolution"
    exit 0
else
    echo "MPI consistency test FAILED for Spin time evolution"
    echo "Expected difference < ${TOLERANCE}, got ${diff}"
    echo "First few lines of comparison:"
    head -5 compare.dat
    exit 1
fi
