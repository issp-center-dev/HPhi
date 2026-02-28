#!/bin/sh -e
# Test MPI consistency for Hubbard model
# Compares eigenvalues from MPI run vs non-MPI run

TOLERANCE="0.000001"

mkdir -p mpi_consistency_hubbard/
cd mpi_consistency_hubbard

cat > stan.in <<EOF
model = "Hubbard"
method = "lanczos"
lattice = "square"
W = 4
L = 2
t = 1.0
U = 4.0
nelec = 8
2Sz = 0
lanczos_max = 1000
EOF

# Run without MPI
echo "Running without MPI..."
../../src/HPhi -s stan.in
cp output/zvo_energy.dat energy_nompi.dat

# Clean output for MPI run
rm -rf output

# Run with MPI (uses MPIRUN env variable set by ctest)
echo "Running with MPI..."
${MPIRUN} ../../src/HPhi -s stan.in
cp output/zvo_energy.dat energy_mpi.dat

# Compare eigenvalues
paste energy_nompi.dat energy_mpi.dat > compare.dat
diff=$(awk -v tol=${TOLERANCE} '
BEGIN { maxdiff = 0.0 }
{
    d = sqrt(($2 - $4) * ($2 - $4))
    if (d > maxdiff) maxdiff = d
}
END { printf "%.10f", maxdiff }
' compare.dat)

echo "Max eigenvalue difference: ${diff}"

# Check if difference is within tolerance
result=$(awk -v diff=${diff} -v tol=${TOLERANCE} 'BEGIN { print (diff < tol) ? "PASS" : "FAIL" }')

if [ "${result}" = "PASS" ]; then
    echo "MPI consistency test PASSED for Hubbard model"
    exit 0
else
    echo "MPI consistency test FAILED for Hubbard model"
    echo "Expected difference < ${TOLERANCE}, got ${diff}"
    cat compare.dat
    exit 1
fi
