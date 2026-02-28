#!/bin/sh -e
# Test MPI consistency for SpinGC model (includes Exchange and PairLift)
# Compares eigenvalues from MPI run vs non-MPI run

TOLERANCE="0.000001"

mkdir -p mpi_consistency_spingc/
cd mpi_consistency_spingc

cat > stan.in <<EOF
model = "SpinGC"
method = "lanczos"
lattice = "Honeycomb"
W = 2
L = 3
J = 1.0
2S = 1
lanczos_max = 1000
EOF

# Run without MPI
echo "Running without MPI..."
../../src/HPhi -s stan.in
cp output/zvo_energy.dat energy_nompi.dat

# Clean output for MPI run
rm -rf output

# Run with MPI
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

result=$(awk -v diff=${diff} -v tol=${TOLERANCE} 'BEGIN { print (diff < tol) ? "PASS" : "FAIL" }')

if [ "${result}" = "PASS" ]; then
    echo "MPI consistency test PASSED for SpinGC model"
    exit 0
else
    echo "MPI consistency test FAILED for SpinGC model"
    echo "Expected difference < ${TOLERANCE}, got ${diff}"
    cat compare.dat
    exit 1
fi
