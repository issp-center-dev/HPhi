#!/bin/sh -e
# Test MPI consistency for SpinlessFermion model
# Compares eigenvalues from MPI run vs non-MPI run

TOLERANCE="0.000001"

mkdir -p mpi_consistency_spinless/
cd mpi_consistency_spinless

# Generate input files using Python script
python3 "$1/test/testSpinlessCalc.py" -p "../../src/HPhi" -m "SpinlessFermion" -s 8

# Run without MPI
echo "Running without MPI..."
../../src/HPhi -e namelist.def
cp output/zvo_energy.dat energy_nompi.dat

# Clean output for MPI run
rm -rf output

# Run with MPI
echo "Running with MPI..."
${MPIRUN} ../../src/HPhi -e namelist.def
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
    echo "MPI consistency test PASSED for SpinlessFermion model"
    exit 0
else
    echo "MPI consistency test FAILED for SpinlessFermion model"
    echo "Expected difference < ${TOLERANCE}, got ${diff}"
    cat compare.dat
    exit 1
fi
