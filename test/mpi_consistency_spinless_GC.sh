#!/bin/sh -e
# Test MPI consistency for SpinlessFermionGC model
# Compares eigenvalues from MPI run vs non-MPI run

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


mkdir -p mpi_consistency_spinless_GC/
cd mpi_consistency_spinless_GC

# Generate input files using Python script
python3 "$1/test/testSpinlessCalc.py" -p "../../src/HPhi" -m "SpinlessFermionGC" -s 8

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
    echo "MPI consistency test PASSED for SpinlessFermionGC model"
    exit 0
else
    echo "MPI consistency test FAILED for SpinlessFermionGC model"
    echo "Expected difference < ${TOLERANCE}, got ${diff}"
    cat compare.dat
    exit 1
fi
