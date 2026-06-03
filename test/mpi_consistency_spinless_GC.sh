#!/bin/sh -e
# Test MPI consistency for SpinlessFermionGC model
# Compares eigenvalues and one-body Green's function from MPI run vs non-MPI run

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
python3 "$1/test/testSpinlessCalc.py" -p "../../src/HPhi" -m "SpinlessFermionGC" -s 8 -V 0.5 --onebody-offdiag

# Run without MPI
echo "Running without MPI..."
../../src/HPhi -e namelist.def
cp output/zvo_energy.dat energy_nompi.dat
cp output/zvo_cisajs.dat onebody_nompi.dat

# Clean output for MPI run
rm -rf output

# Run with MPI
echo "Running with MPI..."
${MPIRUN} ../../src/HPhi -e namelist.def
cp output/zvo_energy.dat energy_mpi.dat
cp output/zvo_cisajs.dat onebody_mpi.dat

# Compare eigenvalues
paste energy_nompi.dat energy_mpi.dat > compare_energy.dat
energy_diff=$(awk -v tol=${TOLERANCE} '
BEGIN { maxdiff = 0.0 }
{
    d = sqrt(($2 - $4) * ($2 - $4))
    if (d > maxdiff) maxdiff = d
}
END { printf "%.10f", maxdiff }
' compare_energy.dat)

echo "Max eigenvalue difference: ${energy_diff}"

# Compare one-body Green's function
paste onebody_nompi.dat onebody_mpi.dat > compare_onebody.dat
onebody_diff=$(awk '
BEGIN { maxdiff = 0.0 }
{
    dre = ($5 - $11)
    dim = ($6 - $12)
    d = sqrt(dre * dre + dim * dim)
    if (d > maxdiff) maxdiff = d
}
END { printf "%.10f", maxdiff }
' compare_onebody.dat)

echo "Max one-body Green function difference: ${onebody_diff}"

result=$(awk -v de=${energy_diff} -v dg=${onebody_diff} -v tol=${TOLERANCE} \
    'BEGIN { print (de < tol && dg < tol) ? "PASS" : "FAIL" }')

if [ "${result}" = "PASS" ]; then
    echo "MPI consistency test PASSED for SpinlessFermionGC model"
    exit 0
else
    echo "MPI consistency test FAILED for SpinlessFermionGC model"
    echo "Expected both differences < ${TOLERANCE}"
    echo "  Eigenvalue difference: ${energy_diff}"
    echo "  One-body Green function difference: ${onebody_diff}"
    cat compare_energy.dat
    cat compare_onebody.dat
    exit 1
fi
