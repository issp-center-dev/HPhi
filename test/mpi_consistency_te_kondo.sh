#!/bin/sh -e
# Test MPI consistency for Kondo time evolution
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


mkdir -p mpi_consistency_te_kondo/
cd mpi_consistency_te_kondo

# First run Lanczos to get ground state
cat > stan.in <<EOF
model = "Kondo"
method = "lanczos"
lattice = "Chain"
L = 4
t = 1.0
J = 2.0
nelec = 4
2Sz = 0
lanczos_max = 1000
initial_iv = 1
EigenvecIO = "out"
EOF

echo "Running Lanczos without MPI..."
../../src/HPhi -s stan.in

# Time evolution without MPI
cat > stan_te.in <<EOF
model = "Kondo"
method = "Time-Evolution"
lattice = "Chain"
L = 4
t = 1.0
J = 2.0
nelec = 4
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

# Compare Flct values
paste flct_nompi.dat flct_mpi.dat > compare.dat
diff=$(awk -v tol=${TOLERANCE} '
BEGIN { maxdiff = 0.0 }
NR > 1 {
    # Compare columns 2-5 (physical quantities)
    for (i = 2; i <= 5; i++) {
        d = sqrt(($(i) - $(i+8)) * ($(i) - $(i+8)))
        if (d > maxdiff) maxdiff = d
    }
}
END { printf "%.10f", maxdiff }
' compare.dat)

echo "Max Flct difference: ${diff}"

result=$(awk -v diff=${diff} -v tol=${TOLERANCE} 'BEGIN { print (diff < tol) ? "PASS" : "FAIL" }')

if [ "${result}" = "PASS" ]; then
    echo "MPI consistency test PASSED for Kondo time evolution"
    exit 0
else
    echo "MPI consistency test FAILED for Kondo time evolution"
    echo "Expected difference < ${TOLERANCE}, got ${diff}"
    echo "First few lines of comparison:"
    head -5 compare.dat
    exit 1
fi
