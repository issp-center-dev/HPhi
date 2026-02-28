#!/bin/sh -e

mkdir -p lanczos_spinless/
cd lanczos_spinless
python3 "$1/test/testSpinlessCalc.py" -p "../../src/HPhi" -mpi "${MPIRUN}" -m "SpinlessFermion" -s 8

# Check energy value only
# Note: Green's functions are not compared because the ground state may be degenerate,
# and MPI/non-MPI runs may find different linear combinations of degenerate eigenvectors.
cat > reference.dat <<EOF
   -4.8284271247461898
    0.0000000000000000
    0.0000000000000000
EOF

paste output/zvo_energy.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste1.dat`

test "${diff}" = "0.000000"
exit $?
