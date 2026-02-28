#!/bin/sh -e
# Test script for HubbardGC MPIdouble batching
# This test uses 16 MPI ranks with a 4-site chain to trigger MPIdouble transfers
# (both sites in inter-process region)

# Check that MPIRUN is set and non-empty
if [ -z "${MPIRUN}" ]; then
    echo "Error: MPIRUN is not set. Please set MPIRUN to run MPI tests."
    echo "This test requires 16 MPI ranks to trigger MPIdouble transfers."
    echo "Example: MPIRUN=\"mpirun -np 16\" ./lanczos_hubbardgc_mpidouble.sh"
    exit 1
fi

# Note: This test requires 16 MPI ranks to properly test MPIdouble transfers.
# If running with fewer ranks, MPIdouble code paths may not be exercised.

mkdir -p lanczos_hubbardgc_mpidouble/
cd lanczos_hubbardgc_mpidouble

cat > stan.in <<EOF
L = 4
model = "HubbardGC"
method = "Lanczos"
lattice = "Chain"
t = 1.0
U = 4.0
outputmode = "none"
EOF

# Run with 16 MPI ranks to ensure MPIdouble transfers
# For 4 sites: Tpow[8] = 256 = 4^4, with 16 ranks we get 4^2 = 16 states per rank
# Sites 0,1 are intra-process (4 bits each), sites 2,3 are inter-process (4 bits each)
${MPIRUN} ../../src/HPhi -s stan.in

# Check energy value (reference: exact diagonalization)
cat > reference.dat <<EOF
    -3.4185507188738526
    0.0000000000000000
    0.0000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste1.dat`

# Tolerance check
test "${diff}" = "0.000000"
exit $?
