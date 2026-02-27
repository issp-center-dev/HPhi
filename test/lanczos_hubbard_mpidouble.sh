#!/bin/sh -e
# Test script for Hubbard (canonical) MPIdouble batching
# This test uses 16 MPI ranks with a 4-site chain to trigger MPIdouble transfers
# (both sites in inter-process region)

mkdir -p lanczos_hubbard_mpidouble/
cd lanczos_hubbard_mpidouble

cat > stan.in <<EOF
L = 4
model = "Hubbard"
method = "Lanczos"
lattice = "Chain"
t = 1.0
U = 4.0
nelec = 4
2Sz = 0
outputmode = "none"
EOF

# Run with 16 MPI ranks to ensure MPIdouble transfers
# For 4 sites: Tpow[8] = 256 = 4^4, with 16 ranks we get variable states per rank
# Sites 0,1 are intra-process, sites 2,3 are inter-process (triggering MPIdouble)
${MPIRUN} ../../src/HPhi -s stan.in

# Check energy value (reference: exact diagonalization for half-filled 4-site Hubbard)
cat > reference.dat <<EOF
    -2.1027484834620682
    0.0000000000000000
    0.0000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste1.dat`

# Tolerance check
test "${diff}" = "0.000000"
exit $?
