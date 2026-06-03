#!/bin/sh -e
# Test script for Hubbard (canonical) MPIdouble batching
# This test uses 16 MPI ranks with a 4-site chain to trigger MPIdouble transfers
# (both sites in inter-process region)

# Check that MPIRUN is set and non-empty
if [ -z "${MPIRUN}" ]; then
    echo "Error: MPIRUN is not set. Please set MPIRUN to run MPI tests."
    echo "This test requires 16 MPI ranks to trigger MPIdouble transfers."
    echo "Example: MPIRUN=\"mpirun -np 16\" ./lanczos_hubbard_mpidouble.sh"
    exit 1
fi

# This test specifically targets MPIdouble paths that require 16 ranks.
MPI_NP=$(printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}')
if ! printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$"; then
    echo "Warning: Could not parse -np/-n from MPIRUN='${MPIRUN}'. Skipping MPIdouble test."
    exit 0
fi
if [ "${MPI_NP}" -ne 16 ]; then
    echo "MPIdouble test requires 16 MPI ranks (current: ${MPI_NP}). Skipping."
    exit 0
fi

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

# Check ground-state energy only (Doublon/Sz are model observables, not zero)
ref_energy="-2.1027484834620682"
energy=$(awk 'NR==1 {print $2}' output/zvo_energy.dat)
diff=$(awk -v a="${energy}" -v b="${ref_energy}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%.12e", d}')

# Tolerance check
awk -v d="${diff}" 'BEGIN{exit (d < 1e-10) ? 0 : 1}'
exit $?
