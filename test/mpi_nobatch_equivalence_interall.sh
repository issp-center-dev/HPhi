#!/bin/sh -e
# Assert HPHI_MPI_NOBATCH=1 (per-term MPI) reproduces the batched result for the
# HubbardGC off-diagonal InterAll path (X_child_GC_InterAll_Hubbard_MPI_batched),
# the most complex batched routine and the one PR #216 review H-1 found a stale-
# state bug in. No StdFace standard-mode input produces an inter-process
# off-diagonal InterAll, so this builds the base defs with -sdry and injects an
# explicit off-diagonal InterAll term coupling site 0 to the last site (site 3),
# which straddles the inter-process boundary at np=4.
if [ -z "${MPIRUN}" ]; then echo "MPIRUN not set. Skipping."; exit 0; fi

mkdir -p mpi_nobatch_equivalence_interall
cd mpi_nobatch_equivalence_interall

cat > stan.in <<EOF
model = "HubbardGC"
method = "Lanczos"
lattice = "chain"
L = 4
t = 1.0
U = 4.0
Lanczos_max = 2000
initial_iv = 1
EOF

# Generate the standard-mode definition files (serial -sdry; np-independent).
../../src/HPhi -sdry stan.in > log_sdry.txt 2>&1

# Inject an explicit off-diagonal InterAll term + its Hermitian conjugate:
#   c^dag_{3,up} c_{0,up} n_{0,down}   and its h.c.
# Couples site 0 and site 3 -> inter-process at np=4 -> batched InterAll path.
cat > interall.def <<EOF
========================
NInterAll 2
========================
========zInterAll=======
========================
    3 0 0 0 0 1 0 1 0.5000000000000000 0.0000000000000000
    0 1 0 1 0 0 3 0 0.5000000000000000 0.0000000000000000
EOF
# Register the InterAll file in the StdFace-generated namelist.
printf '    InterAll  interall.def\n' >> namelist.def

run_side() {  # $1 = tag, $2 = HPHI_MPI_NOBATCH value (empty or 1)
  tag="$1"; nb="$2"
  rm -rf output
  if [ -n "$nb" ]; then env_pfx="HPHI_MPI_NOBATCH=$nb"; else env_pfx=""; fi
  env $env_pfx ${MPIRUN} ../../src/HPhi -e namelist.def > "log_${tag}.txt" 2>&1
  cp output/zvo_energy.dat "energy_${tag}.dat"
}

run_side batched ""
# Confirm the batched off-diagonal InterAll path actually fired.
if ! grep -q "InterAll:.*MPI terms -> .* groups" log_batched.txt; then
  echo "ERROR: batched HubbardGC InterAll path did not fire"; grep -i "InterAll" log_batched.txt || true; exit 1
fi
run_side nobatch 1

d=$(paste energy_batched.dat energy_nobatch.dat \
    | awk '$2 ~ /^[-+0-9.eE]+$/ && $4 ~ /^[-+0-9.eE]+$/ { x=$2-$4; if(x<0)x=-x; if(x>m)m=x } END{ printf "%.3e", m+0 }')
echo "[interall] max |batched - nobatch| energy = ${d}"
awk -v d="${d}" 'BEGIN{ exit (d < 1e-10) ? 0 : 1 }' || { echo "[interall] MISMATCH"; exit 1; }
echo "HubbardGC off-diagonal InterAll: batched == no-batch within tolerance."
