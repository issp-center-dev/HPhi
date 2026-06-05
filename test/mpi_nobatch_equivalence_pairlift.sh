#!/bin/sh -e
# Assert HPHI_MPI_NOBATCH=1 (per-term MPI) reproduces the batched result for the
# SpinGC PairLift MPIsingle path. The SpinGC batched init scans Exchange AND
# PairLift together, but every existing SpinGC test input is pure Heisenberg
# (zero PairLift terms), so the PairLift branch is never exercised. This injects
# an explicit PairLift coupling site 0 to an inter-process site (site 7).
# J is set to 0 so there are NO Exchange terms: the single inter-process
# MPIsingle term the batched group reports then comes ONLY from the PairLift,
# so the path-fired guard below proves the PairLift branch (not Exchange) fired.
if [ -z "${MPIRUN}" ]; then echo "MPIRUN not set. Skipping."; exit 0; fi

mkdir -p mpi_nobatch_equivalence_pairlift
cd mpi_nobatch_equivalence_pairlift

cat > stan.in <<EOF
model = "SpinGC"
method = "Lanczos"
lattice = "chain"
L = 8
2S = 1
J = 0.0
Lanczos_max = 2000
initial_iv = 1
EOF

# Generate the standard-mode definition files (serial -sdry).
../../src/HPhi -sdry stan.in > log_sdry.txt 2>&1

# Inject a PairLift term coupling site 0 and the inter-process site 7.
cat > pairlift.def <<EOF
=============================================
NPairLift          1
=============================================
====== Pair-Lift term ============
=============================================
    0     7         0.300000000000000
EOF
printf '    PairLift  pairlift.def\n' >> namelist.def

run_side() {  # $1 = tag, $2 = HPHI_MPI_NOBATCH value (empty or 1)
  tag="$1"; nb="$2"
  rm -rf output
  if [ -n "$nb" ]; then env_pfx="HPHI_MPI_NOBATCH=$nb"; else env_pfx=""; fi
  env $env_pfx ${MPIRUN} ../../src/HPhi -e namelist.def > "log_${tag}.txt" 2>&1
  cp output/zvo_energy.dat "energy_${tag}.dat"
}

run_side batched ""
# With J=0 the only inter-process MPIsingle term is the injected PairLift, so the
# batched group must report exactly "1 MPIsingle terms" -> this proves the
# PairLift branch fired (without the PairLift it would be "No MPIsingle terms").
if ! grep -q "SpinGC Exchange: 1 MPIsingle terms -> 1 groups" log_batched.txt; then
  echo "ERROR: batched SpinGC PairLift path did not fire as expected"; grep -i "SpinGC Exchange" log_batched.txt || true; exit 1
fi
run_side nobatch 1

d=$(paste energy_batched.dat energy_nobatch.dat \
    | awk '$2 ~ /^[-+0-9.eE]+$/ && $4 ~ /^[-+0-9.eE]+$/ { x=$2-$4; if(x<0)x=-x; if(x>m)m=x } END{ printf "%.3e", m+0 }')
echo "[pairlift] max |batched - nobatch| energy = ${d}"
awk -v d="${d}" 'BEGIN{ exit (d < 1e-10) ? 0 : 1 }' || { echo "[pairlift] MISMATCH"; exit 1; }
echo "SpinGC PairLift: batched == no-batch within tolerance."
