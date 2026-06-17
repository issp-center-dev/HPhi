#!/bin/sh -e
# Assert HPHI_MPI_NOBATCH=1 (per-term MPI) reproduces the batched result for
# HubbardGC off-diagonal InterAll terms whose inter-process factor is a NUMBER
# operator, i.e. terms whose communication partner is the local rank
# (origin==myrank).
#
# This is the complement of mpi_nobatch_equivalence_interall.sh, which only
# injects c^dag_{3,up} c_{0,up} n_{0,down} -- there the inter-process operator
# (c^dag on site 3) flips the rank bit, so origin!=myrank and the term is
# correctly batched. The batched initializer DROPS terms with origin==myrank
# (mltplyMPIBatched.c), and mltplyHubbardGC used to `continue` past every
# inter-PE term under batching, so a term like c^dag_{0,up} c_{1,up} n_{3,up}
# (local hop modulated by an inter-process density) was applied by NEITHER
# path -- silently corrupting H with no warning. See
# docs/2026-06-15-v3.6-prerelease-review.md finding H-1.
#
# The terms below all have site 3 (inter-process at np=4) appearing only as a
# number operator, so every one maps to origin==myrank. With the bug, batched
# mode drops them all and reports the no-InterAll energy; the per-term path
# keeps them. The fix lets origin==myrank terms fall through to the per-term
# path, restoring batched == NOBATCH.
if [ -z "${MPIRUN}" ]; then echo "MPIRUN not set. Skipping."; exit 0; fi

mkdir -p mpi_nobatch_equivalence_interall_selfterm
cd mpi_nobatch_equivalence_interall_selfterm

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

# Off-diagonal InterAll terms whose inter-process factor (site 3) is a number
# operator n_3. Each pair is the term and its Hermitian conjugate (operator
# order reversed). Encoding: i1 s1 i2 s2 i3 s3 i4 s4 = c^dag_{i1,s1} c_{i2,s2}
# c^dag_{i3,s3} c_{i4,s4}.
# The first four pairs put the number operator LAST (c^dag_i c_j n_k ->
# child_GC_CisAjtCkuAku_Hubbard_MPI); the last pair puts it FIRST
# (n_k c^dag_i c_j -> child_GC_CisAisCjtAku_Hubbard_MPI, which delegates to the
# former), so both off-diagonal dispatch targets that route origin==myrank
# terms are exercised.
cat > interall.def <<EOF
========================
NInterAll 10
========================
========zInterAll=======
========================
    0 0 1 0 3 0 3 0 1.0000000000000000 0.0000000000000000
    3 0 3 0 1 0 0 0 1.0000000000000000 0.0000000000000000
    1 0 2 0 3 0 3 0 1.0000000000000000 0.0000000000000000
    3 0 3 0 2 0 1 0 1.0000000000000000 0.0000000000000000
    0 1 1 1 3 0 3 0 1.0000000000000000 0.0000000000000000
    3 0 3 0 1 1 0 1 1.0000000000000000 0.0000000000000000
    0 0 2 0 3 1 3 1 1.0000000000000000 0.0000000000000000
    3 1 3 1 2 0 0 0 1.0000000000000000 0.0000000000000000
    3 0 3 0 0 0 2 0 1.0000000000000000 0.0000000000000000
    2 0 0 0 3 0 3 0 1.0000000000000000 0.0000000000000000
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
run_side nobatch 1

d=$(paste energy_batched.dat energy_nobatch.dat \
    | awk '$2 ~ /^[-+0-9.eE]+$/ && $4 ~ /^[-+0-9.eE]+$/ { x=$2-$4; if(x<0)x=-x; if(x>m)m=x } END{ printf "%.3e", m+0 }')
echo "[interall-selfterm] max |batched - nobatch| energy = ${d}"
awk -v d="${d}" 'BEGIN{ exit (d < 1e-10) ? 0 : 1 }' || { echo "[interall-selfterm] MISMATCH (origin==myrank InterAll terms dropped under batching)"; exit 1; }
echo "HubbardGC origin==myrank InterAll (inter-PE number operator): batched == no-batch within tolerance."
