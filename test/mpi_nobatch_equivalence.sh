#!/bin/sh -e
# Assert HPHI_MPI_NOBATCH=1 (per-term MPI) reproduces the batched result exactly,
# across the model families whose MPI communication is batched.
# Same input, same MPIRUN, batched vs no-batch -> identical energy.
# (Spinless models are Expert-mode only and are covered by the
#  mpi_consistency_spinless / mpi_consistency_spinless_GC tests, which pass in
#  both modes.)
if [ -z "${MPIRUN}" ]; then echo "MPIRUN not set. Skipping."; exit 0; fi

mkdir -p mpi_nobatch_equivalence
cd mpi_nobatch_equivalence

run_one() {  # $1 = tag; stan.in body on stdin
  tag="$1"
  cat > "stan_${tag}.in"
  rm -rf output
  ${MPIRUN} ../../src/HPhi -s "stan_${tag}.in" > "log_${tag}_batched.txt" 2>&1
  cp output/zvo_energy.dat "energy_${tag}_batched.dat"
  rm -rf output
  HPHI_MPI_NOBATCH=1 ${MPIRUN} ../../src/HPhi -s "stan_${tag}.in" > "log_${tag}_nobatch.txt" 2>&1
  cp output/zvo_energy.dat "energy_${tag}_nobatch.dat"
  # zvo_energy.dat rows are "<label>  <value>"; compare the numeric column.
  d=$(paste "energy_${tag}_batched.dat" "energy_${tag}_nobatch.dat" \
      | awk '$2 ~ /^[-+0-9.eE]+$/ && $4 ~ /^[-+0-9.eE]+$/ { x=$2-$4; if(x<0)x=-x; if(x>m)m=x } END{ printf "%.3e", m+0 }')
  echo "[${tag}] max |batched - nobatch| = ${d}"
  awk -v d="${d}" 'BEGIN{ exit (d < 1e-10) ? 0 : 1 }' || { echo "[${tag}] MISMATCH"; exit 1; }
}

run_one hubbard <<EOF
model = "Hubbard"
method = "Lanczos"
lattice = "chain"
L = 8
t = 1.0
U = 4.0
nelec = 8
2Sz = 0
Lanczos_max = 2000
initial_iv = 1
EOF

run_one hubbardgc <<EOF
model = "HubbardGC"
method = "Lanczos"
lattice = "chain"
L = 4
t = 1.0
U = 4.0
Lanczos_max = 2000
initial_iv = 1
EOF

run_one spin <<EOF
model = "Spin"
method = "Lanczos"
lattice = "chain"
L = 8
J = 1.0
2Sz = 0
Lanczos_max = 2000
initial_iv = 1
EOF

run_one spingc <<EOF
model = "SpinGC"
method = "Lanczos"
lattice = "chain"
L = 8
J = 1.0
Lanczos_max = 2000
initial_iv = 1
EOF

echo "All models: batched == no-batch within tolerance."
