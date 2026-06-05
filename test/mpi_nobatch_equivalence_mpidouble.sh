#!/bin/sh -e
# Assert HPHI_MPI_NOBATCH=1 (per-term MPI) reproduces the batched result for the
# MPIdouble transfer path (both transfer sites inter-process). This needs 16 MPI
# ranks so that a 4-site chain puts sites 2,3 inter-process, exercising
# X_child_general_hopp_MPIdouble_batched / X_child_GC_general_hopp_MPIdouble_batched.
# The existing mpi_nobatch_equivalence test only reaches MPIsingle at np=4.
if [ -z "${MPIRUN}" ]; then echo "MPIRUN not set. Skipping."; exit 0; fi

mkdir -p mpi_nobatch_equivalence_mpidouble
cd mpi_nobatch_equivalence_mpidouble

run_one() {  # $1 = tag; stan.in body on stdin
  tag="$1"
  cat > "stan_${tag}.in"
  rm -rf output
  ${MPIRUN} ../../src/HPhi -s "stan_${tag}.in" > "log_${tag}_batched.txt" 2>&1
  cp output/zvo_energy.dat "energy_${tag}_batched.dat"
  # Confirm the batched MPIdouble path actually fired (not "No inter-process transfers").
  if ! grep -q "MPIdouble:.*-> .* groups" "log_${tag}_batched.txt"; then
    echo "[${tag}] ERROR: batched MPIdouble path did not fire"; grep -i "MPIdouble" "log_${tag}_batched.txt" || true; exit 1
  fi
  rm -rf output
  HPHI_MPI_NOBATCH=1 ${MPIRUN} ../../src/HPhi -s "stan_${tag}.in" > "log_${tag}_nobatch.txt" 2>&1
  cp output/zvo_energy.dat "energy_${tag}_nobatch.dat"
  d=$(paste "energy_${tag}_batched.dat" "energy_${tag}_nobatch.dat" \
      | awk '$2 ~ /^[-+0-9.eE]+$/ && $4 ~ /^[-+0-9.eE]+$/ { x=$2-$4; if(x<0)x=-x; if(x>m)m=x } END{ printf "%.3e", m+0 }')
  echo "[${tag}] max |batched - nobatch| energy = ${d}"
  awk -v d="${d}" 'BEGIN{ exit (d < 1e-10) ? 0 : 1 }' || { echo "[${tag}] MISMATCH"; exit 1; }
}

run_one hubbard <<EOF
model = "Hubbard"
method = "Lanczos"
lattice = "chain"
L = 4
t = 1.0
U = 4.0
nelec = 4
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

echo "MPIdouble: batched == no-batch within tolerance."
