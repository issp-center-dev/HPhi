#!/bin/sh -e

mkdir -p fulldiag_elpa_hubbard_chain/
cd fulldiag_elpa_hubbard_chain

cat > stan.in <<EOF
L = 4
model = "FermionHubbard"
method = "FullDiag"
lattice = "chain"
t = 1.0
U = 4.0
nelec = 4
2Sz = 0
EOF

../../src/HPhi -sdry stan.in
echo "Solver  3" >> calcmod.def
echo "NGPU    0" >> calcmod.def
# NOTE: this test is registered via add_hphi_mpi_test(... min:2), whose
# run_with_mpi_precheck.sh precheck validates and parses -np/-n out of the
# ${MPIRUN} env var (not ${MPIRUNFC}, which the plain add_hphi_test-registered
# fulldiag_*.sh scripts use and which is normally unset/empty, i.e. serial).
# Launching with ${MPIRUNFC} here would silently run serially even when the
# precheck confirmed MPIRUN has >=2 ranks, defeating the multi-rank gate this
# test exists to exercise. Use the same variable the precheck validated.
${MPIRUN} ../../src/HPhi -e namelist.def

# Full <H> <N> <Sz> <S2> <D> column comparison (all 5 columns, not just
# energy and doublon): Task 6 unified the distributed FullDiag observable
# path (Mode 0 too) so that S2/Sz are computed identically to the serial
# path instead of being zero-filled, so there is no longer a reason to
# exclude those columns here. These are simply the first 7 (lowest-energy)
# rows of the exact same reference used by fulldiag_hubbard_chain.sh
# (non-ELPA, non-distributed) for the identical stan.in -- reused verbatim
# to assert the ELPA/distributed path reproduces the serial physics
# exactly, S2/Sz included.
cat > reference_ed.dat <<EOF
  -2.102748   4.000000  -0.000000   0.000000   0.287325
  -1.806424   4.000000  -0.000000   2.000000   0.335409
  -1.068140   4.000000   0.000000   0.000000   0.277708
  -0.828427   4.000000  -0.000000   2.000000   0.146447
  -0.828427   4.000000  -0.000000   2.000000   0.146447
   0.000000   4.000000   0.000000   6.000000   0.000000
   0.581449   4.000000  -0.000000   0.000000   1.079437
EOF
awk 'NR>1 && NR<=8 {printf "%11.6f %10.6f %10.6f %10.6f %10.6f\n", $1, $2, $3, $4, $5}' output/zvo_phys_Nup2_Ndown2.dat > ed.dat
paste ed.dat reference_ed.dat > paste_ed.dat
diff=`awk 'BEGIN{max=0}{for(i=1;i<=5;i++){d=$i-$(i+5); if(d<0)d=-d; if(d>max)max=d}}END{print max}' paste_ed.dat`
test "`echo "$diff < 0.000001" | bc`" = "1"

# OutputHam is incompatible with distributed generation (Solver 3, nproc>1):
# must be rejected at startup with a clear message.
cd ..
mkdir -p fulldiag_elpa_hamio_reject/
cd fulldiag_elpa_hamio_reject
cp ../fulldiag_elpa_hubbard_chain/stan.in .
../../src/HPhi -sdry stan.in
printf "Solver  3\nNGPU    0\nOutputHam  1\n" >> calcmod.def
if ${MPIRUN} ../../src/HPhi -e namelist.def > reject.log 2>&1; then
  echo "ERROR: Solver 3 + OutputHam + nproc>1 must fail at startup"
  exit 1
fi
grep -Eqi "OutputHam|InputHam" reject.log

echo "fulldiag_elpa_hubbard_chain: OK"
