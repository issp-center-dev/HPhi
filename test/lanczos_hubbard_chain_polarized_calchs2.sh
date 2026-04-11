#!/bin/sh -e

# Regression test for CalcHS=2 on the fully polarized Hubbard sector
# (Nup=1, Ndown=0). Exercises the Ndown==0 else branch of
# sz_hacker_for_large_systems, which is where the boundary check
# tmp_i_up > tmp_i_max was off-by-one for Nup=1 (snoob(2^(Nsite-1))
# returns tmp_i_max exactly and used to fall through to a spurious
# list_1_[idim_max+1] write, a 1-element heap OOB). Allocator slack
# means the energy check alone will still pass if the bug returns,
# so this test is primarily a ground-state regression guard; ASan
# builds should be used to catch the OOB directly.

mkdir -p lanczos_hubbard_chain_polarized_calchs2/
cd lanczos_hubbard_chain_polarized_calchs2

cat > stan.in <<EOF
L = 4
model = "FermionHubbard"
method = "Lanczos"
lattice = "chain"
t = 1.0
U = 4.0
nelec = 1
2Sz = 1
outputmode = "all"
EOF

${MPIRUN} ../../src/HPhi -s stan.in
echo "CalcHS         2" >> modpara.def

rm -rf output
mkdir -p output
${MPIRUN} ../../src/HPhi -e namelist.def

cat > reference.dat <<EOF
   -2.0000000000000004
    0.0000000000000000
    0.5000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste1.dat`

test "${diff}" = "0.000000"
exit $?
