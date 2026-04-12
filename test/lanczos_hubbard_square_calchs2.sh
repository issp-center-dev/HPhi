#!/bin/sh -e

# Regression test for CalcHS=2 (sz_hacker_for_large_systems) on the
# balanced Hubbard model. Exercises the Hubbard branch of the new path
# (Nup>0, Ndown>0) and verifies that the ground-state energy matches
# the CalcHS=0 reference.

mkdir -p lanczos_hubbard_square_calchs2/
cd lanczos_hubbard_square_calchs2

cat > stan.in <<EOF
W = 4
L = 2
model = "FermionHubbard"
method = "Lanczos"
lattice = "Tetragonal"
t = 1.0
U = 4.0
nelec = 8
2Sz = 0
outputmode = "all"
EOF

${MPIRUN} ../../src/HPhi -s stan.in
echo "CalcHS         2" >> modpara.def

rm -rf output
mkdir -p output
${MPIRUN} ../../src/HPhi -e namelist.def

cat > reference.dat <<EOF
  -10.2529529552637637
    1.0444442739280932
    0.0000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste1.dat`

test "${diff}" = "0.000000"
exit $?
