#!/bin/sh

mkdir -p lobcg_hubbard_square/
cd lobcg_hubbard_square

cat > stan.in <<EOF
W = 4
L = 2
model = "FermionHubbard"
method = "CG"
lattice = "Tetragonal"
t = 1.0
U = 4.0
nelec = 8
2Sz = 0
exct = 3
EOF

${MPIRUN} ../../src/HPhi -s stan.in
if test $? -ne 0 ; then
    exit 1
fi

echo "\nCheck value\n"

cat > reference.dat <<EOF
   0
 -10.2529529552637033 
   1.0444442739280635 
   0.0000000000000000 

   1
  -9.8328986753266356 
   1.1830984216215232 
   0.0000000000000000 

   2
  -9.2310650024220635 
   1.2843076003322857 
   0.0000000000000000 
EOF

diff=`paste output/zvo_energy.dat reference.dat | awk '
BEGIN{diff=0.0} {diff+=sqrt(($2-$3)**2)} END{printf "%8.6f", diff}'`

test "${diff}" = "0.000000"

exit $?
