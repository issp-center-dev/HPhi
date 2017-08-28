#!/bin/sh

mkdir -p lobcg_spin_kagome/
cd lobcg_spin_kagome

cat > stan.in <<EOF
a0w = 1
a0l = 1
a1w = -1
a1l = 2
model = "Spin"
method = "CG"
lattice = "kagome"
J = 1.0
2Sz = 1
exct = 3
EOF

${MPIRUN} ../../src/HPhi -s stan.in
if test $? -ne 0 ; then
    exit 1
fi

echo "\nCheck value\n"

cat > reference.dat <<EOF
 0
  -3.9690017499285175 
  0.0000000000000000 
  0.0000000000000000 

 1
  -3.9690017499285180 
  0.0000000000000000 
  0.0000000000000000 

 2
  -3.9690017499284811 
  0.0000000000000000 
  0.0000000000000000 
EOF

diff=`paste output/zvo_energy.dat reference.dat | awk '
BEGIN{diff=0.0} {diff+=sqrt(($2-$3)**2)} END{printf "%8.6f", diff}'`

echo ${diff}
test "${diff}" = "0.000000"

exit $?
