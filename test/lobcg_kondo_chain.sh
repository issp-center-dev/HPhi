#!/bin/sh

mkdir -p lobcg_kondo_chain/
cd lobcg_kondo_chain

cat > stan.in <<EOF
L = 4
model = "Kondo"
method = "CG"
lattice = "chain"
t = 1.0
J = 4.0
nelec = 4
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
  -12.6776213781764220 
  0.1072445462409922 
  0.0000000000000000 

 1
  -9.8347989641820277 
  0.4289751043491778 
  0.0000000000000000 

 2
  -9.1829527358164000 
  0.3092994621664528 
  0.0000000000000000 
EOF

diff=`paste output/zvo_energy.dat reference.dat | awk '
BEGIN{diff=0.0} {diff+=sqrt(($2-$3)**2)} END{printf "%8.6f", diff}'`

echo ${diff}
test "${diff}" = "0.000000"

exit $?
