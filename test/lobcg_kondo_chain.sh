#!/bin/sh -e

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
exct = 2
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
 0
  -12.6776213781764220 
  0.1072445462409922 
  0.0000000000000000 

 1
  -9.8347989641820277 
  0.4289751043491778 
  0.0000000000000000 
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
