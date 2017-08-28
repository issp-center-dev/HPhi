#!/bin/sh -e

mkdir -p lobcg_hubbardgc_tri/
cd lobcg_hubbardgc_tri

cat > stan.in <<EOF
a0w = 3
a0l = 0
a1w = -1
a1l = 2
model = "HubbardGC"
method = "CG"
lattice = "Triangular"
t = 1.0
U = 4.0
exct = 3
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
 0
  -12.3140744040600154 
  0.3472086205682052 
  -0.5677676388658205 

 1
  -12.3140744040600367 
  0.3472086205821575 
  0.8181548814798043 

 2
  -12.3140744040597578 
  0.3472086201828553 
  -0.2503872426138326 
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)**2)} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
