#!/bin/sh -e

mkdir -p lobcg_spingc_honey/
cd lobcg_spingc_honey

cat > stan.in <<EOF
W = 2
L = 3
model = "SpinGC"
method = "CG"
lattice = "Honeycomb"
J0x = -1.0
J0y =  0.0
J0z =  0.0
J1x =  0.0
J1y = -1.0
J1z =  0.0
J2x =  0.0
J2y =  0.0
J2z = -1.0
2S=1
exct = 3
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
 0
  -2.4500706750607772 
  0.0000000000000000 
  -0.0000000029100357 

 1
  -2.4500706750607781 
  0.0000000000000000 
  0.0000000006622576 

 2
  -2.4500706750606782 
  0.0000000000000000 
  0.0000003816727221 
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)**2)} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
