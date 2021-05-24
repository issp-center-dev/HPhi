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
J0y =  2.0
J0z =  3.3
J1x =  1.0
J1y = 1.0
J1z =  0.0
J2x =  3.0
J2y =  2.0
J2z = 1.0
2S=1
h=-0.1
exct = 3
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
 0
  -10.6884732749673663
  0.0000000000000000
  0.0141890312794323

 1
  -9.9987732950471511
  0.0000000000000000
  0.2032909253553725

 2
  -9.4973428762800314
  0.0000000000000000
  0.0040230516306462
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
