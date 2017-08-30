#!/bin/sh -e

mkdir -p lanczos_spingc_honey/
cd lanczos_spingc_honey

cat > stan.in <<EOF
W = 2
L = 3
model = "SpinGC"
method = "Lanczos"
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
h=1.0
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
  -7.7479167794065980
   0.0000000000000000
   5.8364467634304367
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)**2)} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
