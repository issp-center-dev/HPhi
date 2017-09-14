#!/bin/sh -e

mkdir -p spectrum_spingc_honey/
cd spectrum_spingc_honey
#
# Ground state
#
cat > stan1.in <<EOF
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
h=0.1
EigenVecIO = out
SpectrumQW = 0.5
SpectrumQL = 0.5
NOmega = 5
OmegaIm = 1.0
EOF

${MPIRUN} ../../src/HPhi -s stan1.in
#
# Sz-Sz spectrum
#
cp stan1.in stan2.in
cat >> stan2.in <<EOF
CalcSpec = "Normal"
SpectrumType = "SzSz"
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -69.4500000000 1.0000000000 -0.0345775841 -0.0005618651
  -41.6700000000 1.0000000000 -0.0630142521 -0.0018695444
  -13.8900000000 1.0000000000 -0.3594890516 -0.0657728604
  13.8900000000 1.0000000000 0.0982841459 -0.0045694753
  41.6700000000 1.0000000000 0.0430380901 -0.0008709165
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($3-$7)**2)+sqrt(($4-$8)**2)} END{printf "%8.6f", diff}' paste1.dat`
#
# S+S- spectrum
#
cp stan1.in stan2.in
cat >> stan2.in <<EOF
CalcSpec = "Normal"
SpectrumType = "S+S-"
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -69.4500000000 1.0000000000 -0.0781153561 -0.0012671333
  -41.6700000000 1.0000000000 -0.1421604813 -0.0042050255
  -13.8900000000 1.0000000000 -0.8016131992 -0.1446854133
  13.8900000000 1.0000000000 0.2237319517 -0.0104699567
  41.6700000000 1.0000000000 0.0976362171 -0.0019807240
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste2.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)**2)+sqrt(($4-$8)**2)} END{printf "%8.6f", diff}' paste2.dat`
#
# Density-Density spectrum
#
cp stan1.in stan2.in
cat >> stan2.in <<EOF
CalcSpec = "Normal"
SpectrumType = "Density"
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -69.4500000000 1.0000000000 0.0000000000 0.0000000000
  -41.6700000000 1.0000000000 0.0000000000 0.0000000000
  -13.8900000000 1.0000000000 0.0000000000 0.0000000000
  13.8900000000 1.0000000000 0.0000000000 0.0000000000
  41.6700000000 1.0000000000 0.0000000000 0.0000000000
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste3.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)**2)+sqrt(($4-$8)**2)} END{printf "%8.6f", diff}' paste3.dat`

test "${diff}" = "0.000000"

exit $?
