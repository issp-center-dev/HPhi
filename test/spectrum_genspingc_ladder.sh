#!/bin/sh -e

mkdir -p spectrum_genspingc_ladder/
cd spectrum_genspingc_ladder
#
# Ground state
#
cat > stan1.in <<EOF
L = 2
W = 2
model = "SpinGC"
method = "CG"
lattice = "ladder"
J0 = 1.0
J1 = 1.0
2S = 3
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
-275.1384387633 1.0000000000 -0.0363484913 -0.0001403438
-165.0830632580 1.0000000000 -0.0632064339 -0.0004244105
-55.0276877527 1.0000000000 -0.2421941384 -0.0062432118
55.0276877527 1.0000000000 0.1323558007 -0.0018620664
165.0830632580 1.0000000000 0.0519570272 -0.0002867699
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste1.dat`
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
-275.1384387633 1.0000000000 -0.0726969826 -0.0002806877
-165.0830632580 1.0000000000 -0.1264128678 -0.0008488211
-55.0276877527 1.0000000000 -0.4843882767 -0.0124864235
55.0276877527 1.0000000000 0.2647116014 -0.0037241328
165.0830632580 1.0000000000 0.1039140545 -0.0005735397
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste2.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste2.dat`

test "${diff}" = "0.000000"

exit $?
