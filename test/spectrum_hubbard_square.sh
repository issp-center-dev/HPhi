#!/bin/sh -e

mkdir -p spectrum_hubbard_square/
cd spectrum_hubbard_square
#
# Ground state
#
cat > stan1.in <<EOF
W = 4
L = 2
model = "FermionHubbard"
method = "CG"
lattice = "Tetragonal"
t = 1.0
U = 4.0
nelec = 8
2Sz = 0
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
  -96.0000000000 1.0000000000 -0.0440951047 -0.0005039641
  -57.6000000000 1.0000000000 -0.0785907381 -0.0016034728
  -19.2000000000 1.0000000000 -0.3668649422 -0.0362074986
  19.2000000000 1.0000000000 0.1409881461 -0.0052067704
  57.6000000000 1.0000000000 0.0585833763 -0.0008903186
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
-96.0000000000 1.0000000000 -0.0881902138 -0.0010079283
-57.6000000000 1.0000000000 -0.1571814841 -0.0032069457
-19.2000000000 1.0000000000 -0.7337299221 -0.0724150010
19.2000000000 1.0000000000 0.2819763059 -0.0104135413
57.6000000000 1.0000000000 0.1171667584 -0.0017806372
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste2.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste2.dat`
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
  -96.0000000000 1.0000000000 -0.0413706460 -0.0004441518
  -57.6000000000 1.0000000000 -0.0703904433 -0.0012872601
  -19.2000000000 1.0000000000 -0.2375333432 -0.0149570492
  19.2000000000 1.0000000000 0.1765475276 -0.0081635110
  57.6000000000 1.0000000000 0.0639222321 -0.0010611510
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste3.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste3.dat`
#
# Up-Up spectrum
#
cp stan1.in stan2.in
cat >> stan2.in <<EOF
CalcSpec = "Normal"
SpectrumType = "Up"
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -96.0000000000 1.0000000000 -0.0030684753 -0.0000341177
  -57.6000000000 1.0000000000 -0.0053542539 -0.0001039682
  -19.2000000000 1.0000000000 -0.0210764244 -0.0016381762
  19.2000000000 1.0000000000 0.0109868431 -0.0004402674
  57.6000000000 1.0000000000 0.0043424898 -0.0000683658
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste4.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste4.dat`
#
# Down-Down spectrum
#
cp stan1.in stan2.in
cat >> stan2.in <<EOF
CalcSpec = "Normal"
SpectrumType = "Down"
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -96.0000000000 1.0000000000 -0.0030684753 -0.0000341177
  -57.6000000000 1.0000000000 -0.0053542539 -0.0001039682
  -19.2000000000 1.0000000000 -0.0210764244 -0.0016381762
  19.2000000000 1.0000000000 0.0109868431 -0.0004402674
  57.6000000000 1.0000000000 0.0043424898 -0.0000683658
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste5.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste5.dat`

test "${diff}" = "0.000000"

exit $?
