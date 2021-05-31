#!/bin/sh -e

mkdir -p spectrum_kondogc_chain/
cd spectrum_kondogc_chain
#
# Ground state
#
cat > stan1.in <<EOF
L = 4
model = "KondoGC"
method = "CG"
lattice = "chain"
t = 1.0
J = 4.0
mu = 0.6
h = 0.2
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
  -57.6000000000 1.0000000000 -0.0025722820 -0.0000537686
  -34.5600000000 1.0000000000 -0.0049663255 -0.0002019943
  -11.5200000000 1.0000000000 0.0019939275 -0.0395683493
  11.5200000000 1.0000000000 0.0059245271 -0.0002887485
  34.5600000000 1.0000000000 0.0028058341 -0.0000640101
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
-57.6000000000 1.0000000000 -0.0051231650 -0.0001066420
-34.5600000000 1.0000000000 -0.0098527794 -0.0003974502
-11.5200000000 1.0000000000 -0.0004597280 -0.0893610870
11.5200000000 1.0000000000 0.0119651223 -0.0005890413
34.5600000000 1.0000000000 0.0056373628 -0.0001291993
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
  -57.6000000000 1.0000000000 -0.0084352044 -0.0001719484
  -34.5600000000 1.0000000000 -0.0158985384 -0.0006127767
  -11.5200000000 1.0000000000 -0.1348448908 -0.0724351786
  11.5200000000 1.0000000000 0.0207744673 -0.0010490382
  34.5600000000 1.0000000000 0.0096340777 -0.0002243707
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste3.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste3.dat`
#
# Up-Up spectrum
#
cp stan1.in stan2.in
cat >> stan2.in <<EOF
CalcSpec = "Normal"
SpectrumType = "Down"
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -57.6000000000 1.0000000000 -0.0290903903 -0.0006438572
  -34.5600000000 1.0000000000 -0.0592720194 -0.0026773975
  -11.5200000000 1.0000000000 0.6510110721 -0.7114741311
  11.5200000000 1.0000000000 0.0548650640 -0.0022938688
  34.5600000000 1.0000000000 0.0279845511 -0.0005958293
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste4.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste4.dat`
#
# Down-Down spectrum
#
cp stan1.in stan2.in
cat >> stan2.in <<EOF
CalcSpec = "Normal"
SpectrumType = "Up"
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -57.6000000000 1.0000000000 -0.0289623095 -0.0006381971
  -34.5600000000 1.0000000000 -0.0587434648 -0.0026297578
  -11.5200000000 1.0000000000 0.6189027112 -0.8649587694
  11.5200000000 1.0000000000 0.0553260674 -0.0023326676
  34.5600000000 1.0000000000 0.0281041159 -0.0006009343
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste5.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste5.dat`

test "${diff}" = "0.000000"

exit $?
