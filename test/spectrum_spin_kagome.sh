#!/bin/sh -e

mkdir -p spectrum_spin_kagome/
cd spectrum_spin_kagome
#
# Ground state
#
cat > stan1.in <<EOF
a0w = 1
a0l = 1
a1w = -1
a1l = 2
model = "Spin"
method = "CG"
lattice = "kagome"
J0 = 1.0
J1 = 0.5
J2 = 0.7
J'=0.2
2Sz = 1
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
  -37.8000000000 1.0000000000 -0.0411219716 -0.0011578391
  -22.6800000000 1.0000000000 -0.0715068358 -0.0035095854
  -7.5600000000 1.0000000000 -0.2703429083 -0.0528702157
  7.5600000000 1.0000000000 0.1482108549 -0.0152643621
  22.6800000000 1.0000000000 0.0585967259 -0.0023538947
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
-37.8000000000 1.0000000000 -0.0673363783 -0.0018720807
-22.6800000000 1.0000000000 -0.1160086954 -0.0055679433
-7.5600000000 1.0000000000 -0.4132288825 -0.0735196147
7.5600000000 1.0000000000 0.2566227982 -0.0275716409
22.6800000000 1.0000000000 0.0988745528 -0.0040411036
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
  -37.8000000000 1.0000000000 -0.1462975121 -0.0042842634
  -22.6800000000 1.0000000000 -0.2620517374 -0.0137721571
  -7.5600000000 1.0000000000 -1.2008968068 -0.3073197088
  7.5600000000 1.0000000000 0.4424175519 -0.0394580455
  22.6800000000 1.0000000000 0.1896070353 -0.0072005351
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste3.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%8.6f", diff}' paste3.dat`

test "${diff}" = "0.000000"

exit $?
