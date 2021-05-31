#!/bin/sh -e

mkdir -p spectrum_genspin_ladder/
cd spectrum_genspin_ladder
#
# Ground state
#
cat > stan1.in <<EOF
L = 3
W = 2
model = "Spin"
method = "CG"
lattice = "ladder"
J0 = 1.0
J1 = 1.0
2Sz = 0
2S = 3
EigenVecIO = out
SpectrumQW = 0.5
SpectrumQL = 0.3333333333333333333333
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
  -412.7076581450 1.0000000000 -0.0263199541 -0.0000664587
  -247.6245948870 1.0000000000 -0.0451329726 -0.0001954260
  -82.5415316290 1.0000000000 -0.1582470979 -0.0024036769
  82.5415316290 1.0000000000 0.1050755217 -0.0010594720
  247.6245948870 1.0000000000 0.0394408400 -0.0001492391
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
  -412.7076581450 1.0000000000 -0.0526399081 -0.0001329174
  -247.6245948870 1.0000000000 -0.0902659451 -0.0003908520
  -82.5415316290 1.0000000000 -0.3164941958 -0.0048073538
  82.5415316290 1.0000000000 0.2101510434 -0.0021189441
  247.6245948870 1.0000000000 0.0788816799 -0.0002984782
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste2.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} END{printf "%7.5f", diff}' paste2.dat`
echo ${diff}
test "${diff}" = "0.00000"

exit $?
