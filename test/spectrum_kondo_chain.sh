#!/bin/sh -e

mkdir -p spectrum_kondo_chain/
cd spectrum_kondo_chain
#
# Ground state
#
cat > stan1.in <<EOF
L = 4
model = "Kondo"
method = "CG"
lattice = "chain"
t = 1.0
J = 4.0
nelec = 4
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
  -52.0000000000 1.0000000000 -0.0027565296 -0.0000617734
  -31.2000000000 1.0000000000 -0.0051674668 -0.0002188790
  -10.4000000000 1.0000000000 -0.0359739617 -0.0399362291
  10.4000000000 1.0000000000 0.0071470861 -0.0004231172
  31.2000000000 1.0000000000 0.0032299691 -0.0000849182
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
  -52.0000000000 1.0000000000 -0.0032805430 -0.0000666850
  -31.2000000000 1.0000000000 -0.0057050223 -0.0002054512
  -10.4000000000 1.0000000000 -0.0261701586 -0.0177415792
  10.4000000000 1.0000000000 0.0149006071 -0.0016360066
  31.2000000000 1.0000000000 0.0049769844 -0.0001548985
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
  -52.0000000000 1.0000000000 -0.0090233292 -0.0001967967
  -31.2000000000 1.0000000000 -0.0165075426 -0.0006608521
  -10.4000000000 1.0000000000 -0.0989249107 -0.0289848456
  10.4000000000 1.0000000000 0.0252390332 -0.0015534635
  31.2000000000 1.0000000000 0.0111249204 -0.0002993371
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste3.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)**2)+sqrt(($4-$8)**2)} END{printf "%8.6f", diff}' paste3.dat`
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
  -52.0000000000 1.0000000000 -0.0316822725 -0.0007637764
  -31.2000000000 1.0000000000 -0.0634381695 -0.0030679778
  -10.4000000000 1.0000000000 0.1751530258 -1.2854110796
  10.4000000000 1.0000000000 0.0627602932 -0.0030036702
  31.2000000000 1.0000000000 0.0315100835 -0.0007555217
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste4.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)**2)+sqrt(($4-$8)**2)} END{printf "%8.6f", diff}' paste4.dat`
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
  -52.0000000000 1.0000000000 -0.0316822703 -0.0007637764
  -31.2000000000 1.0000000000 -0.0634381650 -0.0030679776
  -10.4000000000 1.0000000000 0.1751530125 -1.2854109874
  10.4000000000 1.0000000000 0.0627602890 -0.0030036700
  31.2000000000 1.0000000000 0.0315100813 -0.0007555216
EOF
paste output/zvo_DynamicalGreen.dat reference.dat > paste5.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($3-$7)**2)+sqrt(($4-$8)**2)} END{printf "%8.6f", diff}' paste5.dat`

test "${diff}" = "0.000000"

exit $?
