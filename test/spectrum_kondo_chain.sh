#!/bin/sh -e

mkdir -p spectrum_kondo_chain/
cd spectrum_kondo_chain
#
# Sz-Sz spectrum
#
cat > stan2.in <<EOF
L = 4
model = "Kondo"
method = "CG"
lattice = "chain"
t = 1.0
J = 4.0
nelec = 4
2Sz = 0
SpectrumQW = 0.5
SpectrumQL = 0.5
NOmega = 5
OmegaIm = 1.0
CalcSpec = "Scratch"
SpectrumType = "SzSz"
OmegaMax = 64.6776213781762834
Omegamin = -39.3223786218237095
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -52.0000000000 1.0000000000 -0.0027565296 -0.0000617734
  -31.2000000000 1.0000000000 -0.0051674668 -0.0002188790
  -10.4000000000 1.0000000000 -0.0359739617 -0.0399362291
  10.4000000000 1.0000000000 0.0071470861 -0.0004231172
  31.2000000000 1.0000000000 0.0032299691 -0.0000849182
EOF
paste output/zvo_DynamicalGreen_0.dat reference.dat > paste1.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} 
END{printf "%8.6f", diff}
' paste1.dat`
echo "Diff output/zvo_DynamicalGreen_0.dat (SzSz) : " ${diff}
test "${diff}" = "0.000000"
#
# S+S- spectrum
#
cat > stan2.in <<EOF
L = 4
model = "Kondo"
method = "CG"
lattice = "chain"
t = 1.0
J = 4.0
nelec = 4
2Sz = 0
SpectrumQW = 0.5
SpectrumQL = 0.5
NOmega = 5
OmegaIm = 1.0
CalcSpec = "Scratch"
SpectrumType = "S+S-"
OmegaMax = 64.6776213781762834
Omegamin = -39.3223786218237095
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
-52.0000000000 1.0000000000 -0.0055130593 -0.0001235468
-31.2000000000 1.0000000000 -0.0103349337 -0.0004377579
-10.4000000000 1.0000000000 -0.0719479233 -0.0798724508
10.4000000000 1.0000000000 0.0142941726 -0.0008462344
31.2000000000 1.0000000000 0.0064599384 -0.0001698363
EOF
paste output/zvo_DynamicalGreen_0.dat reference.dat > paste2.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} 
END{printf "%8.6f", diff}
' paste2.dat`
echo "Diff output/zvo_DynamicalGreen_0.dat (S+S-) : " ${diff}
test "${diff}" = "0.000000"
#
# Density-Density spectrum
#
cat > stan2.in <<EOF
L = 4
model = "Kondo"
method = "CG"
lattice = "chain"
t = 1.0
J = 4.0
nelec = 4
2Sz = 0
SpectrumQW = 0.5
SpectrumQL = 0.5
NOmega = 5
OmegaIm = 1.0
CalcSpec = "Scratch"
SpectrumType = "Density"
OmegaMax = 64.6776213781762834
Omegamin = -39.3223786218237095
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -52.0000000000 1.0000000000 -0.0090233292 -0.0001967967
  -31.2000000000 1.0000000000 -0.0165075426 -0.0006608521
  -10.4000000000 1.0000000000 -0.0989249107 -0.0289848456
  10.4000000000 1.0000000000 0.0252390332 -0.0015534635
  31.2000000000 1.0000000000 0.0111249204 -0.0002993371
EOF
paste output/zvo_DynamicalGreen_0.dat reference.dat > paste3.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} 
END{printf "%8.6f", diff}
' paste3.dat`
echo "Diff output/zvo_DynamicalGreen_0.dat (Density) : " ${diff}
test "${diff}" = "0.000000"
#
# Up-Up spectrum
#
cat > stan2.in <<EOF
L = 4
model = "Kondo"
method = "CG"
lattice = "chain"
t = 1.0
J = 4.0
nelec = 4
2Sz = 0
SpectrumQW = 0.5
SpectrumQL = 0.5
NOmega = 5
OmegaIm = 1.0
CalcSpec = "Scratch"
SpectrumType = "Up"
OmegaMax = 64.6776213781762834
Omegamin = -39.3223786218237095
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -52.0000000000 1.0000000000 -0.0316822725 -0.0007637764
  -31.2000000000 1.0000000000 -0.0634381695 -0.0030679778
  -10.4000000000 1.0000000000 0.1751530258 -1.2854110796
  10.4000000000 1.0000000000 0.0627602932 -0.0030036702
  31.2000000000 1.0000000000 0.0315100835 -0.0007555217
EOF
paste output/zvo_DynamicalGreen_0.dat reference.dat > paste4.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} 
END{printf "%8.6f", diff}
' paste4.dat`
echo "Diff output/zvo_DynamicalGreen_0.dat (Up) : " ${diff}
test "${diff}" = "0.000000"
#
# Down-Down spectrum
#
cat > stan2.in <<EOF
L = 4
model = "Kondo"
method = "CG"
lattice = "chain"
t = 1.0
J = 4.0
nelec = 4
2Sz = 0
SpectrumQW = 0.5
SpectrumQL = 0.5
NOmega = 5
OmegaIm = 1.0
CalcSpec = "Scratch"
SpectrumType = "Down"
OmegaMax = 64.6776213781762834
Omegamin = -39.3223786218237095
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -52.0000000000 1.0000000000 -0.0316822703 -0.0007637764
  -31.2000000000 1.0000000000 -0.0634381650 -0.0030679776
  -10.4000000000 1.0000000000 0.1751530125 -1.2854109874
  10.4000000000 1.0000000000 0.0627602890 -0.0030036700
  31.2000000000 1.0000000000 0.0315100813 -0.0007555216
EOF
paste output/zvo_DynamicalGreen_0.dat reference.dat > paste5.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($3-$7)*($3-$7))+sqrt(($4-$8)*($4-$8))} 
END{printf "%8.6f", diff}
' paste5.dat`
echo "Diff output/zvo_DynamicalGreen_0.dat (Down) : " ${diff}
test "${diff}" = "0.000000"

exit $?
