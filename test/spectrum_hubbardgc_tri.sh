#!/bin/sh -e

mkdir -p spectrum_hubbardgc_tri/
cd spectrum_hubbardgc_tri
#
# Sz-Sz spectrum
#
cat > stan2.in <<EOF
a0w = 3
a0l = 0
a1w = -1
a1l = 2
model = "HubbardGC"
method = "CG"
lattice = "Triangular"
t = 1.0
U = 4.0
h = 3.0
SpectrumQW = 0.5
SpectrumQL = 0.5
NOmega = 5
OmegaIm = 1.0
CalcSpec = "Scratch"
SpectrumType = "SzSz"
OmegaMax = 129.3140744040596246
Omegamin = -98.6859255959403754
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -114.0000000000 1.0000000000 -0.0077918862 -0.0000747870
  -68.4000000000 1.0000000000 -0.0138584520 -0.0002369117
  -22.8000000000 1.0000000000 -0.0637907099 -0.0052122099
  22.8000000000 1.0000000000 0.0251369120 -0.0007837300
  68.4000000000 1.0000000000 0.0104120987 -0.0001336193
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
a0w = 3
a0l = 0
a1w = -1
a1l = 2
model = "HubbardGC"
method = "CG"
lattice = "Triangular"
t = 1.0
U = 4.0
h = 3.0
SpectrumQW = 0.5
SpectrumQL = 0.5
NOmega = 5
OmegaIm = 1.0
CalcSpec = "Scratch"
SpectrumType = "S+S-"
OmegaMax = 129.3140744040596246
Omegamin = -98.6859255959403754
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
-114.0000000000 1.0000000000 -0.0287050987 -0.0002730922
-68.4000000000 1.0000000000 -0.0507169678 -0.0008542731
-22.8000000000 1.0000000000 -0.2233867683 -0.0173438683
22.8000000000 1.0000000000 0.0970120615 -0.0031631292
68.4000000000 1.0000000000 0.0392484221 -0.0005111302
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
a0w = 3
a0l = 0
a1w = -1
a1l = 2
model = "HubbardGC"
method = "CG"
lattice = "Triangular"
t = 1.0
U = 4.0
h = 3.0
SpectrumQW = 0.5
SpectrumQL = 0.5
NOmega = 5
OmegaIm = 1.0
CalcSpec = "Scratch"
SpectrumType = "Density"
OmegaMax = 129.3140744040596246
Omegamin = -98.6859255959403754
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -114.0000000000 1.0000000000 -0.0238695671 -0.0002257232
  -68.4000000000 1.0000000000 -0.0419778722 -0.0006994446
  -22.8000000000 1.0000000000 -0.1786004023 -0.0133148387
  22.8000000000 1.0000000000 0.0824900641 -0.0027218475
  68.4000000000 1.0000000000 0.0330773262 -0.0004338026
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
a0w = 3
a0l = 0
a1w = -1
a1l = 2
model = "HubbardGC"
method = "CG"
lattice = "Triangular"
t = 1.0
U = 4.0
h = 3.0
SpectrumQW = 0.5
SpectrumQL = 0.5
NOmega = 5
OmegaIm = 1.0
CalcSpec = "Scratch"
SpectrumType = "Down"
OmegaMax = 129.3140744040596246
Omegamin = -98.6859255959403754
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -114.0000000000 1.0000000000 -0.0002609286 -0.0000024571
  -68.4000000000 1.0000000000 -0.0004573225 -0.0000075562
  -22.8000000000 1.0000000000 -0.0018667417 -0.0001285394
  22.8000000000 1.0000000000 0.0009136847 -0.0000303575
  68.4000000000 1.0000000000 0.0003645766 -0.0000047999
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
a0w = 3
a0l = 0
a1w = -1
a1l = 2
model = "HubbardGC"
method = "CG"
lattice = "Triangular"
t = 1.0
U = 4.0
h = 3.0
SpectrumQW = 0.5
SpectrumQL = 0.5
NOmega = 5
OmegaIm = 1.0
CalcSpec = "Scratch"
SpectrumType = "Up"
OmegaMax = 129.3140744040596246
Omegamin = -98.6859255959403754
EOF

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -114.0000000000 1.0000000000 -0.0134927482 -0.0001338471
  -68.4000000000 1.0000000000 -0.0246344010 -0.0004464180
  -22.8000000000 1.0000000000 -0.1412488722 -0.0149296486
  22.8000000000 1.0000000000 0.0378925642 -0.0010593148
  68.4000000000 1.0000000000 0.0166857670 -0.0002047483
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
