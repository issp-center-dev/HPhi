#!/bin/sh -e

mkdir -p lanczos_spingc_hcor/
cd lanczos_spingc_hcor

cat > stan.in <<EOF
W = 2
L = 2
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
h=-1.0
EOF

cat > green3.def <<EOF
===================
num       2
===================
===================
===================
        0         0           0           1           2           1           2           0           3           0           3           0   
        0         0           0           1           2           1           2           0           3           1           3           1   
EOF

cat > green4.def <<EOF
===================
num       2
===================
===================
===================
        0  0  0  1  2  1  2  0  3  0  3  0  4  0  4  0  
        0  0  0  1  2  1  2  0  3  1  3  1  4  1  4  1 
EOF

../../src/HPhi -sdry stan.in

echo   "ThreeBodyG  green3.def" >> namelist.def
echo   "FourBodyG   green4.def" >> namelist.def

${MPIRUN} ../../src/HPhi -e namelist.def

# Check value

cat > reference.dat <<EOF
  -5.1653788071251903 
   0.0000000000000000 
   3.8905991010945784 
EOF

paste output/zvo_energy.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste1.dat`

cat > reference.dat <<EOF
    0    0    0    1    2    1    2    0    3    0    3    0 -0.0067383864 0.0000000000 
    0    0    0    1    2    1    2    0    3    1    3    1 -0.0067383864 0.0000000000 
EOF
paste output/zvo_ThreeBody.dat reference.dat > paste2.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($13-$27)*($13-$27)+($14-$28)*($14-$28))} 
END{printf "%8.6f", diff}' paste2.dat`

cat > reference.dat <<EOF
    0    0    0    1    2    1    2    0    3    0    3    0     4    0    4    0 -0.0001083268 -0.0000000000 
    0    0    0    1    2    1    2    0    3    1    3    1     4    1    4    1 -0.0065792883 0.0000000000
EOF
paste output/zvo_FourBody.dat reference.dat > paste3.dat
diff=`awk 'BEGIN{diff='${diff}'} {diff+=sqrt(($17-$35)*($17-$35)+($18-$36)*($18-$36))} 
END{printf "%8.6f", diff}' paste3.dat`

test "${diff}" = "0.000000"
exit $?
