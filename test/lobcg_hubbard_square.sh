#!/bin/sh -e

mkdir -p lobcg_hubbard_square/
cd lobcg_hubbard_square

cat > stan.in <<EOF
W = 4
L = 2
model = "FermionHubbard"
method = "CG"
lattice = "Tetragonal"
t = 1.0
U = 4.0
nelec = 8
2Sz = 0
exct = 3
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check Energy

cat > reference.dat <<EOF
   0
 -10.2529529552637033 
   1.0444442739280635 
   0.0000000000000000 

   1
  -9.8328986753266356 
   1.1830984216215232 
   0.0000000000000000 

   2
  -9.2310650024220635 
   1.2843076003322857 
   0.0000000000000000 
EOF
paste output/zvo_energy.dat reference.dat > paste1.dat
diff=`awk '
BEGIN{diff=0.0}
{diff+=sqrt(($2-$3)*($2-$3))}
END{printf "%8.6f", diff/NR}
' paste1.dat`
echo "Diff output/zvo_energy.dat : " ${diff}
test "${diff}" = "0.000000"

# Check one-body G

cat > reference.dat <<EOF
   0.5000000000 0.0000000000
   0.0601424447 -0.0000000000
   0.0000000000 0.0000000000
   0.0601424447 0.0000000000
   0.3908178694 0.0000000000
   -0.0000000000 -0.0000000000
   -0.0456141696 -0.0000000000
   -0.0000000000 0.0000000000
   0.5000000000 0.0000000000
   0.0601424447 -0.0000000000
   0.0000000000 0.0000000000
   0.0601424447 -0.0000000000
   0.3908178694 0.0000000000
   -0.0000000000 -0.0000000000
   -0.0456141696 0.0000000000
   -0.0000000000 -0.0000000000
EOF
paste output/zvo_cisajs_eigen0.dat reference.dat > paste2.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($5-$7)*($5-$7)+($6-$8)*($6-$8))}
END{printf "%8.6f", diff/NR}
' paste2.dat`
echo "Diff output/zvo_cisajs_eigen0.dat : " ${diff}
test "${diff}" = "0.000000"

cat > reference.dat <<EOF
   0.5000000000 0.0000000000
   0.1160851096 0.0000000000
   -0.0000000000 -0.0000000000
   0.1160851096 0.0000000000
   0.3390802767 -0.0000000000
   -0.0000000000 -0.0000000000
   -0.1005369505 0.0000000000
   0.0000000000 -0.0000000000
   0.5000000000 0.0000000000
   0.1160851096 -0.0000000000
   0.0000000000 0.0000000000
   0.1160851096 -0.0000000000
   0.3390802767 0.0000000000
   0.0000000000 -0.0000000000
   -0.1005369505 -0.0000000000
   0.0000000000 -0.0000000000
EOF
paste output/zvo_cisajs_eigen1.dat reference.dat > paste3.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($5-$7)*($5-$7)+($6-$8)*($6-$8))}
END{printf "%8.6f", diff/NR}
' paste3.dat`
echo "Diff output/zvo_cisajs_eigen1.dat : " ${diff}
test "${diff}" = "0.000000"

cat > reference.dat <<EOF
   0.4999999066 0.0000000000
   0.1715395136 -0.0000000131
   0.0000000029 0.0000000082
   0.1715395159 -0.0000000029
   0.2774697226 -0.0000000072
   0.0000000613 -0.0000000123
   -0.1530151560 0.0000000195
   -0.0000000038 0.0000000011
   0.5000000944 0.0000000000
   0.1715395159 0.0000000171
   0.0000000020 0.0000000109
   0.1715395044 0.0000000109
   0.2774697296 0.0000000056
   -0.0000000666 -0.0000000105
   -0.1530151585 -0.0000000253
   0.0000000026 -0.0000000023
EOF
paste output/zvo_cisajs_eigen2.dat reference.dat > paste4.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($5-$7)*($5-$7)+($6-$8)*($6-$8))}
END{printf "%8.6f", diff/NR}
' paste4.dat`
echo "Diff output/zvo_cisajs_eigen2.dat : " ${diff}
test "${diff}" = "0.000000"

# Check two-body G

cat > reference.dat <<EOF
   0.5000000000 0.0000000000
   0.1305555342 0.0000000000
   0.2199468210 0.0000000000
   0.2753184555 0.0000000000
   0.2661496850 0.0000000000
   0.2324888052 0.0000000000
   0.2199468210 0.0000000000
   0.2753184555 0.0000000000
   0.0319401468 0.0000000000
   0.3594749660 0.0000000000
   0.2677120490 0.0000000000
   0.2279368239 0.0000000000
   0.2265924281 0.0000000000
   0.2709701358 0.0000000000
   0.2677120490 0.0000000000
   0.2279368239 0.0000000000
   0.3694444658 0.0000000000
   -0.0553716345 -0.0000000000
   0.0336608798 0.0000000000
   -0.0553716345 -0.0000000000
   -0.3275348192 0.0000000000
   0.0397752252 -0.0000000000
   -0.0443777076 -0.0000000000
   0.0397752252 0.0000000000
   0.3694444658 0.0000000000
   -0.0553716345 0.0000000000
   0.0336608798 -0.0000000000
   -0.0553716345 0.0000000000
   -0.3275348192 -0.0000000000
   0.0397752252 0.0000000000
   -0.0443777076 0.0000000000
   0.0397752252 -0.0000000000
   0.1305555342 0.0000000000
   0.5000000000 0.0000000000
   0.2753184555 0.0000000000
   0.2199468210 0.0000000000
   0.2324888052 0.0000000000
   0.2661496850 0.0000000000
   0.2753184555 0.0000000000
   0.2199468210 0.0000000000
   0.3594749660 0.0000000000
   0.0319401468 0.0000000000
   0.2279368239 0.0000000000
   0.2677120490 0.0000000000
   0.2709701358 0.0000000000
   0.2265924281 0.0000000000
   0.2279368239 0.0000000000
   0.2677120490 0.0000000000
EOF
paste output/zvo_cisajscktalt_eigen0.dat reference.dat > paste5.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($9-$11)*($9-$11)+($10-$12)*($10-$12))}
END{printf "%8.6f", diff/NR}
' paste5.dat`
echo "Diff output/zvo_cisajscktalt_eigen0.dat : " ${diff}
test "${diff}" = "0.000000"

cat > reference.dat <<EOF
   0.5000000000 0.0000000000
   0.1478873027 0.0000000000
   0.1821147384 0.0000000000
   0.3040076167 0.0000000000
   0.3030897837 0.0000000000
   0.1910305808 0.0000000000
   0.1821147384 0.0000000000
   0.3040076167 0.0000000000
   0.0604424762 0.0000000000
   0.3518228573 0.0000000000
   0.2945170858 0.0000000000
   0.1961291187 0.0000000000
   0.1832040916 0.0000000000
   0.3089857884 0.0000000000
   0.2945170858 0.0000000000
   0.1961291187 0.0000000000
   0.3521126973 0.0000000000
   0.0098479115 0.0000000000
   0.0051803032 -0.0000000000
   0.0098479115 -0.0000000000
   -0.1586394158 0.0000000000
   0.0115241417 0.0000000000
   0.0086023088 0.0000000000
   0.0115241417 -0.0000000000
   0.3521126973 0.0000000000
   0.0098479115 -0.0000000000
   0.0051803032 0.0000000000
   0.0098479115 0.0000000000
   -0.1586394158 -0.0000000000
   0.0115241417 -0.0000000000
   0.0086023088 -0.0000000000
   0.0115241417 0.0000000000
   0.1478873027 0.0000000000
   0.5000000000 0.0000000000
   0.3040076167 0.0000000000
   0.1821147384 0.0000000000
   0.1910305808 0.0000000000
   0.3030897837 0.0000000000
   0.3040076167 0.0000000000
   0.1821147384 0.0000000000
   0.3518228573 0.0000000000
   0.0604424762 0.0000000000
   0.1961291187 0.0000000000
   0.2945170858 0.0000000000
   0.3089857884 0.0000000000
   0.1832040916 0.0000000000
   0.1961291187 0.0000000000
   0.2945170858 0.0000000000
EOF
paste output/zvo_cisajscktalt_eigen1.dat reference.dat > paste6.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($9-$11)*($9-$11)+($10-$12)*($10-$12))}
END{printf "%8.6f", diff/NR}
' paste6.dat`
echo "Diff output/zvo_cisajscktalt_eigen1.dat : " ${diff}
test "${diff}" = "0.000000"

cat > reference.dat <<EOF
   0.4999999066 0.0000000000
   0.1605384504 0.0000000000
   0.2124152216 0.0000000000
   0.2586285289 0.0000000000
   0.2238715694 0.0000000000
   0.2741542003 0.0000000000
   0.2124152236 0.0000000000
   0.2586285243 0.0000000000
   0.1577326825 0.0000000000
   0.2724843505 0.0000000000
   0.2391462242 0.0000000000
   0.2542761967 0.0000000000
   0.2152724814 0.0000000000
   0.2670132753 0.0000000000
   0.2391463172 0.0000000000
   0.2542761001 0.0000000000
   0.3394614562 0.0000000000
   -0.0462133180 0.0000000055
   -0.0502827226 0.0000000264
   -0.0462133136 -0.0000000006
   -0.1147517021 0.0000000060
   -0.0151298660 0.0000000191
   -0.0517407609 0.0000000087
   -0.0151298670 0.0000000142
   0.3394616440 0.0000000000
   -0.0462133180 -0.0000000055
   -0.0502827226 -0.0000000264
   -0.0462133136 0.0000000006
   -0.1147517021 -0.0000000060
   -0.0151298660 -0.0000000191
   -0.0517407609 -0.0000000087
   -0.0151298670 -0.0000000142
   0.1605384504 0.0000000000
   0.5000000944 0.0000000000
   0.2586286345 0.0000000000
   0.2124153023 0.0000000000
   0.2741543870 0.0000000000
   0.2238715708 0.0000000000
   0.2586286316 0.0000000000
   0.2124153116 0.0000000000
   0.2724844670 0.0000000000
   0.1577327486 0.0000000000
   0.2542761847 0.0000000000
   0.2391464216 0.0000000000
   0.2670133388 0.0000000000
   0.2152726041 0.0000000000
   0.2542762837 0.0000000000
   0.2391463242 0.0000000000
EOF
paste output/zvo_cisajscktalt_eigen2.dat reference.dat > paste7.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($9-$11)*($9-$11)+($10-$12)*($10-$12))}
END{printf "%8.6f", diff/NR}
' paste7.dat`
echo "Diff output/zvo_cisajscktalt_eigen2.dat : " ${diff}
test "${diff}" = "0.000000"

exit $?