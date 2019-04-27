#!/bin/sh -e

mkdir -p lanczos_hubbard_square_restart/
cd lanczos_hubbard_square_restart

cat > stan.in <<EOF
W = 4
L = 2
model = "FermionHubbard"
method = "Lanczos"
Lanczos_max = 10
lattice = "Tetragonal"
t = 1.0
U = 4.0
nelec = 8
2Sz = 0
outputmode = "None"
Restart = "Restart_out"
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value for Restart_in
cp ./stan.in stan2.in
sed -i -e "s/Restart = \"Restart_out\"/Restart = \"Restart_in\"/g" stan2.in
sed -i -e "s/Lanczos_max = 10/Lanczos_max = 1000/g" stan2.in 

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
  -10.2529529552637637
    1.0444442739280932
    0.0000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste1.dat`

# Check value for Restart_inout
rm output/zvo_energy.dat
cp ./stan.in stan3.in
sed -i -e "s/Restart = \"Restart_out\"/Restart = \"Restart\"/g" stan3.in
sed -i -e "s/Lanczos_max = 10/Lanczos_max = 1000/g" stan3.in 

${MPIRUN} ../../src/HPhi -s stan3.in

cat > reference.dat <<EOF
  -10.2529529552637637
    1.0444442739280932
    0.0000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste2.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste2.dat`

test "${diff}" = "0.000000"
exit $?
