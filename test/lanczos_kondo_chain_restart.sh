#!/bin/sh -e

mkdir -p lanczos_kondo_chain_restart/
cd lanczos_kondo_chain_restart

cat > stan.in <<EOF
L = 4
model = "Kondo"
method = "Lanczos"
lattice = "chain"
Lanczos_max = 10
t = 1.0
J = 4.0
nelec = 4
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
# Check value

cat > reference.dat <<EOF
  -12.6776213781764451
    0.1072445462409946
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
  -12.6776213781764451
    0.1072445462409946
    0.0000000000000000
EOF

paste output/zvo_energy.dat reference.dat > paste2.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste2.dat`

test "${diff}" = "0.000000"
exit $?
