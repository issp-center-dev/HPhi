#!/bin/sh -e

mkdir -p lobcg_spinless/
cd lobcg_spinless
python3 "$1/test/testSpinlessCalc.py" -p "../../src/HPhi" -mpi "${MPIRUN}" -m "SpinlessFermion" -s 8 -t "LOBCG"

# Check value: flct
cat > reference.dat <<EOF
   0
   -4.8284271247461898
   0.0000000000000000
   0.0000000000000000

   1
   -4.8284271247461907
   0.0000000000000000
   0.0000000000000000

   2
   -3.4142135623730936
   0.0000000000000000
   0.0000000000000000

   3
   -3.4142135623730958
   0.0000000000000000
   0.0000000000000000

   4
   -3.4142135623730483
   0.0000000000000000
   0.0000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"
exit $?
