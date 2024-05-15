#!/bin/sh -e

mkdir -p lobcg_genspin_chain_threebody/
cp  $1/test/GenSpinThreeBody.py ./lobcg_genspin_chain_threebody
cd  ./lobcg_genspin_chain_threebody
  python3 GenSpinThreeBody.py generate
  ../../src/HPhi -sdry stan.in
  ${MPIRUN} ../../src/HPhi -e open_namelist.def
  python3 GenSpinThreeBody.py aft

cat > reference.dat <<EOF
      1.000000000000
      -9.537535583289
      -0.658624428716
      -0.607832852032
      -0.735737997800
      -0.292658073300 
EOF
paste result_hphi_size6.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($1-$2)*($1-$2))} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"
exit $?
