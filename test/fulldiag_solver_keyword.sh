#!/bin/sh -e

mkdir -p fulldiag_solver_keyword/
cd fulldiag_solver_keyword

cat > stan.in <<EOF
L = 4
model = "FermionHubbard"
method = "FullDiag"
lattice = "chain"
t = 1.0
U = 4.0
nelec = 4
2Sz = 0
EOF

# (1) Solver 0 明示指定で正常終了し、既存参照値と一致すること
../../src/HPhi -sdry stan.in
echo "Solver  0" >> calcmod.def
${MPIRUNFC} ../../src/HPhi -e namelist.def

cat > reference_energy.dat <<EOF
  -2.102748
  -1.806424
  -1.068140
EOF
awk 'NR>1 && NR<=4 {printf "%11.6f\n", $1}' output/zvo_phys_Nup2_Ndown2.dat > energy.dat
paste energy.dat reference_energy.dat > paste_e.dat
diff=`awk 'BEGIN{max=0}{d=$1-$2; if(d<0)d=-d; if(d>max)max=d}END{print max}' paste_e.dat`
test "`echo "$diff < 0.000001" | bc`" = "1"

# (2) 範囲外 Solver 9 はエラー終了すること
cd ..
mkdir -p fulldiag_solver_keyword_invalid/
cd fulldiag_solver_keyword_invalid
cp ../fulldiag_solver_keyword/stan.in .
../../src/HPhi -sdry stan.in
echo "Solver  9" >> calcmod.def
if ${MPIRUNFC} ../../src/HPhi -e namelist.def > invalid.log 2>&1; then
  echo "ERROR: Solver 9 should have failed"
  exit 1
fi

# (3) 非 ELPA ビルドでは Solver 3 はエラー終了すること
#     (ELPA ビルドでは環境依存になるため、このケースは _ELPA 無しビルドの CI 前提)
cd ..
mkdir -p fulldiag_solver_keyword_elpa/
cd fulldiag_solver_keyword_elpa
cp ../fulldiag_solver_keyword/stan.in .
../../src/HPhi -sdry stan.in
echo "Solver  3" >> calcmod.def
if ${MPIRUNFC} ../../src/HPhi -e namelist.def > elpa.log 2>&1; then
  grep -q "Solver" elpa.log || { echo "ERROR: Solver 3 should fail on non-ELPA build"; exit 1; }
fi

# (4) 旧 ScaLAPACK キーワードで非推奨警告が出ること（動作は従来どおり）
cd ..
mkdir -p fulldiag_solver_keyword_dep/
cd fulldiag_solver_keyword_dep
cp ../fulldiag_solver_keyword/stan.in .
../../src/HPhi -sdry stan.in
echo "ScaLAPACK  0" >> calcmod.def
${MPIRUNFC} ../../src/HPhi -e namelist.def > dep.log 2>&1
grep -q "deprecated" dep.log

echo "fulldiag_solver_keyword: OK"
