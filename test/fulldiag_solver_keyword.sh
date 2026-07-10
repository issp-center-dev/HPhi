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
#     ここでは「失敗すること」自体を assert する: 以前の実装は run が成功しても
#     ログに "Solver" が無い場合しか失敗にしなかったため、Solver 3 が非 ELPA
#     ビルドで誤って成功しても (ログに偶然 "Solver" が出れば) テストが green
#     のまま通ってしまうバグがあった (レビュー指摘)。
cd ..
mkdir -p fulldiag_solver_keyword_elpa/
cd fulldiag_solver_keyword_elpa
cp ../fulldiag_solver_keyword/stan.in .
../../src/HPhi -sdry stan.in
echo "Solver  3" >> calcmod.def
if ${MPIRUNFC} ../../src/HPhi -e namelist.def > elpa.log 2>&1; then
  echo "ERROR: Solver 3 should fail on non-ELPA build (HPhi succeeded instead)"
  exit 1
fi
grep -Eq "Solver|ELPA" elpa.log || {
  echo "ERROR: Solver 3 failure log should mention Solver/ELPA"
  exit 1
}

# (3b) 非 ScaLAPACK ビルドでは Solver 1 も同様にエラー終了すること
#     (readdef.c ResolveSolver() の #ifndef _SCALAPACK ガード, 既存動作の確認)
cd ..
mkdir -p fulldiag_solver_keyword_scalapack/
cd fulldiag_solver_keyword_scalapack
cp ../fulldiag_solver_keyword/stan.in .
../../src/HPhi -sdry stan.in
echo "Solver  1" >> calcmod.def
if ${MPIRUNFC} ../../src/HPhi -e namelist.def > scalapack.log 2>&1; then
  echo "ERROR: Solver 1 should fail on non-ScaLAPACK build (HPhi succeeded instead)"
  exit 1
fi
grep -Eq "Solver|ScaLAPACK" scalapack.log || {
  echo "ERROR: Solver 1 failure log should mention Solver/ScaLAPACK"
  exit 1
}

# NOTE (multi-rank gate coverage gap, documented per whole-branch review FIX 1):
# ResolveSolver() in src/readdef.c derives X->iFlgScaLAPACK=1 for Solver 1/3 so
# that the pre-existing multi-process FullDiag gates key on it correctly:
#   - src/HPhiMain.c:~720 (iCalcType==FullDiag && iFlgScaLAPACK==0 && nproc!=1
#     -> error exit)
#   - src/check.c:~100 and src/CheckMPI.c:~542 (iFlgScaLAPACK==1 -> replicated
#     Hilbert-space treatment, NsiteMPI=Nsite, no site separation)
# This machine's regression build (build_noMPI, ENABLE_MPI=OFF) always has
# nproc==1, so the nproc!=1 branch of the HPhiMain.c gate can never be
# exercised here, and Solver 3 additionally requires _ELPA which this build
# does not have (case 3 above). Actual multi-rank gate passage for Solver 1/3
# is covered by fulldiag_elpa_hubbard_chain.sh (registered min:2) on an
# ELPA+MPI build; do not attempt to fake a multi-rank run in this serial-only
# environment.

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
