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
  cat invalid.log
  exit 1
fi

# (2b) Solver is a FullDiag backend selector. An explicit distributed solver
#      on Lanczos must be rejected before capability checks or MPI setup;
#      otherwise iFlgScaLAPACK disables site decomposition even though the
#      requested backend is never called.
cd ..
mkdir -p solver_keyword_nonfulldiag/
cd solver_keyword_nonfulldiag
cat > stan.in <<EOF
L = 4
model = "Spin"
method = "Lanczos"
lattice = "chain"
J = 1.0
2Sz = 0
EOF
../../src/HPhi -sdry stan.in
echo "Solver  1" >> calcmod.def
if ../../src/HPhi -e namelist.def > nonfulldiag.log 2>&1; then
  echo "ERROR: Solver 1 with Lanczos should have failed"
  cat nonfulldiag.log
  exit 1
fi
grep -q "only valid with CalcType=2" nonfulldiag.log || {
  echo "ERROR: non-FullDiag Solver failure should explain FullDiag eligibility"
  cat nonfulldiag.log
  exit 1
}

# (2c) The deprecated ScaLAPACK keyword must remain backward compatible.
#      Outside FullDiag, a ScaLAPACK-enabled build warns and ignores it so the
#      normal MPI site decomposition remains enabled; a non-ScaLAPACK build
#      keeps its longstanding behavior of ignoring the unavailable flag.
cd ..
mkdir -p scalapack_keyword_nonfulldiag/
cd scalapack_keyword_nonfulldiag
cp ../solver_keyword_nonfulldiag/stan.in .
../../src/HPhi -sdry stan.in
echo "ScaLAPACK  1" >> calcmod.def
if [ "${HPHI_HAS_SCALAPACK:-0}" = "1" ]; then
  # Use the configured CI launcher when available so this also proves that
  # clearing iFlgScaLAPACK restores multi-rank site decomposition.
  if ! ${MPIRUN} ../../src/HPhi -e namelist.def > nonfulldiag.log 2>&1; then
    echo "ERROR: legacy ScaLAPACK 1 with Lanczos should remain runnable"
    cat nonfulldiag.log
    exit 1
  fi
  grep -q 'legacy "ScaLAPACK" keyword applies only to CalcType=2 and is ignored' nonfulldiag.log || {
    echo "ERROR: legacy ScaLAPACK fallback should explain that the keyword was ignored"
    cat nonfulldiag.log
    exit 1
  }
else
  if ! ../../src/HPhi -e namelist.def > nonfulldiag.log 2>&1; then
    echo "ERROR: a non-ScaLAPACK build should preserve legacy keyword fallback"
    cat nonfulldiag.log
    exit 1
  fi
fi

# (3) capability-aware: 非 ELPA ビルドでは Solver 3 はエラー終了すること。
#     ELPA ビルド (HPHI_HAS_ELPA=1, test/CMakeLists.txt が USE_ELPA から設定)
#     では Solver 3 は正当な指定になるため、この serial (nproc==1) 環境では
#     CPU-ELPA 経路 (NGPU 0 を明示) が正常終了することを assert する。
#     非 ELPA ビルドでは、従来どおり「失敗すること」自体を assert する: 以前の
#     実装は run が成功してもログに "Solver" が無い場合しか失敗にしなかった
#     ため、Solver 3 が非 ELPA ビルドで誤って成功しても (ログに偶然 "Solver"
#     が出れば) テストが green のまま通ってしまうバグがあった (レビュー指摘)。
cd ..
mkdir -p fulldiag_solver_keyword_elpa/
cd fulldiag_solver_keyword_elpa
cp ../fulldiag_solver_keyword/stan.in .
../../src/HPhi -sdry stan.in
echo "Solver  3" >> calcmod.def
if [ "${HPHI_HAS_ELPA:-0}" = "1" ]; then
  echo "NGPU    0" >> calcmod.def
  if ! ${MPIRUNFC} ../../src/HPhi -e namelist.def > elpa.log 2>&1; then
    echo "ERROR: Solver 3 (NGPU 0) should succeed on an ELPA build"
    cat elpa.log
    exit 1
  fi
else
  if ${MPIRUNFC} ../../src/HPhi -e namelist.def > elpa.log 2>&1; then
    echo "ERROR: Solver 3 should fail on non-ELPA build (HPhi succeeded instead)"
    cat elpa.log
    exit 1
  fi
  grep -Eq "Solver|ELPA" elpa.log || {
    echo "ERROR: Solver 3 failure log should mention Solver/ELPA"
    cat elpa.log
    exit 1
  }
fi

# (3b) capability-aware: 非 ScaLAPACK ビルドでは Solver 1 も同様にエラー終了
#     すること (readdef.c ResolveSolver() の #ifndef _SCALAPACK ガード, 既存
#     動作の確認)。ScaLAPACK ビルド (HPHI_HAS_SCALAPACK=1) では Solver 1 は
#     正当な指定であり、この serial (nproc==1) 環境では lapack_diag() が
#     ZHEEVall にフォールバックするため正常終了することを assert する。
cd ..
mkdir -p fulldiag_solver_keyword_scalapack/
cd fulldiag_solver_keyword_scalapack
cp ../fulldiag_solver_keyword/stan.in .
../../src/HPhi -sdry stan.in
echo "Solver  1" >> calcmod.def
if [ "${HPHI_HAS_SCALAPACK:-0}" = "1" ]; then
  if ! ${MPIRUNFC} ../../src/HPhi -e namelist.def > scalapack.log 2>&1; then
    echo "ERROR: Solver 1 should succeed serially on a ScaLAPACK build (falls back to ZHEEVall)"
    cat scalapack.log
    exit 1
  fi
else
  if ${MPIRUNFC} ../../src/HPhi -e namelist.def > scalapack.log 2>&1; then
    echo "ERROR: Solver 1 should fail on non-ScaLAPACK build (HPhi succeeded instead)"
    cat scalapack.log
    exit 1
  fi
  grep -Eq "Solver|ScaLAPACK" scalapack.log || {
    echo "ERROR: Solver 1 failure log should mention Solver/ScaLAPACK"
    cat scalapack.log
    exit 1
  }
fi

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

# (5) ELPAビルドでは Solver 3 + OutputHam はマルチプロセス時に拒否される
#     （シリアル実行では nproc==1 なので拒否されない = ここでは
#     「calcmod のパースが通り実行が成功する」ことだけを確認し、
#     マルチプロセス拒否そのものは fulldiag_elpa_hubbard_chain 側の
#     ケースで検証する）
if [ "${HPHI_HAS_ELPA:-0}" = "1" ]; then
  cd ..
  mkdir -p fulldiag_solver_keyword_hamio/
  cd fulldiag_solver_keyword_hamio
  cp ../fulldiag_solver_keyword/stan.in .
  ../../src/HPhi -sdry stan.in
  printf "Solver  3\nNGPU  0\nOutputHam  1\n" >> calcmod.def
  ${MPIRUNFC} ../../src/HPhi -e namelist.def > hamio.log 2>&1
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

# (6) ExpecMode 1 は非分散ソルバー（Solver 0）ではエラー終了すること
cd ..
mkdir -p fulldiag_solver_keyword_expecmode_invalid/
cd fulldiag_solver_keyword_expecmode_invalid
cp ../fulldiag_solver_keyword/stan.in .
../../src/HPhi -sdry stan.in
printf "Solver  0\nExpecMode  1\n" >> calcmod.def
if ${MPIRUNFC} ../../src/HPhi -e namelist.def > expecmode_invalid.log 2>&1; then
  echo "ERROR: ExpecMode 1 with Solver 0 should have failed"
  cat expecmode_invalid.log
  exit 1
fi
grep -q "requires CalcType = FullDiag" expecmode_invalid.log || {
  echo "ERROR: ExpecMode failure log should mention the CalcType/Solver eligibility requirement"
  cat expecmode_invalid.log
  exit 1
}

# (7) capability-aware: ELPA ビルド（HPHI_HAS_ELPA=1）でのみ、Solver 3 +
#     ExpecMode 1 をシリアル実行（MPIRUNFC 空 = nproc 1）すると、readdef.c の
#     nproc==1 降格規則により INFO を出して ExpecMode 0 として正常終了し、
#     結果が ExpecMode 0（ケース(1)の参照値）と一致することを確認する。
#     非 ELPA ビルドではこのブロックはスキップされる（Solver 3 がケース(3)の
#     とおりビルドエラーになるため）。
if [ "${HPHI_HAS_ELPA:-0}" = "1" ]; then
  cd ..
  mkdir -p fulldiag_solver_keyword_expecmode/
  cd fulldiag_solver_keyword_expecmode
  cp ../fulldiag_solver_keyword/stan.in .
  ../../src/HPhi -sdry stan.in
  printf "Solver  3\nNGPU  0\nExpecMode  1\n" >> calcmod.def
  if ! ${MPIRUNFC} ../../src/HPhi -e namelist.def > expecmode.log 2>&1; then
    echo "ERROR: Solver 3 + ExpecMode 1 should succeed serially (downgrades to ExpecMode 0)"
    cat expecmode.log
    exit 1
  fi
  grep -q "ExpecMode reverts to 0 for a single process (results are identical)." expecmode.log || {
    echo "ERROR: expected the ExpecMode single-process downgrade message"
    cat expecmode.log
    exit 1
  }
  awk 'NR>1 && NR<=4 {printf "%11.6f\n", $1}' output/zvo_phys_Nup2_Ndown2.dat > energy.dat
  paste energy.dat ../fulldiag_solver_keyword/reference_energy.dat > paste_e.dat
  diff=`awk 'BEGIN{max=0}{d=$1-$2; if(d<0)d=-d; if(d>max)max=d}END{print max}' paste_e.dat`
  test "`echo "$diff < 0.000001" | bc`" = "1"
fi

echo "fulldiag_solver_keyword: OK"
