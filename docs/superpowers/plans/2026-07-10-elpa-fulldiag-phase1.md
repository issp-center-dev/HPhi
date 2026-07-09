# ELPA FullDiag フェーズ1（ソルバー接続）実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** HPhi の FullDiag に `Solver` キーワードと ELPA バックエンド（CPU/GPU、2Dブロックサイクリック分散、フェーズ1では入力行列は既存の全複製 `Ham`）を追加する。

**Architecture:** 設計文書 `docs/superpowers/specs/2026-07-10-elpa-fulldiag-design.md`（v4、承認済み）のフェーズ1。`readdef.c` に `Solver` キーワードと優先順位解決を実装し、`lapack_diag.c` を `iSolver` ベースの分岐に書き換え、新モジュール `matrixlapack_elpa.c` に ELPA 呼び出しを隔離、`matrixscalapack.c` に固有ベクトルのブロック回収 `GetEigenVectorBlock` を追加する。既存経路（LAPACK/ScaLAPACK/MAGMA）の動作は不変。

**Tech Stack:** C99, MPI, BLACS/ScaLAPACK (MKL), ELPA C API (>= 2021.11; GPU は >= 2023.11.001), CMake, ctest + shell テスト。

## Global Constraints

- 設計文書: `docs/superpowers/specs/2026-07-10-elpa-fulldiag-design.md`（判断に迷ったらこれが正）
- ELPA 最低バージョン: CPU = API 20211125（`elpa_init(20211125)`）、GPU = 2023.11.001（`elpa_setup_gpu` の存在で判定、CMake マクロ `_ELPA_GPU`）
- スレッド版 `elpa_openmp` はサポートしない（検出したら configure エラー）
- GPU 要求（`NGPU >= 1`）失敗時に黙って CPU にフォールバックしない（明示エラーで停止）
- 全ての ELPA 集団操作（`elpa_setup` / `elpa_setup_gpu` / `elpa_eigenvectors` / `elpa_deallocate`）の直前に `MPI_Allreduce(MAX)` エラー同期
- 警告・エラー表示は rank 0 のみ（`stdoutMPI` 慣例）
- 既存経路 `Solver 0/1/2` の挙動・出力を一切変えない（`Solver 1` の `GetEigenVector` も変更しない）
- ELPA ブロックサイズは固定 64、ELPA ソルバーは GPU=`ELPA_SOLVER_1STAGE` / CPU=`ELPA_SOLVER_2STAGE`
- C から呼ぶ ScaLAPACK/BLACS ルーチンは全引数ポインタ渡し（既存 `matrixscalapack.c` の慣例）
- コミットメッセージ末尾に既存セッションの Co-Authored-By 行を付ける（このリポジトリでの本セッションの慣例）

## 開発環境の注意

ローカル開発機（macOS）には ELPA が無い前提で進める。各タスクは
(a) `USE_ELPA=OFF`（既定）でのビルドとテストが常に通ること、
(b) ELPA 依存コードはコンパイルガード（`_ELPA`）内に閉じることを守る。
ELPA 実機テスト（clavius）はタスク8のチェックリストとして分離してある。
ローカルでのビルド確認コマンド（以後「標準ビルド確認」と呼ぶ）:

```bash
cd /Users/k-yoshimi/Dropbox/CLionProjects/HPhi-box/HPhi
mkdir -p build && cd build && cmake .. > cmake.log 2>&1 && make HPhi -j4 2>&1 | tail -5
```

Expected: `[100%] Built target HPhi`（警告は既存分のみ）

---

### Task 1: `Solver` キーワード（パース・優先順位解決・検証）

**Files:**
- Modify: `src/include/DefCommon.h`（SOLVER_* 定数）
- Modify: `src/include/struct.h:322-336`（`iSolver` ほかフィールド追加）
- Modify: `src/readdef.c:247-420`（初期化・パース・解決・検証）
- Modify: `src/include/ErrorMessage.h` / `src/ErrorMessage.c`（メッセージ追加）
- Test: `test/fulldiag_solver_keyword.sh`、`test/CMakeLists.txt`

**Interfaces:**
- Produces: `X->Def.iSolver`（`SOLVER_LAPACK=0 / SOLVER_SCALAPACK=1 / SOLVER_MAGMA=2 / SOLVER_ELPA=3`、`ReadcalcmodFile` 完了後は必ず解決済み）。後続タスクは `iSolver` のみを見る（`iNGPU`/`iFlgScaLAPACK` を直接分岐に使わない）。

- [ ] **Step 1: 失敗するテストを書く**

`test/fulldiag_solver_keyword.sh` を新規作成（実行権限 `chmod +x` を忘れない）:

```sh
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
```

`test/CMakeLists.txt` の FullDiag テスト群（`add_hphi_test(fulldiag_hubbard_chain)` の近く）に追加:

```cmake
add_hphi_test(fulldiag_solver_keyword)
```

- [ ] **Step 2: テストが失敗することを確認**

```bash
cd build && cmake .. > /dev/null && make HPhi -j4 > /dev/null 2>&1 && ctest -R fulldiag_solver_keyword --output-on-failure
```

Expected: FAIL（現状 `calcmod.def` の未知キーワード `Solver` で HPhi がエラー終了するため、ケース (1) が落ちる）

- [ ] **Step 3: 定数とフィールドを追加**

`src/include/DefCommon.h` の `#define NUM_CALCTYPE 6` の下に追加:

```c
/*!< FullDiag solver backend (CalcMod keyword "Solver") */
#define NUM_SOLVER 4
#define SOLVER_LAPACK 0
#define SOLVER_SCALAPACK 1
#define SOLVER_MAGMA 2
#define SOLVER_ELPA 3
```

`src/include/struct.h` の `iFlgScaLAPACK` フィールド定義の直後に追加:

```c
    int iSolver;/**<@brief FullDiag solver backend (CalcMod keyword "Solver")
    - 0: LAPACK zheev (serial)
    - 1: ScaLAPACK pzheev
    - 2: MAGMA (single node multi-GPU)
    - 3: ELPA (multi-node CPU/GPU)
    Resolved from legacy keywords when not explicitly given. */

    int iFlgSolverSpec;/**<@brief 1 if the Solver keyword was explicitly given */

    int iFlgNGPUSpec;/**<@brief 1 if the NGPU keyword was explicitly given */
```

- [ ] **Step 4: エラーメッセージを追加**

`src/ErrorMessage.c` の `cErrCUDA` 定義の直後に追加:

```c
char *cErrSolver="Error in %s\n Solver: must be one of 0 (LAPACK), 1 (ScaLAPACK), 2 (MAGMA), 3 (ELPA).\n";
char *cErrSolverBuild="Error in %s\n Solver %d requires HPhi built with %s.\n";
char *cErrElpaGPUBuild="Error in %s\n Solver 3 with NGPU >= 1 requires ELPA >= 2023.11.001 (elpa_setup_gpu).\n Rebuild HPhi against a newer ELPA, or set \"NGPU 0\" for CPU execution.\n";
char *cWarnScaLAPACKDep="Warning in %s\n The \"ScaLAPACK\" keyword is deprecated. Use \"Solver 1\" instead.\n";
char *cWarnSolverConflict="Warning in %s\n Legacy keyword \"%s\" conflicts with the explicit Solver value and is ignored.\n";
```

`src/include/ErrorMessage.h` の `extern char *cErrCUDA;` の直後に対応する extern 宣言を追加:

```c
extern char *cErrSolver;
extern char *cErrSolverBuild;
extern char *cErrElpaGPUBuild;
extern char *cWarnScaLAPACKDep;
extern char *cWarnSolverConflict;
```

- [ ] **Step 5: パースと解決を実装**

`src/readdef.c` の `ReadcalcmodFile` 冒頭の初期化ブロック（`X->iNGPU` 初期化の直後、readdef.c:278 付近）に追加:

```c
  X->iSolver = -1;      /* unresolved; fixed up by ResolveSolver() below */
  X->iFlgSolverSpec = 0;
  X->iFlgNGPUSpec = 0;
```

キーワードパースの `NGPU` 節（readdef.c:333-335）を次に変更:

```c
    else if(CheckWords(ctmp, "NGPU")==0){
        X->iNGPU=itmp;
        X->iFlgNGPUSpec=1;
    }
```

その直後に `Solver` 節を追加:

```c
    else if(CheckWords(ctmp, "Solver")==0){
        X->iSolver=itmp;
        X->iFlgSolverSpec=1;
    }
```

`ScaLAPACK` 節（readdef.c:336-340）を次に変更（警告は全ビルドで出す）:

```c
    else if(CheckWords(ctmp, "ScaLAPACK")==0){
      fprintf(stdoutMPI, cWarnScaLAPACKDep, defname);
#ifdef _SCALAPACK
      X->iFlgScaLAPACK=itmp;
#endif
    }
```

`ReadcalcmodFile` の直前（ファイルスコープ）に解決関数を追加:

```c
/**
 * @brief Resolve the FullDiag solver backend from the Solver keyword and
 * legacy keywords (ScaLAPACK, NGPU), following the precedence table in
 * docs/superpowers/specs/2026-07-10-elpa-fulldiag-design.md section 2.
 * On return X->iSolver is one of SOLVER_*, and X->iNGPU has its
 * solver-dependent default applied when NGPU was not explicitly given.
 */
static int ResolveSolver(struct DefineList *X, const char *defname) {
  if (X->iFlgSolverSpec == 0) {
    /* Legacy resolution: preserve current behavior exactly.
       Compile-time NGPU default (2 on _MAGMA builds) applies here. */
    if (X->iNGPU > 0) {
#ifdef _MAGMA
      X->iSolver = SOLVER_MAGMA;
#else
      fprintf(stdoutMPI, "Warning: MAGMA is not used in this calculation.");
      X->iSolver = SOLVER_LAPACK;
#endif
    }
    else if (X->iFlgScaLAPACK == 1) {
      X->iSolver = SOLVER_SCALAPACK;
    }
    else {
      X->iSolver = SOLVER_LAPACK;
    }
    return 0;
  }

  /* Explicit Solver always wins; warn about conflicting legacy keywords. */
  if (ValidateValue(X->iSolver, 0, NUM_SOLVER - 1)) {
    fprintf(stdoutMPI, cErrSolver, defname);
    return -1;
  }
  if (X->iFlgScaLAPACK == 1 && X->iSolver != SOLVER_SCALAPACK) {
    fprintf(stdoutMPI, cWarnSolverConflict, defname, "ScaLAPACK");
    X->iFlgScaLAPACK = 0;
  }
  if (X->iFlgNGPUSpec == 0) {
    X->iNGPU = (X->iSolver == SOLVER_MAGMA) ? 2 : 0;
  }
#ifndef _SCALAPACK
  if (X->iSolver == SOLVER_SCALAPACK) {
    fprintf(stdoutMPI, cErrSolverBuild, defname, X->iSolver, "ScaLAPACK");
    return -1;
  }
#endif
#ifndef _MAGMA
  if (X->iSolver == SOLVER_MAGMA) {
    fprintf(stdoutMPI, cErrSolverBuild, defname, X->iSolver, "MAGMA");
    return -1;
  }
#endif
#ifndef _ELPA
  if (X->iSolver == SOLVER_ELPA) {
    fprintf(stdoutMPI, cErrSolverBuild, defname, X->iSolver, "ELPA (USE_ELPA=ON)");
    return -1;
  }
#endif
#ifndef _ELPA_GPU
  if (X->iSolver == SOLVER_ELPA && X->iNGPU >= 1) {
    fprintf(stdoutMPI, cErrElpaGPUBuild, defname);
    return -1;
  }
#endif
  return 0;
}
```

検証ブロック（readdef.c:401 の `if(X->iNGPU < 0)` チェックの直前）に解決呼び出しを追加:

```c
  if (ResolveSolver(X, defname) != 0) {
    return (-1);
  }
```

注意: `readdef.c` が `DefCommon.h` を include 済みであることを確認する
（`HubbardNConserved` を使っているので通常は済み。無ければ追加）。

- [ ] **Step 6: テストが通ることを確認**

```bash
cd build && make HPhi -j4 > /dev/null 2>&1 && ctest -R fulldiag_solver_keyword --output-on-failure
```

Expected: PASS

- [ ] **Step 7: 既存 FullDiag 回帰テストを実行**

```bash
cd build && ctest -R fulldiag --output-on-failure
```

Expected: 全 PASS（レガシー入力の挙動が不変であること）

- [ ] **Step 8: コミット**

```bash
git add src/include/DefCommon.h src/include/struct.h src/readdef.c \
        src/include/ErrorMessage.h src/ErrorMessage.c \
        test/fulldiag_solver_keyword.sh test/CMakeLists.txt
git commit -m "Add Solver keyword with legacy-keyword precedence resolution"
```

---

### Task 2: CMake — FindELPA と `USE_ELPA` 配線

**Files:**
- Create: `cmake/FindELPA.cmake`
- Modify: `CMakeLists.txt`（ルート。`find_package(LAPACK)` ブロックの後）
- Modify: `src/CMakeLists.txt`（MAGMA ブロック `if(MAGMA_FOUND)` の直後）
- Create: `config/elpa.cmake`（設定例）

**Interfaces:**
- Produces: CMake 変数 `ELPA_FOUND` / `ELPA_INCLUDE_DIRS` / `ELPA_LIBRARIES`、コンパイル定義 `_ELPA`（検出時）と `_ELPA_GPU`（`elpa_setup_gpu` 存在時）。キャッシュ変数 `ELPA_INCLUDE_DIR` / `ELPA_LIBRARY`（手動指定）、`ELPA_ROOT`（探索ヒント）。

- [ ] **Step 1: FindELPA.cmake を作成**

`cmake/FindELPA.cmake`:

```cmake
# FindELPA.cmake — locate the (non-threaded) ELPA library.
#
# Search order:
#   1. Explicit cache variables ELPA_INCLUDE_DIR / ELPA_LIBRARY
#   2. pkg-config: module "elpa" or versioned "elpa-<version>"
#   3. ELPA_ROOT hint with versioned include dirs (include/elpa-*/elpa/elpa.h)
#
# The threaded variant (elpa_openmp) is NOT supported: HPhi initializes MPI
# with MPI_Init (MPI_THREAD_SINGLE), which is insufficient for it.
#
# Result variables: ELPA_FOUND, ELPA_INCLUDE_DIRS, ELPA_LIBRARIES

set(ELPA_FOUND FALSE)

if(ELPA_INCLUDE_DIR AND ELPA_LIBRARY)
  set(ELPA_INCLUDE_DIRS ${ELPA_INCLUDE_DIR})
  set(ELPA_LIBRARIES ${ELPA_LIBRARY})
  set(ELPA_FOUND TRUE)
endif()

if(NOT ELPA_FOUND)
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_ELPA QUIET elpa)
    if(NOT PC_ELPA_FOUND)
      # Versioned .pc files (elpa-YYYY.MM.SSS.pc). Reject elpa_openmp-*.
      execute_process(
        COMMAND ${PKG_CONFIG_EXECUTABLE} --list-all
        OUTPUT_VARIABLE _elpa_pkg_list ERROR_QUIET)
      string(REGEX MATCH "elpa-[0-9][0-9.]*" _elpa_pkg_name "${_elpa_pkg_list}")
      if(_elpa_pkg_name)
        pkg_check_modules(PC_ELPA QUIET ${_elpa_pkg_name})
      endif()
      string(REGEX MATCH "elpa_openmp-[0-9][0-9.]*" _elpa_omp_name "${_elpa_pkg_list}")
      if(NOT PC_ELPA_FOUND AND _elpa_omp_name)
        message(FATAL_ERROR
          "Only the threaded ELPA variant (${_elpa_omp_name}) was found. "
          "HPhi requires the non-threaded ELPA (MPI_THREAD_SINGLE); "
          "build ELPA without --enable-openmp.")
      endif()
    endif()
    if(PC_ELPA_FOUND)
      set(ELPA_INCLUDE_DIRS ${PC_ELPA_INCLUDE_DIRS})
      set(ELPA_LIBRARIES ${PC_ELPA_LINK_LIBRARIES})
      set(ELPA_FOUND TRUE)
    endif()
  endif()
endif()

if(NOT ELPA_FOUND AND ELPA_ROOT)
  file(GLOB _elpa_inc_candidates "${ELPA_ROOT}/include/elpa-*" "${ELPA_ROOT}/include")
  find_path(ELPA_INCLUDE_DIR_AUTO elpa/elpa.h PATHS ${_elpa_inc_candidates} NO_DEFAULT_PATH)
  find_library(ELPA_LIBRARY_AUTO NAMES elpa PATHS "${ELPA_ROOT}/lib" "${ELPA_ROOT}/lib64" NO_DEFAULT_PATH)
  if(ELPA_INCLUDE_DIR_AUTO AND ELPA_LIBRARY_AUTO)
    set(ELPA_INCLUDE_DIRS ${ELPA_INCLUDE_DIR_AUTO})
    set(ELPA_LIBRARIES ${ELPA_LIBRARY_AUTO})
    set(ELPA_FOUND TRUE)
  endif()
endif()

if(ELPA_FOUND)
  message(STATUS "ELPA include: ${ELPA_INCLUDE_DIRS}")
  message(STATUS "ELPA library: ${ELPA_LIBRARIES}")
  # GPU API (ELPA >= 2023.11.001): detect elpa_setup_gpu
  include(CheckSymbolExists)
  set(CMAKE_REQUIRED_INCLUDES ${ELPA_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${ELPA_LIBRARIES})
  check_symbol_exists(elpa_setup_gpu "elpa/elpa.h" ELPA_HAVE_SETUP_GPU)
  unset(CMAKE_REQUIRED_INCLUDES)
  unset(CMAKE_REQUIRED_LIBRARIES)
else()
  message(FATAL_ERROR
    "USE_ELPA=ON but ELPA was not found. Set ELPA_ROOT=<prefix>, or set "
    "ELPA_INCLUDE_DIR and ELPA_LIBRARY explicitly.")
endif()
```

- [ ] **Step 2: ルート CMakeLists.txt に配線**

`option(USE_SCALAPACK "Use Scalapack" OFF)` の直後に:

```cmake
option(USE_ELPA "Use ELPA for FullDiag" OFF)
```

`find_package(LAPACK)` ブロック（`add_definitions(-D_SCALAPACK)` を含む if/else）の直後に:

```cmake
if(USE_ELPA)
  if(USE_SCALAPACK MATCHES OFF)
    message(STATUS "USE_ELPA=ON requires ScaLAPACK: enabling USE_SCALAPACK.")
    set(USE_SCALAPACK ON CACHE BOOL "Use Scalapack" FORCE)
    add_definitions(-D_SCALAPACK)
  endif()
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
  find_package(ELPA REQUIRED)
endif(USE_ELPA)
```

注意: この if ブロックは既存の `if(USE_SCALAPACK MATCHES OFF) ... -D_lapack`
判定より**後**に置くこと（順序が逆だと `-D_lapack` と `-D_SCALAPACK` が
両方付く）。既存判定より後に置けない場合は `USE_ELPA` の処理を既存判定の
直前に移し、`USE_SCALAPACK` を先に FORCE してから既存判定に流す。

- [ ] **Step 3: src/CMakeLists.txt に配線**

`if(MAGMA_FOUND) ... endif(MAGMA_FOUND)` ブロック（src/CMakeLists.txt:101-104）の直後に:

```cmake
if(ELPA_FOUND)
  include_directories(${ELPA_INCLUDE_DIRS})
  target_link_libraries(HPhi ${ELPA_LIBRARIES})
  add_definitions(-D_ELPA)
  if(ELPA_HAVE_SETUP_GPU)
    add_definitions(-D_ELPA_GPU)
  endif(ELPA_HAVE_SETUP_GPU)
endif(ELPA_FOUND)
```

- [ ] **Step 4: 設定例を作成**

`config/elpa.cmake`:

```cmake
# Example toolchain fragment for building HPhi with ELPA.
#   cmake -DCONFIG=elpa -DUSE_ELPA=ON -DELPA_ROOT=/path/to/elpa/prefix ..
# For GPU execution, ELPA >= 2023.11.001 built with CUDA
# (--enable-nvidia-gpu-kernels) is required; HPhi detects elpa_setup_gpu
# automatically and enables the GPU path (_ELPA_GPU).
set(USE_ELPA ON CACHE BOOL "" FORCE)
set(USE_SCALAPACK ON CACHE BOOL "" FORCE)
# MKL-provided ScaLAPACK (adjust BLACS layer to your MPI):
# set(SCALAPACK_LIBRARIES "-L$ENV{MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64")
```

- [ ] **Step 5: 検証（非 ELPA 環境）**

```bash
cd /Users/k-yoshimi/Dropbox/CLionProjects/HPhi-box/HPhi
rm -rf build_elpa_check && mkdir build_elpa_check && cd build_elpa_check
cmake .. > default.log 2>&1 && echo "default: OK"
cmake -DUSE_ELPA=ON .. > elpa_on.log 2>&1 && echo "UNEXPECTED: should fail" || grep -q "ELPA was not found" elpa_on.log && echo "USE_ELPA=ON fails as expected"
cd .. && rm -rf build_elpa_check
```

Expected: `default: OK` と `USE_ELPA=ON fails as expected` の両方が出る。
その後、標準ビルド確認も実行して既定ビルドが壊れていないことを確認。

- [ ] **Step 6: コミット**

```bash
git add cmake/FindELPA.cmake CMakeLists.txt src/CMakeLists.txt config/elpa.cmake
git commit -m "Add USE_ELPA build option with FindELPA module"
```

---

### Task 3: `matrixlapack_elpa.c` — `diag_elpa_cmp()`

**Files:**
- Create: `src/matrixlapack_elpa.c`
- Create: `src/include/matrixlapack_elpa.h`
- Modify: `src/CMakeLists.txt:53`（ソースリストに `matrixlapack_elpa.c` を追加）

**Interfaces:**
- Consumes: `_ELPA` / `_ELPA_GPU` コンパイル定義（Task 2）
- Produces:
  ```c
  int diag_elpa_cmp(int xNsize, double complex *A_distr,
                    double complex *Z_distr, double *w,
                    int local_nrows, int local_ncols,
                    int myrow, int mycol, int ngpu);
  ```
  戻り値 0=成功 / -1=失敗（全ランクで同一値、内部でエラー同期済み）。
  `A_distr`/`Z_distr` は 2D ブロックサイクリック（nblk=64）のローカル配列
  （カラムメジャー連続、`A_distr` は破壊される）。`w` は長さ `xNsize` の
  固有値配列（全ランクに返る）。

- [ ] **Step 1: ヘッダを作成**

`src/include/matrixlapack_elpa.h`:

```c
/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
#ifndef HPHI_MATRIXLAPACK_ELPA_H
#define HPHI_MATRIXLAPACK_ELPA_H

#ifdef _ELPA
#include <complex.h>

/* ELPA block size used for all 2D block-cyclic descriptors on the ELPA
   path (design doc section 3). */
#define ELPA_NBLK 64

int diag_elpa_cmp(int xNsize, double complex *A_distr,
                  double complex *Z_distr, double *w,
                  int local_nrows, int local_ncols,
                  int myrow, int mycol, int ngpu);
#endif /* _ELPA */

#endif /* HPHI_MATRIXLAPACK_ELPA_H */
```

- [ ] **Step 2: 実装を書く**

`src/matrixlapack_elpa.c`（設計文書 §3「ELPA API 契約」の手順そのまま）:

```c
/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#ifdef _ELPA
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpi.h>
#include <elpa/elpa.h>
#include "matrixlapack_elpa.h"

/* GPU enable option name; single point of change for future amd-gpu etc.
   (design doc section 3, GPU policy) */
static const char *ELPA_GPU_OPTION = "nvidia-gpu";

/* Share a local error across all ranks so every rank takes the same
   branch before each collective ELPA call (design doc section 4). */
static int SyncError(int ierr) {
  int gerr = 0;
  MPI_Allreduce(&ierr, &gerr, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return gerr;
}

/* elpa_set for int values with mandatory error check. */
static int SetElpaInt(elpa_t handle, const char *name, int value) {
  int error = ELPA_OK;
  elpa_set(handle, name, value, &error);
  if (error != ELPA_OK) {
    fprintf(stdout, "  Error: elpa_set(\"%s\", %d) failed: %s\n",
            name, value, elpa_strerr(error));
    return -1;
  }
  return 0;
}

/**
 * @brief Diagonalize a 2D block-cyclic distributed Hermitian matrix with
 * ELPA (design doc section 3). Eigenvalues are returned on all ranks in w;
 * eigenvectors stay distributed in Z_distr. A_distr is destroyed.
 * @return 0 on success, -1 on failure (same value on all ranks).
 */
int diag_elpa_cmp(int xNsize, double complex *A_distr,
                  double complex *Z_distr, double *w,
                  int local_nrows, int local_ncols,
                  int myrow, int mycol, int ngpu) {
  elpa_t handle = NULL;
  int error = ELPA_OK;
  int ierr = 0;

  if (elpa_init(20211125) != ELPA_OK) {
    fprintf(stdout, "  Error: the linked ELPA is older than API 20211125 (2021.11).\n");
    ierr = -1;
  }
  if (SyncError(ierr) != 0) return -1;

  handle = elpa_allocate(&error);
  if (error != ELPA_OK) {
    fprintf(stdout, "  Error: elpa_allocate failed: %s\n", elpa_strerr(error));
    ierr = -1;
  }

  /* Mandatory parameters: set BEFORE elpa_setup (ELPA manual sec. 2). */
  if (ierr == 0) {
    if (SetElpaInt(handle, "na", xNsize) != 0 ||
        SetElpaInt(handle, "nev", xNsize) != 0 ||
        SetElpaInt(handle, "local_nrows", local_nrows) != 0 ||
        SetElpaInt(handle, "local_ncols", local_ncols) != 0 ||
        SetElpaInt(handle, "nblk", ELPA_NBLK) != 0 ||
        SetElpaInt(handle, "mpi_comm_parent",
                   (int)MPI_Comm_c2f(MPI_COMM_WORLD)) != 0 ||
        SetElpaInt(handle, "process_row", myrow) != 0 ||
        SetElpaInt(handle, "process_col", mycol) != 0) {
      ierr = -1;
    }
  }
  if (SyncError(ierr) != 0) goto cleanup_fail;

  error = elpa_setup(handle);
  if (error != ELPA_OK) {
    fprintf(stdout, "  Error: elpa_setup failed: %s\n", elpa_strerr(error));
    ierr = -1;
  }
  if (SyncError(ierr) != 0) goto cleanup_fail;

  /* Tunable runtime options: set AFTER elpa_setup (ELPA manual sec. 2).
     GPU: 1stage is usually faster on GPU; 2stage on CPU. */
  if (ngpu >= 1) {
    if (SetElpaInt(handle, "solver", ELPA_SOLVER_1STAGE) != 0) ierr = -1;
    if (ierr == 0 && SetElpaInt(handle, ELPA_GPU_OPTION, 1) != 0) {
      fprintf(stdout,
              "  Error: the linked ELPA has no NVIDIA GPU support.\n"
              "         Set \"NGPU 0\" in calcmod.def for CPU execution.\n");
      ierr = -1;
    }
  } else {
    if (SetElpaInt(handle, "solver", ELPA_SOLVER_2STAGE) != 0) ierr = -1;
  }
  if (SyncError(ierr) != 0) goto cleanup_fail;

#ifdef _ELPA_GPU
  if (ngpu >= 1) {
    error = elpa_setup_gpu(handle);
    if (error != ELPA_OK) {
      fprintf(stdout,
              "  Error: elpa_setup_gpu failed: %s\n"
              "         Check the GPU environment, or set \"NGPU 0\" for CPU execution.\n",
              elpa_strerr(error));
      ierr = -1;
    }
    if (SyncError(ierr) != 0) goto cleanup_fail;
  }
#else
  if (ngpu >= 1) {
    /* readdef.c rejects this combination at startup; defense in depth. */
    fprintf(stdout, "  Error: this HPhi build has no ELPA GPU API (_ELPA_GPU).\n");
    ierr = -1;
    if (SyncError(ierr) != 0) goto cleanup_fail;
  }
#endif

  /* Type-generic macro resolves to the double-complex solver. */
  elpa_eigenvectors(handle, A_distr, w, Z_distr, &error);
  if (error != ELPA_OK) {
    fprintf(stdout, "  Error: elpa_eigenvectors failed: %s\n", elpa_strerr(error));
    ierr = -1;
  }
  if (SyncError(ierr) != 0) goto cleanup_fail;

  elpa_deallocate(handle, &error);
  elpa_uninit(&error);
  return 0;

cleanup_fail:
  if (handle != NULL) {
    elpa_deallocate(handle, &error);
  }
  elpa_uninit(&error);
  return -1;
}
#endif /* _ELPA */
```

- [ ] **Step 3: ビルドに追加**

`src/CMakeLists.txt:53` のソースリストで `matrixlapack_magma.c` の直後に
`matrixlapack_elpa.c` を追加:

```cmake
  matrixlapack.c matrixlapack_magma.c matrixlapack_elpa.c matrixscalapack.c
```

- [ ] **Step 4: 非 ELPA ビルドが壊れていないことを確認**

標準ビルド確認を実行（`_ELPA` 未定義ではファイル全体が空になるので通る）:

Expected: `[100%] Built target HPhi`

- [ ] **Step 5: コミット**

```bash
git add src/matrixlapack_elpa.c src/include/matrixlapack_elpa.h src/CMakeLists.txt
git commit -m "Add ELPA diagonalization module diag_elpa_cmp"
```

---

### Task 4: `GetEigenVectorBlock()` — rank 0 へのブロック回収

**Files:**
- Modify: `src/matrixscalapack.c`（末尾、`diag_scalapack_cmp` の後）
- Modify: `src/include/matrixscalapack.h`（宣言追加）

**Interfaces:**
- Consumes: `Z_vec`（2D ブロックサイクリック分散固有ベクトル行列）と
  `descZ_vec`（既存グローバル、`descZ[1]` = 2D コンテキスト）
- Produces:
  ```c
  int GetEigenVectorBlock(long int idx, long int xNsize,
                          double complex *Z, int *descZ,
                          double complex *vec);
  void FreeEigenVectorGatherContext(void);
  ```
  `idx` は 0 始まりの固有状態番号。rank 0 の `vec[0..xNsize-1]` に固有
  ベクトルが入る（他ランクの `vec` は作業バッファで内容不定、非 NULL 必須）。
  1×1 宛先グリッドは初回呼び出しで生成・キャッシュし、
  `FreeEigenVectorGatherContext()` で解放（phys ループ後に呼ぶ）。
  戻り値 0=成功。

- [ ] **Step 1: 実装を書く**

`src/matrixscalapack.c` 末尾（`#endif` の直前）に追加:

```c
/* Cached destination grid for GetEigenVectorBlock: a 1x1 BLACS grid
   containing rank 0 only, built once per run (design doc section 3). */
static int ictxt_gather = -100;   /* -100: not initialized */
static int desc_gather[9];

/**
 * @brief Initialize the rank-0-only destination grid and descriptor.
 * All ranks must call this (blacs_gridmap_ is collective). Follows the
 * p?gemr2d contract: ranks outside the destination grid keep
 * desc[CTXT_] = -1 while all other descriptor fields stay valid.
 */
static void InitEigenVectorGatherContext(long int xNsize) {
  int i_negone = -1, i_zero = 0;
  int imap[1] = {0};
  int ld = 1, np_gather = 1;
  int myrow_g, mycol_g, nprow_g, npcol_g;

  blacs_get_(&i_negone, &i_zero, &ictxt_gather);
  blacs_gridmap_(&ictxt_gather, imap, &ld, &np_gather, &np_gather);

  /* Fully initialize the descriptor on ALL ranks (some implementations
     inspect fields other than CTXT_), then mark non-participants. */
  desc_gather[0] = 1;               /* DTYPE_: dense */
  desc_gather[1] = ictxt_gather;    /* CTXT_ */
  desc_gather[2] = (int)xNsize;     /* M_ */
  desc_gather[3] = 1;               /* N_ */
  desc_gather[4] = (int)xNsize;     /* MB_ */
  desc_gather[5] = 1;               /* NB_ */
  desc_gather[6] = 0;               /* RSRC_ */
  desc_gather[7] = 0;               /* CSRC_ */
  desc_gather[8] = (int)xNsize;     /* LLD_ */

  blacs_gridinfo_(&ictxt_gather, &nprow_g, &npcol_g, &myrow_g, &mycol_g);
  if (myrow_g < 0) {
    desc_gather[1] = -1;            /* not in the destination grid */
  }
}

/**
 * @brief Gather one distributed eigenvector (column idx of Z) into
 * vec[0..xNsize-1] on rank 0 with a single pzgemr2d_ block transfer,
 * replacing the element-wise pzelget_ loop of GetEigenVector.
 * @return 0 on success.
 */
int GetEigenVectorBlock(long int idx, long int xNsize,
                        double complex *Z, int *descZ,
                        double complex *vec) {
  const long int i_one = 1;
  long int icol = idx + 1;
  long int m = xNsize, n = 1;

  if (ictxt_gather == -100) {
    InitEigenVectorGatherContext(xNsize);
  }
  /* Last argument: a context containing the union of both grids
     = the all-rank 2D context of Z. */
  pzgemr2d_(&m, &n, Z, (long int *)&i_one, &icol, descZ,
            vec, (long int *)&i_one, (long int *)&i_one, desc_gather,
            &descZ[1]);
  return 0;
}

/**
 * @brief Release the cached destination grid (call after the phys loop).
 */
void FreeEigenVectorGatherContext(void) {
  int myrow_g, mycol_g, nprow_g, npcol_g;
  if (ictxt_gather == -100) return;
  blacs_gridinfo_(&ictxt_gather, &nprow_g, &npcol_g, &myrow_g, &mycol_g);
  if (myrow_g >= 0) {
    blacs_gridexit_(&ictxt_gather);
  }
  ictxt_gather = -100;
}
```

注意（実装時に必ず確認）:
- `matrixscalapack.h` にある既存の BLACS/ScaLAPACK 外部宣言の**引数型
  （`int` か `long int` か）に合わせる**こと。上のコードは既存
  `diag_scalapack_cmp` が `long int` と `int` を混用している慣例に従うが、
  ヘッダの `pzgemr2d_` / `blacs_gridmap_` プロトタイプが無ければ追加し、
  既存宣言と同じ整数幅にそろえる（不一致は実行時クラッシュの元）。
- `descZ[1]`（CTXT_）をコンテキスト引数に渡す。

`src/include/matrixscalapack.h` の `GetEigenVector` 宣言の近くに追加:

```c
int GetEigenVectorBlock(long int idx, long int xNsize,
                        double complex *Z, int *descZ,
                        double complex *vec);
void FreeEigenVectorGatherContext(void);
void pzgemr2d_(long int *m, long int *n,
               double complex *A, long int *ia, long int *ja, int *desca,
               double complex *B, long int *ib, long int *jb, int *descb,
               int *ictxt);
void blacs_gridmap_(int *ictxt, int *usermap, int *ldumap,
                    int *nprow, int *npcol);
void blacs_gridexit_(int *ictxt);
```

（既に同名宣言があれば重複させない。整数幅は既存宣言に合わせる。）

- [ ] **Step 2: ビルド確認（_SCALAPACK 無しの既定ビルド）**

標準ビルド確認を実行。`matrixscalapack.c` は `#ifdef _SCALAPACK` ガード内
なので既定ビルドに影響しないこと。

Expected: `[100%] Built target HPhi`

- [ ] **Step 3: コミット**

```bash
git add src/matrixscalapack.c src/include/matrixscalapack.h
git commit -m "Add block-transfer eigenvector gather GetEigenVectorBlock"
```

---

### Task 5: `lapack_diag.c` の Solver 分岐と ELPA 経路、`phys.c` の回収切替

**Files:**
- Modify: `src/lapack_diag.c:49-92`
- Modify: `src/phys.c:96-132`（回収関数の切替）と 210 付近（コンテキスト解放）

**Interfaces:**
- Consumes: `X->Def.iSolver`（Task 1）、`diag_elpa_cmp`（Task 3）、
  `GetEigenVectorBlock` / `FreeEigenVectorGatherContext`（Task 4）、
  グローバル `Ham` / `v0` / `Z_vec` / `descZ_vec` / `use_scalapack`
- Produces: `Solver 3` での FullDiag 実行経路。固有値は `v0[0..N-1]`
  （全ランク）、固有ベクトルは `Z_vec` 分散のまま、`use_scalapack=1`。

- [ ] **Step 1: lapack_diag.c を iSolver 分岐に書き換え**

`src/lapack_diag.c` の `if (X->Def.iNGPU == 0) { ... } else { ... }`
（lapack_diag.c:56-92）を次に置き換える（`#include "DefCommon.h"` と
`#ifdef _ELPA` 用 include を先頭に追加。既存の SCALAPACK 用ローカル変数
宣言ブロックはそのまま使う）:

```c
  switch (X->Def.iSolver) {
  case SOLVER_SCALAPACK:
#ifdef _SCALAPACK
    if (nproc > 1) {
      fprintf(stdoutMPI, "Using SCALAPACK\n\n");
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      MPI_Dims_create(size, 2, dims);
      nprow = dims[0]; npcol = dims[1];

      blacs_pinfo_(&iam, &nprocs);
      blacs_get_(&i_negone, &i_zero, &ictxt);
      blacs_gridinit_(&ictxt, "R", &nprow, &npcol);
      blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);

      mb = GetBlockSize(xMsize, size);
      mp = numroc_(&xMsize, &mb, &myrow, &i_zero, &nprow);
      nq = numroc_(&xMsize, &mb, &mycol, &i_zero, &npcol);
      Z_vec = malloc(mp * nq * sizeof(complex double));
      diag_scalapack_cmp(xMsize, Ham, v0, Z_vec, descZ_vec);
    } else {
      ZHEEVall(xMsize, Ham, v0, L_vec);
    }
#endif
    break;

  case SOLVER_MAGMA:
#ifdef _MAGMA
    if (myrank == 0) {
      if (diag_magma_cmp(xMsize, Ham, v0, L_vec, X->Def.iNGPU) != 0) {
        return -1;
      }
    }
#endif
    break;

  case SOLVER_ELPA:
#ifdef _ELPA
    if (lapack_diag_elpa(X, xMsize) != 0) {
      return -1;
    }
#endif
    break;

  default: /* SOLVER_LAPACK */
    ZHEEVall(xMsize, Ham, v0, L_vec);
    break;
  }
```

注意: `ResolveSolver`（Task 1）がビルドフラグと不整合な `iSolver` を
起動時に弾いているので、各 case の `#ifdef` 外に到達することはない
（ガードは空分岐になるだけで安全）。既存の SCALAPACK/MAGMA case の中身は
**現行コードの字句をそのまま移す**こと（挙動維持）。

- [ ] **Step 2: ELPA 経路のヘルパを同ファイルに実装**

`lapack_diag()` の直前に追加（`#ifdef _ELPA` ガード内）:

```c
#ifdef _ELPA
#include "matrixlapack_elpa.h"

/**
 * @brief FullDiag via ELPA (phase 1: fill the 2D block-cyclic matrix from
 * the replicated Ham with pzelset, then call diag_elpa_cmp).
 * Eigenvalues land in v0 on all ranks; eigenvectors stay in Z_vec.
 * NOTE (phase 1): the replicated Ham plus A_distr plus Z_vec coexist in
 * memory, so verification runs must stay at small N (design doc sec. 3).
 */
static int lapack_diag_elpa(struct BindStruct *X, long int xMsize) {
  int i_negone = -1, i_zero_i = 0;
  const long int i_zero = 0;
  int rank, size;
  int nprow, npcol, myrow, mycol;
  int ictxt;
  long int mb = ELPA_NBLK, mp, nq, i, j;
  int lld, dims[2] = {0, 0};
  int iam, nprocs, info;
  double complex *A_distr;
  double *w;
  int descA[9];
  int ierr;

  fprintf(stdoutMPI, "Using ELPA (%s)\n\n",
          X->Def.iNGPU >= 1 ? "GPU" : "CPU");

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size, 2, dims);
  nprow = dims[0]; npcol = dims[1];

#ifdef _ELPA_GPU
  /* Startup consistency warning (design doc sec. 2): ranks per node
     should be a multiple of NGPU (ideally equal: 1 rank per GPU). */
  if (X->Def.iNGPU >= 1) {
    MPI_Comm comm_node;
    int nrank_node;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                        MPI_INFO_NULL, &comm_node);
    MPI_Comm_size(comm_node, &nrank_node);
    MPI_Comm_free(&comm_node);
    if (nrank_node % X->Def.iNGPU != 0) {
      fprintf(stdoutMPI,
              "Warning: ranks per node (%d) is not a multiple of NGPU (%d):\n"
              "         GPUs may idle or be shared unevenly. Recommended: 1 rank per GPU.\n",
              nrank_node, X->Def.iNGPU);
    }
  }
#endif

  blacs_pinfo_(&iam, &nprocs);
  blacs_get_(&i_negone, &i_zero_i, &ictxt);
  blacs_gridinit_(&ictxt, "R", &nprow, &npcol);
  blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);

  mp = numroc_(&xMsize, &mb, &myrow, &i_zero, &nprow);
  nq = numroc_(&xMsize, &mb, &mycol, &i_zero, &npcol);
  lld = (mp > 0) ? mp : 1;

  descinit_(descA, &xMsize, &xMsize, &mb, &mb, &i_zero, &i_zero, &ictxt, &lld, &info);
  descinit_(descZ_vec, &xMsize, &xMsize, &mb, &mb, &i_zero, &i_zero, &ictxt, &lld, &info);

  A_distr = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));
  Z_vec = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));
  w = malloc(xMsize * sizeof(double));

  for (i = 0; i < xMsize; i++) {
    for (j = 0; j < xMsize; j++) {
      DivMat(i, j, Ham[i][j], A_distr, descA);
    }
  }

  ierr = diag_elpa_cmp((int)xMsize, A_distr, Z_vec, w,
                       (int)mp, (int)nq, (int)myrow, (int)mycol,
                       X->Def.iNGPU);
  free(A_distr);
  if (ierr != 0) {
    free(w);
    return -1;
  }

  for (i = 0; i < xMsize; i++) {
    v0[i] = w[i];
  }
  free(w);
  use_scalapack = 1;
  return 0;
}
#endif /* _ELPA */
```

注意（実装時に必ず確認）:
- `DivMat` / `descinit_` / `numroc_` / BLACS 各関数のプロトタイプ整数幅は
  `matrixscalapack.h` の既存宣言に合わせる（`diag_scalapack_cmp` が
  同じ呼び方をしているのでそれを写す）。
- `descZ_vec` はグローバル（`global.h:76-77`）。ELPA 経路では
  `GetBlockSize` でなく `ELPA_NBLK` を使うので、`descZ_vec` もここで
  nblk=64 により初期化し直している。

- [ ] **Step 3: phys.c の回収切替**

`src/phys.c` の `GetEigenVector(i, i_max, Z_vec, descZ_vec, vec_tmp);`
（phys.c:103）を次に置き換える:

```c
#ifdef _ELPA
      if (X->Def.iSolver == SOLVER_ELPA) {
        GetEigenVectorBlock(i, i_max, Z_vec, descZ_vec, vec_tmp);
      } else {
        GetEigenVector(i, i_max, Z_vec, descZ_vec, vec_tmp);
      }
#else
      GetEigenVector(i, i_max, Z_vec, descZ_vec, vec_tmp);
#endif
```

`phys()` 末尾の `if(use_scalapack) free(vec_tmp);`（phys.c:212）の直後に追加:

```c
#ifdef _ELPA
  if (use_scalapack && X->Def.iSolver == SOLVER_ELPA) {
    FreeEigenVectorGatherContext();
  }
#endif
```

`phys.c` の include に `DefCommon.h` が無ければ追加する。

- [ ] **Step 4: ビルドと回帰テスト**

標準ビルド確認＋既存 FullDiag 回帰:

```bash
cd build && make HPhi -j4 2>&1 | tail -3 && ctest -R fulldiag --output-on-failure
```

Expected: ビルド成功、全 FullDiag テスト PASS（既定ビルドでは ELPA 経路は
コンパイルされないが、lapack_diag.c の分岐書き換えの回帰を検出する）

- [ ] **Step 5: ELPA 統合テストを追加**

`test/fulldiag_elpa_hubbard_chain.sh` を新規作成（`chmod +x`）:

```sh
#!/bin/sh -e

mkdir -p fulldiag_elpa_hubbard_chain/
cd fulldiag_elpa_hubbard_chain

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

../../src/HPhi -sdry stan.in
echo "Solver  3" >> calcmod.def
echo "NGPU    0" >> calcmod.def
${MPIRUNFC} ../../src/HPhi -e namelist.def

# エネルギーと二重占有率のみ比較する（ELPA 経路は use_scalapack 扱いで
# S2/Sz を計算しないため、既存 fulldiag_hubbard_chain の参照から
# 該当列だけを使う）
cat > reference_ed.dat <<EOF
  -2.102748   0.287325
  -1.806424   0.335409
  -1.068140   0.277708
  -0.828427   0.146447
  -0.828427   0.146447
   0.000000   0.000000
   0.581449   1.079437
EOF
awk 'NR>1 && NR<=8 {printf "%11.6f %10.6f\n", $1, $5}' output/zvo_phys_Nup2_Ndown2.dat > ed.dat
paste ed.dat reference_ed.dat > paste_ed.dat
diff=`awk 'BEGIN{max=0}{d=$1-$3; if(d<0)d=-d; if(d>max)max=d; d=$2-$4; if(d<0)d=-d; if(d>max)max=d}END{print max}' paste_ed.dat`
test "`echo "$diff < 0.000001" | bc`" = "1"

echo "fulldiag_elpa_hubbard_chain: OK"
```

`test/CMakeLists.txt` の末尾付近に登録（ELPA ビルド時のみ）:

```cmake
if(USE_ELPA)
  add_hphi_mpi_test(fulldiag_elpa_hubbard_chain min:2)
endif(USE_ELPA)
```

注意: `zvo_phys_*.dat` の列構成（$1=`<H>`, $2=`<N>`, $3=`<Sz>`, $4=`<S2>`,
$5=`<D>`）は `fulldiag_hubbard_chain.sh` の reference と同一。実装時に
実出力のヘッダ行で列を確認し、awk の列番号を合わせること。

- [ ] **Step 6: 非 ELPA 環境での確認**

```bash
cd build && cmake .. > /dev/null && ctest -N | grep elpa
```

Expected: 出力なし（`USE_ELPA=OFF` では ELPA テストが登録されない）

- [ ] **Step 7: コミット**

```bash
git add src/lapack_diag.c src/phys.c test/fulldiag_elpa_hubbard_chain.sh test/CMakeLists.txt
git commit -m "Wire Solver 3 (ELPA) path into FullDiag with block eigenvector gather"
```

---

### Task 6: 固有対の数学的検証プログラム（ELPA ビルド時のみ）

**Files:**
- Create: `test/unit/elpa_eigen_check.c`
- Modify: `test/CMakeLists.txt`（実行ターゲットとテスト登録）

**Interfaces:**
- Consumes: `diag_elpa_cmp`（Task 3 のシグネチャ）、
  `GetEigenVectorBlock` / `FreeEigenVectorGatherContext`（Task 4）、
  `DivMat`（`matrixscalapack.c` 既存）
- Produces: 実行ファイル `elpa_eigen_check`（引数なし、MPI 実行、
  成功で exit 0 / 失敗で exit 1）

- [ ] **Step 1: 検証プログラムを書く**

`test/unit/elpa_eigen_check.c`（N=97 の決定的エルミート行列で、設計文書 §6
のスケール則閾値により残差・直交性・LAPACK 固有値一致を検査する）:

```c
/* Mathematical verification of the ELPA FullDiag path:
   residual ||A z - w z||, orthogonality ||Z^H Z - I||, and eigenvalue
   agreement with LAPACK zheev, with thresholds c*N*eps*||A|| and c*N*eps
   (c = 50) per the design doc section 6. Run with any rank count. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <mpi.h>
#include "matrixscalapack.h"
#include "matrixlapack_elpa.h"

#define NDIM 97
#define TOL_C 50.0

extern void zheev_(char *jobz, char *uplo, int *n, double complex *a,
                   int *lda, double *w, double complex *work, int *lwork,
                   double *rwork, int *info);

static double complex MatElem(int i, int j) {
  double re = 1.0 / (1.0 + fabs((double)(i - j)));
  double im = (double)(i - j) / (double)(NDIM * NDIM);
  return re + im * I;   /* A[j][i] = conj(A[i][j]) by construction */
}

int main(int argc, char **argv) {
  int i_negone = -1, i_zero_i = 0;
  const long int i_zero = 0;
  int rank, size, ictxt, iam, nprocs, info, ok = 1;
  int nprow, npcol, myrow, mycol;
  long int n = NDIM, mb = ELPA_NBLK, mp, nq, i, j, k;
  int lld, dims[2] = {0, 0};
  int descA[9], descZ[9];
  double complex *A_distr, *Z_distr, *vecs, *vec_tmp;
  double *w, anorm = 0.0, eps = DBL_EPSILON;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size, 2, dims);
  nprow = dims[0]; npcol = dims[1];
  blacs_pinfo_(&iam, &nprocs);
  blacs_get_(&i_negone, &i_zero_i, &ictxt);
  blacs_gridinit_(&ictxt, "R", &nprow, &npcol);
  blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);

  mp = numroc_(&n, &mb, &myrow, &i_zero, &nprow);
  nq = numroc_(&n, &mb, &mycol, &i_zero, &npcol);
  lld = (mp > 0) ? mp : 1;
  descinit_(descA, &n, &n, &mb, &mb, &i_zero, &i_zero, &ictxt, &lld, &info);
  descinit_(descZ, &n, &n, &mb, &mb, &i_zero, &i_zero, &ictxt, &lld, &info);

  A_distr = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));
  Z_distr = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));
  w = malloc(n * sizeof(double));
  vec_tmp = malloc(n * sizeof(double complex));
  vecs = malloc(n * n * sizeof(double complex)); /* rank 0: gathered Z */

  for (i = 0; i < n; i++) {
    double colsum = 0.0;
    for (j = 0; j < n; j++) {
      DivMat(i, j, MatElem(i, j), A_distr, descA);
      colsum += cabs(MatElem(i, j));
    }
    if (colsum > anorm) anorm = colsum;  /* ||A||_1 */
  }

  if (diag_elpa_cmp((int)n, A_distr, Z_distr, w,
                    (int)mp, (int)nq, (int)myrow, (int)mycol, 0) != 0) {
    if (rank == 0) fprintf(stderr, "diag_elpa_cmp failed\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  for (k = 0; k < n; k++) {
    GetEigenVectorBlock(k, n, Z_distr, descZ, vec_tmp);
    if (rank == 0) for (i = 0; i < n; i++) vecs[k * n + i] = vec_tmp[i];
  }
  FreeEigenVectorGatherContext();

  if (rank == 0) {
    /* (a) eigenvalues vs LAPACK zheev */
    double complex *a_full = malloc(n * n * sizeof(double complex));
    double *w_ref = malloc(n * sizeof(double));
    double *rwork = malloc((3 * n - 2) * sizeof(double));
    int lwork = 4 * NDIM, n_int = NDIM;
    double complex *work = malloc(lwork * sizeof(double complex));
    for (j = 0; j < n; j++) for (i = 0; i < n; i++)
      a_full[j * n + i] = MatElem(i, j);
    zheev_("N", "U", &n_int, a_full, &n_int, w_ref, work, &lwork, rwork, &info);
    for (k = 0; k < n; k++) {
      if (fabs(w[k] - w_ref[k]) > TOL_C * n * eps * anorm) {
        fprintf(stderr, "eigenvalue %ld mismatch: %e vs %e\n", k, w[k], w_ref[k]);
        ok = 0;
      }
    }
    /* (b) residual ||A z_k - w_k z_k||_inf */
    for (k = 0; k < n && ok; k++) {
      for (i = 0; i < n; i++) {
        double complex r = -w[k] * vecs[k * n + i];
        for (j = 0; j < n; j++) r += MatElem(i, j) * vecs[k * n + j];
        if (cabs(r) > TOL_C * n * eps * anorm) {
          fprintf(stderr, "residual too large: state %ld row %ld: %e\n", k, i, cabs(r));
          ok = 0; break;
        }
      }
    }
    /* (c) orthogonality |z_k^H z_l - delta_kl| */
    for (k = 0; k < n && ok; k++) {
      for (j = k; j < n; j++) {
        double complex dot = 0.0;
        for (i = 0; i < n; i++) dot += conj(vecs[k * n + i]) * vecs[j * n + i];
        if (cabs(dot - (k == j ? 1.0 : 0.0)) > TOL_C * n * eps) {
          fprintf(stderr, "orthogonality violated: (%ld,%ld) %e\n", k, j, cabs(dot));
          ok = 0; break;
        }
      }
    }
    printf("elpa_eigen_check: %s\n", ok ? "OK" : "FAILED");
    free(a_full); free(w_ref); free(rwork); free(work);
  }
  MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  free(A_distr); free(Z_distr); free(w); free(vec_tmp); free(vecs);
  MPI_Finalize();
  return ok ? 0 : 1;
}
```

- [ ] **Step 2: ビルドとテストを登録**

`test/CMakeLists.txt` 末尾（Task 5 で追加した `if(USE_ELPA)` ブロック内）を
拡張:

```cmake
if(USE_ELPA)
  add_hphi_mpi_test(fulldiag_elpa_hubbard_chain min:2)

  add_executable(elpa_eigen_check unit/elpa_eigen_check.c
    ${CMAKE_SOURCE_DIR}/src/matrixscalapack.c
    ${CMAKE_SOURCE_DIR}/src/matrixlapack_elpa.c)
  target_include_directories(elpa_eigen_check PRIVATE
    ${CMAKE_SOURCE_DIR}/src/include ${ELPA_INCLUDE_DIRS})
  target_compile_definitions(elpa_eigen_check PRIVATE -D_ELPA -D_SCALAPACK -DMPI)
  target_link_libraries(elpa_eigen_check ${ELPA_LIBRARIES}
    ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_C_LIBRARIES} m)
  add_test(NAME elpa_eigen_check
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3
            $<TARGET_FILE:elpa_eigen_check>)
endif(USE_ELPA)
```

注意: `matrixscalapack.c` が `elpa_eigen_check` 単体でリンクできない依存
（HPhi グローバル等）を持っていた場合は、その関数群（`DivMat`/`GetEigenVectorBlock`/
`GetBlockSize` 系）だけを使う形なので通常は自己完結。リンクエラーが出たら
欠けた依存シンボルを確認し、必要最小のスタブを `unit/` に置くこと。

- [ ] **Step 3: 非 ELPA 環境での確認とコミット**

標準ビルド確認（既定ビルドに影響なし）ののち:

```bash
git add test/unit/elpa_eigen_check.c test/CMakeLists.txt
git commit -m "Add ELPA eigenpair verification test program"
```

---

### Task 7: ドキュメント更新

**Files:**
- Modify: `doc/ja/source/filespecification/expertmode_ja/CalcMod_file_ja.rst:205-231`
- Modify: `doc/en/source/filespecification/expertmode_en/CalcMod_file_en.rst:220-245`
- Create: `test/manual/elpa_gpu_check.md`（手動 GPU 検証プロトコル）

**Interfaces:**
- Consumes: `Solver`/`NGPU` の確定仕様（設計文書 §2、Task 1 の実装）
- Produces: ユーザー向け仕様記述（フェーズ1スコープ: `Solver` キーワード、
  `NGPU` の ELPA 時の意味、`ScaLAPACK` 非推奨注記、手動検証手順）

- [ ] **Step 1: 日本語 CalcMod 仕様に Solver を追加**

`CalcMod_file_ja.rst` の `Scalapack` 項の直前に追加し、`Scalapack` 項に
非推奨注記、`NGPU` 項に ELPA の意味を追記:

```rst
-  ``Solver``

   **形式 :** int型 (デフォルト値: 旧キーワードから自動決定)

   | **説明 :** (FullDiag)
     全対角化計算の対角化バックエンドを指定します。
   | 0: LAPACK (逐次)
   | 1: ScaLAPACK
   | 2: MAGMA (シングルノード・マルチGPU)
   | 3: ELPA (マルチノード対応、CPU/GPU)
   | 未指定の場合は従来の ``Scalapack``/``NGPU`` キーワードから
     従来どおりの動作になるよう自動決定されます。
     ``Solver 3`` は ELPA を有効にしたビルド (``USE_ELPA=ON``) が必要です。
```

`Scalapack` 項の説明末尾に追加:

```rst
   | (非推奨) 本キーワードは後方互換のために残されています。
     今後は ``Solver 1`` を使用してください。
```

`NGPU` 項の説明を次で置き換え:

```rst
   **説明 :** (FullDiag)
   ノードあたりの使用 GPU 枚数を指定します。
   ``Solver 2`` (MAGMA) では単一プロセスから使う GPU 枚数です。
   ``Solver 3`` (ELPA) では 0 で CPU 実行、1 以上で GPU 実行となり、
   プロセスと GPU の対応は ELPA が自動割当します (1 プロセス 1 GPU)。
   ノードあたりの MPI プロセス数を GPU 枚数に合わせる実行を推奨します。
   なお ``NGPU`` は使用 GPU 枚数を物理的には制限しません。枚数を厳密に
   制限する場合はジョブスケジューラ側 (``CUDA_VISIBLE_DEVICES`` 等) で
   行ってください。GPU 実行には ELPA 2023.11.001 以降が必要です。
```

- [ ] **Step 2: 英語版に同内容を追加**

`CalcMod_file_en.rst` の対応箇所（`NGPU` 項 232 行付近）に Step 1 と同内容の
英語版を追加する（既存英語項の文体に合わせる）:

```rst
*  ``Solver``

   **Type :** int (default: resolved from legacy keywords)

   **Description :** (FullDiag)
   Diagonalization backend for the full diagonalization method:
   0 (LAPACK, serial), 1 (ScaLAPACK), 2 (MAGMA, single-node multi-GPU),
   3 (ELPA, multi-node CPU/GPU; requires a build with ``USE_ELPA=ON``).
   When omitted, the backend is resolved from the legacy ``Scalapack``
   and ``NGPU`` keywords so that existing inputs behave as before.
```

および `NGPU`/`Scalapack` 項へ Step 1 と対応する英語の追記。

- [ ] **Step 3: 手動 GPU 検証プロトコルを作成**

`test/manual/elpa_gpu_check.md`:

```markdown
# ELPA GPU manual verification protocol (no GPU CI)

Run before releases that touch the FullDiag/ELPA path.
Target: clavius (single node multi-GPU) and one multi-node GPU system.

## Build
    cmake -DUSE_ELPA=ON -DELPA_ROOT=<prefix> ..   # ELPA >= 2023.11.001 (CUDA build)
    make HPhi elpa_eigen_check

## Checklist
1. `ctest -R elpa_eigen_check` ... expect PASS (CPU path sanity)
2. `ctest -R fulldiag_elpa_hubbard_chain` ... expect PASS
3. GPU smoke (1 rank / 1 GPU):
   `Solver 3`, `NGPU 1`, 12-site Hubbard chain FullDiag;
   compare zvo_phys energies against a `Solver 0` run (tol 1e-8);
   confirm "Using ELPA (GPU)" in stdout and GPU utilization in nvidia-smi.
4. GPU multi-rank (ranks = GPUs per node): same comparison at 14 sites.
5. Failure-path check: run with `NGPU 1` against a CPU-only ELPA build;
   expect a clear error mentioning NGPU 0 fallback instruction, no silent
   CPU execution.
Record results (date, host, ELPA version, commit) at the bottom of this file.
```

- [ ] **Step 4: ドキュメントのビルド確認（任意）とコミット**

rst はビルド必須ではない（Sphinx 環境があれば `make -C doc/ja html` 相当で
警告確認）。

```bash
git add doc/ja/source/filespecification/expertmode_ja/CalcMod_file_ja.rst \
        doc/en/source/filespecification/expertmode_en/CalcMod_file_en.rst \
        test/manual/elpa_gpu_check.md
git commit -m "Document Solver keyword, ELPA NGPU semantics, and manual GPU protocol"
```

---

### Task 8: ELPA 実機検証（clavius、チェックリスト）

**Files:** なし（実行のみ。結果は `test/manual/elpa_gpu_check.md` 末尾に記録して追記コミット）

**Interfaces:**
- Consumes: Task 1–7 の全成果物

このタスクはローカルでは実行できない。clavius に ELPA（CPU 版で開始、
CUDA 版は別途）を導入した上で以下を実施する:

- [ ] **Step 1: clavius へ転送し ELPA CPU 版でビルド**（`USE_ELPA=ON`）
- [ ] **Step 2: `ctest -R "elpa"` を ranks 1/2/4 相当で実行し PASS を確認**
- [ ] **Step 3: `fulldiag_hubbard_chain` の参照値と `Solver 3` の全固有値一致を確認**（tol 1e-8）
- [ ] **Step 4: CUDA 版 ELPA（>= 2023.11.001）でビルドし、`test/manual/elpa_gpu_check.md` のチェックリスト 3–5 を実施**
- [ ] **Step 5: 結果を `elpa_gpu_check.md` に記録してコミット**

GPU ジョブ実行時は `~/.claude/CLAUDE.md` の GPU リソース管理ルール
（`nvidia-smi` での事前確認、`CUDA_VISIBLE_DEVICES` の明示、1 ランク
1 GPU）に従うこと。

---

## 完了条件（フェーズ1）

- `Solver 0/1/2` の全既存テストが無変更で PASS
- 非 ELPA ビルド（既定・CI）が従来どおりビルド・テスト可能
- ELPA ビルドで `elpa_eigen_check`・`fulldiag_elpa_hubbard_chain`・
  `fulldiag_solver_keyword` が PASS
- clavius での実機チェックリスト完了（CPU 必須、GPU は ELPA CUDA 版導入後）
- ベンチマークゲート: clavius で `Solver 1` 比の対角化時間を記録
  （設計文書 §8、記録先は `test/manual/elpa_gpu_check.md`）
