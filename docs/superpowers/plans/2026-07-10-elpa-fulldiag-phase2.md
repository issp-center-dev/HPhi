# ELPA FullDiag フェーズ2（ハミルトニアン分散生成）実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `Solver 3`（ELPA）かつ `nproc > 1` のとき、全複製 `Ham`（16N²バイト/ランク）を廃し、各ランクが担当列だけの1D列パネルを生成して `pzgemr2d` で2Dブロックサイクリックへ再分散する — 行列格納を O(N²/P) に落とし N~10⁵ を可能にする。

**Architecture:** 設計文書 `docs/superpowers/specs/2026-07-10-elpa-fulldiag-design.md` §3「フェーズ2」（承認済みv4）。ハミルトニアン書き込みを格納マクロ（`hamstore.h`）に一本化し、分散モードでは連続カラムメジャーの1Dパネル `Ham_local` に格納、生成ループの列範囲を担当分に制限、`RedistPanelToBlockCyclic()`（`pzgemr2d` 1回）でELPAの2D分散行列へ変換する。既存の複製経路（`Solver 0/1/2` と `Solver 3`＋`nproc==1`）は実行時分岐で完全温存。

**Tech Stack:** C99, MPI, BLACS/ScaLAPACK（`pzgemr2d_`）, ELPA（フェーズ1で導入済み）, CMake/ctest。

## Global Constraints

- 設計文書 §3「フェーズ 2」が正。判断に迷ったらそこへ戻る。
- 分散パネル有効条件は **`X->Def.iSolver == SOLVER_ELPA && nproc > 1` のみ**。それ以外の経路（`Solver 0/1/2`、`Solver 3`+`nproc==1`）は挙動・出力とも一切不変。
- 列所有は1Dブロック（連続範囲）: `NC = ceil(N/P)`、rank p は 1-based 列 `[p*NC+1, min((p+1)*NC, N)]` を所有。
- `Ham_local` は**連続カラムメジャー1次元配列**（`double complex*`）、要素 (i,j)（1-based）→ `Ham_local[(j - HamColBegin) * HamPanelLd + (i - 1)]`、`HamPanelLd = idim_max`。所有列0のランクも最低1要素を確保。
- 書き込みマクロはデバッグビルド（`NDEBUG`未定義）で列所有をアサート。
- `Solver 3` かつ `nproc > 1` かつ `OutputHam`/`InputHam` は readdef 起動時エラー（スペック§2で「フェーズ2マージ時に有効化」と定めた検証）。
- 分散モードでは `L_vec` を確保しない（固有ベクトルは `Z_vec` 分散のまま、回収は既存 `GetEigenVectorBlock`）。
- OpenMP `default(none)` プラグマの shared リストに、マクロが参照する新グローバルを必ず追加（漏れはコンパイルエラーになるので検出可能）。
- 警告・エラー表示は rank 0 のみ（`stdoutMPI`）。
- コミットメッセージ末尾に本セッションの Co-Authored-By/Claude-Session トレーラを付ける（フェーズ1と同じ）。
- ローカル（macOS, ELPAなし）の合格基準: 既定ビルドと `build_noMPI` の `ctest -R fulldiag` 16/16 維持。ELPA実行系の検証は最終タスクで clavius（`~/HPhi-elpa`, conda env `hphi_elpa`, `~/opt/elpa-2025.06-cuda`）にて実施。

## 事前監査結果（計画時点で確認済みの事実）

- `Ham` への書き込み箇所: `src/makeHam.c`（ゼロ初期化 :90-93、対角 :97、散布 :122-249 = 全て `Ham[...][j]` 列局所）、`src/nbody_interall.c:1851,1879,1914,1918`（列局所）、`src/anomalous_pair.c:445`（列局所）、`src/input.c:52-53`（**非局所** — `Ham[i][j]` と `Ham[j][i]` の両方 = InputHam 経路。分散モードでは起動時に禁止するので変換不要）。
- `Ham` の読み出し: `src/output.c:87-103`（OutputHam — 同じく禁止）、`src/lapack_diag.c:150-153`（1-based→0-based シフト。分散モードではシフト不要でパネルを直接使う）。
- `makeHam.c:94` の対角ループは `v0[j]/v1[j]` の全域初期化と `Ham[j][j]` 書き込みが同居 — ループ範囲は全域のまま、`Ham` 書き込みだけ所有ガードで間引く。
- スペック§3の「全ての書き込みが makeHam.c:92-249」という記述は不正確（nbody_interall/anomalous_pair にもある。ただし全て列局所なので設計は成立）— Task 1 でスペックを訂正する。

---

### Task 1: 格納マクロとパネルグローバル（`hamstore.h`）＋スペック訂正

**Files:**
- Create: `src/include/hamstore.h`
- Modify: `src/include/global.h`（`Ham`/`L_vec` 宣言の近く）
- Modify: `src/global.c`（定義追加）
- Modify: `docs/superpowers/specs/2026-07-10-elpa-fulldiag-design.md`（§3フェーズ2の書き込み箇所記述の訂正）

**Interfaces:**
- Produces（後続タスク全てが使用）:
  ```c
  /* global.h / global.c */
  extern double complex *Ham_local; /* 1D column panel (distributed mode) */
  extern long int HamColBegin;      /* first owned column, 1-based (0 = panel inactive) */
  extern long int HamColEnd;        /* last owned column, 1-based inclusive */
  extern long int HamPanelLd;       /* leading dimension = idim_max */
  extern int iHamPanelActive;       /* 1: distributed panel mode */
  ```
  ```c
  /* hamstore.h（マクロ、二重評価に注意して j を1回だけ評価する形にする） */
  HAM_OWNED_COL(j)          /* 真: このランクが列 j (1-based) を所有（非分散時は常に真） */
  AddHamElem(i, j, val)     /* Ham(i,j) += val（1-based。分散時はパネルへ、所有アサート付き） */
  ```

- [ ] **Step 1: hamstore.h を書く**

```c
/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */
/* （GPL v3 ヘッダ: 他ソースと同一の15行） */
#ifndef HPHI_HAMSTORE_H
#define HPHI_HAMSTORE_H

#include <assert.h>
#include "global.h"

/**
 * Storage abstraction for the dense FullDiag Hamiltonian
 * (design doc section 3, phase 2).
 *
 * Replicated mode (iHamPanelActive == 0): writes go to the global
 * Ham[i][j] (1-based, as before).
 * Distributed-panel mode (Solver 3, nproc > 1): each rank stores only
 * its owned column range [HamColBegin, HamColEnd] (1-based, inclusive)
 * in the contiguous column-major panel Ham_local; element (i, j) maps
 * to Ham_local[(j - HamColBegin) * HamPanelLd + (i - 1)].
 *
 * Every Hamiltonian write in the generation code MUST go through
 * AddHamElem (or be guarded by HAM_OWNED_COL); a direct Ham[i][j]
 * write silently corrupts nothing in panel mode (Ham is not allocated)
 * but crashes on NULL — the debug assert below catches ownership bugs
 * before that.
 */

#define HAM_OWNED_COL(jcol)                                            \
  (!iHamPanelActive ||                                                 \
   ((long int)(jcol) >= HamColBegin && (long int)(jcol) <= HamColEnd))

#define AddHamElem(irow, jcol, val)                                    \
  do {                                                                 \
    long int hs_i_ = (long int)(irow);                                 \
    long int hs_j_ = (long int)(jcol);                                 \
    if (iHamPanelActive) {                                             \
      assert(hs_j_ >= HamColBegin && hs_j_ <= HamColEnd);              \
      Ham_local[(hs_j_ - HamColBegin) * HamPanelLd + (hs_i_ - 1)]      \
        += (val);                                                      \
    } else {                                                           \
      Ham[hs_i_][hs_j_] += (val);                                      \
    }                                                                  \
  } while (0)

#endif /* HPHI_HAMSTORE_H */
```

- [ ] **Step 2: グローバルを追加**

`src/include/global.h` の `Ham`/`L_vec` の extern 宣言（`grep -n "L_vec" src/include/global.h` で位置確認）の直後に上記 Interfaces の extern 5行を追加。`src/global.c` の対応箇所（`Ham` の定義近く）に:

```c
double complex *Ham_local = NULL;
long int HamColBegin = 0;
long int HamColEnd = -1;
long int HamPanelLd = 0;
int iHamPanelActive = 0;
```

- [ ] **Step 3: スペック訂正**

設計文書 §3 フェーズ2の「`makeHam.c:92-249` で全書き込み箇所を確認済み」の段落に、
「（実装時の再監査で `nbody_interall.c` / `anomalous_pair.c` にも列局所の書き込みがあることを確認。
`input.c` の InputHam 読み込みのみ非列局所であり、分散モードでは起動時に禁止される）」と追記。

- [ ] **Step 4: ビルド確認とコミット**

標準ビルド確認（`cd build && cmake .. && make HPhi -j4` → `[100%] Built target HPhi`。
未使用の新グローバルが増えるだけなので通る）。

```bash
git add src/include/hamstore.h src/include/global.h src/global.c \
        docs/superpowers/specs/2026-07-10-elpa-fulldiag-design.md
git commit -m "Add Hamiltonian storage abstraction for distributed panel mode"
```

---

### Task 2: readdef — OutputHam/InputHam × 分散モードの起動時拒否

**Files:**
- Modify: `src/readdef.c`（`ResolveSolver()` 呼び出し後の検証ブロック、フェーズ1で追加したスペクトル検証の近く）
- Modify: `src/ErrorMessage.c` / `src/include/ErrorMessage.h`
- Test: `test/fulldiag_solver_keyword.sh`（ケース追加）

**Interfaces:**
- Consumes: `X->iSolver`（SOLVER_ELPA=3, DefCommon.h）、`nproc`（global.h 経由）、`X->iOutputHam`/`X->iInputHam`（既存フィールド）
- Produces: `Solver 3 && nproc > 1 && (iOutputHam || iInputHam)` → エラー終了という起動時保証（Task 3以降のパネルモードは Ham 全体を持たない前提をこれで守る）

- [ ] **Step 1: 失敗するテストを書く**

`test/fulldiag_solver_keyword.sh` に、ケース(3b)の後に追加（このビルドが ELPA 対応か否かに関わらず動くよう、非ELPAビルドでは Solver 3 自体が落ちることに注意 — ELPAビルドでのみ意味を持つ分岐にする）:

```sh
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
```

`test/fulldiag_elpa_hubbard_chain.sh` の末尾（比較の後）にマルチプロセス拒否ケースを追加:

```sh
# OutputHam is incompatible with distributed generation (Solver 3, nproc>1):
# must be rejected at startup with a clear message.
cd ..
mkdir -p fulldiag_elpa_hamio_reject/
cd fulldiag_elpa_hamio_reject
cp ../fulldiag_elpa_hubbard_chain/stan.in .
../../src/HPhi -sdry stan.in
printf "Solver  3\nNGPU    0\nOutputHam  1\n" >> calcmod.def
if ${MPIRUN} ../../src/HPhi -e namelist.def > reject.log 2>&1; then
  echo "ERROR: Solver 3 + OutputHam + nproc>1 must fail at startup"
  exit 1
fi
grep -qi "OutputHam\|InputHam" reject.log

echo "fulldiag_elpa_hubbard_chain: OK"
```
（既存の最後の `echo ... OK` はこの新ケースの後ろに移す。）

- [ ] **Step 2: ローカルで確認**

`build_noMPI` で `ctest -R fulldiag_solver_keyword` → PASS のまま（非ELPAビルドでは新ケースはスキップ分岐）。

- [ ] **Step 3: エラーメッセージと検証を実装**

`src/ErrorMessage.c`（フェーズ1で追加した Solver 系メッセージの直後）:

```c
char *cErrElpaHamIO="Error in %s\n Solver 3 (ELPA) with more than one MPI process generates the Hamiltonian distributed,\n which is incompatible with OutputHam/InputHam. Run with 1 process, or use another solver.\n";
```

`ErrorMessage.h` に `extern char *cErrElpaHamIO;` を追加。

`src/readdef.c` の、フェーズ1で追加したスペクトル×Solver検証ブロックの直後に:

```c
  if (X->iSolver == SOLVER_ELPA && nproc > 1
      && (X->iOutputHam == TRUE || X->iInputHam == TRUE)) {
    fprintf(stdoutMPI, cErrElpaHamIO, defname);
    return (-1);
  }
```
（`TRUE` マクロと `iOutputHam`/`iInputHam` の実際のフィールド名・値は readdef.c 内の既存検証
`if(X->iInputHam == 1 && X->iOutputHam==1)` に合わせること。）

- [ ] **Step 4: 回帰＋コミット**

`build_noMPI`: `ctest -R fulldiag` 16/16。

```bash
git add src/readdef.c src/ErrorMessage.c src/include/ErrorMessage.h \
        test/fulldiag_solver_keyword.sh test/fulldiag_elpa_hubbard_chain.sh
git commit -m "Reject OutputHam/InputHam with distributed ELPA FullDiag at startup"
```

---

### Task 3: xsetmem — 分散モードのパネル確保（`Ham`/`L_vec` を確保しない）

**Files:**
- Modify: `src/xsetmem.c`（FullDiag 用確保ブロック、`Ham = cd_2d_allocate(...)` / `L_vec = cd_2d_allocate(...)` の箇所 = xsetmem.c:286-287 付近）

**Interfaces:**
- Consumes: Task 1 のグローバル、`X->Def.iSolver`、`nproc`、`X->Check.idim_max`
- Produces: 分散モードで `Ham_local`（`idim_max * ncols` 要素、最低1）確保・`iHamPanelActive=1`・
  `HamColBegin/HamColEnd/HamPanelLd` 設定済み、`Ham`/`L_vec` は **NULL のまま**。
  複製モードでは従来どおり `Ham`/`L_vec` を確保し `iHamPanelActive=0`。

- [ ] **Step 1: 実装**

`xsetmem.c` の該当2行を次のブロックに置き換える（`c_1d_allocate` 相当が無ければ `malloc` 直書きで既存スタイルに合わせる。`#include "hamstore.h"` は不要 — グローバルは global.h 経由。`DefCommon.h` を include して SOLVER_ELPA を参照。`nproc` は global.h 経由で参照可能なことを `grep -n "extern int nproc" src/include/*.h` で確認して合わせる）:

```c
#ifdef _ELPA
      if (X->Def.iSolver == SOLVER_ELPA && nproc > 1) {
        /* Distributed-panel mode (design doc sec. 3 phase 2):
           each rank stores only its owned 1D column block.
           NC = ceil(N/P); rank p owns 1-based columns
           [p*NC+1, min((p+1)*NC, N)]. Ham/L_vec stay unallocated. */
        long int NN = X->Check.idim_max;
        long int NC = (NN + nproc - 1) / nproc;
        long int jb = (long int)myrank * NC + 1;
        long int je = ((long int)myrank + 1) * NC;
        if (je > NN) je = NN;
        if (jb > NN) { jb = 1; je = 0; } /* rank owns no column */
        HamColBegin = jb;
        HamColEnd = je;
        HamPanelLd = NN;
        iHamPanelActive = 1;
        {
          long int ncols = (je >= jb) ? (je - jb + 1) : 0;
          long int nelem = NN * ncols;
          if (nelem < 1) nelem = 1;
          Ham_local = (double complex *)malloc(nelem * sizeof(double complex));
          if (Ham_local == NULL) return -1;
          for (long int k = 0; k < nelem; k++) Ham_local[k] = 0.0;
        }
      } else
#endif
      {
        Ham = cd_2d_allocate(X->Check.idim_max + 1, X->Check.idim_max + 1);
        L_vec = cd_2d_allocate(X->Check.idim_max + 1, X->Check.idim_max + 1);
      }
```

直後の NULL チェックループ（`if (Ham[j] == NULL || L_vec[j] == NULL)` xsetmem.c:298 付近）を
`if (!iHamPanelActive)` で囲む（パネルモードでは Ham/L_vec は NULL が正）。
`myrank` が xsetmem.c から見えるか確認（global.h の extern。見えなければ include 追加）。

**解放経路の監査（必須）**: `grep -rn "cd_2d_free\|free.*Ham\|free.*L_vec" src/*.c` で
`Ham`/`L_vec` を解放している箇所（xsetmem のメモリ解放関数や終了処理）を特定し、
パネルモードで NULL の `Ham`/`L_vec` を解放しないよう `if (!iHamPanelActive)` ガード
（または `cd_2d_free` が NULL 安全ならその旨をレポートに記載）を入れる。
`Ham_local` は lapack_diag_elpa 内で消費時に解放される（Task 5）が、
対角化前に異常終了した場合に備え、同じ解放箇所に `free(Ham_local)`（NULL 安全）も追加する。

- [ ] **Step 2: ビルド＋回帰＋コミット**

標準ビルド確認 + `build_noMPI` `ctest -R fulldiag` 16/16（非ELPAビルドではブロックはコンパイルされない）。

```bash
git add src/xsetmem.c
git commit -m "Allocate distributed Hamiltonian panel instead of replicated Ham/L_vec"
```

---

### Task 4: 生成コードのマクロ化と担当列制限（makeHam / nbody_interall / anomalous_pair）

**Files:**
- Modify: `src/makeHam.c`（全 `Ham[` 書き込み + j ループ範囲 + OpenMP shared リスト）
- Modify: `src/nbody_interall.c:1851,1879,1914,1918` と当該関数の j ループ
- Modify: `src/anomalous_pair.c:445` と当該関数の j ループ

**Interfaces:**
- Consumes: Task 1 の `AddHamElem`/`HAM_OWNED_COL` と Task 3 のパネル確保
- Produces: 分散モードで各ランクが担当列のみ O(N²/P) 生成。複製モードは従来と bit 同一の `Ham`。

- [ ] **Step 1: 書き込み箇所の完全列挙（監査）**

```bash
grep -rn "Ham\[" src/*.c | grep -v "Ham_local\|lapack_diag\|output.c\|input.c\|xsetmem\|global.c"
grep -rn "L_vec" src/*.c src/include/*.h
```
`Ham[` の全行がこのタスクの変換対象（計画時点の既知: makeHam.c 多数 + nbody_interall.c 4 + anomalous_pair.c 1）。
`L_vec` の各参照は「パネルモードで到達し得ないこと」を確認してレポートに列挙する
（計画時点の監査: phys.c は `use_scalapack` 分岐で回避済み、CalcSpectrumByFullDiag は
readdef ゲートで到達不能、lapack_diag の非ELPA case は iSolver 分岐で到達不能、
xsetmem は Task 3 で分岐済み — 新顔があれば同様に確認）。

- [ ] **Step 2: makeHam.c の変換**

1. `#include "hamstore.h"` を追加。
2. ゼロ初期化ループ（makeHam.c:90-93 `Ham[i][j] = 0;`）を分岐化:
   ```c
   if (!iHamPanelActive) {
     for (i = 0; i <= i_max; i++)
       for (j = 0; j <= i_max; j++)
         Ham[i][j] = 0;
   }
   /* panel mode: Ham_local was zero-initialized at allocation (xsetmem) */
   ```
3. 対角ループ（makeHam.c:94-97）: ループは全域のまま `Ham[j][j] += list_Diagonal[j];` を
   ```c
   if (HAM_OWNED_COL(j)) AddHamElem(j, j, list_Diagonal[j]);
   ```
   に置換。OpenMP プラグマの shared リストに `iHamPanelActive, HamColBegin, HamColEnd, HamPanelLd, Ham_local` を追加（`default(none)` なので漏れはコンパイルエラー）。
4. 主要項ループ: `Ham[X + 1][j] += dmv;` 型の全書き込みを `AddHamElem(X + 1, j, dmv)` に置換。
   各書き込みを囲む状態ループ `for (j = 1; j <= i_max; j++)` は、ループ本体が**列 j への書き込みだけ**を行う場合、
   ```c
   long int hs_jb = iHamPanelActive ? HamColBegin : 1;
   long int hs_je = iHamPanelActive ? HamColEnd : (long int)i_max;
   for (j = hs_jb; j <= hs_je; j++) { ... }
   ```
   に置換して O(N²/P) 化する（`hs_jb > hs_je` = 所有0列ランクはループが自然にスキップされる）。
   ループ本体に `v0`/`v1` 等の全域副作用が同居する箇所（対角ループのみ、確認済み）は範囲を変えずガード方式にする。
5. 変換後、`grep -n "Ham\[" src/makeHam.c` が 2（ゼロ初期化の分岐内のみ）になることを確認。

- [ ] **Step 3: nbody_interall.c / anomalous_pair.c の変換**

同じパターン: `#include "hamstore.h"`、4+1 箇所の `Ham[...][j] +=` を `AddHamElem(..., j, ...)` に、
それぞれを囲む j 状態ループ（各書き込み行から上に遡って特定する）を `hs_jb..hs_je` 範囲に制限。
ループが `M_Ham` モード以外（Lanczos の行列ベクトル積等）と共用されている場合は、
**範囲制限は `X->Large.mode == M_Ham && iHamPanelActive` のときだけ**適用する条件付き境界にする
（`hs_jb`/`hs_je` の初期化式に mode 判定を含める）。共用の有無は関数を読んで判断し、レポートに記載。

- [ ] **Step 4: 回帰（複製経路の bit 同一性）**

`build_noMPI`: `ctest -R fulldiag` 16/16（これが「複製モードで出力不変」の確認。
FullDiag テストは固有値・グリーン関数を参照値と照合するので変換ミスを検出する）。
標準ビルド確認も実施。

- [ ] **Step 5: コミット**

```bash
git add src/makeHam.c src/nbody_interall.c src/anomalous_pair.c
git commit -m "Route Hamiltonian writes through storage macros with owned-column loops"
```

---

### Task 5: 再分散 `RedistPanelToBlockCyclic()` と ELPA 経路のパネル対応

**Files:**
- Modify: `src/matrixscalapack.c`（`GetEigenVectorBlock` 群の近くに追加）
- Modify: `src/include/matrixscalapack.h`（宣言）
- Modify: `src/lapack_diag.c`（`lapack_diag_elpa()` のパネル分岐）

**Interfaces:**
- Consumes: Task 1/3/4 の成果、フェーズ1の `diag_elpa_cmp(..., int nblk, int ngpu)` と `ElpaBlockSize()`
- Produces:
  ```c
  /* 1D column panel (this rank owns ncols_panel columns starting at
     global 1-based column jbegin) -> 2D block-cyclic A_distr(descA_2d).
     All ranks call collectively. Returns 0. */
  int RedistPanelToBlockCyclic(long int xNsize, long int jbegin,
                               long int ncols_panel, long int panel_ld,
                               double complex *panel,
                               double complex *A_distr, int *descA_2d);
  ```

- [ ] **Step 1: RedistPanelToBlockCyclic を実装**

`src/matrixscalapack.c` 末尾（`#endif` 前）に追加。整数幅は既存宣言（`descinit_` 等）に厳密に合わせる:

```c
/**
 * @brief Redistribute a 1D column-block panel into an existing 2D
 * block-cyclic matrix with a single pzgemr2d_ call (design doc sec. 3
 * phase 2). The 1D source: a 1 x P grid ('R'), MB = N (all rows in one
 * block), NB = NC = ceil(N/P) (one column block per rank), RSRC=CSRC=0,
 * LLD = panel_ld. Ranks owning zero columns still participate (their
 * local part is an unused 1-element buffer).
 */
int RedistPanelToBlockCyclic(long int xNsize, long int jbegin,
                             long int ncols_panel, long int panel_ld,
                             double complex *panel,
                             double complex *A_distr, int *descA_2d) {
  int i_negone = -1, i_zero_i = 0, info;
  int ictxt_1d, nprow_1, npcol_1, myrow_1, mycol_1;
  int desc1d[9];
  long int NC, lld;
  const long int i_one = 1;
  int size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  NC = (xNsize + size - 1) / size;

  blacs_get_(&i_negone, &i_zero_i, &ictxt_1d);
  nprow_1 = 1; npcol_1 = size;
  blacs_gridinit_(&ictxt_1d, "R", &nprow_1, &npcol_1);
  blacs_gridinfo_(&ictxt_1d, &nprow_1, &npcol_1, &myrow_1, &mycol_1);

  lld = (panel_ld > 0) ? panel_ld : 1;
  {
    long int mb1 = xNsize, nb1 = NC;
    descinit_(desc1d, &xNsize, &xNsize, &mb1, &nb1, &i_zero_i, &i_zero_i,
              &ictxt_1d, &lld, &info);
  }

  pzgemr2d_((long int *)&xNsize, (long int *)&xNsize,
            panel, (long int *)&i_one, (long int *)&i_one, desc1d,
            A_distr, (long int *)&i_one, (long int *)&i_one, descA_2d,
            &descA_2d[1]);

  blacs_gridexit_(&ictxt_1d);
  (void)jbegin; (void)ncols_panel;
  return 0;
}
```

注意（実装時に必ず確認・調整）:
- `descinit_` の引数幅は既存呼び出し（`diag_scalapack_cmp`）と同一に（`&xNsize` は long int、
  `RSRC/CSRC/ictxt/lld/info` は既存宣言の幅）。合わなければ既存に従いローカル変数の型を変える。
- `pzgemr2d_` の最終引数は両グリッドを包含するコンテキスト = 2D側 `descA_2d[1]`（フェーズ1の
  `GetEigenVectorBlock` と同じ規約）。
- rank の所有列は 1D グリッドの `mycol` と一致する設計（xsetmem の `myrank*NC+1` と
  BLACS 'R' 1×P グリッドの列番号が一致すること）— `assert(mycol_1 == myrank)` 相当の
  デバッグチェックを入れる（`myrank` は global）。

`matrixscalapack.h` に宣言を追加（`GetEigenVectorBlock` の隣、同じ整数幅表記で）。

- [ ] **Step 2: lapack_diag_elpa のパネル分岐**

`src/lapack_diag.c` の `lapack_diag_elpa()` で、`DivMat` 充填ループ（現行）を分岐化:

```c
  if (iHamPanelActive) {
    /* Distributed generation (phase 2): redistribute the 1D panel and
       free it before ELPA to lower the memory peak. */
    RedistPanelToBlockCyclic(xMsize, HamColBegin,
                             (HamColEnd >= HamColBegin)
                               ? (HamColEnd - HamColBegin + 1) : 0,
                             HamPanelLd, Ham_local, A_distr, descA);
    free(Ham_local);
    Ham_local = NULL;
    iHamPanelActive = 0; /* panel consumed; phys.c uses Z_vec only */
  } else {
    for (i = 0; i < xMsize; i++) {
      for (j = 0; j < xMsize; j++) {
        DivMat(i, j, Ham[i][j], A_distr, descA);
      }
    }
  }
```

重要: 複製モードの `lapack_diag()` 冒頭には 1-based→0-based シフト（`Ham[i][j] = Ham[i+1][j+1]`,
lapack_diag.c:150-153）がある。パネルは 1-based 生成だが `Ham_local[(j-jb)*ld + (i-1)]` で行は
既に 0-based 詰めになっており、列も先頭詰めなのでシフト不要。ただし `lapack_diag()` の
シフトループ自体が `Ham` を触るので、**シフトループを `if (!iHamPanelActive)` で囲む**こと
（パネルモードで Ham は NULL）。

- [ ] **Step 3: ビルド＋回帰＋コミット**

標準ビルド確認 + `build_noMPI` 16/16。

```bash
git add src/matrixscalapack.c src/include/matrixscalapack.h src/lapack_diag.c
git commit -m "Redistribute the generated 1D panel into the ELPA matrix"
```

---

### Task 6: 再分散の検証テスト（複製経路とパネル経路の全要素一致）

**Files:**
- Create: `test/unit/elpa_redist_check.c`
- Modify: `test/CMakeLists.txt`（`if(USE_ELPA)` ブロック内に登録）

**Interfaces:**
- Consumes: `RedistPanelToBlockCyclic`（Task 5）、`DivMat`、`ElpaBlockSize`、BLACS 基盤
- Produces: 実行ファイル `elpa_redist_check`（MPI 実行、exit 0/1）— N=97（P・nblk 非整除）で
  (a) 複製→`pzelset` 詰めの 2D 行列と (b) 1Dパネル→`pzgemr2d` の 2D 行列を全ローカル要素比較。

- [ ] **Step 1: テストプログラムを書く**

```c
/* Verify RedistPanelToBlockCyclic: build the same deterministic matrix
   (i) replicated + DivMat/pzelset and (ii) as a 1D column panel +
   pzgemr2d, then compare the resulting 2D block-cyclic local arrays
   element-wise. N=97 exercises non-divisible N vs both the process
   count and the (capped) block size. Run with any np. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>
#include "matrixscalapack.h"
#include "matrixlapack_elpa.h"

#define NDIM 97

static double complex MatElem(long int i, long int j) { /* 0-based */
  return (1.0 / (1.0 + labs(i - j))) + ((double)(i - j) / (NDIM * NDIM)) * I;
}

int main(int argc, char **argv) {
  int i_negone = -1, i_zero_i = 0;
  int rank, size, ictxt, iam, nprocs, info, ok = 1;
  int nprow, npcol, myrow, mycol;
  long int n = NDIM, mb, mp, nq, i, j, k;
  int lld, dims[2] = {0, 0};
  int descA[9], descB[9];
  double complex *A_ref, *B_panel2d, *panel;
  long int NC, jb, je, ncols, lde;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size, 2, dims);
  nprow = dims[0]; npcol = dims[1];
  mb = ElpaBlockSize(n, nprow, npcol);
  blacs_pinfo_(&iam, &nprocs);
  blacs_get_(&i_negone, &i_zero_i, &ictxt);
  blacs_gridinit_(&ictxt, "R", &nprow, &npcol);
  blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);

  mp = numroc_(&n, &mb, &myrow, &i_zero_i, &nprow);
  nq = numroc_(&n, &mb, &mycol, &i_zero_i, &npcol);
  lld = (mp > 0) ? (int)mp : 1;
  descinit_(descA, &n, &n, &mb, &mb, &i_zero_i, &i_zero_i, &ictxt, &lld, &info);
  descinit_(descB, &n, &n, &mb, &mb, &i_zero_i, &i_zero_i, &ictxt, &lld, &info);

  A_ref = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));
  B_panel2d = malloc(((mp * nq > 0) ? mp * nq : 1) * sizeof(double complex));

  /* (i) replicated fill */
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      DivMat(i, j, MatElem(i, j), A_ref, descA);

  /* (ii) 1D column panel fill: rank owns 0-based columns [jb, je) */
  NC = (n + size - 1) / size;
  jb = (long int)rank * NC;
  je = jb + NC; if (je > n) je = n;
  ncols = (je > jb) ? (je - jb) : 0;
  lde = n;
  panel = malloc((((n * ncols) > 0) ? n * ncols : 1) * sizeof(double complex));
  for (j = jb; j < je; j++)
    for (i = 0; i < n; i++)
      panel[(j - jb) * lde + i] = MatElem(i, j);

  RedistPanelToBlockCyclic(n, jb + 1, ncols, lde, panel, B_panel2d, descB);

  for (k = 0; k < mp * nq; k++) {
    if (cabs(A_ref[k] - B_panel2d[k]) > 1e-14) {
      fprintf(stderr, "rank %d: mismatch at local index %ld: %g\n",
              rank, k, cabs(A_ref[k] - B_panel2d[k]));
      ok = 0; break;
    }
  }
  { int gok; MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); ok = gok; }
  if (rank == 0) printf("elpa_redist_check: %s\n", ok ? "OK" : "FAILED");
  free(A_ref); free(B_panel2d); free(panel);
  blacs_gridexit_(&ictxt);
  MPI_Finalize();
  return ok ? 0 : 1;
}
```

- [ ] **Step 2: 登録**

`test/CMakeLists.txt` の `if(USE_ELPA)` ブロック内、`elpa_eigen_check` の登録に倣って
（ソースは `unit/elpa_redist_check.c` + `matrixscalapack.c` + `matrixlapack_elpa.c`、
同じ include/definitions/正規化済み `_ec_sc_libs` リンク、`run_with_mpi_precheck.sh min:2` +
`SKIP_RETURN_CODE 77`）登録する。フェーズ1と同じ変数を再利用。

- [ ] **Step 3: ローカル確認＋コミット**

既定ビルドで `ctest -N | grep -i elpa` → 空（USE_ELPA=OFF）。標準ビルド確認。

```bash
git add test/unit/elpa_redist_check.c test/CMakeLists.txt
git commit -m "Add redistribution verification test for the 1D panel path"
```

---

### Task 7: clavius 実機検証（分散生成のエンドツーエンド）

**Files:** なし（実行と記録。結果は `test/manual/elpa_gpu_check.md` に追記してコミット）

- [ ] **Step 1**: rsync（**リポジトリルートから**。フェーズ1の教訓: cwd に注意）→ `build_elpa`（CPU, conda env `hphi_elpa`。CPU版elpaはenvから削除済みなので `-DELPA_ROOT=$HOME/opt/elpa-2025.06-cuda` を使い、実行時 `LD_LIBRARY_PATH=$HOME/opt/elpa-2025.06-cuda/lib:$CONDA_PREFIX/lib`）で再ビルド。
- [ ] **Step 2**: `elpa_redist_check` を np=1,2,3,4,6,8 で実行 → 全 OK。
- [ ] **Step 3**: `elpa_eigen_check`・`fulldiag_elpa_hubbard_chain`（新 OutputHam 拒否ケース込み）・`fulldiag_solver_keyword` を再実行 → PASS。
- [ ] **Step 4**: 分散生成の物理検証: L=6 Hubbard 鎖（N=400）を `Solver 3`/np=4 で実行し `Solver 0` と全固有値比較（フェーズ1と同じ手順、今回はパネル経路が使われる）。L=8（N=4900）np=2,4 でも energy/doublon 比較。
- [ ] **Step 5**: **メモリスケーリング確認**（スペック§6）: L=8（N=4900）で np=1（複製: Ham+L_vec ≈ 2×16N² ≈ 770MB）と np=4 分散（パネル+分散行列）を `/usr/bin/time -v`（Linux）の Maximum resident set size で比較し、分散側の1ランクあたり常駐が大きく下がることを記録。
- [ ] **Step 6**: GPU（両GPU空きを nvidia-smi で確認の上、ルール遵守）: L=8 `Solver 3`/`NGPU 1`/np=2（分散生成→GPU対角化）で energy 一致確認。
- [ ] **Step 7**: 結果を `test/manual/elpa_gpu_check.md` に「Phase 2 validation」節として追記しコミット。

---

## 完了条件（フェーズ2）

- `Solver 0/1/2` と `Solver 3`+`nproc==1` の既存テストが無変更で PASS（複製経路 bit 同一）。
- ELPA ビルドで `elpa_redist_check`（np 複数）・既存 ELPA テスト・OutputHam 拒否ケースが PASS。
- clavius で分散生成の固有値一致（複数 np）とメモリ削減を記録。
- スペック§3フェーズ2の記述（書き込み箇所の訂正含む）と実装が一致。
