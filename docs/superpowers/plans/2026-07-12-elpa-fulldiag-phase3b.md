# ELPA FullDiag フェーズ3b（ExpecMode 2 / トレースカーネル）実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `ExpecMode 2` を実キーワードにする — 演算子ごとに基底写像 (k→k′, 振幅) を 1 回前計算し、所有全状態を密ループでストリーミング評価するトレースカーネル（`src/expec_trace.c`）を、モデル別の段階的有効化とフォールバック付きで実装する。

**Architecture:** 設計文書 `docs/superpowers/specs/2026-07-11-elpa-fulldiag-phase3-design.md`（v5.1、§3「Mode 2」が正）のフェーズ3b。フェーズ3a の Mode 1 基盤（状態パネル、ExpecLocal、パーシャル+マニフェスト出力）は完成済みで、3b は `phys_stateparallel()` 内のカーネル選択として載る。①スケルトン+ケイパビリティ表+INFO（全量フォールバック=挙動は Mode 1 と同一）→ ②写像抽出エンジン+純粋性テスト → ③対角量・エネルギー系カーネル → ④一体GFカーネル → ⑤二体GFカーネル → ⑥ゴールデンテストで合格モデルのみ表を TRUE 化 → ⑦equivテスト拡張 → ⑧docs → ⑨clavius検証+ベンチ。

**Tech Stack:** C99, MPI（オーケストレーション層のみ）, 既存 mltply*Core 要素関数, フェーズ3a の green_output/ExpecLocal/状態パネル基盤, CMake/ctest。

## Global Constraints

- スペック §3「Mode 2: トレースカーネル」が正。判断に迷ったらスペック。**演算子代数を再実装しない** — 写像は「汎用経路が呼ぶのと同一の要素関数を k ごとに 1 回呼んで記録」する（規約を構造ごと継承）。
- **保証**: `ExpecMode` は速度のみを変える。0/1/2 の結果は丸め誤差の範囲で一致（1e-8、bit 一致ではない）。カーネル担当量とフォールバック量の書き手は**常に一意**（二重書き込み・部分上書きの構造的排除、スペック §3 所有表）。
- モデル別段階的有効化: ゴールデンテスト（Task 6）合格モデルのみケイパビリティ表で有効。初期対象は **Hubbard（正準+GC）と Spin/SpinGC の half**。汎用スピン・spinless・Kondo は全量フォールバック。`all_s2` は 3b では**全モデルでフォールバック**（expec_totalspin 経由。表の「トレースカーネル（対象モデル）」化は 3c 以降の拡張余地として明記）。
- `src/expec_trace.c` は **MPI フリー**（mpi.h include 禁止、生 MPI・exitMPI 禁止、許可 wrapperMPI は SumMPI_dc/d/li/i, fopenMPI, childfopenMPI, stdoutMPI のみ）。作成と同一コミットで `test/check_expec_local_calls.sh` の `FILES` に追加。集団操作はすべて `src/phys_distributed.c`（恒久的にスキャン対象外）。
- 要素関数呼び出しは **`X->Large.mode = M_CORR` のまま**行う（tmp_v0 への副作用が構造的に無効になる。mltplyCommon.h: M_MLTPLY=0, M_ENERGY=1, M_CORR=3, M_CALCSPEC=4, H_CORR=5）。H_CORR は使わない。
- 既定（`ExpecMode 0/1`）の挙動は不変。build_noMPI（`cmake -DENABLE_MPI=OFF ..`）の `ctest -R "fulldiag|check_expec"` 18/18 基準。ELPA 実行系はローカルでは登録のみ、実行は Task 9（clavius, `~/HPhi-elpa`, conda env `hphi_elpa`, `LD_LIBRARY_PATH=$HOME/opt/elpa-2025.06-cuda/lib:$CONDA_PREFIX/lib`, **`export CUDA_MPS_PIPE_DIRECTORY=/tmp/nonexistent-mps-hphi` 必須**、rsync はリポジトリルートから）。
- GC 系要素関数は list_1 間接参照なし（基底添字=ビットパターン）→ **GC モデルの単体テストはローカル MPI ビルドで実行可能**。正準系（list_1/GetOffComp 要）は HPhi 実行経由のゴールデンテスト（clavius）で検証。
- コミットメッセージ末尾に本セッションの Co-Authored-By/Claude-Session トレーラ。

## File Structure

- Create: `src/expec_trace.c` — ケイパビリティ表・写像抽出・ストリーミング・カーネル本体（MPI フリー）
- Create: `src/include/expec_trace.h` — 公開 API（下記 Interfaces）
- Modify: `src/phys.c` — Mode 2 ダウングレード INFO の撤去（分岐は phys_stateparallel 内へ）
- Modify: `src/phys_distributed.c` — Mode 2 のカーネル選択+INFO 集約表示
- Modify: `src/phys_distributed_local.c` — 状態ループの量別スキップ（カーネル担当量はフォールバック実行しない）
- Modify: `src/CMakeLists.txt`, `test/CMakeLists.txt`, `test/check_expec_local_calls.sh`（FILES に expec_trace.c）
- Test: `test/unit/expec_trace_map_check.c`（GC 写像・純粋性、ローカル実行可）, `test/fulldiag_expecmode_equiv.sh` 拡張
- Docs: CalcMod ja/en, `docs/superpowers/specs/2026-07-11-phase3a-migration-note.md` の 3b 追記（新規ファイル `...-phase3b-migration-note.md`）

---

### Task 1: expec_trace スケルトン + ケイパビリティ表 + ディスパッチ配線（全量フォールバック）

**Files:**
- Create: `src/expec_trace.c`, `src/include/expec_trace.h`
- Modify: `src/phys.c:90-96`（ダウングレード INFO 撤去）, `src/phys_distributed.c`, `src/phys_distributed_local.c`, `src/CMakeLists.txt`, `test/check_expec_local_calls.sh`（FILES へ追加）
- Test: `test/fulldiag_solver_keyword.sh` は不変（nproc==1 降格は readdef のまま）。ローカル受け入れ = ビルド 2 系 + ガード + 回帰 18/18

**Interfaces:**
- Produces（Task 2-7 が依存する公開 API、`src/include/expec_trace.h`）:
  ```c
  typedef enum {
    TRACE_Q_ENERGY_FLCT = 0,  /* energy+var+doublon+num+sz (expec_energy_flct 相当一式) */
    TRACE_Q_ONEBODY,          /* expec_cisajs 相当 */
    TRACE_Q_TWOBODY,          /* expec_cisajscktaltdc 相当 */
    TRACE_Q_NQUANT
  } TraceQuantity;
  /* 1=カーネル担当 / 0=フォールバック（Mode 1 経路）。モデル×量の定数表を参照 */
  int TraceKernelAvailable(const struct BindStruct *X, TraceQuantity q);
  /* rank 0 用: 量ごとの担当（kernel/fallback）を "  INFO: ExpecMode 2 ..." 形式で
     fp に出力（呼ぶのはオーケストレーション層。expec_trace.c 自身は表示しない） */
  void TraceKernelReportPlan(const struct BindStruct *X, FILE *fp);
  /* 所有状態一括評価: panel は列優先 ld=NN の所有状態パネル（phys_stateparallel と同じ）。
     担当量のみ評価し X->Phys.all_* と GF 出力へ書く。戻り値 0/-1（ローカル rc） */
  int expec_trace_owned_states(struct BindStruct *X, const double complex *panel,
                               long int jb, long int je, long int NN);
  ```
- Consumes: 3a の `phys_stateparallel()`（phys_distributed.c:43-187 の構造）と `phys_stateparallel_local_loop()`（phys_distributed_local.c:67-92 の per-state ループ）。

- [ ] **Step 1: ケイパビリティ表とスケルトンを書く**

`src/expec_trace.c` に定数表（Task 6 で TRUE 化するまで全 FALSE）:

```c
/* model × quantity capability matrix. A row is enabled ONLY after the
   model passes the golden cross-check tests (see plan Task 6 / spec §3.5).
   NBodyG / AnomalousG / totalspin(S2) are ALWAYS fallback in phase 3b and
   deliberately have no row here. */
typedef struct { int calc_model; int flg_general_spin; int q[TRACE_Q_NQUANT]; } TraceCap;
static const TraceCap kTraceCap[] = {
  /* model,      genspin, {ENERGY_FLCT, ONEBODY, TWOBODY} */
  { Hubbard,     0,       {0, 0, 0} },
  { HubbardGC,   0,       {0, 0, 0} },
  { Spin,        0,       {0, 0, 0} },   /* half のみ; genspin=1 は非対応 */
  { SpinGC,      0,       {0, 0, 0} },
  /* Kondo/KondoGC/tJ/tJGC/Spinless*/
};
```

`TraceKernelAvailable` は iCalcModel と iFlgGeneralSpin で表を引き、行が無い/genspin 不一致なら 0。`expec_trace_owned_states` はこの段階では担当量なし（即 return 0）。`TraceKernelReportPlan` は 3 量それぞれについて
`"  INFO: ExpecMode 2: <quantity> uses the <trace kernel|ExpecMode-1 fallback> path.\n"` を出力（quantity 名は `energy/fluctuation`, `one-body Green functions`, `two-body Green functions`。NBodyG/AnomalousG/S2 は常時フォールバックなので `"  INFO: ExpecMode 2: NBodyG/AnomalousG/S2 always use the ExpecMode-1 fallback path.\n"` を 1 行固定で出す）。

- [ ] **Step 2: ディスパッチ配線**

`src/phys.c:91-92` の
```c
    if (X->Def.iExpecMode == EXPECMODE_TRACE)
      fprintf(stdoutMPI, "  INFO: ExpecMode 2 kernels are not available in this build; running as ExpecMode 1.\n");
```
を撤去（`phys_stateparallel(X, neig)` 呼び出しと assert は不変）。`src/phys_distributed.c` の `phys_stateparallel()` に、ローカルループ呼び出しの**直前**（パネル確保・再分散の後）で:

```c
  int use_trace = (X->Def.iExpecMode == EXPECMODE_TRACE);
  if (use_trace && myrank == 0) TraceKernelReportPlan(X, stdoutMPI);
  if (use_trace) {
    rc_local = expec_trace_owned_states(X, panel, jb, je, NN);   /* 担当量 */
    if (rc_local == 0)
      rc_local = phys_stateparallel_local_loop_fallback(X, panel, jb, je, NN);
  } else {
    rc_local = phys_stateparallel_local_loop(X, panel, jb, je, NN);
  }
```

`phys_distributed_local.c` に `phys_stateparallel_local_loop_fallback()` を追加: 既存 `phys_stateparallel_local_loop()` と同一の状態ループだが、**カーネル担当量の evaluator 呼び出しと all_* 代入をスキップ**する（`TraceKernelAvailable(X, TRACE_Q_ENERGY_FLCT)` なら `expec_energy_flct` と all_energy/doublon/num_up/num_down/sz 代入を飛ばす、ONEBODY なら `expec_cisajs` を、TWOBODY なら `expec_cisajscktaltdc` を飛ばす。NBodyG/AnomalousG/totalspin と all_s2 代入は常に実行）。重複を避けるため、両ループは共通の static 関数 `state_loop_impl(X, panel, jb, je, NN, skip_mask)` に統合し、既存 `phys_stateparallel_local_loop` は `skip_mask=0` の薄いラッパにする。

**注意（v0/v1 規約, phys_distributed_local.c:70-81 の事実）**: `expec_energy_flct` をスキップする場合、後続 evaluator が読む `v1` は誰も詰めない。fallback ループでは ENERGY_FLCT スキップ時に**明示的に `v1[j+1]=panel[...]` を直接詰め、`v0` はゼロクリア**する（後続の expec_cisajs/expec_cisajscktaltdc/expec_totalspin は vec=v1 しか読まない — 3a Task 6 で検証済みの事実。expec_nbodyg/expec_anomalousg も同様に v1 を受け取る）。この詰め替えの検証は Task 6 のゴールデンテストが担う。

- [ ] **Step 3: ガード FILES 追加 + ビルド + 回帰**

`test/check_expec_local_calls.sh` の `FILES` に `src/expec_trace.c` を追加。`src/CMakeLists.txt` のソースリストに追加（`#ifdef _SCALAPACK` で空 TU 化は**しない** — expec_trace.c は MPI 非依存の純ローカルコードなので全ビルドでコンパイルし、呼び出し側だけが `_SCALAPACK` 内。ただし list_1 等 extern 参照のみ）。

Run: `cd build_noMPI && cmake .. && make HPhi -j8 && ctest -R "fulldiag|check_expec"`
Expected: 18/18 PASS（この段階の Mode 2 実挙動 = 全量フォールバック = Mode 1 と同一、INFO のみ変化）

- [ ] **Step 4: equiv スクリプトの Mode-2 INFO アサーション更新**

`test/fulldiag_expecmode_equiv.sh` の Mode 2 サブケース（3a で追加、"ExpecMode 2 kernels are not available in this build" を grep している箇所）を、新 INFO 形式（`"ExpecMode 2:"` を含む行が 1 行以上 + 出力一致は従来どおり Mode 1 と比較）に更新。`sh -n` で構文確認。

- [ ] **Step 5: コミット** `git commit -m "Wire ExpecMode 2 dispatch through a trace-kernel capability table"`

---

### Task 2: 写像抽出エンジン + 純粋性テスト（GC 系はローカル実行）

**Files:**
- Modify: `src/expec_trace.c`, `src/include/expec_trace.h`
- Create: `test/unit/expec_trace_map_check.c`
- Modify: `test/CMakeLists.txt`（green_partial_merge_check の登録様式を踏襲: `-DMPI`, min:2 は不要 — **シリアル/np=1 で十分**なので通常の add_test + ラベル unit）

**Interfaces:**
- Produces:
  ```c
  typedef struct {
    long int n;            /* 基底次元 */
    long int *kprime;      /* [n] 0-based 行き先。基底外(消滅)は -1 */
    double complex *amp;   /* [n] 振幅 s_k×結合定数（dam_pr プローブ値そのもの） */
    int is_diagonal;       /* kprime[k]==k 恒等のとき 1 */
  } TraceMap;
  /* ones プローブで要素関数を k=1..n に 1 回ずつ呼び、(k', amp) を記録する。
     戻り値 0/-1。probe_v1/probe_v0 は呼び出し側が確保した n+1 要素の全1配列
     （抽出中に書き換えられないことは純粋性テストの検査対象） */
  int TraceMapExtractOneBody(struct BindStruct *X, int ipair, TraceMap *map,
                             double complex *probe_v0, double complex *probe_v1);
  int TraceMapExtractTwoBody(struct BindStruct *X, int ipair, TraceMap *map,
                             double complex *probe_v0, double complex *probe_v1);
  void TraceMapFree(TraceMap *map);
  ```
- Consumes: 要素関数（下記）と Task 1 のスケルトン。

- [ ] **Step 1: 抽出原理をコメントで固定し、一体 Hubbard/GC から実装**

原理（スペック §3.1 の実装形。**プローブ = 全1ベクトル**）: `X->Large.mode==M_CORR` では全要素関数の `tmp_v0` 書き込みが無効（mltplyHubbardCore.c:387-389 等）で、戻り値は常に `dam_pr = conj(tmp_v1[off]) * tmp_V * sgn * tmp_v1[j]`。`tmp_v1[*]=1` を渡せば **`dam_pr = tmp_V×sgn`（=amp）** が得られ、off 引数（または対角なら k 自身）が k′。消滅（正準系で GetOffComp 失敗、または PauliBlock）は dam_pr==0 → `kprime[k]=-1`。

一体 Hubbard 正準（expec_cisajs.c:409-499 の LOCAL 分岐と同一の要素列）:
```c
  /* pair: X->Def.CisAjt[ipair] = {isite1-1, sigma1, isite2-1, sigma2} */
  general_hopp_GetInfo(X, org_isite1, org_isite2, org_sigma1, org_sigma2);
  is1 = X->Large.is1_spin; is2 = X->Large.is2_spin;   /* GetInfo 直後に snapshot（Large は次の GetInfo で上書きされる） */
  Asum = X->Large.A_spin;  Adiff = X->Large.isA_spin;
  for (k = 1; k <= n; k++) {
    if (diagonal) { map->kprime[k-1] = k-1; map->amp[k-1] = (list_1[k] & is1) ? 1.0 : 0.0; }
    else {
      dam = CisAjt(k, probe_v0, probe_v1, X, is1, is2, Asum, Adiff, 1.0, &off);
      map->kprime[k-1] = (dam != 0.0) ? (long int)off - 1 : -1;   /* off は 1-based */
      map->amp[k-1] = dam;
    }
  }
```
（**注意**: 現行 `CisAjt`（mltplyHubbardCore.c:354-396）は off を out-param で返さない — 内部 `GetOffComp` の結果 `off` はローカル。**要素関数のシグネチャは変更しない**。代わりに正準一体は `child_CisAjt` 系（off を返す薄い関数）を直接使うか、`CisAjt` と同じ列（`list_1[k]` → ビット判定 → `GetOffComp(list_2_1,list_2_2,...,&off)`）を expec_trace.c 内の static ヘルパで**同一 API 呼び出しにより**再構成する。どちらを採るかは実装時に `src/mltplyHubbardCore.c` の child 層 API（`child_CisAjt` が存在するか、off を返すか）を確認して決め、**採った方式と根拠を実装コミットのメッセージに記録**する。GC は `GC_CisAjt(j,...,&tmp_off)` が off を out-param で返す（mltplyHubbardCore.c:403-441）のでそのまま使える。）

一体 Spin-half（対角のみ: child_Spin_CisAis / child_SpinGC_CisAis, mltplySpinCore.c:210-239）と SpinGC-half（対角+横磁場: child_SpinGC_CisAit が off を返す, mltplySpinCore.c:247-271）も同じ形で。

- [ ] **Step 2: 二体の抽出**（expec_cisajscktaltdc.c:899-940 [Hubbard] / :1994-2024 [SpinGC-half] の LOCAL 分岐と同一の要素列）

`Rearray_Interactions(ipair, ..., X, 2)` → `general_int_GetInfo(...)` → 4 分岐（diag/diag, diag/off, off/diag, off/off）の `*_element` 関数を ones プローブで k ごとに呼ぶ。`*_element` は off を out-param（`&tmp_off`）で返すのでそのまま記録。diag/diag（CisAisCisAis 系, off-param なし）は `kprime[k]=k-1`。SpinGC-half の同一サイト縮約 4 分岐（GC_CisAisCisAis_spin_element 等, mltplySpinCore.c:471-620）も同様。Rearray が非ゼロを返すペア（規約外）は `map->n=0` で返し、呼び出し側（Task 5）が Mode 1 と同じ「0.0 行を書く」フォールバック処理をする。

- [ ] **Step 3: 純粋性テストを書く（RED→GREEN）** — `test/unit/expec_trace_map_check.c`

MPI 不要（シリアル）。**GC モデルのみ**（list_1 不要 — HubbardGC と SpinGC-half）。スタブは green_partial_merge_check.c の流儀（myrank=0, stdoutMPI=stdout を file-scope 定義）。手組みの小さい BindStruct: L=4 サイト HubbardGC（n=256）と L=4 SpinGC-half（n=16）、Tpow/Large.i_max 等の最小初期化（`X.Large.mode=M_CORR` 固定）。検査:
1. **写像正当性**: ランダム複素ベクトル z について、`Σ_k conj(z[k′])·amp[k]·z[k]`（kprime≥0 のみ）が、同じ要素関数を M_CORR で z 直接評価した `Σ_k dam_pr(z)` と 1e-13 で一致（一体: 対角 1 個・非対角 1 個・ゼロ結果[常に消滅する組]1 個。二体: 4 分岐各 1 個）。
2. **純粋性（スペック §3.2b）**: 同一演算子で抽出を 2 回実行し kprime/amp が bit 一致、かつ抽出前後で `X->Large` の snapshot（memcmp）と probe_v0/probe_v1 の内容が不変。

Run: `cd build_mpi && make expec_trace_map_check && ./test/expec_trace_map_check`
Expected: 最初は FAIL（関数未実装）→ 実装後 PASS。build_noMPI でも同テストを登録・実行（MPI 非依存）。

- [ ] **Step 4: 回帰 + コミット** `git commit -m "Add trace-map extraction engine with purity checks"`

---

### Task 3: エネルギー系カーネル（対角量 + H 写像による energy/var）

**Files:**
- Modify: `src/expec_trace.c`
- Test: `test/unit/expec_trace_map_check.c` にエネルギー系ケース追加

**Interfaces:**
- Produces: `expec_trace_owned_states()` が ENERGY_FLCT 担当時に `X->Phys.all_energy/all_doublon/all_num_up/all_num_down/all_sz`（インデックス n-1）を所有状態分埋める。
- Consumes: Task 2 の TraceMap。

- [ ] **Step 1: 対角量の重み前計算**

doublon/num_up/num_down/Sz は対角写像: `w_k` を k=1..n で 1 回計算（expec_energy_flct.c の各モデル実装 `expec_energy_flct_Hubbard`(:358-474, list_1 使用)/`_HubbardGC`(:238-351)/`_HalfSpinGC`(:481-546) と**同一のビット演算式**を per-k ヘルパに抽出して使う。Spin 正準は expec_energy_flct.c:140-155 のとおり定数[doublon=0, num=NsiteMPI, Sz=0.5*Total2SzMPI]なのでストリーミング不要）。状態ストリーミング: `q_n = Σ_k w_k·|z_n[k]|²`。

- [ ] **Step 2: energy/var は H 写像ストリーミング**

`<H>` と `<H²>` は **mltply を per-state で呼ぶ代わりに**、ハミルトニアンの全項（Trans/InterAll/CoulombIntra/... — `expec_energy_flct` が呼ぶ `mltply(X,v0,v1)` の項構成）を… **実装最小化の決定**: 3b の ENERGY_FLCT カーネルは `<H>` を**対角化の固有値から採用**し（スペック §3.3 第一文）、`var = <H²>-<H>²` は**フォールバック時のみ**出力する（カーネル担当時は var 検査を実施しない旨を INFO で明示し、`X->Phys.var` は 0 埋めではなく**固有値² を代入**して従来出力の列を壊さない — つまり var 列は定義上 0 になる）。固有値の取得: `lapack_diag.c` の ELPA/ScaLAPACK 経路が `X->Phys.all_energy`… **ではなく** 固有値配列をどのグローバルに置くかを実装時に `src/lapack_diag.c` で確認し（Eigenvalue.dat を書いているコードが読んでいる配列）、`phys_stateparallel()` から `expec_trace_owned_states` へ引数で渡す形に整える（グローバル追加はしない。`double *eigenvalues` 引数を API に追加: `expec_trace_owned_states(X, panel, jb, je, NN, eigenvalues)`。Task 1 の API をこの形に更新するのはこのタスクの冒頭で行い、ヘッダ・呼び出し側・スケルトンを同時に直す）。
**この決定の含意**（レビュー用に明記): Mode 2 の `<H>` は固有値そのもの（Mode 0/1 は `conj(v1)·H·v1` の再計算）— 両者は残差 ~1e-13 で一致し 1e-8 等価性を満たす。var 列は Mode 0/1 の「固有ベクトル品質検査」の意味を失う（0 になる）ため、docs（Task 8）と INFO で「Mode 2 は var 検査を行わない」ことを明示する。

- [ ] **Step 3: ケイパビリティ表で Hubbard/HubbardGC/Spin/SpinGC の ENERGY_FLCT を暫定 TRUE 化**（Task 6 のゴールデン前だが、equiv テストが 3 モードで全量比較するため、TRUE 化はローカル回帰+単体テストの範囲で行い、最終確定は Task 6/9）

- [ ] **Step 4: 単体テストケース追加**（GC 2 モデル: ランダムベクトルで w_k ストリーミング結果が expec_energy_flct_HubbardGC/_HalfSpinGC の per-state 実装と 1e-13 一致）→ RED→GREEN → 回帰 → コミット

---

### Task 4: 一体GFカーネル（バッファリング + 出力）

**Files:**
- Modify: `src/expec_trace.c`
- Test: 単体テストに一体 GF ストリーミングケース追加

**Interfaces:**
- Produces: ONEBODY 担当時、演算子外側×状態内側で `g[ipair][n]` を評価し、**Mode 1 と同一の出力機構・同一のファイル内容**（状態別 `zvo_cisajs_eigen%d.dat` は per-state ローカル書き、集約は GreenOutput パーシャル）で書く。
- Consumes: Task 2 の TraceMapExtractOneBody。

- [ ] **Step 1: バッファ設計を実装**

`double complex *gbuf` サイズ `NCisAjt × ncols`。**メモリゲート**: `NCisAjt*ncols*16 > TRACE_GBUF_MAX_BYTES`（`#define TRACE_GBUF_MAX_BYTES (1UL<<30)` /* 1 GiB per rank */）なら ONEBODY をこの実行に限りフォールバックへ降格し、rank 0 INFO で理由を表示（`TraceKernelReportPlan` 時点で判定できるよう、判定関数は X->Def.NCisAjt と ncols から静的に計算）。演算子ループ: pair ごとに TraceMapExtract → 全所有状態をストリーミング（`g = Σ conj(z[k′])·amp·z[k]`）→ gbuf[ipair][*] に格納 → TraceMapFree（同時保持は 1 演算子分の O(N)、スペック §3.1）。

- [ ] **Step 2: 出力**

全 pair 終了後、状態順に: `X->Phys.eigen_num = n-1` をセットし、expec_cisajs.c が使うのと**同一のファイル名規約・行フォーマット**（expec_cisajs.c の fprintf 書式を関数化するか、書式文字列を 1 箇所の #define に共通化 — 実装時に expec_cisajs.c の該当 fprintf を確認し、**書式の二重定義を避ける**方を選ぶ）で per-state ファイル/パーシャル集約に書く。ExpecLocal は**オーケストレーション層が Mode 2 でも Enter/SetPartialSuffix 済み**の区間で走るため fopenMPI はランクローカル（Task 1 の配線で `expec_trace_owned_states` 呼び出しを ExpecLocalEnter/Leave で囲む — Task 1 Step 2 のコードに含める）。

- [ ] **Step 3: 単体テストケース**（GC: ランダム 3 状態パネルで gbuf の値が expec_cisajs_HubbardGC 直接実行と 1e-13 一致）→ RED→GREEN → 回帰 → コミット

---

### Task 5: 二体GFカーネル

**Files:**
- Modify: `src/expec_trace.c`
- Test: 単体テストに二体ケース追加

一体（Task 4）と同一構造: `NCisAjtCkuAlvDC × ncols` バッファ + 同じメモリゲート（一体・二体は別々に判定 — 片方だけ降格可）。演算子ごとに `Rearray_Interactions(...,2)` → TraceMapExtractTwoBody → ストリーミング。Rearray 失敗ペアは Mode 1 と同一の 0.0 行を書く（expec_cisajscktaltdc.c の該当 fprintf と同一書式）。出力書式は expec_cisajscktaltdc.c の per-state 書式と同一。単体テスト（GC 4 分岐 + 同一添字 + ゼロ結果）→ RED→GREEN → 回帰 → コミット。

---

### Task 6: ゴールデンテストとケイパビリティ最終確定

**Files:**
- Modify: `test/fulldiag_expecmode_equiv.sh`（Mode 2 の実カーネル検証を強化）
- Modify: `src/expec_trace.c`（表の最終値）

- [ ] **Step 1: equiv スクリプトに Mode 2 実カーネル検証を追加**

既存 4 ケース（Hubbard+NBodyG / SpinGC Gamma / 正準 Spin / SpinGC honeycomb 多体）は既に mode0/mode1/mode2 を回して全出力比較している — 3b では mode2 が実カーネルになるため**追加変更は INFO アサーションのみ**（Task 1 Step 4 で対応済み）。ここでは新規ケースを 1 つ追加: **Hubbard L=4 正準で一体+二体 GF を持つケースの mode2 出力が mode0 と 1e-8 一致**（正準系のゴールデン — ローカルでは登録のみ、実行は Task 9）。スペック §6「非対角項・密度-密度対角項・同一添字・ゼロ結果演算子」は greentwo.def の内容として盛り込む（サイト範囲は **Nsite 未満**を厳守 — 3a で範囲外サイトが既存 SIGFPE を踏んだ教訓）。

- [ ] **Step 2: ケイパビリティ表の最終確定**

単体テスト（GC）+ equiv（正準含む、Task 9 で実行）に合格した組だけ TRUE を残す。合格しなかった組は FALSE に戻し、理由をコメントで表に記録。

- [ ] **Step 3: 回帰 + コミット**

---

### Task 7: docs + 移行ノート

**Files:**
- Modify: `doc/en/.../CalcMod_file_en.rst`, `doc/ja/.../CalcMod_file_ja.rst`（ExpecMode 2 の実態化: 対象量・対象モデル・フォールバック規則・「var 検査は行わない」・INFO の読み方。3a で入れた「2 は 1 として動作」の記述を置換）
- Create: `docs/superpowers/specs/2026-07-12-phase3b-migration-note.md`（PR 転記用: ExpecMode 2 実装、利用指針「フェーズ3b 以降は通常 2 を推奨」[スペック §2]、var 列の意味変更）
- Modify: `test/manual/elpa_gpu_check.md`（Mode 2 検証項目 + ベンチ項目追加）

検証: INFO 文字列は実装から verbatim 引用（字下げ注記の流儀は 3a と同じ）。sphinx/rst2html 両言語レンダー確認。コミット。

---

### Task 8: ローカル最終回帰 + ガード確認

- build_noMPI 18/18 / build_mpi の unit（expec_trace_map_check 含む）green / `check_expec_local_calls` PASS（expec_trace.c スキャン対象で違反ゼロ）/ `sh -n` 各スクリプト。コミット（残作業があれば）。

---

### Task 9: clavius 実機検証 + ベンチマークゲート（コントローラ直接実行）

- [ ] rsync（リポジトリルートから）→ build_elpa 再構成・ビルド（+ expec_trace_map_check / equiv 用バイナリ）
- [ ] `expec_trace_map_check`（np=1）/ equiv np=2,3（全ケース: mode0 vs 1 vs 2）/ 既存 ELPA テスト全部（statepanel sweep, zero_owner, merge check, hubbard_chain, solver_keyword）
- [ ] **ベンチマークゲート**（スペック §6）: N=4900 全状態の一体+二体 GF 実時間を Mode 0/1/2 で比較（np=4）。目標: Mode 2 ≥ Mode 1（実測値と内訳を記録することが完了条件 — 「≥」未達なら損益分岐の分析を記録し、docs の利用指針を実測に合わせて修正）
- [ ] 可能なら L=10 Hubbard（N≈63504, 対角化 GPU）で Mode 0/1/2 の実時間 1 点（スペック §6 の外挿根拠。GPU 使用前に **nvidia-smi で空き確認**、他ユーザーのジョブと競合しないこと）
- [ ] 結果を `test/manual/elpa_gpu_check.md` に「Phase 3b validation」節として追記・コミット

## 完了条件（フェーズ3b）

- 既定（ExpecMode 0/1）全既存テスト無変更 PASS、ガード PASS（expec_trace.c 含む）
- GC 単体テスト（写像正当性・純粋性・一体・二体・エネルギー系）PASS（ローカル+clavius）
- equiv np=2/3 で mode0/1/2 全出力一致（1e-8）— 正準 Hubbard の一体+二体ゴールデン含む
- ベンチ記録（Mode 2 vs 1 vs 0、損益分岐コメント付き）
- ケイパビリティ表の最終値がテスト結果と一致（合格モデルのみ TRUE）
- docs/移行ノートが実装（INFO 文字列・var 列の扱い・フォールバック規則）と一致
