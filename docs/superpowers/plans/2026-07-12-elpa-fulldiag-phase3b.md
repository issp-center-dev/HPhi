# ELPA FullDiag フェーズ3b（ExpecMode 2 / トレースカーネル）実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `ExpecMode 2` を実キーワードにする — 演算子ごとに基底写像 (k→k′, 振幅) を 1 回前計算し、所有全状態を密ループでストリーミング評価するトレースカーネル（`src/expec_trace.c`）を、**一体GF・二体GF に限定して**モデル別の段階的有効化とフォールバック付きで実装する。

**Architecture:** 設計文書 `docs/superpowers/specs/2026-07-11-elpa-fulldiag-phase3-design.md`（v5.1、§3「Mode 2」が正）のフェーズ3b。フェーズ3a の Mode 1 基盤（状態パネル、ExpecLocal、パーシャル+マニフェスト出力）は完成済みで、3b は `phys_stateparallel()` 内のカーネル選択として載る。**唯一の真実源は呼び出しごとに 1 回構築される不変の `TraceExecutionPlan`**（静的ケイパビリティ ∧ 実行時メモリゲート）で、INFO 表示・カーネル実行・フォールバックループの 3 者が同じ plan を消費する。①スケルトン+plan+INFO（全量フォールバック=挙動は Mode 1 と同一）→ ②写像プローブアダプタ（mltply*Core.c 側）+家系別監査+純粋性テスト → ③一体GFカーネル → ④二体GFカーネル → ⑤ゴールデンテストでケイパビリティ確定 → ⑥docs → ⑦ローカル回帰 → ⑧clavius検証+ベンチ（正準系の早期チェックポイント込み）。

**Tech Stack:** C99, MPI（オーケストレーション層のみ）, 既存 mltply*Core 要素関数, フェーズ3a の green_output/ExpecLocal/状態パネル基盤, CMake/ctest。

## Global Constraints

- スペック §3「Mode 2: トレースカーネル」が正。**演算子代数を再実装しない** — 写像抽出は既存要素関数と同一コード経路を共有する**私設プローブアダプタ**（各要素関数の定義ファイル内に併設、公開 API 変更なし）経由で行う。
- **3b のカーネル対象は ONEBODY（expec_cisajs 相当）と TWOBODY（expec_cisajscktaltdc 相当）のみ。** エネルギー系（energy/var/doublon/num/sz = expec_energy_flct 一式）・S²（expec_totalspin）・NBodyG・AnomalousG は **3b では常にフォールバック**（Mode 1 経路）。理由: (a) var は固有ベクトル品質検査であり、固有値代用は「ExpecMode は速度のみを変える」保証に違反する（レビュー却下済みの旧案）。(b) expec_energy_flct は energy/var/対角量を 1 回の mltply で同時計算するため部分カーネル化は二重書き手を生む。エネルギー系のトレース化（H 全項写像）は 3c 拡張として移行ノートに明記。**この方針により expec_energy_flct はフォールバックループで常に実行され、v0→v1 の受け渡し規約（phys_distributed_local.c:70-81）は Mode 1 と完全に同一のまま**（v1 詰め替えの特殊処理は不要）。
- **保証**: `ExpecMode` は速度のみを変える。0/1/2 の全出力列（var 含む）は丸め誤差の範囲で一致（1e-8、bit 一致ではない）。各量の書き手は**常に一意**（plan のマスクが唯一の判定。ケイパビリティ表・メモリゲート・フォールバックループが独立に判断することを構造的に禁止する）。
- モデル別段階的有効化: ゴールデンテスト（Task 5）合格の（モデル×量×**演算子家系ブランチ全網羅**）のみケイパビリティ表で TRUE。**実装中は表は全 FALSE のまま**（単体テストはカーネル関数を直接呼ぶ — 本番ディスパッチを経由しない）。初期候補は Hubbard（正準+GC）と Spin/SpinGC の half。汎用スピン・spinless・Kondo は対象外。
- `src/expec_trace.c` は **MPI フリー**（mpi.h include 禁止、生 MPI・exitMPI 禁止、許可 wrapperMPI は SumMPI_dc/d/li/i, fopenMPI, childfopenMPI, stdoutMPI のみ）。作成と同一コミットで `test/check_expec_local_calls.sh` の `FILES` に追加。集団操作はすべて `src/phys_distributed.c`（恒久的にスキャン対象外）。プローブアダプタを併設する `mltply*Core.c` は現状どおりスキャン対象外（MPI 面なし — 凍結インベントリ §2 の記録どおり）。
- 要素関数プローブは **`X->Large.mode = M_CORR` のまま**行う（tmp_v0 書き込みが構造的に無効。M_MLTPLY=0, M_ENERGY=1, M_CORR=3, M_CALCSPEC=4, H_CORR=5）。H_CORR は使わない。**プローブの結合定数は常に tmp_V=1.0 を渡す**（GF には結合定数がなく、汎用経路も tmp_OneGreen=1.0 / Rearray の tmp_V を使う。これにより「戻り値==0 ⟺ 遷移消滅」が成立 — 振幅は ±1 または ±位相×1 のため。二体は Rearray が返す tmp_V(±1/共役符号) をそのまま渡し、`assert(tmp_V != 0.0)` を置く）。
- 既定（`ExpecMode 0/1`）の挙動は不変。build_noMPI（`cmake -DENABLE_MPI=OFF ..`）の `ctest -R "fulldiag|check_expec"` 18/18 基準。ELPA 実行系はローカルでは登録のみ、実行は Task 8（clavius, `~/HPhi-elpa`, conda env `hphi_elpa`, `LD_LIBRARY_PATH=$HOME/opt/elpa-2025.06-cuda/lib:$CONDA_PREFIX/lib`, **`export CUDA_MPS_PIPE_DIRECTORY=/tmp/nonexistent-mps-hphi` 必須**、rsync はリポジトリルートから）。**正準系（GetOffComp/list_1）の写像はローカル単体テスト不能のため、Task 3/4 完了直後にそれぞれ clavius 早期チェックポイント（コントローラ実行）を置く** — 最終ゲートまで正準系の検証を遅らせない。
- GC 系要素関数は list_1 間接参照なし（基底添字=ビットパターン）→ GC モデルの単体テストはローカル MPI ビルドで実行可能。
- コミットメッセージ末尾に本セッションの Co-Authored-By/Claude-Session トレーラ。

## File Structure

- Create: `src/expec_trace.c` — TraceExecutionPlan・ストリーミング・カーネル本体（MPI フリー）
- Create: `src/include/expec_trace.h` — 公開 API（Task 1 で確定、以後変更しない）
- Modify: `src/mltplyHubbardCore.c`, `src/mltplySpinCore.c` — 私設写像プローブアダプタ（各要素関数の直後に併設。公開ヘッダには載せず、`src/include/expec_trace_probe.h`（新規、src 内部用）に宣言）
- Modify: `src/phys.c`（ダウングレード INFO 撤去）, `src/phys_distributed.c`（plan 構築+カーネル呼び出し）, `src/phys_distributed_local.c`（plan マスク消費）
- Modify: `src/CMakeLists.txt`, `test/CMakeLists.txt`, `test/check_expec_local_calls.sh`
- Modify: `docs/superpowers/specs/2026-07-11-expec-call-inventory.md` — §2c として家系別プローブ監査表を追加
- Test: `test/unit/expec_trace_map_check.c`（GC 写像・純粋性・境界、ローカル実行可）, `test/fulldiag_expecmode_equiv.sh` 拡張
- Docs: CalcMod ja/en, `docs/superpowers/specs/2026-07-12-phase3b-migration-note.md`

---

### Task 1: TraceExecutionPlan + スケルトン + ディスパッチ配線（全量フォールバック）

**Files:**
- Create: `src/expec_trace.c`, `src/include/expec_trace.h`
- Modify: `src/phys.c:90-96`, `src/phys_distributed.c`, `src/phys_distributed_local.c`, `src/CMakeLists.txt`, `test/check_expec_local_calls.sh`
- Modify: `test/fulldiag_expecmode_equiv.sh`（INFO アサーション更新）

**Interfaces:**
- Produces（`src/include/expec_trace.h` — **この形が最終。以後のタスクはこのヘッダを変更しない**。写像抽出などカーネル内部の型・関数は Task 2 の `expec_trace_internal.h`（別ヘッダ）に置き、こちらへは追加しない。HPhi はヘッダを外部インストールしないため両者とも内部ヘッダだが、オーケストレーション境界の凍結として区別する）:
  ```c
  typedef enum {
    TRACE_Q_ONEBODY = 0,   /* expec_cisajs 相当 */
    TRACE_Q_TWOBODY,       /* expec_cisajscktaltdc 相当 */
    TRACE_Q_NQUANT
  } TraceQuantity;
  typedef struct {
    /* kernel[q]==1: トレースカーネルが担当。0: フォールバック（Mode 1 経路）。
       構築後は不変。INFO 表示・カーネル・フォールバックループの全員がこの
       同一インスタンスを消費する（他の判定源を持たない） */
    int kernel[TRACE_Q_NQUANT];
    /* demoted_memory[q]==1: 静的には対応モデルだがメモリゲートで降格した */
    int demoted_memory[TRACE_Q_NQUANT];
    long int nc_uniform;   /* 構築に使った一様ブロック幅 NC=ceil(neig/nproc)
                              （rank 局所の所有数ではない。ゼロ所有判定は
                              je<jb で行い、この値は使わない） */
  } TraceExecutionPlan;
  /* 静的ケイパビリティ表 ∧ 実行時メモリゲート（チェック付きサイズ計算）で
     plan を 1 回構築する。nc_uniform は全ランク同値の NC=ceil(neig/nproc) を
     呼び出し側が渡す（rank 局所 ncols は渡してはならない — plan の全ランク
     一致が崩れる）。ExpecMode!=2 なら全量 fallback の plan を返す */
  void TraceBuildPlan(const struct BindStruct *X, long int nc_uniform,
                      size_t gbuf_max_bytes, TraceExecutionPlan *plan);
  /* gbuf_max_bytes は MPI オーケストレーション層（phys_distributed.c）が
     rank 0 で環境変数を解析し MPI_Bcast した値を渡す — 環境変数がノード間で
     不一致でも plan は全ランク一致（expec_trace.c は MPI フリーのまま） */
  /* rank 0 用: plan の内容を量ごとに INFO 表示（降格理由込み）。
     エネルギー系/S2/NBodyG/AnomalousG が常時フォールバックである旨の固定行も出す */
  void TraceReportPlan(const TraceExecutionPlan *plan, FILE *fp);
  /* plan->kernel[q]==1 の量だけを所有状態一括評価し GF 出力へ書く。
     戻り値 0/-1（ローカル rc）。書き込み開始前に量単位で完結性を保証:
     写像抽出やバッファ確保の失敗は「その量を書かずに rc=-1」（部分出力なし） */
  int expec_trace_owned_states(struct BindStruct *X, const TraceExecutionPlan *plan,
                               const double complex *panel,
                               long int jb, long int je, long int NN);
  ```
- Consumes: 3a の `phys_stateparallel()`（phys_distributed.c:43-187）と `phys_stateparallel_local_loop()`（phys_distributed_local.c:67-92）。

- [ ] **Step 1: plan・表・スケルトン**

`src/expec_trace.c` に静的表（**全 FALSE で出荷。Task 5 まで誰も TRUE にしない**）:

```c
typedef struct { int calc_model; int flg_general_spin; int q[TRACE_Q_NQUANT]; } TraceCap;
static const TraceCap kTraceCap[] = {
  /* Rows are flipped to 1 ONLY by plan Task 5, after the golden
     cross-checks for EVERY reachable operator-family branch of that
     (model, quantity) pass. See the branch-coverage table in
     docs/superpowers/specs/2026-07-11-expec-call-inventory.md §2c. */
  { Hubbard,   0, {0, 0} },
  { HubbardGC, 0, {0, 0} },
  { Spin,      0, {0, 0} },   /* half のみ; iFlgGeneralSpin==1 は行に一致させない */
  { SpinGC,    0, {0, 0} },
};
```

`TraceBuildPlan`: 表を引き、TRUE の量についてのみメモリゲートを評価する。**判定は rank 依存の ncols ではなく、全ランクで同一の一様ブロック幅 `NC = ceil(neig/nproc)` を用いる**（各 rank の ncols ≤ NC、かつ NC は通信なしで全ランク同値 → plan が構造的に全ランク一致し、rank 間での kernel/fallback 混在実行が起こらない。API の引数・フィールドは `nc_uniform` と命名済み — 上記 Interfaces 参照）。サイズ判定と確保は**同一のチェック付きヘルパ**を共用する:

```c
/* 0 を返したら「収まらない/表現不能」。plan 構築（判定）と Task 3/4 の
   malloc（確保サイズ計算）の両方がこの 1 関数を使う — 2 判定の乖離を構造的に禁止 */
static size_t TraceGbufBytes(long int nops, long int nc_uniform, size_t max_bytes) {
  size_t sn, sc;
  if (nops <= 0 || nc_uniform <= 0) return 0;
  /* long int が size_t より広い環境での切り詰めを先に排除（uintmax_t 経由の
     表現可能性チェック — 64bit 前提にしない） */
  if ((uintmax_t)nops > (uintmax_t)SIZE_MAX ||
      (uintmax_t)nc_uniform > (uintmax_t)SIZE_MAX) return 0;
  sn = (size_t)nops; sc = (size_t)nc_uniform;
  if (sc > SIZE_MAX / sizeof(double complex)) return 0;      /* ncols*16 が overflow */
  if (sn > SIZE_MAX / (sc * sizeof(double complex))) return 0; /* nops*(ncols*16) が overflow */
  if (sn * sc * sizeof(double complex) > max_bytes) return 0; /* キャップ超過 */
  return sn * sc * sizeof(double complex);
}
/* HPHI_TRACE_BUF_MAX_MB の MiB→byte 変換も同様にチェック付きで行う
   （value_mb > SIZE_MAX >> 20 なら既定値へフォールバックし警告） */
```
`TraceGbufBytes(...)==0` なら demoted_memory[q]=1, kernel[q]=0。`TraceGbufMaxBytes()`: 環境変数 `HPHI_TRACE_BUF_MAX_MB`（1..1048576 の整数のみ受理、不正値は既定にフォールバックして stderr に 1 行警告）×2^20、未設定は既定 1024 MiB。**このゲートは量ごと・ランクごとの結果バッファ 1 本のキャップであり、プロセス総メモリの上限ではない**（コメントで明記。写像 O(N) 1 本と panel は別勘定）。

`TraceReportPlan` の出力（equiv テストがこの全行を検証する — 量ごと 1 行+固定行 1 行）:
```
  INFO: ExpecMode 2: one-body Green functions use the trace kernel.
  INFO: ExpecMode 2: two-body Green functions use the ExpecMode-1 fallback (unsupported model).
  INFO: ExpecMode 2: two-body Green functions use the ExpecMode-1 fallback (result buffer would exceed HPHI_TRACE_BUF_MAX_MB).
  INFO: ExpecMode 2: energy/fluctuation, S2, NBodyG, and AnomalousG always use the ExpecMode-1 path in this version.
```
（2 行目と 3 行目は排他 — 理由テキストは "unsupported model" / "result buffer would exceed HPHI_TRACE_BUF_MAX_MB" の 2 種。）

`expec_trace_owned_states` はこの段階では担当量なし（plan が全 FALSE）で即 return 0。

**開発用強制フック（Task 5 で削除）**: `TraceBuildPlan` は環境変数 `HPHI_TRACE_FORCE` を読む。値は**量名のカンマ区切りリスト**（`onebody` / `twobody`。例: `HPHI_TRACE_FORCE=onebody,twobody`）。列挙された量だけをケイパビリティ表の値に関係なく kernel[q]=1 にする（メモリゲートは通常どおり適用）。**列挙されない量は絶対に強制されない**ため、未実装カーネルが選択されて出力が欠落する事故は構造的に起こらない（Task 3 のチェックポイントは `=onebody`、Task 4 は `=onebody,twobody` を使う）。認識できないトークンは stderr 警告のうえ無視。実装箇所には `/* development-only hook: plan Task 5 REMOVES this */` コメントを付す。

- [ ] **Step 2: ディスパッチ配線（plan が唯一の判定源）**

`src/phys.c:91-92` のダウングレード INFO 2 行を撤去。`src/phys_distributed.c` の `phys_stateparallel()`、ローカルループ呼び出し部を:

```c
  TraceExecutionPlan tplan;
  long int nc_uniform = (NN + (long int)nproc - 1) / (long int)nproc; /* 全ランク同値 */
  unsigned long gbuf_max = 0;
  if (myrank == 0) gbuf_max = (unsigned long)TraceGbufMaxBytesFromEnv(); /* env は rank 0 のみ解析 */
  MPI_Bcast(&gbuf_max, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  TraceBuildPlan(X, nc_uniform, (size_t)gbuf_max, &tplan); /* ExpecMode!=2 なら全量 fallback */
  if (X->Def.iExpecMode == EXPECMODE_TRACE && myrank == 0)
    TraceReportPlan(&tplan, stdoutMPI);
  ExpecLocalEnter();
  GreenOutputSetPartialSuffix(myrank);
  rc_local = expec_trace_owned_states(X, &tplan, panel, jb, je, NN);
  if (rc_local == 0)
    rc_local = phys_stateparallel_local_loop(X, panel, jb, je, NN, &tplan);
  GreenOutputClearPartialSuffix();
  ExpecLocalLeave();
```

**注意（3a からの構造変更を最小化）**: 3a では Enter/SetPartialSuffix はローカルループ内部にあった。plan 消費とカーネルを同じ ExpecLocal セッションに入れるため、Enter/Set/Clear/Leave を**オーケストレーション層へ引き上げ**、`phys_stateparallel_local_loop` からは除去する（同関数のシグネチャに `const TraceExecutionPlan *plan` を追加し、per-state ループで `plan->kernel[TRACE_Q_ONEBODY]` なら `expec_cisajs` 呼び出しをスキップ、`plan->kernel[TRACE_Q_TWOBODY]` なら `expec_cisajscktaltdc` をスキップ。**それ以外の evaluator と all_* 代入群は無条件に従来どおり** — エネルギー系は常にフォールバックなので v0/v1 受け渡しは Mode 1 と同一）。`phys_distributed_local.c` は MPI フリーのまま（plan は値渡しの読み取り専用構造体）。ガードの許可リストに影響なし。3a の既存呼び出し（Mode 1）も同じ経路を通る（plan 全 FALSE）ため挙動不変。**単一出口不変条件**: `ExpecLocalEnter()` 成功後は、カーネル/フォールバックのいずれが失敗しても `GreenOutputClearPartialSuffix()` と `ExpecLocalLeave()` を**ちょうど 1 回ずつ**通ってから return する（上記コード形を崩す早期 return を将来追加することを禁ずるコメントを添える）。ランク局所の失敗の集約は Mode 1 と同一 — rc_local が既存の単一ランデブー `MPI_Allreduce(MIN)` に入り、失敗時は Gatherv/Merge に到達しない。ゼロ所有ランク（ncols==0）はカーネルでも即 return 0（Mode 1 のループ不実行と同型）。

- [ ] **Step 3: ガード FILES 追加 + ビルド + 回帰**

`test/check_expec_local_calls.sh` の `FILES` に `src/expec_trace.c` を追加。`src/CMakeLists.txt` に expec_trace.c を追加（全ビルドでコンパイル。MPI 非依存）。

Run: `cd build_noMPI && cmake .. && make HPhi -j8 && ctest -R "fulldiag|check_expec"`
Expected: 18/18 PASS（Mode 2 実挙動 = 全量フォールバック = Mode 1 と同一、INFO のみ変化）。build_mpi でも `make HPhi` + `ctest -R "green_partial|check_expec"` green。

- [ ] **Step 4: equiv スクリプトの Mode-2 INFO アサーション更新**

3a の "ExpecMode 2 kernels are not available in this build" grep を撤去し、**plan の全行検証**に置換: mode2 実行ログに (a) `"ExpecMode 2: one-body Green functions use the"` 行、(b) `"ExpecMode 2: two-body Green functions use the"` 行、(c) `"always use the ExpecMode-1 path"` 行が**すべて**存在すること（この段階では (a)(b) とも fallback 理由付き）。出力一致検証は従来どおり。`sh -n` 確認。

- [ ] **Step 5: コミット** `git commit -m "Introduce the TraceExecutionPlan dispatch for ExpecMode 2"`

---

### Task 2: 写像プローブアダプタ（mltply*Core 併設）+ 家系別監査 + 純粋性テスト

**Files:**
- Modify: `src/mltplyHubbardCore.c`, `src/mltplySpinCore.c`（プローブアダプタ併設）
- Create: `src/include/expec_trace_probe.h`（src 内部宣言）
- Modify: `src/expec_trace.c`（TraceMap 型と抽出ドライバ）
- Modify: `docs/superpowers/specs/2026-07-11-expec-call-inventory.md`（§2c 家系別監査表）
- Create: `test/unit/expec_trace_map_check.c` / Modify: `test/CMakeLists.txt`

**Interfaces:**
- Produces（**`src/include/expec_trace_internal.h`** — Task 1 で確定済みの `expec_trace.h`（オーケストレーション API のみ、以後不変）とは別の src 内部ヘッダ。単体テストはこちらを include してカーネル内部を直接呼ぶ。HPhi はヘッダをインストールしないが、公開面の規律として区別する）:
  ```c
  typedef struct {
    long int n;            /* = X->Check.idim_max（ヒルベルト次元。正準/GC とも。
                              kprime の値域は [-1, n-1] — 抽出後に範囲アサート） */
    long int *kprime;      /* [n] 0-based 行き先; 遷移消滅は -1 */
    double complex *amp;   /* [n] 振幅（kprime>=0 のときのみ有意） */
    int is_diagonal;
  } TraceMap;
  int TraceMapExtractOneBody(struct BindStruct *X, int ipair, TraceMap *map);
  int TraceMapExtractTwoBody(struct BindStruct *X, int ipair, TraceMap *map);
  void TraceMapFree(TraceMap *map);
  ```
  （probe ベクトルは抽出側の内部実装詳細 — アダプタが値そのものを返すため ones 配列は不要。`expec_trace_probe.h` のアダプタ宣言もこの内部ヘッダに統合してよい[実装時にファイル数を最小化する側を選ぶ]。）
- Produces（`src/include/expec_trace_probe.h`、アダプタ群 — **既存要素関数の公開シグネチャは一切変更しない**）:
  ```c
  /* 各既存要素関数の定義ファイル内、当該関数の直後に併設。命名規約:
     <元関数名>_TraceProbe。戻り値: 1=遷移あり(*kprime_out,*amp_out 有効) /
     0=消滅。実装は元関数と同一のビット判定・GetOffComp 呼び出し列を共有する
     （元関数の本体を「写像計算部」と「ベクトル適用部」に分ける static 関数
     抽出リファクタで共有し、コピーは作らない — 元関数の数値挙動が変わらない
     ことは既存回帰 18/18 が担保） */
  int CisAjt_TraceProbe(long unsigned int j, struct BindStruct *X,
      long unsigned int is1_spin, long unsigned int is2_spin,
      long unsigned int sum_spin, long unsigned int diff_spin,
      long int *kprime_out, double complex *amp_out);
  /* 同様に: GC_CisAjt_TraceProbe, GC_CisAis_TraceProbe,
     child_Spin_CisAis_TraceProbe, child_SpinGC_CisAis_TraceProbe,
     child_SpinGC_CisAit_TraceProbe,
     CisAisCisAis_element_TraceProbe, CisAisCjtAku_element_TraceProbe,
     CisAjtCkuAku_element_TraceProbe, CisAjtCkuAlv_element_TraceProbe,
     （GC_ 版 4 種）, GC_CisAisCisAis_spin_element_TraceProbe,
     GC_CisAisCitAiu_spin_element_TraceProbe,
     GC_CisAitCiuAiu_spin_element_TraceProbe,
     GC_CisAitCiuAiv_spin_element_TraceProbe,
     （Spin-half 正準二体の到達家系 — Step 0 の監査で確定した分） */
  ```

- [ ] **Step 0: 家系別プローブ監査（成果物 = インベントリ §2c 表）**

対象 4 モデル×2 量について、汎用経路（expec_cisajs.c / expec_cisajscktaltdc.c の LOCAL 分岐）が**到達し得る要素関数ブランチを全数列挙**し、各行に: (i) 元関数名と定義位置 (ii) kprime の取得法（out-param / 対角恒等 / GetOffComp 内部）(iii) 振幅の構成（符号×tmp_V、共役の有無、Rearray 前処理で吸収済みの係数）(iv) M_CORR で tmp_v0 書き込みが無効であることのソース根拠（行番号）(v) プローブアダプタ名 (vi) **当該経路の到達ヘルパ（GetOffComp/SgnBit 等含む）が書き込む可能性のある X のフィールド・グローバル配列の全列挙**（純粋性テストの snapshot 対象リストはこの列から導出する — 場当たりで選ばない）、を記録する。**Spin 正準（half）の二体家系はフェーズ3a の探索で未踏なので、ここで expec_cisajscktalt_SpinHalf（expec_cisajscktaltdc.c:961 以降のディスパッチ先）を読み、到達家系を列挙して表に含める**。表にない家系が汎用経路に存在した場合はその（モデル×量）を対象から外す（表が根拠）。この表は Task 5 のケイパビリティ TRUE 化の前提条件リストになる。

- [ ] **Step 1: アダプタ実装（一体: Hubbard 正準/GC, Spin-half, SpinGC-half）**

方式は関数抽出リファクタ: 例として `CisAjt`（mltplyHubbardCore.c:354-396）を
```c
static int CisAjt_map(long unsigned int j, struct BindStruct *X,
    long unsigned int is1_spin, long unsigned int is2_spin,
    long unsigned int sum_spin, long unsigned int diff_spin,
    long unsigned int *off_out, double complex *sgn_out);  /* 既存本体の写像計算部 */
```
に抽出し、既存 `CisAjt` は `CisAjt_map` を呼んでから従来どおり `tmp_v0`/`dam_pr` を処理、`CisAjt_TraceProbe` も `CisAjt_map` を呼ぶだけ、と 3 層にする（**同一コード経路の共有 = 代数の再実装なし**。対角系（CisAis 系）はプローブが `*kprime_out=j-1` を返し振幅=ビット判定結果）。`expec_trace.c` の `TraceMapExtractOneBody` は expec_cisajs.c の LOCAL 分岐と同一の前処理（`general_hopp_GetInfo` → snapshot）を行い、k=1..n でアダプタを呼んで詰める。

- [ ] **Step 2: アダプタ実装（二体: Hubbard 正準/GC 4 分岐, SpinGC-half 4 分岐, Spin-half 正準の到達家系）**

`Rearray_Interactions(ipair,...,2)` → `general_int_GetInfo` →（監査表の家系ごとに）`*_element_TraceProbe`。`*_element` 系も同じ関数抽出方式（写像計算部 `*_element_map` を共有）。プローブの結合定数規約を精密化: **tmp_V は 1.0 で初期化し、`Rearray_Interactions` が並べ替えに伴い変換した値（±1/共役因子。正規の GF ペアでは常に非ゼロ — `assert`）をそのまま共有 `*_map` ルーチンへ渡す**。共役・符号反転を伴う Rearray 分岐は単体テストで明示的に踏む。Rearray 非ゼロ返却ペアは `map->n=0` で返し、呼び出し側（Task 4）が Mode 1 と同一の 0.0 行フォールバック（**Rearray 非ゼロ = 「規約外ペア → 0.0 行を書く」が Mode 1 の実挙動であることを expec_cisajscktaltdc.c:1007 付近で確認し、エラー流用でないことを監査表に記録**）。

- [ ] **Step 3: 純粋性・写像正当性テスト（RED→GREEN）** — `test/unit/expec_trace_map_check.c`

シリアル・MPI 不要・**GC モデルのみ**（HubbardGC L=4: n=256 / SpinGC-half L=4: n=16。green_partial_merge_check.c のスタブ流儀）。検査:
1. **写像正当性**: ランダム複素ベクトル z で `Σ_{k:kprime≥0} conj(z[kprime])·amp·z[k]` が、同じ要素関数（元関数、M_CORR、vec=z 直接）の `Σ dam_pr` と 1e-13 一致。一体: 対角/非対角/ゼロ結果（常時消滅の組）各 1。二体: 4 分岐各 1 + 同一添字 + ゼロ結果。
2. **純粋性（スペック §3.2b）**: 同一演算子で抽出 2 回 → kprime/amp が bit 一致。抽出前後の snapshot 対象は**Step 0 監査表 §2c の (vi) 書き込み集合列に記録された全フィールド・全配列**（実装時に確定した名前付きリストをテスト内に列挙する — 「等」で省略しない。memcmp 全域比較はパディングで無効なので使わない）。
3. **境界**: 演算子 0 個（NCisAjt=0）で抽出ドライバが何もしないこと。

Run: `cd build_mpi && make expec_trace_map_check && ./test/expec_trace_map_check` → RED（未実装）→ 実装 → GREEN。build_noMPI にも登録・実行。

- [ ] **Step 4: 回帰（アダプタ抽出リファクタが既存経路を壊していないこと = 18/18）+ コミット**

---

### Task 3: 一体GFカーネル + clavius 早期チェックポイント（正準）

**Files:**
- Modify: `src/expec_trace.c`
- Test: `test/unit/expec_trace_map_check.c` にストリーミング+バッファケース追加

- [ ] **Step 1: ストリーミングとバッファ**

`expec_trace_owned_states` の ONEBODY 部: `gbuf`（**確保バイト数は `TraceGbufBytes(nops, plan->nc_uniform, TraceGbufMaxBytes())` の戻り値のみを使う** — サイズ式の再記述禁止。0 が返る事態は plan 構築時に排除済みだが、0 なら防御的に rc=-1。malloc 失敗は**書き込み前なので**その量を rc=-1 で報告 — 部分出力なし）。演算子外側ループ: pair → `TraceMapExtractOneBody` → 全所有状態ストリーミング → gbuf → `TraceMapFree`（写像の同時保持は 1 本、スペック §3.1）。**全 pair 完了後に**出力フェーズ: 状態順に `X->Phys.eigen_num = n-1` を設定し、expec_cisajs.c と同一のファイル名規約・行書式で per-state ファイル/パーシャル集約へ書く（書式文字列は expec_cisajs.c の該当 fprintf と共通の #define へ抽出し二重定義を避ける。eigen_num の設定は Mode 1 の per-state 慣行と同じで、後続フォールバックループが状態ごとに再設定するため干渉しない）。書き込み中の失敗は Mode 1 の書き込み失敗と同じ扱い（sticky manifest エラー → 集団 rc=-1。**フォールバックへの再試行はしない** — 二重出力防止）。**量またぎの原子性は保証しない**（一体を書き終えた後に二体の準備で失敗した場合、一体の part は残るが、集団 rc=-1 により Merge は公開せず（マニフェスト規則）実行全体が失敗として終わる — Mode 1 のループ途中失敗と同じ回復モデル[再実行]。これは意図した仕様として docs に記載不要[内部挙動]、コード内コメントに記す）。

- [ ] **Step 2: 単体テストケース**（GC: ランダム 3 状態パネルで gbuf の中身が expec_cisajs_HubbardGC / expec_cisajs_SpinGCHalf の直接実行と 1e-13 一致。メモリゲート境界: `HPHI_TRACE_BUF_MAX_MB=1` で nops×ncols がゲートを跨ぐ 2 ケース — 降格した plan では kernel[q]==0 になること）→ RED→GREEN → 回帰 → コミット

- [ ] **Step 3（clavius 早期チェックポイント — コントローラ実行）**: rsync → build_elpa 再構成・ビルド → **一時的に**単体テストレベルで正準 Hubbard の写像正当性を検証: `expec_trace_map_check` に正準ケースを追加するのではなく、equiv の Hubbard ケースを `ExpecMode 2` + ケイパビリティ強制 ON（開発用フック `HPHI_TRACE_FORCE=onebody` — Task 1 で定義済みのセマンティクス）で np=2 実行し、mode0 と比較。**ログで一体が "use the trace kernel" 行になっていることを必ず確認**（一体のみ強制有効の状態）。不一致ならここで修正してから Task 4 へ進む。

---

### Task 4: 二体GFカーネル + clavius 早期チェックポイント（正準）

Task 3 と同一構造: `nops = X->Def.NCisAjtCkuAlvDC`、別バッファ・別ゲート判定（片方だけ降格可）。Rearray 失敗ペアは Mode 1 と同一の 0.0 行（expec_cisajscktaltdc.c の該当 fprintf と同一書式、GreenOutputWriteIndexPrefix 込み）。単体テスト（GC 4 分岐+同一添字+ゼロ結果+ゲート境界）→ RED→GREEN → 回帰 → コミット → **clavius 早期チェックポイント**（`HPHI_TRACE_FORCE=onebody,twobody` で equiv の SpinGC honeycomb 多体ケース[green6 含む — NBodyG 系がフォールバックであることも同時に確認できる]と Hubbard ケースを np=2/3 実行、mode0 比較。**ログで一体・二体の両方が "use the trace kernel" 行になっていることを必ず確認**する — フック値の打ち間違いでフォールバック経路だけを検証してしまう事故の防止）。

---

### Task 5: ゴールデンテストとケイパビリティ確定

**Files:**
- Modify: `test/fulldiag_expecmode_equiv.sh`, `src/expec_trace.c`（表の最終値のみ）

- [ ] **Step 1: equiv に正準 Hubbard 一体+二体のゴールデンケースを追加**（L=4、greentwo.def はスペック §6 の 4 種 — 非対角・密度-密度対角・同一添字・ゼロ結果 — を**サイト範囲 < Nsite 厳守**で含める。既存 4 ケースはそのまま）。**GetOffComp の分岐網羅**: 上向き/下向きホップ・順/逆順ペア・境界サイト対を greenone.def に含め、1 ケース内で正準写像の主要分岐を踏む。
- [ ] **Step 2: 表の TRUE 化**（監査表 §2c の全到達家系にアダプタがあり、単体テスト+早期チェックポイントに合格した (モデル×量) のみ。**各 TRUE 行のコメントに、根拠となる §2c 監査行の範囲と検証テスト名（expec_trace_map_check のケース名 / equiv ケース名）を機械可読に記録**する — 後日ディスパッチ分岐が変わったとき TRUE 行の妥当性を追跡できるようにする。合格しなかった組は FALSE のまま理由をコメントで記録）。`HPHI_TRACE_FORCE` フックはこの時点で**削除**（開発専用のため出荷しない — Task 1 で導入時に「Task 5 で削除」とコメントしておく）。
- [ ] **Step 3: equiv の INFO アサーションを最終 plan（kernel 行）に更新** → 回帰 → コミット

---

### Task 6: docs + 移行ノート

- CalcMod ja/en: ExpecMode 2 の実態化（対象量 = 一体・二体 GF、対象モデル、フォールバック規則、メモリゲートと `HPHI_TRACE_BUF_MAX_MB`、INFO の読み方。「2 は 1 として動作」の 3a 記述を置換。**エネルギー系・S² のトレース化は将来拡張**である旨）。INFO 文字列は実装から verbatim（字下げ注記は 3a の流儀）。
- Create: `docs/superpowers/specs/2026-07-12-phase3b-migration-note.md`（PR 転記用: ExpecMode 2 実装、利用指針は**ベンチ結果を見てから確定**[スペック §2 は「3b 以降は通常 2 を推奨」だが、Task 8 の実測が Mode 1 未満なら推奨文言を実測に合わせる]、var 列は Mode 2 でも従来どおり計算される[エネルギー系フォールバック]こと）。
- `test/manual/elpa_gpu_check.md` に Mode 2 検証項目+ベンチ項目追加。
- 両言語 rst レンダー確認 → コミット。

---

### Task 7: ローカル最終回帰

- build_noMPI 18/18 / build_mpi unit（expec_trace_map_check 含む）green / ガード PASS（expec_trace.c 対象、違反ゼロ）/ `sh -n` 各スクリプト。必要なら修正コミット。

---

### Task 8: clavius 実機検証 + ベンチマークゲート（コントローラ直接実行）

- [ ] rsync → 再構成・ビルド → `expec_trace_map_check`（np=1）/ equiv np=2,3（全ケース mode0/1/2）/ 既存 ELPA テスト全部
- [ ] **ベンチマークゲート**（スペック §6）: N=4900 全状態の一体+二体 GF 実時間を Mode 0/1/2 で比較（np=4）。目標 Mode 2 ≥ Mode 1 — **実測値と内訳（写像抽出時間 vs ストリーミング時間 vs 出力時間）の記録が完了条件**。未達なら損益分岐の分析を記録し Task 6 の利用指針を修正
- [ ] 可能なら L=10 Hubbard（N≈63504、対角化 GPU — **nvidia-smi で空き確認、他ユーザーと競合しない**）で Mode 0/1/2 実時間 1 点
- [ ] 結果を `test/manual/elpa_gpu_check.md` に「Phase 3b validation」節として追記・コミット

## 完了条件（フェーズ3b）

- 既定（ExpecMode 0/1）全既存テスト無変更 PASS、ガード PASS（expec_trace.c 含む）
- 家系別監査表（§2c）が存在し、TRUE 化された全 (モデル×量) の到達家系を網羅
- GC 単体テスト（写像正当性・純粋性・一体・二体・ゲート境界・境界条件）PASS（ローカル+clavius）
- equiv np=2/3 で mode0/1/2 全出力一致（1e-8、var 列含む）— 正準 Hubbard の一体+二体ゴールデン含む
- ベンチ記録（Mode 2 vs 1 vs 0、内訳と損益分岐コメント付き）と docs 利用指針の整合
- ケイパビリティ表の最終値がテスト結果と一致、`HPHI_TRACE_FORCE` フックが削除済み
- docs/移行ノートが実装（INFO 文字列・メモリゲート・フォールバック規則・var 列が従来どおりであること）と一致
