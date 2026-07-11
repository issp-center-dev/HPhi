# ELPA FullDiag フェーズ3a（ExpecMode / Mode 1 状態タスク並列）実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** FullDiag の物理量計算に `ExpecMode`（0=従来/1=状態タスク並列/2=3aでは1と同動作）を導入し、分散固有ベクトルを状態パネルに再分散して各ランクが担当状態の `expec_*` を独立実行する — 物理量計算の実時間を ~P 倍化し、状態ごとの集団転送を排除する。

**Architecture:** 設計文書 `docs/superpowers/specs/2026-07-11-elpa-fulldiag-phase3-design.md`（v5.1、収束済み。判断に迷ったらこれが正）のフェーズ3a。①コールインベントリ凍結 → ②`ExpecMode` キーワード → ③`ExpecLocal` フック（wrapperMPI/FileIO） → ④green_output のパーシャル出力＋マニフェスト → ⑤状態パネル再分散 → ⑥Mode 1 ドライバ＋Mode 0 の S²/Sz 統一 → ⑦モード等価性テスト → ⑧docs。

**Tech Stack:** C99, MPI, ScaLAPACK（`pzgemr2d_`）, 既存 expec/green_output 基盤, CMake/ctest。

## Global Constraints

- スペック §2/§3 が正。`ExpecMode` 有効条件: `iCalcType==FullDiag && iSolver∈{SOLVER_SCALAPACK,SOLVER_ELPA}`、`nproc==1` は INFO を出して 0 に降格（文言: "ExpecMode reverts to 0 for a single process (results are identical)."）。3a では `ExpecMode 2` は INFO を出して 1 と同動作。
- **保証**: `ExpecMode` は速度のみを変える。0/1/2 の結果は丸め誤差の範囲で一致（bit 一致ではない）。
- `ExpecLocal` 中の禁止事項: `exitMPI` 呼び出し禁止、生 `MPI_*` 到達禁止（nbody/anomalous の partner≠myrank 分岐には防御ガード `if (iExpecLocal) return エラー` を入れる）。ON/OFF は `ExpecLocalEnter()/ExpecLocalLeave()` のみ（入れ子不可アサート、`phys()` 出口で OFF アサート）。
- 集約 Green 出力はマニフェスト方式（spec §3 の {attempted, opened, open_error, bytes, closed_ok, part_path, final_path}。サイズでの成否推測禁止、連結・公開・削除はグローバル rc 成功時のみ、当該実行のマニフェストのみ読む、書き込み前に part を unlink）。
- 状態所有はフェーズ2 と同一の連続ブロック: `NC=ceil(N/P)`, `first_state(r)=r*NC+1`(1-based), `ncols_local(r)=max(0, min((r+1)*NC,N)-r*NC)`。ゼロ所有ランクは 1 要素ダミー＋有効ディスクリプタで集団参加（フェーズ2 np=50 実機検証済みパターン）。
- **Mode 0 も変更あり（唯一の例外）**: 分散 Mode 0 の状態ループに `ExpecLocalEnter/Leave` で囲んだ `expec_totalspin` を追加し、S²/Sz のゼロ埋めと stdout の短縮形式を廃止（シリアル形式に統一）。意図的な挙動修正としてテスト・docs・PR 移行ノートに反映。
- 既定（`ExpecMode 0` かつ非分散）の全既存テストは無変更 PASS。ローカル回帰は `build_noMPI`（`cmake -DENABLE_MPI=OFF ..`）の `ctest -R fulldiag` 17/17 基準（Mode 0 分散の挙動修正は ELPA テスト側で検証）。
- コミットメッセージ末尾に本セッションの Co-Authored-By/Claude-Session トレーラ。
- ローカル（macOS, ELPAなし）では ELPA 実行系テストは登録のみ。実行はタスク9（clavius, `~/HPhi-elpa`, conda env `hphi_elpa`, `LD_LIBRARY_PATH=$HOME/opt/elpa-2025.06-cuda/lib:$CONDA_PREFIX/lib`）。**rsync は必ずリポジトリルートから**（フェーズ2の教訓）。

---

### Task 1: ExpecLocal コールインベントリ凍結 + ガードスクリプト

**Files:**
- Create: `docs/superpowers/specs/2026-07-11-expec-call-inventory.md`
- Create: `test/check_expec_local_calls.sh`（+x）
- Modify: `test/CMakeLists.txt`（`add_hphi_test_with_srcdir(check_expec_local_calls)` を validation テスト群の近くに追加）

**Interfaces:**
- Produces: 凍結済み許可マトリクス（Task 3 のフック実装対象リスト）と、防御ガードが必要な生 MPI 分岐の完全な一覧（Task 3 で使用）。

- [ ] **Step 1: 到達コールグラフのインベントリを機械的に生成**

対象ファイル集合（expec 到達層）: `expec_energy_flct.c expec_cisajs.c expec_cisajscktaltdc.c expec_totalspin.c nbody_correlation.c anomalous_pair.c` と、それらが呼ぶ要素関数層 `mltplyHubbardCore.c mltplySpinCore.c mltplyMPIHubbardCore.c mltplyMPISpinCore.c`（実在するファイル名は `ls src/mltply*` で確認して調整）。各ファイルで:

```bash
cd /Users/k-yoshimi/Dropbox/CLionProjects/HPhi-box/HPhi
# 注意: gcc -fpreprocessed は Apple clang に存在しない（このマシンで確認済み）。
# コメント除去は Python で行い、失敗・空出力は必ずエラーにする（黙って空の
# インベントリで PASS する事故を構造的に禁止）。
python3 - "$f" <<'PY' などではなく、次の共有ストリッパを使う:
cat > /tmp/strip_c_comments.py <<'PY'
import re, sys
src = open(sys.argv[1]).read()
out = re.sub(r'/\*.*?\*/', lambda m: '\n'*m.group(0).count('\n'), src, flags=re.S)
out = re.sub(r'//[^\n]*', '', out)
assert out.strip(), f"stripper produced empty output for {sys.argv[1]}"
sys.stdout.write(out)
PY
for f in src/expec_energy_flct.c src/expec_cisajs.c src/expec_cisajscktaltdc.c \
         src/expec_totalspin.c src/nbody_correlation.c src/anomalous_pair.c; do
  echo "== $f =="
  python3 /tmp/strip_c_comments.py "$f" | \
    grep -noE "(SumMPI_[a-z]+|MaxMPI_[a-z]+|BcastMPI_[a-z]+|BarrierMPI|NormMPI_dc|VecProdMPI|MPI_[A-Za-z_]+|exitMPI|fopenMPI|childfopenMPI)\(" | sort | uniq -c
done
```

（ストリッパは最終的に `test/strip_c_comments.py` としてコミットし、ガード
スクリプトから共用する。stderr は捨てない。）

間接層（各 expec ファイルが呼ぶ `child_*`/`GC_child_*`/ヘルパ関数の定義ファイル）も同じ抽出を行い、結果を全て `docs/superpowers/specs/2026-07-11-expec-call-inventory.md` に「関数 × ローカルモード挙動（return-input / local-open / 防御ガード / 到達不能）」の表として記録する。

- [ ] **Step 2: partner_rank 監査（spec の証明義務）**

`nbody_correlation.c` / `anomalous_pair.c` の `*_partner_rank` 系関数を読み、複製 FullDiag モード（`iFlgScaLAPACK=1`、MPI 分離サイトなし = `Nsite == NsiteMPI`）で partner/origin が常に `myrank` になる根拠（分離サイトのビットが存在しない）をインベントリ文書に記録する。結論が「証明できない」場合はその分岐の具体的な到達条件を書き、Task 3 の防御ガードを「必須の安全網」として明記する（どちらでもガードは入れる）。

- [ ] **Step 3: ガードスクリプトを書く**

`test/check_expec_local_calls.sh`: 検査対象ファイル集合は **Step 1 と同じ既存ファイルのみ**
（`phys_distributed.c` はまだ存在しない — **Task 6 が作成と同時にこのスクリプトの
リストへ追加する**。スクリプト冒頭に `FILES="..."` 変数と
`# Task 6 adds src/phys_distributed.c here` コメントを置く）。
処理: `test/strip_c_comments.py`（Task 1 でコミット、失敗・空出力で即エラー、
stderr 温存）でコメント除去 → 許可リスト（インベントリ文書の表と一致:
`SumMPI_dc SumMPI_d SumMPI_li SumMPI_i fopenMPI childfopenMPI stdoutMPI`）に
無い `MPI_[A-Z]` / `exitMPI` / wrapperMPI 呼び出しが現れたら exit 1。
**防御ガード済み例外のマーカー規約（Task 3 と共通、両タスクの正）**:
ガード済み領域は元ソース上で `/* EXPEC_LOCAL_GUARDED_BEGIN */` と
`/* EXPEC_LOCAL_GUARDED_END */` の**行マーカーで囲む**（領域方式）。
スクリプトは**コメント除去前の元ソース**で BEGIN/END 行の行番号範囲を先に
収集し、コメント除去後のマッチのうちその範囲内のものを除外する（2 パス）。
行番号ハードコード禁止。

- [ ] **Step 4: 現状で実行し、未ガードの生 MPI 箇所（個数はインベントリが正 —
  概数を受け入れ基準にしない）が検出されることを確認 → Task 3 完了までは
  スクリプト冒頭の名前付き変数 `TEMP_UNGUARDED_FILES="src/nbody_correlation.c src/anomalous_pair.c"`
  （コメント `# Task 3 must empty this variable`）で当該ファイルの生 MPI 検査
  のみ一時スキップして PASS させる。Task 3 がこの変数を空にする。**

- [ ] **Step 5: ctest 登録・ローカル確認・コミット**

```bash
cd build_noMPI && cmake .. > /dev/null && ctest -R check_expec_local_calls --output-on-failure
git add docs/superpowers/specs/2026-07-11-expec-call-inventory.md test/check_expec_local_calls.sh test/CMakeLists.txt
git commit -m "Freeze the ExpecLocal call inventory with a guard test"
```

---

### Task 2: `ExpecMode` キーワード（パース・検証・降格規則）

**Files:**
- Modify: `src/include/DefCommon.h`（`#define NUM_EXPECMODE 3` と `EXPECMODE_SERIAL 0 / EXPECMODE_STATEPARALLEL 1 / EXPECMODE_TRACE 2`）
- Modify: `src/include/struct.h`（`int iExpecMode;` を `iSolver` 群の直後に）
- Modify: `src/readdef.c`（初期化 `X->iExpecMode=EXPECMODE_SERIAL`、パース節 `CheckWords(ctmp,"ExpecMode")`、検証を Solver 検証群の**後**に追加）
- Modify: `src/ErrorMessage.c`/`.h`（`cErrExpecMode`: 値域と有効条件の説明）
- Test: `test/fulldiag_solver_keyword.sh`（ケース追加）

**Interfaces:**
- Produces: `X->Def.iExpecMode`（ReadcalcmodFile 完了後は検証・降格済み）。降格規則: `nproc==1 → 0`（INFO）、3a では `2 → 1`（INFO "ExpecMode 2 kernels are not available in this build; running as ExpecMode 1."。この降格は readdef でなく **phys 側**で行う — 3b 導入時に readdef を触らないため。readdef は 0-2 を受理して検証のみ）。

- [ ] **Step 1: テスト追加（失敗確認→実装→通過）**

`fulldiag_solver_keyword.sh` に、**既存の最後のケースの後ろへ、その時点の
採番規則に続けて** 2 ケースを追加する（以下の (6)/(7) は説明用の仮番号）:
(6) `ExpecMode 1` + `Solver 0` → エラー終了すること（非分散ソルバー）; (7) ELPA ビルド（`HPHI_HAS_ELPA=1`）でのみ: `Solver 3` + `ExpecMode 1` がシリアル実行（`MPIRUNFC` 空 = nproc 1）で INFO を出して正常終了し、結果が `ExpecMode 0` と一致すること。検証コード（readdef.c、Solver 検証群の直後）:

```c
  if (ValidateValue(X->iExpecMode, 0, NUM_EXPECMODE - 1)) {
    fprintf(stdoutMPI, cErrExpecMode, defname);
    return (-1);
  }
  if (X->iExpecMode != EXPECMODE_SERIAL) {
    if (X->iCalcType != FullDiag ||
        (X->iSolver != SOLVER_SCALAPACK && X->iSolver != SOLVER_ELPA)) {
      fprintf(stdoutMPI, cErrExpecMode, defname);
      return (-1);
    }
    if (nproc == 1) {
      fprintf(stdoutMPI,
        "  INFO: ExpecMode reverts to 0 for a single process (results are identical).\n");
      X->iExpecMode = EXPECMODE_SERIAL;
    }
  }
```

（`cErrExpecMode = "Error in %s\n ExpecMode: must be 0 (serial), 1 (state-parallel), or 2 (trace kernels),\n and requires CalcType = FullDiag with Solver 1 (ScaLAPACK) or 3 (ELPA).\n"`）

- [ ] **Step 2: 回帰（build_noMPI 17/17 + 新ケース）→ コミット** `git commit -m "Add ExpecMode keyword with validation and single-process downgrade"`

---

### Task 3: ExpecLocal フック（wrapperMPI / FileIO / 防御ガード）

**Files:**
- Modify: `src/wrapperMPI.c` / `src/include/wrapperMPI.h`
- Modify: `src/FileIO.c`（`fopenMPI` のローカルモード）
- Modify: `src/nbody_correlation.c` / `src/anomalous_pair.c`（partner≠myrank 分岐の防御ガード）
- Modify: `test/check_expec_local_calls.sh`（Task 1 の一時除外を空にする）

**Interfaces:**
- Produces:
  ```c
  void ExpecLocalEnter(void);   /* 入れ子不可（アサート）。エラーフラグもクリア */
  void ExpecLocalLeave(void);
  int  ExpecLocalActive(void);  /* iExpecLocal の読み取り */
  void ExpecLocalSetError(void);/* ローカルループ中の遅延エラー通知（常に定義） */
  int  ExpecLocalError(void);   /* 蓄積エラーの読み取り（常に定義） */
  ```
  フラグ実体は wrapperMPI.c の static とし、外部はアクセサのみ使用。
  フック対象（Task 1 のインベントリが確定させた表に従う。既定案）:
  `SumMPI_dc/_d/_li/_i` → ローカル時は入力を返す（MPI を呼ばない）;
  `MaxMPI_li/_d`, `BcastMPI_li`, `NormMPI_dc`, `VecProdMPI`, `BarrierMPI` →
  インベントリで expec 到達が確認されたもののみ同様に無通信化、到達しない
  ものは「ローカル中の呼び出しをデバッグアサートで禁止」;
  `fopenMPI` → ローカル時は呼び出しランクで開く。

- [ ] **Step 1: フック実装**（各関数の先頭に `if (iExpecLocal) return <無通信結果>;`。`NormMPI_dc`/`VecProdMPI` のローカル値はローカル和をそのまま返す — 関数本体のローカル計算部を通ってリダクションだけ飛ばす形にする）
- [ ] **Step 2: 防御ガード**: nbody_correlation.c / anomalous_pair.c の `#ifdef MPI` 生 Sendrecv ブロックの直前（partner/origin≠myrank が確定した位置）に:

```c
    if (ExpecLocalActive()) {
      /* EXPEC_LOCAL_GUARDED: cross-rank exchange is unreachable in the
         replicated FullDiag basis (see the call-inventory audit); if we
         ever get here in local mode, fail this state instead of touching
         raw MPI, which would deadlock the state-parallel loop. */
      fprintf(stdout, "  Error: cross-rank term reached in ExpecMode local loop.\n");
      return /* 各関数のエラー慣例に従う（dam_pr 経路なら 0 を返して上位 rc を立てるのではなく、関数シグネチャに応じてエラー値。実装時に各関数の戻り値契約を確認して統一） */;
    }
```

ガードの統一契約（全箇所共通）: `ExpecLocalSetError()` を呼んでから
`return 0.0;`（double complex 返しの関数）または関数のエラー慣例値を返す。
NaN センチネルは使わない（0.0 + SetError で十分かつ一貫）。ガード領域は
`/* EXPEC_LOCAL_GUARDED_BEGIN */` / `_END` 行マーカーで囲む（Task 1 の
スクリプト規約と同一 — 規約の正は両タスクに同文で記載済み）。
対象箇所は Task 1 インベントリの生 MPI 一覧**全件**（概数でなくリストが正）。
この契約をインベントリ文書に追記。

- [ ] **Step 3: アサート**: `ExpecLocalEnter` で `assert(!iExpecLocal)`、`Leave` で `assert(iExpecLocal)`。
- [ ] **Step 4: ガードスクリプトの一時除外を撤去し PASS 確認 → 回帰 17/17 → コミット** `git commit -m "Add ExpecLocal hook with no-comm reductions and raw-MPI guards"`

---

### Task 4: green_output パーシャル出力＋マニフェスト

**Files:**
- Modify: `src/green_output.c` / `src/include/green_output.h`
- Modify: 集約ファイルを開いている expec 側の呼び出し箇所（`grep -rn "GreenOutputFileName" src/*.c` で全列挙し、`childfopenMPI(GreenOutputFileName...)` 型のオープン/クローズを新ヘルパに置換）

**Interfaces:**
- Produces:
  ```c
  void GreenOutputSetPartialSuffix(int rank);   /* ローカルモード開始時 */
  void GreenOutputClearPartialSuffix(void);
  int  GreenOutputOpenAggregate(struct BindStruct *X, GreenOutputKind kind, FILE **fp);
      /* パーシャル時: final 名から part_path を導出(".part%d")、unlink→open、
         マニフェスト {attempted=1, opened, open_error} を更新。
         非パーシャル時: 従来の childfopenMPI と同動作 */
  int  GreenOutputCloseAggregate(GreenOutputKind kind, FILE *fp);
      /* bytes=ftell, closed_ok を記録して fclose */
  int  GreenOutputMergePartials(struct BindStruct *X);
      /* rank 0: マニフェストを Gather（固定長レコード×kind数）、
         全ランク成功時のみ final 名を初期化(トランケート)→ランク順に
         part を連結→part 削除。attempted&&open_error はエラー。
         正当な空(attempted&&!open_error&&bytes==0)は正常 */
  ```
- 注意: `GreenOutputInitializeAggregateFiles()` の呼び出しは Mode 1 経路では**行わない**（Merge が公開時に初期化する）。Mode 0 経路は従来どおり。

- [ ] **Step 1: 実装**（マニフェストは `static` 配列 `[GreenOutputAnomalous+1]` の
  レコード構造体 {attempted, opened, open_error, bytes, closed_ok, part_path,
  final_path}。**bytes は spec の rows の実装形**（ftell による正当な空との区別
  という意図は同一 — spec 側の字句も bytes に合わせて 1 行修正すること）。
  part_path/final_path は別フィールド。Gather は `MPI_Gather` 固定長。
  **green_output.c は build_noMPI でもコンパイルされるため、MPI を使う
  Gather/Merge 部は `#ifdef MPI` で包み、非 MPI ビルドでは Merge は
  非パーシャル動作の no-op にする**）
- [ ] **Step 2: expec 側の集約オープン箇所を新ヘルパに置換**（機械的置換。置換一覧をレポートに）
- [ ] **Step 3: 回帰 17/17（非パーシャル経路の等価性はこれが担保）→ コミット** `git commit -m "Add manifest-based partial aggregate output to green_output"`

---

### Task 5: 状態パネル再分散 `RedistBlockCyclicToStatePanel`

**Files:**
- Modify: `src/matrixscalapack.c` / `src/include/matrixscalapack.h`
- Create: `test/unit/elpa_statepanel_check.c`
- Modify: `test/CMakeLists.txt`（`if(USE_ELPA)` 内に elpa_redist_check と同型で登録、`min:2`）

**Interfaces:**
- Produces:
  ```c
  /* Z(2D block-cyclic, descZ) -> state panel (this rank owns full
     eigenvectors for 1-based states [jbegin, jbegin+ncols-1]).
     panel: N x ncols column-major, ld = N. Collective. Returns 0/-1. */
  int RedistBlockCyclicToStatePanel(long int xNsize,
                                    double complex *Z, int *descZ,
                                    long int jbegin, long int ncols,
                                    long int panel_ld, double complex *panel);
  ```
  実装は `RedistPanelToBlockCyclic` の逆向き（1D 側が**宛先**）: 同じ 1×P 'R' グリッド、`desc1d(M=N,N=N,MB=N,NB=NC,RSRC=CSRC=0,LLD=panel_ld)`、`pzgemr2d_(N, N, Z, 1,1, descZ, panel, 1,1, desc1d, &descZ[1])`
  — **この呼び出し表記は模式図（SCHEMATIC）**。実引数は全てポインタ渡し・
  整数幅は既存宣言厳守で、`src/matrixscalapack.c` の `RedistPanelToBlockCyclic`
  の実装をソース/宛先の役割だけ入れ替えて**そのまま写す**こと（フェーズ1/2 で
  幅バグが 2 度出た箇所）。所有権整合ガード（mycol==myrank + 呼び出し引数と
  導出値の照合、Allreduce(MIN) 同期、フェーズ2 と同一パターン）。

- [ ] **Step 1: 実装**（`RedistPanelToBlockCyclic` を鏡映しにする。整数幅は既存宣言厳守）
- [ ] **Step 2: 単体テスト** `elpa_statepanel_check.c`: N=97 の決定行列を `DivMat` で 2D 分散 → `RedistBlockCyclicToStatePanel` → 各所有状態 n についてパネル列を `MatElem(i,n)` と直接比較（元行列の列 = この行列を「固有ベクトル行列」と見なした検証。対角化不要で純粋にデータ移動を検査）。任意 np、`MPI_Allreduce(MIN)` で ok 共有。elpa_redist_check と同じ CMake 材料（`_ec_sc_libs` 等）で登録。
- [ ] **Step 3: 既定ビルド無変更確認（`ctest -N | grep -i elpa` 空）→ コミット** `git commit -m "Add block-cyclic to state-panel redistribution with unit check"`

---

### Task 6: Mode 1 ドライバ + phys.c 分岐 + Mode 0 S²/Sz 統一

**Files:**
- Create: `src/phys_distributed.c` / `src/include/phys_distributed.h`
- Modify: `src/phys.c`（冒頭分岐＋Mode 0 例外＋stdout 形式統一）
- Modify: `src/CMakeLists.txt`（ソースリストに `phys_distributed.c`）

**Interfaces:**
- Consumes: Task 2 の `iExpecMode`、Task 3 の `ExpecLocalEnter/Leave/Active` + `iExpecLocalError`、Task 4 の GreenOutput 系、Task 5 の `RedistBlockCyclicToStatePanel`、グローバル `Z_vec`/`descZ_vec`/`v0`/`use_scalapack`
- Produces: `int phys_stateparallel(struct BindStruct *X, unsigned long int neig);`

- [ ] **Step 1: ドライバ実装**（spec §3 の擬似コードを忠実に）:

```c
int phys_stateparallel(struct BindStruct *X, unsigned long int neig) {
  long int NN = (long int)neig, P = nproc;
  long int NC = (NN + P - 1) / P;
  long int jb = (long int)myrank * NC + 1;
  long int je = ((long int)myrank + 1) * NC; if (je > NN) je = NN;
  long int ncols = (je >= jb) ? (je - jb + 1) : 0;
  int rc_local = 0;
  double complex *panel;
  long int n, j;

  panel = malloc((((NN * ncols) > 0) ? NN * ncols : 1) * sizeof(double complex));
  /* malloc NULL チェックを Allreduce(MIN) で同期（フェーズ2 lapack_diag_elpa と同型）→ 失敗時は全ランク -1 */
  if (RedistBlockCyclicToStatePanel(NN, Z_vec, descZ_vec, jb, ncols, NN, panel) != 0) { free(panel); return -1; }
  free(Z_vec); Z_vec = NULL;

  if (GreenOutputInitializeAggregateFiles 相当の初期化はここでは呼ばない /* Merge が公開時に行う */);
  ExpecLocalEnter();
  GreenOutputSetPartialSuffix(myrank);
  for (n = jb; n <= je && rc_local == 0; n++) {
    X->Phys.eigen_num = n - 1;                      /* 既存 phys.c と同じ 0-based */
    for (j = 0; j < NN; j++) v0[j + 1] = panel[(n - jb) * NN + j];
    if (expec_energy_flct(X) != 0) { rc_local = -1; break; }
    if (expec_cisajs(X, v1) != 0)  { rc_local = -1; break; }
    if (expec_cisajscktaltdc(X, v1) != 0) { rc_local = -1; break; }
    if (expec_nbodyg(X, v1) != 0)  { rc_local = -1; break; }
    if (expec_anomalousg(X, v1) != 0) { rc_local = -1; break; }
    if (X->Def.iCalcType == FullDiag && expec_totalspin(X, v1) != 0) { rc_local = -1; break; }
    if (iExpecLocalError) { rc_local = -1; break; }
    /* all_* への記録（既存 phys.c 末尾と同じ代入群、インデックスは n-1） */
  }
  GreenOutputClearPartialSuffix();
  ExpecLocalLeave();
  free(panel);
  /* 単一ランデブー */
  { int g; MPI_Allreduce(&rc_local, &g, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (g != 0) return -1; }
  /* all_* の Gatherv: 一時受信バッファに gather し rank 0 で X->Phys.all_* へ
     コピー（MPI_IN_PLACE のエイリアス問題を避ける）。
     rank 0 が状態順に i=... 行を再レンダリング出力（シリアル形式・S2 列あり） */
  if (GreenOutputMergePartials(X) != 0) {
    return -1;   /* マニフェストがエラーを報告: 公開しない（集団的に失敗） */
  }
  return 0;
}
```
（`GreenOutputMergePartials` は内部で全ランクのマニフェスト Gather を行う
集団関数なので、全ランクが同順で呼ぶ。戻り値も全ランク一致で返す設計に
する — rank 0 の判定を Bcast/Allreduce で共有。）

（`expec_*` の実シグネチャ・`all_*` 代入群・進捗行のフォーマットは `src/phys.c` の現物から正確に写す。Gatherv の recvcounts/displs は全ランクの ncols から導出。）

- [ ] **Step 2: phys.c の分岐と Mode 0 例外**: `phys()` 冒頭に

```c
#ifdef _SCALAPACK
  if (use_scalapack && X->Def.iExpecMode != EXPECMODE_SERIAL) {
    if (X->Def.iExpecMode == EXPECMODE_TRACE)
      fprintf(stdoutMPI, "  INFO: ExpecMode 2 kernels are not available in this build; running as ExpecMode 1.\n");
    if (phys_stateparallel(X, neig) != 0) exitMPI(-1);
    return;
  }
#endif
```

Mode 0 分散経路の変更（唯一の例外）: 既存の `use_scalapack` 分岐内、S²/Sz をゼロ埋めしている箇所（phys.c:168-169 付近）と S2 列なし printf（phys.c:197 付近）を、

```c
      ExpecLocalEnter();
      if (myrank == 0) {
        if (expec_totalspin(X, v1) != 0) { ExpecLocalLeave(); exitMPI(-1); }
      } else { X->Phys.s2 = 0.0; X->Phys.Sz = 0.0; }
      ExpecLocalLeave();
      /* rank 0 の値が正; 表示・all_* は rank 0 の値を使う（従来から表示は
         stdoutMPI=rank 0。all_* の Gather は Mode 0 では不要 — rank 0 が正しい値を持つ） */
```

に置換し、printf をシリアル形式（S2 列あり）に統一する。**注意**: 非分散（`!use_scalapack`）経路は一切触らない。

- [ ] **Step 3: `phys_distributed.c` の本体全体を `#ifdef _SCALAPACK` で包む**
  （既定 build_noMPI には mpi.h も `_SCALAPACK` も無い。matrixscalapack.c と
  同じ「空の翻訳単位」パターン。ヘッダの宣言も同様にガード）。phys.c 側の
  分岐は既に `#ifdef _SCALAPACK` 内。**Task 1 のガードスクリプトの `FILES` に
  `src/phys_distributed.c` を追加する**（作成と同一コミットで）。
- [ ] **Step 4: `phys()` の全出口（分散分岐の早期 return 含む）で
  `assert(!ExpecLocalActive())` をデバッグアサート**（spec の出口保証）。
- [ ] **Step 5: ビルド＋回帰 17/17（既定ビルドでは新ファイルは空 TU）＋
  check_expec_local_calls PASS → コミット** `git commit -m "Add state-parallel observables driver and unify distributed Mode 0 S2/Sz"`
- [ ] **Step 6（早期チェックポイント）: clavius に rsync（リポジトリルートから）し
  `build_elpa` で **ビルドのみ**通す**（`make HPhi -j32`。Task 5/6 の新コードは
  ローカルでは一切コンパイルされないため、構文・幅エラーをここで前倒し検出。
  実行は Task 9）。

---

### Task 7: モード等価性テスト

**Files:**
- Create: `test/fulldiag_expecmode_equiv.sh`（+x）
- Modify: `test/CMakeLists.txt`（`if(USE_ELPA)` 内、`add_hphi_mpi_test(fulldiag_expecmode_equiv min:2)`）
- Modify: `test/fulldiag_elpa_hubbard_chain.sh`（S2/Sz 除外の回避コメントと列選別を撤去し、S² 列込み比較へ強化 — spec §6 の棚卸し）

**Interfaces:**
- Consumes: Task 2-6 の全成果物

- [ ] **Step 1: 等価性テストスクリプト**: 3 ケース（Hubbard 鎖 L=4 一体+二体GF・集約形式 ON・**さらに NBodyG 定義を追加してフォールバック経路と NBody 集約 kind も演習**（既存 `fulldiag_hubbard_nbody_interall` テストの def を流用）、SpinGC Gamma=0.5 L=6、Spin 鎖 L=8）× {ExpecMode 0, 1}（3a では 2 は 1 と同動作なので 2 も 1 ケースだけ回して INFO と一致を確認）で実行し、`zvo_phys_*` 全列（S²/Sz 込み）と全 Green ファイル（状態別・集約とも）を `paste`+awk 1e-8 比較。集約形式は `OutputGreenFormat` の集約値を calcmod に指定（既存 green_output_format テストの指定方法を流用）。np は `${MPIRUN}`（min:2）で、追加で np=3（**非整除の検証。ゼロ所有状態
ランクはこれらの N では発生しない — その経路は Task 5 の
`elpa_statepanel_check` を小さな N 引数で回して担保**する: 単体テストに
`argv[1]` で N を渡せるようにし、ctest 登録に N=4 np>4 相当のケースを追加）
を script 内の 2 回目の mpirun で実行。
**集約 kind の網羅**: OneBody/TwoBody/NBody は上記ケースで演習される。
ThreeBody/FourBody/SixBody は `expec_cisajscktaltdc` の多体出力を持つ入力
（既存テストの def を流用できるものがあれば追加）で 1 ケース演習し、
入力で到達させられない kind と AnomalousG は「分散 FullDiag での到達可否」を
Task 1 インベントリ監査の結論として文書に記録する（spec §2 の監査項目）。
**失敗注入**: 等価性スクリプトの最後に「part ファイルを 1 つ書き込み後に
削除 → Merge が非ゼロで失敗すること」を確認するシナリオを追加する。
これがスクリプト単体で困難なら Task 9（clavius）での手動確認項目として
明記し、スキップ理由をスクリプトのコメントに残す。
- [ ] **Step 2: 既定ビルドで未登録確認 → sh -n → コミット** `git commit -m "Add ExpecMode equivalence tests and strengthen the ELPA chain test"`

---

### Task 8: ドキュメント

**Files:**
- Modify: `doc/ja/source/filespecification/expertmode_ja/CalcMod_file_ja.rst` / `doc/en/.../CalcMod_file_en.rst`（`ExpecMode` エントリ新設: 値表・有効条件・「速度のみ」保証と丸め注記・S²/Sz 挙動修正・一時メモリ 2×O(N²/P) 注記・利用指針。`Solver` エントリの近く、既存様式厳守）
- Modify: `test/manual/elpa_gpu_check.md`（チェックリストにフェーズ3 項目: equiv テスト、Mode 1 ベンチ、S²/Sz 確認）

- [ ] **Step 1: ja/en 追記（同内容・各言語様式）**
- [ ] **Step 2: PR 移行ノートの草稿**を `docs/superpowers/specs/2026-07-11-phase3a-migration-note.md` に作成
  （内容: 分散 Mode 0 で S²/Sz がゼロ埋めでなくなる・stdout がシリアル形式に
  統一される・ExpecMode は結果を変えない[丸め除く]。push 時に PR 説明文へ
  転記するための原稿 — spec §8 の成果物）。
- [ ] **Step 3: コミット** `git commit -m "Document ExpecMode and the distributed S2/Sz unification"`

---

### Task 9: clavius 実機検証（コントローラ直接実行）

- [ ] rsync（リポジトリルートから）→ `build_elpa` 再構成（`-DELPA_ROOT=$HOME/opt/elpa-2025.06-cuda`）→ ビルド
- [ ] `elpa_statepanel_check` np=1,2,3,4,8 / `fulldiag_expecmode_equiv`（np=2,3,4）/ 既存 ELPA テスト全部
- [ ] **ベンチマークゲート**: L=8 Hubbard（N=4900）全状態の一体+二体GF 実時間を `ExpecMode 0` vs `1`（np=4）で比較 — 目標 ~P 倍（`CalcTimer.dat` の expec 区間 + wall）
- [ ] S²/Sz 統一の確認: 分散 Mode 0 の `zvo_phys` がシリアル実行と一致
- [ ] 結果を `test/manual/elpa_gpu_check.md` に「Phase 3a validation」節として追記・コミット

---

## 完了条件（フェーズ3a）

- 既定（非分散）全既存テスト無変更 PASS、`check_expec_local_calls` PASS（一時除外なし）
- ELPA ビルドで equiv テスト・statepanel 単体・既存 ELPA テスト全 PASS
- clavius で Mode 0/1 等価（S² 込み）とベンチ ~P 倍を記録
- spec v5.1 §2/§3/§5 の全項目に対応する実装・テスト・docs が存在
- PR 移行ノート草稿（S²/Sz・stdout 形式変更）が作成済み（push 時に PR 説明文へ転記）
- Mode 0 の変更で `zvo_phys` を書くのが rank 0 であることを実装時に output.c/phys.c で確認済み（all_* の Gather 不要判断の根拠）
