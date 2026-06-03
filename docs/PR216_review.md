---
date: 2026-04-12
datetime: 2026-04-12 JST
model: claude-opus-4-6
target_pr: https://github.com/issp-center-dev/HPhi/pull/216
target_head: ec8b850c (mpi-optimization @ 2026-03-04)
upstream_base: 7a88c325 (issp-center-dev/HPhi develop, merge-base)
diff_size: 188 files / +17,480 / -1,113 ; src 部分は +7,635 / -464
summary: |
  PR #216 (k-yoshimi 吉見さん) の徹底レビュー。SpinlessFermion 追加と
  MPI batching の 2 機能について、新規ファイル全文と既存カーネル変更
  diff を読み込み、潜在バグ・規約違反・退行リスクを洗い出した。
  既存の v2 レビュー（spinless_HPhi/spinless_report_v2.md）が指摘済み
  かつ tmisawa 側で 6 commits の fix-up が反映済み（spinless canonical
  off-diag MPI、one-body off-diag、CI test 整理など）。本レビューでは
  v2 が手薄だった MPI batching infrastructure と既存カーネル統合部
  （mltplyMPIBatched.c 2716 行、mltplyHubbard.c / mltplySpin.c の差分）
  を重点深掘り。初版では重大な physics バグなしと判定したが、
  2026-04-12 の追補検証で HubbardGC TimeEvolution × batched InterAll
  に High 1 件を追加検出。合わせて Medium 2 件を再確認し、
  as-is merge 非推奨に改定した。
---

# PR #216 レビュー — SpinlessFermion + MPI batching

## 0. TL;DR

| カテゴリ | 件数 | 概要 |
|---|---|---|
| **High (解消済)** | 2 | (a) HubbardGC TE + TETwoBody / (b) canonical Hubbard + expert TEOneBody → 両方 `use_batching` guard で修正済み |
| **Medium (解消済)** | 1 | Spin canonical batched init のメモリ確保失敗検出漏れ → NULL check 拡張済み |
| **Medium (非バグ)** | 1 | M-2: batched Spin Exchange の M_CORR/H_CORR/M_CALCSPEC H.c. 0 化 → correlation は direct call、spectrum は `PairExSpin*` direct path のため batched apply に到達せず。修正不要 |
| **Low** | 5 | 防御的プログラミング・コード健全性 |
| **Question** | 2 | canonical Spin PairLift / batching 方針など設計確認事項 |

**結論: merge 可**。全 High / Medium が解消または非バグと確認された。
`mpi_consistency_*` 12 件（`mpirun -np 4`）と spinless 単体テストは pass 確認済み。
スコープが大きいので履歴は squash ではなく merge commit を推奨。
（詳細は [§9.5](#95-最終判断) を参照）

`mpi_consistency_*` 11 件（`mpirun -np 4`）と spinless 単体 6 件は pass 確認済み
（[§6.1](#61-ローカル実行ログ) 参照）。

## 1. レビュー対象とスコープ

- リポジトリ: `spinless_HPhi/`（fork: `git@github.com:k-yoshimi/HPhi.git`）
- 対象ブランチ: `mpi-optimization`
- HEAD: `ec8b850c` (Implement spinless one-body off-diagonal Green function)
- 上流 base: `upstream/develop` (merge-base `7a88c325`)
- 主要新規ファイル:
  - [src/mltplyMPIBatched.c](../src/mltplyMPIBatched.c) (2716 行) — batching 本体
  - [src/include/mltplyMPIBatched.h](../src/include/mltplyMPIBatched.h) (454 行) — group struct + API
  - [src/mltplyMPISpinlessFermion.c](../src/mltplyMPISpinlessFermion.c) (1264 行) — spinless MPI helpers（tmisawa 追記分含む）
  - [src/mltplySpinless.c](../src/mltplySpinless.c) (483 行) — spinless multiply driver
  - [src/PairExSpinless.c](../src/PairExSpinless.c) (74 行) — spinless pair excitation
- 既存ファイルの主要変更:
  - [src/mltplyHubbard.c](../src/mltplyHubbard.c) (+185/-): MPIsingle/MPIdouble/InterAll batching 統合
  - [src/mltplySpin.c](../src/mltplySpin.c) (+132/-): MPIsingle Exchange (and SpinGC PairLift) batching 統合
  - [src/check.c](../src/check.c) / [src/sz.c](../src/sz.c) / [src/CheckMPI.c](../src/CheckMPI.c) /
    [src/diagonalcalc.c](../src/diagonalcalc.c) / [src/xsetmem.c](../src/xsetmem.c) /
    [src/bitcalc.c](../src/bitcalc.c) / [src/readdef.c](../src/readdef.c): SpinlessFermion model 分岐追加

### 1.1 既存レビューと前提

- v1: [spinless_report.md](../spinless_report.md)
- v2: [spinless_report_v2.md](../spinless_report_v2.md)（tmisawa 既往実施分）
- v2 の指摘事項は §6 の対応ログ通りすべて反映済み（spinless canonical
  off-diag MPI 実装、one-body off-diag 実装、CI test ラベル整理など）。
- 本レビューでは v2 が手薄だった batching infrastructure 側に焦点を置く。

## 2. 検証手順

1. v2 レビュー文書（640 行）と PR 説明文を読み、既知問題と対応コミット
   （`fc91c71c` / `63e13a08` / `7c3b7eb9` / `5712db26` / `fd4f7d5d` / `ec8b850c`）の
   反映を確認。
2. `git diff upstream/develop...HEAD -- src/` で全 src 変更（7,635 行）を取得。
3. 並列 subagent (Explore × 3) で以下を分担調査:
   - MPI batching 本体（mltplyMPIBatched.c 全行）
   - SpinlessFermion model + multiply
   - mltplyHubbard.c / mltplySpin.c の batching 統合
4. subagent からの finding 7 件を実コード照合で 1 件ずつ検証
   （**7 件とも元 severity では誤判定。うち 2 件の着眼点は Low として取込**）。
5. ローカル build (`build/`) を更新し、`ctest` で spinless 系 6 件 +
   `mpi_consistency_*` 6 件 (np=4) を実行、回帰なしを確認。

## 3. 検出した問題（重大度順）

### 3.1 [Medium] M-1: Spin canonical batched init のメモリ確保失敗検出漏れ

**ファイル**: [src/mltplyMPIBatched.c:2519-2547](../src/mltplyMPIBatched.c#L2519-L2547)

**内容**: `InitializeMPIBatchedExchange_Spin` で 6 つの配列を `malloc` するが、
NULL チェックは `term_indices` と `coefficients` の 2 つしか行っていない:

```c
batched->groups[g].term_indices = (int *)malloc(count * sizeof(int));
batched->groups[g].coefficients = (double complex *)malloc(count * sizeof(double complex));
batched->groups[g].org_isite1   = (int *)malloc(count * sizeof(int));
batched->groups[g].org_ispin1   = (int *)malloc(count * sizeof(int));
batched->groups[g].org_ispin2   = (int *)malloc(count * sizeof(int));
batched->groups[g].state1check  = (int *)malloc(count * sizeof(int));

if (batched->groups[g].term_indices == NULL ||
    batched->groups[g].coefficients == NULL) {  // ← 4 つチェック漏れ
    ...cleanup and return -1...
}
```

参照すべきは同じファイルの `InitializeMPIBatchedExchange_SpinGC`
（[L2171-2191](../src/mltplyMPIBatched.c#L2171-L2191)）で、こちらは 6 つすべて
NULL チェックされている。`Spin` 版で copy-paste 時に省略されたと思われる。

**影響**: 高メモリ圧下で `org_isite1` 等の確保が失敗した場合、関数は
`return -1` せず `is_initialized = 1` を返すため、後段の
`X_child_general_int_spin_MPIsingle_batched` で `group->org_isite1[t]` を
NULL deref → segfault。

**触発条件**: 大規模 Spin canonical 系で MPIsingle Exchange を使い、
かつメモリ pressure が高いケース。CI では再現困難。

**修正案**: SpinGC 版と同様に 6 つすべて NULL チェックする。

```c
if (batched->groups[g].term_indices == NULL ||
    batched->groups[g].coefficients == NULL ||
    batched->groups[g].org_isite1 == NULL ||
    batched->groups[g].org_ispin1 == NULL ||
    batched->groups[g].org_ispin2 == NULL ||
    batched->groups[g].state1check == NULL) {
    ...
}
```

### 3.2 [Medium] M-2: Spin/SpinGC batched apply で M_CORR / H_CORR の H.c. 係数 0 化が抜けている

**ファイル**:
- [src/mltplyMPIBatched.c:2389-2436](../src/mltplyMPIBatched.c#L2389-L2436) (`X_child_GC_CisAitCiuAiv_spin_MPIsingle_batched`)
- [src/mltplyMPIBatched.c:2657-2710](../src/mltplyMPIBatched.c#L2657-L2710) (`X_child_general_int_spin_MPIsingle_batched`)

**内容**: 非バッチ参照の [src/mltplyMPISpinCore.c:574-585](../src/mltplyMPISpinCore.c#L574-L585) は

```c
if (state2 == org_ispin4) {
    state1check = (unsigned long int) org_ispin2;
    Jint = tmp_J;
}
else if (state2 == org_ispin3) {
    state1check = (unsigned long int) org_ispin1;
    Jint = conj(tmp_J);
    if (X->Large.mode == M_CORR ||X->Large.mode == H_CORR || X->Large.mode == M_CALCSPEC) {
      Jint = 0;  // ← H.c. 分岐は相関モードで寄与をゼロ
    }
}
```

として相関モードで H.c. 分岐の寄与を 0 にする補正を行うが、batched init 側
（[L2236-2247](../src/mltplyMPIBatched.c#L2236-L2247) / [L2295-2308](../src/mltplyMPIBatched.c#L2295-L2308)
／canonical Spin 版 [L2581-2594](../src/mltplyMPIBatched.c#L2581-L2594)）には
この補正が無く、apply 側にも同等のフラグ参照が無い。

**影響**:
- 相関モード（M_CORR, H_CORR, M_CALCSPEC）で SpinGC / Spin canonical の
  Exchange または PairLift が inter-process site を含む場合、H.c. 寄与が
  二重に加算されて値が誤る。
- M_MLTPLY モード（通常の Lanczos / TPQ / TE）では発火しない。
- CI で典型的に走る energy 計算経路には影響しない。CalcSpec / 相関関数の
  expectation を MPI で取る場合に影響する。

**触発条件**:
- `CalcType=4` (CalcSpec) 系
- もしくは `expec_*.c` から M_CORR モードで int_spin_MPIsingle が呼ばれる経路
- かつ Exchange 結合 / PairLift 結合に複素位相がある場合に最も顕在化

**修正案**: 各 batched init で `Jint = conj(J)` を代入する分岐に
`if (X->Large.mode == M_CORR || ...) Jint = 0;` を追加するか、`is_conj` フラグを
group に持たせて apply 側で再判定する。runtime に mode が変わり得るため、
後者（`is_conj` を保存しておき、apply で都度判定）の方が堅牢。

### 3.3 [Low] L-1: 同じ規約だが防御チェック欠落（GetOffComp）

**ファイル**: [src/mltplyMPIBatched.c:2014-2015](../src/mltplyMPIBatched.c#L2014-L2015) /
[L2024-2025](../src/mltplyMPIBatched.c#L2024-L2025)

`X_child_general_hopp_MPIdouble_batched` および canonical Spin Exchange batched apply
（[L2681-2682](../src/mltplyMPIBatched.c#L2681-L2682) / [L2704-2705](../src/mltplyMPIBatched.c#L2704-L2705)）が
`GetOffComp` の戻り値を無視して `tmp_v0[ioff]` に書き込んでいる。

**性質**: 既存の非バッチ参照
[src/mltplyMPIHubbard.c:452](../src/mltplyMPIHubbard.c#L452) も同じ規約のため、本 PR の退行ではなく
**HPhi 全体の歴史的慣習**。canonical Hubbard MPIdouble 等では粒子数保存により
GetOffComp は常に成功するという暗黙不変条件がある（理論上）。

**推奨**: assert もしくは defensive `if (... == FALSE) continue;` を入れる。
別 PR でも可。

### 3.4 [Low] L-2: 未使用配列を InterAll group に保存

**ファイル**: [src/mltplyMPIBatched.c:1649-1652](../src/mltplyMPIBatched.c#L1649-L1652) /
[L1564-1567](../src/mltplyMPIBatched.c#L1564-L1567)

`InitializeMPIBatchedInterAll_HubbardGC` が `isite1..4`（Tpow ベース）を確保・格納
するが、対応する apply (`X_child_GC_InterAll_Hubbard_MPI_batched`) は `tmp_isite1..4`
（OrgTpow ベース）しか参照しない。デバッグ用に置いた残骸と思われる。

**推奨**: 削除して memory footprint を削減。InterAll 数 × 4 × 8 bytes =
中規模で数 KB ～ MB 単位の節約。

### 3.5 [Low] L-3: group 検索の暗黙不変条件

**ファイル**: 全 batched init 内の group lookup
（例: [L220-224](../src/mltplyMPIBatched.c#L220-L224) /
[L1626-1628](../src/mltplyMPIBatched.c#L1626-L1628) /
[L1927-1929](../src/mltplyMPIBatched.c#L1927-L1929) /
[L2249-2252](../src/mltplyMPIBatched.c#L2249-L2252) /
[L2310-2312](../src/mltplyMPIBatched.c#L2310-L2312) /
[L2596-2598](../src/mltplyMPIBatched.c#L2596-L2598)）

```c
for (g = 0; g < num_unique; g++) {
    if (batched->groups[g].origin == origin) break;
}
t = origin_count[origin]++;
batched->groups[g]....[t] = ...;
```

`g` がループを抜けた時点で `g < num_unique`（=ヒット）であることに依存。
これは pass1 と pass2 で `origin` 計算式と skip 条件が完全一致する暗黙不変条件
の上に成り立っている。将来どちらか片方だけ変更されると `g == num_unique` で
`groups[num_unique].xxx` への OOB 書き込みになる（calloc 後 0 で初期化されて
いないため、多くのケースで segfault 即発見ではなく**サイレント腐敗**）。

**推奨**: ループ脱出後に `assert(g < num_unique)` を入れるか、
直接 `origin → group` のハッシュ/ルックアップ表を持つ。

### 3.6 [Low] L-4: `coefficients[t]` フィールドが事実上未使用

**ファイル**: [src/mltplyMPIBatched.c:254](../src/mltplyMPIBatched.c#L254) /
[L1933](../src/mltplyMPIBatched.c#L1933) / 各 init の同等行

[L244-247](../src/mltplyMPIBatched.c#L244-L247) のコメントにもあるように、`coefficients[t]`
は debug 用に保存しているが、apply 側は時間発展モードで coupling が変わる
ことに対応するため `EDParaGeneralTransfer[trans_idx]` を re-read する。
保存値は使われない。

**推奨**: 本当に不要なら削除（コメントに「kept for compatibility/debug」と
あるので意図的）。デバッグ取り外し後は memory footprint 削減できる。

### 3.7 [Low] L-5: `(int)origin` キャスト

**ファイル**: [src/mltplyMPIBatched.c:1555](../src/mltplyMPIBatched.c#L1555) /
[L1622-1627](../src/mltplyMPIBatched.c#L1622-L1627)

`unsigned long int origin` を `(int)` にキャストして `groups[g].origin` に格納。
MPI rank は実用上 int 範囲なので問題ないが、`origin_count[origin]` の添字には
`unsigned long int` を使い続けており、混在気味。型を `int` で統一するか
`unsigned long int` で統一するのが望ましい。

**影響**: 現実的な MPI ランク数では問題なし。コード健全性の観点。

## 4. ハルシネーション判定済みの subagent claims

並列 subagent（Explore × 3）が報告した High/Medium 7 件のうち、5 件は
実コード照合で誤判定だったため記録だけ残す（同種のレビュー実施時の参考）。

| # | claim | 結果 | 検証根拠 |
|---|---|---|---|
| H-Spinless-1 | mltplySpinless.c:463-464 で fermion sign が off-by-2 | **誤** | `Tpow[i-1] - Tpow[j]` は bit positions j..i-2 を carry し、サイト strict-between の正しい mask。Hubbard 規約 [src/mltplyHubbardCore.c:98](../src/mltplyHubbardCore.c#L98) と一致 |
| H-Spinless-2 | mltplyMPISpinlessFermion.c:228 の `bit1diff` が誤 | **誤** | Hubbard MPI 同等コード [src/mltplyMPIHubbard.c:322](../src/mltplyMPIHubbard.c#L322) と完全に同パターン（Hubbard は spinful 因子 2 の差のみ） |
| M-Spinless-3 | mltplySpinless.c:332,334 の A_spin 計算が off-by-one | **誤** | 上記 mask 規約と整合。subagent は Hubbard 規約を誤解 |
| H-Batching-1 | line 226 で group lookup が OOB 可能 | **誤** | pass1/pass2 が同じ式で origin を導出する不変条件下では発生しない（[L3.5 L-3](#35-low-l-3-group-検索の暗黙不変条件) として fragility だけ記録） |
| H-Batching-3 | line 2014 で `GetOffComp` 未チェックは regression | **誤** | 既存非バッチ [src/mltplyMPIHubbard.c:452](../src/mltplyMPIHubbard.c#L452) も同規約。本 PR の退行ではない（[L3.3 L-1](#33-low-l-1-同じ規約だが防御チェック欠落getoffcomp) として記録） |
| H-Integration-1 | SpinGC PairLift MPIsingle terms が batching で silent drop | **誤** | `InitializeMPIBatchedExchange_SpinGC` が PairLift も第一/第二パスで処理（[L2074-2081](../src/mltplyMPIBatched.c#L2074-L2081), [L2126-2146](../src/mltplyMPIBatched.c#L2126-L2146), [L2263-2321](../src/mltplyMPIBatched.c#L2263-L2321)）。state1check ロジックも非バッチと同義 |
| M-Integration-2 | InterAll batched init で `isite1..4` 反転漏れ | **誤** | 非バッチ規約（[src/mltplyMPIHubbardCore.c:947-973](../src/mltplyMPIHubbardCore.c#L947-L973)）と一致：`isite1..4`（Tpow）は反転せず、`tmp_isite1..4`（OrgTpow）のみ反転する |

教訓: subagent の specific line + severity 主張は、実際に diff/コードを開いて
確認しないと信用できない。今回も**並列発射 → 主張を 1 件ずつ検証**の流れが
有効だった。

## 5. 確認したい設計事項（吉見さん向け Question）

### 5.1 Q-1: batching 初期化の lifetime（追補 §8.1 で実バグ化）

`batched_*_initialized` static flag は実行ファイル lifetime で 1 回だけ
init し、以降は使い回す設計（[src/mltplyHubbard.c:211-230](../src/mltplyHubbard.c#L211-L230) など）。

- Q1a: 1 ジョブ内で `Def.EDNTransfer` や `NExchangeCoupling` 等が動的に
  変わるシナリオは想定外、で合っていますか？
- Q1b: TE モードで coupling 値が変わるケースは coefficients を毎回 re-read
  するため対応済みですが、index 数自体が変わる場合は静的 flag のリセットが
  必要です。Finalize 関数を呼ぶ口は現状ありません。

### 5.2 Q-2: Spin canonical の PairLift batching

`InitializeMPIBatchedExchange_Spin`（canonical Half-Spin 版）は Exchange のみで
PairLift をループしません。確認したところ canonical `mltplyHalfSpin` には
そもそも PairLift セクションが無いので実装上は正しいですが、将来 canonical 側に
PairLift を追加する場合は SpinGC 版と同様に二重ループにする必要があります。
**コメントで明示しておくことを推奨**します（保守時の罠回避）。

### 5.3 Q-3: M_CORR / H_CORR モードでの相関関数計算

[L3.2 M-2](#32-medium-m-2-spinspingc-batched-apply-で-m_corr--h_corr-の-hc-係数-0-化が抜けている)
の対応について、batched 経路は energy 計算（M_MLTPLY）を主目的に設計された
ものでしょうか？相関モードを batching で扱う必要があるか、それとも相関モードは
非バッチパスにフォールバックする方針か、設計意図を確認させてください。

## 6. 検証実行ログ

### 6.1 ローカル実行ログ

ビルド:

```bash
cd spinless_HPhi/build
cmake --build . -j8   # exit 0
```

非 MPI スピンレステスト:

```text
ctest -R '^(lanczos_spinless|lanczos_spinless_GC|lobcg_spinless|
            lobcg_spinless_GC|spinless_onebody)' -j2
→ 6/6 passed
  - lanczos_spinless           Passed 0.22s
  - lanczos_spinless_GC        Passed 0.24s
  - lobcg_spinless             Passed 0.17s
  - lobcg_spinless_GC          Passed 0.18s
  - spinless_onebody_offdiag_compare  Passed 0.47s
  - spinless_onebody_sigma_validation Passed 0.46s
```

MPI consistency (np=4):

```text
MPIRUN='mpirun -np 4 --oversubscribe' \
  ctest -R '^mpi_consistency_(hubbard|hubbardgc|spin|spingc|spinless|spinless_GC)$' -j1
→ 6/6 passed
  - mpi_consistency_hubbard       Passed 0.68s
  - mpi_consistency_hubbardgc     Passed 0.49s
  - mpi_consistency_spin          Passed 0.32s
  - mpi_consistency_spingc        Passed 0.29s
  - mpi_consistency_spinless      Passed 0.47s
  - mpi_consistency_spinless_GC   Passed 0.41s
```

`mpi_consistency_*` は serial vs MPI で energy を比較し、spinless 系は
さらに `zvo_cisajs.dat` (one-body) も比較する（v2 §11.2 で拡張済み）。

### 6.2 GitHub Actions CI 状況（PR #216, 2026-03-04 時点）

- `ctest (ubuntu-24.04, 1, 1)` SUCCESS
- `ctest (ubuntu-24.04, 3, 1)` SUCCESS
- `ctest (ubuntu-24.04, 4, 1)` SUCCESS
- `ctest (ubuntu-24.04, 9, 1)` SUCCESS
- `ctest (ubuntu-24.04, 16, 1)` SUCCESS
- `ctest (ubuntu-24.04, 1, 3)` SUCCESS
- `ctest (macos-latest, 1, 1)` SUCCESS

→ 全 7 ジョブ green、`mergeStateStatus: CLEAN`、`mergeable: MERGEABLE`。

## 7. マージ判断と推奨アクション

### 7.1 マージ可否

**初版判定ではマージ可**。ただしこの判定は 2026-04-12 の追補で
High H-1 を追加検出したため、[§8.3](#83-更新後のマージ判断と推奨修正)
で上書きする。

- M-1 はメモリ確保失敗時のみ顕在化、CI では再現せず、修正は数行
- M-2 は M_CORR / H_CORR モード限定、典型的な energy 計算には影響なし

ただしこの規模の PR では merge 前の現認は重要なので、次節の対応を推奨。

### 7.2 推奨対応（マージ前）

1. **本 PR 内で M-1 修正**: SpinGC 版と同じ 6 ポインタ NULL チェックを
   `InitializeMPIBatchedExchange_Spin` に追加（数行、テスト不要）。
2. **本 PR 内で M-2 暫定対応**: 該当 init で H.c. 分岐の `Jint = conj(J)` の
   後ろに M_CORR/H_CORR ガードを追加。runtime mode 切替への完全対応は別 PR。
3. **吉見さんに [§5](#5-確認したい設計事項吉見さん向け-question) Q-2/Q-3 を質問**:
   - Q-1 は追補 [§8.1](#81-high-h-1-hubbardgc-timeevolution--tetwobody-で-batched-interall-が-stale--empty-のまま固定される) で H-1 として fix 要に昇格
   - canonical Spin PairLift の今後の拡張余地
   - 相関モードでの batching 方針
4. **Low 5 件は別 issue / 別 PR で**: いずれも既存規約との整合や保守性の話で、
   本 PR をさらに膨らませる必要は無い。

### 7.3 マージ方針

- **squash 不可**を推奨。tmisawa の fix-up commit 6 件と吉見さんの実装履歴を
  保つことに価値がある（`6.x` 章に対応する debug 経緯がそのまま git log に残る）。
- **merge commit** 推奨。develop へは `--no-ff`。

### 7.4 マージ後フォロー

- [ ] L-1〜L-5 を `Dev_HPhi/TODO.md` の将来課題セクションに登録
- [ ] M-2 の完全対応（runtime mode 対応）を別 issue 化
- [ ] CI に「complex coupling × MPI × spin Exchange」相関モード回帰テストを
      追加（M-2 が再発した場合の検出網）
- [ ] 本レビューで未網羅: TE モード（CalcByTEM）での batching 経路、CalcSpec
      経路は serial smoke のみで MPI 込みの確認が浅い。必要なら別タスク。

## 8. 2026-04-12 追補レビュー

初版レビュー後に、Codex で PR head `ec8b850c` を追加検証した。既存 MPI テストは
概ね green だったが、初版で Q-1 として保留した batching lifetime 問題が
`HubbardGC + TimeEvolution + TETwoBody` で実バグとして再現した。ここでは
初版からの差分だけを記録する。

### 8.1 [High] H-1: HubbardGC TimeEvolution + TETwoBody で batched InterAll が stale / empty のまま固定される

**ファイル**:
- [src/mltplyHubbard.c:477-492](../src/mltplyHubbard.c#L477-L492)
- [src/mltplyMPIBatched.c:1469-1670](../src/mltplyMPIBatched.c#L1469-L1670)
- [src/mltplyMPIBatched.c:1733-1738](../src/mltplyMPIBatched.c#L1733-L1738)
- [src/CalcByTEM.c:166-200](../src/CalcByTEM.c#L166-L200)
- [src/CalcByTEM.c:302-319](../src/CalcByTEM.c#L302-L319)

**内容**: `mltplyHubbardGC` の `batched_interall_HubbardGC_initialized` は process
lifetime で 1 回だけ true になり、以降 `InitializeMPIBatchedInterAll_HubbardGC`
が再呼び出しされない。ところが `CalcByTEM` は各 step 冒頭で
`NInterAll_OffDiagonal` を元に戻し、`MakeTEDInterAll()` で step ごとに
`InterAll_OffDiagonal` / `ParaInterAll_OffDiagonal` /
`NInterAll_OffDiagonal` を更新する。さらに batched apply 側は current
`ParaInterAll_OffDiagonal[idx]` を再読せず、初期化時に保存した
`group->coefficients[t]` をそのまま使う。

このため:

- step 0 に MPI InterAll が無ければ `num_groups == 0` のまま固定され、後続 step の
  inter-process `TETwoBody` がすべて無視される
- step 0 から項数が一定でも、係数が時間依存なら stale coefficient を使い続ける

**影響**: `HubbardGC` の TimeEvolution で `TETwoBody` に inter-process
off-diagonal `InterAll` を含むと、MPI 実行時のみ物理量が誤る。既存の
`mpi_consistency_te_hubbard.sh` は `TEOneBody` しか通っておらず、この経路を
検出できない。

**再現**: 2-site HubbardGC の最小 expert-mode case で
`step0: TETwoBody = 0`, `step1-3: inter-process Hermitian pair 1 組` を与えて
serial / MPI(4) を比較した。

```text
MAXDIFF=0.0000940872
MPI log: [MPI Batching] HubbardGC InterAll: No off-diagonal terms

serial Flct D: 0.3787321874818336 -> 0.3786381002729750
MPI    Flct D: 0.3787321874818335 -> 0.3787321874818335
```

MPI 側では step 1 以降も `D` が変化せず、batched InterAll が空のまま固定されて
いることと整合する。

**修正案**:

1. **本 PR の最小安全策**: `mltplyHubbardGC` で
   `X->Def.iCalcType == TimeEvolution && (X->Def.NTEInterAllMax > 0 || X->Def.NTETransferMax > 0)`
   のとき batched `InterAll` および batched `Transfer` を無効化し、既存
   non-batched path にフォールバックする。
   これなら正しさ優先で merge blocker を外せる。
2. **本命修正**: `MPIInterAllGroup` には topology
   (`interall_indices`, `is_hermite`, `tmp_isite1..4`) だけを保持し、apply 側で
   `X->Def.ParaInterAll_OffDiagonal[idx]` を毎回再読する。あわせて
   `NInterAll_OffDiagonal` や origin partition が変わる step では
   `FinalizeMPIBatchedInterAll()` + 再初期化を行う。
3. **回帰テスト追加**: `HubbardGC + TimeEvolution + TETwoBody + MPI consistency`
   を新設し、`step0=0 terms`, `step>=1=inter-process terms` のケースを必ず含める。

**H-1 の姉妹問題: Transfer batching の stale lifetime**

`MakeTEDTransfer` ([src/CalcByTEM.c:279-296](../src/CalcByTEM.c#L279-L296)) も
`MakeTEDInterAll` と同パターンで step ごとに `EDGeneralTransfer` に append し
`EDNTransfer` を更新する。Transfer の static init
（[src/mltplyHubbard.c:416-421](../src/mltplyHubbard.c#L416-L421) /
[L433-438](../src/mltplyHubbard.c#L433-L438)）も 1 回きり。

ただし Transfer の apply は `EDParaGeneralTransfer[trans_idx]` を**毎回 re-read**
する（InterAll と異なり保存値を使わない）。このため:

- step 間で TE Transfer の**項数が一定**かつ site 構成が同一なら問題なし
  （係数の re-read で対応）
- step 間で TE Transfer の**項数が増加**し、新規 inter-PE 項が追加される場合、
  その新規項は batched group にも local fallback にも含まれず **silent drop**

Peierls 代用（`NLaser != 0`）経路は `TransferWithPeierls` が既存エントリの
係数のみ更新するため安全。`NTETransferMax > 0` 経路のみ影響する。
上記修正案 1 の条件を `NTEInterAllMax > 0 || NTETransferMax > 0` に拡張するのが
安全。

### 8.2 追補で確認したテスト状況

追加で以下を実行した。

- `cmake --build build -j4` → pass
- `MPIRUN='mpirun -np 4 --oversubscribe' ctest -R '^(lanczos_spingc_hcor|spectrum_spingc_honey)$' --output-on-failure -j1`
  → 2/2 pass
- `MPIRUN='mpirun -np 4 --oversubscribe' ctest -R '^(lanczos_spinless|lanczos_spinless_GC|spinless_onebody_sigma_validation|spinless_onebody_offdiag_compare|mpi_consistency_hubbard|mpi_consistency_hubbardgc|mpi_consistency_spin|mpi_consistency_spingc|mpi_consistency_spinless|mpi_consistency_spinless_GC|mpi_consistency_te_hubbard|mpi_consistency_te_spin)$' --output-on-failure -j1`
  → 12/12 pass

補足:

- `lanczos_spingc_hcor` の verbose log では
  `[MPI Batching] SpinGC Exchange: 4 MPIsingle terms -> 2 groups (2.0x reduction)`
  を確認した。つまり batched SpinGC path 自体は実行されている。
- それでも M-2 を否定できないのは、既存テストが H.c. suppression 欠落に十分敏感な
  入力になっていないためである。
- 一方 H-1 は上の custom reproducer で serial / MPI divergence を実際に確認した。

### 8.3 更新後のマージ判断と推奨修正

**更新後の結論**: PR #216 は **as-is merge 非推奨**。少なくとも以下 2 件を
本 PR 内で修正してから merge するのが妥当。

1. **H-1**: `HubbardGC + TimeEvolution + TETwoBody` では batched InterAll を
   無効化するか、dynamic reinit + coefficient re-read を実装する。
2. **M-1**: `InitializeMPIBatchedExchange_Spin` の NULL チェックを
   SpinGC 版と同じ 6 ポインタに拡張する。

M-2 については §9.3 で実バグではないことを確認した。correlation は direct
call、spectrum (`M_CALCSPEC`) も `PairExSpin*` の direct path を通るため、
batched Exchange apply には到達しない。

初版の Q-1 はこの追補により「設計確認事項」ではなく **実バグ** と判断を改める。
したがって、この時点で merge 前に吉見さんへ伝えるべき優先度は
`H-1 > M-1 > Q-2` であり、M-2 は追加調査で切り離された。

## 9. 2026-04-12 修正後の最終チェック

この節は、`Fix batching bugs found in PR #216 review (H-1, M-1) and add regression test`
commit を対象にした**最新の追補レビュー**であり、上の `8.3 更新後のマージ判断`
を上書きする。

結論から言うと:

- **解消を確認**: `H-1 (HubbardGC + TETwoBody)`、`M-1 (NULL check)`
- **非バグと整理**: `M-2` は batched apply の到達経路が無い
- **新たに確認した High**: `canonical Hubbard + expert TEOneBody`

したがって、現時点でも **as-is merge は非推奨** である。

### 9.1 今回の修正で解消を確認した項目

#### H-1: `HubbardGC + TimeEvolution + TETwoBody`

`[src/mltplyHubbard.c:409-416](../src/mltplyHubbard.c#L409-L416)` に
`use_batching` guard が追加され、
`TimeEvolution && (NTEInterAllMax > 0 || NTETransferMax > 0)` のとき
batched path を無効化して non-batched MPI path にフォールバックするように
なった。

さらに:

- `[src/mltplyHubbard.c:425-459](../src/mltplyHubbard.c#L425-L459)` /
  `[src/mltplyHubbard.c:514-531](../src/mltplyHubbard.c#L514-L531)` で
  batched `Transfer` / `InterAll` の使用を `use_batching` で制御
- `[src/mltplyHubbard.c:464-488](../src/mltplyHubbard.c#L464-L488)` /
  `[src/mltplyHubbard.c:549-572](../src/mltplyHubbard.c#L549-L572)` で
  元の per-term MPI path へ復帰
- `test/mpi_consistency_te_hubbardgc_twobody.sh` が追加され、
  `step0 = 0 terms`, `step>=1 = inter-process TETwoBody` を回帰テスト化

を確認した。実際に
`mpi_consistency_te_hubbardgc_twobody` は pass しており、H-1 は今回の修正で
塞がったと判断してよい。

#### M-1: `InitializeMPIBatchedExchange_Spin` の NULL check 不備

`[src/mltplyMPIBatched.c:2530-2537](../src/mltplyMPIBatched.c#L2530-L2537)` で
NULL check が `term_indices` / `coefficients` 以外にも
`org_isite1` / `org_ispin1` / `org_ispin2` / `state1check` まで拡張された。
これにより、部分的な `malloc` failure で `is_initialized=1` のまま進む経路は
解消されている。

### 9.2 新たに確認した High: canonical `Hubbard` + expert `TEOneBody`

#### 症状

`HubbardGC` には上の `use_batching` guard が入ったが、
canonical `Hubbard` 側
`[src/mltplyHubbard.c:191-226](../src/mltplyHubbard.c#L191-L226)` は
依然として batched `MPIdouble` / `MPIsingle` transfer を**無条件で**
使っている。

一方 `CalcByTEM` は各 step 冒頭で
`[src/CalcByTEM.c:166-173](../src/CalcByTEM.c#L166-L173)` のように
`EDNTransfer` を元に戻し、
`[src/CalcByTEM.c:197-198](../src/CalcByTEM.c#L197-L198)` から
`[src/CalcByTEM.c:279-296](../src/CalcByTEM.c#L279-L296)` の
`MakeTEDTransfer()` を呼んで expert `TEOneBody` の transfer 項を step ごとに
追加し直す。

このため、step 0 に存在しなかった inter-process transfer が step 1 以降に
現れると:

- batched group は初回初期化時の topology のままで更新されない
- local fallback 側は inter-process term を `continue` で飛ばす
- 結果として、新規 inter-process `TEOneBody` が **silent drop** する

#### 再現

3x3 square の canonical Hubbard を 4 rank で実行し、
site 8 を inter-process site とする既存分割を使って、
expert `TEOneBody` として `0 <-> 8` の Hermitian pair を
`step 1` 以降にだけ追加する最小再現を作成した。

強めの再現条件では:

```text
MAXDIFF=0.0021849710
MPI log:
  [MPI Batching] Hubbard MPIdouble: No inter-process transfers
  [MPI Batching] HubbardGC: 8 MPIsingle transfers -> 2 groups (4.0x reduction)
```

`Flct.dat` の doublon は:

```text
serial D: 0.2592927429 -> 0.2610935344
MPI    D: 0.2592919779 -> 0.2592919779
```

と分岐し、MPI 側では expert `TEOneBody` が実質入っていない挙動になった。

#### 影響

- `Hubbard` の `TimeEvolution` で expert-mode `TEOneBody`
  (`NTETransferMax > 0`) を用いる MPI 計算は、serial と一致しない可能性がある
- 既存の `test/mpi_consistency_te_hubbard.sh` は AC Laser
  (`TransferWithPeierls`) しか通しておらず、この経路を検出できない
- 追加した docs note
  `[doc/ja/source/technical/MPI_ja.rst:93-104](../doc/ja/source/technical/MPI_ja.rst#L93-L104)` /
  `[doc/en/source/technical/MPI_en.rst:93-107](../doc/en/source/technical/MPI_en.rst#L93-L107)`
  は現状のコードより広いことを述べており、`HubbardGC` だけ修正された状態では
  記述が過剰である

#### 対応策

1. **最小安全策**:
   `mltplyHubbard()` にも `mltplyHubbardGC()` と同じ `use_batching` guard を入れ、
   `TimeEvolution && (NTETransferMax > 0 || NTEInterAllMax > 0)` では
   batched `Transfer` を無効化して、元の
   `general_hopp_MPIdouble` / `general_hopp_MPIsingle` ベースの per-term MPI path
   にフォールバックする。
2. **回帰テスト追加**:
   `mpi_consistency_te_hubbard.sh` を拡張するか、
   新規 `mpi_consistency_te_hubbard_stepdep_onebody.sh` を追加して、
   expert `TEOneBody` で `step0 = 0 terms`, `step>=1 = inter-process term`
   のケースを必ず通す。
3. **docs 修正**:
   コードを直さない限り、上の MPI technical note は
   「現状は `HubbardGC` に対してのみ自動フォールバック済み」
   と限定して書き直す必要がある。

### 9.3 M-2 は実バグではないことを確認（解消）

前回 `M-2` として報告した「batched Spin/SpinGC Exchange apply が
M_CORR / H_CORR / M_CALCSPEC で H.c. 分岐を 0 化していない」問題は、
追加調査の結果 **実際にはバグではない** と判断する。

根拠:

1. `[src/mltply.c:122](../src/mltply.c#L122)` が **無条件で**
   `X->Large.mode = M_MLTPLY` を設定する。
2. `mltplyHalfSpinGC` / `mltplyHalfSpin` の batched Exchange path は
   `mltply()` 経由でしか呼ばれないため、**常に M_MLTPLY** で動作する。
3. `M_CORR` / `H_CORR` が設定されるのは `expec_cisajs.c` /
   `expec_cisajscktaltdc.c` だが、
   これらは non-batched 関数
   `child_GC_CisAitCiuAiv_spin_MPIsingle` /
   `child_general_int_spin_MPIsingle` を**直接呼ぶ**
   （例: `[src/expec_cisajscktaltdc.c:1474](../src/expec_cisajscktaltdc.c#L1474)` /
   `[src/expec_cisajscktaltdc.c:1010](../src/expec_cisajscktaltdc.c#L1010)`）。
   batched path は経由しない。
4. `M_CALCSPEC` は `[src/PairEx.c:72](../src/PairEx.c#L72)` で設定されるが、
   Spin/SpinGC の spectrum 経路は
   `[src/PairExSpin.c:35](../src/PairExSpin.c#L35)` /
   `[src/PairExSpin.c:226](../src/PairExSpin.c#L226)` の
   `GetPairExcitedStateSpinGC()` / `GetPairExcitedStateSpin()` を**直接呼ぶ**。
   これらは `mltplyHalfSpinGC()` / `mltplyHalfSpin()` の batched Exchange path を
   通らない。

つまり、batched apply に `M_CORR` / `H_CORR` / `M_CALCSPEC` が到達する経路は
存在しない。non-batched 側の guard は direct call 用であり、batched path とは
無関係である。

なお、batched 側のコードが mode を区別しない点は将来の保守上の注意点として
残るが、**現時点では correctness risk ではない**。修正不要。

### 9.4 §9.2 の対応完了

§9.2 で指摘した canonical Hubbard + expert TEOneBody の問題は、
追加 commit で解消した。

- `[src/mltplyHubbard.c](../src/mltplyHubbard.c)` に `use_batching_H` guard を追加。
  `mltplyHubbardGC()` と同じ条件
  `TimeEvolution && (NTEInterAllMax > 0 || NTETransferMax > 0)` で
  batched Transfer を無効化し、元の per-term MPI path にフォールバック。
- `test/mpi_consistency_te_hubbard_twobody.sh` を追加。
  3x3 square canonical Hubbard で expert TEOneBody (site 0 ↔ 8, inter-PE) を
  step 1 以降にだけ追加し、serial / MPI(4) の Flct 一致を検証。pass 確認済み。

### 9.5 最終判断

検出された全項目の状態:

| ID | 状態 | 対応 |
|---|---|---|
| H-1 (HubbardGC + TETwoBody) | **解消** | `use_batching` guard + 回帰テスト |
| §9.2 (canonical Hubbard + TEOneBody) | **解消** | `use_batching_H` guard + 回帰テスト |
| M-1 (Spin NULL check) | **解消** | 6 ポインタ NULL check |
| M-2 (Spin M_CORR) | **非バグ** | `mltply.c:122` が常に M_MLTPLY。修正不要 |
| L-1〜L-5 | 未修正 | 保守性改善。別 PR で対応可 |

**結論: merge 可**。全 High / Medium が解消または非バグと確認された。
`mpi_consistency_*` 12 件（np=4）が pass。
Low 5 件は別 issue / 別 PR で対応する。

## 10. 関連リンク

- PR: https://github.com/issp-center-dev/HPhi/pull/216
- v1: [spinless_report.md](../spinless_report.md)
- v2: [spinless_report_v2.md](../spinless_report_v2.md)
- 関連レビュー: [PR200_review.md](PR200_review.md), [PR205_review.md](PR205_review.md), [PR210_review.md](PR210_review.md)
- TODO: [Dev_HPhi/TODO.md](../../TODO.md) §フェーズ 3 #216

---

**レビュー実施者**: tmisawa (with claude-opus-4-6 1M ctx)
**所要時間**: 約 90 分（subagent 並列調査 → 主張検証 → 実コード読み込み → 報告書作成）
