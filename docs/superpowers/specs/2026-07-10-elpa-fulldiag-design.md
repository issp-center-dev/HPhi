# ELPA による FullDiag マルチノード GPU 対応 — 設計文書

- 日付: 2026-07-10
- 対象: HPhi FullDiag (`CalcType 2`) の対角化バックエンド拡張
- ステータス: ユーザー承認済み設計（実装前）

## 1. 背景と目的

現行の FullDiag の GPU 対応は MAGMA (`magma_zheevd_m`) によるシングルノード
マルチ GPU のみで、マルチノード GPU 計算はできない（`src/lapack_diag.c:82-91`）。
一方 CPU 側には ScaLAPACK (`pzheev`) による分散対角化が存在する
（`src/matrixscalapack.c:200`）。

本設計では、ScaLAPACK と同じ 2 次元ブロックサイクリック分散を前提とする
マルチノード対応固有値ソルバー **ELPA**（GPU カーネル搭載）を導入し、
FullDiag をマルチノード GPU 対応にする。

### 目標

- ヒルベルト空間次元 **N ~ 10^5** の全対角化を複数ノード（CPU/GPU）で実行可能にする。
- 汎用実装とし特定環境に依存しない（環境固有設定は `config/*.cmake` と
  ドキュメントに分離）。
- 検証環境: clavius（シングルノードマルチ GPU）、bowmore（CPU マルチノード）、
  ISSP スパコン ohtaka/kugui（本番）。

### 非スコープ

- 期待値計算 (`expec_*`) の分散化。固有ベクトルは従来どおり 1 本ずつ rank 0 の
  `v0` に集約してから期待値計算する（`src/phys.c` の既存構造を維持）。
- 既存経路（LAPACK / ScaLAPACK `pzheev` / MAGMA）の機能変更。
- Lanczos / LOBCG / TPQ 等、FullDiag 以外の計算タイプ。

## 2. ユーザーインターフェース

`calcmod.def` に新キーワード **`Solver`**（int 型）を追加する。

| 値 | バックエンド | 必要ビルドフラグ | 備考 |
|---|---|---|---|
| 0 | LAPACK `zheev`（逐次） | なし | 既定値 |
| 1 | ScaLAPACK `pzheev` | `_SCALAPACK` | 既存経路 |
| 2 | MAGMA `magma_zheevd_m` | `_MAGMA` | 既存経路、シングルノード限定 |
| 3 | ELPA（2stage） | `_SCALAPACK` + `_ELPA` | 本設計で追加 |

- `NGPU` は「ノードあたり GPU 枚数」として存続する。
  - `Solver 3` のとき: `NGPU 0` = ELPA を CPU 実行、`NGPU >= 1` = GPU 実行。
  - `Solver 2` のとき: 従来どおり MAGMA の GPU 枚数。
- 後方互換マッピング（`Solver` 未指定時、既存入力の挙動を不変に保つ。
  上から順に評価し、最初に一致したものを採用）:
  1. MAGMA ビルドかつ `NGPU > 0`（既定値 2）→ `Solver 2`
     （現行 `lapack_diag.c` は `iNGPU > 0` が ScaLAPACK 指定より優先のため）
  2. `ScaLAPACK 1` 指定 → `Solver 1`
  3. それ以外 → `Solver 0`
- `struct.h` の `DefineList` に `int iSolver` を追加。

### 起動時検証（`readdef.c`）

- `Solver` の値域チェック（0–3）。
- `Solver 1` かつ非 `_SCALAPACK` ビルド → エラー。
- `Solver 2` かつ非 `_MAGMA` ビルド → エラー。
- `Solver 3` かつ非 `_ELPA` ビルド → エラー。
- `Solver 3` かつ `nproc > 1` かつ `OutputHam`/`InputHam` 指定 → エラー
  （フェーズ 2 の分散生成と非互換のため）。
- `NGPU < 0` → エラー（既存チェックを流用）。
- `Solver 3` は `nproc == 1` でも許可する（clavius での単一ランク GPU 検証用）。

## 3. アーキテクチャ

新モジュール **`src/matrixlapack_elpa.c`**（+ `src/include/matrixlapack_elpa.h`）
に ELPA 呼び出しを隔離する。ELPA **C API**（`elpa.h`）を使用し、
**ELPA >= 2021.11** を要求する。ソルバーは `ELPA_SOLVER_2STAGE`、
GPU は `elpa_set("nvidia-gpu", 1)` で有効化する。

実装は 2 フェーズに分割し、各フェーズを個別に検証・マージする。

### フェーズ 1: ソルバー接続

入力行列は当面既存の全複製 `Ham` のまま、対角化のみ ELPA 化する。

- `lapack_diag.c` に `Solver == 3` 分岐を追加:
  1. BLACS グリッド生成（既存 ScaLAPACK 経路と同様、`MPI_Dims_create` で
     nprow × npcol）。
  2. ブロックサイズは **固定値 64**（ELPA 推奨域。現行 `GetBlockSize` は使わない）。
  3. 全複製 `Ham` から `pzelset_`（既存 `DivMat`）で 2D ブロックサイクリック
     分散行列 `A_distr` を構成。
  4. `diag_elpa_cmp()` を呼ぶ。
- `diag_elpa_cmp(xNsize, A_distr, desca, r, Z_distr, descZ, ngpu)`:
  - `elpa_init` → `elpa_allocate` → `elpa_set`（`na`, `nev = na`,
    `local_nrows`, `local_ncols`, `nblk`, `mpi_comm_parent`,
    `process_row`, `process_col`）→ `elpa_setup`。
  - `ngpu >= 1` のとき `nvidia-gpu = 1` を設定（枚数の明示割当は行わず、
    ランク→GPU 対応は ELPA/CUDA ランタイムに委ねる。ノードあたりランク数を
    GPU 枚数に合わせる運用をドキュメントに記載）。
  - `elpa_eigenvectors_double_complex` を実行。固有値 `r` は全ランクへ、
    固有ベクトルは分散行列 `Z_distr` のまま返す。
  - `elpa_deallocate` → `elpa_uninit`。
- **固有ベクトル回収の置換**: 新関数 `GetEigenVectorBlock(i, m, Z, descZ, vec)`
  を `matrixscalapack.c` に追加。固有状態 1 本（分散行列 Z の第 i 列）を
  `pzgemr2d_` で rank 0 に転送する。`phys.c` の分散経路
  （`use_scalapack` フラグ）でこれを使用する。既存の要素単位
  `GetEigenVector`（`pzelget_` を N 回/状態）は ELPA 経路では使わない。
  ScaLAPACK 経路 (`Solver 1`) も同関数に切り替える（同一の分散表現のため
  リスクなしで高速化される）。
- `use_scalapack` グローバルフラグは「固有ベクトルが分散格納されている」印
  として ELPA 経路でも 1 にセットする（`phys.c` の分岐を共用）。

### フェーズ 2: ハミルトニアンの分散生成（案 A: 1D 列パネル → 再分散）

全複製 `Ham`（16N² バイト/ランク）を排し、メモリを O(N²/P) に落とす。

- 根拠: `makeHam.c` の生成ループは「列 j を固定し、その列内の行に散布する」
  構造であり、全ての書き込みが `Ham[*][j]`（現在の列）で閉じている
  （`makeHam.c:92-249` で確認済み）。よって列方向の 1D ブロック分割が
  生成コード改修最小で成立する。
- `xsetmem.c`: `Solver == 3` かつ `nproc > 1` のとき、全複製 `Ham` / `L_vec`
  の代わりに 1D 列パネル `Ham_local`（N 行 × ceil(N/P) 列）のみ確保する。
- `makeHam.c`: 列ループを自ランク担当範囲 `[j_start, j_end)` に制限し、
  格納先を `Ham_local[i][j - j_start]` とする。物理項ごとの散布ロジックは
  変更しない。
- `diag_elpa_cmp()` の入口で 1D 列パネル → 2D ブロックサイクリックを
  `pzgemr2d_` 一回で再分散する（1D 側は nprow=1, npcol=P,
  列ブロックサイズ ceil(N/P) のディスクリプタで表現できるため、
  ScaLAPACK の標準機能で完結する）。
- 一時メモリ: パネルと分散行列の 2 枚が同時に存在する瞬間があり、
  ピークは約 2 × 16N²/P バイト/ランク + ELPA ワーク領域。
  N = 10^5、P = 32 ランクなら 1 ランクあたり約 20 GB で現実的。
- `Solver 0/1/2` および `nproc == 1` では従来どおり全複製 `Ham` を使う
  （フェーズ 2 の変更はビルドではなく実行時分岐）。

### データフロー（フェーズ 2 完成形、Solver 3・P ランク）

```
makeHam (各ランクが担当列を生成)
  → Ham_local [N × N/P, 1D 列パネル]
  → pzgemr2d (再分散)
  → A_distr [2D ブロックサイクリック, nblk=64]
  → ELPA 2stage (CPU or nvidia-gpu)
  → 固有値 (全ランク) + Z_distr (分散)
  → 固有状態ごとに GetEigenVectorBlock → rank 0 の v0
  → expec_* (既存の rank 0 集約計算)
```

## 4. エラー処理

- ELPA の全ステップ（`elpa_init` / `elpa_allocate` / `elpa_setup` /
  `elpa_set` / `elpa_eigenvectors`）の戻り値を検査し、`ELPA_OK` 以外なら
  確保済みリソースを解放して `-1` を返し、HPhi 標準のエラー終了に乗せる。
- GPU 要求時（`NGPU >= 1`）に ELPA の GPU 設定が失敗した場合は、
  「`NGPU 0` で CPU 実行せよ」という明示メッセージを出して停止する。
  **黙って CPU にフォールバックしない**（性能の取り違え防止。
  現行 MAGMA 経路の警告つき LAPACK フォールバックとは方針を変える）。
- `Solver` × ビルドフラグ × `nproc` × `OutputHam`/`InputHam` の整合性は
  すべて `readdef.c` の起動時検証で弾く（計算開始後に落とさない）。

## 5. ビルドシステム

- `CMakeLists.txt` に `option(USE_ELPA "Use ELPA" OFF)` を追加。
- 検出は pkg-config（`pkg_check_modules(ELPA elpa)`）を第一とし、
  見つからない場合は `ELPA_ROOT` ヒントによるパス探索をフォールバックにする。
- 検出成功時に `-D_ELPA` を定義し、`include_directories` / リンクを設定。
- **依存関係: `USE_ELPA=ON` は `USE_SCALAPACK=ON` を要求する**。
  ELPA 自体が ScaLAPACK/BLACS に依存しており、かつ本設計の再分散
  （`pzgemr2d`）・詰め込み（`pzelset`）も ScaLAPACK を使うため。
  `USE_ELPA=ON` かつ `USE_SCALAPACK=OFF` は CMake の configure 時に
  明示エラーにする。
- ScaLAPACK は MKL 同梱のもので可（bowmore / ISSP は MKL 環境）。
- `config/` に ELPA 有効ビルドの設定例を追加（GPU 版はドキュメントで
  ELPA 側のビルド要件 — CUDA 対応 configure — を注記）。

## 6. テスト計画

### 正しさ（各フェーズ完了時）

- 既存 `test/` の FullDiag サンプル（小規模スピン鎖・Hubbard、N <= 10^3）で
  `Solver 0` と `Solver 3` の固有値（許容 1e-8）とグリーン関数出力を比較。
- ランク数 1 / 2 / 4、`NGPU 0 / 1 / 2` の組み合わせ。
- 検証環境の順序: clavius（CPU → GPU、フェーズ 1）→ bowmore
  （マルチノード CPU、フェーズ 2）→ ISSP（マルチノード GPU、最終）。

### メモリ（フェーズ 2）

- N ~ 10^4 級で 1 ランクあたり常駐メモリが O(N²/P) に落ちていることを確認
  （全複製時代の 16N² と比較）。

### 回帰

- `Solver 0/1/2` の既存経路が無変更で通ることを既存テストで確認。
- ELPA 非搭載環境ではビルド・テストとも自動スキップ（ScaLAPACK と同じ扱い）。

## 7. ドキュメント更新

- `doc/{ja,en}/source/filespecification/expertmode_*/CalcMod_file_*.rst`:
  `Solver` キーワード追加、`NGPU` の説明更新（ELPA 時の意味を追記）。
- 「HPhi はマルチノード GPU 計算に対応していません」の記述を改訂。
- ビルドドキュメントに ELPA の導入手順（CPU 版 / CUDA 版）を追加。

## 8. 実装フェーズとマージ計画

1. **フェーズ 1**（ソルバー接続 + `Solver` キーワード + CMake + 固有ベクトル
   ブロック回収）: 単体で価値があり（既存全複製のまま ELPA の速度向上）、
   小規模で正しさを確定させる。
2. **フェーズ 2**（1D 列パネル分散生成 + 再分散）: フェーズ 1 の検証完了後に
   着手。メモリスケーリングを解放し N ~ 10^5 を可能にする。

各フェーズは独立の実装計画・PR とする。
