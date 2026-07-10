# ELPA による FullDiag マルチノード GPU 対応 — 設計文書

- 日付: 2026-07-10（v4: AI設計レビュー3巡 + ELPA公式マニュアル照合を反映）
- 対象: HPhi FullDiag (`CalcType 2`) の対角化バックエンド拡張
- ステータス: レビュー反映済み設計（実装前）
- 一次資料: ELPA Manual — User's Guide and Best Practices
  (https://elpa.mpcdf.mpg.de/userguide/elpa_userguide.pdf, Version 2026.02.001)。
  以下「ELPA マニュアル」。API 手順・GPU 割当・ソルバー選択の記述は
  本マニュアルに照合済み。

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
- NVIDIA 以外の GPU バックエンド（AMD/Intel）。ただし GPU 有効化オプション名
  （`"nvidia-gpu"`）は `matrixlapack_elpa.c` 内の定数 1 箇所に集約し、
  将来 `amd-gpu` / `intel-gpu` へ拡張可能にする（ELPA はいずれも
  同一 API 形態でサポート）。

### 認識済みの限界（承知の上で進める）

- N ~ 10^5 の全対角化は演算量 O(N^3)・固有ベクトル格納 O(N^2) であり、
  行列格納を分散しても計算時間と GPU メモリが先に律速になり得る。
  本設計の主目的は (a) メモリ的に実行可能にすること、(b) 対角化カーネルの
  GPU 加速であり、実効性能はフェーズごとのベンチマークゲートで確認する（§8）。
- フェーズ 1 は全複製 `Ham` を残すため N ~ 10^5 のメモリ挙動は検証できない。
  フェーズ 1 の目的は ELPA の正しさと対角化カーネル性能の確認に限定し、
  スケーリング結論はフェーズ 2 のベンチマークで出す。
- 固有ベクトルを 1 本ずつ rank 0 に集約する構造は O(N^2) の対 rank 0 通信を
  生み、大きな N では律速になり得る。これは「期待値計算は rank 0 集約」という
  既存構造（非スコープ）に由来する制約であり、`GetEigenVectorBlock` の
  インターフェースは複数列の一括転送（バッチ幅 B）に拡張可能な形にしておき、
  実測で問題になれば フォローアップで最適化する。

## 2. ユーザーインターフェース

`calcmod.def` に新キーワード **`Solver`**（int 型）を追加する。

| 値 | バックエンド | 必要ビルドフラグ | 備考 |
|---|---|---|---|
| 0 | LAPACK `zheev`（逐次） | なし | 既定値 |
| 1 | ScaLAPACK `pzheev` | `_SCALAPACK` | 既存経路 |
| 2 | MAGMA `magma_zheevd_m` | `_MAGMA` | 既存経路、シングルノード限定 |
| 3 | ELPA | `_SCALAPACK` + `_ELPA` | 本設計で追加 |

### `NGPU` の意味

`NGPU` は「**ノードあたりの使用 GPU 枚数**」として存続する。

- `Solver 2`（MAGMA）: 従来どおり単一ランクから使う GPU 枚数。
- `Solver 3`（ELPA）: `NGPU 0` = CPU 実行。`NGPU >= 1` = GPU 実行。
  ランク→GPU の割当は **ELPA の自動ラウンドロビン割当**に任せる
  （ELPA マニュアル §「GPU usage」: "By default, ELPA will automatically
  assign each MPI task to a certain GPU device in a round-robin fashion"）。
  HPhi からの `use_gpu_id` 手動指定は行わない。これによりスケジューラが
  `CUDA_VISIBLE_DEVICES` でランクごとに単一 GPU を見せる構成
  （Slurm `--gpus-per-task=1` 等）とも衝突しない。
  `NGPU` の数値（>= 1 の部分）は §「起動時検証」の整合性警告にのみ使う。
- MAGMA と異なり、ELPA では 1 ランクが複数 GPU を使うことはない
  （1 ランク 1 GPU）。「GPU 1 枚あたりの MPI ランク数は全 GPU で一定」
  という ELPA の要件（マニュアル同節）と併せてドキュメントに明記する。

### `NGPU` の既定値

`NGPU` の既定値の扱い（明示指定の有無をパーサで記録する）:

- **`Solver` 未指定（レガシー入力）**: 現行のコンパイル時既定値を維持する
  （`_MAGMA` ビルドで 2、それ以外で 0。`readdef.c:274-278`）。
  レガシー解決（下表）はこの値を使う。既存入力の挙動を厳密に保つため。
- **`Solver` 明示指定時**: `NGPU` 未指定なら既定値はソルバー依存 —
  `Solver 2` → 2、`Solver 0/1/3` → 0（ELPA の GPU は明示オプトイン）。

### 優先順位表（`Solver` × 旧キーワード × ビルドフラグ）

パース後に `readdef.c` 内の単一関数で解決する（実行時経路に互換ロジックを
分散させない）。

| `Solver` | `ScaLAPACK` | `NGPU`（解決後） | ビルド | 結果 |
|---|---|---|---|---|
| 指定あり | 任意 | 任意 | 対応ビルド | **`Solver` が常に優先**。矛盾する旧キーワード（例: `Solver 0` + `ScaLAPACK 1`）は警告を出して無視 |
| 指定あり | — | — | 非対応ビルド（例: `Solver 3` で `_ELPA` なし） | エラー終了 |
| なし | 任意 | >0（コンパイル時既定含む） | `_MAGMA` | `Solver 2`（現行 `lapack_diag.c` の `iNGPU>0` 優先を踏襲） |
| なし | 1 | 0 | `_SCALAPACK` | `Solver 1` |
| なし | 1 | 任意 | 非 `_SCALAPACK` | 現行同様無視（警告を追加） |
| なし | 0/なし | 0 | 任意 | `Solver 0` |

- 旧 `ScaLAPACK` キーワード検出時は「`Solver 1` を使え」という
  非推奨（deprecation）警告を常に出す（動作は維持）。

### 起動時検証（`readdef.c`）

警告・エラーの表示は **rank 0 のみ**が行う（既存 `stdoutMPI` 慣例に従い、
全ランクからの重複出力を避ける）。

- `Solver` の値域チェック（0–3）。
- `Solver 1` かつ非 `_SCALAPACK` ビルド → エラー。
- `Solver 2` かつ非 `_MAGMA` ビルド → エラー。
- `Solver 3` かつ非 `_ELPA` ビルド → エラー。
- `Solver 3` かつ `NGPU >= 1` かつ GPU 対応 ELPA なしビルド
  （§5 の `_ELPA_GPU` 未定義）→ エラー（「`NGPU 0` で CPU 実行せよ」を案内）。
- `Solver 3` かつ `nproc > 1` かつ `OutputHam`/`InputHam` 指定 → エラー
  （フェーズ 2 の分散生成と非互換。フェーズ 1 の間は全複製 `Ham` が存在する
  ため許可し、フェーズ 2 マージ時にこの検証を有効化する。分散 I/O 対応は
  恒久的に非スコープとし、需要があれば別設計とする）。
- `NGPU < 0` → エラー（既存チェックを流用）。
- `Solver 3` は `nproc == 1` でも許可する（clavius での単一ランク GPU 検証用）。
- `Solver 3` かつ `NGPU >= 1` のとき、ノード内ランク数
  （`MPI_Comm_split_type(MPI_COMM_TYPE_SHARED)` で取得）が `NGPU` の
  整数倍でない場合は rank 0 から警告（GPU 遊休または GPU あたりランク数
  不均一が起きる旨）。実行は継続する（テスト構成の柔軟性を優先）。

## 3. アーキテクチャ

新モジュール群に ELPA 呼び出しと分散ユーティリティを隔離する:

- **`src/matrixlapack_elpa.c`**（+ `src/include/matrixlapack_elpa.h`）:
  ELPA ハンドル管理と対角化のみを担う（再分散は含まない）。
- **`matrixscalapack.c` への追加**: 再分散ユーティリティ
  `RedistPanelToBlockCyclic()` と固有ベクトル回収 `GetEigenVectorBlock()`。
  ソルバーとデータ移動の責務を分離する。

実装は 2 フェーズに分割し、各フェーズを個別に検証・マージする。

### ELPA API 契約

ELPA マニュアル §2（Fortran/C 使用例）および §「GPU usage」の手順に従う。

- ELPA **C API**（`elpa.h`）。最低バージョン:
  - CPU 実行: **API バージョン 20211125**（ELPA 2021.11 以降）。
    `elpa_init(20211125)` の戻り値で非対応版を検出しエラー終了。
  - GPU 実行: **ELPA 2023.11.001 以降**（`elpa_setup_gpu` が必要。§5 で
    configure 時に検出し `_ELPA_GPU` を定義。無ければ GPU 実行は起動時
    エラー）。
- 呼び出し手順。**集団操作（`elpa_setup` / `elpa_setup_gpu` /
  `elpa_eigenvectors`）の直前には必ず `MPI_Allreduce`（0=成功/-1=失敗の符号化では `MPI_MIN`。`MPI_MAX` は部分失敗を隠すため不可）による
  エラーフラグ同期を置き**、どこかのランクの局所的失敗（allocate/set）で
  他ランクが集団操作に進んでハングする事態を防ぐ:
  1. `handle = elpa_allocate(&error)`
  2. 必須パラメータを `elpa_set`（**setup 前**）:
     `na = N`, `nev = N`,
     `local_nrows` / `local_ncols`（**2D グリッド・`nblk = 64` に対する
     `numroc_` の値**。フェーズ 2 の 1D パネルのブロックサイズ NC とは無関係）,
     `nblk = 64`,
     `mpi_comm_parent = MPI_Comm_c2f(MPI_COMM_WORLD)`（**Fortran ハンドル
     変換必須**）, `process_row = myrow`, `process_col = mycol`
     （BLACS グリッド値）
  3. ★エラー同期 → `error = elpa_setup(handle)`
  4. ランタイムオプション（**setup 後**。ELPA マニュアルの規定順）:
     - `solver`: GPU 実行時は `ELPA_SOLVER_1STAGE`、CPU 実行時は
       `ELPA_SOLVER_2STAGE`（マニュアル: "ELPA1 is usually the better
       choice than ELPA2 for the performance on GPU"。局所行列が小さい場合の
       2stage 逆転はチューニング課題としてドキュメントに記載）。
     - GPU 実行時のみ: `elpa_set(handle, "nvidia-gpu", 1, &error)` —
       error が `ELPA_OK` でなければ「ELPA が CUDA 対応でビルドされて
       いない」と診断して停止（ELPA は未対応オプションにエラーを返す）。
  5. GPU 実行時のみ: ★エラー同期 → `status = elpa_setup_gpu(handle)`
     （マニュアル: "To finalize the GPU setup"。失敗時は GPU 環境不備として
     診断メッセージ付きで停止）。
  6. ★エラー同期 → `elpa_eigenvectors(handle, a, w, z, &error)`
     （型総称マクロ。ELPA C API の標準的な呼び方で、`a` が
     `double complex*` なら複素倍精度実体に解決される。入力 `a` は
     破壊される。出力固有ベクトル `z` は別配列）
  7. ★エラー同期 → `elpa_deallocate(handle, &error)` →
     `elpa_uninit(&error)`（`elpa_deallocate` も集団的に扱い、全ランクが
     同じ順序で呼ぶ）。
- **ライフサイクル**: `elpa_init` / `elpa_uninit` は `diag_elpa_cmp` 内で
  完結させる。FullDiag では `lapack_diag`（→ `diag_elpa_cmp`）は 1 回の
  実行で 1 度しか呼ばれないため、複数回 init/uninit の互換性問題は
  発生しない（将来複数回呼ぶ場合は init を初回のみに変える）。
- コミュニケータは `MPI_COMM_WORLD` 固定でよい（HPhi は FullDiag 経路で
  独自コミュニケータを使っておらず、既存 ScaLAPACK 経路も
  `MPI_COMM_WORLD` 前提）。
- `elpa_set` は全て int 値。各呼び出し直後に error を検査するヘルパ
  （設定名を含むメッセージ付き）を `matrixlapack_elpa.c` に置く。
- 複素 2stage カーネルの明示選択（`complex_kernel`）は行わない
  （GPU では 1stage を使い、CPU では ELPA の自動選択に任せる。
  必要になればチューニングオプションとして後日追加）。

### GPU 割当ポリシー

- HPhi は CUDA API を直接呼ばず、`use_gpu_id` も設定しない。
  ランク→GPU 対応は ELPA の自動ラウンドロビン（マニュアル明記）と、
  スケジューラによる `CUDA_VISIBLE_DEVICES` 束縛のいずれでも正しく動く。
- 黙った CPU フォールバックはしない: GPU 要求（`NGPU >= 1`）で
  `"nvidia-gpu"` 設定または `elpa_setup_gpu` が失敗したら停止し、
  「`NGPU 0` で CPU 実行せよ」と案内する。
- 運用上の推奨（ドキュメント記載）: **1 ランク 1 GPU を強く推奨**
  （ノードあたりランク数 = `NGPU`）。GPU あたり複数ランクは各ランクの
  ローカル行列＋ELPA ワークが同一デバイスに同居するため GPU OOM リスクが
  高い。GPU あたりランク数は全 GPU で一定にする（ELPA 要件）。
  クラスタはノード構成均一を前提とする。
- `NGPU` は使用 GPU 枚数を**物理的に制限しない**（割当は ELPA / スケジューラ
  側の可視デバイスで決まる）。枚数を厳密に制限したい場合はスケジューラの
  GPU 割当（`CUDA_VISIBLE_DEVICES` 等）で行う旨をドキュメントに明記する。

### フェーズ 1: ソルバー接続

入力行列は当面既存の全複製 `Ham` のまま、対角化のみ ELPA 化する。

- `lapack_diag.c` に `Solver == 3` 分岐を追加:
  1. BLACS グリッド生成（既存 ScaLAPACK 経路と同様、`MPI_Dims_create` で
     nprow × npcol）。
  2. ブロックサイズは **固定値 64**（ELPA 推奨域。現行 `GetBlockSize` は
     使わない。`pzgemr2d`/ELPA とも N が 64 で割り切れない場合を扱えるため
     制約はない）。
  3. 全複製 `Ham` から `pzelset_`（既存 `DivMat`）で 2D ブロックサイクリック
     分散行列 `A_distr` を構成（フェーズ 1 のみの暫定経路）。
     注意: フェーズ 1 では全複製 `Ham`（16N²/ランク）に加えて
     `A_distr`+`Z_distr` が同時に載るため、メモリはむしろ増える。
     フェーズ 1 の検証は小規模 N に限る（開発時の注意として明記）。
  4. `diag_elpa_cmp(xNsize, A_distr, desca, r, Z_distr, descZ, ngpu)` を呼ぶ。
     固有値 `r` は全ランクに、固有ベクトルは分散行列 `Z_distr` のまま返す。
     `A_distr` と `Z_distr` は別配列（ELPA が `a` を破壊するため）。
- **固有ベクトル回収** `GetEigenVectorBlock(i, m, Z, descZ, vec)`:
  固有状態 1 本（分散行列 Z の第 i+1 列、N×1 部分行列）を `pzgemr2d_` で
  rank 0 に転送する。
  - **ScaLAPACK 再分散ルーチンの規定パターンに従う**（p?gemr2d の仕様:
    転送に関与する全プロセスが呼び出し、宛先グリッドに属さないプロセスは
    宛先ディスクリプタの `CTXT_` を -1 とし、最終引数のコンテキストは
    両グリッドの全プロセスを包含するものを渡す）:
    - 宛先グリッド: **全ランクが** 同一の親（システム）BLACS コンテキストと
      同一の usermap（rank 0 のみを含む 1×1）を引数に `blacs_gridmap_` を
      呼ぶ（グリッド生成は集団操作のため。グリッドに含まれるのは rank 0
      のみ）。
    - 宛先ディスクリプタ `descV` は **全ランクで全フィールドを完全に初期化**
      する（`DTYPE_=1, M=N, N=1, MB=N, NB=1, RSRC=CSRC=0, LLD=N` を
      rank 0 と同値で設定）。その上で rank 0 以外は `descV[CTXT_] = -1` に
      上書きする（実装によっては CTXT 以外のフィールドも参照するため）。
    - C から呼ぶ全ての ScaLAPACK/BLACS ルーチン（`descinit_` /
      `pzgemr2d_` / `pzelset_` / `blacs_gridmap_` 等）は Fortran 呼出規約に
      従い、**スカラーを含む全引数をポインタ渡し**し、`descinit_` の末尾
      `info` のような出力引数も省略しない（既存 `matrixscalapack.c` の
      慣例に一致）。
    - 全ランクが `pzgemr2d_(N, 1, Z, 1, i+1, descZ, vec, 1, 1, descV,
      ictxt_2d)` を呼ぶ。`ictxt_2d`（全ランク参加の 2D コンテキスト）を
      最終引数に渡す。
    - コンテキストは phys ループ前に 1 度生成し、ループ後に
      `blacs_gridexit_` で明示的に解放する（固有状態ごとの生成・破棄は
      しない。2D コンテキスト等、本設計で生成する全 BLACS コンテキストも
      同様に使用終了時に解放する）。**回収ループ内では
      `MPI_Allreduce` エラー同期は行わない**（ディスクリプタ・バッファは
      ループ前に検証済みで、`pzgemr2d_` に局所失敗モードがないため。
      §4 のエラー同期は ELPA 呼び出し周りに限る）。
  - この呼び出しパターン自体を小規模単体テストで検証する（§6）。
    万一プラットフォーム固有の問題が出た場合の代替は、2D ディスクリプタ
    から所有ランクを計算する明示的 MPI 送受信（設計済みの代替案として
    記録のみ）。
  - `phys.c` の ELPA 経路でのみ使用する。**既存 ScaLAPACK 経路
    (`Solver 1`) の `GetEigenVector` は本件では変更しない**（共有経路の
    挙動変更を避ける）。`Solver 1` の回収高速化は、`GetEigenVectorBlock` が
    ELPA 経路で検証済みになった後の独立フォローアップとする。
  - インターフェースは列範囲 `[i, i+B)` の一括転送に拡張可能な形
    （実装は当面 B=1）にし、将来の通信レイテンシ最適化に備える。
- `use_scalapack` グローバルフラグは「固有ベクトルが分散格納されている」印
  として ELPA 経路でも 1 にセットする（`phys.c` の分岐を共用）。ELPA 経路か
  否かは `iSolver` で判別し、回収関数を切り替える。

### フェーズ 2: ハミルトニアンの分散生成（1D 列パネル → 再分散）

全複製 `Ham`（16N² バイト/ランク）を排し、行列格納を O(N²/P) に落とす。

- 根拠: `makeHam.c` の生成ループは「列 j を固定し、その列内の行に散布する」
  構造であり、全ての書き込みが `Ham[*][j]`（現在の列）で閉じている
  （`makeHam.c:92-249` で全書き込み箇所を確認済み。行側への鏡像書き込みや
  生成後の対称化は存在しない）。
- **書き込みの不変条件を格納抽象で強制する**: フェーズ 2 では `Ham` への
  書き込みを格納マクロ/インライン関数（例 `SetHamElem(i, j, dmv)`）経由に
  統一し、デバッグビルドでは「j が自ランク所有列である」ことをアサートする。
  将来の項追加が列所有を破った場合に即座に検出できる。
- **1D 列パネルの正確な定義**（P = `nproc`）:
  - 専用 BLACS コンテキスト `ictxt_1d` を 1×P グリッド（行優先 'R'、
    全 P ランク参加）で作る。
  - `NC = iceil(N, P)`（= ceil(N/P)）とし、ソースディスクリプタは
    `descinit_(desc1d, M=N, N=N, MB=N, NB=NC, RSRC=0, CSRC=0,
    ictxt=ictxt_1d, LLD=N)`。
  - この配置ではランク p（グリッド列 p）がグローバル列
    `[p*NC, min((p+1)*NC, N))` を所有する（NB=NC なので 1 ランク 1 ブロック。
    N が P で割り切れない場合、末尾ランクの所有列数は減る。所有列数は
    `numroc_(N, NC, mycol, 0, P)` と一致する）。
  - **所有列数 0 のランク**（P > ceil(N/NC) となる極小 N のとき）も
    ディスクリプタ上は有効とし、ローカル配列は最低 1 要素をダミー確保して
    非 NULL ポインタを渡す（`LLD` は常に N >= 1 で有効）。
  - グローバル列 j → 所有ランク `j / NC`、ローカル列番号 `j % NC`。
    この 2 式を上記格納マクロに実装する。
- `xsetmem.c`: `Solver == 3` かつ `nproc > 1` のとき、全複製 `Ham` / `L_vec`
  の代わりに 1D 列パネル `Ham_local`（N 行 × `numroc_` 列）のみ確保する。
  `L_vec` は確保しない（固有ベクトルは `Z_distr` に分散のまま）。
  フェーズ 2 実装時に `Ham` / `L_vec` グローバルの全参照箇所を grep 監査し、
  FullDiag 経路以外からの参照がないことを確認する。
- `makeHam.c`: 列ループを自ランク担当範囲 `[j_start, j_end)` に制限し、
  格納先を格納マクロ経由にする。物理項ごとの散布ロジックは変更しない。
- **`Ham_local` のメモリレイアウト**: ScaLAPACK/ELPA が要求する連続
  カラムメジャー 1 次元配列（`double complex*`）で確保し、格納マクロは
  `Ham_local[(j - j_start) * N + i]` に展開する（既存 `Ham` のような
  ポインタ配列 `double complex**` にはしない）。
- **再分散** `RedistPanelToBlockCyclic()`: `pzgemr2d_(N, N, Ham_local, 1, 1,
  desc1d, A_distr, 1, 1, desc2d, ictxt_2d)` の 1 回で行う。1D/2D とも同一の
  P ランク上のグリッドなので、包含コンテキストとして 2D 側を渡せる。
  再分散完了後、ELPA 呼び出し前に `Ham_local` を解放してピークメモリを
  下げる。
- `Solver 0/1/2` および `nproc == 1` では従来どおり全複製 `Ham` を使う
  （フェーズ 2 の変更はビルドではなく実行時分岐）。

### メモリバジェット（フェーズ 2 完成形、複素倍精度 16 バイト/要素）

| 領域 | サイズ/ランク | N=10^5, P=32 |
|---|---|---|
| `Ham_local`（1D パネル、再分散後に解放） | 16N²/P | 5 GB |
| `A_distr`（ELPA 入力、破壊される） | 16N²/P | 5 GB |
| `Z_distr`（固有ベクトル出力、`A_distr` と別配列） | 16N²/P | 5 GB |
| 固有値・回収バッファ等 | O(N) | ~2 MB |
| ELPA ワーク領域 | O(N²/P) 未満（実測で確認） | < 5 GB |

- ホストピーク: 再分散時 `Ham_local`+`A_distr` = 2 枚、対角化時
  `A_distr`+`Z_distr` = 2 枚 → **約 2×16N²/P = 10 GB/ランク**（P=32）。
- GPU（デバイス）側: ELPA はランクのローカル部分行列＋ワークを
  デバイスに置くため約 2–3×16N²/P。P=32 で 10–15 GB/ランクとなり
  A100 40GB 級で成立。**P（総ランク数）がデバイスメモリで決まる**。
- ホスト側の注意: コア数の多いノードでランク数を増やすと
  ランクあたり 10 GB × ランク数/ノード がノードの RAM を超え得る。
  実行ガイド（§7）に「必要総メモリ ≈ 2×16N² を全ノードの RAM 合計・
  全 GPU のデバイスメモリ合計がそれぞれ上回ること」という目安式と、
  ランク数/スレッド数のバランス指針を記載する。
- v1 の「2×16N²/P ≈ 20 GB」という記載は誤りで、本表に訂正
  （2 枚合計で 10 GB/ランク @ P=32）。

### データフロー（フェーズ 2 完成形、Solver 3・P ランク）

```
makeHam (各ランクが担当列を生成)
  → Ham_local [N × ~N/P, 1D 列パネル, ictxt_1d]
  → RedistPanelToBlockCyclic (pzgemr2d)  → Ham_local 解放
  → A_distr [2D ブロックサイクリック, nblk=64, ictxt_2d]
  → diag_elpa_cmp (ELPA 1stage(GPU)/2stage(CPU), GPU割当はELPA自動)
  → 固有値 (全ランク) + Z_distr (分散)
  → 固有状態ごとに GetEigenVectorBlock (pzgemr2d) → rank 0 の v0
  → expec_* (既存の rank 0 集約計算)
```

## 4. エラー処理

- ELPA の全ステップ（`elpa_init` / `elpa_allocate` / `elpa_set`（各回）/
  `elpa_setup` / `elpa_setup_gpu` / `elpa_eigenvectors`）の戻り値を検査する。
- **集団的エラー同期**: §3「ELPA API 契約」の手順どおり、**全ての集団操作
  （`elpa_setup`、`elpa_setup_gpu`、`elpa_eigenvectors`、および後続の
  `pzgemr2d_`）の直前**に `MPI_Allreduce`（0=成功/-1=失敗なら `MPI_MIN`）でエラーフラグを全ランクで
  共有する。どこかのランクの局所的失敗（`elpa_allocate` / `elpa_set` の
  失敗、メモリ確保失敗を含む）を全ランクが同期的に検知し、リソース解放 →
  `-1` を返して HPhi 標準のエラー終了（`exitMPI`）に乗せる。
  一部ランクだけが集団操作に進んでハングする事態を構造的に防ぐ。
- GPU 要求時（`NGPU >= 1`）に GPU 設定が失敗した場合は、
  「ELPA が GPU 対応でビルドされていない。`NGPU 0` で CPU 実行せよ」
  という明示メッセージを出して停止する。**黙って CPU にフォールバック
  しない**（性能の取り違え防止。現行 MAGMA 経路の警告つき LAPACK
  フォールバックとは方針を変える）。
- `Solver` × ビルドフラグ × `nproc` × `OutputHam`/`InputHam` の整合性は
  すべて `readdef.c` の起動時検証で弾く（計算開始後に落とさない）。

## 5. ビルドシステム

- `CMakeLists.txt` に `option(USE_ELPA "Use ELPA" OFF)` を追加。
- **依存関係: `USE_ELPA=ON` は ScaLAPACK を要求する**（ELPA 自体が
  ScaLAPACK/BLACS に依存し、再分散・回収も ScaLAPACK を使うため）。
  `USE_ELPA=ON` の場合は `USE_SCALAPACK` を自動的に ON にする。
  ユーザーが明示的に `USE_SCALAPACK=OFF` を指定していた場合のみ
  configure エラーにする。
- **検出（`cmake/FindELPA.cmake` を新規作成）**、優先順位順に:
  1. pkg-config: モジュール名は `elpa` に加え、バージョン付き
     （`elpa-<version>`）を探索して拾う。
  2. `ELPA_ROOT` ヒント: `find_path` + `file(GLOB)` で
     `${ELPA_ROOT}/include/elpa*-*/elpa/elpa.h` のバージョン付き
     ディレクトリを自動検出する（ライブラリ名は `elpa`）。
  3. 手動指定フォールバック: `ELPA_INCLUDE_DIR` / `ELPA_LIBRARY`
     キャッシュ変数の明示指定を常にサポートする（クラスタ固有の
     非標準配置対策）。
  4. いずれも失敗なら configure エラー（見つからない旨と `ELPA_ROOT` /
     手動変数の指定方法を表示）。
- 検出成功時に `-D_ELPA` を定義し、include/リンクを設定。さらに
  `CMAKE_REQUIRED_INCLUDES` / `CMAKE_REQUIRED_LIBRARIES` に検出済み
  ELPA を設定した上で `check_symbol_exists(elpa_setup_gpu elpa/elpa.h)`
  を実行し、成功なら `-D_ELPA_GPU` を定義する。これは **GPU 用 API
  （ELPA >= 2023.11.001）の存在検出**であり、リンク先 ELPA が実際に
  CUDA 対応でビルドされているかはランタイムの `"nvidia-gpu"` 設定
  エラー検査（§3）が担う。`_ELPA_GPU` 無しビルドでの `NGPU >= 1` は
  起動時エラー。
- **スレッド版の扱い**: **非スレッド版 `elpa` のみをサポートする**。
  `elpa_openmp` は HPhi が `MPI_Init`（= `MPI_THREAD_SINGLE`）で初期化
  しているため要件（`MPI_THREAD_SERIALIZED` 以上）を満たせず、リンクすると
  未定義動作になり得る。FindELPA は `elpa_openmp` しか見つからない場合、
  その旨と非スレッド版のビルド方法を示して configure エラーにする。
  スレッド版対応（`MPI_Init_thread` への移行）は将来の別課題とする。
- ScaLAPACK は MKL 同梱のもので可（bowmore / ISSP は MKL 環境）。
- `config/` に ELPA 有効ビルドの設定例を追加（GPU 版はドキュメントで
  ELPA 側のビルド要件 — CUDA 対応 configure — を注記）。

## 6. テスト計画

### 単体・部品テスト（フェーズごとに実装と同時）

- **再分散検証**: 決定的に生成した小行列（N が P・nblk=64 で割り切れない
  ケースを含む: 例 N=97, P=3）で、(a) 全複製→`pzelset` 経路と
  (b) 1D パネル→`pzgemr2d` 経路の 2D 分散行列を全要素比較（フェーズ 2）。
- **固有ベクトル回収検証**: `GetEigenVectorBlock` の rank 0 宛先グリッド
  パターンを、既知行列（対角行列等）の固有ベクトル回収で単体検証
  （フェーズ 1。所有列 0 ランクが生じる構成を含む）。
- **固有対の数学的検証**: 小規模エルミート行列で残差 `||A Z − Z D||` と
  直交性 `||Z^H Z − I||` を閾値検査し、LAPACK 結果と固有値を比較
  （グリーン関数比較だけでは位相・列順序の問題を見逃すため）。
  閾値は固定値でなくスケール則で定める: 残差は `c·N·ε·||A||`、直交性は
  `c·N·ε`（ε = 倍精度マシンイプシロン、c は O(10) の定数）。

### パーサ・後方互換テスト

- `Solver` × 旧 `ScaLAPACK`/`NGPU` × ビルドフラグの優先順位表（§2）の
  各行を calcmod.def の組み合わせテストで検証（エラー終了すべきものが
  エラーになること、非推奨警告が出ること、レガシー入力の解決結果が
  現行動作と一致することを含む）。

### 統合テスト（各フェーズ完了時）

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
- **GPU CI は存在しない**前提で、クラスタでの手動検証プロトコル
  （実行コマンド・期待値・チェックリスト）を `test/` 配下に文書として置き、
  リリース前チェックに組み込む。

## 7. ドキュメント更新

- `doc/{ja,en}/source/filespecification/expertmode_*/CalcMod_file_*.rst`:
  `Solver` キーワード追加、`NGPU` の説明更新（ELPA 時の意味 =
  ノードあたり GPU 枚数、割当は ELPA 自動ラウンドロビン、1 ランク 1 GPU、
  GPU あたりランク数一定の要件、MAGMA との挙動差を明記）。
  旧 `ScaLAPACK` キーワードに非推奨注記。
- 「HPhi はマルチノード GPU 計算に対応していません」の記述を改訂。
- ビルドドキュメントに ELPA の導入手順（CPU 版 / CUDA 版、
  GPU には ELPA >= 2023.11.001 が必要な旨）を追加。
- **実行ガイド**: ノードあたりランク数と `NGPU` の合わせ方、
  メモリ目安式（総メモリ ≈ 2×16N² をホスト RAM 合計・GPU デバイスメモリ
  合計が上回ること、GPU 側 2–3×16N²/P/ランク）、ランク/スレッド配分、
  ノード構成均一の前提、ISSP（kugui）向けジョブスクリプト例
  （srun/mpirun のバインディング指定込み）を ja/en 両方に追加。

## 8. 実装フェーズとマージ計画

1. **フェーズ 1**（ソルバー接続 + `Solver` キーワード + CMake/FindELPA +
   `GetEigenVectorBlock`（ELPA 経路のみ））: 単体で価値があり
   （既存全複製のまま ELPA の速度向上）、小規模で正しさを確定させる。
   ベンチマークゲート: clavius で `Solver 1` 比の対角化時間を記録。
2. **フェーズ 2**（1D 列パネル分散生成 + 再分散 + 格納マクロ/アサート +
   `Ham`/`L_vec` 参照監査）: フェーズ 1 の検証完了後に着手。メモリ
   スケーリングを解放し N ~ 10^5 を可能にする。ベンチマークゲート:
   bowmore マルチノードで N ~ 10^4 級の強スケーリングを記録。
3. **フォローアップ（別 PR）**: `Solver 1`（ScaLAPACK）経路の固有ベクトル
   回収を検証済み `GetEigenVectorBlock` に切り替え。必要に応じ回収の
   バッチ化（列範囲一括転送）。

各フェーズは独立の実装計画・PR とする。

## 9. レビュー対応記録

### v1 → v2（AI 設計レビュー 1 巡目）

| 指摘 | 対応 |
|---|---|
| ELPA C API 手順が不完全（solver 設定漏れ・`MPI_Comm_c2f`・エラー検査） | §3「ELPA API 契約」に手順を明文化 |
| GPU 設定契約が不十分（CUDA 非対応ビルド検出・カーネル選択） | §3 に検出・診断方針を明文化 |
| GPU バインディングが環境任せ | v2 で `use_gpu_id` 明示割当 → **v3 で撤回**（下記） |
| 1D→2D 再分散のディスクリプタ仕様不足 | §3 フェーズ 2 に descinit 引数・所有計算・非整除時挙動・検証テストを明記 |
| メモリ見積もりの不整合 | §3「メモリバジェット」表に訂正 |
| `Solver`×旧キーワードの優先順位が不完全 | §2 に優先順位表・`NGPU` 既定値規則・非推奨警告を追加 |
| CMake のバージョン付き ELPA パス検出 | §5 FindELPA.cmake の探索仕様を明記 |
| `Solver 1` の回収変更が共有経路リスク | フェーズ 1 から除外しフォローアップ PR に分離 |
| 1 ランク失敗時のハング | §4 集団的エラー同期を追加 |

### v2 → v3（AI 設計レビュー 2 巡目 + ELPA マニュアル照合）

| 指摘 | 対応 |
|---|---|
| Codex: solver/GPU オプションは setup **前**に設定すべき | **一次資料で棄却**: ELPA マニュアル §2 の使用例は必須パラメータのみ setup 前、solver/GPU/カーネルは setup **後**と明記。v3 はマニュアルの手順（setup 後設定 + `elpa_setup_gpu` で確定）を採用し、出典を明記 |
| （マニュアル照合で判明）`elpa_setup_gpu` 呼び出しが欠落 | §3 手順 5 に追加。GPU 実行の最低版を ELPA 2023.11.001 に設定、configure 時検出（`_ELPA_GPU`） |
| （マニュアル照合で判明）GPU では 1stage が通常高速 | ソルバー既定を GPU=1stage / CPU=2stage に変更 |
| Codex: rank-0 宛先 `pzgemr2d` の契約が未定義 | §3 に p?gemr2d 規定パターン（`blacs_gridmap_` による 1×1 グリッド、非参加ランクは `CTXT_=-1`、包含コンテキストを最終引数に）と単体テスト・代替案を明記 |
| Antigravity: setup 前の局所失敗で集団ハング | §3/§4: 全集団操作の直前にエラー同期を配置 |
| Codex: 優先順位表の「既定含む」が曖昧 | §2 `NGPU` 既定値: レガシー解決はコンパイル時既定（現行踏襲）、ソルバー依存既定は `Solver` 明示時のみ、と規定 |
| Codex: 所有列 0 ランクの割付規則 | §3 フェーズ 2 にダミー 1 要素確保を明記 |
| Codex: パーサ/後方互換テスト欠落 | §6 に追加 |
| Antigravity: `numroc_` のブロックサイズ曖昧 | §3 手順 2 で 2D・nblk=64 と明記 |
| Antigravity: 警告の全ランク重複出力 | §2 起動時検証を rank 0 のみに |
| Antigravity: `use_gpu_id` 強制はスケジューラ束縛（`CUDA_VISIBLE_DEVICES`）と衝突 | **v2 の明示割当を撤回**。ELPA の自動ラウンドロビン（マニュアル明記）に委ね、`use_gpu_id` は設定しない。スケジューラ束縛と両立 |
| Antigravity: 高コア数ノードでのホスト OOM | §3 メモリバジェット注記 + §7 実行ガイドに目安式 |
| Codex: elpa_openmp のスレッド制御 | §5 に非スレッド版優先 + スレッド指針を追加 |
| Codex: 列所有不変条件の恒久化 | §3 フェーズ 2 に格納マクロ + デバッグアサートを追加 |
| Codex/Antigravity: 回収の O(N²) 通信・バッチ化 | §1 認識済みの限界 + §8 フォローアップに記録（B 列一括転送に拡張可能な IF） |

### v3 → v4（AI 設計レビュー 3 巡目）

3 巡目で Codex は「must_fix なし・実装開始可」。Antigravity の must_fix
1 件と両者の should_fix / open_questions への対応:

| 指摘 | 対応 |
|---|---|
| Antigravity(must): `elpa_openmp` は `MPI_THREAD_SINGLE`（HPhi の `MPI_Init`）と非互換 | §5: 非スレッド版 `elpa` のみサポート。`elpa_openmp` のみ検出時は configure エラー。スレッド版対応は将来課題 |
| Codex: `elpa_deallocate`/`elpa_uninit` の集団性が曖昧 | §3 手順 7 で集団的に扱うことを明記、ライフサイクル（1 回呼び出し）も規定 |
| Codex: `blacs_gridmap_` は全ランクが呼ぶことを明確化 | §3 GetEigenVectorBlock を修正 |
| Codex: `_ELPA_GPU` は GPU-API 検出であって GPU 対応検出ではない | §5 の記述を訂正（実 CUDA 対応はランタイム検査が担う） |
| Antigravity: 非参加ランクの `descV` は全フィールド初期化すべき | §3 に反映（全初期化 + `CTXT_=-1` 上書き） |
| Antigravity: C からの `descinit_` は `info` 引数必須 | §3 に明記 |
| Antigravity: `check_symbol_exists` 前の `CMAKE_REQUIRED_*` 設定 | §5 に明記 |
| Antigravity: 回収ループ内の Allreduce は不要 | §3: 回収ループ内はエラー同期なしと明記 |
| Antigravity: `NGPU` は GPU 枚数を物理制限しない | §2/§3/§7 にドキュメント方針を明記 |
| Antigravity: GPU あたり複数ランクの GPU OOM リスク | §3 GPU 割当ポリシー: 1 ランク 1 GPU を強い推奨に格上げ |
| Codex: 検証閾値のスケール則 | §6 に `c·N·ε·||A||` 形式で規定 |
| Codex: CPU の 2stage 固定は妥当か | 固定で開始し、チューニングオプション化はフォローアップ（§8）— 実測が出るまで選択肢を増やさない |
| Codex: FindELPA のクラスタ固有問題 | §5 に `ELPA_INCLUDE_DIR`/`ELPA_LIBRARY` 手動フォールバックを追加、クラスタ検証は §6 の環境順テストで実施 |
| Antigravity: elpa_init ライフサイクル・独自コミュニケータ | §3 に明記（1 回呼び出しで完結、`MPI_COMM_WORLD` 固定で可） |

### v4 確定（AI 設計レビュー 4 巡目 = Antigravity 再確認）

4 巡目で **Codex（3 巡目）・Antigravity（4 巡目）とも must_fix ゼロ**。
最終 should_fix / open_questions への対応:

| 指摘 | 対応 |
|---|---|
| `Ham_local` はポインタ配列でなく連続カラムメジャー 1 次元配列に | §3 フェーズ 2 にレイアウトと格納マクロ展開式を明記 |
| BLACS コンテキストの明示解放 | §3: `blacs_gridexit_` による解放を明記 |
| Fortran 呼出規約（全引数ポインタ渡し）の一般化 | §3 に一般規則として明記 |
| `elpa_eigenvectors` 型総称マクロの採用 | §3 手順 6 を総称マクロに変更 |
| ノード不均一警告のロジック | rank 0 の属するノードの情報のみで判定する簡易版とし、ノード構成均一の前提を §3/§7 に明記（全ノード集約はしない） |
| CPU ソルバー（1stage/2stage）のパラメータ化 | 固定で開始、フォローアップ課題（§8）。実測が出るまで選択肢を増やさない |
| フェーズ 1 の一時的メモリ増（全複製+分散 2 枚） | §3 フェーズ 1 に開発時注意として明記（検証は小規模 N に限る） |
| GPU あたり複数ランクの OOM | §3/§7: 1 ランク 1 GPU の強い推奨と実行ガイド記載（対応済み、再掲） |
