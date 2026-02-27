.. highlight:: none

.. _Sec:troubleshooting:

トラブルシューティング
======================

本節では :math:`{\mathcal H}\Phi` の実行中に発生する可能性のあるエラーメッセージとその対処法について説明します。

クイックリファレンス
--------------------

**コマンドライン・ファイルエラー**

.. list-table::
   :header-rows: 1
   :widths: 40 30 30
   :align: left

   * - エラーメッセージ
     - 原因
     - 対処法
   * - ``argument count is incorrect``
     - コマンドライン引数の数が不正
     - ``HPhi -e namelist.def`` または ``HPhi -s stan.in`` の形式で実行
   * - ``Error: The file (ファイル名) does not exist.``
     - 指定したファイルが存在しない
     - ファイルパスとファイル名を確認
   * - ``ERROR: Open input file``
     - 入力ファイルを開けない
     - ファイルの存在と読み取り権限を確認

**キーワード・パラメータエラー**

.. list-table::
   :header-rows: 1
   :widths: 40 30 30
   :align: left

   * - エラーメッセージ
     - 原因
     - 対処法
   * - ``ERROR ! Unsupported Keyword !``
     - 存在しないキーワードを指定
     - マニュアルでサポートされているキーワードを確認
   * - ``ERROR ! Keyword (名前) is duplicated !``
     - 同じキーワードを2回指定
     - 重複したキーワードを1つに統一
   * - ``ERROR ! (名前) is NOT specified !``
     - 必須キーワードが指定されていない
     - 必須キーワードを追加

**モデル・ソルバーエラー**

.. list-table::
   :header-rows: 1
   :widths: 40 30 30
   :align: left

   * - エラーメッセージ
     - 原因
     - 対処法
   * - ``ERROR ! Unsupported Solver : (名前)``
     - サポートされていないソルバーを指定
     - Lanczos, TPQ, FullDiag, CG, TimeEvolution のいずれかを指定
   * - ``ERROR ! Unsupported Model : (名前)``
     - サポートされていないモデルを指定
     - Hubbard, Spin, Kondo, SpinlessFermion のいずれかを指定
   * - ``Sorry, this system is unsupported in the STANDARD MODE...``
     - モデルと格子の組み合わせが未サポート
     - エキスパートモードを使用するか、サポートされている組み合わせを選択

**メモリエラー**

.. list-table::
   :header-rows: 1
   :widths: 40 30 30
   :align: left

   * - エラーメッセージ
     - 原因
     - 対処法
   * - ``Malloc ERROR for XXX``
     - メモリ割り当て失敗
     - システムサイズを小さくするか、メモリの大きいマシンを使用
   * - ``There is not enough memory. Please reduce MPI process number.``
     - 1プロセスあたりのメモリ不足
     - MPIプロセス数を減らして1プロセスあたりのメモリを増やす

**MPIエラー**

.. list-table::
   :header-rows: 1
   :widths: 40 30 30
   :align: left

   * - エラーメッセージ
     - 原因
     - 対処法
   * - ``You can not use OpenMP and MPI at the same time.``
     - OpenMPとMPIの併用が制限されている
     - どちらか一方を使用
   * - ``ERROR: MPI process count is not a power of 2``
     - MPIプロセス数が2の冪乗でない
     - MPIプロセス数を1, 2, 4, 8, 16などに設定
   * - ``ERROR: Not enough MPI processes.``
     - MPIプロセス数が不足
     - プロセス数を増やす

**ハミルトニアン検証エラー（カノニカル集団）**

以下のエラーはカノニカル集団の計算で、サイト数・電子数・スピンの組み合わせが物理的に実現できない場合に発生します。

.. list-table::
   :header-rows: 1
   :widths: 50 50
   :align: left

   * - エラーメッセージ
     - 意味
   * - ``ERROR ! abs(2 * Sz) > nsite in Hubbard model !``
     - Szの絶対値がサイト数の半分を超えている
   * - ``ERROR ! Nelec > 2 * nsite in Hubbard model !``
     - 電子数がサイト数の2倍を超えている
   * - ``ERROR ! (nelec + 2 * Sz) % 2 != 0 in Hubbard model !``
     - 電子数とSzのパリティが合わない
   * - ``ERROR ! nelec <= nsite && 2 * |Sz| > nelec in Hubbard model !``
     - 電子数に対してSzが大きすぎる
   * - ``ERROR ! nelec > nsite && 2 * |Sz| > 2 * nsite - nelec in Hubbard model !``
     - 電子数が多い場合にSzが大きすぎる
   * - ``ERROR ! abs(2 * Sz) > nsite in Spin model !``
     - スピン模型でSzがサイト数の半分を超えている
   * - ``ERROR ! (nsite + 2 * Sz) % 2 != 0 in Spin model !``
     - スピン模型でサイト数とSzのパリティが合わない
   * - ``ERROR ! Nelec_cond / 2 + Nelec_loc > nsite in Kondo model !``
     - 近藤模型で電子数がサイト数を超えている
   * - ``ERROR ! (nelec_cond + nelec_loc + 2 * Sz) % 2 != 0 in Kondo model !``
     - 近藤模型で電子数とSzのパリティが合わない

**対処法:** 入力パラメータ ``Nsite``, ``Nelec``, ``2Sz`` の値を確認し、物理的に実現可能な組み合わせに修正してください。

**警告メッセージ**

.. list-table::
   :header-rows: 1
   :widths: 50 50
   :align: left

   * - メッセージ
     - 説明
   * - ``Check ! (名前) is SPECIFIED but will NOT be USED.``
     - 指定したパラメータが現在の計算では使用されない。不要な場合はコメントアウト
   * - ``(名前) = (値) ###### DEFAULT VALUE IS USED ######``
     - デフォルト値が使用されている（エラーではない）

----

よくあるエラーの詳細な対処法
----------------------------

メモリ不足エラー
^^^^^^^^^^^^^^^^

**症状:** ``Malloc ERROR`` が表示され、計算が開始できない。

**原因:** ヒルベルト空間の次元が大きすぎて必要なメモリが確保できない。

**対処法:**

1. **システムサイズを小さくする**

   サイト数を減らすことで、ヒルベルト空間の次元が指数関数的に小さくなります。

2. **MPIプロセス数を増やす**

   ``mpirun -np 16 HPhi ...`` のようにプロセス数を増やすと、ヒルベルト空間がプロセス間で分割されます。

3. **使用可能メモリを確認**

   ``CHECK_Memory.dat`` を確認して、必要なメモリ量を把握してください。

MPIプロセス数のエラー
^^^^^^^^^^^^^^^^^^^^^

**症状:** ``ERROR: MPI process count is not a power of 2``

**原因:** :math:`{\mathcal H}\Phi` ではMPIプロセス数は2の冪乗（1, 2, 4, 8, 16, 32, ...）である必要があります。

**対処法:**

.. code-block:: bash

   # 正しい例
   mpirun -np 4 HPhi -e namelist.def
   mpirun -np 8 HPhi -e namelist.def
   mpirun -np 16 HPhi -e namelist.def

   # 間違った例（3は2の冪乗ではない）
   mpirun -np 3 HPhi -e namelist.def  # エラー

ファイルが見つからないエラー
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**症状:** ``Error: The file (ファイル名) does not exist.``

**対処法:**

1. namelist.def 内のファイル名を確認
2. 相対パスの場合、実行ディレクトリからの相対位置を確認
3. ファイル名の大文字・小文字を確認（Linux/macOSは大文字・小文字を区別）

収束しない場合
^^^^^^^^^^^^^^

**症状:** Lanczos法やCG法で収束しない。

**対処法:**

1. **Lanczosステップ数を増やす**

   ModParaファイルの ``LanczosEps`` や ``Lanczos_max`` を調整

2. **初期ベクトルを変更**

   ``initial_iv`` パラメータで異なる乱数シードを試す

3. **縮退がある場合**

   縮退した固有値がある場合、収束が遅くなることがあります
