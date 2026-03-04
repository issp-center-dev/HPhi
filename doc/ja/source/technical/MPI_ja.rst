MPI並列化
=========

本節では :math:`{\mathcal H}\Phi` で使用されるMPI並列化戦略について説明します。

概要
----

:math:`{\mathcal H}\Phi` はMPIを使用して複数プロセス間で計算を並列化します。
ヒルベルト空間の基底状態はMPIプロセス間に分散され、
ハミルトニアンの行列ベクトル積 :math:`\hat{H}|\psi\rangle` が並列に実行されます。

サイトの分類
------------

:math:`N` サイト系で :math:`N_{\rm proc}` MPIプロセスを使用する場合、
サイトは2つのカテゴリに分類されます：

**ローカルサイト** (:math:`N_{\rm local}` サイト)
  量子状態が各プロセス内に格納されるサイト。
  これらのサイトに対する基底状態操作はローカルで実行可能。

**プロセス間サイト** (:math:`N_{\rm inter}` サイト)
  量子状態が複数プロセスにまたがるサイト。
  これらのサイトに対する操作にはMPI通信が必要。

関係式：

.. math::
   N_{\rm local} = N - \log_2(N_{\rm proc})

.. math::
   N_{\rm inter} = \log_2(N_{\rm proc})

例えば、8サイト系で4 MPIプロセス（:math:`N_{\rm proc}=4`）を使用する場合、
ローカルサイト :math:`N_{\rm local}=6`、プロセス間サイト :math:`N_{\rm inter}=2` となります。

MPI通信パターン
---------------

ハミルトニアン中の演算は、必要なMPI通信によって分類されます：

**ローカル演算**
  演算に関わる両サイトがローカルサイト。
  MPI通信は不要。

**MPIsingle演算**
  一方がローカルサイト、他方がプロセス間サイト。
  1つのパートナープロセスとのMPI通信が必要。
  通信相手（origin）は以下で決定：

  .. math::
     \text{origin} = \text{myrank} \oplus T_{\rm pow}[\text{site}]

  ここで :math:`\oplus` はXOR演算、:math:`T_{\rm pow}` は各サイトの2のべき乗です。

**MPIdouble演算**
  両サイトがプロセス間サイト。
  MPI通信が必要だが、受信データは一様に処理可能。
  通信相手は：

  .. math::
     \text{origin} = \text{myrank} \oplus (T_{\rm pow}[\text{site}_1] + T_{\rm pow}[\text{site}_2])

バッチ処理MPI通信
-----------------

MPI通信オーバーヘッドを削減するため、:math:`{\mathcal H}\Phi` は
同じ通信相手（origin）を持つ複数の演算をバッチにグループ化します：

**最適化前：**
  プロセス間サイトを持つ各転送項/相互作用項が個別に ``MPI_Sendrecv`` を呼び出す。
  同じoriginを持つ :math:`N` 項がある場合、:math:`N` 回の独立した通信が発生。

**最適化後：**
  同じoriginを持つ項をグループ化。
  1回の ``MPI_Sendrecv`` で必要なデータを取得し、
  その後 :math:`N` 項すべてをローカルで処理。

この最適化は以下の場合に特に効果的：

* 多くのホッピング項が同じ通信相手を共有する2D/3D格子モデル
* 隣接サイト間に多くの転送項を持つ系

バッチ通信は以下のモデルに実装済み：

* SpinlessFermionモデル（MPIsingle）
* Hubbardモデル（MPIsingle、MPIdouble、InterAll）
* Spinモデル（Exchange相互作用）

SpinlessFermionのoff-diagonal two-body Green関数
-----------------------------------------------

SpinlessFermionにおける
:math:`\langle c^\dagger_i c_j c^\dagger_k c_l \rangle`
のMPI経路は、グランドカノニカル（GC）に加えてカノニカルでも
プロセス間サイトを含む場合に対応しています。

この経路では、グローバルビット表現で演算子列を適用し、
ランク間通信を伴う場合でもフェルミオン符号を一貫して評価します。

MPI前提テストの挙動（ctest）
----------------------------

以下のテストはMPI実行環境を前提とします：

* ``mpi_consistency_*`` （目安: ``-np >= 2``）
* ``lanczos_*_mpidouble`` （目安: ``-np = 16``）

``MPIRUN`` が未設定、または必要なMPIランク条件を満たさない場合、
これらは ``ctest`` 上で ``Skipped`` として報告されます。

例：

.. code-block:: bash

   MPIRUN='mpirun -np 4 --oversubscribe' ctest -L consistency
   MPIRUN='mpirun -np 16 --oversubscribe' ctest -L batching

フェルミオン符号
----------------

フェルミオン系（Hubbard、SpinlessFermion）では、
演算子の反交換から生じるフェルミオン符号をMPI通信時に注意深く扱う必要があります。

符号は ``SgnBit`` 関数を使用して計算され、
生成演算子と消滅演算子の間の占有状態数をカウントします。
この符号は初期化時に事前計算され、各転送グループと共に格納されます。

プロセス数の要件
----------------

MPIプロセス数は特定の制約を満たす必要があります：

**Hubbard/近藤モデル**
  プロセス数 = :math:`4^n` （1サイトあたり4状態: 空、上、下、二重占有）

**スピン1/2モデル**
  プロセス数 = :math:`2^n` （1サイトあたり2状態: 上、下）

**一般スピンモデル**
  プロセス数 = プロセス間サイトの :math:`(2S_i+1)` の積

プロセス数設定の詳細は :ref:`Subsec:process` を参照してください。
