.. highlight:: none

.. _Sec:outputfile:

出力ファイル
============

本節では :math:`{\mathcal H}\Phi` が生成する出力ファイルについて説明します。

クイックリファレンス
--------------------

**チェック・検証ファイル**

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :align: left

   * - ファイル
     - 説明
   * - :doc:`CHECK_Chemi.dat <CHECK_Chemi_ja>`
     - 化学ポテンシャル項
   * - :doc:`CHECK_InterAll.dat <CHECK_InterAll_ja>`
     - 一般二体相互作用
   * - :doc:`CHECK_CoulombIntra.dat <CHECK_CoulombIntra_ja>`
     - オンサイトクーロン相互作用
   * - :doc:`CHECK_Hund.dat <CHECK_Hund_ja>`
     - フント結合項
   * - :doc:`CHECK_INTER_Sr.dat <CHECK_INTER_U_ja>`
     - サイト間相互作用
   * - :doc:`CHECK_Memory.dat <CHECK_Memory_ja>`
     - メモリ使用量の見積もり
   * - :doc:`WarningOnTransfer.dat <WarningOnTransfer_ja>`
     - 転送積分に関する警告

**計算時間・進捗ファイル**

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :align: left

   * - ファイル
     - 説明
   * - :doc:`CalcTimer.dat <CalcTimer_ja>`
     - 各計算ステップの詳細時間
   * - :doc:`TimeKeeper.dat <TimeKeeper_ja>`
     - 全体の計算時間ログ
   * - :doc:`sz_TimeKeeper.dat <sz_TimeKeeper_ja>`
     - Sz計算の時間
   * - :doc:`Time_CG_EigenVector.dat <Time_CG_EigenVector_ja>`
     - CG固有ベクトル計算時間

**Lanczos法出力**

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :align: left

   * - ファイル
     - 説明
   * - :doc:`zvo_energy.dat <energy_ja>`
     - Lanczos反復中のエネルギー収束
   * - :doc:`zvo_Lanczos_Step.dat <Lanczos_Step_ja>`
     - Lanczosステップ情報

**TPQ（有限温度）出力**

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :align: left

   * - ファイル
     - 説明
   * - :doc:`Time_TPQ_Step.dat <Time_TPQ_Step_ja>`
     - TPQステップ時間
   * - :doc:`Norm_rand*.dat / Norm_tpq.dat <Norm_rand_ja>`
     - TPQ状態のノルム
   * - :doc:`SS_rand*.dat / SS_tpq.dat <SS_rand_ja>`
     - TPQのスピン-スピン相関
   * - :doc:`Flct_rand*.dat / Flct_tpq.dat <Flct_rand_ja>`
     - TPQの揺らぎ

**時間発展出力**

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :align: left

   * - ファイル
     - 説明
   * - :doc:`Time_TE_Step.dat <Time_TE_Step_ja>`
     - 時間発展ステップ時間
   * - :doc:`Norm.dat <Norm_ja>`
     - 時間発展中のノルム
   * - :doc:`SS.dat <SS_ja>`
     - 時間発展中のスピン-スピン相関
   * - :doc:`Flct.dat <Flct_ja>`
     - 時間発展中の揺らぎ

**物理量**

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :align: left

   * - ファイル
     - 説明
   * - :doc:`zvo_Eigenvalue.dat <Eigenvalue_ja>`
     - 固有値
   * - :doc:`zvo_phys*.dat <phys_ja>`
     - 物理量（エネルギー、ダブロン、Sz等）
   * - :doc:`zvo_cisajs*.dat <cisajs_ja>`
     - 一体グリーン関数
   * - :doc:`zvo_AnomalousG*.dat <AnomalousG_ja>`
     - 異常ペアグリーン関数
   * - :doc:`zvo_cisajscktalt*.dat <cisajscktalt_ja>`
     - 二体グリーン関数
   * - :doc:`zvo_NBodyG*.dat <NBodyG_ja>`
     - N体グリーン関数

**波動関数・ベクトル**

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :align: left

   * - ファイル
     - 説明
   * - :doc:`zvo_Ham.dat <ham_ja>`
     - ハミルトニアン行列要素（小規模系）
   * - :doc:`zvo_eigenvec*.dat <eigenvec_ja>`
     - 固有ベクトルデータ
   * - :doc:`zvo_tmpvec*.dat <tmpvec_ja>`
     - リスタート用一時ベクトル
   * - :doc:`zvo_recalcvec*.dat <recalcvec_ja>`
     - 再計算ベクトル
   * - :doc:`zvo_excited*.dat <excitedvec_ja>`
     - 励起状態ベクトル

**動的グリーン関数**

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :align: left

   * - ファイル
     - 説明
   * - :doc:`zvo_DynamicalGreen.dat <DynamicalGreen_ja>`
     - 動的グリーン関数データ
   * - :doc:`zvo_TMcomponents.dat <TMcomponents_ja>`
     - 三重対角行列成分
   * - :doc:`residual.dat <residual_ja>`
     - スペクトル計算の残差

----

詳細仕様
--------

チェック・検証ファイル
^^^^^^^^^^^^^^^^^^^^^^

入力の検証とメモリ要件の確認のために生成されるファイルです。

.. toctree::
   :maxdepth: 1

   CHECK_Chemi_ja
   CHECK_InterAll_ja
   CHECK_CoulombIntra_ja
   CHECK_Hund_ja
   CHECK_INTER_U_ja
   CHECK_Memory_ja
   WarningOnTransfer_ja

計算時間・進捗ファイル
^^^^^^^^^^^^^^^^^^^^^^

計算の進捗と時間を記録するファイルです。

.. toctree::
   :maxdepth: 1

   CalcTimer_ja
   TimeKeeper_ja
   sz_TimeKeeper_ja
   Time_CG_EigenVector_ja

Lanczos法出力
^^^^^^^^^^^^^

Lanczos計算固有の出力ファイルです。

.. toctree::
   :maxdepth: 1

   energy_ja
   Lanczos_Step_ja

TPQ（有限温度）出力
^^^^^^^^^^^^^^^^^^^

熱的純粋量子（TPQ）状態計算の出力ファイルです。

.. toctree::
   :maxdepth: 1

   Time_TPQ_Step_ja
   Norm_rand_ja
   SS_rand_ja
   Flct_rand_ja

時間発展出力
^^^^^^^^^^^^

時間発展計算の出力ファイルです。

.. toctree::
   :maxdepth: 1

   Time_TE_Step_ja
   Norm_ja
   SS_ja
   Flct_ja

物理量
^^^^^^

計算された物理量を含むファイルです。

.. toctree::
   :maxdepth: 1

   Eigenvalue_ja
   phys_ja
   cisajs_ja
   AnomalousG_ja
   cisajscktalt_ja
   NBodyG_ja

波動関数・ベクトル
^^^^^^^^^^^^^^^^^^

波動関数とベクトルデータを含むファイルです。

.. toctree::
   :maxdepth: 1

   ham_ja
   eigenvec_ja
   tmpvec_ja
   recalcvec_ja
   excitedvec_ja

動的グリーン関数
^^^^^^^^^^^^^^^^

動的グリーン関数計算の出力ファイルです。

.. toctree::
   :maxdepth: 1

   DynamicalGreen_ja
   TMcomponents_ja
   residual_ja
