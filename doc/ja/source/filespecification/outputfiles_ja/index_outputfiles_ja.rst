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

   * - ファイル
     - 説明
   * - CHECK_Chemi.dat
     - 化学ポテンシャル項
   * - CHECK_InterAll.dat
     - 一般二体相互作用
   * - CHECK_CoulombIntra.dat
     - オンサイトクーロン相互作用
   * - CHECK_Hund.dat
     - フント結合項
   * - CHECK_INTER_Sr.dat
     - サイト間相互作用
   * - CHECK_Memory.dat
     - メモリ使用量の見積もり
   * - WarningOnTransfer.dat
     - 転送積分に関する警告

**計算時間・進捗ファイル**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ファイル
     - 説明
   * - CalcTimer.dat
     - 各計算ステップの詳細時間
   * - TimeKeeper.dat
     - 全体の計算時間ログ
   * - sz_TimeKeeper.dat
     - Sz計算の時間
   * - Time_CG_EigenVector.dat
     - CG固有ベクトル計算時間

**Lanczos法出力**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ファイル
     - 説明
   * - zvo_energy.dat
     - Lanczos反復中のエネルギー収束
   * - zvo_Lanczos_Step.dat
     - Lanczosステップ情報

**TPQ（有限温度）出力**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ファイル
     - 説明
   * - Time_TPQ_Step.dat
     - TPQステップ時間
   * - Norm_rand*.dat
     - TPQ状態のノルム
   * - SS_rand*.dat
     - TPQのスピン-スピン相関
   * - Flct_rand*.dat
     - TPQの揺らぎ

**時間発展出力**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ファイル
     - 説明
   * - Time_TE_Step.dat
     - 時間発展ステップ時間
   * - Norm.dat
     - 時間発展中のノルム
   * - SS.dat
     - 時間発展中のスピン-スピン相関
   * - Flct.dat
     - 時間発展中の揺らぎ

**物理量**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ファイル
     - 説明
   * - zvo_Eigenvalue.dat
     - 固有値
   * - zvo_phys*.dat
     - 物理量（エネルギー、ダブロン、Sz等）
   * - zvo_cisajs*.dat
     - 一体グリーン関数
   * - zvo_cisajscktalt*.dat
     - 二体グリーン関数

**波動関数・ベクトル**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ファイル
     - 説明
   * - zvo_Ham.dat
     - ハミルトニアン行列要素（小規模系）
   * - zvo_eigenvec*.dat
     - 固有ベクトルデータ
   * - zvo_tmpvec*.dat
     - リスタート用一時ベクトル
   * - zvo_recalcvec*.dat
     - 再計算ベクトル
   * - zvo_excited*.dat
     - 励起状態ベクトル

**動的グリーン関数**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ファイル
     - 説明
   * - zvo_DynamicalGreen.dat
     - 動的グリーン関数データ
   * - zvo_TMcomponents.dat
     - 三重対角行列成分
   * - residual.dat
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
   cisajscktalt_ja

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
