.. highlight:: none

.. _Ch:HowToExpert:

エキスパートモード用入力ファイル
================================

本節では :math:`{\mathcal H}\Phi` のエキスパートモードで使用する入力ファイル(\*def)について説明します。

クイックリファレンス
--------------------

以下の表はエキスパートモードの全入力ファイルの一覧です。

**基本設定ファイル**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65
   :align: left

   * - ファイル
     - 必須
     - 説明
   * - :doc:`List <List_file_for_the_input_files_ja>`
     - Yes
     - 入力ファイル名のリスト
   * - :doc:`CalcMod <CalcMod_file_ja>`
     - Yes
     - 計算モードの設定
   * - :doc:`ModPara <ModPara_file_ja>`
     - Yes
     - 基本パラメータ（サイト数、電子数、Lanczosステップ等）
   * - :doc:`LocSpin <LocSpin_file_ja>`
     - 近藤のみ
     - 局在スピンの位置

**ハミルトニアン定義**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65
   :align: left

   * - ファイル
     - 必須
     - 説明
   * - :doc:`Trans <Trans_file_ja>`
     - No
     - 一体項: :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}`
   * - :doc:`InterAll <InterAll_file_ja>`
     - No
     - 一般二体相互作用: :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}`
   * - :doc:`NBodyInterAll <NBodyInterAll_file_ja>`
     - No
     - 一般N体相互作用: :math:`\prod_p c_{i_p\sigma'_p}^{\dagger}c_{j_p\sigma_p}`
   * - :doc:`CoulombIntra <CoulombIntra_file_ja>`
     - No
     - オンサイトクーロン: :math:`n_{i\uparrow}n_{i\downarrow}`
   * - :doc:`CoulombInter <CoulombInter_file_ja>`
     - No
     - サイト間クーロン: :math:`n_i n_j`
   * - :doc:`Hund <Hund_file_ja>`
     - No
     - フント結合: :math:`n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow}`
   * - :doc:`PairHop <PairHop_file_ja>`
     - No
     - ペアホッピング: :math:`c_{i\uparrow}^{\dagger}c_{j\uparrow}c_{i\downarrow}^{\dagger}c_{j\downarrow}`
   * - :doc:`Exchange <Exchange_file_ja>`
     - No
     - 交換相互作用: :math:`c_{i\uparrow}^{\dagger}c_{j\uparrow}c_{j\downarrow}^{\dagger}c_{i\downarrow}`
   * - :doc:`Ising <Ising_file_ja>`
     - No
     - イジング相互作用: :math:`S_i^z S_j^z`
   * - :doc:`PairLift <PairLift_file_ja>`
     - No
     - ペアリフト: :math:`c_{i\uparrow}^{\dagger}c_{i\downarrow}c_{j\uparrow}^{\dagger}c_{j\downarrow}`

**出力指定**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65
   :align: left

   * - ファイル
     - 必須
     - 説明
   * - :doc:`OneBodyG <OneBodyG_file_ja>`
     - No
     - 一体グリーン関数: :math:`\langle c^{\dagger}_{i\sigma_1}c_{j\sigma_2}\rangle`
   * - :doc:`TwoBodyG <TwoBodyG_file_ja>`
     - No
     - 二体グリーン関数: :math:`\langle c^{\dagger}_{i\sigma_1}c_{j\sigma_2}c^{\dagger}_{k\sigma_3}c_{l\sigma_4}\rangle`
   * - :doc:`NBodyG <NBodyG_file_ja>`
     - No
     - N体グリーン関数: :math:`\left\langle \prod_p c^{\dagger}_{i_p\sigma'_p}c_{j_p\sigma_p}\right\rangle`

**スペクトル・時間発展**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65
   :align: left

   * - ファイル
     - 必須
     - 説明
   * - :doc:`SingleExcitation <SingleExcitation_file_ja>`
     - スペクトル
     - 動的グリーン関数用の一粒子励起演算子
   * - :doc:`PairExcitation <PairExcitation_file_ja>`
     - スペクトル
     - 動的グリーン関数用のペア励起演算子
   * - :doc:`SingleExcitationBra <SingleExcitationBra_file_ja>`
     - オプション
     - 非対角グリーン関数用のブラ側一粒子励起演算子
   * - :doc:`PairExcitationBra <PairExcitationBra_file_ja>`
     - オプション
     - 非対角グリーン関数用のブラ側ペア励起演算子
   * - :doc:`SpectrumVec <SpectrumVec_File_ja>`
     - スペクトル
     - スペクトル計算用入力ベクトル
   * - :doc:`OneBodyTE <OneBodyTE_File_ja>`
     - 時間発展
     - 時間依存一体項
   * - :doc:`TwoBodyTE <TwoBodyTE_File_ja>`
     - 時間発展
     - 時間依存二体相互作用

----

詳細仕様
--------

基本設定ファイル
^^^^^^^^^^^^^^^^

計算の基本パラメータを定義するファイルです。

.. toctree::
   :maxdepth: 1

   List_file_for_the_input_files_ja
   CalcMod_file_ja
   ModPara_file_ja
   LocSpin_file_ja

ハミルトニアン定義
^^^^^^^^^^^^^^^^^^

ハミルトニアンの各項を指定するファイルです。

**一体項:**

.. toctree::
   :maxdepth: 1

   Trans_file_ja

**二体相互作用:**

.. toctree::
   :maxdepth: 1

   InterAll_file_ja
   NBodyInterAll_file_ja
   CoulombIntra_file_ja
   CoulombInter_file_ja
   Hund_file_ja
   PairHop_file_ja
   Exchange_file_ja
   Ising_file_ja
   PairLift_file_ja

出力指定
^^^^^^^^

計算・出力する物理量を指定するファイルです。

.. toctree::
   :maxdepth: 1

   OneBodyG_file_ja
   TwoBodyG_file_ja
   NBodyG_file_ja

スペクトル・時間発展
^^^^^^^^^^^^^^^^^^^^

動的グリーン関数計算および時間発展計算に使用するファイルです。

.. toctree::
   :maxdepth: 1

   SingleExcitation_file_ja
   PairExcitation_file_ja
   SingleExcitationBra_file_ja
   PairExcitationBra_file_ja
   SpectrumVec_File_ja
   OneBodyTE_File_ja
   TwoBodyTE_File_ja
