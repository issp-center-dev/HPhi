.. highlight:: none

.. _Subsec:pairexcitationbra:

PairExcitationBra指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

非対角動的グリーン関数

.. math::

   G_{BA}(z) = \langle \phi | B^{\dagger} \frac{1}{z-H} A | \phi \rangle

を計算するための\ **ブラ側**\ 二体励起演算子\ :math:`B`\ を定義します。
ケット側演算子\ :math:`A`\ は\ :doc:`PairExcitation <PairExcitation_file_ja>`\ 指定ファイルで定義します。
本ファイルを省略した場合は、通常の対角グリーン関数（\ :math:`B=A`\ ）が計算されます。

ファイルフォーマットおよびパラメータは\ :doc:`PairExcitation <PairExcitation_file_ja>`\ 指定ファイルと同一で、
演算子数のキーワードのみ ``NPairExcitationBra`` です。以下にファイル例（サイト 1 の密度演算子\ :math:`n_{1\uparrow}`\ ）を記載します。

::

    ===============================
    NPairExcitationBra    1
    ===============================
    ======= Pair Excitation Bra ====
    ===============================
        1     0     1     0     1    1.0    0.0

.. _use_rules_pairexcitationbra_ja:

使用ルール
^^^^^^^^^^

*  ヘッダの省略はできません。フォーマット・パラメータは\ :doc:`PairExcitation <PairExcitation_file_ja>`\ と同一です。

*  本ファイルに記載する演算子は\ :math:`B`\ であり\ :math:`B^{\dagger}`\ ではありません
   （計算される量は\ :math:`\langle\phi|B^{\dagger}(z-H)^{-1}A|\phi\rangle`\ です）。

*  非対角グリーン関数はシフト型 BiCG 法（\ ``method="CG"``\ ）かつ ``CalcSpec="Normal"`` のときのみ利用可能です
   （ブラ演算子を指定した場合、リスタート・保存系のスペクトルモードは利用できません）。

*  ケット（\ :doc:`PairExcitation <PairExcitation_file_ja>`\ ）とブラの演算子は、\ :math:`|\phi\rangle`\ を
   \ **同一**\ の励起ヒルベルト空間（粒子数・\ :math:`S_z`\ の変化が同じ）へ写す必要があります。
   1 つのファイル内の全演算子も同一のセクター変化をもつ必要があります。これを満たさない場合はプログラムが終了します。

*  対応モデル：Hubbard, Spin, SpinGC 模型。（その他の模型は非対角グリーン関数に未対応のため、
   これらの場合はプログラムが終了します。）``SingleExcitationBra`` と ``PairExcitationBra`` を同時に指定することはできません。

.. raw:: latex

   \newpage
