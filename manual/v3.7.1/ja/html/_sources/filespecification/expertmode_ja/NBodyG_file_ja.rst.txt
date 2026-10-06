.. highlight:: none

.. _Subsec:nbodyg:

NBodyG指定ファイル
~~~~~~~~~~~~~~~~~~

N体グリーン関数の計算対象成分を指定します。1行で :math:`N` 個の因子を指定し、
以下の期待値を表します。

.. math::

   \left\langle \prod_{p=1}^{N}
   c_{i_p\sigma'_p}^{\dagger} c_{j_p\sigma_p} \right\rangle .

Spin/SpinGC の場合、各因子はフェルミオンの一体演算子ではなく同一サイト上の
局所スピン行列要素を表すため、\ :math:`i_p=j_p`\ を満たす必要があります。
以下にファイル例を記載します。

::

    ==========================
    NNBodyG 3
    ==========================
    ============NBodyG========
    ==========================
    1 0 0 1 0
    2 0 0 0 0 1 1 1 1
    3 0 0 0 0 1 1 1 1 2 0 2 0

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降:
   [int02] ([int03] [int04] [int05] [int06]) ...

   括弧で示した4整数の組を [int02] 回繰り返します。

パラメータ
^^^^^^^^^^

-  :math:`[`\ string01\ :math:`]`

   **形式 :** string型 (空白不可)

   **説明 :** N体グリーン関数成分総数のキーワード名を指定します(任意)。

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** N体グリーン関数成分の総数を指定します。

-  :math:`[`\ int02\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** この成分に含まれる因子数 :math:`N` を指定します。
   正の整数でなければなりません。

-  各因子の :math:`[`\ int03\ :math:`]`, :math:`[`\ int05\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** サイト番号を指定する整数。0以上\ ``Nsite``\ 未満で指定します。
   :math:`[`\ int03\ :math:`]`\ が\ :math:`i_p`\ 、
   :math:`[`\ int05\ :math:`]`\ が\ :math:`j_p`\ に対応します。

-  各因子の :math:`[`\ int04\ :math:`]`, :math:`[`\ int06\ :math:`]`

   **形式 :** int型 (空白不可)

   | **説明 :** スピンまたは局所状態を指定する整数。
   | Hubbard、tJ、Kondo およびそれらの GC/NConserved 版では
   | 0: アップスピン
   | 1: ダウンスピン
   | を表します。
   | Spin/SpinGC では OneBodyG/TwoBodyG と同じ局所状態インデックスを用います。
   | SpinlessFermion/SpinlessFermionGC では 0 のみ指定可能です
     (0以外を指定した場合はエラー終了します)。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  :math:`[`\ int01\ :math:`]`\ と定義されているN体グリーン関数成分の総数が異なる場合はエラー終了します。

-  サイト番号またはスピン番号が許容範囲外の場合はエラー終了します。

-  Spin/SpinGC では、すべての因子が\ :math:`i_p=j_p`\ を満たす必要があります。

-  Kondo/KondoGC/KondoNConserved では、局在スピンサイトを含む因子は
   :math:`i_p=j_p`\ を満たす必要があります。一般スピンのKondo局在スピンには対応していません。

-  NBodyG は SpinGC, Spin, HubbardGC, Hubbard, HubbardNConserved,
   SpinlessFermionGC, SpinlessFermion, tJGC, tJ, tJNConserved,
   KondoGC, Kondo, KondoNConserved で使用できます。

-  NBodyG は等時刻期待値を出力するための指定であり、CalcSpec の
   励起演算子指定には使用されません。CalcSpec で動的相関関数を計算する場合は、
   SingleExcitation または PairExcitation を使用してください。

-  HubbardNConserved, tJNConserved, KondoNConserved では、NBodyG と
   CalcSpec の同時指定はエラーになります。

-  出力値は :ref:`Subsec:nbodygdat` で説明するファイルに出力されます。

.. raw:: latex

   \newpage
