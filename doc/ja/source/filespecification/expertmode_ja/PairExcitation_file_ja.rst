.. highlight:: none

.. _Subsec:pairexcitation:

PairExcitation指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~

動的関数

.. math:: G_n^{O_l,O_r}(z) = \langle \Phi_n | \hat{O}_l (z + E_n - \hat{\cal H})^{-1} \hat{O}_r| \Phi_n \rangle

において、\ :math:`\hat{O}_{l,r}`\ として二体励起状態を作成するための演算子

.. math:: \sum_{i, j, \sigma_1, \sigma_2} A_{i \sigma_1 j \sigma_2} c_{i \sigma_1}c_{j \sigma_2}^{\dagger} \quad \textrm{or} \quad
          \sum_{i, j, \sigma_1, \sigma_2} A_{i \sigma_1 j \sigma_2} c_{i\sigma_1}^{\dagger}c_{j\sigma_2}

を定義します。
ひとつの\ :math:`\hat{O}_r`\ と複数(1個以上)の\ :math:`\hat{O}_l`\ を指定することにより、効率よく計算を行うことが可能です。
なお、\ :math:`c_{i\sigma_1}c_{j\sigma_2}^{\dagger}`\ と\ :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}`\ を混同することは出来ません。
以下にファイル例を記載します。
この例では、

.. math::

    \hat{O}_r = \hat{S}_{\textbf{R}=\textbf{0}}^z = \frac{1}{2} (c_{\textbf{0}\uparrow}^{\dagger}c_{\textbf{0}\uparrow}-c_{\textbf{0}\downarrow}^{\dagger}c_{\textbf{0}\downarrow})
    \\
    \hat{O}_l = \hat{S}_{\textbf{R}}^z = \frac{1}{2} (c_{\textbf{R}\uparrow}^{\dagger}c_{\textbf{R}\uparrow}-c_{\textbf{R}\downarrow}^{\dagger}c_{\textbf{R}\downarrow})

としています。

::

    =============================================
    NPair 9
    =============================================
    =============== Pair Excitation =============
    =============================================
    2
    0 0 0 0 1        -0.500000000000000 0.0
    0 1 0 1 1         0.500000000000000 0.0
    2
    0 0 0 0 1        -0.500000000000000 0.0
    0 1 0 1 1         0.500000000000000 0.0
    2
    1 0 1 0 1        -0.500000000000000 0.0
    1 1 1 1 1         0.500000000000000 0.0
    2
    2 0 2 0 1        -0.500000000000000 0.0
    2 1 2 1 1         0.500000000000000 0.0
    2
    3 0 3 0 1        -0.500000000000000 0.0
    3 1 3 1 1         0.500000000000000 0.0
    :

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降:

    ::
      
        [int02]
        [int03]  [int04]  [int05]  [int06]  [int07]  [double01]  [double02]
        :
        [int02]個繰り返し

    というブロックを[int01]個繰り返す。
    1個目のブロックが\ :math:`\hat{O}_{r}`\ 、その後が\ :math:`\hat{O}_{l}`\ を表す。

パラメータ
^^^^^^^^^^

-  :math:`[`\ string01\ :math:`]`

   **形式 :** string型 (空白不可)

   **説明 :** 二体励起演算子のキーワード名を指定します(任意)。

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** 演算子\ :math:`\hat{O}_{r}`\ と\ :math:`\hat{O}_{l}`\ の数を合わせた総数を指定します。
   上の場合は1個の\ :math:`\hat{O}_{r}`\ と8個の\ :math:`\hat{O}_{l}`\ を合わせた9個となります。

-  :math:`[`\ int02\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** 各演算子\ :math:`\hat{O}_{r,l}`\ に含まれる項数を指定します。

-  :math:`[`\ int03\ :math:`]`, :math:`[`\ int05\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上\ ``Nsite``\ 未満で指定します。

-  :math:`[`\ int04\ :math:`]`, :math:`[`\ int06\ :math:`]`

   **形式 :** int型 (空白不可)

   | **説明 :** スピンを指定する整数。電子系・近藤格子系の場合は
   | 0: アップスピン、
   | 1: ダウンスピン、
   | スピン系の場合には、
   | :math:`0, 1, \cdots, 2S+1`
     (:math:`-S-0.5, -S+0.5, \cdots, S+0.5`\ に対応\ :math:`)`
   | を選択することが出来ます。

-  :math:`[`\ int07\ :math:`]`

   **形式 :** int型 (空白不可)

   | **説明 :** 二体励起演算子のタイプを指定する整数。
   | 0: :math:`c_{i\sigma_1}c_{j\sigma_2}^{\dagger}`
   | 1: :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}`
   | を選択することが出来ます。

-  :math:`[`\ double01\ :math:`]`, :math:`[`\ double02\ :math:`]`

   **形式 :** double型 (空白不可)

   **説明 :**
   :math:`A_{i \sigma_1 j \sigma_2}`\ の実部を\ :math:`[`\ double01\ :math:`]`\ 、虚部を\ :math:`[`\ double02\ :math:`]`\ でそれぞれ指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  成分が重複して指定された場合にはエラー終了します。

-  :math:`[`\ int01\ :math:`]`\ と定義されている二体励起演算子の総数が異なる場合はエラー終了します。

-  :math:`[`\ int02\ :math:`]`-:math:`[`\ int06\ :math:`]`\ を指定する際、範囲外の整数を指定した場合はエラー終了します。


.. raw:: latex

   \newpage
