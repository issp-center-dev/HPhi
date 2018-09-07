.. highlight:: none

Exchange指定ファイル
~~~~~~~~~~~~~~~~~~~~

Exchangeカップリングをハミルトニアンに付け加えます
(:math:`S=1/2`\ の系でのみ使用可能)。 電子系の場合には

.. math:: H+=\sum_{i,j}J_{ij}^{\rm Ex} (c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{j \downarrow}^{\dagger}c_{i  \downarrow}+c_ {i \downarrow}^{\dagger}c_{j\downarrow}c_{j \uparrow}^{\dagger}c_{i  \uparrow})

が付け加えられ、スピン系の場合には

.. math:: H+=\sum_{i,j}J_{ij}^{\rm Ex} (S_i^+S_j^-+S_i^-S_j^+)

が付け加えられます。 **スピン系の**\ :math:`(S_i^+S_j^-+S_i^-S_j^+)`\ **を
電子系の演算子で書き直すと、** 
:math:`-(c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{j \downarrow}^{\dagger}c_{i  \downarrow}+c_ {i \downarrow}^{\dagger}c_{j\downarrow}c_{j \uparrow}^{\dagger}c_{i  \uparrow})`
**となることに注意して下さい。** 
以下にファイル例を記載します。

::

    ====================== 
    NExchange 6  
    ====================== 
    ========Exchange ====== 
    ====================== 
       0     1  0.50000
       1     2  0.50000
       2     3  0.50000
       3     4  0.50000
       4     5  0.50000
       5     0  0.50000

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降: [int02] [int03] [double01]

パラメータ
^^^^^^^^^^

-  :math:`[`\ string01\ :math:`]`

   **形式 :** string型 (空白不可)

   **説明 :**
   Exchangeカップリングの総数のキーワード名を指定します(任意)。

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** Exchangeカップリングの総数を指定します。

-  :math:`[`\ int02\ :math:`]`, :math:`[`\ int03\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上\ ``Nsite``\ 未満で指定します。

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型 (空白不可)

   **説明 :** :math:`J_{ij}^{\rm Ex}`\ を指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  成分が重複して指定された場合にはエラー終了します。

-  :math:`[`\ int01\ :math:`]`\ と定義されているExchangeカップリングの総数が異なる場合はエラー終了します。

-  :math:`[`\ int02\ :math:`]`-:math:`[`\ int03\ :math:`]`\ を指定する際、範囲外の整数を指定した場合はエラー終了します。

.. raw:: latex

   \newpage
