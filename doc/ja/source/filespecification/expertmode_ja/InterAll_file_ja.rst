.. highlight:: none

.. _Subsec:interall:

InterAll指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~~

ここでは一般二体相互作用をハミルトニアンに付け加えます。付け加える項は以下で与えられます。

.. math::

   \mathcal{H}+=\sum_{i,j,k,l}\sum_{\sigma_1,\sigma_2, \sigma_3, \sigma_4}
   I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}


なお、スピンに関して計算する場合には、\ :math:`i=j, k=l`\ となるよう設定してください。
また、オフダイアゴナル成分については、i1 sigma1 i2 sigma2 i3 sigma3 i4 sigma4 とi4 sigma4 i3 sigma3 i2 sigma2 i1 sigma1 のペアを必ず用意するようにしてください。

フェルミオンの項は、上式に書かれた演算子の順序で解釈します。
:math:`\sigma\ne\tau` のスピン交換型の積を並べ替える際には、縮約による一体項も含まれます。

.. math::

   c^\dagger_{i\sigma}c_{j\tau}c^\dagger_{k\tau}c_{l\sigma}
   = \delta_{jk}c^\dagger_{i\sigma}c_{l\sigma}
     - c^\dagger_{i\sigma}c_{l\sigma}c^\dagger_{k\tau}c_{j\tau}.

:math:`j=k` かつ :math:`i\ne l` の場合、この一体項を内部で ``Transfer`` に加算します。
Transferのハミルトニアンには全体にマイナス符号が付くため、入力係数としては :math:`-I` に相当します。
元のInterAll項を残したまま、この縮約一体項を手動で重ねて追加しないでください。
``TETwoBody`` の各時刻にも同じ規則を適用します。

.. note::

   以前のバージョンでは、Hubbard・Kondo・tJ系でこのサイト間の縮約一体項が欠落していました。
   該当する計算は再計算が必要です。古いsolver checkpointから継続しても、過去の計算結果は補正されません。
   スピン保存のbilinearの積と縮約一体項へ既に書き換えた入力は、そのまま使用できます。

以下にファイル例を記載します。

::

    ====================== 
    NInterAll      36  
    ====================== 
    ========zInterAll===== 
    ====================== 
    0    0    0    1    1    1    1    0   0.50  0.0
    0    1    0    0    1    0    1    1   0.50  0.0
    0    0    0    0    1    0    1    0   0.25  0.0
    0    0    0    0    1    1    1    1  -0.25  0.0
    0    1    0    1    1    0    1    0  -0.25  0.0
    0    1    0    1    1    1    1    1   0.25  0.0
    2    0    2    1    3    1    3    0   0.50  0.0
    2    1    2    0    3    0    3    1   0.50  0.0
    2    0    2    0    3    0    3    0   0.25  0.0
    2    0    2    0    3    1    3    1  -0.25  0.0
    2    1    2    1    3    0    3    0  -0.25  0.0
    2    1    2    1    3    1    3    1   0.25  0.0
    4    0    4    1    5    1    5    0   0.50  0.0
    4    1    4    0    5    0    5    1   0.50  0.0
    4    0    4    0    5    0    5    0   0.25  0.0
    4    0    4    0    5    1    5    1  -0.25  0.0
    4    1    4    1    5    0    5    0  -0.25  0.0
    4    1    4    1    5    1    5    1   0.25  0.0
    …

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降:
   [int02] [int03] [int04] [int05] [int06] [int07] [int08] [int09] [double01] [double02]

パラメータ
^^^^^^^^^^

-  :math:`[`\ string01\ :math:`]`

   **形式 :** string型 (空白不可)

   **説明 :** 二体相互作用の総数のキーワード名を指定します(任意)。

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** 二体相互作用の総数を指定します。

-  :math:`[`\ int02\ :math:`]`, :math:`[`\ int04\ :math:`]`,
   :math:`[`\ int06\ :math:`]`, :math:`[`\ int08\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上\ ``Nsite``\ 未満で指定します。

-  :math:`[`\ int03\ :math:`]`, :math:`[`\ int05\ :math:`]`,
   :math:`[`\ int07\ :math:`]`, :math:`[`\ int09\ :math:`]`

   **形式 :** int型 (空白不可)

   | **説明 :** スピンを指定する整数。
   | 0: アップスピン
   | 1: ダウンスピン
   | を選択することが出来ます。

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型 (空白不可)

   **説明 :**
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}`\ の実部を指定します。

-  :math:`[`\ double02\ :math:`]`

   **形式 :** double型 (空白不可)

   **説明 :**
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}`\ の虚部を指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  Hamiltonianがエルミートという制限から\ :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}=I_{lkji\sigma_4\sigma_3\sigma_2\sigma_1}^{\dagger}`\ の関係を満たす必要があります。上記の関係が成立しない場合にはエラー終了します。
   また、エルミート共役の形式は\ :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}`\ に対して、\ :math:`I_{lkji\sigma_4\sigma_3\sigma_2\sigma_1}`
   :math:`c_{l\sigma_4}^{\dagger}c_{k\sigma_3}c_{j\sigma_2}^{\dagger}c_{i\sigma_1}`\ を満たすように入力してください。

-  スピンに関して計算する場合、\ :math:`i=j, k=l`\ を満たさないペアが存在するとエラー終了します。

-  成分が重複して指定された場合にはエラー終了します。

-  :math:`[`\ int01\ :math:`]`\ と定義されているInterAllの総数が異なる場合はエラー終了します。

-  :math:`[`\ int02\ :math:`]`-:math:`[`\ int09\ :math:`]`\ を指定する際、範囲外の整数を指定した場合はエラー終了します。

.. raw:: latex

   \newpage
