.. highlight:: none

.. _Subsec:ham:

ham.dat
~~~~~~~

(FullDiagでのみ出力)
``CalcMod``\ ファイルで\ ``OutputHam=1``\ の場合に、HΦ内部で計算されたハミルトニアンをMatrixMarket形式で出力します。\ ``CalcMod``\ ファイルで\ ``InputHam=1``\ とすると、定義ファイル一式と本ファイルを読み込み、再計算することができます。以下にファイル例を記載します。

::

    %%%%MatrixMarket matrix coordinate complex hermitian
    28 28 56
    1 1 1.000000 0.000000
    2 1 0.500000 0.000000
    3 2 0.500000 0.000000
    4 3 0.500000 0.000000
    5 4 0.500000 0.000000
    6 5 0.500000 0.000000
    7 6 0.500000 0.000000
    7 7 1.000000 0.000000
        …

ファイル名
^^^^^^^^^^

-  ##\_ham.dat

##はModParaファイル内の[string02]で指定されるヘッダを表します。

ファイル形式
^^^^^^^^^^^^

-  1行：ヘッダ

-  2行：\ :math:`[`\ int01\ :math:`]`  :math:`[`\ int02\ :math:`]`  :math:`[`\ int03\ :math:`]`

-  3行-：\ :math:`[`\ int04\ :math:`]`  :math:`[`\ int05\ :math:`]`  :math:`[`\ double01\ :math:`]`  :math:`[`\ double02\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :** Hamiltonianの行数。

-  :math:`[`\ int02\ :math:`]`

   **形式 :** int型

   **説明 :** Hamiltonianの列数。

-  :math:`[`\ int03\ :math:`]`

   **形式 :** int型

   **説明 :** Hamiltonianの非零の要素数。

-  :math:`[`\ double01\ :math:`]`, :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   | **説明 :** Hamiltonianの値を表します。
   | :math:`[`\ double01\ :math:`]`\ が実部、\ :math:`[`\ double02\ :math:`]`\ が虚部を表します。

.. raw:: latex

   \newpage
