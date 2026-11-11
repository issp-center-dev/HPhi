.. highlight:: none

.. _Subsec:residual:

residual.dat
~~~~~~~~~~~~~

``method``\ =\ ``CG``\ (``CalcMod``\ ファイルで\ ``CalcType``\ が3の場合)を用いた動的グリーン関数の計算のステップ数に依存した途中結果が出力されます。
ファイル名およびファイル形式は以下の通りです。

ファイル名
^^^^^^^^^^

-  residual.dat

ファイル形式
^^^^^^^^^^^^

-  1行目-\ ``NOmega``\ 行目：\ :math:`[`\ int01\ :math:`]` 
   :math:`[`\ double01\ :math:`]`  :math:`[`\ double02\ :math:`]`
   :math:`[`\ double03\ :math:`]`  :math:`[`\ double04\ :math:`]`

-  \ ``NOmega``\ +1行目：空行

-  以下繰り返し

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :** 反復ステップ数を表します。10の倍数のみ出力されます。

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型

   **説明 :** 振動数の値を表します。

-  :math:`[`\ double02\ :math:`]`, :math:`[`\ double03\ :math:`]`

   **形式 :** double型

   | **説明 :** 反復数が\ :math:`[`\ int01\ :math:`]`\ の時の動的グリーン関数の値を出力します。
     :math:`[`\ double02\ :math:`]`\ が実部、\ :math:`[`\ double03\ :math:`]`\ が虚部を表します。

-  :math:`[`\ double04\ :math:`]`

   **形式 :** double型

   | **説明 :** 反復数が\ :math:`[`\ int01\ :math:`]`\ の時の残差ベクトルのノルムを出力します。

.. raw:: latex

   \newpage
