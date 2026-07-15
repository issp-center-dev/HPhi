.. highlight:: none

.. _Subsec:dynamicalG:

DynamicalGreen.dat
~~~~~~~~~~~~~~~~~~

動的グリーン関数の計算結果を出力します。ファイル名およびファイル形式は以下の通りです。

ファイル名
^^^^^^^^^^

-  ##\_DynamicalGreen.dat

##はModParaファイル内の[string02]で指定されるヘッダを表します。

有限温度の固有状態ループおよび multiple-operator / multiple-bra モード
(\ :ref:`ModPara <Subsec:modpara>`\ ファイルの\ ``SpectrumLoopExct``\ 、
``SpectrumNumOp``\ 、\ ``SpectrumNumBra``\ )を有効にすると、
``(固有状態, ket 演算子, bra)``\ の組み合わせごとに1つのファイルが出力され、
0始まりの添字がファイル名に付与されます。

-  ##\_DynamicalGreen\_\ *idx*\ .dat — 固有状態 *idx*\ (\ ``SpectrumLoopExct``\ )。

-  ##\_DynamicalGreen\_\ *idx*\ \_\ *op*\ .dat — 固有状態 *idx*\ 、ket 演算子セット
   *op*\ (\ ``SpectrumNumOp``\ )。

-  ##\_DynamicalGreen\_\ *idx*\ \_\ *op*\ \_\ *bra*\ .dat — 固有状態 *idx*\ 、
   ket 演算子セット *op*\ 、bra 演算子セット *bra*\ (\ ``SpectrumNumBra``\ )。

以下のファイル形式はこれらすべての変種で共通です。

ファイル形式
^^^^^^^^^^^^

-  1行目-:
   :math:`[`\ double01\ :math:`]`  :math:`[`\ double02\ :math:`]`  :math:`[`\ double03\ :math:`]`  :math:`[`\ double04\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ double01\ :math:`]`, :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   | **説明 :** 振動数の値を表します。
   | :math:`[`\ double01\ :math:`]`\ が実部、\ :math:`[`\ double02\ :math:`]`\ が虚部を表します。

-  :math:`[`\ double03\ :math:`]`, :math:`[`\ double04\ :math:`]`

   **形式 :** double型

   | **説明 :** 動的グリーン関数の値を表します。
   | :math:`[`\ double03\ :math:`]`\ が実部、\ :math:`[`\ double04\ :math:`]`\ が虚部を表します。

.. raw:: latex

   \newpage
