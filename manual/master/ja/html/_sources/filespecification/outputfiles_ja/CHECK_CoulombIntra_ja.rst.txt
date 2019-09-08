.. highlight:: none

CHECK\_CoulombIntra.dat
~~~~~~~~~~~~~~~~~~~~~~~

Hamiltonianのオンサイトクーロン相互作用

.. math:: H+=\sum_{i} U_i n_{i\uparrow} n_{j \downarrow}

に関する入力確認を行います。\ :math:`U_i`\ が出力されます。
以下にファイル例を記載します。

::

    i=0 isite1=1 tmp_V=4.000000 
    i=1 isite1=2 tmp_V=4.000000 
    i=2 isite1=3 tmp_V=4.000000 
    i=3 isite1=4 tmp_V=4.000000 
    i=4 isite1=5 tmp_V=4.000000 
    i=5 isite1=6 tmp_V=4.000000 

ファイル形式
^^^^^^^^^^^^

以下のようなファイル形式をとります。

-  i=\ :math:`[`\ int01\ :math:`]` isite1=\ :math:`[`\ int02\ :math:`]`
   tmp\_V=\ :math:`[`\ double01\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :** 入力番号。

-  :math:`[`\ int02\ :math:`]`

   **形式 :** int型

   **説明 :** :math:`U_i`\ のサイト番号\ :math:`i`\ を表す整数。

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型

   **説明 :** :math:`U_i`\ の値。

.. raw:: latex

   \newpage