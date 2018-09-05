.. highlight:: none

CHECK\_Hund.dat
~~~~~~~~~~~~~~~

HamiltonianのHundカップリング

.. math::

   \begin{aligned}
   H += -\sum_{i,j}J_{ij}^{\rm Hund} (n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow})\end{aligned}


に関する入力確認を行います。\ :math:`J_{ij}^{\rm Hund}`\ が出力されます。
以下にファイル例を記載します。

::

    i=0 isite1=1 isite2=2 tmp_V=0.250000 
    i=1 isite1=2 isite2=3 tmp_V=0.250000 
    i=2 isite1=3 isite2=4 tmp_V=0.250000 
    i=3 isite1=4 isite2=5 tmp_V=0.250000 
    i=4 isite1=5 isite2=6 tmp_V=0.250000 
    i=5 isite1=6 isite2=1 tmp_V=0.250000 

ファイル形式
^^^^^^^^^^^^

以下のようなファイル形式をとります。

-  i=\ :math:`[`\ int01\ :math:`]` isite1=\ :math:`[`\ int02\ :math:`]`
   isite2=\ :math:`[`\ int03\ :math:`]`
   tmp\_V=\ :math:`[`\ double01\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :** 入力番号。

-  :math:`[`\ int02\ :math:`]`, :math:`[`\ int03\ :math:`]`

   **形式 :** int型

   | **説明 :** :math:`J_{ij}^{\rm Hund}`\ のサイト番号を表す整数。
   | :math:`[`\ int02\ :math:`]`\ が\ :math:`i`,
     :math:`[`\ int03\ :math:`]`\ が\ :math:`j`\ を表します。

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型

   **説明 :** :math:`J_{ij}^{\rm Hund}`\ の値。

.. raw:: latex

   \newpage