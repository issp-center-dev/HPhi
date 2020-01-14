.. highlight:: none

CHECK\_INTER\_U.dat
~~~~~~~~~~~~~~~~~~~

Hamiltonianのオフサイトクーロン相互作用

.. math:: H+=\sum_{i} V_{ij} n_{i} n_{j}

に関する入力確認を行います。\ :math:`V_{ij}`\ が出力されます。
以下にファイル例を記載します。

::

    i=0 isite1=1 isite2=2 tmp_V=-0.125000 
    i=1 isite1=2 isite2=3 tmp_V=-0.125000 
    i=2 isite1=3 isite2=4 tmp_V=-0.125000 
    i=3 isite1=4 isite2=5 tmp_V=-0.125000 
    i=4 isite1=5 isite2=6 tmp_V=-0.125000 
    i=5 isite1=6 isite2=1 tmp_V=-0.125000 

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

   | **説明 :** :math:`V_{ij}`\ のサイト番号を表す整数。
   | :math:`[`\ int02\ :math:`]`\ が\ :math:`i`,
     :math:`[`\ int03\ :math:`]`\ が\ :math:`j`\ を表します。

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型

   **説明 :** :math:`V_{ij}`\ の値。


.. raw:: latex

   \newpage