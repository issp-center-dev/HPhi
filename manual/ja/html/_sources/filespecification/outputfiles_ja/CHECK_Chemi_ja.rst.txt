.. highlight:: none

.. _Subsec:checkchemi:

CHECK\_Chemi.dat
~~~~~~~~~~~~~~~~

Hamiltonianのうち化学ポテンシャル

.. math:: H+=\sum_{i,\sigma} \mu_{i\sigma} c_{i\sigma}^{\dagger}c_{i\sigma}

に関する入力確認を行います。\ :math:`\mu_{i\sigma}`\ が出力されます。
以下にファイル例を記載します。

::

    i=0 spin=0 isite1=1 tmp_V=0.000000 
    i=1 spin=0 isite1=2 tmp_V=0.000000 
    i=2 spin=0 isite1=3 tmp_V=0.000000 
    i=3 spin=0 isite1=4 tmp_V=0.000000 
    i=4 spin=0 isite1=5 tmp_V=0.000000 
    i=5 spin=0 isite1=6 tmp_V=0.000000 
    …

ファイル形式
^^^^^^^^^^^^

以下のようなファイル形式をとります。

-  i=\ :math:`[`\ int01\ :math:`]` spin=\ :math:`[`\ int02\ :math:`]`
   isite1=\ :math:`[`\ int03\ :math:`]`
   tmp\_V=\ :math:`[`\ double01\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :** 入力番号。

-  :math:`[`\ int02\ :math:`]`

   **形式 :** int型

   | **説明 :**
     :math:`\mu_{i\sigma}`\ のスピン番号\ :math:`\sigma`\ を表す整数。
   | 0: アップスピン
   | 1: ダウンスピン
   | を表します。

-  :math:`[`\ int03\ :math:`]`

   **形式 :** int型

   **説明 :**
   :math:`\mu_{i\sigma}`\ のサイト番号\ :math:`i`\ を表す整数。

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型

   **説明 :** :math:`\mu_{i\sigma}`\ の値。


.. raw:: latex

   \newpage