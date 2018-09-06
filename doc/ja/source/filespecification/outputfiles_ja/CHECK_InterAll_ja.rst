.. highlight:: none

CHECK\_InterAll.dat
~~~~~~~~~~~~~~~~~~~

Hamiltonianの一般二体相互作用のうち対角成分

.. math:: H+=\sum_{i,j, \sigma} I_{iijj\sigma_1\sigma_1\sigma_2\sigma_2} c_{i\sigma_1}^{\dagger}c_{i\sigma_1}c_{i\sigma_2}^{\dagger}c_{i\sigma_2}


に関する入力確認を行います。\ :math:`I_{iijj\sigma_1\sigma_1\sigma_2\sigma_2}`\ が出力されます。
以下にファイル例を記載します。

::

    i=0 isite1=1 A_spin=0 isite2=2 B_spin=0 tmp_V=0.500000 
    i=1 isite1=1 A_spin=0 isite2=2 B_spin=1 tmp_V=-0.500000 
    i=2 isite1=1 A_spin=1 isite2=2 B_spin=0 tmp_V=-0.500000 
    i=3 isite1=1 A_spin=1 isite2=2 B_spin=1 tmp_V=0.500000 
    i=4 isite1=2 A_spin=0 isite2=3 B_spin=0 tmp_V=0.500000 
    i=5 isite1=2 A_spin=0 isite2=3 B_spin=1 tmp_V=-0.500000 
    …

ファイル形式
^^^^^^^^^^^^

以下のようなファイル形式をとります。

-  i=\ :math:`[`\ int01\ :math:`]` isite1=\ :math:`[`\ int02\ :math:`]`
   A\_spin=\ :math:`[`\ int03\ :math:`]`
   isite2=\ :math:`[`\ int04\ :math:`]`
   B\_spin=\ :math:`[`\ int05\ :math:`]`
   tmp\_V=\ :math:`[`\ double01\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :** 入力番号。

-  :math:`[`\ int02\ :math:`]`, :math:`[`\ int04\ :math:`]`

   **形式 :** int型

   | **説明 :**
     :math:`I_{iijj\sigma_1\sigma_1\sigma_2\sigma_2}`\ のサイト番号を表す整数。
   | :math:`[`\ int02\ :math:`]`\ が\ :math:`i`,
     :math:`[`\ int04\ :math:`]`\ が\ :math:`j`\ を表します。

-  :math:`[`\ int03\ :math:`]`, :math:`[`\ int05\ :math:`]`

   **形式 :** int型

   | **説明 :**
     :math:`I_{iijj\sigma_1\sigma_1\sigma_2\sigma_2}`\ のスピン番号を表す整数。
   | :math:`[`\ int03\ :math:`]`\ が\ :math:`\sigma_1`,
     :math:`[`\ int05\ :math:`]`\ が\ :math:`\sigma_2`\ に対応し、それぞれ
   | 0: アップスピン
   | 1: ダウンスピン
   | を表します。

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型

   **説明 :** :math:`I_{iijj\sigma_1\sigma_1\sigma_2\sigma_2}`\ の値。

.. raw:: latex

   \newpage