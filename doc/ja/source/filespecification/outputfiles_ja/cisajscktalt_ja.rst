.. highlight:: none

.. _Subsec:cisajscktalt:

cisajscktalt.dat
~~~~~~~~~~~~~~~~

TwoBodyGで指定された二体グリーン関数\ :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle`\ の計算結果を出力します。以下にファイル例を記載します。

::

        0    0    0    0    0    0    0    0 0.4452776740 0.0000000000
        0    0    0    0    0    1    0    1 0.1843355815 0.0000000000
        0    0    0    0    1    0    1    0 0.1812412105 0.0000000000
        0    0    0    0    1    1    1    1 0.2640364635 0.0000000000
        0    0    0    0    2    0    2    0 0.0279690007 0.0000000000
        0    0    0    0    2    1    2    1 0.2009271524 0.0000000000
        0    0    0    0    3    0    3    0 0.2512810778 0.0000000000
        0    0    0    0    3    1    3    1 0.1939965962 0.0000000000
        …

ファイル名
^^^^^^^^^^

Lanczos法: ##\_cisajscktalt.dat

TPQ法: ##\_cisajscktalt\_set??step%%.dat

実時間発展法: ##\_cisajscktalt\_step%%.dat

全対角化法、LOBCG法: ##\_cisajscktalt\_eigen&&.dat

##はModParaファイル内の[string02]で指定されるヘッダ、??はTPQ法計算時のrunの番号、%%はTPQ法でのステップ数、&&は固有値の番号を表します。

ファイル形式
^^^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`  :math:`[`\ int02\ :math:`]`  :math:`[`\ int03\ :math:`]`  :math:`[`\ int04\ :math:`]`  :math:`[`\ int05\ :math:`]`  :math:`[`\ int06\ :math:`]`  :math:`[`\ int07\ :math:`]`  :math:`[`\ int08\ :math:`]`  :math:`[`\ double01\ :math:`]`  :math:`[`\ double02\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`,
   :math:`[`\ int03\ :math:`]`,\ :math:`[`\ int05\ :math:`]`,
   :math:`[`\ int07\ :math:`]`

   **形式 :** int型

   **説明 :** サイト番号を指定する整数。
   :math:`[`\ int01\ :math:`]`\ が\ :math:`i`\ サイト、\ :math:`[`\ int03\ :math:`]`\ が\ :math:`j`\ サイト、\ :math:`[`\ int05\ :math:`]`\ が\ :math:`k`\ サイト、\ :math:`[`\ int07\ :math:`]`\ が\ :math:`l`\ サイトを表します。

-  :math:`[`\ int02\ :math:`]`,
   :math:`[`\ int04\ :math:`]`,\ :math:`[`\ int06\ :math:`]`,
   :math:`[`\ int08\ :math:`]`

   **形式 :** int型

   | **説明 :** スピンを指定する整数。
     :math:`[`\ int02\ :math:`]`\ が\ :math:`\sigma_1`\ 、\ :math:`[`\ int04\ :math:`]`\ が\ :math:`\sigma_2`\ 、\ :math:`[`\ int06\ :math:`]`\ が\ :math:`\sigma_3`\ 、\ :math:`[`\ int08\ :math:`]`\ が\ :math:`\sigma_4`\ に対応します。
   | 0: アップスピン
   | 1: ダウンスピン
   | を表します。

-  :math:`[`\ double01\ :math:`]`, :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   | **説明 :**
     :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle`\ の値を表します。
   | :math:`[`\ double01\ :math:`]`\ が実部、\ :math:`[`\ double02\ :math:`]`\ が虚部を表します。

.. raw:: latex

   \newpage
