.. highlight:: none

.. _Subsec:cgcisajs:

cisajs.dat
~~~~~~~~~~

OneBodyGで指定された一体グリーン関数\ :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`\ の計算結果を出力します。以下にファイル例を記載します。

::

        0    0    0    0 0.4452776740 0.0000000000
        0    1    0    1 0.4452776740 0.0000000000
        1    0    1    0 0.5000000000 0.0000000000
        1    1    1    1 0.5000000000 0.0000000000
        2    0    2    0 0.4452776740 0.0000000000
        2    1    2    1 0.4452776740 0.0000000000
        3    0    3    0 0.5000000000 0.0000000000
        3    1    3    1 0.5000000000 0.0000000000
        …

ファイル名
^^^^^^^^^^

Lanczos法: ##\_cisajs.dat

TPQ法: ##\_cisajs\_set??step%%.dat

実時間発展法: ##\_cisajs\_step%%.dat

全対角化法、LOBCG法: ##\_cisajs\_eigen&&.dat

##はModParaファイル内の[string02]で指定されるヘッダ、??はTPQ法計算時のrunの番号、%%はTPQ法でのステップ数、&&は固有値の番号を表します。

ファイル形式
^^^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`  :math:`[`\ int02\ :math:`]`  :math:`[`\ int03\ :math:`]`  :math:`[`\ int04\ :math:`]`  :math:`[`\ double01\ :math:`]`  :math:`[`\ double02\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`, :math:`[`\ int03\ :math:`]`

   **形式 :** int型

   **説明 :**
   サイト番号を指定する整数。\ :math:`[`\ int01\ :math:`]`\ が\ :math:`i`\ サイト、\ :math:`[`\ int03\ :math:`]`\ が\ :math:`j`\ サイトを表します。

-  :math:`[`\ int02\ :math:`]`, :math:`[`\ int04\ :math:`]`

   **形式 :** int型

   | **説明 :**
     スピンを指定する整数。\ :math:`[`\ int02\ :math:`]`\ が\ :math:`\sigma_1`\ 、\ :math:`[`\ int03\ :math:`]`\ が\ :math:`\sigma_2`\ に対応します。
   | 0: アップスピン
   | 1: ダウンスピン
   | を表します。

-  :math:`[`\ double01\ :math:`]`, :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   | **説明 :**
     :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`\ の値を表します。
   | :math:`[`\ double01\ :math:`]`\ が実部、\ :math:`[`\ double02\ :math:`]`\ が虚部を表します。

.. raw:: latex

   \newpage
