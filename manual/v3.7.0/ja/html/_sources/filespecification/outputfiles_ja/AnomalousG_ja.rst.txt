.. highlight:: none

.. _Subsec:anomalousgoutput:

AnomalousG.dat
~~~~~~~~~~~~~~

``AnomalousG`` キーワードで指定された異常ペアグリーン関数の計算結果を出力します。
typeが ``0`` の行は :math:`c_{i\sigma_1}c_{j\sigma_2}` の期待値、typeが
``1`` の行は :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\dagger}` の期待値を
表します。以下にファイル例を記載します。

::

    0 0 0 0 1 -0.5000000000 0.0000000000
    1 0 1 0 0 -0.5000000000 0.0000000000

ファイル名
^^^^^^^^^^

Lanczos法、LOBCG法: ##\_AnomalousG.dat

TPQ法、cTPQ法: ##\_AnomalousG\_set??step%%.dat

``OutputGreenFormat=1`` のTPQ/cTPQ法: ##\_AnomalousG\_tpq.dat

全対角化法: ##\_AnomalousG\_eigen&&.dat

``OutputGreenFormat=1`` の全対角化法、LOBCG法: ##\_AnomalousG\_eigen.dat

実時間発展法: ##\_AnomalousG\_step%%.dat

``OutputGreenFormat=1`` の実時間発展法: ##\_AnomalousG\_te.dat

##はModParaファイル内の[string02]で指定されるヘッダ、??はTPQ法計算時のrunの番号、%%はTPQ法または時間発展計算でのステップ数、&&は固有値の番号を表します。集約形式では、各行の先頭にTPQ/cTPQでは ``set step``、実時間発展では ``step``、全対角化/LOBCGでは ``eigen`` が追加されます。

ファイル形式
^^^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`  :math:`[`\ int02\ :math:`]`  :math:`[`\ int03\ :math:`]`  :math:`[`\ int04\ :math:`]`  :math:`[`\ int05\ :math:`]`  :math:`[`\ double01\ :math:`]`  :math:`[`\ double02\ :math:`]`

``OutputGreenFormat=1`` の集約ファイルでは、上記の列の前に計算手法に応じた ``set step``、``step``、または ``eigen`` のindex列が追加されます。

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   | **説明 :** 異常ペアのtypeを表します。
   | 0: :math:`\langle c_{i\sigma_1}c_{j\sigma_2}\rangle`
   | 1: :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\dagger}\rangle`

-  :math:`[`\ int02\ :math:`]`, :math:`[`\ int04\ :math:`]`

   **形式 :** int型

   **説明 :**
   サイト番号を指定する整数。\ :math:`[`\ int02\ :math:`]`\ が\ :math:`i`\ サイト、\ :math:`[`\ int04\ :math:`]`\ が\ :math:`j`\ サイトを表します。

-  :math:`[`\ int03\ :math:`]`, :math:`[`\ int05\ :math:`]`

   **形式 :** int型

   | **説明 :**
     スピンを指定する整数。\ :math:`[`\ int03\ :math:`]`\ が\ :math:`\sigma_1`\ 、\ :math:`[`\ int05\ :math:`]`\ が\ :math:`\sigma_2`\ に対応します。
   | 0: アップスピン
   | 1: ダウンスピン。

-  :math:`[`\ double01\ :math:`]`, :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   | **説明 :** 異常ペアグリーン関数の値を表します。
   | :math:`[`\ double01\ :math:`]`\ が実部、\ :math:`[`\ double02\ :math:`]`\ が虚部を表します。

.. raw:: latex

   \newpage
