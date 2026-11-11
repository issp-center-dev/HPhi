.. highlight:: none

.. _Subsec:nbodygdat:

NBodyG.dat
~~~~~~~~~~

NBodyG指定ファイルで指定したN体グリーン関数の計算結果を出力します。1行で
:math:`N` 個の因子を持つ

.. math::

   \left\langle \prod_{p=1}^{N}
   c_{i_p\sigma'_p}^{\dagger} c_{j_p\sigma_p} \right\rangle

の値を表します。以下にファイル例を記載します。

::

    1    0    0    1    0 0.1250000000 0.0000000000
    2    0    0    0    0    1    1    1    1 0.2500000000 0.0000000000
    3    0    0    0    0    1    1    1    1    2    0    2    0 0.0312500000 0.0000000000

ファイル名
^^^^^^^^^^

Lanczos法: ##\_NBodyG.dat

TPQ法: ##\_NBodyG\_set??step%%.dat

``OutputGreenFormat=1`` のTPQ/cTPQ法: ##\_NBodyG\_tpq.dat

全対角化法、LOBCG法: ##\_NBodyG\_eigen&&.dat

``OutputGreenFormat=1`` の全対角化法、LOBCG法: ##\_NBodyG\_eigen.dat

実時間発展法: ##\_NBodyG\_step%%.dat

``OutputGreenFormat=1`` の実時間発展法: ##\_NBodyG\_te.dat

##はModParaファイル内の[string02]で指定されるヘッダ、??はTPQ法計算時のrunの番号、%%はTPQ法または実時間発展法でのステップ数、&&は固有値の番号を表します。集約形式では、各行の先頭にTPQ/cTPQでは ``set step``、実時間発展では ``step``、全対角化/LOBCGでは ``eigen`` が追加されます。

ファイル形式
^^^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]` (:math:`[`\ int02\ :math:`]` :math:`[`\ int03\ :math:`]` :math:`[`\ int04\ :math:`]` :math:`[`\ int05\ :math:`]`) ... :math:`[`\ double01\ :math:`]` :math:`[`\ double02\ :math:`]`

   括弧で示した4整数の組を :math:`[`\ int01\ :math:`]` 回繰り返します。

``OutputGreenFormat=1`` の集約ファイルでは、上記の列の前に計算手法に応じた ``set step``、``step``、または ``eigen`` のindex列が追加されます。

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :** この成分に含まれる因子数 :math:`N` を表します。

-  各因子の :math:`[`\ int02\ :math:`]`, :math:`[`\ int04\ :math:`]`

   **形式 :** int型

   **説明 :** サイト番号を指定する整数。
   :math:`[`\ int02\ :math:`]`\ が\ :math:`i_p`\ 、
   :math:`[`\ int04\ :math:`]`\ が\ :math:`j_p`\ に対応します。

-  各因子の :math:`[`\ int03\ :math:`]`, :math:`[`\ int05\ :math:`]`

   **形式 :** int型

   | **説明 :** スピンまたは局所状態を指定する整数。
   | Hubbard、tJ、Kondo およびそれらの GC/NConserved 版では
   | 0: アップスピン
   | 1: ダウンスピン
   | を表します。
   | SpinlessFermion/SpinlessFermionGC では 0 のみ出力されます。

-  :math:`[`\ double01\ :math:`]`, :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   | **説明 :** 指定したN体グリーン関数の値を表します。
   | :math:`[`\ double01\ :math:`]`\ が実部、\ :math:`[`\ double02\ :math:`]`\ が虚部を表します。

.. raw:: latex

   \newpage
