.. highlight:: none

.. _Subsec:ssrand:


SS\_rand.dat
~~~~~~~~~~~~

| (TPQ法でのみ出力)TPQ法での有限温度計算結果を出力します。
  再計算の場合は値が追記されます。 以下にファイル例を示します。

::

    # inv_tmp, energy, phys_var, phys_doublon, phys_num, step_i
    0.017471  5.526334 45.390269 1.464589 6.000000 1
    0.034863  5.266718 42.655559 1.434679 6.000000 2
    …
    31.999572  -4.814170 23.176231 0.590568 6.000000 1997
    32.015596  -4.814170 23.176231 0.590568 6.000000 1998
    32.031620  -4.814170 23.176231 0.590568 6.000000 1999

ファイル名
^^^^^^^^^^

-  SS\_rand??.dat

??はTPQ法計算時のrunの番号を表します。

ファイル形式
^^^^^^^^^^^^

1行目はヘッダで、2行目以降は以下のファイル形式で記載されます。

-  :math:`[`\ double01\ :math:`]` :math:`[`\ double02\ :math:`]`
   :math:`[`\ double03\ :math:`]` :math:`[`\ double04\ :math:`]`
   :math:`[`\ double05\ :math:`]` :math:`[`\ int01\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型

   **説明 :** 逆温度\ :math:`1/{k_{\rm B}T} ~(`\ ただし,
   :math:`k_{\rm B} = 1)`\ 。

-  :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   **説明 :** エネルギーの期待値\ :math:`\langle H \rangle`\ 。

-  :math:`[`\ double03\ :math:`]`

   **形式 :** double型

   **説明 :**
   ハミルトニアンの2乗の期待値\ :math:`\langle H^2 \rangle`\ 。

-  :math:`[`\ double04\ :math:`]`

   **形式 :** double型

   **説明 :** ダブロン
   ダブロンの期待値\ :math:`\sum_{i}\langle n_{i\uparrow}n_{i\downarrow}\rangle`
   。

-  :math:`[`\ double05\ :math:`]`

   **形式 :** double型

   **説明 :** 粒子数\ :math:`\langle {\hat n} \rangle`\ 。

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :**
   初期ランダムベクトルに\ :math:`(l-\hat{H}/N_{s})`\ (:math:`l`\ はModParaファイルの\ ``LargeValue``\ 、\ :math:`N_{s}`\ はサイト数)を作用させた回数。


.. raw:: latex

   \newpage