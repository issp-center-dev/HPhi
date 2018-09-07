.. highlight:: none

.. _Subsec:flctrand:

Flct\_rand.dat
~~~~~~~~~~~~~~

| (TPQ法でのみ出力)
  TPQ法での揺らぎに関する有限温度計算結果を出力します。
  再計算の場合は値が追記されます。 以下にファイル例を示します。

::

     # inv_temp, N, N^2, D, D^2, Sz, Sz^2, step_i
    0.0826564 12.00 144.00 0.00 0.00 0.0009345626081113 0.2500 1
    0.1639935 12.00 144.00 0.00 0.00 0.0023147006319775 0.2500 2
    0.2440168 12.00 144.00 0.00 0.00 0.0037424057659867 0.2500 3
    …
    135.97669 12.00 144.00 0.00 0.00 -0.0000000000167368 0.2500 1998
    136.04474 12.00 144.00 0.00 0.00 -0.0000000000165344 0.2500 1999

ファイル名
^^^^^^^^^^

-  Flct\_rand??.dat

??はTPQ法計算時のrunの番号を表します。

ファイル形式
^^^^^^^^^^^^

1行目はヘッダで、2行目以降は以下のファイル形式で記載されます。

-  :math:`[`\ double01\ :math:`]` :math:`[`\ double02\ :math:`]`
   :math:`[`\ double03\ :math:`]` :math:`[`\ double04\ :math:`]`
   :math:`[`\ double05\ :math:`]` :math:`[`\ double06\ :math:`]`
   :math:`[`\ double07\ :math:`]` :math:`[`\ int01\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型

   **説明 :** 逆温度\ :math:`1/{k_{\rm B}T}`\ 。

-  :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   **説明 :** 全粒子数\ :math:`\sum_{i} \langle \hat{n}_i \rangle`\ 。

-  :math:`[`\ double03\ :math:`]`

   **形式 :** double型

   **説明 :**
   粒子数の2乗の期待値\ :math:` \langle (\sum_{i}\hat{n}_i)^2 \rangle`\ 。

-  :math:`[`\ double04\ :math:`]`

   **形式 :** double型

   **説明 :** ダブロン
   :math:`\frac{1}{N_s} \sum_{i}\langle n_{i\uparrow}n_{i\downarrow}\rangle`
   の期待値(ただし\ :math:`N_s`\ はサイト数)。

-  :math:`[`\ double05\ :math:`]`

   **形式 :** double型

   **説明 :** ダブロンの二乗
   :math:`\frac{1}{N_s} \langle (\sum_{i}n_{i\uparrow} n_{i\downarrow})^2\rangle`
   の期待値(ただし\ :math:`N_s`\ はサイト数)。

-  :math:`[`\ double06\ :math:`]`

   **形式 :** double型

   **説明 :** スピンの\ :math:`Sz`\ 成分
   :math:`\frac{1}{N_s} \sum_{i}\langle \hat{S}_i^z\rangle`
   の期待値(ただし\ :math:`N_s`\ はサイト数)。

-  :math:`[`\ double07\ :math:`]`

   **形式 :** double型

   **説明 :** スピンの\ :math:`Sz`\ 成分の二乗
   :math:`\frac{1}{N_s} \langle (\sum_{i} \hat{S}_i^z)^2\rangle`
   の期待値(ただし\ :math:`N_s`\ はサイト数)。

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :**
   初期ランダムベクトルに\ :math:`(l-\hat{H}/N_{s})`\ (:math:`l`\ はModParaファイルの\ ``LargeValue``\ 、\ :math:`N_{s}`\ はサイト数)を作用させた回数。


.. raw:: latex

   \newpage