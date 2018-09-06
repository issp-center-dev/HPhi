.. highlight:: none

.. _Subsec:normrand:

Norm\_rand.dat
~~~~~~~~~~~~~~

| (TPQ法でのみ出力) TPQ法での有限温度計算時のログを出力します。
  再計算の場合は値が追記されます。 以下にファイル例を示します。

::

     # inv_temp, global_norm, global_1st_norm, step_i 
    0.017471 19.046586 11.288975 1
    0.034863 19.089752 11.288975 2
    …
    31.999572 20.802362 11.288975 1997
    32.015596 20.802362 11.288975 1998
    32.031620 20.802362 11.288975 1999

ファイル名
^^^^^^^^^^

-  Norm\_rand??.dat

??はTPQ法計算時のrunの番号を表します。

ファイル形式
^^^^^^^^^^^^

1行目はヘッダで、2行目以降は以下のファイル形式で記載されます。

-  :math:`[`\ double01\ :math:`]` :math:`[`\ double02\ :math:`]`
   :math:`[`\ double03\ :math:`]` :math:`[`\ int01\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型

   **説明 :** 逆温度\ :math:`1/{k_{\rm B}T}`\ 。

-  :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   **説明 :** TPQ法で計算される規格化前の波動関数(ベクトル)のノルム:
   :math:`\langle \tilde{\psi}_{k} |\tilde{\psi}_{k}\rangle`\ 。ただし、\ :math:`|\tilde{\psi}_{k}\rangle \equiv(l-\hat{H}/N_{s})|\psi_{k-1}\rangle`\ 。

-  :math:`[`\ double03\ :math:`]`

   **形式 :** double型

   **説明 :**
   規格化前の初期波動関数(ランダムベクトル)のノルム：\ :math:`\langle \tilde{\psi}_{0} |\tilde{\psi}_{0}\rangle`\ 。ただし、\ :math:`|\tilde{\psi}_{0}\rangle`\ は規格化前の初期波動関数。

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :**
   初期ランダムベクトルに\ :math:`(l-\hat{H}/N_{s})`\ (:math:`l`\ はModParaファイルの\ ``LargeValue``\ 、\ :math:`N_{s}`\ はサイト数)を作用させた回数。

.. raw:: latex

   \newpage