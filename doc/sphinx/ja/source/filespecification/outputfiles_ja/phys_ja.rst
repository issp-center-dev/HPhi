.. highlight:: none

.. _Subsec:phys:

phys.dat
~~~~~~~~

| (FullDiagでのみ出力)全対角法で計算したエネルギーと物理量を出力します。エネルギーの低い基底エネルギーから順に出力されます。以下にファイル例を記載します。

::

     <H>         <N>        <Sz>       <S2>       <D> 
      -4.814170   0.000000   0.000000  -0.000000   0.590568
      -3.796850   0.000000   0.000000   1.333333   0.423804
     …
     14.489622   0.000000   0.000000   0.000000   2.550240
     14.852520   0.000000   0.000000   0.000000   2.329157

ファイル名
^^^^^^^^^^

-  カノニカル: ##\_phys\_Nup\_$$Ndown%%.dat

-  グランドカノニカル: ##\_phys.dat

##はModParaファイル内の[string02]で指定されるヘッダ、$$はNup、%%はNdownを表します。

ファイル形式
^^^^^^^^^^^^

1行目はヘッダで、2行目以降は以下のファイル形式で記載されます。

-  :math:`[`\ double01\ :math:`]` :math:`[`\ double02\ :math:`]`
   :math:`[`\ double03\ :math:`]` :math:`[`\ double04\ :math:`]`
   :math:`[`\ double05\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型

   **説明 :** エネルギーの期待値\ :math:`\langle H\rangle`\ 。

-  :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   **説明 :** 粒子数の期待値\ :math:`\langle \hat{n}\rangle`\ 。

-  :math:`[`\ double03\ :math:`]`

   **形式 :** double型

   **説明 :** スピンのz成分の期待値\ :math:`\langle S_z\rangle`\ 。

-  :math:`[`\ double04\ :math:`]`

   **形式 :** double型

   **説明 :** スピンの2乗の期待値\ :math:`\langle  S^2\rangle`\ 。

-  :math:`[`\ double05\ :math:`]`

   **形式 :** double型

   **説明 :** ダブロン
   :math:`\frac{1}{N_s} \sum_{i}\langle n_{i\uparrow}n_{i\downarrow}\rangle`
   (ただし:math:`N_s`\ はサイト数)。

.. raw:: latex

   \newpage