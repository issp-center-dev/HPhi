.. highlight:: none

.. _subsec:energy.dat:

energy.dat
~~~~~~~~~~

| (Lanczos法、LOBCG法のみ)
  Lanczos法もしくはCG法で求めた固有ベクトルを用いて計算したエネルギー、ダブロンと\ :math:`\langle S_z \rangle` を出力します。
  以下にファイル例を記載します。
| ``method="Lanczos"``\ の場合

::

    Energy  -7.1043675920 
    Doublon  0.4164356536 
    Sz  0.0000000000 

``method="CG"``\ の場合

::

    State 0
      Energy  -7.1043675920 
      Doublon  0.4164356536 
      Sz  0.0000000000 

    State 1
    :

ファイル名
^^^^^^^^^^

-  ##\_energy.dat

##はModParaファイル内の[string02]で指定されるヘッダを表します。

ファイル形式
^^^^^^^^^^^^

-  Energy :math:`[`\ double01\ :math:`]`

-  Doublon :math:`[`\ double02\ :math:`]`

-  Sz :math:`[`\ double03\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型

   **説明 :**
   Lanczos法・CG法で求めた固有ベクトルを用い計算したエネルギー。

-  :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   **説明 :** Lanczos法・CG法で求めた固有ベクトルを用い計算したダブロン
   :math:`\frac{1}{N_s} \sum_{i}\langle n_{i\uparrow}n_{i\downarrow}\rangle`
   (:math:`N_s`\ はサイト数)。

-  :math:`[`\ double03\ :math:`]`

   **形式 :** double型

   **説明 :** 固有ベクトルを用い計算した\ :math:`\langle S_z\rangle`\ 。

.. raw:: latex

   \newpage


