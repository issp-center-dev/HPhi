.. highlight:: none

TMcomponents.dat
~~~~~~~~~~~~~~~~

Lanczos法を用いた動的グリーン関数の再計算に必要な三重対角行列の要素と励起状態のノルムが出力されます。
ファイル名およびファイル形式は以下の通りです。

ファイル名
^^^^^^^^^^

-  ##\_TMcomponents.dat

##はModParaファイル内の[string02]で指定されるヘッダを表します。

ファイル形式
^^^^^^^^^^^^

-  1行目：\ :math:`[`\ int01\ :math:`]`

-  2行目：\ :math:`[`\ double01\ :math:`]`

-  3行目 - :math:`[`\ int02\ :math:`]`\ +3行:
   :math:`[`\ double02\ :math:`]`  :math:`[`\ double03\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :**
   動的グリーン関数の計算に要したステップ数\ :math:`N_d`\ を表します。

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型

   **説明 :** 励起状態のノルムを表します。

-  :math:`[`\ double02\ :math:`]`, :math:`[`\ double03\ :math:`]`

   **形式 :** double型

   | **説明 :**
     Lanczos法を用いた動的グリーン関数の再計算に必要な三重対角行列の要素\ :math:`\alpha_i,~\beta_i~(i =1,\cdots N_d)`\ の値を表します。一番最初の成分に計算に使用しない値が入っています。
     :math:`[`\ double02\ :math:`]`\ が\ :math:`\alpha_i`\ 、\ :math:`[`\ double03\ :math:`]`\ が\ :math:`\beta_i`\ に対応します。


.. raw:: latex

   \newpage

