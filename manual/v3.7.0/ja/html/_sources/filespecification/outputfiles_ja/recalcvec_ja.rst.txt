.. highlight:: none

recalcvec.dat
~~~~~~~~~~~~~

Lanczso法を用いた動的グリーン関数の再計算に必要な2つのベクトルが出力されます。
ファイルはバイナリ形式で出力されます。
ファイル名およびファイル形式は以下の通りです。

ファイル名
^^^^^^^^^^

-  ##\_recalcvec\_rank\_$$.dat

##はModParaファイル内の[string02]で指定されるヘッダ、$$はランク番号を表します。

ファイル形式
^^^^^^^^^^^^

-  1行目：\ :math:`[`\ int01\ :math:`]`

-  2行目：\ :math:`[`\ int02\ :math:`]`

-  3行目 - :math:`[`\ int02\ :math:`]`\ +3行:
   :math:`[`\ double01\ :math:`]`  :math:`[`\ double02\ :math:`]`

-  4行目 - 2\ :math:`\times`\ :math:`[`\ int02\ :math:`]`\ +4行:
   :math:`[`\ double03\ :math:`]`  :math:`[`\ double04\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :** 動的グリーン関数の計算に要したステップ数を表します。

-  :math:`[`\ int02\ :math:`]`

   **形式 :** Long int型

   **説明 :** 計算対象のヒルベルト空間数を指定する整数。

-  :math:`[`\ double01\ :math:`]`, :math:`[`\ double02\ :math:`]`

   **形式 :** double

   | **説明 :** Lanczos法での\ :math:`\vec{v}_{k+1}`\ を出力します。
     :math:`[`\ double01\ :math:`]`\ が実部、\ :math:`[`\ double02\ :math:`]`\ が虚部を表します。一番最初の成分に計算に使用しない値が入っています。

-  :math:`[`\ double03\ :math:`]`, :math:`[`\ double04\ :math:`]`

   **形式 :** double

   | **説明 :** Lanczos法での\ :math:`\vec{v}_{k+1}`\ を出力します。
     :math:`[`\ double01\ :math:`]`\ が実部、\ :math:`[`\ double02\ :math:`]`\ が虚部を表します。一番最初の成分に計算に使用しない値が入っています。

.. raw:: latex

   \newpage
