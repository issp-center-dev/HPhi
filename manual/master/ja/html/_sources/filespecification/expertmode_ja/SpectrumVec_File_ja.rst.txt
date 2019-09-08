.. highlight:: none

.. _Subsec:spectrumvec:

SpectrumVec指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~

スペクトル計算用の入力ファイルのヘッダを指定します。
ファイル名およびファイル形式は以下の通りです(ファイル形式はeigenvec.datと同様です)。
ファイルデータはバイナリ形式です。

ファイル名
^^^^^^^^^^

-  ##\_rank\_$$.dat

##はSpectrumVecで指定されるヘッダ、$$はランク番号を表します。また、&&はTPQ時のサンプリングの番号を表します。

ファイル形式
^^^^^^^^^^^^

-  1行目：\ :math:`[`\ int01\ :math:`]`

-  2行目：\ :math:`[`\ int02\ :math:`]`

-  2行目-:
   :math:`[`\ double01\ :math:`]`  :math:`[`\ double02\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :** 計算対象のヒルベルト空間数を指定する整数。

-  :math:`[`\ int02\ :math:`]`

   **形式 :** int型

   **説明 :**
   計算に要したステップ数を表します。Lanczos法ではLanczosステップ数、TPQ法ではハミルトニアンを乗算した回数を記載します。

-  :math:`[`\ double01\ :math:`]`, :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   | **説明 :** 入力するベクトルの値を表します。
   | :math:`[`\ double01\ :math:`]`\ が実部、\ :math:`[`\ double02\ :math:`]`\ が虚部を表します。


.. raw:: latex

   \newpage
