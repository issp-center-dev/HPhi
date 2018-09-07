.. highlight:: none

.. _Subsec:eigenvec:

eigenvec.dat
~~~~~~~~~~~~

CalcModファイルのOutputEigenVec=1の場合に、Lanczos法で計算された固有ベクトルを出力します。IntputEigenVec=1の場合には、出力されたファイルの読み込みを行います。ファイルはバイナリ形式で出力されます。ファイル名およびファイル形式は以下の通りです。

ファイル名
^^^^^^^^^^

-  ##\_eigenvec\_&&\_rank\_$$.dat

##はModParaファイル内の[string02]で指定されるヘッダ、&&は固有値の番号、$$はランク番号を表します。

ファイル形式
^^^^^^^^^^^^

**Note:**  ``eigen_vector``\ の一番最初の成分に計算に使用しない値が入っています。

.. raw:: latex

   \newpage