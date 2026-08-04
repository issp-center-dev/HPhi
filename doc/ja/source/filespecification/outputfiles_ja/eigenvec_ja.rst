.. highlight:: none

.. _Subsec:eigenvec:

eigenvec.dat
~~~~~~~~~~~~

``CalcMod`` ファイルで ``OutputEigenVec=1`` とした場合に、固有ベクトルを
バイナリ形式で出力します。反復法では指定された状態を出力し、FullDiagでは
全固有ベクトルを出力します。

``InputEigenVec=1`` は、対応するリスタート・スペクトル計算で反復法の
出力ファイルを読み込みます。FullDiag出力も同じバイナリ配置ですが、
FullDiagのリスタート機能ではありません。FullDiagでは常に対角化を実行します。

ファイル名
^^^^^^^^^^

-  ##\_eigenvec\_&&\_rank\_$$.dat

##はModParaファイル内の[string02]で指定されるヘッダ、&&は0始まりの
固有状態番号、$$はランク番号を表します。

反復法では、各MPIランクが局所Hilbert空間の成分を自身のランク番号の
ファイルへ出力します。FullDiagでは各ファイルに分割されていない1本の
固有ベクトルを格納するため、ScaLAPACK/ELPAや ``ExpecMode=1`` / ``2`` を
用いたMPI実行でも ``##_eigenvec_&&_rank_0.dat`` だけを出力します。

ファイル形式
^^^^^^^^^^^^

ファイルには、順に反復回数（``int``）、ベクトル次元
（``unsigned long int``）、``complex double`` 型の
``local_size+1`` 成分を格納します。

**Note:**  ``eigen_vector``\ の一番最初の成分は計算に使用しません。

FullDiagでは反復回数は0、ベクトル次元はHilbert空間の全次元です。この指定を
有効にすると、:math:`N_{\rm H}` 本のファイルへ合計
:math:`N_{\rm H}^2` 個の複素数成分を出力するため、ディスク使用量は
:math:`O(N_{\rm H}^2)` になります。

.. raw:: latex

   \newpage
