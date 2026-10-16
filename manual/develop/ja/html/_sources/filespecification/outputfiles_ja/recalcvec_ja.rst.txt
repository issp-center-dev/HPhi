.. highlight:: none

recalcvec.dat
~~~~~~~~~~~~~

基底状態Lanczos計算またはshifted BiCGによるスペクトル計算を継続するための
ベクトルを、MPIランクごとのバイナリファイルに保存します。

ファイル名
^^^^^^^^^^

-  ``##_recalcvec_rank_$$.dat``

``##`` はModParaファイルの ``CDataFileHead``、``$$`` はMPIランク番号です。

ファイル形式
^^^^^^^^^^^^

基底状態Lanczos法では次の内容を保存します。

::

    unsigned long step;
    unsigned long local_size;
    complex double v0[local_size + 1];
    complex double v1[local_size + 1];

2個のヘッダ値を個別に書き、その後にLanczosベクトル2本を書きます。
各ベクトルの添字0は計算には使用しません。

shifted BiCGでは別のlegacy payloadを保存します。

::

    int step_or_status;
    unsigned long local_size;
    complex double residual[local_size + 1];
    complex double old_residual[local_size + 1];
    complex double shadow_residual[local_size + 1];
    complex double old_shadow_residual[local_size + 1];

バイナリ表現は環境に依存します。リスタート時のMPIプロセス数と各ランクの
ベクトル次元は、保存時と一致させる必要があります。

.. raw:: latex

   \newpage
