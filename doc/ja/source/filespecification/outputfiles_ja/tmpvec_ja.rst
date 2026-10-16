.. highlight:: none

.. _Subsec:restart:

tmpvecファイル
~~~~~~~~~~~~~~

TPQ/cTPQ法およびLOBPCG（CG）法のlegacy ``ReStart`` で使用する
バイナリファイルです。基底状態Lanczos法はこの形式を使用せず、
``recalcvec`` と ``TMComponents`` を使用します。
新しいcheckpointは :ref:`Subsec:restart_v3` のv3 manifest形式を使い、
``tmpvec`` を出力しません。legacy CGの ``Lanczos_max`` は追加反復数で、表示する
反復labelを1から開始し、以前の ``Lanczos_Step`` prefixを復旧しません。

ファイル名
^^^^^^^^^^

-  ``tmpvec_set&&_rank_$$.dat``

``&&`` は0始まりのTPQサンプル番号またはLOBPCG固有状態番号、
``$$`` はMPIランク番号です。ファイル名に ``CDataFileHead`` は付きません。

ファイル形式
^^^^^^^^^^^^

対応する64 bit環境では、先頭に12 byteのヘッダがあり、その後に
複素ベクトルが1本続きます。

::

    fwrite(&step, sizeof(int), 1, fp);
    fwrite(&local_size, sizeof(unsigned long int), 1, fp);
    fwrite(&vector[0], sizeof(complex double), local_size + 1, fp);

``step`` は保存した反復番号、``local_size`` は各ランクが保持する
ベクトル次元です。互換性のため ``local_size + 1`` 個の複素数を保存し、
``vector[0]`` は計算には使用しません。バイナリ表現は環境に依存します。

.. raw:: latex

   \newpage
