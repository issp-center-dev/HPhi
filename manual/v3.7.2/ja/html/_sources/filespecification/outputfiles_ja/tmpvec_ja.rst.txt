.. highlight:: none

.. _Subsec:restart:

tmpvec.dat
~~~~~~~~~~

CalcModファイルのReStart=1, 2の場合に、計算途中のベクトルを出力します。
ファイルはバイナリ形式で出力されます。
ファイル名およびファイル形式は以下の通りです(ファイル形式はeigenvec.datと同様です)。

ファイル名
^^^^^^^^^^

-  Lanczos法：##\_tmpvec\_rank\_$$.dat

-  TPQ法, LOBPCG法：##\_tmpvec\_set\_&&\_rank\_$$.dat

##はModParaファイル内の[string02]で指定されるヘッダ、
$$はランク番号を表します。また、&&はTPQ時のサンプリングの番号
もしくはLOBPCG法での固有値の番号を表します。

ファイル形式
^^^^^^^^^^^^

次のようなソースコードを用いて、バイナリファイルとして出力されます
(実際の :math:`{\mathcal H}\Phi`\ のソースコードとは変数名等が異なります)。

| Lanczos

::

    fp = fopen("zvo_tmpvec_rank_0.dat", "wb");
    fwrite(&number_of_interations, sizeof(int), 1,fp);
    fwrite(&local_size, sizeof(unsigned long int),1,fp);
    fwrite(&last_vector[0], sizeof(complex double),local_size+1, fp);
    fwrite(&second_last_vector[0], sizeof(complex double),local_size+1, fp);
    fclose(fp);

| TPQおよびLOBPCG

::

    fp = fopen("zvo_tmpvec_set_0_rank_0.dat", "wb");
    fwrite(&number_of_interations, sizeof(int), 1,fp);
    fwrite(&local_size, sizeof(unsigned long int),1,fp);
    fwrite(&last_vector[0], sizeof(complex double),local_size+1, fp);
    fclose(fp);

ただし、\ ``number_of_interations``\ は反復回数、
``local_size``\ は固有ベクトルのサイズ(MPI並列を使っている場合は全ヒルベルト次元とは異なります)、
``last_vector``\ は最後の反復でのベクトル、
``second_last_vector``\ は最後から2番目の反復でのベクトルをそれぞれ表します。

※
``last_vector``\ 、\ ``second_last_vector``\ の一番最初の成分に計算に使用しない値が入っています。


.. raw:: latex

   \newpage
