インストール方法
================

:math:`{\mathcal H}\Phi` は次の場所からダウンロードできます。

https://github.com/issp-center-dev/HPhi/releases

ダウンロードしたファイルを次のように展開してください。

``$ tar xzvf HPhi-xxx.tar.gz``

:math:`{\mathcal H}\Phi`\ はcmakeを利用してインストールできます。
:math:`{\mathcal H}\Phi`\ を展開したディレクトリのパスを$PathTohphi
、ビルドディレクトリを$HOME/build/hphi
(任意の場所を指定可能)とした場合に、 ::

 cd $HOME/build/hphi
 cmake -DCONFIG=gcc $PathTohphi
 make

でコンパイルすることができます。全対角化の計算に\ ``ScaLAPACK``\ を使用する場合には、\ ``-DUSE_SCALAPACK=ON``\
のオプションをcmake時につけてください。
全対角化の計算に\ ``ELPA``\ を使用する場合には、\ ``-DUSE_ELPA=ON``\
のオプションをcmake時につけてください（ScaLAPACKとMPIが必要です。
``USE_SCALAPACK`` は自動的に有効化されます）。ELPAが標準外のパスに
インストールされている場合は\ ``-DELPA_ROOT=/path/to/elpa``\ で場所を
指定してください（``ELPA_INCLUDE_DIR``/``ELPA_LIBRARY`` での個別指定も
可能です）。

.. note::

   ELPAビルド時のプラットフォーム注意（AMD EPYC + Mellanox InfiniBand
   クラスタ、Intel oneAPI 2022 / Intel MPI 2021.7 で実機確認済み。
   いずれも物性研スパコン kugui で実際に遭遇・診断した事例です）:

   * **AVX-512 を持たない CPU**\ （AMD EPYC Rome/Milan など）では、ELPA
     本体の configure に ``--disable-avx512 --disable-avx512-kernels``
     を付けてください。ELPA はビルドホストの能力に関係なく AVX-512
     カーネルをコンパイルし、ランタイムのカーネル選択がそれを選ぶ
     ことがあります。その場合、2段ソルバ経路
     （\ ``single_hh_trafo_*_AVX512_*``\ ）で不正命令（SIGILL）により
     クラッシュします。小さい行列（1段ソルバ経路）では発症しない
     ため、見落としやすい点に注意してください。
   * **Intel MPI** を Mellanox（mlx/UCX プロバイダ）環境で使う場合、
     ``ExpecMode`` 1 または 2 の複数ランク実行が物理量収集
     （\ ``MPI_Gatherv``\ ）で既定チューニングのアルゴリズムにより
     デッドロックすることがあります。回避策として
     ``I_MPI_ADJUST_GATHERV=3``\ （マルチノード可）または
     ``I_MPI_FABRICS=shm``\ （単一ノードのみ）を設定してください。

コンパイル後、$HOME/build/hphi直下にsrcフォルダが作成され、
実行ファイルであるHPhiがそのフォルダ内に作成されます。
MPIライブラリがない場合には、MPI非対応の実行ファイルが作成されます。

なお、上の例ではgccコンパイラを前提としたコンパイルになっていますが、

 *  ``sekirei`` : 物性研究所システムB ”sekirei”
 *  ``sekirei_acc`` : 物性研究所システムB ”sekirei” (MAGMAライブラリを使用する場合)
 * ``fujitsu`` : 富士通コンパイラ
 * ``intel`` : intelコンパイラ + Linux PC
 * ``gcc`` : GCC + Linux PC

のオプションが用意されています。
以下、HPhiを展開したディレクトリでビルドする例を示します(intelコンパイラの場合)。

::

 mkdir ./build
 cd ./build
 cmake -DCONFIG=intel ../
 make

実行後、buildフォルダ直下にsrcフォルダが作成され、HPhiがsrcフォルダ内に作成されます。
なお、コンパイラを変更しコンパイルし直したい場合には、都度buildフォルダごと削除を行った上で、
新規に上記作業を行うことをお薦めします。
また、SSE2が使用出来る場合には、cmakeでのコンパイル時\ ``-DHAVE_SSE2``\ を付け加えてください.


.. tip::

 | CMake 中に StdFace (https://github.com/issp-center-dev/StdFace, スタンダードモードのパーサー) が自動でダウンロードされます。
 | もしもこのダウンロードに失敗した場合は、 ``sh src/StdFace/download.sh`` を実行してダウンロードをしてください。
 | 
 | もしくは、別の場所にダウンロードした StdFace を使いたい場合は、
 | ``-DSTDFACE_DIR=<path_to_stdface>`` を cmake コマンドに渡してください。
