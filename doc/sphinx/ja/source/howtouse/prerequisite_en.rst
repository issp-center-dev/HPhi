.. highlight:: none

.. _Ch:Prerequisite:

要件
====

:math:`{\mathcal H}\Phi` のコンパイルおよび使用には次のものが必要です。

 * C/fortran コンパイラ (インテル、富士通、GNUなど)
 * BLAS/LAPACKライブラリ (インテルMKL, 富士通, ATLASなど)
 * MPIライブラリ (MPI並列を行わない場合は必要ありません)
 * ScaLAPACKライブラリ (全対角化で使用しない場合は必要ありません)
 * MAGMAライブラリ (全対角化で使用しない場合は必要ありません)

.. tip::

 | **例/ intelコンパイラーでの設定**
 | intelコンパイラを使用する場合には、コンパイラに付属の設定用スクリプトを使用するのが簡単です。
 | 64ビットOSでbashを使っている場合には
 
 | ``source /opt/intel/bin/compilervars.sh intel64``
 | または
 | ``source /opt/intel/bin/iccvars.sh intel64``
 | ``source /opt/intel/mkl/bin/mklvars.sh``
 
等を\ ``~/.bashrc``\ に記載してください。
詳しくはお手持ちのコンパイラ、ライブラリのマニュアルをお読みください。
