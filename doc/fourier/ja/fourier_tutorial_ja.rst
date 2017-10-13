.. _tutorial:

チュートリアル
==============

このチュートリアルは ``sample/Standard/Spin/HeisenbergSquare/`` (mVMC)
および ``sample/CG/Heisenberg/`` (HPhi)
にあるインプットファイルを用いて行う.

HPhi/vmc.out の実行
-------------------

- :math:`{\mathcal H}\Phi` の場合

  基底状態および相関関数の計算を行う.
  
  .. code-block:: bash

     $ ../../../../src/HPhi -s stan.in

- mVMC の場合

  変分波動関数の最適化を行う.
  
  .. code-block:: bash

     $ ../../../../src/vmc.out -s StdFace.def

  相関関数を計算するために, ``StdFace.def`` に以下の行を付け加える.

  ::

     NVMCCalMode = 1

  相関関数を計算する.
  
  .. code-block:: bash

     $ ../../../../src/vmc.out -s StdFace.def output/zqp_opt.dat
         
これにより, カレントディレクトリの ``output/`` 以下に
1体および2体の相関関数が出力される.

関連するファイル

- StdFace.def (mVMC/:math:`{\mathcal H}\Phi` のマニュアル参照)
- zqp_opt.dat (mVMCのマニュアル参照)
- greenone.def (:ref:`greenindex`)
- greentwo.def (:ref:`greenindex`)

相関関数のフーリエ変換
----------------------

ユーティリティプログラム ``fourier`` を使って,
相関関数をフーリエ変関する.

.. code-block:: bash

   $ ../../../../tool/fourier namelist.def geometry.dat
     
これにより, カレントディレクトリの ``output/`` 以下に
フーリエ変換された相関関数が出力される.

関連するファイル

- output/zvo_cisajs_001.dat (:ref:`zvocisajs`)
- output/zvo_cisajs.dat (:ref:`zvocisajs`)
- output/zvo_cisajscktalt_001.dat (:ref:`zvocisajs`)
- output/zvo_cisajscktalt.dat (:ref:`zvocisajs`)
- geometry.dat (:ref:`geometry`)
- output/zvo_corr.dat (:ref:`zvocorr`)

相関関数のプロット
------------------

ユーティリティプログラム ``corplot`` を使って,
相関関数を :math:`k` 空間でプロットする.

.. code-block:: bash

   $ ../../../../tool/corplot output/zvo_corr.dat
   or
   $ ../../../../tool/corplot output/zvo_corr_eigen0.dat

この時, ターミナルには次のように標準入力を促すメッセージが現れる.

::

    #####  Plot Start  #####

       Please specify target number from below (0 or Ctrl-C to exit):

       Real Part Without ErrorBar
         [ 1] Up-Up [ 2] Down-Down [ 3] Density-Density [ 4] SzSz [ 5] S+S- [ 6] S-S+
       Imaginary Part Without ErrorBar
         [11] Up-Up [12] Down-Down [13] Density-Density [14] SzSz [15] S+S- [16] S-S+
       Real Part With ErrorBar
         [21] Up-Up [22] Down-Down [23] Density-Density [24] SzSz [25] S+S- [26] S-S+
       Imaginary Part With ErrorBar
         [31] Up-Up [32] Down-Down [33] Density-Density [34] SzSz [35] S+S- [36] S-S+

      Target : 

プロットしたい量に対応する数字(例えば4)を入力し,
``Enter`` キーを押すと gnuplot が起動して3Dグラフが表示される(図 :num:`corplotpng` ).

.. _corplotpng:
     
.. figure:: ../figs/corplot.png

            Target : 4 としてプロットした図.
            黒線は第一ブリルアンゾーンを表す.

関連するファイル

- kpoint.dat (:ref:`kpoint`)
- correlation.gp (:ref:`gnuplot`)
- correlation.dat (:ref:`correlation`)
