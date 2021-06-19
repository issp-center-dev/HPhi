.. _fileformat:

ファイルフォーマット
====================

.. _geometry:

サイトの位置と軌道のインデックス, *k* 点
----------------------------------------

:ref:`tutorial` でのファイル名は ``geometry.dat`` .
各サイトの位置と軌道の情報は
mVMC/:math:`{\mathcal H}\Phi` のスタンンダードモードを用いた場合には
自動的に生成される.

::

   1.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00 (1)
   0.000000000000000e+00     1.000000000000000e+00     0.000000000000000e+00 (1)
   0.000000000000000e+00     0.000000000000000e+00     1.000000000000000e+00 (1)
   0.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00 (2)
   2 2 0        (3)
   -2 2 0       (3)
   0 0 1        (3)
   0 0 0 0      (4)
   -1 1 0 0     (4)
   0 1 0 0      (4)
   1 1 0 0      (4)
   -1 2 0 0     (4)
   0 2 0 0      (4)
   1 2 0 0      (4)
   0 3 0 0      (4)
   4 20         (5)
   G 0 0 0      (6)
   X 0.5 0 0    (6)
   M 0.5 0.5 0  (6)
   G 0 0 0      (6)
   16 16 1      (7)

#. 単位格子ベクトル. 任意の単位 (スタンダードモードで自動生成).
#. 1体項がシミュレーションセルの境界を跨いだときに付く位相(単位degree)
   (スタンダードモードで自動生成)   
#. シミュレーションセルの形状を指定する三本の整数ベクトル.
   スタンダードモードの入力パラメーター ``a0W``, ``a0L``, ``a0H``, ``a1W``...
   に対応する(スタンダードモードで自動生成).
#. 各サイトの座標結晶並進ベクトル(指数)および内部座標(軌道)のインデックス
   (スタンダードモードで自動生成).
#. *k* パスのノード(対称性の高い点)の数と, ノード間の *k* 点の分割数.
#. *k* ノードのラベルとフラクショナル座標
#. 運動量分布関数のFermiSurferファイルを作成する時の *k* グリッド

サイト表示の1体および2体相関関数
--------------------------------

.. _greenindex:

計算する相関関数のインデックスの指定
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mVMC/:math:`{\mathcal H}\Phi` で計算する相関関数を指定する.
スタンダードモードを使った場合には自動的に生成される.
総合的な説明はmVMC/:math:`{\mathcal H}\Phi` のマニュアルを参照.
:ref:`tutorial` でのファイル名は ``greenone.def`` (1体)および ``greentwo.def`` (2体)である.

:ref:`supported` にある相関関数を計算するためには, 
以下のようにインデックスを指定する必要がある.

- :math:`\langle {\hat c}_{{\bf k}\alpha\uparrow}^{\dagger} {\hat c}_{{\bf k}\beta\uparrow}\rangle`

  :math:`\langle {\hat c}_{{\bf 0}\alpha\uparrow}^{\dagger} {\hat c}_{{\bf R}\beta\uparrow}\rangle`
  に対して, :math:`{\bf R}` が全ての単位胞,
  :math:`\alpha` および :math:`\beta` がそれぞれ全ての軌道を網羅するようにする.
  
- :math:`\langle {\hat c}_{{\bf k}\alpha\downarrow}^{\dagger} {\hat c}_{{\bf k}\beta\downarrow}\rangle`

  :math:`\langle {\hat c}_{{\bf 0}\alpha\downarrow}^{\dagger} {\hat c}_{{\bf R}\beta\downarrow}\rangle`
  に対して, :math:`{\bf R}` が全ての単位胞,
  :math:`\alpha` および :math:`\beta` がそれぞれ全ての軌道を網羅するようにする.
  
- :math:`\langle {\hat \rho}_{{\bf k}\alpha} {\hat \rho}_{{\bf k}\beta}\rangle` および
  :math:`\langle {\hat S}_{{\bf k}\alpha}^{z} {\hat S}_{{\bf k}\beta}^{z} \rangle`

  :math:`\langle {\hat c}_{{\bf 0}\alpha\sigma}^{\dagger} {\hat c}_{{\bf 0}\alpha\sigma} {\hat c}_{{\bf R}\beta \sigma'}^{\dagger} {\hat c}_{{\bf R}\beta \sigma'}\rangle`
  に対して, :math:`{\bf R}` が全ての単位胞,
  :math:`\alpha` および :math:`\beta` がそれぞれ全ての軌道を網羅し,
  :math:`\sigma` および :math:`\sigma'` が :math:`\uparrow`, :math:`\downarrow` を網羅するようにする.

- :math:`\langle {\hat S}_{{\bf k}\alpha}^{+} {\hat S}_{{\bf k}\beta}^{-} \rangle` および
  :math:`\langle {\hat {\bf S}}_{{\bf k}\alpha} \cdot {\hat {\bf S}}_{{\bf k}\beta} \rangle`

  :math:`{\mathcal H}\Phi` の場合は
  :math:`\langle {\hat c}_{{\bf 0}\alpha\sigma}^{\dagger} {\hat c}_{{\bf 0}\alpha-\sigma} {\hat c}_{{\bf R}\beta -\sigma}^{\dagger} {\hat c}_{{\bf R}\beta \sigma}\rangle`
  に対して, :math:`{\bf R}` が全ての単位胞,
  :math:`\alpha` および :math:`\beta` がそれぞれ全ての軌道を網羅し,
  :math:`\sigma` が :math:`\uparrow`, :math:`\downarrow` を網羅するようにする.
  mVMC の場合は
  :math:`\langle {\hat c}_{{\bf 0}\alpha\sigma}^{\dagger} {\hat c}_{{\bf R}\beta \sigma} {\hat c}_{{\bf R}\beta -\sigma}^{\dagger} {\hat c}_{{\bf 0}\alpha-\sigma}\rangle`
  に対して, :math:`{\bf R}` が全ての単位胞,
  :math:`\alpha` および :math:`\beta` がそれぞれ全ての軌道を網羅し,
  :math:`\sigma` が :math:`\uparrow`, :math:`\downarrow` を網羅するようにする.
  いずれの場合も演算子の順番に注意にすること.
  
スタンダードモードのデフォルト(``outputmode="corr"``)では,
自動的に上記のインデックスが指定されるため, 特に気にする必要はない.

.. _zvocisajs:

サイト表示の1体および2体相関関数の計算結果
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:ref:`greenindex` で指定したインデックスを持つ相関関数が
mVMC/:math:`{\mathcal H}\Phi` によって計算され,
ファイルに出力される.
総合的な説明はmVMC/:math:`{\mathcal H}\Phi` のマニュアルを参照.
:ref:`tutorial` でのファイル名は
``output/zvo_cisajs_001.dat`` および ``output/zvo_cisajscktalt_001.dat`` (mVMC), 
``output/zvo_cisajs.dat`` および ``output/zvo_cisajscktalt.dat`` (:math:`{\mathcal H}\Phi`).

``greenr2k`` ユーティリティはこのファイルを読み込んで計算を行う.
この時, (スタンダードモードを使わず自分でインデックスを指定するなどにより)
:ref:`greenindex` で挙げたインデックスの相関関数のなかで欠けているものがある場合,
それを 0 として扱う.

.. _zvocorr:

*k* パス上での相関関数
----------------------

Fourier変換された相関関数(波数表示)が入っている.
ユーティリイティ ``greenr2k`` によって生成される.
:ref:`tutorial` でのファイル名は ``output/zvo_corr_eigen0.dat`` である.

::
   
   # k-length[1]
   # Orbital  1 to Orbital  1
   #  UpUp[   2,   3] (Re. Im.) DownDown[   4,   5]
   #  Density[   6,   7] SzSz[   8,   9] S+S-[  10,  11] S.S[  12,  13]
   0.00000E+00    0.88211E+00   -0.50000E-09    0.88211E+00    0.40000E-09 ... 
   0.25000E-01    0.87976E+00   -0.46625E-09    0.87976E+00    0.42882E-09 ...
   0.50000E-01    0.87276E+00   -0.42841E-09    0.87276E+00    0.45201E-09 ...
   :                                                               :

はじめに各カラムに出力されている量の説明がコメントとして書かれ,
それに続いて *k* 点の距離とそれぞれの相関関数の実部と虚部が書かれている.
      
.. _gnuplot:

gnuplot スクリプト
------------------

``greenr2k`` にて作成される.
gnuplotでこれを読み込むことでグラフ中に *k* 点のラベルを表示する.
ファイル名は ``kpath.gp`` である.

.. code-block:: gnuplot

   set xtics ('G'     0.00000, 'X'     0.50000, 'M'     1.00000, 'G'    1.70711)
   set ylabel 'Correlation function'
   set grid xtics lt 1 lc 0

.. _correlation:

運動量分布関数の等値面をプロットするためのFermiSurferファイル
-------------------------------------------------------------

``greenr2k`` にて作成される.
ファイル名は ``output/zvo_corr_eigen0.dat.frmsf``
