各ユーティリティの動作について
==============================

``fourier`` ユーティリティ
--------------------------

このユーティリティーは, 次のようにして使う.

.. code-block:: bash

   $ ${PATH}/fourier ${NAMELIST} ${GEOMETRY}

ここで, ``${PATH}`` は ``fourier`` ユーティリティのバイナリのあるディレクトリのパス,
${NAMELIST}は :math:`{\mathcal H}\Phi`/mVMC の NameList インプットファイル名,
${GEOMETRY}は :ref:`geometry` ファイルへのパスである.

:math:`{\mathcal H}\Phi` の各モード
(Lanczos, TPQ, 全対角化, LOBCG)および mVMC のどの計算で得られた
相関関数のFourier変換を行うかによって, 動作が若干異なる.
以下では ModPara インプットファイルの ``CDataFileHead`` が
``"zvo"`` (デフォルト値)であるとする.

HPhi-Lanczos
~~~~~~~~~~~~

この場合に ``HPhi`` が ``output/`` ディレクトリに出力するサイト表示の相関関数は,
``zvo_cisajs.dat`` (1体), ``zvo_cisajscktalt.dat`` (2体)である.
``fourier`` ユーティリティーは, これらを読み込みFourier変換を行った後,
単一のファイル ``zvo_corr.dat`` を ``output/`` ディレクトリに出力する.

HPhi-TPQ
~~~~~~~~

この場合に ``HPhi`` は, 各試行/TPQステップ毎に
``zvo_cisajs_run*step*.dat`` (1体), ``zvo_cisajscktalt_run*step*.dat`` (2体)というファイルを
``output/`` ディレクトリに出力する.
``fourier`` ユーティリティーは, 各試行/TPQステップ毎に
1体および2体の相関関数を読み込みFourier変換を行った後,
``zvo_corr_run*step*.dat`` という名前のファイルとして ``output/`` ディレクトリに出力する.

HPhi-全対角化およびLOBCG
~~~~~~~~~~~~~~~~~~~~~~~~

この場合に ``HPhi`` は, 各波動関数ごとに
``zvo_cisajs_eigen*.dat`` (1体), ``zvo_cisajscktalt_eigen*.dat`` (2体)というファイルを
``output/`` ディレクトリに出力する.
``fourier`` ユーティリティーは, 各波動関数ごとに
1体および2体の相関関数を読み込みFourier変換を行った後,
``zvo_corr_eigen*.dat`` という名前のファイルとして ``output/`` ディレクトリに出力する.

mVMC
~~~~

この場合に ``vmc.out`` は, ``ModPara`` インプットファイルで指定された
``NDataIdxStart`` および ``NDataQtySmp`` というパラメーターに応じて
試行を行いインデックスをつけられた
``zvo_cisajs_???.dat`` (1体), ``zvo_cisajscktalt_???.dat`` (2体)というファイルを
``output/`` ディレクトリに出力する.
``fourier`` ユーティリティーはそれらのファイルを読み込み, 
各試行に対してFourier変換を行った後,
それらの実部, 虚部ごとに平均値

.. math::

   \begin{align}
   \langle A \rangle = \frac{1}{N_{\rm Try}} \sum_{i=1}^{N_{\rm Try}} A_i
   \end{align}

および標準誤差

.. math::
   
   \begin{align}
   \delta A = \frac{1}{N_{\rm Try} - 1}
   \sqrt{\frac{1}{N_{\rm Try}} \sum_{i=1}^{N_{\rm Try}} (A_i - \langle A \rangle)^2}
   \end{align}

を計算し, 平均値と誤差を含んだ単一のファイル
``zvo_corr_eigen*.dat`` を ``output/`` ディレクトリに出力する.

``corplot`` ユーティリティ
--------------------------

このユーティリティーは, 次のようにして使う.

.. code-block:: bash

   $ ${PATH}/corplot ${CORR1} ${CORR2} ${CORR3} ...

ここで, ``${PATH}`` は ``corplot`` ユーティリティのバイナリのあるディレクトリのパス,
${CORR1}, ${CORR2}, ${CORR3}, ...は ``fourier`` ルーチンによって生成された
:ref:`zvocorr` ファイルへのパスである.
すなわち, このユーティリティーでは, (TPQ計算による温度依存性を調べる等の用途で)
複数の相関関数ファイルを読み込み, それらを同時にプロットすることができる.
