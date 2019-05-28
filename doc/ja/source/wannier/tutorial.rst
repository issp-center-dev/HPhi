.. _tutorialwannier:

チュートリアル
==============

このチュートリアルでは Sr\ :sub:`2`\ VO\ :sub:`4`
を2次元3軌道Hubbardモデルにダウンフォールドして,
それをHPhi/mVMCで計算する.
DFT計算はQuantumESPRESSOで行う.

電荷密度のSCF計算
-----------------

まず, DFTによる電荷密度のSCF計算を行う.

:download:`scf.in <../../../../samples/Wannier/Sr2VO4/scf.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/scf.in

擬ポテンシャル (UPF ファイル) は
`The SG15 Optimized Norm-Conserving Vanderbilt (ONCV) pseudopotentials <www.quantum-simulation.org/potentials/sg15_oncv/>`_ のものを使う.

http://www.quantum-simulation.org/potentials/sg15_oncv/sg15_oncv_upf_2015-10-07.tar.gz

QuantumESPRESSO内のプログラム ``pw.x`` を使う.

.. code-block:: bash

   $ pw.x -in scf.in

(Optional) バンド計算と描画
---------------------------

:download:`band.in <../../../../samples/Wannier/Sr2VO4/band.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/band.in

ここでも ``pw.x`` を使う.
                    
.. code-block:: bash

   $ pw.x -in band.in

:download:`bands.in <../../../../samples/Wannier/Sr2VO4/bands.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/bands.in

QuantumESPRESSOの ``bands.x`` を使う.
                    
.. code-block:: bash

   $ bands.x -in bands.in

出力された ``bands.out.gnu`` をGnuPlotなどで読み込んでバンドを描く.
   
Kohn-Sham軌道の計算
-------------------

:download:`nscf.in <../../../../samples/Wannier/Sr2VO4/nscf.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/nscf.in

ここでも ``pw.x`` を使う.
                    
.. code-block:: bash

   $ pw.x -in nscf.in

次にRESPACKに付属のユーティリティー ``qe2respack.py`` を使う.
引数は ``[prefix].save`` ディレクトリ名.

.. code-block:: bash

   $ qe2respack.py sr2vo4.save
                
Wannier関数, 誘電関数, 有効相互作用の計算
-----------------------------------------

:download:`respack.in <../../../../samples/Wannier/Sr2VO4/respack.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/respack.in

RESPACKのプログラム ``calc_wannier``, ``calc_chiqw``, ``calc_j3d``,
``calc_w3d`` を使う.                    
                    
.. code-block:: bash

   $ calc_wannier < respack.in
   $ calc_chiqw < respack.in
   $ calc_w3d < respack.in
   $ calc_j3d < respack.in

これにより, RESPACKによって出力されたホッピング等のファイルが,
Wannier90の形式で ``dir-mvmc`` フォルダに格納される.
(RESPACKの次のバージョンでは,  ``dir-model`` にディレクトリ名が変更予定)

HPhi/mVMCによるモデル計算
-------------------------

HPhi/mVMCのスタンダードモードを利用することで,
``dir-mvmc`` のファイルを読み込み該当したモデルの計算ができる.
最初に ``dir-mvmc`` 以下のファイル一式を, 実行するディレクトリに移したあとに,
スタンダードモードで計算実行を行えばよい.
例えば, mVMCの場合は以下のコマンドを打つことで計算が実行される(HPhiでもほぼ同様).
                    
:download:`respack.in <../../../../samples/Wannier/Sr2VO4/stan.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/stan.in

.. code-block:: bash

   $ cp ./dir-mvmc/* .
   $ vmc.out -s stan.in
