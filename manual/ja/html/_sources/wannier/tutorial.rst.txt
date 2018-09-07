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

次にRESPACKに付属のユーティリティー ``qe2respack.sh`` を使う.
引数は ``[prefix].save`` ディレクトリ名.

.. code-block:: bash

   $ qe2respack.sh sr2cuo3.save
                
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

HPhi/mVMCによるモデル計算
-------------------------

まず, RESPACKによって出力されたホッピング等のファイルを
Wannier90の形式に変換する.
これにはHPhi/mVMCに付属のユーティリティー
respack2wan90.pyを使う.
引数はHPhi/mVMCのスタンダードモードのパラメーター ``CDataFileHead`` と同じにする.
引数を指定しない場合は ``zvo`` (HPhi/mVMCの ``CDataFileHead`` のデフォルト)
が指定されたものとする.

.. code-block:: bash

   $ respack2wan90.py zvo

これで, HPhi/mVMCで実行する準備が出来たので,
スタンダードモードを用いて計算する.   
                    
:download:`respack.in <../../../../samples/Wannier/Sr2VO4/stan.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/stan.in

.. code-block:: bash

   $ vmc.out -s stan.in
