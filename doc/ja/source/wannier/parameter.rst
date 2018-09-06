スタンダードモードの入力パラメーター
====================================

以下に入力ファイルの例を示す.

:download:`respack.in <../../../../samples/Wannier/Sr2VO4/stan.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/stan.in

Wannier関数を用いたダウンフォールディングに特有のパラメーター設定は次の通りである.
                    
* ``lattice = "wannier90"``

* ``cutoff_t``, ``cutoff_u``, ``cutoff_j``

   **形式 :** 実数

   **デフォルト値 :** ``1.0e-8``

   ホッピング, Coulomb積分, 交換積分に対して,
   これより小さい値を無視する.

*  ``W``, ``L``, ``Height``
*  ``a0W``, ``a0L``, ``a0H``, ``a1W``, ``a1L``, ``a1H``, ``a2W``, ``a2L``, ``a2H``
*  ``Wsub``, ``Lsub``, ``Hsub``
*  ``a0Wsub``, ``a0Lsub``, ``a0Hsub``, ``a1Wsub``, ``a1Lsub``, ``a1Hsub``, ``a2Wsub``, ``a2Lsub``, ``a2Hsub``

   このように3番目の次元が現れる.
