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

* ``cutoff_tW``, ``cutoff_tL``, ``cutoff_tH``
* ``cutoff_UW``, ``cutoff_UL``, ``cutoff_UH``
* ``cutoff_JW``, ``cutoff_JL``, ``cutoff_JH``

   **形式 :** 整数

   **デフォルト値 :** すべてのレンジの項を含む

   ホッピング, Coulomb積分, 交換積分に対して, これらの値を越える並進ベクトル
   :math:`{\bf R}` を持つものを無視するようにする.

* ``cutoff_length_t``, ``cutoff_length_U``, ``cutoff_length_J``

   **形式 :** 実数

   **デフォルト値 :** -1.0 (すべてのレンジの項を含む)

   ホッピング, Coulomb積分, 交換積分に対して,
   この距離を超えるものを無視する.
   距離はワニエ関数の中心座標と単位格子ベクトルから求める.

*  ``W``, ``L``, ``Height``
*  ``a0W``, ``a0L``, ``a0H``, ``a1W``, ``a1L``, ``a1H``, ``a2W``, ``a2L``, ``a2H``
*  ``Wsub``, ``Lsub``, ``Hsub``
*  ``a0Wsub``, ``a0Lsub``, ``a0Hsub``, ``a1Wsub``, ``a1Lsub``, ``a1Hsub``, ``a2Wsub``, ``a2Lsub``, ``a2Hsub``

   このように3番目の次元が現れる.
