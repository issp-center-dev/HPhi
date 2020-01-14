スタンダードモードの入力パラメーター
====================================

以下に入力ファイルの例を示す.

:download:`stan.in <../../../../samples/Wannier/Sr2VO4/stan.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/stan.in

Wannier関数を用いたダウンフォールディングに特有のパラメーター設定は次の通りである.
                    
- 格子

  * ``lattice = "wannier90"``

- 格子サイズ関連のパラメータ
   
  *  ``W``, ``L``, ``Height``

     **形式 :** 自然数.

     **説明 :** 標準の単位胞の並び方を指定する.
     
  
  *  ``a0W``, ``a0L``, ``a0H``, ``a1W``, ``a1L``, ``a1H``, ``a2W``, ``a2L``, ``a2H``
   
     **形式 :** 整数.

     **説明 :**
     格子を指定する3本のベクトル( :math:`{\vec a}_0, {\vec a}_1, {\vec a}_2`)を指定する.
     これらのベクトルは標準の並進ベクトルを基底とした座標(Fractional coordinate) で指定される.
   
- 副格子サイズ関連のパラメータ
  
  *  ``Wsub``, ``Lsub``, ``Hsub``
   
     **形式 :** 自然数.デフォルトでは ``Wsub=W``, ``Lsub=L``, ``Hsub=Height`` となる.

     **説明 :**
     mVMC でのみ利用可能.
     変分波動関数のペア軌道部分に副格子を持たせるためのパラメータで,副格子のサイズを与える.
     元の計算セルが副格子に整合しない場合にはプログラムを終了する.
   
   
  *  ``a0Wsub``, ``a0Lsub``, ``a0Hsub``, ``a1Wsub``, ``a1Lsub``, ``a1Hsub``, ``a2Wsub``, ``a2Lsub``, ``a2Hsub``
   
     **形式 :** 自然数. デフォルトでは ``a0Wsub=a0W``, ``a0Lsub=a0L``, ``a0Hsub=a0H``,
     ``a1Wsub=a1W``, ``a1Lsub=a1L``, ``a1Hsub=a1H``, ``a2Wsub=a2W``, ``a2Lsub=a2L``, ``a2Hsub=a2H`` となる.

     **説明 :** これらのパラメーターの指定の仕方は ``a0W``, ``a0L``, ``a0H``,
     ``a1W``, ``a1L``, ``a1H``,   ``a2W``, ``a2L``, ``a2H`` と同様である.
     ただし,元の計算セルが副格子に整合しない場合にはプログラムを終了する.


- 相互作用の制御関連パラメータ
     
  * ``lambda_U``

     **形式 :** 実数 (0以上)

     **デフォルト値 :** ``1.0``

     クーロン積分の大きさを :math:`\lambda_U` 倍にして調整するパラメータ.

  * ``lambda_J``
	 
     **形式 :** 実数 (0以上)

     **デフォルト値 :** ``1.0``

     交換積分の大きさを :math:`\lambda_J` 倍にして調整するパラメータ.

  * ``lambda``

     **形式 :** 実数 (0以上)

     **デフォルト値 :** ``1.0``

     クーロン積分, 交換積分の大きさを :math:`\lambda` 倍にして調整するパラメータ.
     :math:`\lambda_U` , :math:`\lambda_J` が定義されている場合には, そちらの値を優先する. 
     
     
  * ``cutoff_t``, ``cutoff_u``, ``cutoff_j``

     **形式 :** 実数

     **デフォルト値 :** ``1.0e-8``

     ホッピング, クーロン積分, 交換積分に対して,
     これより小さい値を無視する.

  * ``cutoff_tW``, ``cutoff_tL``, ``cutoff_tH``
  * ``cutoff_UW``, ``cutoff_UL``, ``cutoff_UH``
  * ``cutoff_JW``, ``cutoff_JL``, ``cutoff_JH``

     **形式 :** 整数.

     **デフォルト値 :** ``cutoff_tW = int((W-1)/2)``, ``cutoff_tL=int((L-1)/2)``, ``cutoff_tH=int((Height-1)/2)`` に指定される(ただし, ``W``, ``L``, ``Height`` が指定されていない場合は0). それ以外は0.

     ホッピング, Coulomb積分, 交換積分に対して,
     これらの値を越える並進ベクトル :math:`{\bf R}` を持つものを無視するようにする.

  * ``cutoff_length_t``, ``cutoff_length_U``, ``cutoff_length_J``

     **形式 :** 実数.

     **デフォルト値 :** ``cutoff_length_t = -1.0`` (すべてのレンジの項を含む), それ以外は0.3.

     ホッピング, Coulomb積分, 交換積分に対して,この距離を超えるものを無視する.
     距離はワニエ関数の中心座標と単位格子ベクトルから算出される.

- 一体補正に関するパラメータ
  
  一体補正では一体項に対して下記の項を差し引くことで、模型を解く際のダブルカウンティングを避けることが可能となる.
  
       .. math::
	  \begin{aligned}
	  t_{mm}^{\rm DC}({\bf 0}) &\equiv \alpha U_{mm}({\bf 0}) D_{mm}({\bf 0})
	  + \sum_{({\bf R}, n) \neq ({\bf 0}, m)} U_{m n} ({\bf R})D_{nn}({\bf 0})\\
	  & - (1-\alpha) \sum_{({\bf R}, n) \neq ({\bf 0}, 0)} J_{m n}({\bf R}) D_{nn}({\bf R}),\\
	  t_{mn}^{\rm DC}({\bf R}_{ij}) &\equiv \frac{1}{2} J_{mn}({\bf R}_{ij}) \left(D_{nm}({\bf R}_{ji}) + 2 {\rm Re} [D_{nm}({\bf R}_{ji})]\right)\\
	  &-\frac{1}{2}  U_{mn}({\bf R}_{ij}) D_{nm}({\bf R}_{ji}),
	  \quad ({\bf R}_{ij}, m) \neq ({\bf 0}, n),
	  \\
	  D_{mn}({\bf R}_{ij}) &\equiv \sum_{\sigma}
	  \left\langle c_{im \sigma}^{\dagger} c_{jn \sigma}\right\rangle_{\rm KS},
	  \end{aligned}

  ここで, 第一項はHartree補正、第二項はFock補正を表す. :math:`\alpha` はオンサイト相互作用の寄与を調整するパラメータを表す.
	

  * ``doublecounting``

     **形式 :** char型

     **デフォルト値 :** ``none``

     ``none``: 一体補正を行わない.  ``Hartree_U``: クーロン積分 :math:`U_{Rii}` のみを考慮したHartree補正を行う. ``Hartree``: 通常のHartree補正を行う. ``full``: Fock項も含んだ一体補正を行う. 
     電子密度 :math:`D_{Rij}` に対してはRESPACKで出力した密度に関するファイル ``[CDataFileHead]_dr.dat`` に記載された値を採用する.ただし, ``[CDataFileHead]_dr.dat`` では, サイトあたりの電荷数が出力されているため,
     電荷密度にスピン依存性がないと仮定し処理している. 
     
  * ``alpha``

     **形式 :** 実数

     **デフォルト値 :** ``0.5``

     一体補正のうち, オンサイト相互作用の寄与を調整するパラメータ (:math:`0\le \alpha \le 1`). 

    
