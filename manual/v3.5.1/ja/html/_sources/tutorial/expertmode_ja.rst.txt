.. highlight:: none

エキスパートモード
============================

エキスパートモードでは、入力ファイルとして

1. 入力ファイルリスト
2. 基本パラメータ用ファイル
3. Hamiltonian作成用ファイル
4. 出力結果指定用ファイル

を用意した後、計算を行います。計算開始後に関しては、スタンダードモードと同様です。
ここでは前節のスタンダードモードでのチュートリアルを行った後の状況を例に入力ファイルの作成に関する説明を行います。

入力ファイルリストファイル
------------------------------------------

入力ファイルの種類と名前を指定するファイルnamelist.defには、下記の内容が記載されています。
各入力ファイルリストファイルでは、行毎にKeywordとファイル名を記載し、ファイルの種類の区別を行います。
詳細はセクション :ref:`Subsec:InputFileList` をご覧ください。 ::

        ModPara  modpara.def
        LocSpin  locspn.def
   CoulombInter  coulombinter.def
           Hund  hund.def
       Exchange  exchange.def
       OneBodyG  greenone.def
       TwoBodyG  greentwo.def
        CalcMod  calcmod.def
 PairExcitation  pair.def
    SpectrumVec  zvo_eigenvec_0
    
基本パラメータの指定
--------------------------

計算モード、計算用パラメータ、局在スピンの位置を以下のファイルで指定します。

**計算モードの指定**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CalcModでひも付けられるファイル(ここではcalcmod.def)で計算モードを指定します。
ファイルの中身は下記の通りです。
計算モード、計算用パラメータ、局在スピンの位置を以下のファイルで指定します。::

 #CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 3:CG, ...
 #CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, ...
 #Restart = 0:None, 1:Save, 2:Restart&Save, 3:Restart
 #CalcSpec = 0:None, 1:Normal, 2:No H*Phi, 3:Save, ...
 CalcType   3
 CalcModel   1
 ReStart   0
 CalcSpec   0
 CalcEigenVec   0
 InitialVecType   0
 InputEigenVec   0

CalcTypeで計算手法の選択、CalcModelで対象モデルの選択を行います。
ここでは、計算手法としてLOBCG法、対象モデルとしてスピン系(カノニカル)を選択しています。
CalcModファイルでは固有ベクトルの入出力機能も指定することができます。
CalcModファイルの詳細は :ref:`Subsec:calcmod` をご覧ください。


**計算用パラメータの指定**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ModParaでひも付けられるファイル(ここではmodpara.def)で計算用パラメータを指定します。ファイルの中身は下記の通りです。::

 --------------------
 Model_Parameters   0
 --------------------
 HPhi_Cal_Parameters
 --------------------
 CDataFileHead  zvo
 CParaFileHead  zqp
 --------------------
 Nsite          16   
 2Sz            0    
 Lanczos_max    2000 
 initial_iv     -1   
 exct           1    
 LanczosEps     14   
 LanczosTarget  2    
 LargeValue     4.500000000000000e+00    
 NumAve         5    
 ExpecInterval  20   
 NOmega         200  
 OmegaMax       7.200000000000000e+01     4.000000000000000e-02    
 OmegaMin       -7.200000000000000e+01    4.000000000000000e-02    
 OmegaOrg       0.0 0.0

このファイルでは、サイト数、伝導電子の数、トータル :math:`S_z` やLanczosステップの最大数などを指定します。ModParaファイルの詳細はセクション :ref:`Subsec:modpara` をご覧ください。
  
**局在スピンの位置の指定**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

LocSpinでひも付けられるファイル(ここではlocspn.def)で局在スピンの位置と :math:`S` の値を指定します。ファイルの中身は下記の通りです。::

 ================================
 NlocalSpin    16  
 ================================ 
 ========i_0IteElc_1LocSpn ====== 
 ================================ 
     0      1
     1      1
     2      1
     3      1
     4      1
     5      1
 ...
 
LocSpinファイルの詳細は :ref:`Subsec:locspn` をご覧ください。
 

Hamiltonianの指定
----------------------------------

基本パラメータを設定した後は、Hamiltonianを構築するためのファイルを作成します。
 :math:`{\mathcal H}\Phi` では、電子系の表現で計算を行うため、スピン系では以下の関係式

.. math::

    S^z_{i}&=(c_{i\uparrow}^{\dagger}c_{i\uparrow}-c_{i\downarrow}^{\dagger}c_{i\downarrow})/2,\\
    S^+_{i}&=c_{i\uparrow}^{\dagger}c_{i\downarrow},\\
    S^-_{i}&=c_{i\downarrow}^{\dagger}c_{i\uparrow}.

を用い、電子系の演算子に変換しHamiltonianの作成をする必要があります。


**Transfer部の指定**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Transでひも付けられるファイル(ここではzTrans.def)で電子系のTransferに相当するHamiltonian

.. math::

   \mathcal{H} +=-\sum_{ij\sigma_1\sigma2}
   t_{ij\sigma_1\sigma2}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}.
   
を指定します。ファイルの中身は下記の通りです。::

 ======================== 
 NTransfer       0  
 ======================== 
 ========i_j_s_tijs====== 
 ======================== 

スピン系では外場を掛ける場合などに使用することができます。
例えば、サイト1に :math:`-0.5 S_z^{(1)}` ( :math:`S=1/2` )の外場を掛けたい場合には、
電子系の表現 :math:`-0.5/2(c_{1\uparrow}^{\dagger}c_{1\uparrow}-c_{1\downarrow}^{\dagger}c_{1\downarrow})` に書き換えた以下のファイルを作成することで計算することが出来ます。 ::

 ======================== 
 NTransfer      1   
 ======================== 
 ========i_j_s_tijs====== 
 ======================== 
 1 0 1 0 -0.25 0
 1 1 1 1 0.25 0
 
Transファイルの詳細はセクション :ref:`Subsec:Trans` をご覧ください。

**二体相互作用部の指定**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

InterAllでひも付けられるファイル(ここではzInterAll.def)で電子系の二体相互作用部に相当するHamiltonian

.. math::

   \mathcal{H}+=\sum_{i,j,k,l}\sum_{\sigma_1,\sigma_2, \sigma_3, \sigma_4}I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}.

を指定します。ファイルの中身は下記の通りです。 ::

 ====================== 
 NInterAll      96  
 ====================== 
 ========zInterAll===== 
 ====================== 
     0     0     0     0     1     0     1     0   0.500000  0.000000
     0     0     0     0     1     1     1     1  -0.500000  0.000000
     0     1     0     1     1     0     1     0  -0.500000  0.000000
     0     1     0     1     1     1     1     1   0.500000  0.000000
     0     0     0     1     1     1     1     0   1.000000  0.000000
     0     1     0     0     1     0     1     1   1.000000  0.000000
 ...

ここでは、簡単のためサイト :math:`i` とサイト :math:`i+1` 間の相互作用に着目して説明します。
:math:`S = 1/2` の場合、相互作用の項をフェルミオン演算子で書き換えると、 

.. math::
   \mathcal{H}_{i,i+1}&=J(S^x_{i}S^x_{i+1}+S^y_{i}S^y_{i+1}+S^z_{i}S^z_{i+1}) \nonumber\\
   &=J \left( \frac{1}{2}S^+_{i}S^-_{i+1}+\frac{1}{2}S^-_{i}S^+_{i+1}+S^z_{i}S^z_{i+1} \right) \nonumber\\
   &=J \left[ \frac{1}{2}c_{i\uparrow}^{\dagger}c_{i\downarrow}c_{i+1\downarrow}^{\dagger}c_{i+1\uparrow}+\frac{1}{2}c_{i\downarrow}^{\dagger}c_{i\uparrow}c_{i+1\uparrow}^{\dagger}c_{i+1\downarrow}+\frac{1}{4}(c_{i\uparrow}^{\dagger}c_{i\uparrow}-c_{i\downarrow}^{\dagger}c_{i\downarrow})(c_{i+1\uparrow}^{\dagger}c_{i+1\uparrow}-c_{i+1\downarrow}^{\dagger}c_{i+1\downarrow}) \right]. \nonumber 

となります。したがって、 :math:`J=2` に対してInterAllファイルのフォーマットを参考に相互作用を記載すると、
 :math:`S^z_{i}S^z_{i+1}` の相互作用は ::

    i     0     i     0    i+1     0    i+1     0   0.500000  0.000000
    i     0     i     0    i+1     1    i+1     1  -0.500000  0.000000
    i     1     i     1    i+1     0    i+1     0  -0.500000  0.000000
    i     1     i     1    i+1     1    i+1     1   0.500000  0.000000
  
となり、それ以外の項は ::

    i     0     i     1    i+1     1    i+1     0   1.000000  0.000000
    i     1     i     0    i+1     0    i+1     1   1.000000  0.000000
  
と表せばよいことがわかります。なお、InterAll以外にも、Hamiltonianを簡易的に記載するための下記のファイル形式に対応しています。
詳細はセクション :ref:`Subsec:interall` - :ref:`Subsec:pairlift` をご覧ください。

出力ファイルの指定
-------------------------

一体Green関数および二体Green関数の計算する成分を、それぞれOneBodyG, TwoBodyGでひも付けられるファイルで指定します。 

**一体Green関数の計算対象の指定**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OneBodyGでひも付けられるファイル(ここではgreenone.def)で計算する一体Green関数  :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2} \rangle` の成分を指定します。ファイルの中身は下記の通りです ::

 ===============================
 NCisAjs         32
 ===============================
 ======== Green functions ======
 ===============================
    0     0     0     0
    0     1     0     1
    1     0     1     0
    1     1     1     1
    2     0     2     0
 ...
 
一体Green関数計算対象成分の指定に関するファイル入力形式の詳細はセクション :ref:`Subsec:onebodyg` をご覧ください。

**二体Green関数の計算対象の指定**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TwoBodyGでひも付けられるファイル(ここではgreentwo.def)で計算する二体Green関数 :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4} \rangle` の成分を指定します。ファイルの中身は下記の通りです。 ::

 =============================================
 NCisAjsCktAltDC       1024
 =============================================
 ======== Green functions for Sq AND Nq ======
 =============================================
    0     0     0     0     0     0     0     0
    0     0     0     0     0     1     0     1
    0     0     0     0     1     0     1     0
    0     0     0     0     1     1     1     1
    0     0     0     0     2     0     2     0
 ...

二体Green関数計算対象成分の指定に関するファイル入力形式の詳細はセクション :ref:`Subsec:twobodyg` をご覧ください。

計算の実行
--------------------------

全ての入力ファイルが準備できた後、計算実行します。実行時はエキスパートモードを指定する \"-e\" をオプションとして指定の上、入力ファイルリストファイル(ここではnamelist.def)を引数とし、ターミナルから :math:`{\mathcal H}\Phi` を実行します。 ::

 $ Path/HPhi -e namelist.def

計算開始後のプロセスは、スタンダードモードと同様になります。
