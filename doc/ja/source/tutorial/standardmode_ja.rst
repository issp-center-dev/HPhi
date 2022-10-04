スタンダードモード
==============================

Heisenberg模型
------------------------------

以下のチュートリアルはディレクトリ ::

 samples/old/CG/Heisenberg/
 
内で行います。
Heisenberg模型におけるサンプル入力ファイルは ::

 samples/old/CG/Heisenberg/stan.in

にあります。
この例では2次元正方格子のHeisenberg模型(最近接サイト間の反強磁性的スピン結合のみを持つ)
を考察します。

.. math::

   \mathcal H=J \sum_{i,j=1}^{4} ( S_{i j} \cdot  S_{i+1 j} +  S_{i j} \cdot  S_{i j+1})
..   \hat{\mathcal H}=J \sum_{i,j=1}^{4} (\hat{ S }_{i j} \cdot \hat{ S }_{i+1 j} + \hat{ S }_{i j} \cdot \hat{ S }_{i j+1},)

ただし、周期境界条件 :math:`(S_{15}=S_{51}= S_{11})` を採用します。
インプットファイルの中身は次のとおりです。 ::

 model = "Spin"
 method = "CG"
 lattice = "square"
 W = 4
 L = 4
 J = 1.0
 2Sz = 0

この例ではスピン結合 :math:`J=1` (任意単位)とし、サイト数は16としました。

**Log 出力**
^^^^^^^^^^^^^^^^^^^^^^^

標準出力でログが出力されます。
また、 \"output\" ディレクトリが自動生成され、
そこにも計算経過を示すログが出力されます。
例えば、サンプル実行時には以下のファイルが出力されます。 ::

 CHECK_InterAll.dat     Time_CG_EigenVector.dat  zvo_Lanczos_Step.dat  
 CHECK_Memory.dat       WarningOnTransfer.dat    zvo_sz_TimeKeeper.dat
 CHECK_Sdim.dat         zvo_TimeKeeper.dat
 
ログ出力されるファイルの詳細は :ref:`Subsec:checkchemi` 等をご覧ください。

実行コマンドと標準出力(MPI並列/ハイブリッド並列でコンパイルした場合の結果)は次のとおりです。::

 $ Path/HPhi -s stan.in
 
::


       ,ammmmmmmmmmmmmmb,,       Welcome to the
     ,@@` dm          mb  ===m
   ,@@` d@@@@@@@@@@@@@@@@b Pm,   @@          @@       @@
  d@  d@@@ @@@ @@@@@@ @@@@b ~@a  @@          @@    @@@@@@@@
 d@   @@@@ ^^^ @@@@ m m @@@   @, @@          @@  @@@  @@  @@@
 @    @@@@_@@@_@@@@mm mm@@@   @| @@mmmmmmmmmm@@ @@    @@    @@
 P@    9@@@@@@@@@@@@@@@@@P    @~ @@@@@@@@@@@@@@ @@    @@    @@
  @@      ~~9@@@@@@PPP~      @P  @@          @@  @@@  @@  @@@
   ~@@b      @@@@@@@      ,@@~   @@          @@    @@@@@@@@
     ~@@@m,,@@@@@@@@@  ,m@~`     @@          @@       @@
         ~~9@@@@@@@@@  ~
            9@P~~~9@P            Version 2.0.3


 #####  Parallelization Info.  #####

   OpenMP threads : 1
   MPI PEs : 1


 ######  Standard Intarface Mode STARTS  ######

   Open Standard-Mode Inputfile stan.in

   KEYWORD : model                | VALUE : Spin
   KEYWORD : method               | VALUE : CG
   KEYWORD : lattice              | VALUE : square
   KEYWORD : w                    | VALUE : 4
   KEYWORD : l                    | VALUE : 4
   KEYWORD : j                    | VALUE : 1.0
   KEYWORD : 2sz                  | VALUE : 0

 #######  Parameter Summary  #######

   @ Lattice Size & Shape

                 a = 1.00000     ######  DEFAULT VALUE IS USED  ######
           Wlength = 1.00000     ######  DEFAULT VALUE IS USED  ######
           Llength = 1.00000     ######  DEFAULT VALUE IS USED  ######
                Wx = 1.00000     ######  DEFAULT VALUE IS USED  ######
                Wy = 0.00000     ######  DEFAULT VALUE IS USED  ######
                Lx = 0.00000     ######  DEFAULT VALUE IS USED  ######
                Ly = 1.00000     ######  DEFAULT VALUE IS USED  ######
            phase0 = 0.00000     ######  DEFAULT VALUE IS USED  ######
            phase1 = 0.00000     ######  DEFAULT VALUE IS USED  ######

   @ Super-Lattice setting

                 L = 4
                 W = 4
            Height = 1           ######  DEFAULT VALUE IS USED  ######
          Number of Cell = 16

   @ Hamiltonian

                 h = 0.00000     ######  DEFAULT VALUE IS USED  ######
             Gamma = 0.00000     ######  DEFAULT VALUE IS USED  ######
                2S = 1           ######  DEFAULT VALUE IS USED  ######
                 D = 0.00000     ######  DEFAULT VALUE IS USED  ######
               J0x = 1.00000
               J0y = 1.00000
               J0z = 1.00000
               J1x = 1.00000
               J1y = 1.00000
               J1z = 1.00000

   @ Numerical conditions

        LargeValue = 4.50000     ######  DEFAULT VALUE IS USED  ######

 ######  Print Expert input files  ######

     locspn.def is written.
     coulombinter.def is written.
     hund.def is written.
     exchange.def is written.
     CDataFileHead = zvo         ######  DEFAULT VALUE IS USED  ######
       Lanczos_max = 2000        ######  DEFAULT VALUE IS USED  ######
        initial_iv = -1          ######  DEFAULT VALUE IS USED  ######
              exct = 1           ######  DEFAULT VALUE IS USED  ######
        LanczosEps = 14          ######  DEFAULT VALUE IS USED  ######
     LanczosTarget = 2           ######  DEFAULT VALUE IS USED  ######
            NumAve = 5           ######  DEFAULT VALUE IS USED  ######
     ExpecInterval = 20          ######  DEFAULT VALUE IS USED  ######
            NOmega = 200         ######  DEFAULT VALUE IS USED  ######
          OmegaMax = 72.00000    ######  DEFAULT VALUE IS USED  ######
          OmegaMin = -72.00000   ######  DEFAULT VALUE IS USED  ######
           OmegaIm = 0.04000     ######  DEFAULT VALUE IS USED  ######
               2Sz = 0
      modpara.def is written.

   @ Spectrum

        SpectrumQW = 0.00000     ######  DEFAULT VALUE IS USED  ######
        SpectrumQL = 0.00000     ######  DEFAULT VALUE IS USED  ######
        SpectrumQH = 0.00000     ######  DEFAULT VALUE IS USED  ######
      SpectrumType = szsz        ######  DEFAULT VALUE IS USED  ######
         pair.def is written.


   @ CalcMod

           Restart = none        ######  DEFAULT VALUE IS USED  ######
    InitialVecType = c           ######  DEFAULT VALUE IS USED  ######
        EigenVecIO = none        ######  DEFAULT VALUE IS USED  ######
          CalcSpec = none        ######  DEFAULT VALUE IS USED  ######
      calcmod.def is written.

       ioutputmode = 1           ######  DEFAULT VALUE IS USED  ######
     greenone.def is written.
     greentwo.def is written.
     namelist.def is written.

 ######  Input files are generated.  ######

   Read File 'namelist.def'.
   Read File 'calcmod.def' for CalcMod.
   Read File 'modpara.def' for ModPara.
   Read File 'locspn.def' for LocSpin.
   Read File 'coulombinter.def' for CoulombInter.
   Read File 'hund.def' for Hund.
   Read File 'exchange.def' for Exchange.
   Read File 'greenone.def' for OneBodyG.
   Read File 'greentwo.def' for TwoBodyG.
   Read File 'pair.def' for PairExcitation.

 ######  Definition files are correct.  ######

   Read File 'locspn.def'.
   Read File 'coulombinter.def'.
   Read File 'hund.def'.
   Read File 'exchange.def'.
   Read File 'greenone.def'.
   Read File 'greentwo.def'.
   Read File 'pair.def'.

 ######  Indices and Parameters of Definition files(*.def) are complete.  ######

   MAX DIMENSION idim_max=12870
   APPROXIMATE REQUIRED MEMORY  max_mem=0.001647 GB


 ######  MPI site separation summary  ######

   INTRA process site
     Site    Bit
        0       2
        1       2
        2       2
        3       2
        4       2
        5       2
        6       2
        7       2
        8       2
        9       2
       10       2
       11       2
       12       2
       13       2
       14       2
       15       2

   INTER process site
     Site    Bit

   Process element info
     Process       Dimension   Nup  Ndown  Nelec  Total2Sz   State
           0           12870     8      8      8         0

    Total dimension : 12870


 ######  LARGE ALLOCATE FINISH !  ######

   Start: Calculate HilbertNum for fixed Sz.
   End  : Calculate HilbertNum for fixed Sz.

   Start: Calculate diagaonal components of Hamiltonian.
   End  : Calculate diagaonal components of Hamiltonian.

 ######  Eigenvalue with LOBPCG  #######

   initial_mode=1 (random): iv = -1 i_max=12870 k_exct =1

     Step   Residual-2-norm     Threshold      Energy
         1     2.44343e+00     1.00000e-07          -5.27456e-01
         2     2.76604e+00     1.87217e-07          -1.87217e+00
         3     2.61923e+00     4.19088e-07          -4.19088e+00
         4     2.57106e+00     5.97098e-07          -5.97098e+00

 ( snip )

        40     7.39431e-06     1.12285e-06          -1.12285e+01
        41     4.15948e-06     1.12285e-06          -1.12285e+01
        42     2.04898e-06     1.12285e-06          -1.12285e+01
        43     9.92048e-07     1.12285e-06          -1.12285e+01

 ######  End  : Calculate Lanczos EigenValue.  ######


 ######  End  : Calculate Lanczos EigenVec.  ######

 i=    0 Energy=-11.228483 N= 16.000000 Sz=  0.000000 Doublon=  0.000000

この実行では、はじめにハミルトニアンの詳細を記述するファイル
(``locspin.def`` , ``trans.def`` , ``exchange.def`` , ``coulombintra.def`` , ``hund.def`` , ``namelist.def`` , ``calcmod.def`` , ``modpara.def`` ) と、結果として出力する相関関数の要素を指定するファイル( ``greenone.def`` , ``greentwo.def`` ) が生成されます。これらのファイルはエキスパートモードと共通です。

**計算結果出力**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**局所最適ブロック共役勾配(LOBCG)法**
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

入力ファイルで ``method = "CG"`` を選択すると、LOBCG法での計算が行われます。
LOBCG法での計算が正常終了すると、固有エネルギーおよび一体グリーン関数、二体グリーン関数が計算され、ファイル出力されます。
以下に、このサンプルでの出力ファイル例を記載します。
(xxには0から始まる固有値番号が入ります)。 ::
 
 zvo_energy.dat
 zvo_cisajscktalt_eigen_xx.dat  zvo_phys_Nup4_Ndown4.dat


スタンダードモードの場合は、\"greenone.def\"、\"greentwo.def\"に基づき、::

 zvo_cisajs_eigen_xx.dat、zvo_cisajscktalt_eigen_xx.dat

に固有値番号に対応した一体グリーン関数および二体グリーン関数の値が出力されます。
 

**Lanczos法**
""""""""""""""""""

Lanczos法での計算が正常終了すると、固有エネルギーおよび一体グリーン関数、二体グリーン関数が計算され、ファイル出力されます。 ::
 
 zvo_energy.dat zvo_cisajs.dat
 zvo_cisajscktalt.dat


スタンダードモードの場合は、``greenone.def`` 、``greentwo.def`` に基づき、
一体グリーン関数には :math:`\langle n_{i\sigma} \rangle` 、
二体グリーン関数には :math:`\langle n_{i\sigma} n_{j\sigma'} \rangle` が自動出力されます。なお、Lanczos法で求めた固有ベクトルが十分な精度を持つ場合には
その固有ベクトルで計算されます。
一方、Lanczos法で求めた固有ベクトルが十分な精度を持たない場合には、
ログ出力に「Accuracy of Lanczos vetor is not enough」が表示され、
CG法で固有ベクトルが求められます。
各ファイルの詳細は、セクション :ref:`Subsec:energy.dat` , :ref:`Subsec:cgcisajs` , :ref:`Subsec:cisajscktalt` に記載がありますので、ご参照ください。

**TPQ法**
""""""""""""""

入力ファイルで ``method = "TPQ"`` を選択すると、TPQ法での計算が行われます。
TPQ法での計算が正常終了すると、以下のファイルが出力されます
(\%\%にはrunの回数、\&\&にはTPQのステップ数が入ります)。::
  
 Norm_rand%%.dat SS_rand%%.dat
 zvo_cisajs_set%%step&&.dat
 zvo_cisajscktalt_set%%step&&.dat
 
Norm\_rand\%\%.datには、逆温度や波動関数の規格前の大きさなどの基礎情報が、各run回数に応じステップ数とともに出力されます。また、SS\_rand\%\%.datには、逆温度、エネルギー、ハミルトニアンの二乗の期待値などの物理量が各run回数に応じステップ数とともに出力されます。zvo\_cisajs\_set\%\%step\&\&.datとzvo\_cisajscktalt\_set\%\%step\&\&.datには各run回数でのステップ数に応じた一体グリーン関数および二体グリーン関数が出力されます。各ファイルの詳細はそれぞれ、セクション :ref:`Subsec:normrand`, :ref:`Subsec:ssrand`, :ref:`Subsec:cgcisajs`, :ref:`Subsec:cisajscktalt` に記載がありますので、ご参照ください。

**全対角化法**
"""""""""""""""""""""""""""""""

入力ファイルで ``method = "fulldiag"`` を選択すると、全対角化法での計算が行われます。
全対角化法での計算が正常終了すると、下記のファイルが出力されます(xxには0から始まる固有値番号が入ります)。::
 
 Eigenvalue.dat zvo_cisajs_eigen_xx.dat
 zvo_cisajscktalt_eigen_xx.dat  zvo_phys_Nup4_Ndown4.dat

Eigenvalue.datには固有値番号およびエネルギー固有値が出力されます。また、zvo\_cisajs\_eigen\_xx.dat、zvo\_cisajscktalt\_eigen\_xx.datには固有値番号に対応した一体グリーン関数および二体グリーン関数の値が出力されます。また、zvo\_phys\_Nup4\_Ndown4.dat, physical quantitiesには、エネルギーやダブロンの期待値などの物理量が出力されます。各ファイルの詳細は、それぞれ :ref:`Subsec:eigenvalue` - :ref:`Subsec:cisajscktalt` に記載がありますので、ご参照ください。

その他の系でのチュートリアル
---------------------------------

 ``samples`` 以下にはこの他にも様々なチュートリアルが置いてあります。それぞれのチュートリアルの内容や手順については、各フォルダにある ``README.md`` または `チュートリアル向けの説明資料 <https://issp-center-dev.github.io/HPhi/manual/develop/tutorial/en/html/index.html>`_ をご覧ください。
