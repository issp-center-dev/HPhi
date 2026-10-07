.. highlight:: none

.. _Subsec:TransSym:

TransSym指定ファイル
~~~~~~~~~~~~~~~~~~~~

計算を一つの対称性セクターに制限します。サイト置換の群 :math:`G` と、
各要素 :math:`g` に対する1次元指標 :math:`\chi(g)` （絶対値1の複素数）を
指定します。 :math:`{\mathcal H}\Phi` は対称化された状態

.. math:: |r;\chi\rangle \propto \sum_{g\in G}\chi(g)^{*}\,T_g|r\rangle

を基底として計算します。ここで :math:`T_g` は各サイト :math:`i` の内容を
サイト :math:`g(i)` へ移す演算子で、 :math:`|r\rangle` は代表配置を走ります。
セクター内のすべての状態は :math:`T_g|\psi\rangle=\chi(g)|\psi\rangle`
を満たし、ヒルベルト空間の次元はおよそ群の位数分の1に縮小されます。
セクターの次元はログに ``Symmetry basis: raw_dim=... sector_dim=...``
として出力されます。

群は並進に限りません。ハミルトニアンがすべての操作で不変で、表現が
1次元である限り（可換群、またはより大きな群の1次元既約表現）、鏡映、
回転、およびそれらと並進の積も同じ書式で指定できます。2次元既約表現が
必要な場合は、可換部分群に制限してください。

スタンダードモードで ``MomentumIndex`` を指定すると、このファイルが
``qptransidx.def`` として生成されます
（:doc:`MomentumIndex パラメータ <../standardmode_ja/Parameters_for_conserved_quantities_ja>` を参照）。
長さ :math:`L` の chain では :math:`L` 個の並進と指標
:math:`\chi(T^{g})=\exp(-2\pi i m g/L)` が書かれ、 ``MomentumIndex`` の節では
これを運動量 :math:`k=2\pi m/L` と呼んでいます。

以下にファイル例を記載します（6サイトのリングの並進、
``MomentumIndex = 1``）。

::

    # MomentumIndex 1
    =============================================
    NQPTrans          6
    =============================================
    ======== TrIdx_TrWeight_and_TrIdx_i_xi ======
    =============================================
    0  1.000000000000000  0.000000000000000
    1  0.500000000000000 -0.866025403784439
    2 -0.500000000000000 -0.866025403784439
    3 -1.000000000000000  0.000000000000000
    4 -0.500000000000000  0.866025403784439
    5  0.500000000000000  0.866025403784439
    0 0 0 1
    0 1 1 1
    ...
    1 0 1 1
    1 1 2 1
    ...
    5 5 4 1

次の例は同じリングの鏡映 :math:`i\to 5-i` で、指標 :math:`-1`
（奇パリティ）のセクターです。

::

    =============================================
    NQPTrans          2
    =============================================
    ======== TrIdx_TrWeight_and_TrIdx_i_xi ======
    =============================================
    0  1.0
    1 -1.0
    0 0 0 1
    0 1 1 1
    0 2 2 1
    0 3 3 1
    0 4 4 1
    0 5 5 1
    1 0 5 1
    1 1 4 1
    1 2 3 1
    1 3 2 1
    1 4 1 1
    1 5 0 1

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  続く [int01] 行: [int02] [double01] [double02]

-  続く [int01] :math:`\times` ``Nsite`` 行: [int02] [int03] [int04] [int05]

``#`` で始まる行と空行はどこにあっても読み飛ばされ、上の行数には
数えません。

パラメータ
^^^^^^^^^^

-  :math:`[`\ string01\ :math:`]`

   **形式 :** string型 (空白不可)

   **説明 :** 対称操作の総数のキーワード名として ``NQPTrans`` を指定します。
   大文字小文字は区別しません。

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** 対称操作の総数、すなわち群の位数を指定します。

-  :math:`[`\ int02\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** 対称操作の番号を指定します。0以上\ :math:`[`\ int01\ :math:`]`\ 未満です。

-  :math:`[`\ double01\ :math:`]`, :math:`[`\ double02\ :math:`]`

   **形式 :** double型 (:math:`[`\ double02\ :math:`]` は省略可)

   **説明 :** 操作\ :math:`[`\ int02\ :math:`]`\ の指標 :math:`\chi(g)`
   の実部と虚部を指定します。\ :math:`[`\ double02\ :math:`]`\ を省略すると
   虚部は0になります。

-  :math:`[`\ int03\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** サイト番号を指定します。0以上\ ``Nsite``\ 未満で指定します。

-  :math:`[`\ int04\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** 操作\ :math:`[`\ int02\ :math:`]`\ がサイト\ :math:`[`\ int03\ :math:`]`\ を
   移す先のサイト番号を指定します。0以上\ ``Nsite``\ 未満で指定します。

-  :math:`[`\ int05\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** 反周期境界条件用に予約されています。本バージョンでは
   ``1`` を指定してください。

メタデータ
^^^^^^^^^^

``# MomentumIndex`` :math:`[`\ int06\ :math:`]` の形のコメント行は、
このファイルを生成した ``MomentumIndex`` を記録します。
:math:`[`\ int06\ :math:`]` は0以上でC言語の ``int`` 型に収まる整数
（``INT_MAX`` 以下）です。スタンダードモードは
これを1行目に書き、 :math:`{\mathcal H}\Phi` はログに
``TransSym metadata: MomentumIndex=``\ :math:`[`\ int06\ :math:`]`
と出力します。本バージョンでは計算には使われません。
その他のコメント行は無視されます。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  ヘッダの省略はできません。

-  各操作はサイトの全単射でなければならず、\ :math:`[`\ int02\ :math:`]`\ と
   :math:`[`\ int03\ :math:`]`\ の各組はちょうど1回ずつ指定します。

-  操作の集合は群でなければなりません。恒等操作を含み、任意の2つの
   操作の合成がいずれかの操作に一致する必要があります。指標は絶対値1で
   乗法的（:math:`\chi(gh)=\chi(g)\chi(h)`）、恒等操作の指標は
   :math:`\chi(e)=1` でなければなりません。

-  ハミルトニアンはすべての操作で不変でなければなりません。各項を操作で
   写した結果が元の項と一致しない場合、
   ``TransSym Hamiltonian invariance failed`` のエラーで終了します。

-  ``SpinlessFermion``、 ``Hubbard``、 ``tJ`` では、フェルミオンの置換に伴う符号は
   自動的に考慮されます。例えば占有軌道を入れ替える鏡映は :math:`-1` の
   因子を与えるため、偶・奇セクターの次元はスピン系の数え方とは異なります。

-  対応手法は ``Lanczos``、 ``CG``、 ``TPQ`` （microcanonical TPQ）、 ``cTPQ``、 ``FullDiag``、 ``TimeEvolution`` です。模型は :math:`S=1/2` で ``2Sz`` を
   固定した ``Spin``、 ``Ncond`` を固定した ``SpinlessFermion``、 ``Nup`` と
   ``Ndown`` を固定した ``Hubbard`` / ``tJ`` です。固定Szを持たないスピン1/2の
   ``SpinGC`` にも対応します（下記参照）。expert mode では次の項に対応します。

   - ``Spin``: 縦磁場の ``Trans``、 ``Exchange``、 ``Ising``、 ``CoulombInter``、
     ``Hund``、固定Szを保存する ``InterAll``。
   - ``SpinlessFermion``: サイト内ポテンシャルを含む ``Trans``、
     ``CoulombInter``、 ``InterAll``。
   - ``Hubbard``: スピンを保存する ``Trans``、 ``CoulombIntra``、 ``CoulombInter``、
     ``Hund``、 ``Ising``、 ``Exchange``、 ``PairHop``、固定スピンを保存する ``InterAll``。
   - ``tJ``: ``Hubbard`` と同じ項を二重占有のない配置へ射影します。
     ``CoulombIntra`` と ``PairHop`` の寄与は0になります。raw次元は
     :math:`\binom{N_{\rm site}}{N_\uparrow}\binom{N_{\rm site}-N_\uparrow}{N_\downarrow}` です。
     replicated / distributedの両layoutに対応し、rawのサイト分割で使えない
     MPIプロセス数でも実行できます。

   拡張項ではフェルミオンの正規順序化、または局所スピン行列の積の簡約後に
   係数を集約します。置換符号、縮約、重複項、family間の相殺を含めて、
   不変性と固定量子数の保存を検査します。係数の許容誤差は :math:`10^{-10}` です。
   これらcanonical模型では ``PairLift``、 ``NBodyInterAll``、異常項は非対応です。
   SpinGCは下記のとおり ``PairLift`` に対応します。
   spinlessのraw solverでは非対角 ``InterAll`` は引き続き非対応で、今回の拡張は
   ``TransSym`` に適用されます。Standard modeの入力生成は変更していません。

   相関関数（``OneBodyG``、``TwoBodyG``、``ThreeBodyG``、``FourBodyG``、
   ``SixBodyG``、``NBodyG``）は ``Lanczos``、``CG``、``TPQ``、``cTPQ``、
   ``TimeEvolution`` でsector内計算されます。
   値はsector内の厳密な期待値で、対称操作で互いに移り合う行は同じ値になり、
   sectorの外へ出る演算子は0になります。canonical ``Spin`` の
   ``ThreeBodyG``/``FourBodyG``/``SixBodyG`` はraw基底と同様に起動時に拒否されます。
   同じ積は ``NBodyG`` で指定できます。``AnomalousG``、スペクトル計算、リスタート、
   ハミルトニアンの入出力は本ファイルと併用できません。固有ベクトルの入出力は、
   CGとTimeEvolutionで下記のsector checkpoint形式により利用できます。
   未対応の組み合わせはエラーで終了します。


SpinGCセクター（スピン1/2、expert mode）
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``CalcModel=4`` と ``TransSym`` の組み合わせは、射影前の次元が
:math:`2^{N_{\rm site}}` の全スピン1/2空間から空間対称性のセクターを選びます。
全 :math:`S_z` は固定しません。``2Sz``、``Nup``、``Ndown``、``Ncond`` は
省略してください。値0の明示指定も拒否します。サイト数は
``0 < Nsite < CHAR_BIT * sizeof(unsigned long)`` を満たす必要があります。
bit 0が下向き、bit 1が上向きで、site 0が最下位bitです。各軌道の最小bit列を
代表とし、規格化した射影状態でその係数を正実にします。
セクターsolverは全空間ベクトルを必要としません。

対応するHamiltonian項は、横磁場・複素磁場を含むサイト内 ``Trans``、
``Ising``、``Exchange``、``CoulombInter``、``Hund``、``PairLift`` と、
サイト内スピン行列単位の積で表す ``InterAll`` です。全Szを変える項も使えますが、
Hermitian性と群の各操作に対する不変性は必要です。
:math:`E_i^{ab}=|a\rangle_i\langle b|` とすると、実数係数の ``PairLift`` の1行
``i j J`` は次の項を表します。

.. math::

   J(E_i^{10}E_j^{10}+E_i^{01}E_j^{01}).

追加の1/2はありません。逆向きの行と重複行は加算され、同一サイトの行は0です。
canonical ``Spin`` の固定Sz条件は変更しません。

.. list-table:: SpinGCセクターの手法と機能
   :header-rows: 1
   :widths: 17 24 22 19 18

   * - 手法 (CalcType)
     - 結果
     - 相関関数
     - ベクトル読込/出力
     - MPI layout
   * - Lanczos (0)
     - 低エネルギー状態
     - 全6形式
     - 不可/不可
     - distributed / replicated
   * - mTPQ (1)
     - 単一sectorの標本
     - 全6形式
     - 不可/不可
     - distributed / replicated
   * - FullDiag (2)
     - sectorの全固有値
     - 不可
     - 不可/不可
     - metadataはreplicated
   * - CG / LOBCG (3)
     - 低エネルギー状態
     - 全6形式
     - 可/可
     - distributed / replicated
   * - TimeEvolution (4)
     - 静的・動的時間発展
     - 全6形式
     - 必須/任意
     - distributed / replicated
   * - cTPQ (5)
     - 単一sectorの標本
     - 全6形式
     - 不可/不可
     - distributed / replicated

相関の6形式は ``OneBodyG``、``TwoBodyG``、``ThreeBodyG``、``FourBodyG``、
``SixBodyG``、``NBodyG`` です。``OneBodyG`` はサイト内演算子のみ指定でき、
サイト間の行は拒否します。aggregate形式と従来形式の両方に対応します。
FullDiagは下記の制限のもとでLAPACK（Solver 0、1rank）、ScaLAPACK（1）、
ELPA（3）を使えます。

:math:`S_z=\sum_i S_i^z` と :math:`S_z^2` は固定量子数から代入せず、実際の
sectorベクトルで評価します。既存の出力列は変更しません。CGのenergyファイルは
``Sz`` を出力し、TPQ/cTPQ/TEのFlctファイルは磁化の1次・2次モーメントを出力します。
CGに ``Sz2`` 列は追加しません。mTPQとcTPQは **単一の空間対称性sector** の
熱標本であり、全SpinGC空間の熱平均ではありません。sector間の和には適切な統計的
重みが必要で、単一sectorの標本はその和を表しません。

CGのベクトル読込は入力状態を評価し、最適化を再開しません。TEの読込は新しい時間列を
開始し、同一sector内のHamiltonian変更（quench）も許します。rank数、所有範囲、模型、
sectorと位相規約は一致させます。どちらもsolver restartではなく ``ReStart`` は拒否します。
SpinGC checkpointもversion 1で、headerは ``model=4``、``nup/ndown/ne`` は0です。
Hamiltonian fingerprintは解析済みPairLift行を含む
``hphi-parsed-hamiltonian-fnv1a64-v3`` です。canonical模型はv2と既存checkpoint形式を
維持します。manifestは ``fixed_quantities=none``、``full_dim=2^Nsite`` を記録します。

SpinGCは ``TEOneBody`` / ``TETwoBody`` の駆動に対応し、下記の右端Taylor規約を使います。
使用する全時刻の入力を伝播前に検査します。``Laser`` は拒否されます。
代わりに、群不変なサイト内 ``TEOneBody`` またはスピン積の ``TETwoBody`` を指定してください。
general spin、Boost、Kondo、新しいStandard-modeのSpinGC運動量入力、スピン軸回転、
全スピン反転、反ユニタリ操作、多次元既約表現、spectrum、Hamiltonian入出力、
solver restart、rank数を変更したcheckpoint再分配は対象外です。
このSpinGC sector経路では ``CoulombIntra``、``PairHop``、``NBodyInterAll``、
``AnomalousG`` も非対応です。FullDiagは固有ベクトル・相関出力、distributed基底metadata、
MAGMAと非serial ``ExpecMode`` も拒否します。

8サイト横磁場のCG例
^^^^^^^^^^^^^^^^^^^

空のディレクトリで次の標準ライブラリだけを使うPythonコードを ``make_input.py`` として
保存し、``python3 make_input.py`` を実行します。
:math:`H=-\sum_{i=0}^7 S_i^x`、characterが1の並進（運動量0）、CGのportableな
expert入力を生成します。正の ``Trans`` 係数はHPhiのtransfer規約の負号に対応します。

.. code-block:: python

   from pathlib import Path

   def definition(name, rows, count=None, keyword="NData"):
       rows = list(rows)
       header = "====\n{} {}\n====\n====\n====\n".format(
           keyword, len(rows) if count is None else count)
       Path(name).write_text(header + "".join(
           " ".join(map(str, row)) + "\n" for row in rows))

   Path("sym.def").write_text(
       "CalcMod calc.def\nModPara mod.def\nLocSpin loc.def\n"
       "TransSym group.def\nTrans trans.def\n")
   Path("calc.def").write_text(
       "CalcType 3\nCalcModel 4\nOutputMode 0\nOutputDataHead 1\n")
   Path("mod.def").write_text(
       "====\nModel_Parameters 0\n====\n====\n====\n"
       "CDataFileHead zvo\nCParaFileHead zqp\n====\n"
       "Nsite 8\nLanczos_max 400\ninitial_iv -1\nexct 1\n"
       "LanczosEps 18\nLanczosTarget 1\nLargeValue 100\nPreCG 0\n")
   definition("loc.def", ((i, 1) for i in range(8)))
   definition("trans.def", ((i, a, i, b, 0.5, 0)
              for i in range(8) for a, b in [(1, 0), (0, 1)]))
   definition("group.def", [(g, 1.0, 0.0) for g in range(8)] +
              [(g, i, (i+g) % 8, 1) for g in range(8) for i in range(8)],
              count=8, keyword="NQPTrans")

``OMP_NUM_THREADS=1`` として ``HPhi -e sym.def``、または
``mpiexec -np 4 HPhi -e sym.def`` を実行します。収束した基底状態は
:math:`E=-4`、:math:`\langle S_z\rangle=0`、:math:`\langle S_z^2\rangle=2` です。
``output/zvo_energy.dat`` にenergyとSzが出力されます。このCG例でSz2を求めるには、
``TwoBodyG`` で :math:`S_i^z S_j^z` を指定して全 :math:`i,j` について和を取ります
（:math:`S_i^z=(E_i^{11}-E_i^{00})/2`）。新しいStandard-mode keywordは使いません。

セクター内Lanczosの基底layout
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``TransSym`` を指定したセクター内Lanczosは、対称性で縮約した基底上で実行します。
既定の基底layoutは ``distributed`` で、各MPI rankは担当する行と、その行に必要な
ghost要素だけを保持します。行列ベクトル積には既存のdistributed planとhalo exchangeを
使用し、ベクトルnorm、内積、energy、varianceはMPI reductionで集約します。
MPI rank数がセクター次元より多い場合は担当行が0のrankもありますが、そのrankも
collective operationに参加します。

開発者が比較・診断する場合は ``HPHI_SYMMETRY_BASIS_LAYOUT=replicated`` を指定して
rollback用layoutを選択できます。この設定は通常運用の推奨値ではありません。
layoutの選択によって出力ファイル名、出力形式、セクターmanifestのschemaは変わりません。

セクター内TPQ
^^^^^^^^^^^^^

expert modeで ``CalcType=1`` と ``TransSym`` を指定すると、選択したセクター内の
ランダムベクトルへ :math:`l-H_q/N_{\rm site}` を反復して作用させます。
既定の基底layoutはdistributedで、 ``HPHI_SYMMETRY_BASIS_LAYOUT=replicated`` により
参照用layoutを選択できます。要素を持たないMPI rankもglobal normの集約に参加します。
``Lanczos_max``、 ``NumAve``、 ``LargeValue``、 ``initial_iv``、 ``InitialVecType`` は
通常のTPQと同じ意味です。 ``exct`` によるセクター次元の制限はありません。
Standard modeの ``MomentumIndex`` 生成は引き続きLanczosとCGに限られます。

Hubbardのdoublonとその二乗平均を含む、既存SS/Norm/Flctの全列に対応します。
SSの ``phys_var`` 列は従来どおり :math:`\langle H^2\rangle` で、
エネルギー平均の二乗を引いた分散ではありません。
``OutputGreenFormat=1`` によりSS/Norm/Flctと相関関数の集約出力を使えます。
相関関数は ``ExpecInterval`` が選ぶTPQ評価stepで再計算されます。大規模sectorで
全site pairを毎step要求すると高コストになります。restart、ベクトル入出力、
スペクトルは ``TransSym`` と併用できません。

計算結果は **単一の対称性セクター内** のtraceを推定します。
固定量子数空間全体の熱平均ではありません。manifestには
``ensemble=single_symmetry_sector``、 ``num_ave``、 ``large_value``、
``initial_vec_type`` を追記します。既存の乱数生成はMPIの要素分担とOpenMP threadに
依存するため、seedを固定してもprocess/thread数が異なる場合の標本は一致しません。

セクター内cTPQ
^^^^^^^^^^^^^^

``CalcType=5`` はmTPQと同じセクター保存方式・出力を使い、
:math:`\sum_{n=0}^{n_{\max}}(-\Delta\beta H_q/2)^n/n!` を作用させて各stepで正規化します。
Taylor打ち切り次数の収束は利用者が確認してください。
``InvTemp`` を指定しない場合は :math:`\Delta\beta=1/\mathrm{LargeValue}` で、
``ExpandCoef`` の省略時既定値は10です。明示指定した次数は正の整数で
``INT_MAX`` 未満に限ります。normがゼロ・非有限の場合とmatvec失敗時は停止します。

``InvTemp`` の各行は ``beta nmax physcal eigen`` です。betaは0から始まり、
有限値で単調非減少である必要があります。同じbetaの繰り返しは刻み0として扱います。 ``nmax`` は ``INT_MAX`` 未満の正の整数、
両flagは0または1です。行 ``i`` の次数はbeta ``i`` から ``i+1`` への発展に使い、
最終行の次数は使いません。セクター内cTPQでは初期・最終点も含めて
``eigen`` が非ゼロの行を拒否します。ベクトル入出力とrestartは引き続き非対応です。
明示scheduleでは初期点と ``physcal=1`` の行で相関関数を出力します。``InvTemp`` を
指定しない場合は ``ExpecInterval`` が出力stepを選びます。``OutputGreenFormat=1`` は
全sample・stepを添字付きの相関関数集約ファイルへまとめます。
manifestには ``canonical_tpq_steps``、 ``beta_schedule`` と、一定刻み・次数または
全beta・次数行を記録します。mTPQと同様に、全セクターの和ではなく単一セクターの結果です。

静的ハミルトニアンでのセクター時間発展
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

expert modeの ``CalcType=4`` は、保存済みセクターベクトルを静的Hamiltonianで
時間発展させます。 ``InputEigenVec=1``、 ``ReStart=0``、正整数の ``ExpandCoef`` を
指定し、namelistへ ``SpectrumVec seed_eigenvec_0`` を追加します。
``output/seed_eigenvec_0_rank_<rank>.dat`` は下記のsector形式が必要です。
CGの ``OutputEigenVec=1`` で作成できます。基底・sector・MPI数・layoutは一致が
必要ですが、quenchのためHamiltonianは変更できます。既定layoutはdistributed、
``HPHI_SYMMETRY_BASIS_LAYOUT=replicated`` で参照用layoutを選択できます。

TE実行前にCGのseedを全rank分、別prefixで保存するか、TE側の
``CDataFileHead`` を変更してください。入出力のprefixが同じ場合、row 0の出力
``<prefix>_eigenvec_0_rank_<rank>.dat`` がseedを上書きします。
sectorの最終出力が別名であっても、このrow 0での上書きは防げません。

全行の項数が0の ``TEOneBody`` で時刻列を指定し、 ``Lanczos_max`` 行を使用します。
時刻は有限で単調非減少とします。最初の行は入力状態を伝播せず記録し、以降は
:math:`\exp[-iH(t_j-t_{j-1})]` の ``ExpandCoef`` 次Taylor多項式を作用させて
規格化します。同一時刻も許容します。入力checkpointのstep/timeは来歴として
記録し、新しい時刻列の時計には引き継ぎません。時間依存項とPeierls駆動については次節を参照してください。

``SS`` / ``Norm`` / ``Flct`` はsector内の期待値を出力します。 ``Norm`` は各stepの
規格化前のnormで、Taylor打切りにより1からずれる場合があります。次数を増やすか
時間間隔を減らして収束を確認してください。相関関数は伝播後、``ExpecInterval`` が
選ぶstepで評価されます。``OutputGreenFormat=1`` は0始まりの時刻stepを添字として
集約ファイルへまとめます。大規模sectorで全site pairを毎step要求すると高コストに
なります。``ReStart`` は引き続き非対応です。
manifestには実際の時刻列・次数と入力checkpointのmethod/state/step/time/H digestを
記録します。

``OutputEigenVec=1`` と正の ``OutputInterval`` により、周期出力を
``<prefix>_eigenvec_<step>_rank_<rank>.dat``、最終出力を別名
``<prefix>_eigenvec_final_rank_<rank>.dat`` に保存します。sector headerの ``step`` は
完了した0始まりの時刻行、 ``time`` はその物理時刻、 ``state_index`` は入力状態の
ラベルです。最終出力はrow 0を上書きしません。これらのcheckpointを ``SpectrumVec``
で新しい同一sector計算へ渡せますが、restartではありません。raw TEの既存binary形式とは
異なります。

時間依存するセクターHamiltonian
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``TEOneBody`` または ``TETwoBody`` の1種類で駆動できます。
canonical模型では ``Laser`` も選べますが、SpinGCでは拒否します。
one-body/two-body項は静的Hamiltonianへの加算、Peierls駆動は解析済みtransfer係数の
位相変更です。上記4模型で対角・非対角項に対応し、spinless fermionも含みます。
raw spinless solverの対角TE項に対する制限は変更しません。

初期vector読込と時系列出力より前に、使用する全時刻について係数の有限性、
固定量子数の保存、指定群に対する不変性を検査します。後の時刻だけで対称性が
破れる場合も伝播開始前に停止します。 ``Laser`` は既存の9パラメータを使い、
時刻を ``Tinit + step * TimeSlice`` とします。時間刻みは有限・非負とし、最初の行は
``Tinit`` における入力状態を伝播せず記録します。 ``TEOneBody`` / ``TETwoBody`` は
ファイルの有限・単調非減少の時刻列を使います。

区間の終点 :math:`t_j` のHamiltonian :math:`H(t_j)` を固定し、そのTaylor多項式を
前の状態へ作用させます。一次項を含む全次数で同じ行列を使います。
連続的な時間依存Hamiltonianに対しては右端点での区分一定近似なので、Taylor次数だけを
増やしても時間離散化誤差は消えません。時間刻みと次数の両方で収束を確認してください。
このsectorの規約は、従来のraw動的TEの漸化式とは区別されます。

現実装は基底の代表状態・順序・位相・MPI分担を保持し、各時刻で対角成分を更新します。
既存のrank局所基底からlookup directoryを再作成し、行列作用のplanを再構築します。
raw Hilbert空間を列挙し直す必要がなく、項の出現・消滅にも対応します。
共通疎構造のcacheやplanの係数だけの更新は、まだ行いません。manifestには ``te_hamiltonian=time_dependent``、
``te_integrator=right_endpoint_taylor``、 ``te_plan_update=rebuild_with_fixed_basis`` と
各時刻の実効Hamiltonian digest ``te_hamiltonian_<step>`` を記録します。
同じdigestをその行のcheckpoint headerにも記録します。

セクターベクトルcheckpoint
^^^^^^^^^^^^^^^^^^^^^^^^^^

expert modeのCGでは ``OutputEigenVec=1`` と ``InputEigenVec=1`` を使えます。
ファイル名は ``output/<CDataFileHead>_eigenvec_<state>_rank_<rank>.dat`` で、
状態・rankの番号は0始まりです。要素を持たないrankもファイルを出力します。
読込には、模型・固定量子数・群と指標・sector・MPI rank数・layout・局所分担・
global indexと代表状態の対応・基底位相規約の一致が必要です。
従来のraw基底vector、 ``InputEigenVec=2``、 ``ReStart`` は拒否します。
CGの入力実行は与えた状態の物理量を評価し、固有状態の再計算や反復再開は行いません。

生成時のHamiltonian digestは別に記録し、同一sector内のquench用に現在のHamiltonianと
異なる場合も読み込めます。変更の有無をログへ表示します。この状態読込はsolverのrestartとは別です。
全rankでmetadataを照合してからvector本体を検査し、成功した場合だけ計算用vectorへ反映します。
切れたファイル、余分なデータ、非有限値、1から :math:`10^{-8}` を超えてずれたglobal normの二乗、
checksum不一致は拒否します。異なるcheckpoint setのrankファイルも混在できません。

version 1形式はendiannessに依存しません。28個のlittle-endian unsigned 64-bit wordのheaderに、
``local_dim`` 個のIEEE binary64実部・虚部の組が続きます。添字0の未使用成分は保存しません。
headerの順序は次のとおりです。

.. code-block:: text

   magic version phase scalar model nsite nup ndown ne
   raw_dim sector_dim ranks rank layout offset local_dim group_digest
   sector_count sector_xor sector_sum order_digest hamiltonian_digest
   source_method state_index step time_bits payload_xor payload_sum

``magic`` は8 byteの ``HPHISV1\n``、versionとphaseは1、scalarは128、layoutは
0（replicated）または1（distributed）です。 ``offset`` は0始まりです。
``source_method`` はCalcType番号、 ``time_bits`` はIEEE binary64の物理時刻（CGでは0）です。
phase 1は最小代表状態に対する正規化済み
:math:`\sum_g\overline{\chi(g)}T_g|r\rangle` で、代表状態の係数を正の実数に固定します。

群・sector・Hamiltonianのfingerprintはmanifestと同じalgorithmです。order digestは所有成分順に
1始まりのglobal index・代表整数・orbitサイズ・stabilizerサイズを各8 byteのlittle-endianで
FNV-1a-64へ入力します。payload digestはrank・offset・局所次元（各8 byte）に続けて
保存した実部・虚部を入力し、そのrank別hashのXORと :math:`2^{64}` を法とする和を記録します。
これらは整合性検査用fingerprintで、暗号学的hashではありません。
全rankの ``.part`` 一時ファイル書込成功後に正式名へ変更します。正式名への変更失敗も
異常終了とし、readerは全rankが同じcheckpoint setに属することを検査します。

セクター内FullDiag
^^^^^^^^^^^^^^^^^^

expert modeの ``CalcType=2`` で、選択したセクターの全固有値を計算します。
buildに含まれていれば ``Solver=0`` （LAPACK、MPI 1 rank）、 ``Solver=1``
（ScaLAPACK）、 ``Solver=3`` （ELPA）に対応します。
基底metadataはreplicatedに限られ、FullDiagでの
``HPHI_SYMMETRY_BASIS_LAYOUT=distributed`` は拒否します。
solverの作業配列には全rankでセクター全次元を使います。
複数rankのELPAでは所有するcolumn panelだけを生成します。
ScaLAPACKは現段階ではHamiltonianを複製しますが、使わない複製固有ベクトル行列は
確保しません。ELPAではセクター次元がprocess gridより小さければエラーとなるため、
その場合はrank数を減らしてください。

``output/<CDataFileHead>_energy_sector.dat`` へ0始まりの番号と全固有値を出力し、
セクターmanifestも記録します。この固有値ファイル名は常にprefix付きです。
通常のFullDiagは従来どおり ``Eigenvalue.dat`` を出力します。
セクターFullDiagでは固有状態物理量、相関関数、固有ベクトルを出力しません。
0以外の ``ExpecMode`` とMAGMAは非対応です。
Hamiltonian/ベクトル入出力とrestartも非対応です。
manifestにはsolver、行列storage、 ``output_scope=eigenvalues`` を記録します。

セクター情報ファイル
^^^^^^^^^^^^^^^^^^^^

非空の対称化基底を構築し、セクターのオプションを検証した後、solver の開始前に
``output/symmetry_sector.dat`` を出力します。 ``OutputDataHead=1`` の場合は
``output/<CDataFileHead>_symmetry_sector.dat`` です。 ``TransSym`` を指定しない
通常の計算と、 ``-sdry`` による定義ファイル生成では出力しません。
書き込みエラーがあれば全 MPI rank で計算を終了します。

先頭行は ``format=HPhiSymmetrySector version=1`` です。以降は ``key=value`` の
形式で、計算手法、模型、サイト数、固定量子数、固定量子数空間の全次元
（``full_dim``）、セクター次元（``sector_dim``）、群の位数、任意指定の
``momentum_index``、基底 layout、MPI rank 数、OpenMP thread 数の上限、
項の件数、solver のパラメータを記録します。入力セクターを記録するファイルであり、
ファイルの存在は計算の完了・収束を意味しません。

次の3種類の、algorithm version を含む fingerprint を出力します。

- ``group_digest``: ``hphi-group-fnv1a64-v1`` は、操作をサイト置換の辞書順に
  並べて置換と指標を hash 化します。操作番号の変更には依存しません。
  指標の各成分は :math:`10^{-10}` 刻みに丸めます。この量子化は、丸め境界付近の
  浮動小数点入力を一般的に同値判定するものではありません。
- ``sector_digest``: ``hphi-sector-multiset-v1:count:xor:sum`` は、各基底の
  代表状態、軌道サイズ、固定部分群サイズを FNV-1a 64 で hash 化し、列挙順や
  MPI 分割によらず集約します。整数は固定幅 little-endian、和は :math:`2^{64}` を
  法として計算します。ノルムとハミルトニアンの対角値は含めません。
- ``hamiltonian_digest``: canonical模型の ``hphi-parsed-hamiltonian-fnv1a64-v2`` は、対応する
  ハミルトニアンの項を格納順に記録し、係数には binary64 の bit 列を使います。
  version 2ではサイト内ポテンシャル、pair hopping、分離済みの対角・非対角InterAllを
  記録対象に追加しています。物理的に同じハミルトニアンでも、項の順序や分解が違えば値が異なることがあります。

セクターの識別には、模型、固定量子数、セクター次元、 ``group_digest``、
``sector_digest`` の組を使います。共役表現では ``sector_digest`` が同じに
なり得るため、単独では使用できません。結合定数の変更はセクター識別に影響しません。
これらは照合のための fingerprint であり、衝突のない同値性の証明ではありません。
順序に依存しない fingerprint だけでは、checkpoint のベクトル成分の並びも検証できません。
従来の計算結果と ``CalcTimerRankStats.dat`` の形式は変更しません。

.. raw:: latex

   \newpage
