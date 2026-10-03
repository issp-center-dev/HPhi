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

-  対応手法は ``Lanczos``、 ``CG``、 ``TPQ`` （microcanonical TPQ）、 ``FullDiag`` です。模型は :math:`S=1/2` で ``2Sz`` を
   固定した ``Spin``、 ``Ncond`` を固定した ``SpinlessFermion``、 ``Nup`` と
   ``Ndown`` を固定した ``Hubbard`` / ``tJ`` です。expert mode では次の項に対応します。

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
   ``PairLift``、 ``NBodyInterAll``、異常項は非対応です。
   spinlessのraw solverでは非対角 ``InterAll`` は引き続き非対応で、今回の拡張は
   ``TransSym`` に適用されます。Standard modeの入力生成は変更していません。

   相関関数、スペクトル計算、リスタート、ハミルトニアンと固有ベクトルの入出力は
   本ファイルと併用できません。非対応の組み合わせはエラーで終了します。

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
``OutputGreenFormat=1`` によるSS/Norm/Flct集約出力も使えます。
相関関数、restart、ベクトル入出力、スペクトル、cTPQは ``TransSym`` と併用できません。

計算結果は **単一の対称性セクター内** のtraceを推定します。
固定量子数空間全体の熱平均ではありません。manifestには
``ensemble=single_symmetry_sector``、 ``num_ave``、 ``large_value``、
``initial_vec_type`` を追記します。既存の乱数生成はMPIの要素分担とOpenMP threadに
依存するため、seedを固定してもprocess/thread数が異なる場合の標本は一致しません。

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
- ``hamiltonian_digest``: ``hphi-parsed-hamiltonian-fnv1a64-v2`` は、対応する
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
