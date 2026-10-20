.. highlight:: none

.. _Subsec:ham:

ham.dat
~~~~~~~

(FullDiagでのみ出力)
``CalcMod``\ ファイルで\ ``OutputHam=1``\ の場合に、HΦ内部で計算されたハミルトニアンをMatrixMarket形式で出力します。\ ``CalcMod``\ ファイルで\ ``InputHam=1``\ とすると、定義ファイル一式と本ファイルを読み込み、再計算することができます。以下にファイル例を記載します。

::

    %%%%MatrixMarket matrix coordinate complex hermitian
    28 28 56
    1 1 1.000000 0.000000
    2 1 0.500000 0.000000
    3 2 0.500000 0.000000
    4 3 0.500000 0.000000
    5 4 0.500000 0.000000
    6 5 0.500000 0.000000
    7 6 0.500000 0.000000
    7 7 1.000000 0.000000
        …

ファイル名
^^^^^^^^^^

-  ##\_ham.dat

##はModParaファイル内の[string02]で指定されるヘッダを表します。

ファイル形式
^^^^^^^^^^^^

-  1行：ヘッダ

-  2行：\ :math:`[`\ int01\ :math:`]`  :math:`[`\ int02\ :math:`]`  :math:`[`\ int03\ :math:`]`

-  3行-：\ :math:`[`\ int04\ :math:`]`  :math:`[`\ int05\ :math:`]`  :math:`[`\ double01\ :math:`]`  :math:`[`\ double02\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :** Hamiltonianの行数。

-  :math:`[`\ int02\ :math:`]`

   **形式 :** int型

   **説明 :** Hamiltonianの列数。

-  :math:`[`\ int03\ :math:`]`

   **形式 :** int型

   **説明 :** Hamiltonianの非零の要素数。

-  :math:`[`\ double01\ :math:`]`, :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   | **説明 :** Hamiltonianの値を表します。
   | :math:`[`\ double01\ :math:`]`\ が実部、\ :math:`[`\ double02\ :math:`]`\ が虚部を表します。

sector FullDiagのハミルトニアン
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

エキスパートモードで ``TransSym`` （:ref:`Subsec:TransSym`）を指定した
``CalcType=2`` では、 ``OutputHam=1`` が選択したsectorの :math:`D\times D`
ハミルトニアンを ``output/##_ham_sector.mtx`` へ出力して終了します。
固有値・相関・固有状態は計算しません。 ``InputHam=1`` は同じファイルを読み込み、
その行列を対角化します。このときも固有値に加えて、指定した相関と
``OutputEigenVec=1`` の全固有状態を出力できます。 ``OutputHam`` と
``InputHam`` の同時指定はできません。

ファイルはMatrix Marketの ``coordinate complex hermitian`` 形式です。
基底はHPhiのsector代表元の順序と位相規約で、 ``% HPHI key=value`` の行が
その基底を識別します。次はTransSymを指定したスピン鎖（6サイト、運動量 :math:`k=1` 、 :math:`D=3` ）の出力例です。

::

    %%MatrixMarket matrix coordinate complex hermitian
    % HPHI format=sector-hamiltonian-v1
    % HPHI scalar=complex-binary64
    % HPHI model=Spin
    % HPHI nsite=6
    % HPHI nup=3
    % HPHI ndown=3
    % HPHI ne=3
    % HPHI local_site_mask=0000000000000000
    % HPHI fixed_flags=0
    % HPHI phase=1
    % HPHI full_dim=20
    % HPHI sector_dim=3
    % HPHI group_digest=8a08ff3f5c2433e3
    % HPHI sector_count=3
    % HPHI sector_xor=ba8fc1732fc70785
    % HPHI sector_sum=640126bcd43ffa9b
    % HPHI basis_order=110eb94c9a548c63
    3 3 6
    1 1 0.80000000000000004 0
    2 1 0.62 -0.11
    3 1 0.21473720558371184 0.59193575034635193
    2 2 0.48999999999999999 0
    3 2 0.40526279441628799 0.48193575034635183
    3 3 0.48999999999999999 0

-  1行目は ``%%MatrixMarket matrix coordinate complex hermitian`` です。

-  ``% HPHI`` の行は次の17個のキーを1回ずつ持ちます（順序は任意）。
   整数は10進、maskとdigestは16桁の16進です。
   ``format=sector-hamiltonian-v1``、 ``scalar=complex-binary64``、
   ``model`` （ ``Spin``、 ``SpinGC``、 ``SpinlessFermion``、 ``Hubbard``、
   ``tJ``、 ``Kondo``、 ``KondoNConserved``、 ``KondoGC``、 ``HubbardNConserved``、
   ``tJNConserved``、 ``HubbardGC``、 ``SpinlessFermionGC`` ）、
   ``nsite``、 ``nup``、 ``ndown``、 ``ne``、 ``local_site_mask``、
   ``fixed_flags``、 ``phase``、 ``full_dim``、 ``sector_dim``、
   ``group_digest``、 ``sector_count``、 ``sector_xor``、 ``sector_sum``、
   ``basis_order`` 。値はsector checkpointの基底識別と同じです。群とsectorの
   digestなど一部は :ref:`Subsec:TransSym` のmanifestにも記録されます。

-  次元行は ``D D nnz`` 、データ行は ``row column real imag`` です。
   1始まりの下三角（ :math:`{\rm row}\ge{\rm column}` ）だけを書きます。
   HPhiは列の昇順・列内の行の昇順に、0でない要素だけを有効数字17桁で書きます
   （符号付きの0は+0とします）。対角化した行列と同じbinary64値に戻ります。

読込では、行の順序は任意で、明示的な0も書けます。書かれていない要素は0、
上三角は下三角の複素共役です。対角の虚部は厳密に0でなければなりません。
1行目と次元行の間および最後のデータ行の後には、空行と一般の ``%`` コメントを
置けます。改行はLFとCRLFのどちらでも読めます。

識別キーが現在の入力のsectorと一致しない場合（同じ次元の別sectorを含む）、
キーの欠落・重複・未知のキー、値の書式違い、次元・ ``nnz`` の不一致、範囲外や
上三角の添字、座標の重複（0を含む）、非有限値、0でない対角虚部、宣言より
少ない・多いデータ行、1行内の余分なtoken、長すぎる行は、エラーとして
対角化の前に全MPIランクで終了します。

行列の係数を編集しても ``% HPHI`` の行は残してください。基底識別が一致すれば、
外部で作成した行列も同じ規則で読み込めます。定義ファイルのハミルトニアンと
係数が異なることはエラーにしません。例えば、全係数を2倍して固有値を2倍にする
POSIX shellの例は次のとおりです。

::

    awk '/^%/ || NF != 4 { print; next }
         { $3 = sprintf("%.17g", 2 * $3); $4 = sprintf("%.17g", 2 * $4); print }' \
        output/zvo_ham_sector.mtx > edited.mtx
    mv edited.mtx output/zvo_ham_sector.mtx
    # calc.def の OutputHam を 0、InputHam を 1 にして同じ入力で実行します。

同じ基底の行列ファイルは、MPIランク数とsolver（LAPACK、ScaLAPACK、ELPA）を
変えて読み込めます。複数ランクのELPAでも、各ランクは自身の担当列だけを保持し、
ルートが読んだデータ行を一定件数ごとに配ります。一方、固有状態のcheckpointは
保存時と同じランク数とlayoutが必要です。

``InputHam=1`` の実行では、manifestに ``matrix_source=input`` と、読み込んだ値の
``matrix_digest`` を記録し、固有状態のcheckpointもこのdigestを持ちます。
保存した状態をTEやCGへ渡すと、利用側は自身の定義ファイルのハミルトニアンで
時間発展・再評価します（ :ref:`Subsec:eigenvec` ）。編集した行列の固有値とは
一致しないことがあります。

.. raw:: latex

   \newpage
