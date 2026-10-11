.. highlight:: none

.. _Subsec:eigenvec:

eigenvec.dat
~~~~~~~~~~~~

``CalcMod`` ファイルで ``OutputEigenVec=1`` とした場合に、固有ベクトルを
バイナリ形式で出力します。反復法では指定された状態を出力し、FullDiagでは
全固有ベクトルを出力します。

``InputEigenVec=1`` は、対応するリスタート・スペクトル計算で反復法の
出力ファイルを読み込みます。FullDiag出力も同じバイナリ配置ですが、
FullDiagのリスタート機能ではありません。FullDiagでは常に対角化を実行します。

ファイル名
^^^^^^^^^^

-  ##\_eigenvec\_&&\_rank\_$$.dat

##はModParaファイル内の[string02]で指定されるヘッダ、&&は0始まりの
固有状態番号、$$はランク番号を表します。

反復法では、各MPIランクが局所Hilbert空間の成分を自身のランク番号の
ファイルへ出力します。FullDiagでは各ファイルに分割されていない1本の
固有ベクトルを格納するため、ScaLAPACK/ELPAや ``ExpecMode=1`` / ``2`` を
用いたMPI実行でも ``##_eigenvec_&&_rank_0.dat`` だけを出力します。

ファイル形式
^^^^^^^^^^^^

ファイルには、順に反復回数（``int``）、ベクトル次元
（``unsigned long int``）、``complex double`` 型の
``local_size+1`` 成分を格納します。

**Note:**  ``eigen_vector``\ の一番最初の成分は計算に使用しません。

FullDiagでは反復回数は0、ベクトル次元はHilbert空間の全次元です。この指定を
有効にすると、:math:`N_{\rm H}` 本のファイルへ合計
:math:`N_{\rm H}^2` 個の複素数成分を出力するため、ディスク使用量は
:math:`O(N_{\rm H}^2)` になります。

sector FullDiagの固有状態
^^^^^^^^^^^^^^^^^^^^^^^^^

エキスパートモードで ``TransSym`` を指定し、``CalcType=2`` と
``OutputEigenVec=1`` を与えると、sector内の **全** 固有状態を
``##_eigenvec_&&_rank_$$.dat`` に保存します。状態番号&&はエネルギー昇順の
0始まりで、``##_energy_sector.dat`` と同じ順序です。縮退した準位内の状態は
一意ではなく、順序・位相・基底はsolverや実行ごとに異なり得ます。
各MPIランクは所有成分を持たない場合も含めて自身の区間を出力するため、
次元 :math:`D` のsectorを :math:`P` ランクで計算すると
:math:`D\times P` 個のファイルになります
（:math:`D=10{,}000`、:math:`P=16` なら160,000個）。状態範囲の指定は
できません。

これらのファイルは上記のrawバイナリ配置ではなく、:ref:`Subsec:TransSym` に
記載したsector checkpoint形式です。raw FullDiagと意図的に同じファイル名を
使っており、``TransSym`` の有無で形式が決まります。読み込む前に、どちらの
計算の出力かを確認してください。ヘッダには ``source_method=2``
（FullDiag）、状態番号、step 0、time 0と、定義ファイルのHamiltonian
digestではなく対角化した行列の ``hphi-sector-matrix-fnv1a64-v1`` digestを
記録します。

保存した状態は、同じsector・MPIランク数・基底順序・位相・layoutの計算で
再利用できます。

* 時間発展: 例えば状態2のランク別ファイルを、別の実行ディレクトリの
  ``output/seed_rank_$$.dat`` へコピーし、namelistに ``SpectrumVec seed``
  を加えます（``CalcType=4``、``InputEigenVec=1``）。TEは自身の
  ``##_eigenvec_&&`` ファイルを出力するため、同じディレクトリでは
  他のFullDiag状態を上書きします。
* CG: ``InputEigenVec=1`` は ``##_eigenvec_0`` から
  ``##_eigenvec_(exct-1)`` を読み込み、反復せずにCG入力のHamiltonianで
  再評価します。再評価だけなら ``OutputEigenVec=0`` を指定してください。
  ``OutputEigenVec=1`` は入力ファイルをCGのcheckpointで上書きし、CGの
  相関ファイルもFullDiagと同じ名前（ ``##_cisajs_eigen&&.dat`` など）で
  出力します。元のFullDiagの出力を残すには別ディレクトリへコピーしてください。

TEとCGの既定layoutは ``HPHI_SYMMETRY_BASIS_LAYOUT=distributed`` ですが、
FullDiagはreplicated layoutで状態を書きます。利用側の実行では
``HPHI_SYMMETRY_BASIS_LAYOUT=replicated`` を指定してください。利用側は
自身の定義ファイルのHamiltonianを評価し、quenchも可能です。FullDiagの
digestは行列digestなので、読込ログは
``source_hamiltonian_kind=sector_matrix_v1 hamiltonian_changed=unknown``
を表示します。

``fulldiag/`` で2ランクのFullDiagを実行した後、状態2からTEを行う
POSIX shellの例:

::

    mkdir -p te/output
    cp fulldiag/*.def te/
    for rank in 0 1; do
      cp fulldiag/output/zvo_eigenvec_2_rank_${rank}.dat \
         te/output/seed_rank_${rank}.dat
    done
    # te/ で calc.def に CalcType 4 と InputEigenVec 1 を設定し、
    # TEOneBody の時刻表と "SpectrumVec seed" を namelist に加えてから、
    # 同じ2ランクで実行します:
    # HPHI_SYMMETRY_BASIS_LAYOUT=replicated mpiexec -np 2 HPhi -e namelist.def

.. raw:: latex

   \newpage
