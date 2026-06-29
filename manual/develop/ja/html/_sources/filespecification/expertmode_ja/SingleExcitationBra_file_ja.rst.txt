.. highlight:: none

.. _Subsec:singleexcitationbra:

SingleExcitationBra指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

非対角動的グリーン関数

.. math::

   G_{BA}(z) = \langle \phi | B^{\dagger} \frac{1}{z-H} A | \phi \rangle

を計算するための\ **ブラ側**\ 一体励起演算子\ :math:`B=c_{i\sigma_1}(c_{i\sigma_1}^{\dagger})`\ を定義します。
ケット側演算子\ :math:`A`\ は\ :doc:`SingleExcitation <SingleExcitation_file_ja>`\ 指定ファイルで定義します。
本ファイルを省略した場合は、通常の対角グリーン関数（\ :math:`B=A`\ ）が計算されます。

ファイルフォーマットおよびパラメータは\ :doc:`SingleExcitation <SingleExcitation_file_ja>`\ 指定ファイルと同一で、
演算子数のキーワードのみ ``NSingleExcitationBra`` です。以下にファイル例を記載します。

::

    ===============================
    NSingleExcitationBra    1
    ===============================
    ====== Single Excitation Bra ===
    ===============================
        0     0     0    1.0    0.0

.. _use_rules_singleexcitationbra_ja:

使用ルール
^^^^^^^^^^

*  ヘッダの省略はできません。フォーマット・パラメータは\ :doc:`SingleExcitation <SingleExcitation_file_ja>`\ と同一です。

*  本ファイルに記載する演算子は\ :math:`B`\ であり\ :math:`B^{\dagger}`\ ではありません
   （計算される量は\ :math:`\langle\phi|B^{\dagger}(z-H)^{-1}A|\phi\rangle`\ です）。

*  非対角グリーン関数はシフト型 BiCG 法（\ ``method="CG"``\ ）かつ ``CalcSpec="Normal"`` のときのみ利用可能です
   （ブラ演算子を指定した場合、リスタート・保存系のスペクトルモードは利用できません）。

*  ケット（\ :doc:`SingleExcitation <SingleExcitation_file_ja>`\ ）とブラの演算子は、\ :math:`|\phi\rangle`\ を
   \ **同一**\ の励起ヒルベルト空間へ写す必要があり、1 つのファイル内の全演算子も同一のセクター変化をもつ必要があります
   （これを満たさない場合はプログラムが終了します）。「同一セクター」の意味は模型によって異なります。

   *  Hubbard 模型（カノニカル）：粒子数\ **および**\ :math:`S_z`\ の変化が同じであること。

   *  HubbardGC 模型（グランドカノニカル）：一体励起はセクター変化を伴わず、全 Fock 空間が
      1 つのヒルベルト空間となるため、ケットとブラは（どちらのスピンでも）常に同一空間を共有します。

   *  HubbardNConserved 模型（\ :math:`S_z`\ 自由）：粒子数の変化のみが同じであればよく、\ :math:`S_z`\ の変化は異なってもよい。

*  対応モデル：Hubbard 模型、HubbardGC 模型（グランドカノニカル）、および
   HubbardNConserved 模型（カノニカル、電子数固定で\ :math:`S_z`\ は自由）。
   （一体励起はスピン系では定義されず、その他の模型は非対角グリーン関数に未対応のため、
   これらの場合はプログラムが終了します。）

*  HubbardNConserved 模型は Hubbard 模型（\ ``CalcModel=0``\ ）において\ :ref:`ModPara <Subsec:modpara>`\ ファイルで
   ``2Sz``\ を指定せず\ ``Ncond``\ （電子数）のみを与えることで選択されます。このとき\ :math:`S_z`\ は自由となるため、
   スピン交差（2Sz 非保存）の非対角ルート（例えばスピン軌道結合をもつ不純物に用いる）が有効になります。励起ヒルベルト空間は電子数のみで決まり、
   ケットとブラは粒子数の変化が同じであればよく、スピン交差の\ :math:`G_{BA}`\ を 1 回の求解で得られます。

*  \ :ref:`ModPara <Subsec:modpara>`\ ファイルで\ ``SpectrumNumBra``\ :math:`>1`\ を
   指定した場合、追加の bra 演算子セット\ ``1``\ ...\ :math:`(B-1)`\ は実行
   ディレクトリの\ ``single_ex_bra_1.def``\ ...\ ``single_ex_bra_``\ :math:`(B-1)`\ ``.def``\ から
   本ファイルと同じ形式で読み込まれます(bra セット\ ``0``\ は namelist で指定する
   本\ ``SingleExcitationBra``\ ファイルです)。

.. raw:: latex

   \newpage
