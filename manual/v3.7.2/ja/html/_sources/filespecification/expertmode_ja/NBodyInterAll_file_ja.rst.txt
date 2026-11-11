.. highlight:: none

.. _Subsec:nbodyinterall:

NBodyInterAll指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~

N体相互作用をハミルトニアンに付け加えます。1行で :math:`N` 個の因子を指定し、
以下の項を表します。

.. math::

   \mathcal{H} += V \prod_{p=1}^{N}
   c_{i_p\sigma'_p}^{\dagger} c_{j_p\sigma_p}.

Spin/SpinGC の場合、各因子はフェルミオンの一体演算子ではなく同一サイト上の
局所スピン行列要素を表すため、\ :math:`i_p=j_p`\ を満たす必要があります。
因子は演算子積で用いる順に記述します。以下にファイル例を記載します。

::

    ==========================
    NNBodyInterAll 3
    ==========================
    ========NBodyInterAll=====
    ==========================
    1 0 0 1 0 -1.0 0.0
    1 1 0 0 0 -1.0 0.0
    2 0 0 0 0 1 1 1 1 0.5 0.0

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降:
   [int02] ([int03] [int04] [int05] [int06]) ... [double01] [double02]

   括弧で示した4整数の組を [int02] 回繰り返します。

パラメータ
^^^^^^^^^^

-  :math:`[`\ string01\ :math:`]`

   **形式 :** string型 (空白不可)

   **説明 :** N体相互作用項の総数のキーワード名を指定します(任意)。

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** N体相互作用項の総数を指定します。

-  :math:`[`\ int02\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** この項に含まれる因子数 :math:`N` を指定します。
   正の整数でなければなりません。

-  各因子の :math:`[`\ int03\ :math:`]`, :math:`[`\ int05\ :math:`]`

   **形式 :** int型 (空白不可)

   **説明 :** サイト番号を指定する整数。0以上\ ``Nsite``\ 未満で指定します。
   :math:`[`\ int03\ :math:`]`\ が\ :math:`i_p`\ 、
   :math:`[`\ int05\ :math:`]`\ が\ :math:`j_p`\ に対応します。

-  各因子の :math:`[`\ int04\ :math:`]`, :math:`[`\ int06\ :math:`]`

   **形式 :** int型 (空白不可)

   | **説明 :** スピンまたは局所状態を指定する整数。
   | Hubbard、tJ、Kondo およびそれらの GC/NConserved 版では
   | 0: アップスピン
   | 1: ダウンスピン
   | を表します。
   | Spin/SpinGC では OneBodyG/TwoBodyG と同じ局所状態インデックスを用います。
   | SpinlessFermion/SpinlessFermionGC では 0 のみ指定可能です
     (0以外を指定した場合はエラー終了します)。

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型 (空白不可)

   **説明 :** 係数 :math:`V` の実部を指定します。

-  :math:`[`\ double02\ :math:`]`

   **形式 :** double型 (空白不可)

   **説明 :** 係数 :math:`V` の虚部を指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  :math:`[`\ int01\ :math:`]`\ と定義されているN体相互作用項の総数が異なる場合はエラー終了します。

-  サイト番号またはスピン番号が許容範囲外の場合はエラー終了します。

-  Spin/SpinGC では、すべての因子が\ :math:`i_p=j_p`\ を満たす必要があります。

-  Kondo/KondoGC/KondoNConserved では、局在スピンサイトを含む因子は
   :math:`i_p=j_p`\ を満たす必要があります。一般スピンのKondo局在スピンには対応していません。

-  対角 NBodyInterAll 項の係数の虚部は 0 でなければなりません。

-  ハミルトニアンがエルミートであるため、非対角 NBodyInterAll 項は隣接するエルミート共役ペアとして記述する必要があります。
   フェルミオン模型では、因子の順序を逆にし、各因子の生成・消滅側のインデックスを入れ替え、係数を複素共役にした行を続けて入力します。
   Spin/SpinGC でも、因子の順序を逆にし、各因子の局所状態インデックスを入れ替え、係数を複素共役にした厳密なエルミート共役の行を入力できます。
   すべての因子が異なるサイトに作用する場合は、スピン演算子が可換であるため、因子の順序は保ったまま局所状態インデックスを入れ替えた同等な行も許容されます。

   例えば、フェルミオン模型では以下の隣接する2行がエルミート共役ペアになります。

   ::

       2 0 0 1 0 2 1 3 1 0.25  0.10
       2 3 1 2 1 1 0 0 0 0.25 -0.10

   対応する演算子は以下の通りです。

   .. math::

      \begin{aligned}
      &(0.25+0.10\mathrm{i})
      c_{0,0}^{\dagger} c_{1,0}
      c_{2,1}^{\dagger} c_{3,1},\\
      &(0.25-0.10\mathrm{i})
      c_{3,1}^{\dagger} c_{2,1}
      c_{1,0}^{\dagger} c_{0,0}.
      \end{aligned}

   Spin/SpinGC でも、厳密なエルミート共役は同じように記述できます。

   ::

       2 0 1 0 0 2 0 2 1 0.25  0.10
       2 2 1 2 0 0 0 0 1 0.25 -0.10

   :math:`X_i^{ab}=|a\rangle_i\langle b|` とすると、対応する演算子は以下の通りです。
   スピン1/2の場合、同じペアは :math:`S^\pm` を用いても書けます。

   .. math::

      \begin{aligned}
      &(0.25+0.10\mathrm{i}) X_0^{1,0} X_2^{0,1}
        =(0.25+0.10\mathrm{i}) S_0^+ S_2^-,\\
      &(0.25-0.10\mathrm{i}) X_2^{1,0} X_0^{0,1}
        =(0.25-0.10\mathrm{i}) S_2^+ S_0^-.
      \end{aligned}

   この例では2つの因子が異なるサイトに作用するため、2行目は1行目と同じサイト順のままでも記述できます。

   ::

       2 0 1 0 0 2 0 2 1 0.25  0.10
       2 0 0 0 1 2 1 2 0 0.25 -0.10

   これは同じエルミート共役ペアに対応します。

   .. math::

      \begin{aligned}
      &(0.25-0.10\mathrm{i}) X_0^{0,1} X_2^{1,0}
        =(0.25-0.10\mathrm{i}) S_0^- S_2^+
        =(0.25-0.10\mathrm{i}) S_2^+ S_0^-.
      \end{aligned}

-  NBodyInterAll は SpinGC, Spin, HubbardGC, Hubbard, HubbardNConserved,
   SpinlessFermionGC, SpinlessFermion, tJGC, tJ, tJNConserved,
   KondoGC, Kondo, KondoNConserved で使用できます。

-  NBodyInterAll は TimeEvolution では使用できません。

-  NBodyInterAll は HubbardNConserved, tJNConserved, KondoNConserved の
   CalcSpec では使用できません。

-  NBodyInterAll は SpinlessFermion/SpinlessFermionGC の FullDiag では使用できません。

.. raw:: latex

   \newpage
