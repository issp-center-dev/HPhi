.. _Subsec:restart_v3:

v3 solver restartファイル
=========================

CG、mTPQ、cTPQ、TimeEvolutionは ``ReStart`` 継続にこの形式を使います。
保存時と現在のMPI process数は異なってもかまいません。HPhiは保存されたglobal連続区間と
現在の所有区間の交差を読み、vectorを再配置します。``TransSym`` v3継続ではglobal
vector順序を変えずにmetadata layoutも ``replicated`` (0) と ``distributed`` (1) の
間で変更できます。独立して公開するTPQ set間では保存rank数・layoutが異なる場合が
ありますが、同じset内の全rank fileは一致が必要です。旧sector v1/v2状態importも
rank数とlayoutの変更に対応します。

名前、単位、manifest
--------------------

vector名は ``<CDataFileHead>_restart_set<k>_comp<c>_gen<g>_rank_<r>.dat`` です。
``k`` はCG状態、TPQ sample、またはTimeEvolutionの0、component 0は状態、CGの
component 1は探索方向です。generationは0または1です。公開単位は、CGでは全状態の
両component、TPQ/cTPQでは1 sample、TimeEvolutionではset 0です。

``<CDataFileHead>_restart_manifest.dat`` だけが公開済みdataを指します::

   hphi_restart_manifest 1
   method <CG|TPQ|cTPQ|TimeEvolution>
   basis <raw|sector>
   set <k> ranks <P> gen <g> step <STEP> ncomp <1|2>

CG行はrank数、generation、stepが共通で、TPQ行は独立です。``ReStart=3`` は
これらのfileを変更しません。新しい ``ReStart=1`` 公開は無関係な旧行を置換します。
writerは全 ``.part`` fileを作り、完成単位をrenameした後でmanifestを原子的に公開します。
古いfileは公開後だけ削除するため、manifestが不完全な単位を指すことはありません。

binary header
-------------

各vectorはlittle-endian unsigned 64-bitの40語headerと、``LOCAL_DIM`` 個のIEEE
binary64実部・虚部の組です。vector添字0は保存しません。

.. list-table:: v3 header語
   :header-rows: 1
   :widths: 7 22 71

   * - index
     - 名称
     - 意味
   * - 0
     - MAGIC
     - ``0x0a31565349485048`` (``HPHISV1\n``)
   * - 1
     - VERSION
     - ``3``
   * - 2
     - PHASE
     - basis位相規約。importでも検査
   * - 3
     - SCALAR
     - binary64 complex payloadでは ``128``
   * - 4
     - MODEL
     - 正規化した入力模型identity
   * - 5
     - NSITE
     - 保存basisのsite数
   * - 6
     - NUP
     - 固定up数。非固定なら0
   * - 7
     - NDOWN
     - 固定down数。非固定なら0
   * - 8
     - NE
     - 固定粒子数。非固定なら0
   * - 9
     - FULL_DIM
     - sectorのfull dimension。rawは0
   * - 10
     - DIM
     - global vector dimension
   * - 11
     - RANKS
     - 保存時MPI process数。manifest行と一致
   * - 12
     - RANK
     - 保存rank index
   * - 13
     - LAYOUT
     - rawは2、sectorは0 (replicated) または1 (distributed)
   * - 14
     - OFFSET
     - payloadの0始まりglobal先頭位置
   * - 15
     - LOCAL_DIM
     - payload要素数
   * - 16
     - GROUP
     - sector group digest。rawは0
   * - 17
     - SECTOR_COUNT
     - 分割非依存sector count
   * - 18
     - SECTOR_XOR
     - 分割非依存sector XOR digest
   * - 19
     - SECTOR_SUM
     - 分割非依存sector sum digest
   * - 20
     - ORDER
     - sector global順序XOR digest
   * - 21
     - HAMILTONIAN
     - 情報用Hamiltonian digest。rawは0/unknown
   * - 22
     - METHOD
     - source ``CalcType``
   * - 23
     - STATE
     - set index。file名とmanifestに一致
   * - 24
     - STEP
     - 次の反復/step。FLAGS指定時は終端CG行
   * - 25
     - TIME
     - TimeEvolutionの物理時刻bit。他は0
   * - 26
     - DATA_XOR
     - 分割非依存payload XOR digest
   * - 27
     - DATA_SUM
     - 分割非依存payload sum digest
   * - 28
     - LOCAL_SITE_MASK
     - Kondo局在site mask。他は0
   * - 29
     - FIXED_FLAGS
     - raw固定量またはsector Kondo flags
   * - 30
     - BASIS_KIND
     - 0 raw、1 sector
   * - 31
     - COMPONENT
     - 0 state、1 CG direction
   * - 32
     - TOTAL2SZ
     - 固定raw ``2Sz`` のbit表現。非固定なら0
   * - 33
     - SITE_DIGEST
     - raw site/local-spin identity digest
   * - 34
     - ORDER_SUM
     - sector global順序sum digest
   * - 35
     - NSET
     - CG状態数、source TPQ sample数、TEでは1
   * - 36
     - SCALAR_A
     - CG固有値またはTPQ first norm
   * - 37
     - SCALAR_B
     - cTPQ累積逆温度。他は0
   * - 38
     - SCALAR_C
     - cTPQ ``delta_tau``。他は0
   * - 39
     - FLAGS
     - bit 0は収束終端CG checkpoint。他bitは不正

比較単位とdigest
----------------

同じset/componentのrank間では ``RANK``、``OFFSET``、``LOCAL_DIM`` だけが
異なります。CGの状態・探索方向間では ``COMPONENT`` とdata digestも異なって
構いません。CGの状態間では ``STATE``、data digest、``SCALAR_A`` が異なり、
step、flags、set数、generation、rank数は共通です。TPQ sampleはstate、step、
data/scalar、manifest行を独立に持ちますが、basis identityとmethodは一致します。

data digestは各global indexと実部・虚部のbit列をhashし、rank間でXORとunsigned sumを
集約します。sector順序digestはglobal index、代表状態、orbit size、stabilizer sizeを
hashします。この2種類の集約により、分割に依存せず破損と並べ替えを検出します。

rank再配置
----------

v3とsector v1/v2 fileは、保存rank数、offset、local dimensionを記録します。
保存rank 0から保存rank数−1まで、payloadが0のfileも含めて全fileを残してください。
HPhiは保存file一式を検証し、各保存区間と現在区間の交差部分を読みます。
自己記述形式なので ``RestartSourceRanks`` は不要ですが、指定する場合は保存rank数と
一致する必要があります。

sector v1/v2では、保存時の分割に依存するDATA checksumを保存fileから最初に再現し、
検証に成功してからpayloadを再配置します。分割非依存sector digestでglobal基底順序も
検証します。source/destination layoutは0または1のどちらでも構いませんが、同じ
checkpointのsource fileは全て同じ有効なlayoutを宣言する必要があります。

legacy raw fileはnative 12-byteの ``STEP`` / local dimension headerだけを持ち、
保存rank数もidentityも記録しません。このためrank数を変えるlegacy読込では
``RestartSourceRanks=<保存rank数>`` が必須です。宣言した全fileの存在、header、長さ、
共通step、次元和を検証し、模型・sector identityではなく次元だけを検査したことを
表示します。

rank変更専用の回帰試験coverageは次のとおりです。

.. list-table:: rank変更restart/import回帰試験表
   :header-rows: 1
   :widths: 18 18 20 44

   * - 形式
     - 用途
     - 基底・手法
     - 試験済みsourceからdestination
   * - v3
     - 継続
     - raw CG
     - 4→1、4→16、16→4、1→4 rank
   * - v3
     - 継続
     - raw mTPQ/cTPQ/TE
     - mTPQ 4→2、cTPQ 4→2、TE 4→1 rank
   * - v3
     - 継続
     - sector CG
     - 3 distributed→6 replicated・16 distributed rank（空ownerを含む）
   * - v3
     - 状態import
     - raw CG/TE
     - 4→1 rank
   * - sector v1/v2
     - 状態import
     - Spin/Kondo
     - Spin v1 4→1/2/16、Kondo v2 4→1 rank。layout変更とchecksum破損も検査
   * - legacy raw
     - 状態import
     - CG
     - 4→1・4→16 rank。空source rank fileの欠落拒否も検査
   * - legacy raw
     - 継続
     - mTPQ
     - 4→2 rank

sector v3のmTPQ/cTPQ/TimeEvolution継続と状態import、legacy rawのCG/cTPQ/
TimeEvolution継続とTE状態importには、rank変更専用の数値回帰試験がまだありません。
これらも同じ区間readerを使いますが、専用caseが無いことを追加の数値検証と解釈しないで
ください。

復旧とlegacy選択
----------------

必要な各setについて、最初にmanifest行を使います。参照されたv3 fileの欠損・破損は
致命的で、legacyへfallbackしません。行が無い場合だけ完全なlegacy fileを検証して
使います。両方無ければ ``ReStart=2`` は新規setを開始し、``ReStart=3`` はエラーです。
sector v1/v2はimport専用です。legacy CGは追加反復をlabel 1から始め、以前のCG履歴prefixを
復旧しません。legacyのrank数変更には上記の ``RestartSourceRanks`` が必要です。

履歴復旧は ``STEP`` より前の行を保持し、それ以後とEOFの未完断片を削除し、途中の破損を
拒否します。終端CG行は現在の収束判定も満たす場合だけ保持します。TPQ aggregate行は
実行順に追記されるので、複数setの継続は ``(set, step, observable)`` でsortして比較します。

v3 vectorを状態importに使う方法
--------------------------------

manifestからsetとgenerationを選び、全rankの
``..._restart_set<k>_comp0_gen<g>_rank_<r>.dat`` を既存のimport名へcopyします。
TimeEvolutionでは ``<SpectrumVec>_rank_<r>.dat``、CGでは
``<CDataFileHead>_eigenvec_<state>_rank_<r>.dat`` です。全source rank fileを保持し、
``RestartSourceRanks`` を指定する場合はその数にします。v3では指定は任意ですが、指定値は
保存rank数と一致する必要があります。現在のrank数と、sectorではlayoutも変更できます。
importは新しいsolver時計を始め、継続ではありません。
