リスタートと状態import
=======================

HPhiには独立した3つの機能があります。``ReStart`` は保存したalgorithm状態から
solverを継続しますが、後述のlegacy経路には例外があります。``InputEigenVec=1`` はCGでの評価、または新しい時刻を0から
始めるTimeEvolutionの初期状態としてvectorをimportします。``CalcSpec`` 4/5は
``TMComponents`` と ``recalcvec`` によりLanczos spectrum列を継続します。異なる
機能はvectorのlayoutが一致する場合でも契約が異なります。

v3 solver checkpointとsector v1/v2状態形式は、v3.7.2までのreleaseに使われた
コードには含まれません。これらのreleaseでは、ここで説明するCG/TPQ/TimeEvolutionの
vectorはlegacy raw形式です。file形式のversion 3とHPhiのrelease番号は別のものです。

solverの継続
------------

``CalcMod`` の ``ReStart`` は、0が継続なし、1が新規計算とcheckpoint出力、
2がcheckpoint読込と更新、3が読込のみです。v3 manifest/vector形式はraw基底と
``TransSym`` 基底のCG、mTPQ、cTPQ、TimeEvolutionで使われます。保存時と現在の
MPI process数は異なってもかまいません。HPhiは保存済みglobal連続区間と現在の
所有区間の交差を読み、vectorを再配置します。sector v3継続ではmetadata layoutも
``replicated`` (0) と ``distributed`` (1) の間で変更できます。rank所有は変わっても
global vector順序は変わりません。

``ReStart=3`` は物理量出力を読込専用にする指定ではありません。継続時には時系列を
復旧して追記します。またTE継続checkpointが選択されない場合、明示的な
``InputEigenVec=1`` / ``SpectrumVec`` からimportへ進むことがあり、後述のとおり
既存の時系列を作り直す場合があります。

v3継続の ``Lanczos_max`` は絶対的な最終上限です。CGは保存した反復label、
TPQ/cTPQとTimeEvolutionは保存した次stepから再開します。CGでは上限が保存反復label
より小さい場合、TPQ/cTPQとTimeEvolutionでは上限が保存した次step以下の場合、
復元状態を評価し、時系列行を追加せず、checkpointも書き直しません。CGは
``Lanczos_max`` と同じlabelの反復を含み、TPQ/cTPQとTimeEvolutionはそのstepを含みません。
収束終端CG checkpointは現在の ``LanczosEps`` で再判定します。厳しい許容誤差では
反復を再開できますが、上限によるno-opを現在も収束しているとは表示しません。

``CheckpointInterval`` はCG、mTPQ、cTPQの周期checkpointを指定し、0は終了時だけを
意味します。TimeEvolutionは継続checkpointに ``OutputInterval`` を使い、時刻表は
保存された次stepを含む必要があります。詳細は
:doc:`restartファイル仕様 <../filespecification/outputfiles_ja/restart_ja>` を参照してください。

検証範囲と制限
--------------

v3は保存state・基底・solver単位を検証しますが、Hの一致は要求しません。
HAMILTONIAN語は情報用で、rawでは0、sectorではparsed Hのdigestを保存します。
``hamiltonian_changed=no`` も伝播履歴の検証済みという意味ではありません。
時間依存sector TE入力では、保存stepのHと現在のbase Hをscheduleとして比較できないため
``unknown`` と表示します。継続には同じHと伝播設定を用い、別Hや新protocolには
状態importを使ってください。

個別の検査と制限は次のとおりです。

* CGはstate・探索方向・Ritz値・反復labelを復元します。Hx/Hpは再計算しますが、
  Hのidentityは照合しません。
* mTPQはHと ``LargeValue`` を照合しません。cTPQは保存した虚時間刻みと現在の
  ``1/LargeValue`` を照合しますが、Hと ``ExpandCoef`` は照合しません。
  ``InvTemp`` と ``ReStart=2/3`` の組合せは拒否します。
* TEは保存境界時刻と、表が保存next stepを覆うことを照合します。過去・将来のH/時刻行や
  Taylor次数は照合しません。table/sector TEの保存時計は ``time[next_step-1]``、
  raw Laser TEは ``Tinit + next_step*TimeSlice`` です。境界一致だけではprotocol一致とはいえません。
* TPQ sampleは保存stateとset labelで識別します。初期乱数設定と並列数は初期stateの
  生成条件であり、継続時にrank数を変えてもそのstateを生成し直しません。

状態import
----------

v3とsector v1/v2の状態importでは模型、basis identity、sector、phase、payload長、norm、digestを検証します。
source method/state/step/timeとsolver scalarは来歴であり、source solverを再開しません。
CGはimport状態を評価し、TimeEvolutionは現在の時刻表のstep 0から開始します。v3 importは
保存時と現在のprocess数が異なる場合も自動的に再配置します。sector v3と旧sector
v1/v2 importはmetadata layout 0/1の変更も許します。v1/v2では保存時の分割に従って
payload checksumを検証してから再配置し、分割非依存のsector digestで基底集合を
検証します。global順序はbasis builderが代表状態をcanonicalにsortする規則で決まるため、
集合が同じなら順序も同じです。必要なcopy手順は出力ファイル仕様に示します。

CGの ``InputEigenVec=1`` はimport状態の評価であり、その状態からCG反復を開始する
指定ではありません。importでは元Hと現在Hの一致は要求しませんが、現在Hが選択sectorで
有効である必要はあります。

raw CG/TEの ``OutputEigenVec`` とcTPQの ``InvTemp`` 固有vector出力は、現在も
legacy raw layoutです。通常のraw CG→TE importもこの形式を使い、v3 solver manifestを
使いません。このimportは引き続き可能ですが、模型・基底identityとpayload checksumを
fileから検証することはできません。

TE継続checkpointが選択されず、``SpectrumVec`` がv3またはsector v1/v2の場合は、
``ReStart=2/3`` でもstep 0からimportします。SS/Norm/Flct/集約Green時系列を作り直す前に
警告します。意図的に新protocolを開始する場合は、別の出力directoryと ``ReStart=0/1`` を
使ってください。legacy rawの ``SpectrumVec`` は、``ReStart=2/3`` では保存next stepから
継続し、``ReStart=0/1`` ではstep 0からimportします。legacy raw TEには比較できる保存時計がありません。

MPI process数を変更する手順
-----------------------------

要素を持たないownerのfileも含め、全source rank fileを残します。v3継続ではmanifestと
そこから参照される全vector fileをcopyし、希望するprocess数で実行します。追加parameterは
不要です。v3継続は各manifest行のprocess数を使い、 ``RestartSourceRanks`` を無視します。
このため、TPQ/cTPQの移行中断後に、新しいprocess数で書いたv3 sampleと移行前のlegacy
sampleが混在していても再開できます。このparameterは残ったlegacy sampleだけに適用され、
全sampleがv3として再出力された後は削除できます。

v3またはsector v1/v2状態importでは、全source fileをimport用prefixへcopyし、希望する
process数で実行します。source数はheaderから得られます。状態importでの
``RestartSourceRanks`` は任意で、指定値の一致を追加検査するだけです。

legacy raw importまたは継続では、保存済みrank別fileを全てcopyし、process数を変更する前に
``ModPara`` へ ``RestartSourceRanks=<保存process数>`` を追加します。末尾の空rank fileが
欠けた状態と、より少ないrankで完結したfile一式を区別できないため、HPhiはfile名から
この値を推測しません。

legacyとspectrumの挙動
----------------------

legacy raw ``tmpvec_set*_rank_*.dat`` はCGとTPQ/cTPQで引き続き読めます。
この形式はsource rank数もidentity語も保存しないため、rank数を変更する読込では
``RestartSourceRanks`` が必須です。要素数0の末尾rankを含む全保存fileを残してください。
HPhiは再配置前にnative header、file長、step、global次元を検証しますが、fileだけから
模型やsectorのidentityを検証することはできません。legacy CGはstateだけを復元し、
探索方向を0にし、現在HでRitz値を再計算します。従来の追加反復の意味を維持し、
反復labelは1から始まり、以前の ``Lanczos_Step`` prefixを復旧しません。
入力stateの正規化・直交化は既に済んでいる前提で、legacy readerはやり直しません。
この動作はwarm start相当ですが、v3のX/P継続や一般的なwarm start用interfaceではありません。
legacy TPQ/cTPQとraw TimeEvolutionは、可能な場合に検証済みの次stepより
前の出力行を復旧します。旧sector v1/v2は状態import形式で、solver継続checkpointでは
ありません。

legacy mTPQは初期normを復元できず0として出力します。legacy cTPQはβを
``(step-1)/LargeValue`` から再構成するため、``LargeValue`` が変わっていない前提です。
いずれもHのidentityは照合しません。

基底状態Lanczosと ``CalcSpec`` 4/5は、従来の ``TMComponents`` / ``recalcvec`` と
反復回数の意味を維持し、v3 solver manifestを使いません。

出力の復旧
----------

追記前に全checkpointを検証し、復旧境界以後の完全な行と未完の末尾行を削除します。
途中の不正な行はエラーです。履歴ファイルが無い場合は新規作成され、その後の行だけが
書かれます。失われたprefixは再構成できません。set別fileとCG/TimeEvolution履歴は
中断なし実行と同じ行順です。共有TPQ aggregateは実行順なので、複数setの継続結果は
``(set, step, observable)`` keyでsortして比較してください。

``CalcTimer.dat`` の ``Restart write`` (7100)、``Restart read`` (7110)、
``Restart recover outputs`` (7120)は、それぞれv3公開、継続vector読込、時系列の
検証・書換えを計測します。これは機能別の計測値であり、管理されたstorage benchmark
ではありません。
