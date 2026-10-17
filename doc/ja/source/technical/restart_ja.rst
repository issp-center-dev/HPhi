リスタートと状態import
=======================

HPhiには独立した3つの機能があります。``ReStart`` は保存したalgorithm状態から
solverを継続します。``InputEigenVec=1`` はCGでの評価、または新しい時刻を0から
始めるTimeEvolutionの初期状態としてvectorをimportします。``CalcSpec`` 4/5は
``TMComponents`` と ``recalcvec`` によりLanczos spectrum列を継続します。異なる
機能のファイルを互いに置き換えて使用することはできません。

solverの継続
------------

``CalcMod`` の ``ReStart`` は、0が継続なし、1が新規計算とcheckpoint出力、
2がcheckpoint読込と更新、3が読込のみです。v3 manifest/vector形式はraw基底と
``TransSym`` 基底のCG、mTPQ、cTPQ、TimeEvolutionで使われます。保存時と現在の
MPI process数は異なってもかまいません。HPhiは保存済みglobal連続区間と現在の
所有区間の交差を読み、vectorを再配置します。sector v3継続ではmetadata layoutも
``replicated`` (0) と ``distributed`` (1) の間で変更できます。rank所有は変わっても
global vector順序は変わりません。

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

状態import
----------

状態importでは模型、basis identity、sector、phase、payload長、norm、digestを検証します。
source method/state/step/timeとsolver scalarは来歴であり、source solverを再開しません。
CGはimport状態を評価し、TimeEvolutionは現在の時刻表のstep 0から開始します。v3 importは
保存時と現在のprocess数が異なる場合も自動的に再配置します。sector v3と旧sector
v1/v2 importはmetadata layout 0/1の変更も許します。v1/v2では保存時の分割に従って
payload checksumを検証してから再配置し、分割非依存のsector digestでglobal基底順序を
検証します。必要なcopy手順は出力ファイル仕様に示します。

MPI process数を変更する手順
-----------------------------

要素を持たないownerのfileも含め、全source rank fileを残します。v3継続ではmanifestと
そこから参照される全vector fileをcopyし、希望するprocess数で実行します。追加parameterは
不要です。v3またはsector v1/v2状態importでは、全source fileをimport用prefixへcopyし、
希望するprocess数で実行します。source数はheaderから得られます。これらの形式での
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
模型やsectorのidentityを検証することはできません。legacy CGは
従来の追加反復の意味を維持し、反復labelは1から始まり、以前の ``Lanczos_Step`` prefixを
復旧しません。legacy TPQ/cTPQとraw TimeEvolutionは、可能な場合に検証済みの次stepより
前の出力行を復旧します。旧sector v1/v2は状態import形式で、solver継続checkpointでは
ありません。

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
