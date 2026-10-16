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
``TransSym`` 基底のCG、mTPQ、cTPQ、TimeEvolutionで使われます。このreleaseでは
保存時と現在のMPI process数が同じでなければなりません。sector継続ではbasis
layoutも同じである必要があります。rank再配置とlayout変更は拒否します。

v3継続の ``Lanczos_max`` は絶対的な最終上限です。CGは保存した反復label、
TPQ/cTPQとTimeEvolutionは保存した次stepから再開します。保存stepが小さくした
上限以上なら、復元状態を評価し、時系列行を追加せず、checkpointも書き直しません。
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
現在、同じprocess数を必要とし、sector importは同じlayoutも必要です。必要なcopy手順は
出力ファイル仕様に示します。

legacyとspectrumの挙動
----------------------

legacy raw ``tmpvec_set*_rank_*.dat`` はCGとTPQ/cTPQで引き続き読めます。legacy CGは
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
