.. highlight:: none

実時間発展
----------

:math:`{\cal H}\Phi`\ では初期波動関数を\ :math:`|\Phi(t_0)\rangle`\ として、

.. math:: |\Phi (t_n)\rangle = e^{-i {\cal H}  \Delta t_n}|\Phi (t_{n-1})\rangle


の関係を用いて逐次時間発展の計算をしています。
ここで時刻 :math:`t_n = \sum_{j=1}^n  \Delta t_j` です。
実際の計算では\ :math:`e^{-i {\cal H}  \Delta t_n}`\ を

.. math:: e^{-i {\cal H}  \Delta t_n} =\sum_{l=0}^m \frac{1}{l!}(-i {\cal H}  \Delta t_n)^l


と近似しています。展開次数の上限\ :math:`m`\ は\ ``ModPara``\ ファイルの\ ``ExpandCoef``\ で指定することができます。
各ステップの状態は正規化され、``Norm.dat`` には正規化前のノルムが出力されます。
この値と ``ExpandCoef`` を増やした際の収束を確認して、Taylor打切り誤差を評価してください。
時間に依存しないHamiltonianでは、エネルギー保存も確認できます。

時間依存Hamiltonian
~~~~~~~~~~~~~~~~~~~~

各ステップでは、一次項を含む **全てのTaylor項に同じHamiltonian** を使います。
時刻の規約は次のとおりです。

* ``TEOneBody`` / ``TETwoBody``: 行0は入力状態を伝播せず記録します。
  行 :math:`n>0` は :math:`\Delta t_n=t_n-t_{n-1}` の区間を
  :math:`H(t_n)`（右端点）で伝播し、時刻 :math:`t_n` の状態・物理量を出力します。
  同じ時刻が連続する場合は状態を伝播しません。
* raw基底のExpert mode ``Laser``: step :math:`n` は
  :math:`H(T_{\mathrm{init}}+n\Delta t)` を使い、step 0から ``TimeSlice`` だけ伝播します。
  出力時刻は :math:`T_{\mathrm{init}}+(n+1)\Delta t` です。
  エネルギーはそのステップの伝播に使ったHamiltonianで評価します。
* 対称性sectorの ``Laser``: 行0は ``Tinit`` の入力状態を記録し、
  以降は右端点規約を使います。詳細は ``TransSym`` ファイルの仕様を参照してください。

Standard modeの ``PumpType`` が ``TEOneBody`` / ``TETwoBody`` を生成する場合は、
laser駆動を含めて時刻表の規約に従います。raw Expert modeの ``Laser`` とは区別してください。

``ExpandCoef`` を増やすと、各区間で固定したHの指数関数の近似が改善します。
時間離散化の精度次数が上がるわけではありません。
滑らかな時間依存Hに対する端点評価は一般に一次精度なので、最終時刻を固定して
時間刻みを半分ずつにした収束も確認してください。
例えば :math:`H(t)=(1+t)K` の厳密な伝播は :math:`\exp[-iK(T+T^2/2)]` であり、
十分高いTaylor次数でも右端点評価の誤差は時間刻みに比例します。

旧バージョンのTE結果との互換性
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

今回のH更新修正より前（3.6.0、3.7.2を含む）のraw TEは、一次項に前ステップの
:math:`H\psi`、高次項に現在のHを使っていました。
raw Expert mode ``Laser`` の初回では、この一次項自体が欠ける場合もありました。
修正によって時間依存系の状態・物理量が変わります。
静的Hの時刻表による伝播と、上記の時刻ラベルの規約は変わりません。

checkpoint形式は変更していません。再開時も、保存状態と選択したHamiltonianから
:math:`H\psi` を作り直します。ただし旧binaryのcheckpointを読んでも、保存状態に既に
蓄積された誤差は取り消せません。影響する軌道について修正済みの結果が必要な場合は、
元の初期状態から再計算してください。
