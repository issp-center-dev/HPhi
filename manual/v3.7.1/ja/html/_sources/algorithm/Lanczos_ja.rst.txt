.. highlight:: none

Lanczos法
==============

手法概要
-------------------------

Lanczos法の解説については TITPACK [#]_ のマニュアルと 線形計算の数理 [#]_ を参考にしています.

Lanczos法はある初期ベクトルにハミルトニアン
を作用させて最大・最小近傍の固有値・固有ベクトルを
求める方法です。
Lanczos法で固有値を求める際にはHilbert空間の次元の
大きさの波動関数を表すベクトルが2つ [#]_ あれば原理的には実行可能なため、
大規模疎行列の対角化手法として有用であることが知られています。
後述するように、固有ベクトルを求める際にはHilbert空間の次元の大きさのベクトルがもう1本必要です。

Lanczos法の 原理はべき乗法に基づいています。
べき乗法ではある任意のベクトル\ :math:`\vec{x}_{0}`\ に
Hamitonianを逐次的に作用させて, :math:`{\mathcal H}^{n}\vec{x}_{0}`
を作成します。 このとき、生成される空間
:math:`\mathcal{K}_{n+1}({\mathcal H},\vec{x}_{0})=\{\vec{x}_{0},{\mathcal H}^{1}\vec{x}_{0},\dots,{\mathcal H}^{n}\vec{x}_{0}\}`
はKrylov部分空間といわれます。
初期ベクトルを\ :math:`{\mathcal H}`\ の固有ベクトル\ :math:`\vec{e}_{i}`\ (対応する固有値を :math:`E_{i}`\ とする)
の重ね合わせで表すと

.. math::

   \begin{aligned}
   \vec{x}_{0}=\sum_{i}a_{i}\vec{e}_{i}\end{aligned}

となります。 ここで、\ :math:`E_{0}`\ を絶対値最大の固有値としました。
またハミルトニアンはエルミートであるため、
固有値は全て実数であることに注意する必要があります。
これにHamiltonianの\ :math:`{\mathcal H}^{n}`\ を作用させると、

.. math::

   \begin{aligned}
   {\mathcal H}^{n}\vec{x}_{0}=E_{0}^{n}\Big[ a_{0}\vec{e}_{0}+\sum_{i\neq0}\left(\frac{E_{i}}{E_{0}}\right)^na_{i}\vec{e}_{i}\Big]\end{aligned}

となり、絶対値最大固有値\ :math:`E_{0}`\ に対応する固有ベクトルが支配的になります。
適切な変換を行って、この固有ベクトルの成分を抽出するのが Lanczos法です。

Lanczos法では、 :math:`\mathcal{K}_{n}({\mathcal H},\vec{x}_{0})`
から正規直交ベクトル\ :math:`{\vec{v}_{0},\dots,\vec{v}_{n-1}}`\ を次の手続きにしたがって
順次生成していきます。 初期条件を
:math:`\vec{v}_{0} =\vec{x}_{0}/|\vec{x}_{0}|`,
:math:`\beta_{0}=0,\vec{x}_{-1}=0`
とすると、正規直交基底は次の手続きによって逐次的に生成することができます。

.. math::

   \begin{aligned}
   \alpha_{k} &= ({\mathcal H}\vec{v}_{k},\vec{v}_{k}) \\
   \vec{w}   &= {\mathcal H}\vec{v}_{k}-\beta_{k}\vec{v}_{k-1}-\alpha_{k}\vec{v}_{k} \\
   \beta_{k+1} &= |\vec{w}| \\
   \vec{v}_{k+1} &= \frac{\vec{v}_{k}}{|\vec{v}_{k}|}\end{aligned}

この定義から、\ :math:`\alpha_{k}`,\ :math:`\beta_{k}` ともに実数であることに注意する必要があります。

これらの正規直交基底が成すKryrov部分空間の中で
もとのハミルトニアンに対する固有値問題は、

.. math::

   \begin{aligned}
   T_{n}=V_{n}^{\dagger}{\mathcal H} V_{n}\end{aligned}

と変形されます。 ここで、
:math:`V_{n}`\ は\ :math:`\vec{v}_{i}(i=0,1,\dots,n-1)`\ を
並べた行列です。 :math:`T_{n}`\ は三重対角行列であり、 その対角成分は
:math:`\alpha_{i}`, 副対角成分は\ :math:`\beta_{i}`\ で与えられます。
この三重対角行列\ :math:`T_{n}`\ の
固有値はもとのハミルトニアンの固有値の
近似値となっています(\ :math:`V^{\dagger}V=I`,\ :math:`I`\ は単位行列であることに注意)。
:math:`T_{n}`\ の固有ベクトルを\ :math:`\tilde{\vec{e}}_{i}`
とするともとのハミルトニアンの固有ベクトルとの関係は
:math:`\vec{e}_{i}=V\tilde{\vec{e}}_{i}`\ で与えられます。
:math:`V`\ を覚えていれば、Lanczos法を行うと同時に固有ベクトル
を求めることができますが、実際の場合は
(Hilbert空間の次元 :math:`\times` Lanczos法の反復回数)の
大きさの行列を保持することは不可能です。
そこで、Lanczos法で固有値を求めた後、
:math:`\tilde{\vec{e}_{i}}`\ を保存しておき  [4]_、
再び同じLanczos法の計算を行い
求めた\ :math:`V`\ から元の固有ベクトルを再現 するようにしています。

Lanczos法では, 元のHilbert空間の次元より十分小さい反復回数
で最大及び最小に近い固有値を精度よく求めることが
できることが知られています。
すなわち\ :math:`T_{n}`\ の次元\ :math:`n`\ は数百-数千程度ですみます。
定量的には、最大及び最小固有値の評価の誤差は
Lanczosの反復回数に対して指数関数的に減少することが
示されています(詳細はRef. [2]_ を参照して下さい)。

逆反復法
-------------------------

固有値の近似値がわかっているときは
適当なベクトル\ :math:`\vec{y}_{0}`\ に対して
:math:`({\mathcal H}-E_{n})^{-1}`\ を 逐次的に作用させれば、
固有値\ :math:`E_{n}`\ に対応する固有ベクトルの 成分が支配的になり、
固有ベクトルを精度良く求めることができます。

:math:`({\mathcal H}-E_{n})^{-1}`\ を作用させる方程式を書き換えると
以下の連立方程式が得られます。

.. math::

   \begin{aligned}
   \vec{y}_{k}&=({\mathcal H}-E_{n})\vec{y}_{k+1}\end{aligned}

この連立方程式を 共役勾配法(CG法)を用いて解くことで、
固有ベクトルを求めることができます。 その固有ベクトルから
固有値およびその他の相関関数を 求めることができます。 ただし,
CG法の実行にはヒルベルト空間の次元のベクトル を4本確保する必要があり、
大規模系の計算を実行する際にはメモリが足りなくなる恐れがあるので注意が必要です。

実際の実装について
-------------------------

初期ベクトルの設定について
^^^^^^^^^^^^^^^^^^^^^^^^^^

Lanczos法では、スタンダードモード用の入力ファイル、もしくはエキスパートモードではModParaで指定する\ ``initial_iv``\ (:math:`\equiv r_s`)により初期ベクトルの設定方法を指定します。また、初期ベクトルはModParaで指定される入力ファイルの\ ``InitialVecType``\ を用い、実数もしくは複素数の指定をすることができます。

-  カノニカル集団かつ ``initial_iv`` :math:`\geq 0`\ の場合

   ヒルベルト空間のうち、

   .. math::

      \begin{aligned}
      (N_{\rm dim}/2 + r_s ) \% N_{\rm dim}\end{aligned}


   の成分が与えられます。ここで、 :math:`N_{\rm dim}` は対象となるヒルベルト空間の総数で、:math:`N_{\rm dim}/2` はデフォルト値\ ``initial_iv``
   :math:`=1`\ で特殊なヒルベルト空間の選択をさけるために加えられています。なお、選択された成分の係数は実数の場合は\ :math:`1`\ 、複素数の場合には\ :math:`(1+i)/\sqrt{2}`\ が与えられます。

-  グランドカノニカル集団 もしくは ``initial_iv`` :math:`< 0`\ の場合

   初期ベクトルはランダムベクトルとして与えられます。乱数のシードは

   .. math::

      \begin{aligned}
      123432+|r_s|\end{aligned}


   で指定します。ここで、\ :math:`n_{\rm run}`\ はrunの回数であり、runの最大回数はスタンダードモード用入力ファイル、もしくはModParaで指定される入力ファイルの\ ``NumAve``\ で指定します。\ ``initial_iv``\ は入力ファイルで指定します。乱数はSIMD-oriented
   Fast Mersenne Twister
   (dSFMT)を用い発生させています [5]_ 。

収束判定について
^^^^^^^^^^^^^^^^

:math:`{\mathcal H}\Phi` では、\ :math:`T_{n}`\ の対角化にlapackのルーチン
:math:`\rm dsyev`\ を使用しており、
:math:`T_{n}`\ の基底状態の次の固有値（第一励起状態のエネルギー）
を収束判定条件に用いています。 デフォルトの設定では、
最初の5回のLanczosステップの後に、
２回毎に\ :math:`T_{n}`\ の対角化を行い、
前のLanczosステップの第一励起状態のエネルギーと
指定した精度以内で一致すれば、収束したと判定しています。
なお、収束する際の精度は ``CDataFileHead``  (エキスパートモードではModParaファイル内)で指定することが可能です。

その後、Lanczos法を再度行い、 逐次\ :math:`V`\ を求めて、指定した
準位の固有ベクトルを求めます。
得られた固有ベクトル\ :math:`|n\rangle`\ を用い、
エネルギーの期待値\ :math:`E_{n}=\langle n|{\mathcal H}|n\rangle` 
およびバリアンス\ :math:`\Delta=\langle n|{\mathcal H}^{2}|n\rangle -(\langle n|{\mathcal H}|n\rangle)^2`
を求めて、\ :math:`E_{n}`\ がLaczos法で求めた固有値と
指定した精度で一致しているか、
バリアンスが指定した精度以下になっているかを チェックしています。
指定した精度に達していれば、対角化を終了しています。

指定した精度に達していない場合には
逆反復法を用いて再度固有ベクトルを求め直します。
逆反復法の初期ベクトルとしてLanczos法で求めた
固有ベクトルをとった方が一般に収束が早いので、
標準の設定ではそのように取っています。


.. [#] \http://www.qa.iir.titech.ac.jp/~nishimori/titpack2_new/index-e.html
.. [#] \M. Sugihara, K. Murota, Theoretical Numerical Linear Algebra, Iwanami Stud-ies in Advanced Mathematics, Iwanami Shoten, Publishers, 2009.
.. [#] 高速化のために、\In :math:`{\mathcal H}\Phi` ではハミルトニアンの対角成分を表すベクトル1本と,スピン :math:`z` 成分 :math:`S_{z}` 保存, 粒子数保存の場合はその状態を指定するベクトル1本を余計に確保しています。いずれのベクトルの大きさもHilbert空間の次元です。 
.. [#] :math:`\tilde{\vec{e}_{i}}`\ の次元は高々Lanczos法の反復回数であることに注意。
.. [#] \http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/index.html
