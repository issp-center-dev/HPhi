.. highlight:: none

熱的純粋量子状態による有限温度計算
----------------------------------

杉浦・清水によって、 少数個（サイズが大きい場合はほぼ一つ）の
波動関数から有限温度の物理量を計算する方法が提案されました [1]_ [2]_ 。
その状態は熱的純粋量子状態(TPQ)と呼ばれています。
TPQはハミルトニアンを波動関数に順次作用させて得られるので、
Lanczos法の技術がそのまま使うことができます。
TPQ状態は次のように与えられます。

.. math::

   |\Phi(\beta)\rangle\equiv\exp(-\frac{\beta\hat{\mathcal H}}{2})|\Phi_{\rm rand}\rangle.

ここで :math:`\beta` は逆温度、:math:`|\Phi_{\rm rand}\rangle`  はランダムな初期ベクトルです。
:math:`|\Phi(\beta)\rangle` の期待値として、有限温度の物理量が計算できることが
示されています.

実際の実装について
~~~~~~~~~~~~~~~~~~

mTPQ状態の構成について
^^^^^^^^^^^^^^^^^^^^^^^^^^
ここでは、とくに計算が簡単な, micro canonical TPQ(mTPQ)の
概要を述べます [1]_。

:math:`|\Phi_{\rm rand}\rangle`\ をあるランダムベクトルとします。
これに\ :math:`(l-{\mathcal H}/N_{s})`\ (:math:`l`\ はある定数、\ :math:`N_{s}`\ はサイト数)を\ :math:`k`\ 回作用させた
（規格化された）ベクトルは次のように与えられます。

.. math::

   \begin{aligned}
   |\Phi_{k}\rangle \equiv \frac{(l-{\mathcal H}/N_{s})|\Phi_{k-1}\rangle}{|(l-{\mathcal H}/N_{s})|\Phi_{k-1}\rangle|}.\end{aligned}

この\ :math:`|\Phi_{k}\rangle`\ がmTPQ状態で、このmTPQ状態に対応する逆温度\ :math:`\beta_{k}`\ は
以下のように内部エネルギー\ :math:`u_{k}`\ から求めることができます。

.. math::

   \begin{aligned}
   \beta_{k}\sim \frac{2k/N_{s}}{l-u_{k}},~~
   u_{k} = \langle \Phi_{k}|{\mathcal H}|\Phi_{k}\rangle/N_{s}.\end{aligned}

そして、任意 [3]_ の物理量\ :math:`\hat{A}`\ の\ :math:`\beta_{k}`\ での平均値は

.. math::

   \begin{aligned}
   \langle \hat{A}\rangle_{\beta_{k}} =  \langle \Phi_{k}|\hat{A}|\Phi_{k}\rangle/N_{s}\end{aligned}

となります。 有限系では最初の乱数ベクトルによる誤差がありますので、
いくつか独立な計算を行って、\ :math:`|\psi_{0}\rangle`
に関する平均値および標準偏差を見積もっています。

cTPQ状態の構成について
^^^^^^^^^^^^^^^^^^^^^^^^^^
カノニカルTPQ(cTPQ)状態の構成方法について述べます [2]_,
cTPQ法では :math:`\exp[-\beta\hat{\mathcal H}/2]` は
次にように近似されます.

.. math::

  &\exp(-\frac{\beta\hat{\mathcal H}}{2})|\Phi_{\rm rand}\rangle\sim|\Phi_{k}\rangle = U(\Delta\tau)^{k}|\Phi_{\rm rand}\rangle\\
  &U(\Delta\tau)=\sum_{n=0}^{n_{\rm max}}\frac{1}{n!}(-\frac{\Delta\tau}{2}\hat{\mathcal H})^{n}\\
  &\beta_{k}=k\Delta \tau
   
mTPQと同じように物理量はcTPQ状態の期待値として計算できます。

.. math::

   \langle \hat{A}\rangle_{\beta_{k}} =  \frac{\langle\Phi_{k}|\hat{A}|\Phi_{k}\rangle}{\langle\Phi_{k}|\Phi_{k}\rangle}.



初期ベクトルの設定について
^^^^^^^^^^^^^^^^^^^^^^^^^^

熱的純粋量子状態による有限温度計算では、初期ベクトルは全ての成分に対してランダムな係数を与えます。
初期ベクトルの係数の型はModParaで指定される入力ファイルの\ ``InitialVecType``\ を用い、
実数もしくは複素数の指定をすることができます。乱数のシードは\ ``initial_iv``
(:math:`\equiv r_s`)により

.. math::

   \begin{aligned}
   123432+(n_{\rm run}+1)\times  |r_s|+k_{\rm Thread}+N_{\rm Thread} \times k_{\rm Process}\end{aligned}


で与えられます。ここで、\ :math:`n_{\rm run}`\ はrunの回数であり、runの最大回数はスタンダードモード用入力ファイル、
もしくはModParaで指定される入力ファイルの\ ``NumAve``\ で指定します。
``initial_iv``\ はスタンダードモード用の入力ファイル、もしくはエキスパートモードではModParaで指定される入力ファイルで指定します。乱数はSIMD-oriented
Fast Mersenne Twister(dSFMT)を用い発生させています [4]_ 。
また、\ :math:`k_{\rm Thread}, N_{\rm Thread}, k_{\rm Process}`\ 
はそれぞれスレッド番号、スレッド数、プロセス番号を表します。
したがって同じ\ ``initial_iv``\ を用いても、並列数が異なる場合には別の初期波動関数が生成されます。

.. [1] \S. Sugiura, A. Shimizu, Phys. Rev. Lett. **108**, 240401 (2012).
.. [2] \S. Sugiura and A. Shimizu, Phys. Rev. Lett. **111**, 010401 (2013).
.. [3] 局所的にという条件がつきます。
.. [4] \http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/SFMT.
